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


def quadratic_products():
    from sympy.integrals.rubi.constraints import cons47, cons2, cons3, cons8, cons227, cons5, cons228, cons130, cons229, cons230, cons13, cons165, cons231, cons139, cons232, cons233, cons234, cons235, cons236, cons237, cons70, cons71, cons49, cons238, cons29, cons50, cons19, cons239, cons240, cons241, cons242, cons68, cons243, cons244, cons245, cons246, cons148, cons247, cons248, cons249, cons250, cons251, cons252, cons253, cons254, cons168, cons255, cons33, cons170, cons256, cons257, cons96, cons149, cons258, cons40, cons259, cons260, cons43, cons20, cons261, cons262, cons263, cons264, cons265, cons266, cons267, cons268, cons269, cons45, cons270, cons56, cons271, cons272, cons273, cons274, cons275, cons276, cons277, cons278, cons279, cons280, cons281, cons282, cons283, cons284, cons285, cons286, cons287, cons288, cons21, cons289, cons290, cons291, cons292, cons293, cons294, cons295, cons296, cons297, cons298, cons299, cons300, cons301, cons302, cons303, cons304, cons305, cons306, cons307, cons308, cons309, cons310, cons311, cons312, cons313, cons314, cons315, cons86, cons87, cons316, cons317, cons318, cons127, cons210, cons319, cons320, cons321, cons64, cons322, cons323, cons324, cons325, cons326, cons4, cons327, cons328, cons329, cons141, cons330, cons331, cons332, cons333, cons152, cons334, cons150, cons335, cons198, cons336, cons337, cons338, cons339, cons340, cons91, cons341, cons342, cons343, cons90, cons89, cons344, cons345, cons346, cons128, cons347, cons348, cons209, cons349, cons350, cons351, cons352, cons353, cons354, cons355, cons356, cons357, cons358, cons359, cons360, cons361, cons362, cons363, cons364, cons365, cons366, cons367, cons368, cons369, cons370, cons371, cons372, cons373, cons374, cons375, cons376, cons377, cons151, cons378, cons126, cons379, cons95, cons25, cons167, cons75, cons380, cons82, cons381, cons382, cons383, cons384, cons385, cons386, cons387, cons52, cons388, cons389, cons390, cons391, cons392, cons393, cons394, cons395, cons396, cons397, cons398, cons399, cons400, cons401, cons402, cons403, cons404, cons405, cons406, cons407, cons408, cons409, cons410, cons411, cons412, cons413, cons414, cons415, cons416, cons417, cons418, cons211, cons419, cons420, cons421, cons422, cons423, cons424, cons425, cons426, cons427, cons428, cons429, cons430, cons431, cons432, cons433, cons222, cons434, cons435, cons436, cons437, cons438, cons439, cons440, cons441, cons442, cons443, cons444, cons445, cons446, cons447, cons448, cons449, cons450, cons451, cons452, cons453, cons454, cons455, cons226, cons36, cons37, cons38, cons456, cons457, cons458, cons459, cons460


    pattern192 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons47)
    rule192 = ReplacementRule(pattern192, replacement192)

    pattern193 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons47, cons227)
    rule193 = ReplacementRule(pattern193, replacement193)

    pattern194 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons228, cons130, cons229)
    rule194 = ReplacementRule(pattern194, With194)

    pattern195 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons228, cons130, cons230)
    rule195 = ReplacementRule(pattern195, replacement195)

    pattern196 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons228, cons13, cons165, cons231)
    rule196 = ReplacementRule(pattern196, replacement196)

    pattern197 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(-3)/2), x_), cons2, cons3, cons8, cons228)
    rule197 = ReplacementRule(pattern197, replacement197)

    pattern198 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons228, cons13, cons139, cons232, cons231)
    rule198 = ReplacementRule(pattern198, replacement198)

    pattern199 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons228, cons233, cons229)
    rule199 = ReplacementRule(pattern199, With199)

    pattern200 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons228, CustomConstraint(With200))
    rule200 = ReplacementRule(pattern200, replacement200)

    pattern201 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons228)
    rule201 = ReplacementRule(pattern201, replacement201)

    pattern202 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons5, cons234)
    rule202 = ReplacementRule(pattern202, replacement202)

    pattern203 = Pattern(Integral(S(1)/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons3, cons8, cons235)
    rule203 = ReplacementRule(pattern203, replacement203)

    pattern204 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons228)
    rule204 = ReplacementRule(pattern204, replacement204)

    pattern205 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons3, cons8, cons13, cons236)
    rule205 = ReplacementRule(pattern205, replacement205)

    pattern206 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons228, cons13, CustomConstraint(With206))
    rule206 = ReplacementRule(pattern206, replacement206)

    pattern207 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons5, cons228, cons237)
    rule207 = ReplacementRule(pattern207, With207)

    pattern208 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons5, cons70, cons71)
    rule208 = ReplacementRule(pattern208, replacement208)

    pattern209 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons47, cons49, cons238)
    rule209 = ReplacementRule(pattern209, replacement209)

    pattern210 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons47, cons49, cons239)
    rule210 = ReplacementRule(pattern210, replacement210)

    pattern211 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons47, cons49, cons240)
    rule211 = ReplacementRule(pattern211, replacement211)

    pattern212 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons47, cons241, cons242, cons68)
    rule212 = ReplacementRule(pattern212, replacement212)

    pattern213 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons47, cons241)
    rule213 = ReplacementRule(pattern213, replacement213)

    pattern214 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons47, cons241, cons243)
    rule214 = ReplacementRule(pattern214, replacement214)

    pattern215 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons47, cons241)
    rule215 = ReplacementRule(pattern215, replacement215)

    pattern216 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons47, cons241, cons244, cons243)
    rule216 = ReplacementRule(pattern216, replacement216)

    pattern217 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons47, cons241, cons245)
    rule217 = ReplacementRule(pattern217, replacement217)

    pattern218 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons47, cons241, cons246, cons148, cons247, cons248)
    rule218 = ReplacementRule(pattern218, replacement218)

    pattern219 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons47, cons241, cons246, cons148, cons249, cons248, cons250)
    rule219 = ReplacementRule(pattern219, replacement219)

    pattern220 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons47, cons241, cons13, cons165, cons251, cons240, cons250, cons252, cons253, cons248)
    rule220 = ReplacementRule(pattern220, replacement220)

    pattern221 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons47, cons241, cons246, cons139, cons254, cons248)
    rule221 = ReplacementRule(pattern221, replacement221)

    pattern222 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons47, cons241, cons246, cons139, cons168, cons248)
    rule222 = ReplacementRule(pattern222, replacement222)

    pattern223 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons47, cons241, cons246, cons139, cons255, cons248)
    rule223 = ReplacementRule(pattern223, replacement223)

    pattern224 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons47, cons241, cons33, cons170, cons240, cons256, cons257)
    rule224 = ReplacementRule(pattern224, replacement224)

    pattern225 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons47, cons241, cons33, cons96, cons248)
    rule225 = ReplacementRule(pattern225, replacement225)

    pattern226 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons47, cons149, cons241)
    rule226 = ReplacementRule(pattern226, replacement226)

    pattern227 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons228, cons258, cons40)
    rule227 = ReplacementRule(pattern227, replacement227)

    pattern228 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons5, cons259, cons260)
    rule228 = ReplacementRule(pattern228, replacement228)

    pattern229 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons258, cons149, cons43)
    rule229 = ReplacementRule(pattern229, replacement229)

    pattern230 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons19, cons5, cons259, cons149, cons43)
    rule230 = ReplacementRule(pattern230, replacement230)

    pattern231 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons258, cons149, cons242)
    rule231 = ReplacementRule(pattern231, replacement231)

    pattern232 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons5, cons259, cons149, cons242)
    rule232 = ReplacementRule(pattern232, replacement232)

    pattern233 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons228, cons258, cons149, cons13, cons139)
    rule233 = ReplacementRule(pattern233, replacement233)

    pattern234 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**S(2), x_), cons2, cons8, cons29, cons50, cons5, cons259, cons149, cons13, cons139)
    rule234 = ReplacementRule(pattern234, replacement234)

    pattern235 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons228, cons258, cons149, cons20, cons13, cons261, cons262, cons263)
    rule235 = ReplacementRule(pattern235, replacement235)

    pattern236 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons19, cons5, cons259, cons149, cons20, cons13, cons261, cons262, cons263)
    rule236 = ReplacementRule(pattern236, replacement236)

    pattern237 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons258, cons149, cons264)
    rule237 = ReplacementRule(pattern237, replacement237)

    pattern238 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons5, cons259, cons149, cons264)
    rule238 = ReplacementRule(pattern238, replacement238)

    pattern239 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons258, cons149, cons265)
    rule239 = ReplacementRule(pattern239, replacement239)

    pattern240 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons19, cons5, cons259, cons149, cons265)
    rule240 = ReplacementRule(pattern240, replacement240)

    pattern241 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons258)
    rule241 = ReplacementRule(pattern241, replacement241)

    pattern242 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons259)
    rule242 = ReplacementRule(pattern242, replacement242)

    pattern243 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons228, cons258, cons246, cons165, cons266, cons255, cons248)
    rule243 = ReplacementRule(pattern243, replacement243)

    pattern244 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons259, cons246, cons165, cons266, cons255, cons248)
    rule244 = ReplacementRule(pattern244, replacement244)

    pattern245 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons228, cons258, cons246, cons165, cons267, cons240, cons248)
    rule245 = ReplacementRule(pattern245, replacement245)

    pattern246 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons259, cons246, cons165, cons267, cons240, cons248)
    rule246 = ReplacementRule(pattern246, replacement246)

    pattern247 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons228, cons258, cons246, cons139, cons254, cons248)
    rule247 = ReplacementRule(pattern247, replacement247)

    pattern248 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons259, cons246, cons139, cons254, cons248)
    rule248 = ReplacementRule(pattern248, replacement248)

    pattern249 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons228, cons258, cons246, cons139, cons168, cons248)
    rule249 = ReplacementRule(pattern249, replacement249)

    pattern250 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons259, cons246, cons139, cons168, cons248)
    rule250 = ReplacementRule(pattern250, replacement250)

    pattern251 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons228, cons258, cons33, cons268, cons240, cons248)
    rule251 = ReplacementRule(pattern251, replacement251)

    pattern252 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons5, cons259, cons33, cons268, cons240, cons248)
    rule252 = ReplacementRule(pattern252, replacement252)

    pattern253 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons228, cons258, cons33, cons269, cons255, cons248)
    rule253 = ReplacementRule(pattern253, replacement253)

    pattern254 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons5, cons259, cons33, cons269, cons255, cons248)
    rule254 = ReplacementRule(pattern254, replacement254)

    pattern255 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons3, cons8, cons50, cons19, cons149)
    rule255 = ReplacementRule(pattern255, replacement255)

    pattern256 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons5, cons259, cons149, cons45, cons270)
    rule256 = ReplacementRule(pattern256, replacement256)

    pattern257 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons228, cons258, cons149)
    rule257 = ReplacementRule(pattern257, replacement257)

    pattern258 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons259, cons149)
    rule258 = ReplacementRule(pattern258, replacement258)

    pattern259 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons49)
    rule259 = ReplacementRule(pattern259, replacement259)

    pattern260 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons49, cons244, cons56)
    rule260 = ReplacementRule(pattern260, replacement260)

    pattern261 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons228, cons49, cons130, cons271)
    rule261 = ReplacementRule(pattern261, replacement261)

    pattern262 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons49, cons272, cons246, cons165, cons96, cons273, cons248)
    rule262 = ReplacementRule(pattern262, replacement262)

    pattern263 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons228, cons49, cons272, cons13, cons165, cons274, cons275, cons33, cons248)
    rule263 = ReplacementRule(pattern263, replacement263)

    pattern264 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons228, cons49, cons272, cons246, cons139, cons168, cons248)
    rule264 = ReplacementRule(pattern264, replacement264)

    pattern265 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons228, cons49, cons272, cons13, cons139, cons276, cons33, cons248)
    rule265 = ReplacementRule(pattern265, replacement265)

    pattern266 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons49)
    rule266 = ReplacementRule(pattern266, replacement266)

    pattern267 = Pattern(Integral(S(1)/(sqrt(d_ + x_*WC('e', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons49, cons277)
    rule267 = ReplacementRule(pattern267, replacement267)

    pattern268 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons49, cons277)
    rule268 = ReplacementRule(pattern268, replacement268)

    pattern269 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons49, cons278)
    rule269 = ReplacementRule(pattern269, replacement269)

    pattern270 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons228, cons49, cons272, cons33, cons168, cons240, cons279)
    rule270 = ReplacementRule(pattern270, replacement270)

    pattern271 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons228, cons49, cons272, cons33, cons96, cons280)
    rule271 = ReplacementRule(pattern271, replacement271)

    pattern272 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons49)
    rule272 = ReplacementRule(pattern272, replacement272)

    pattern273 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons228, cons281, cons241, cons130)
    rule273 = ReplacementRule(pattern273, replacement273)

    pattern274 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons282, cons130, cons283)
    rule274 = ReplacementRule(pattern274, replacement274)

    pattern275 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241, cons284)
    rule275 = ReplacementRule(pattern275, With275)

    pattern276 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons282, cons285)
    rule276 = ReplacementRule(pattern276, With276)

    pattern277 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241, cons286)
    rule277 = ReplacementRule(pattern277, replacement277)

    pattern278 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons282, cons287)
    rule278 = ReplacementRule(pattern278, replacement278)

    pattern279 = Pattern(Integral(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241)
    rule279 = ReplacementRule(pattern279, replacement279)

    pattern280 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons282)
    rule280 = ReplacementRule(pattern280, replacement280)

    pattern281 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241, cons20, cons168, cons288)
    rule281 = ReplacementRule(pattern281, replacement281)

    pattern282 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons282, cons20, cons168, cons288)
    rule282 = ReplacementRule(pattern282, replacement282)

    pattern283 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241, cons33, cons168)
    rule283 = ReplacementRule(pattern283, replacement283)

    pattern284 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons282, cons33, cons168)
    rule284 = ReplacementRule(pattern284, replacement284)

    pattern285 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241)
    rule285 = ReplacementRule(pattern285, replacement285)

    pattern286 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons282)
    rule286 = ReplacementRule(pattern286, replacement286)

    pattern287 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241)
    rule287 = ReplacementRule(pattern287, replacement287)

    pattern288 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons282)
    rule288 = ReplacementRule(pattern288, replacement288)

    pattern289 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons228, cons281, cons241, cons33, cons96)
    rule289 = ReplacementRule(pattern289, replacement289)

    pattern290 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons19, cons282, cons33, cons96)
    rule290 = ReplacementRule(pattern290, replacement290)

    pattern291 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons228, cons281, cons241, cons21)
    rule291 = ReplacementRule(pattern291, replacement291)

    pattern292 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons19, cons282, cons21)
    rule292 = ReplacementRule(pattern292, replacement292)

    pattern293 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241)
    rule293 = ReplacementRule(pattern293, replacement293)

    pattern294 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2), x_), cons2, cons8, cons29, cons50, cons282)
    rule294 = ReplacementRule(pattern294, replacement294)

    pattern295 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241, cons13, cons139, cons232)
    rule295 = ReplacementRule(pattern295, replacement295)

    pattern296 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons282, cons13, cons139, cons232)
    rule296 = ReplacementRule(pattern296, replacement296)

    pattern297 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons228, cons281, cons241, cons289)
    rule297 = ReplacementRule(pattern297, replacement297)

    pattern298 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons5, cons282, cons289)
    rule298 = ReplacementRule(pattern298, replacement298)

    pattern299 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons290, cons291, cons292, cons149)
    rule299 = ReplacementRule(pattern299, replacement299)

    pattern300 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons3, cons8, cons29, cons50, cons293, cons241, cons33, cons294, cons295, cons296)
    rule300 = ReplacementRule(pattern300, replacement300)

    pattern301 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons3, cons8, cons29, cons50, cons293, cons241, cons33, cons294)
    rule301 = ReplacementRule(pattern301, replacement301)

    pattern302 = Pattern(Integral(x_**m_/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons228, cons297)
    rule302 = ReplacementRule(pattern302, replacement302)

    pattern303 = Pattern(Integral((e_*x_)**m_/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons50, cons228, cons297)
    rule303 = ReplacementRule(pattern303, replacement303)

    pattern304 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241, cons297)
    rule304 = ReplacementRule(pattern304, replacement304)

    pattern305 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/sqrt(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons282, cons297)
    rule305 = ReplacementRule(pattern305, replacement305)

    pattern306 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241, cons246, cons298, cons165)
    rule306 = ReplacementRule(pattern306, replacement306)

    pattern307 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons282, cons246, cons298, cons165)
    rule307 = ReplacementRule(pattern307, replacement307)

    pattern308 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241, cons246, cons298, cons139)
    rule308 = ReplacementRule(pattern308, replacement308)

    pattern309 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons282, cons246, cons298, cons139)
    rule309 = ReplacementRule(pattern309, replacement309)

    pattern310 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons241)
    rule310 = ReplacementRule(pattern310, replacement310)

    pattern311 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons299)
    rule311 = ReplacementRule(pattern311, replacement311)

    pattern312 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons281, cons241, cons149, cons242)
    rule312 = ReplacementRule(pattern312, replacement312)

    pattern313 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons5, cons282, cons149, cons242)
    rule313 = ReplacementRule(pattern313, replacement313)

    pattern314 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons281, cons241, cons244, cons13, cons139)
    rule314 = ReplacementRule(pattern314, replacement314)

    pattern315 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons19, cons5, cons282, cons244, cons13, cons139)
    rule315 = ReplacementRule(pattern315, replacement315)

    pattern316 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons281, cons241, cons244)
    rule316 = ReplacementRule(pattern316, replacement316)

    pattern317 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons19, cons5, cons282, cons244)
    rule317 = ReplacementRule(pattern317, replacement317)

    pattern318 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons228, cons281, cons241, cons13, cons165, cons300, cons68, cons301, cons302)
    rule318 = ReplacementRule(pattern318, replacement318)

    pattern319 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons19, cons282, cons13, cons165, cons300, cons68, cons301, cons303)
    rule319 = ReplacementRule(pattern319, replacement319)

    pattern320 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons228, cons281, cons241, cons13, cons165, cons240, cons304, cons305, cons302)
    rule320 = ReplacementRule(pattern320, replacement320)

    pattern321 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons19, cons282, cons13, cons165, cons240, cons304, cons305, cons303)
    rule321 = ReplacementRule(pattern321, replacement321)

    pattern322 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241, cons246, cons139, cons170, cons306, cons302)
    rule322 = ReplacementRule(pattern322, replacement322)

    pattern323 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons282, cons246, cons139, cons170, cons306, cons303)
    rule323 = ReplacementRule(pattern323, replacement323)

    pattern324 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons241, cons246, cons139, cons168, cons302)
    rule324 = ReplacementRule(pattern324, replacement324)

    pattern325 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons282, cons246, cons139, cons168, cons303)
    rule325 = ReplacementRule(pattern325, replacement325)

    pattern326 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons228, cons281, cons241, cons13, cons139, cons302)
    rule326 = ReplacementRule(pattern326, replacement326)

    pattern327 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons19, cons282, cons13, cons139, cons303)
    rule327 = ReplacementRule(pattern327, replacement327)

    pattern328 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons281, cons241, cons307, cons240, cons302)
    rule328 = ReplacementRule(pattern328, replacement328)

    pattern329 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons19, cons5, cons282, cons307, cons240, cons303)
    rule329 = ReplacementRule(pattern329, replacement329)

    pattern330 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons281, cons241, cons308)
    rule330 = ReplacementRule(pattern330, replacement330)

    pattern331 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons8, cons29, cons50, cons19, cons5, cons282, cons309)
    rule331 = ReplacementRule(pattern331, replacement331)

    pattern332 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons2, cons3, cons8, cons29, cons50, cons241, cons310, cons311)
    rule332 = ReplacementRule(pattern332, With332)

    pattern333 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons312)
    rule333 = ReplacementRule(pattern333, With333)

    pattern334 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons2, cons3, cons8, cons29, cons50, cons241, cons310, cons313)
    rule334 = ReplacementRule(pattern334, With334)

    pattern335 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons314)
    rule335 = ReplacementRule(pattern335, With335)

    pattern336 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/4)*(d_ + x_*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons282)
    rule336 = ReplacementRule(pattern336, replacement336)

    pattern337 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/4)*(d_ + x_*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons282)
    rule337 = ReplacementRule(pattern337, replacement337)

    pattern338 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons234, cons231)
    rule338 = ReplacementRule(pattern338, replacement338)

    pattern339 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons315, cons231)
    rule339 = ReplacementRule(pattern339, replacement339)

    pattern340 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons5, cons282, cons149, cons45, cons295)
    rule340 = ReplacementRule(pattern340, replacement340)

    pattern341 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons228, cons281, cons241, cons149, cons86)
    rule341 = ReplacementRule(pattern341, With341)

    pattern342 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons5, cons282, cons149, cons86)
    rule342 = ReplacementRule(pattern342, With342)

    pattern343 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons228, cons281, cons241, cons149)
    rule343 = ReplacementRule(pattern343, With343)

    pattern344 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons5, cons282, cons149)
    rule344 = ReplacementRule(pattern344, With344)

    pattern345 = Pattern(Integral((u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons70, cons71)
    rule345 = ReplacementRule(pattern345, replacement345)

    pattern346 = Pattern(Integral((a_ + u_**S(2)*WC('c', S(1)))**WC('p', S(1))*(u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons5, cons70, cons71)
    rule346 = ReplacementRule(pattern346, replacement346)

    pattern347 = Pattern(Integral(x_**WC('n', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons5, cons87, cons316)
    rule347 = ReplacementRule(pattern347, replacement347)

    pattern348 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons317, cons47, cons318)
    rule348 = ReplacementRule(pattern348, replacement348)

    pattern349 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons47, cons318, cons149, cons244)
    rule349 = ReplacementRule(pattern349, replacement349)

    pattern350 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons47, cons318, cons149, cons246, cons139, cons170)
    rule350 = ReplacementRule(pattern350, replacement350)

    pattern351 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons317, cons47, cons318, cons149, cons13, cons139, cons319)
    rule351 = ReplacementRule(pattern351, replacement351)

    pattern352 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons47, cons318, cons149, cons33, cons96, cons227, cons320)
    rule352 = ReplacementRule(pattern352, replacement352)

    pattern353 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons47, cons318, cons149, cons33, cons96, cons321)
    rule353 = ReplacementRule(pattern353, replacement353)

    pattern354 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons47, cons318, cons149, cons64, cons321, cons322)
    rule354 = ReplacementRule(pattern354, replacement354)

    pattern355 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons47, cons318, cons149, cons321)
    rule355 = ReplacementRule(pattern355, replacement355)

    pattern356 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons47, cons323, cons149, cons241, cons13, cons324)
    rule356 = ReplacementRule(pattern356, replacement356)

    pattern357 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons47, cons323, cons149, cons325)
    rule357 = ReplacementRule(pattern357, replacement357)

    pattern358 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons47, cons323, cons149, cons241, cons244)
    rule358 = ReplacementRule(pattern358, replacement358)

    pattern359 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons47, cons323, cons149, cons241, cons321, cons272, cons33, cons96)
    rule359 = ReplacementRule(pattern359, replacement359)

    pattern360 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons47, cons323, cons149, cons241, cons321, cons272, cons274, cons326)
    rule360 = ReplacementRule(pattern360, replacement360)

    pattern361 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons317, cons47, cons149)
    rule361 = ReplacementRule(pattern361, replacement361)

    pattern362 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons317, cons228, cons258, cons40)
    rule362 = ReplacementRule(pattern362, replacement362)

    pattern363 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons317, cons259, cons327)
    rule363 = ReplacementRule(pattern363, replacement363)

    pattern364 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons258, cons86, cons316)
    rule364 = ReplacementRule(pattern364, replacement364)

    pattern365 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons259, cons86, cons316)
    rule365 = ReplacementRule(pattern365, replacement365)

    pattern366 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons258, cons328)
    rule366 = ReplacementRule(pattern366, replacement366)

    pattern367 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons259, cons329)
    rule367 = ReplacementRule(pattern367, replacement367)

    pattern368 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons258, cons246, cons139, cons170)
    rule368 = ReplacementRule(pattern368, replacement368)

    pattern369 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons259, cons246, cons139, cons170)
    rule369 = ReplacementRule(pattern369, replacement369)

    pattern370 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons258, cons141, cons330, cons56)
    rule370 = ReplacementRule(pattern370, replacement370)

    pattern371 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons259, cons141, cons330, cons56)
    rule371 = ReplacementRule(pattern371, replacement371)

    pattern372 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons258, cons331, cons255)
    rule372 = ReplacementRule(pattern372, replacement372)

    pattern373 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons259, cons331, cons255)
    rule373 = ReplacementRule(pattern373, replacement373)

    pattern374 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons258, cons321)
    rule374 = ReplacementRule(pattern374, replacement374)

    pattern375 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons259, cons321)
    rule375 = ReplacementRule(pattern375, replacement375)

    pattern376 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons2, cons8, cons127, cons210, cons332, cons13, cons333)
    rule376 = ReplacementRule(pattern376, replacement376)

    pattern377 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons2, cons8, cons127, cons210, cons5, cons332)
    rule377 = ReplacementRule(pattern377, replacement377)

    pattern378 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons258, cons149, cons152, cons13, cons334)
    rule378 = ReplacementRule(pattern378, replacement378)

    pattern379 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons259, cons149, cons152, cons13, cons334)
    rule379 = ReplacementRule(pattern379, replacement379)

    pattern380 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons258, cons149, cons150, cons335)
    rule380 = ReplacementRule(pattern380, replacement380)

    pattern381 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons259, cons149, cons150, cons335)
    rule381 = ReplacementRule(pattern381, replacement381)

    pattern382 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons258, cons149, cons198, cons335)
    rule382 = ReplacementRule(pattern382, replacement382)

    pattern383 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons259, cons149, cons198, cons335)
    rule383 = ReplacementRule(pattern383, replacement383)

    pattern384 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons317, cons228, cons258, cons149, cons43, cons336, cons337)
    rule384 = ReplacementRule(pattern384, replacement384)

    pattern385 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons317, cons259, cons149, cons43, cons338, cons337)
    rule385 = ReplacementRule(pattern385, replacement385)

    pattern386 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons317, cons228, cons258, cons149, cons43, cons339)
    rule386 = ReplacementRule(pattern386, replacement386)

    pattern387 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons317, cons259, cons149, cons43, cons339)
    rule387 = ReplacementRule(pattern387, replacement387)

    pattern388 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons258, cons149, cons43, cons340, cons165, cons91, cons341)
    rule388 = ReplacementRule(pattern388, replacement388)

    pattern389 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons259, cons149, cons43, cons340, cons165, cons91, cons341)
    rule389 = ReplacementRule(pattern389, replacement389)

    pattern390 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons317, cons228, cons258, cons149, cons43, cons340, cons165, cons337, cons342, cons343)
    rule390 = ReplacementRule(pattern390, replacement390)

    pattern391 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons4, cons317, cons259, cons149, cons43, cons340, cons165, cons337, cons342, cons343)
    rule391 = ReplacementRule(pattern391, replacement391)

    pattern392 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons258, cons149, cons43, cons340, cons139, cons90)
    rule392 = ReplacementRule(pattern392, replacement392)

    pattern393 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons259, cons149, cons43, cons340, cons139, cons90)
    rule393 = ReplacementRule(pattern393, replacement393)

    pattern394 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons317, cons228, cons258, cons149, cons43, cons340, cons139)
    rule394 = ReplacementRule(pattern394, replacement394)

    pattern395 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons4, cons317, cons259, cons149, cons43, cons340, cons139)
    rule395 = ReplacementRule(pattern395, replacement395)

    pattern396 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons258, cons149, cons43, cons89, cons90, cons337, cons344)
    rule396 = ReplacementRule(pattern396, replacement396)

    pattern397 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons259, cons149, cons43, cons89, cons90, cons337, cons344)
    rule397 = ReplacementRule(pattern397, replacement397)

    pattern398 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons258, cons149, cons43, cons89, cons91, cons248)
    rule398 = ReplacementRule(pattern398, replacement398)

    pattern399 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons259, cons149, cons43, cons89, cons91, cons248)
    rule399 = ReplacementRule(pattern399, replacement399)

    pattern400 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/((x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons258)
    rule400 = ReplacementRule(pattern400, replacement400)

    pattern401 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons259)
    rule401 = ReplacementRule(pattern401, replacement401)

    pattern402 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons317, cons228, cons258, cons149, cons345, cons346, cons128)
    rule402 = ReplacementRule(pattern402, replacement402)

    pattern403 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons317, cons259, cons149, cons345, cons347, cons128)
    rule403 = ReplacementRule(pattern403, replacement403)

    pattern404 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons258, cons149, cons345, cons89, cons91, cons248)
    rule404 = ReplacementRule(pattern404, replacement404)

    pattern405 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons259, cons149, cons345, cons89, cons91, cons248)
    rule405 = ReplacementRule(pattern405, replacement405)

    pattern406 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons317, cons228, cons258, cons149, cons345, cons348, cons248)
    rule406 = ReplacementRule(pattern406, replacement406)

    pattern407 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons317, cons259, cons149, cons345, cons348, cons248)
    rule407 = ReplacementRule(pattern407, replacement407)

    pattern408 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons317, cons228, cons258, cons149, cons209)
    rule408 = ReplacementRule(pattern408, replacement408)

    pattern409 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons317, cons259, cons349, cons152, cons350, cons351)
    rule409 = ReplacementRule(pattern409, replacement409)

    pattern410 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons317, cons259, cons149, cons209)
    rule410 = ReplacementRule(pattern410, replacement410)

    pattern411 = Pattern(Integral(x_**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons228, cons258)
    rule411 = ReplacementRule(pattern411, replacement411)

    pattern412 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons5, cons259)
    rule412 = ReplacementRule(pattern412, replacement412)

    pattern413 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons258, cons149, cons272)
    rule413 = ReplacementRule(pattern413, replacement413)

    pattern414 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**S(2), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons259, cons149, cons272)
    rule414 = ReplacementRule(pattern414, replacement414)

    pattern415 = Pattern(Integral((x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons3, cons8, cons50, cons127, cons210, cons19, cons4, cons149)
    rule415 = ReplacementRule(pattern415, replacement415)

    pattern416 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons317, cons259, cons149, cons45, cons270)
    rule416 = ReplacementRule(pattern416, replacement416)

    pattern417 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons317, cons228, cons258, cons149)
    rule417 = ReplacementRule(pattern417, replacement417)

    pattern418 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons317, cons259, cons149)
    rule418 = ReplacementRule(pattern418, replacement418)

    pattern419 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons317, cons228, cons281, cons130)
    rule419 = ReplacementRule(pattern419, replacement419)

    pattern420 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons282, cons130)
    rule420 = ReplacementRule(pattern420, replacement420)

    pattern421 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281)
    rule421 = ReplacementRule(pattern421, replacement421)

    pattern422 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282)
    rule422 = ReplacementRule(pattern422, replacement422)

    pattern423 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons281, cons244, cons352)
    rule423 = ReplacementRule(pattern423, replacement423)

    pattern424 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons282, cons244, cons353)
    rule424 = ReplacementRule(pattern424, replacement424)

    pattern425 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons244, cons13, cons139, cons354)
    rule425 = ReplacementRule(pattern425, replacement425)

    pattern426 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons244, cons13, cons139, cons355)
    rule426 = ReplacementRule(pattern426, replacement426)

    pattern427 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons281, cons244)
    rule427 = ReplacementRule(pattern427, replacement427)

    pattern428 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons282, cons244)
    rule428 = ReplacementRule(pattern428, replacement428)

    pattern429 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons228, cons281, cons356)
    rule429 = ReplacementRule(pattern429, replacement429)

    pattern430 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons282, cons357)
    rule430 = ReplacementRule(pattern430, replacement430)

    pattern431 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons13, cons139)
    rule431 = ReplacementRule(pattern431, replacement431)

    pattern432 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons13, cons139)
    rule432 = ReplacementRule(pattern432, replacement432)

    pattern433 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons228, cons281, cons289)
    rule433 = ReplacementRule(pattern433, replacement433)

    pattern434 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons282, cons289)
    rule434 = ReplacementRule(pattern434, replacement434)

    pattern435 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons2, cons8, cons50, cons127, cons210, cons5, cons358, cons359)
    rule435 = ReplacementRule(pattern435, replacement435)

    pattern436 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons360, cons290, cons291)
    rule436 = ReplacementRule(pattern436, replacement436)

    pattern437 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons246, cons165, cons249, cons361, cons362)
    rule437 = ReplacementRule(pattern437, replacement437)

    pattern438 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons246, cons165, cons249, cons361, cons362)
    rule438 = ReplacementRule(pattern438, replacement438)

    pattern439 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons317, cons228, cons281, cons13, cons165, cons363, cons68, cons301, cons364)
    rule439 = ReplacementRule(pattern439, replacement439)

    pattern440 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons317, cons282, cons13, cons165, cons363, cons68, cons301, cons364)
    rule440 = ReplacementRule(pattern440, replacement440)

    pattern441 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons317, cons228, cons281, cons13, cons165, cons365, cons305, cons364)
    rule441 = ReplacementRule(pattern441, replacement441)

    pattern442 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons317, cons282, cons13, cons165, cons365, cons305, cons364)
    rule442 = ReplacementRule(pattern442, replacement442)

    pattern443 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons246, cons139, cons168, cons366)
    rule443 = ReplacementRule(pattern443, replacement443)

    pattern444 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons246, cons139, cons168, cons367)
    rule444 = ReplacementRule(pattern444, replacement444)

    pattern445 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons246, cons139, cons170, cons368, cons364)
    rule445 = ReplacementRule(pattern445, replacement445)

    pattern446 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons246, cons139, cons170, cons368, cons364)
    rule446 = ReplacementRule(pattern446, replacement446)

    pattern447 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons317, cons228, cons281, cons13, cons139, cons364)
    rule447 = ReplacementRule(pattern447, replacement447)

    pattern448 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons13, cons139, cons364)
    rule448 = ReplacementRule(pattern448, replacement448)

    pattern449 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons20)
    rule449 = ReplacementRule(pattern449, replacement449)

    pattern450 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons20)
    rule450 = ReplacementRule(pattern450, replacement450)

    pattern451 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons369, cons170)
    rule451 = ReplacementRule(pattern451, replacement451)

    pattern452 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons369, cons170)
    rule452 = ReplacementRule(pattern452, replacement452)

    pattern453 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281)
    rule453 = ReplacementRule(pattern453, replacement453)

    pattern454 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282)
    rule454 = ReplacementRule(pattern454, replacement454)

    pattern455 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons317, cons228, cons281, cons369, cons96)
    rule455 = ReplacementRule(pattern455, replacement455)

    pattern456 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons317, cons282, cons369, cons96)
    rule456 = ReplacementRule(pattern456, replacement456)

    pattern457 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons358)
    rule457 = ReplacementRule(pattern457, replacement457)

    pattern458 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons358)
    rule458 = ReplacementRule(pattern458, replacement458)

    pattern459 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons228, cons281, cons33, cons170, cons321, cons370, cons364)
    rule459 = ReplacementRule(pattern459, replacement459)

    pattern460 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons282, cons33, cons170, cons321, cons370, cons364)
    rule460 = ReplacementRule(pattern460, replacement460)

    pattern461 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons228, cons281, cons33, cons96, cons364)
    rule461 = ReplacementRule(pattern461, replacement461)

    pattern462 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons5, cons317, cons282, cons33, cons96, cons364)
    rule462 = ReplacementRule(pattern462, replacement462)

    pattern463 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons281, cons371, cons68)
    rule463 = ReplacementRule(pattern463, replacement463)

    pattern464 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons282, cons371, cons68)
    rule464 = ReplacementRule(pattern464, replacement464)

    pattern465 = Pattern(Integral((f_ + x_*WC('g', S(1)))/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons372, cons373, cons374)
    rule465 = ReplacementRule(pattern465, replacement465)

    pattern466 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons127, cons210, cons228)
    rule466 = ReplacementRule(pattern466, replacement466)

    pattern467 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)))), x_), cons2, cons8, cons127, cons210, cons375)
    rule467 = ReplacementRule(pattern467, replacement467)

    pattern468 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(e_*x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons50, cons127, cons210, cons228)
    rule468 = ReplacementRule(pattern468, replacement468)

    pattern469 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(e_*x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)))), x_), cons2, cons8, cons50, cons127, cons210, cons376)
    rule469 = ReplacementRule(pattern469, replacement469)

    pattern470 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons281)
    rule470 = ReplacementRule(pattern470, replacement470)

    pattern471 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons282)
    rule471 = ReplacementRule(pattern471, replacement471)

    pattern472 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons228, cons281, cons377)
    rule472 = ReplacementRule(pattern472, replacement472)

    pattern473 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons8, cons29, cons50, cons127, cons210, cons282, cons377)
    rule473 = ReplacementRule(pattern473, replacement473)

    pattern474 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons151, cons165)
    rule474 = ReplacementRule(pattern474, replacement474)

    pattern475 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons151, cons165)
    rule475 = ReplacementRule(pattern475, replacement475)

    pattern476 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons378, cons369)
    rule476 = ReplacementRule(pattern476, With476)

    pattern477 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons378, cons369)
    rule477 = ReplacementRule(pattern477, With477)

    pattern478 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons126, cons338, cons379)
    rule478 = ReplacementRule(pattern478, replacement478)

    pattern479 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons126, cons338, cons379)
    rule479 = ReplacementRule(pattern479, replacement479)

    pattern480 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons126, cons338)
    rule480 = ReplacementRule(pattern480, replacement480)

    pattern481 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons126, cons338)
    rule481 = ReplacementRule(pattern481, replacement481)

    pattern482 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons228, cons281, cons95)
    rule482 = ReplacementRule(pattern482, replacement482)

    pattern483 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons8, cons29, cons50, cons127, cons210, cons282, cons95)
    rule483 = ReplacementRule(pattern483, replacement483)

    pattern484 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons228, cons281, cons21, cons25, cons95, cons170, cons167)
    rule484 = ReplacementRule(pattern484, replacement484)

    pattern485 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons282, cons21, cons25, cons95, cons170, cons167)
    rule485 = ReplacementRule(pattern485, replacement485)

    pattern486 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons228, cons281, cons21, cons25, cons95, cons170, cons90)
    rule486 = ReplacementRule(pattern486, replacement486)

    pattern487 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons282, cons21, cons25, cons95, cons170, cons90)
    rule487 = ReplacementRule(pattern487, replacement487)

    pattern488 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons228, cons281, cons21, cons25, cons95, cons170, cons91)
    rule488 = ReplacementRule(pattern488, replacement488)

    pattern489 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons282, cons21, cons25, cons95, cons170, cons91)
    rule489 = ReplacementRule(pattern489, replacement489)

    pattern490 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons228, cons281, cons75)
    rule490 = ReplacementRule(pattern490, replacement490)

    pattern491 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons282, cons75)
    rule491 = ReplacementRule(pattern491, replacement491)

    pattern492 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons228, cons281, cons21, cons25)
    rule492 = ReplacementRule(pattern492, replacement492)

    pattern493 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons282, cons21, cons25)
    rule493 = ReplacementRule(pattern493, replacement493)

    pattern494 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons380)
    rule494 = ReplacementRule(pattern494, replacement494)

    pattern495 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons380)
    rule495 = ReplacementRule(pattern495, replacement495)

    pattern496 = Pattern(Integral((x_*WC('g', S(1)))**WC('n', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons19, cons4, cons5, cons360, cons290, cons291)
    rule496 = ReplacementRule(pattern496, replacement496)

    pattern497 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons25, cons149, cons340, cons165, cons91)
    rule497 = ReplacementRule(pattern497, replacement497)

    pattern498 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons25, cons149, cons340, cons165, cons91)
    rule498 = ReplacementRule(pattern498, replacement498)

    pattern499 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons25, cons149, cons340, cons139, cons90)
    rule499 = ReplacementRule(pattern499, replacement499)

    pattern500 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons25, cons149, cons340, cons139, cons90)
    rule500 = ReplacementRule(pattern500, replacement500)

    pattern501 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281)
    rule501 = ReplacementRule(pattern501, With501)

    pattern502 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282)
    rule502 = ReplacementRule(pattern502, With502)

    pattern503 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281, cons82)
    rule503 = ReplacementRule(pattern503, replacement503)

    pattern504 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282, cons82)
    rule504 = ReplacementRule(pattern504, replacement504)

    pattern505 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons317, cons228, cons281)
    rule505 = ReplacementRule(pattern505, replacement505)

    pattern506 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons317, cons282)
    rule506 = ReplacementRule(pattern506, replacement506)

    pattern507 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(f_ + x_*WC('g', S(1)))**S(2), x_), cons2, cons8, cons50, cons127, cons210, cons19, cons5, cons381)
    rule507 = ReplacementRule(pattern507, replacement507)

    pattern508 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(f_ + x_*WC('g', S(1)))**S(3), x_), cons2, cons8, cons50, cons127, cons210, cons19, cons5, cons381)
    rule508 = ReplacementRule(pattern508, replacement508)

    pattern509 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons228, cons281, cons150)
    rule509 = ReplacementRule(pattern509, replacement509)

    pattern510 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons317, cons282, cons150)
    rule510 = ReplacementRule(pattern510, replacement510)

    pattern511 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons382)
    rule511 = ReplacementRule(pattern511, replacement511)

    pattern512 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons383)
    rule512 = ReplacementRule(pattern512, replacement512)

    pattern513 = Pattern(Integral((u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(u_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(a_ + u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons70, cons71)
    rule513 = ReplacementRule(pattern513, replacement513)

    pattern514 = Pattern(Integral((a_ + u_**S(2)*WC('c', S(1)))**WC('p', S(1))*(u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(u_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons70, cons71)
    rule514 = ReplacementRule(pattern514, replacement514)

    pattern515 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons384, cons385, cons386, cons387)
    rule515 = ReplacementRule(pattern515, replacement515)

    pattern516 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons384, cons385, cons149, cons388, cons389)
    rule516 = ReplacementRule(pattern516, replacement516)

    pattern517 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons47, cons149)
    rule517 = ReplacementRule(pattern517, replacement517)

    pattern518 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons127, cons5, cons52, cons47, cons149)
    rule518 = ReplacementRule(pattern518, replacement518)

    pattern519 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons52, cons390, cons391, cons392)
    rule519 = ReplacementRule(pattern519, replacement519)

    pattern520 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons8, cons29, cons50, cons127, cons52, cons393, cons391, cons392)
    rule520 = ReplacementRule(pattern520, replacement520)

    pattern521 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons127, cons52, cons391, cons394)
    rule521 = ReplacementRule(pattern521, replacement521)

    pattern522 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons127, cons52, cons395)
    rule522 = ReplacementRule(pattern522, replacement522)

    pattern523 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons396, cons395)
    rule523 = ReplacementRule(pattern523, replacement523)

    pattern524 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons396, cons395)
    rule524 = ReplacementRule(pattern524, replacement524)

    pattern525 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons396, cons397, cons398, cons399)
    rule525 = ReplacementRule(pattern525, replacement525)

    pattern526 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons8, cons29, cons50, cons127, cons396, cons397, cons398, cons400)
    rule526 = ReplacementRule(pattern526, replacement526)

    pattern527 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons127, cons397, cons398, cons401)
    rule527 = ReplacementRule(pattern527, replacement527)

    pattern528 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons29, cons50, cons127, cons2, cons3, cons8, cons52, cons396, cons402, cons403, cons399)
    rule528 = ReplacementRule(pattern528, replacement528)

    pattern529 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons29, cons50, cons127, cons2, cons8, cons52, cons396, cons402, cons403, cons400)
    rule529 = ReplacementRule(pattern529, replacement529)

    pattern530 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons29, cons127, cons2, cons3, cons8, cons52, cons402, cons403, cons401)
    rule530 = ReplacementRule(pattern530, replacement530)

    pattern531 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons228, cons396, cons404, cons139, cons405)
    rule531 = ReplacementRule(pattern531, replacement531)

    pattern532 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons127, cons228, cons404, cons139, cons405)
    rule532 = ReplacementRule(pattern532, replacement532)

    pattern533 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons396, cons404, cons139, cons405)
    rule533 = ReplacementRule(pattern533, replacement533)

    pattern534 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons52, cons228, cons396, cons13, cons139, cons406, cons407)
    rule534 = ReplacementRule(pattern534, replacement534)

    pattern535 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons127, cons52, cons228, cons13, cons139, cons408, cons407)
    rule535 = ReplacementRule(pattern535, replacement535)

    pattern536 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons52, cons396, cons13, cons139, cons409, cons407)
    rule536 = ReplacementRule(pattern536, replacement536)

    pattern537 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons52, cons228, cons396, cons13, cons148, cons410, cons411)
    rule537 = ReplacementRule(pattern537, replacement537)

    pattern538 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons127, cons52, cons228, cons13, cons148, cons410, cons411)
    rule538 = ReplacementRule(pattern538, replacement538)

    pattern539 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons52, cons396, cons13, cons148, cons410, cons411)
    rule539 = ReplacementRule(pattern539, replacement539)

    pattern540 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons228, cons396, CustomConstraint(With540))
    rule540 = ReplacementRule(pattern540, replacement540)

    pattern541 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons127, cons228, CustomConstraint(With541))
    rule541 = ReplacementRule(pattern541, replacement541)

    pattern542 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons228, cons396, cons412)
    rule542 = ReplacementRule(pattern542, replacement542)

    pattern543 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons228, cons396, cons413, cons233)
    rule543 = ReplacementRule(pattern543, With543)

    pattern544 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons396, cons414)
    rule544 = ReplacementRule(pattern544, replacement544)

    pattern545 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons127, cons228, cons233)
    rule545 = ReplacementRule(pattern545, With545)

    pattern546 = Pattern(Integral(S(1)/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons228, cons396, cons413, cons415)
    rule546 = ReplacementRule(pattern546, With546)

    pattern547 = Pattern(Integral(S(1)/((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons396, cons416)
    rule547 = ReplacementRule(pattern547, With547)

    pattern548 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons127, cons228, cons415)
    rule548 = ReplacementRule(pattern548, With548)

    pattern549 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons228, cons396)
    rule549 = ReplacementRule(pattern549, replacement549)

    pattern550 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(d_ + x_**S(2)*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons127, cons228)
    rule550 = ReplacementRule(pattern550, replacement550)

    pattern551 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)))/(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons396)
    rule551 = ReplacementRule(pattern551, replacement551)

    pattern552 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons228, cons396)
    rule552 = ReplacementRule(pattern552, With552)

    pattern553 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons127, cons228)
    rule553 = ReplacementRule(pattern553, With553)

    pattern554 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons417)
    rule554 = ReplacementRule(pattern554, replacement554)

    pattern555 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons8, cons29, cons50, cons127, cons5, cons52, cons418)
    rule555 = ReplacementRule(pattern555, replacement555)

    pattern556 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons70, cons71)
    rule556 = ReplacementRule(pattern556, replacement556)

    pattern557 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons5, cons52, cons70, cons71)
    rule557 = ReplacementRule(pattern557, replacement557)

    pattern558 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons384, cons385, cons386, cons387)
    rule558 = ReplacementRule(pattern558, replacement558)

    pattern559 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons384, cons385, cons149, cons388, cons389)
    rule559 = ReplacementRule(pattern559, replacement559)

    pattern560 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons52, cons47)
    rule560 = ReplacementRule(pattern560, replacement560)

    pattern561 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons19, cons5, cons52, cons47)
    rule561 = ReplacementRule(pattern561, replacement561)

    pattern562 = Pattern(Integral((g_ + x_*WC('h', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons419, cons420, cons20)
    rule562 = ReplacementRule(pattern562, replacement562)

    pattern563 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons421, cons420, cons20)
    rule563 = ReplacementRule(pattern563, replacement563)

    pattern564 = Pattern(Integral((g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons5, cons419, cons422, cons20)
    rule564 = ReplacementRule(pattern564, replacement564)

    pattern565 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons5, cons421, cons422, cons20)
    rule565 = ReplacementRule(pattern565, replacement565)

    pattern566 = Pattern(Integral(x_**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons50, cons127, cons52, cons228, cons423, cons40)
    rule566 = ReplacementRule(pattern566, replacement566)

    pattern567 = Pattern(Integral(x_**WC('p', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons50, cons127, cons52, cons424, cons40)
    rule567 = ReplacementRule(pattern567, replacement567)

    pattern568 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons425, cons426, cons272)
    rule568 = ReplacementRule(pattern568, replacement568)

    pattern569 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons427, cons428, cons272)
    rule569 = ReplacementRule(pattern569, replacement569)

    pattern570 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons19, cons5, cons429, cons426, cons272)
    rule570 = ReplacementRule(pattern570, replacement570)

    pattern571 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons228, cons396, cons430)
    rule571 = ReplacementRule(pattern571, replacement571)

    pattern572 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons396, cons430)
    rule572 = ReplacementRule(pattern572, replacement572)

    pattern573 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons19, cons228, cons430)
    rule573 = ReplacementRule(pattern573, replacement573)

    pattern574 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons19, cons430)
    rule574 = ReplacementRule(pattern574, replacement574)

    pattern575 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons228, cons396, cons33, cons96, cons431)
    rule575 = ReplacementRule(pattern575, replacement575)

    pattern576 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons396, cons33, cons96, cons432)
    rule576 = ReplacementRule(pattern576, replacement576)

    pattern577 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons5, cons228, cons33, cons96, cons431)
    rule577 = ReplacementRule(pattern577, replacement577)

    pattern578 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons5, cons33, cons96, cons432)
    rule578 = ReplacementRule(pattern578, replacement578)

    pattern579 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, cons431)
    rule579 = ReplacementRule(pattern579, replacement579)

    pattern580 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2)*(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons396, cons432)
    rule580 = ReplacementRule(pattern580, replacement580)

    pattern581 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons228, cons431)
    rule581 = ReplacementRule(pattern581, replacement581)

    pattern582 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2)*(g_ + x_*WC('h', S(1)))), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons432)
    rule582 = ReplacementRule(pattern582, replacement582)

    pattern583 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, cons246, cons139, cons168)
    rule583 = ReplacementRule(pattern583, replacement583)

    pattern584 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons396, cons246, cons139, cons168)
    rule584 = ReplacementRule(pattern584, replacement584)

    pattern585 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons228, cons246, cons139, cons168)
    rule585 = ReplacementRule(pattern585, replacement585)

    pattern586 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(g_ + x_*WC('h', S(1)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons246, cons139, cons168)
    rule586 = ReplacementRule(pattern586, replacement586)

    pattern587 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons228, cons396, cons13, cons139, cons433)
    rule587 = ReplacementRule(pattern587, replacement587)

    pattern588 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons396, cons13, cons139, cons432)
    rule588 = ReplacementRule(pattern588, replacement588)

    pattern589 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons19, cons228, cons13, cons139, cons433)
    rule589 = ReplacementRule(pattern589, replacement589)

    pattern590 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons19, cons13, cons139, cons432)
    rule590 = ReplacementRule(pattern590, replacement590)

    pattern591 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons228, cons396, cons244)
    rule591 = ReplacementRule(pattern591, replacement591)

    pattern592 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons396, cons244)
    rule592 = ReplacementRule(pattern592, replacement592)

    pattern593 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons19, cons5, cons228, cons244)
    rule593 = ReplacementRule(pattern593, replacement593)

    pattern594 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons19, cons5, cons244)
    rule594 = ReplacementRule(pattern594, replacement594)

    pattern595 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons228, cons396, cons272)
    rule595 = ReplacementRule(pattern595, replacement595)

    pattern596 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons396, cons272)
    rule596 = ReplacementRule(pattern596, replacement596)

    pattern597 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons19, cons5, cons228, cons272)
    rule597 = ReplacementRule(pattern597, replacement597)

    pattern598 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)))*(g_ + x_*WC('h', S(1)))**WC('m', S(1)), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons19, cons5, cons272)
    rule598 = ReplacementRule(pattern598, replacement598)

    pattern599 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, cons222, cons165)
    rule599 = ReplacementRule(pattern599, replacement599)

    pattern600 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons396, cons222, cons434)
    rule600 = ReplacementRule(pattern600, replacement600)

    pattern601 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, cons404, cons139, cons405)
    rule601 = ReplacementRule(pattern601, replacement601)

    pattern602 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons396, cons404, cons139, cons405)
    rule602 = ReplacementRule(pattern602, replacement602)

    pattern603 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons228, cons404, cons139, cons405)
    rule603 = ReplacementRule(pattern603, replacement603)

    pattern604 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons52, cons228, cons396, cons13, cons139, cons406, cons407)
    rule604 = ReplacementRule(pattern604, replacement604)

    pattern605 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons52, cons396, cons13, cons139, cons409, cons407)
    rule605 = ReplacementRule(pattern605, replacement605)

    pattern606 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons52, cons228, cons13, cons139, cons408, cons407)
    rule606 = ReplacementRule(pattern606, replacement606)

    pattern607 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons52, cons228, cons396, cons13, cons165, cons435)
    rule607 = ReplacementRule(pattern607, replacement607)

    pattern608 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons52, cons396, cons13, cons165, cons435)
    rule608 = ReplacementRule(pattern608, replacement608)

    pattern609 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons52, cons228, cons13, cons165, cons435)
    rule609 = ReplacementRule(pattern609, replacement609)

    pattern610 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, CustomConstraint(With610))
    rule610 = ReplacementRule(pattern610, replacement610)

    pattern611 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons228, CustomConstraint(With611))
    rule611 = ReplacementRule(pattern611, replacement611)

    pattern612 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons436)
    rule612 = ReplacementRule(pattern612, replacement612)

    pattern613 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons437)
    rule613 = ReplacementRule(pattern613, With613)

    pattern614 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, cons412, cons438)
    rule614 = ReplacementRule(pattern614, replacement614)

    pattern615 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, cons412, cons439)
    rule615 = ReplacementRule(pattern615, replacement615)

    pattern616 = Pattern(Integral(x_/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons228, cons396, cons385)
    rule616 = ReplacementRule(pattern616, replacement616)

    pattern617 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, cons385, cons440)
    rule617 = ReplacementRule(pattern617, replacement617)

    pattern618 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, cons385, cons441)
    rule618 = ReplacementRule(pattern618, replacement618)

    pattern619 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, cons374, cons442)
    rule619 = ReplacementRule(pattern619, replacement619)

    pattern620 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons443)
    rule620 = ReplacementRule(pattern620, replacement620)

    pattern621 = Pattern(Integral((g_ + x_*WC('h', S(1)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons228, cons444)
    rule621 = ReplacementRule(pattern621, replacement621)

    pattern622 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, cons233)
    rule622 = ReplacementRule(pattern622, With622)

    pattern623 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons396, cons414)
    rule623 = ReplacementRule(pattern623, With623)

    pattern624 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons228, cons233)
    rule624 = ReplacementRule(pattern624, With624)

    pattern625 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396, cons374, cons415)
    rule625 = ReplacementRule(pattern625, With625)

    pattern626 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons396, cons416)
    rule626 = ReplacementRule(pattern626, With626)

    pattern627 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons228, cons415)
    rule627 = ReplacementRule(pattern627, With627)

    pattern628 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons396)
    rule628 = ReplacementRule(pattern628, With628)

    pattern629 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons211, cons228)
    rule629 = ReplacementRule(pattern629, With629)

    pattern630 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons412, cons445, cons446, cons447)
    rule630 = ReplacementRule(pattern630, With630)

    pattern631 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons448, cons449, cons45)
    rule631 = ReplacementRule(pattern631, replacement631)

    pattern632 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons412, cons445, cons446, cons315)
    rule632 = ReplacementRule(pattern632, With632)

    pattern633 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons8, cons29, cons127, cons210, cons211, cons448, cons449, cons450)
    rule633 = ReplacementRule(pattern633, replacement633)

    pattern634 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons451)
    rule634 = ReplacementRule(pattern634, replacement634)

    pattern635 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons452)
    rule635 = ReplacementRule(pattern635, replacement635)

    pattern636 = Pattern(Integral((u_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons52, cons70, cons71)
    rule636 = ReplacementRule(pattern636, replacement636)

    pattern637 = Pattern(Integral((u_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons52, cons70, cons71)
    rule637 = ReplacementRule(pattern637, replacement637)

    pattern638 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*z_**WC('m', S(1)), x_), cons19, cons5, cons52, cons453, cons454, cons455)
    rule638 = ReplacementRule(pattern638, replacement638)

    pattern639 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons19, cons4, cons5, cons52, cons338, cons126, cons379)
    rule639 = ReplacementRule(pattern639, replacement639)

    pattern640 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons19, cons4, cons5, cons52, cons130, cons86)
    rule640 = ReplacementRule(pattern640, replacement640)

    pattern641 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons19, cons4, cons5, cons52, cons338, cons126)
    rule641 = ReplacementRule(pattern641, replacement641)

    pattern642 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons5, cons52, cons384, cons385, cons386, cons387)
    rule642 = ReplacementRule(pattern642, replacement642)

    pattern643 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons5, cons52, cons384, cons385, cons386, cons387)
    rule643 = ReplacementRule(pattern643, replacement643)

    pattern644 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons5, cons52, cons384, cons385, cons149, cons388, cons389)
    rule644 = ReplacementRule(pattern644, replacement644)

    pattern645 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons5, cons52, cons384, cons385, cons149, cons388, cons389)
    rule645 = ReplacementRule(pattern645, replacement645)

    pattern646 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons5, cons52, cons47)
    rule646 = ReplacementRule(pattern646, replacement646)

    pattern647 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons5, cons52, cons47)
    rule647 = ReplacementRule(pattern647, replacement647)

    pattern648 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons37, cons38, cons5, cons52, cons47)
    rule648 = ReplacementRule(pattern648, replacement648)

    pattern649 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons38, cons5, cons52, cons47)
    rule649 = ReplacementRule(pattern649, replacement649)

    pattern650 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons228, cons396, cons222, cons165)
    rule650 = ReplacementRule(pattern650, replacement650)

    pattern651 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons228, cons396, cons222, cons165)
    rule651 = ReplacementRule(pattern651, replacement651)

    pattern652 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons396, cons222, cons434)
    rule652 = ReplacementRule(pattern652, replacement652)

    pattern653 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons38, cons396, cons222, cons434)
    rule653 = ReplacementRule(pattern653, replacement653)

    pattern654 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons228, cons396, cons404, cons139, cons405)
    rule654 = ReplacementRule(pattern654, replacement654)

    pattern655 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons228, cons396, cons404, cons139, cons405)
    rule655 = ReplacementRule(pattern655, replacement655)

    pattern656 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons396, cons404, cons139, cons405)
    rule656 = ReplacementRule(pattern656, replacement656)

    pattern657 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons38, cons396, cons404, cons139, cons405)
    rule657 = ReplacementRule(pattern657, replacement657)

    pattern658 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons37, cons38, cons228, cons404, cons139, cons405)
    rule658 = ReplacementRule(pattern658, replacement658)

    pattern659 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons38, cons228, cons404, cons139, cons405)
    rule659 = ReplacementRule(pattern659, replacement659)

    pattern660 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons52, cons228, cons396, cons13, cons139, cons406, cons407)
    rule660 = ReplacementRule(pattern660, replacement660)

    pattern661 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons52, cons228, cons396, cons13, cons139, cons406, cons407)
    rule661 = ReplacementRule(pattern661, replacement661)

    pattern662 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons52, cons396, cons13, cons139, cons409, cons407)
    rule662 = ReplacementRule(pattern662, replacement662)

    pattern663 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons38, cons52, cons396, cons13, cons139, cons409, cons407)
    rule663 = ReplacementRule(pattern663, replacement663)

    pattern664 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons37, cons38, cons52, cons228, cons13, cons139, cons408, cons407)
    rule664 = ReplacementRule(pattern664, replacement664)

    pattern665 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons38, cons52, cons228, cons13, cons139, cons408, cons407)
    rule665 = ReplacementRule(pattern665, replacement665)

    pattern666 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons52, cons228, cons396, cons13, cons165, cons435, cons456)
    rule666 = ReplacementRule(pattern666, replacement666)

    pattern667 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons52, cons228, cons396, cons13, cons165, cons435, cons456)
    rule667 = ReplacementRule(pattern667, replacement667)

    pattern668 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons52, cons396, cons13, cons165, cons435, cons456)
    rule668 = ReplacementRule(pattern668, replacement668)

    pattern669 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons38, cons52, cons396, cons13, cons165, cons435, cons456)
    rule669 = ReplacementRule(pattern669, replacement669)

    pattern670 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons37, cons38, cons52, cons228, cons13, cons165, cons435, cons456)
    rule670 = ReplacementRule(pattern670, replacement670)

    pattern671 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons38, cons52, cons228, cons13, cons165, cons435, cons456)
    rule671 = ReplacementRule(pattern671, replacement671)

    pattern672 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons228, cons396, CustomConstraint(With672))
    rule672 = ReplacementRule(pattern672, replacement672)

    pattern673 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons228, cons396, CustomConstraint(With673))
    rule673 = ReplacementRule(pattern673, replacement673)

    pattern674 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons37, cons38, cons228, CustomConstraint(With674))
    rule674 = ReplacementRule(pattern674, replacement674)

    pattern675 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons38, cons228, CustomConstraint(With675))
    rule675 = ReplacementRule(pattern675, replacement675)

    pattern676 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons228, cons396)
    rule676 = ReplacementRule(pattern676, replacement676)

    pattern677 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons228, cons396)
    rule677 = ReplacementRule(pattern677, replacement677)

    pattern678 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons396)
    rule678 = ReplacementRule(pattern678, replacement678)

    pattern679 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons38, cons396)
    rule679 = ReplacementRule(pattern679, replacement679)

    pattern680 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons37, cons38, cons228)
    rule680 = ReplacementRule(pattern680, replacement680)

    pattern681 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons38, cons228)
    rule681 = ReplacementRule(pattern681, replacement681)

    pattern682 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons5, cons52, cons457)
    rule682 = ReplacementRule(pattern682, replacement682)

    pattern683 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons5, cons52, cons458)
    rule683 = ReplacementRule(pattern683, replacement683)

    pattern684 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons5, cons52, cons459)
    rule684 = ReplacementRule(pattern684, replacement684)

    pattern685 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons8, cons29, cons50, cons127, cons36, cons38, cons5, cons52, cons460)
    rule685 = ReplacementRule(pattern685, replacement685)

    pattern686 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons5, cons52, cons70, cons71)
    rule686 = ReplacementRule(pattern686, replacement686)

    pattern687 = Pattern(Integral((u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons5, cons52, cons70, cons71)
    rule687 = ReplacementRule(pattern687, replacement687)

    pattern688 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons5, cons52, cons70, cons71)
    rule688 = ReplacementRule(pattern688, replacement688)

    pattern689 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('C', S(1)) + u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons5, cons52, cons70, cons71)
    rule689 = ReplacementRule(pattern689, replacement689)

    pattern690 = Pattern(Integral((u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons5, cons52, cons70, cons71)
    rule690 = ReplacementRule(pattern690, replacement690)

    pattern691 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons38, cons5, cons52, cons70, cons71)
    rule691 = ReplacementRule(pattern691, replacement691)
    return [rule192, rule193, rule194, rule195, rule196, rule197, rule198, rule199, rule200, rule201, rule202, rule203, rule204, rule205, rule206, rule207, rule208, rule209, rule210, rule211, rule212, rule213, rule214, rule215, rule216, rule217, rule218, rule219, rule220, rule221, rule222, rule223, rule224, rule225, rule226, rule227, rule228, rule229, rule230, rule231, rule232, rule233, rule234, rule235, rule236, rule237, rule238, rule239, rule240, rule241, rule242, rule243, rule244, rule245, rule246, rule247, rule248, rule249, rule250, rule251, rule252, rule253, rule254, rule255, rule256, rule257, rule258, rule259, rule260, rule261, rule262, rule263, rule264, rule265, rule266, rule267, rule268, rule269, rule270, rule271, rule272, rule273, rule274, rule275, rule276, rule277, rule278, rule279, rule280, rule281, rule282, rule283, rule284, rule285, rule286, rule287, rule288, rule289, rule290, rule291, rule292, rule293, rule294, rule295, rule296, rule297, rule298, rule299, rule300, rule301, rule302, rule303, rule304, rule305, rule306, rule307, rule308, rule309, rule310, rule311, rule312, rule313, rule314, rule315, rule316, rule317, rule318, rule319, rule320, rule321, rule322, rule323, rule324, rule325, rule326, rule327, rule328, rule329, rule330, rule331, rule332, rule333, rule334, rule335, rule336, rule337, rule338, rule339, rule340, rule341, rule342, rule343, rule344, rule345, rule346, rule347, rule348, rule349, rule350, rule351, rule352, rule353, rule354, rule355, rule356, rule357, rule358, rule359, rule360, rule361, rule362, rule363, rule364, rule365, rule366, rule367, rule368, rule369, rule370, rule371, rule372, rule373, rule374, rule375, rule376, rule377, rule378, rule379, rule380, rule381, rule382, rule383, rule384, rule385, rule386, rule387, rule388, rule389, rule390, rule391, rule392, rule393, rule394, rule395, rule396, rule397, rule398, rule399, rule400, rule401, rule402, rule403, rule404, rule405, rule406, rule407, rule408, rule409, rule410, rule411, rule412, rule413, rule414, rule415, rule416, rule417, rule418, rule419, rule420, rule421, rule422, rule423, rule424, rule425, rule426, rule427, rule428, rule429, rule430, rule431, rule432, rule433, rule434, rule435, rule436, rule437, rule438, rule439, rule440, rule441, rule442, rule443, rule444, rule445, rule446, rule447, rule448, rule449, rule450, rule451, rule452, rule453, rule454, rule455, rule456, rule457, rule458, rule459, rule460, rule461, rule462, rule463, rule464, rule465, rule466, rule467, rule468, rule469, rule470, rule471, rule472, rule473, rule474, rule475, rule476, rule477, rule478, rule479, rule480, rule481, rule482, rule483, rule484, rule485, rule486, rule487, rule488, rule489, rule490, rule491, rule492, rule493, rule494, rule495, rule496, rule497, rule498, rule499, rule500, rule501, rule502, rule503, rule504, rule505, rule506, rule507, rule508, rule509, rule510, rule511, rule512, rule513, rule514, rule515, rule516, rule517, rule518, rule519, rule520, rule521, rule522, rule523, rule524, rule525, rule526, rule527, rule528, rule529, rule530, rule531, rule532, rule533, rule534, rule535, rule536, rule537, rule538, rule539, rule540, rule541, rule542, rule543, rule544, rule545, rule546, rule547, rule548, rule549, rule550, rule551, rule552, rule553, rule554, rule555, rule556, rule557, rule558, rule559, rule560, rule561, rule562, rule563, rule564, rule565, rule566, rule567, rule568, rule569, rule570, rule571, rule572, rule573, rule574, rule575, rule576, rule577, rule578, rule579, rule580, rule581, rule582, rule583, rule584, rule585, rule586, rule587, rule588, rule589, rule590, rule591, rule592, rule593, rule594, rule595, rule596, rule597, rule598, rule599, rule600, rule601, rule602, rule603, rule604, rule605, rule606, rule607, rule608, rule609, rule610, rule611, rule612, rule613, rule614, rule615, rule616, rule617, rule618, rule619, rule620, rule621, rule622, rule623, rule624, rule625, rule626, rule627, rule628, rule629, rule630, rule631, rule632, rule633, rule634, rule635, rule636, rule637, rule638, rule639, rule640, rule641, rule642, rule643, rule644, rule645, rule646, rule647, rule648, rule649, rule650, rule651, rule652, rule653, rule654, rule655, rule656, rule657, rule658, rule659, rule660, rule661, rule662, rule663, rule664, rule665, rule666, rule667, rule668, rule669, rule670, rule671, rule672, rule673, rule674, rule675, rule676, rule677, rule678, rule679, rule680, rule681, rule682, rule683, rule684, rule685, rule686, rule687, rule688, rule689, rule690, rule691, ]





def replacement192(a, b, c, x):
    return Dist((b/S(2) + c*x)/sqrt(a + b*x + c*x**S(2)), Int(S(1)/(b/S(2) + c*x), x), x)


def replacement193(a, b, c, p, x):
    return Simp((b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))), x)


def With194(a, b, c, p, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(c**(-p), Int(Simp(b/S(2) + c*x - q/S(2), x)**p*Simp(b/S(2) + c*x + q/S(2), x)**p, x), x)


def replacement195(a, b, c, p, x):
    return Int(ExpandIntegrand((a + b*x + c*x**S(2))**p, x), x)


def replacement196(a, b, c, p, x):
    return -Dist(p*(-S(4)*a*c + b**S(2))/(S(2)*c*(S(2)*p + S(1))), Int((a + b*x + c*x**S(2))**(p + S(-1)), x), x) + Simp((b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))), x)


def replacement197(a, b, c, x):
    return Simp(-S(2)*(b + S(2)*c*x)/((-S(4)*a*c + b**S(2))*sqrt(a + b*x + c*x**S(2))), x)


def replacement198(a, b, c, p, x):
    return -Dist(S(2)*c*(S(2)*p + S(3))/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp((b + S(2)*c*x)*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def With199(a, b, c, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(c/q, Int(S(1)/Simp(b/S(2) + c*x - q/S(2), x), x), x) - Dist(c/q, Int(S(1)/Simp(b/S(2) + c*x + q/S(2), x), x), x)


def With200(a, b, c, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = -S(4)*a*c/b**S(2) + S(1)
    if And(RationalQ(q), Or(EqQ(q**S(2), S(1)), Not(RationalQ(-S(4)*a*c + b**S(2))))):
        return True
    return False


def replacement200(a, b, c, x):

    q = -S(4)*a*c/b**S(2) + S(1)
    return Dist(-S(2)/b, Subst(Int(S(1)/(q - x**S(2)), x), x, S(1) + S(2)*c*x/b), x)


def replacement201(a, b, c, x):
    return Dist(S(-2), Subst(Int(S(1)/Simp(-S(4)*a*c + b**S(2) - x**S(2), x), x), x, b + S(2)*c*x), x)


def replacement202(a, b, c, p, x):
    return Dist((-S(4)*c/(-S(4)*a*c + b**S(2)))**(-p)/(S(2)*c), Subst(Int(Simp(-x**S(2)/(-S(4)*a*c + b**S(2)) + S(1), x)**p, x), x, b + S(2)*c*x), x)


def replacement203(b, c, x):
    return Dist(S(2), Subst(Int(S(1)/(-c*x**S(2) + S(1)), x), x, x/sqrt(b*x + c*x**S(2))), x)


def replacement204(a, b, c, x):
    return Dist(S(2), Subst(Int(S(1)/(S(4)*c - x**S(2)), x), x, (b + S(2)*c*x)/sqrt(a + b*x + c*x**S(2))), x)


def replacement205(b, c, p, x):
    return Dist((-c*(b*x + c*x**S(2))/b**S(2))**(-p)*(b*x + c*x**S(2))**p, Int((-c*x/b - c**S(2)*x**S(2)/b**S(2))**p, x), x)


def With206(a, b, c, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    d = Denominator(p)
    if LessEqual(S(3), d, S(4)):
        return True
    return False


def replacement206(a, b, c, p, x):

    d = Denominator(p)
    return Dist(d*sqrt((b + S(2)*c*x)**S(2))/(b + S(2)*c*x), Subst(Int(x**(d*(p + S(1)) + S(-1))/sqrt(-S(4)*a*c + b**S(2) + S(4)*c*x**d), x), x, (a + b*x + c*x**S(2))**(S(1)/d)), x)


def With207(a, b, c, p, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return -Simp(((-b - S(2)*c*x + q)/(S(2)*q))**(-p + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Hypergeometric2F1(-p, p + S(1), p + S(2), (b + S(2)*c*x + q)/(S(2)*q))/(q*(p + S(1))), x)


def replacement208(a, b, c, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*x + c*x**S(2))**p, x), x, u), x)


def replacement209(a, b, c, d, e, m, p, x):
    return Simp(c**(-m/S(2) + S(-1)/2)*e**m*(a + b*x + c*x**S(2))**(m/S(2) + p + S(1)/2)/(m + S(2)*p + S(1)), x)


def replacement210(a, b, c, d, e, m, p, x):
    return Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*log(RemoveContent(d + e*x, x))/e, x)


def replacement211(a, b, c, d, e, m, p, x):
    return Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))), x)


def replacement212(a, b, c, d, e, m, p, x):
    return -Simp((b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(1))*(-b*e + S(2)*c*d)), x)


def replacement213(a, b, c, d, e, x):
    return Dist(sqrt(a + b*x + c*x**S(2))/(b + S(2)*c*x), Int((b + S(2)*c*x)/(d + e*x)**S(2), x), x)


def replacement214(a, b, c, d, e, m, x):
    return -Dist((-b*e + S(2)*c*d)*sqrt(a + b*x + c*x**S(2))/(e*(b + S(2)*c*x)*(m + S(2))), Int((d + e*x)**m, x), x) + Simp((d + e*x)**(m + S(1))*sqrt(a + b*x + c*x**S(2))/(e*(m + S(2))), x)


def replacement215(a, b, c, d, e, x):
    return Dist(S(2)*c/(-b*e + S(2)*c*d), Int(S(1)/((d + e*x)*sqrt(a + b*x + c*x**S(2))), x), x) + Simp(-S(4)*c*e*sqrt(a + b*x + c*x**S(2))/((d + e*x)*(-b*e + S(2)*c*d)**S(2)), x)


def replacement216(a, b, c, d, e, m, p, x):
    return -Simp((b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(2))*(-b*e + S(2)*c*d)), x) + Simp(-S(2)*c*e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)**S(2)*(m*p + S(-1))), x)


def replacement217(a, b, c, d, e, p, x):
    return Dist((-b*e + S(2)*c*d)/(S(2)*c), Int((a + b*x + c*x**S(2))**p, x), x) + Simp(e*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))), x)


def replacement218(a, b, c, d, e, m, p, x):
    return Dist(p*(S(2)*p + S(-1))*(-b*e + S(2)*c*d)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(1))), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))), x) - Simp(p*(b + S(2)*c*x)*(d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(1))), x)


def replacement219(a, b, c, d, e, m, p, x):
    return Dist(S(2)*c*p*(S(2)*p + S(-1))/(e**S(2)*(m + S(1))*(m + S(2))), Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))), x) - Simp(p*(b + S(2)*c*x)*(d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))/(e**S(2)*(m + S(1))*(m + S(2))), x)


def replacement220(a, b, c, d, e, m, p, x):
    return Dist(p*(S(2)*p + S(-1))*(-b*e + S(2)*c*d)**S(2)/(S(2)*c*e**S(2)*(m + S(2)*p)*(m + S(2)*p + S(1))), Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1)), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))), x) - Simp(p*(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(-b*e + S(2)*c*d)*(a + b*x + c*x**S(2))**(p + S(-1))/(S(2)*c*e**S(2)*(m + S(2)*p)*(m + S(2)*p + S(1))), x)


def replacement221(a, b, c, d, e, m, p, x):
    return Dist(e**S(2)*m*(m + S(2)*p + S(2))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp((b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((S(2)*p + S(1))*(-b*e + S(2)*c*d)), x) - Simp(e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(m + S(2)*p + S(2))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)), x)


def replacement222(a, b, c, d, e, m, p, x):
    return Dist(e**S(2)*m*(m + S(-1))/(S(2)*c*(p + S(1))*(S(2)*p + S(1))), Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp((b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))), x) - Simp(e*m*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))*(S(2)*p + S(1))), x)


def replacement223(a, b, c, d, e, m, p, x):
    return Dist(S(2)*c*e**S(2)*(m + S(2)*p + S(2))*(m + S(2)*p + S(3))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)**S(2)), Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp((b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((S(2)*p + S(1))*(-b*e + S(2)*c*d)), x) + Simp(-S(2)*c*e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(m + S(2)*p + S(2))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)**S(2)), x)


def replacement224(a, b, c, d, e, m, p, x):
    return Dist(m*(-b*e + S(2)*c*d)/(S(2)*c*(m + S(2)*p + S(1))), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x), x) + Simp((b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p/(S(2)*c*(m + S(2)*p + S(1))), x)


def replacement225(a, b, c, d, e, m, p, x):
    return Dist(S(2)*c*(m + S(2)*p + S(2))/((m + S(1))*(-b*e + S(2)*c*d)), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x), x) - Simp((b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(1))*(-b*e + S(2)*c*d)), x)


def replacement226(a, b, c, d, e, m, p, x):
    return Dist(c**(-IntPart(p))*(b/S(2) + c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((b/S(2) + c*x)**(S(2)*p)*(d + e*x)**m, x), x)


def replacement227(a, b, c, d, e, m, p, x):
    return Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x)


def replacement228(a, c, d, e, m, p, x):
    return Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x)


def replacement229(a, b, c, d, e, m, p, x):
    return Simp(e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))), x)


def replacement230(a, c, d, e, m, p, x):
    return Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(p + S(1))), x)


def replacement231(a, b, c, d, e, m, p, x):
    return Simp(e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-b*e + S(2)*c*d)), x)


def replacement232(a, c, d, e, m, p, x):
    return Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(p + S(1))), x)


def replacement233(a, b, c, d, e, p, x):
    return -Dist(e**S(2)*(p + S(2))/(c*(p + S(1))), Int((a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp(e*(d + e*x)*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))), x)


def replacement234(a, c, d, e, p, x):
    return -Dist(e**S(2)*(p + S(2))/(c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1)), x), x) + Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)/(c*(p + S(1))), x)


def replacement235(a, b, c, d, e, m, p, x):
    return Int((a/d + c*x/e)**(-m)*(a + b*x + c*x**S(2))**(m + p), x)


def replacement236(a, c, d, e, m, p, x):
    return Dist(a**(-m)*d**(S(2)*m), Int((a + c*x**S(2))**(m + p)*(d - e*x)**(-m), x), x)


def replacement237(a, b, c, d, e, m, p, x):
    return Dist((m + p)*(-b*e + S(2)*c*d)/(c*(m + S(2)*p + S(1))), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x), x) + Simp(e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))), x)


def replacement238(a, c, d, e, m, p, x):
    return Dist(S(2)*d*(m + p)/(m + S(2)*p + S(1)), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1)), x), x) + Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))), x)


def replacement239(a, b, c, d, e, m, p, x):
    return Dist(c*(m + S(2)*p + S(2))/((-b*e + S(2)*c*d)*(m + p + S(1))), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x), x) - Simp(e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))), x)


def replacement240(a, c, d, e, m, p, x):
    return Dist((m + S(2)*p + S(2))/(S(2)*d*(m + p + S(1))), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x), x) - Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(m + p + S(1))), x)


def replacement241(a, b, c, d, e, x):
    return Dist(S(2)*e, Subst(Int(S(1)/(-b*e + S(2)*c*d + e**S(2)*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))/sqrt(d + e*x)), x)


def replacement242(a, c, d, e, x):
    return Dist(S(2)*e, Subst(Int(S(1)/(S(2)*c*d + e**S(2)*x**S(2)), x), x, sqrt(a + c*x**S(2))/sqrt(d + e*x)), x)


def replacement243(a, b, c, d, e, m, p, x):
    return -Dist(c*p/(e**S(2)*(m + p + S(1))), Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + p + S(1))), x)


def replacement244(a, c, d, e, m, p, x):
    return -Dist(c*p/(e**S(2)*(m + p + S(1))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2)), x), x) + Simp((a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + p + S(1))), x)


def replacement245(a, b, c, d, e, m, p, x):
    return -Dist(p*(-b*e + S(2)*c*d)/(e**S(2)*(m + S(2)*p + S(1))), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))), x)


def replacement246(a, c, d, e, m, p, x):
    return -Dist(S(2)*c*d*p/(e**S(2)*(m + S(2)*p + S(1))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1)), x), x) + Simp((a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(2)*p + S(1))), x)


def replacement247(a, b, c, d, e, m, p, x):
    return -Dist((-b*e + S(2)*c*d)*(m + S(2)*p + S(2))/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp((d + e*x)**m*(-b*e + S(2)*c*d)*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement248(a, c, d, e, m, p, x):
    return Dist(d*(m + S(2)*p + S(2))/(S(2)*a*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x), x) - Simp(d*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*e*(p + S(1))), x)


def replacement249(a, b, c, d, e, m, p, x):
    return -Dist(e**S(2)*(m + p)/(c*(p + S(1))), Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp(e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))), x)


def replacement250(a, c, d, e, m, p, x):
    return -Dist(e**S(2)*(m + p)/(c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2)), x), x) + Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(p + S(1))), x)


def replacement251(a, b, c, d, e, m, p, x):
    return Dist((m + p)*(-b*e + S(2)*c*d)/(c*(m + S(2)*p + S(1))), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x), x) + Simp(e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))), x)


def replacement252(a, c, d, e, m, p, x):
    return Dist(S(2)*d*(m + p)/(m + S(2)*p + S(1)), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1)), x), x) + Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))), x)


def replacement253(a, b, c, d, e, m, p, x):
    return Dist(c*(m + S(2)*p + S(2))/((-b*e + S(2)*c*d)*(m + p + S(1))), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x), x) - Simp(e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))), x)


def replacement254(a, c, d, e, m, p, x):
    return Dist((m + S(2)*p + S(2))/(S(2)*d*(m + p + S(1))), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x), x) - Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(m + p + S(1))), x)


def replacement255(b, c, e, m, p, x):
    return Dist(x**(-m - p)*(e*x)**m*(b + c*x)**(-p)*(b*x + c*x**S(2))**p, Int(x**(m + p)*(b + c*x)**p, x), x)


def replacement256(a, c, d, e, m, p, x):
    return Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x)


def replacement257(a, b, c, d, e, m, p, x):
    return Dist((d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x), x)


def replacement258(a, c, d, e, m, p, x):
    return Dist((a + c*x**S(2))**FracPart(p)*(d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p)), Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x), x)


def replacement259(a, b, c, d, e, x):
    return Dist(b**S(2)/(d**S(2)*(-S(4)*a*c + b**S(2))), Int((d + e*x)/(a + b*x + c*x**S(2)), x), x) + Dist(-S(4)*b*c/(d*(-S(4)*a*c + b**S(2))), Int(S(1)/(b + S(2)*c*x), x), x)


def replacement260(a, b, c, d, e, m, p, x):
    return Simp(S(2)*c*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement261(a, b, c, d, e, m, p, x):
    return Int(ExpandIntegrand((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x)


def replacement262(a, b, c, d, e, m, p, x):
    return -Dist(b*p/(d*e*(m + S(1))), Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))), x)


def replacement263(a, b, c, d, e, m, p, x):
    return -Dist(d*p*(-S(4)*a*c + b**S(2))/(b*e*(m + S(2)*p + S(1))), Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1)), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))), x)


def replacement264(a, b, c, d, e, m, p, x):
    return -Dist(d*e*(m + S(-1))/(b*(p + S(1))), Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp(d*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))), x)


def replacement265(a, b, c, d, e, m, p, x):
    return -Dist(S(2)*c*(m + S(2)*p + S(3))/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp(S(2)*c*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement266(a, b, c, d, e, x):
    return Dist(S(4)*c, Subst(Int(S(1)/(-S(4)*a*c*e + b**S(2)*e + S(4)*c*e*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))), x)


def replacement267(a, b, c, d, e, x):
    return Dist(S(4)*sqrt(-c/(-S(4)*a*c + b**S(2)))/e, Subst(Int(S(1)/sqrt(Simp(-b**S(2)*x**S(4)/(d**S(2)*(-S(4)*a*c + b**S(2))) + S(1), x)), x), x, sqrt(d + e*x)), x)


def replacement268(a, b, c, d, e, x):
    return Dist(S(4)*sqrt(-c/(-S(4)*a*c + b**S(2)))/e, Subst(Int(x**S(2)/sqrt(Simp(-b**S(2)*x**S(4)/(d**S(2)*(-S(4)*a*c + b**S(2))) + S(1), x)), x), x, sqrt(d + e*x)), x)


def replacement269(a, b, c, d, e, m, x):
    return Dist(sqrt(-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))/sqrt(a + b*x + c*x**S(2)), Int((d + e*x)**m/sqrt(-a*c/(-S(4)*a*c + b**S(2)) - b*c*x/(-S(4)*a*c + b**S(2)) - c**S(2)*x**S(2)/(-S(4)*a*c + b**S(2))), x), x)


def replacement270(a, b, c, d, e, m, p, x):
    return Dist(d**S(2)*(m + S(-1))*(-S(4)*a*c + b**S(2))/(b**S(2)*(m + S(2)*p + S(1))), Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**p, x), x) + Simp(S(2)*d*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(m + S(2)*p + S(1))), x)


def replacement271(a, b, c, d, e, m, p, x):
    return Dist(b**S(2)*(m + S(2)*p + S(3))/(d**S(2)*(m + S(1))*(-S(4)*a*c + b**S(2))), Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x), x) + Simp(-S(2)*b*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(d*(m + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement272(a, b, c, d, e, m, p, x):
    return Dist(S(1)/e, Subst(Int(x**m*(a - b**S(2)/(S(4)*c) + c*x**S(2)/e**S(2))**p, x), x, d + e*x), x)


def replacement273(a, b, c, d, e, m, p, x):
    return Int(ExpandIntegrand((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x)


def replacement274(a, c, d, e, m, p, x):
    return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m, x), x)


def With275(a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist((c*d - e*(b/S(2) - q/S(2)))/q, Int(S(1)/(b/S(2) + c*x - q/S(2)), x), x) - Dist((c*d - e*(b/S(2) + q/S(2)))/q, Int(S(1)/(b/S(2) + c*x + q/S(2)), x), x)


def With276(a, c, d, e, x):
    q = Rt(-a*c, S(2))
    return Dist(-c*d/(S(2)*q) + e/S(2), Int(S(1)/(c*x + q), x), x) + Dist(c*d/(S(2)*q) + e/S(2), Int(S(1)/(c*x - q), x), x)


def replacement277(a, b, c, d, e, x):
    return Dist(e/(S(2)*c), Int((b + S(2)*c*x)/(a + b*x + c*x**S(2)), x), x) + Dist((-b*e + S(2)*c*d)/(S(2)*c), Int(S(1)/(a + b*x + c*x**S(2)), x), x)


def replacement278(a, c, d, e, x):
    return Dist(d, Int(S(1)/(a + c*x**S(2)), x), x) + Dist(e, Int(x/(a + c*x**S(2)), x), x)


def replacement279(a, b, c, d, e, x):
    return Dist(S(2)*e, Subst(Int(x**S(2)/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)), x)


def replacement280(a, c, d, e, x):
    return Dist(S(2)*e, Subst(Int(x**S(2)/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)), x)


def replacement281(a, b, c, d, e, m, x):
    return Int(PolynomialDivide((d + e*x)**m, a + b*x + c*x**S(2), x), x)


def replacement282(a, c, d, e, m, x):
    return Int(PolynomialDivide((d + e*x)**m, a + c*x**S(2), x), x)


def replacement283(a, b, c, d, e, m, x):
    return Dist(S(1)/c, Int((d + e*x)**(m + S(-2))*Simp(-a*e**S(2) + c*d**S(2) + e*x*(-b*e + S(2)*c*d), x)/(a + b*x + c*x**S(2)), x), x) + Simp(e*(d + e*x)**(m + S(-1))/(c*(m + S(-1))), x)


def replacement284(a, c, d, e, m, x):
    return Dist(S(1)/c, Int((d + e*x)**(m + S(-2))*Simp(-a*e**S(2) + c*d**S(2) + S(2)*c*d*e*x, x)/(a + c*x**S(2)), x), x) + Simp(e*(d + e*x)**(m + S(-1))/(c*(m + S(-1))), x)


def replacement285(a, b, c, d, e, x):
    return Dist(e**S(2)/(a*e**S(2) - b*d*e + c*d**S(2)), Int(S(1)/(d + e*x), x), x) + Dist(S(1)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((-b*e + c*d - c*e*x)/(a + b*x + c*x**S(2)), x), x)


def replacement286(a, c, d, e, x):
    return Dist(e**S(2)/(a*e**S(2) + c*d**S(2)), Int(S(1)/(d + e*x), x), x) + Dist(S(1)/(a*e**S(2) + c*d**S(2)), Int((c*d - c*e*x)/(a + c*x**S(2)), x), x)


def replacement287(a, b, c, d, e, x):
    return Dist(S(2)*e, Subst(Int(S(1)/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)), x)


def replacement288(a, c, d, e, x):
    return Dist(S(2)*e, Subst(Int(S(1)/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)), x)


def replacement289(a, b, c, d, e, m, x):
    return Dist(S(1)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((d + e*x)**(m + S(1))*Simp(-b*e + c*d - c*e*x, x)/(a + b*x + c*x**S(2)), x), x) + Simp(e*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement290(a, c, d, e, m, x):
    return Dist(c/(a*e**S(2) + c*d**S(2)), Int((d - e*x)*(d + e*x)**(m + S(1))/(a + c*x**S(2)), x), x) + Simp(e*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))), x)


def replacement291(a, b, c, d, e, m, x):
    return Int(ExpandIntegrand((d + e*x)**m, S(1)/(a + b*x + c*x**S(2)), x), x)


def replacement292(a, c, d, e, m, x):
    return Int(ExpandIntegrand((d + e*x)**m, S(1)/(a + c*x**S(2)), x), x)


def replacement293(a, b, c, d, e, x):
    return Simp(-S(2)*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((-S(4)*a*c + b**S(2))*sqrt(a + b*x + c*x**S(2))), x)


def replacement294(a, c, d, e, x):
    return Simp((-a*e + c*d*x)/(a*c*sqrt(a + c*x**S(2))), x)


def replacement295(a, b, c, d, e, p, x):
    return -Dist((S(2)*p + S(3))*(-b*e + S(2)*c*d)/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp((a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement296(a, c, d, e, p, x):
    return Dist(d*(S(2)*p + S(3))/(S(2)*a*(p + S(1))), Int((a + c*x**S(2))**(p + S(1)), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))), x)


def replacement297(a, b, c, d, e, p, x):
    return Dist((-b*e + S(2)*c*d)/(S(2)*c), Int((a + b*x + c*x**S(2))**p, x), x) + Simp(e*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))), x)


def replacement298(a, c, d, e, p, x):
    return Dist(d, Int((a + c*x**S(2))**p, x), x) + Simp(e*(a + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))), x)


def replacement299(a, b, c, d, e, m, p, x):
    return Dist((d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((d + e*x)**(m - p)*(a*d + c*e*x**S(3))**p, x), x)


def replacement300(b, c, d, e, m, x):
    return Int((d + e*x)**m/(sqrt(b*x)*sqrt(S(1) + c*x/b)), x)


def replacement301(b, c, d, e, m, x):
    return Dist(sqrt(x)*sqrt(b + c*x)/sqrt(b*x + c*x**S(2)), Int((d + e*x)**m/(sqrt(x)*sqrt(b + c*x)), x), x)


def replacement302(a, b, c, m, x):
    return Dist(S(2), Subst(Int(x**(S(2)*m + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x, sqrt(x)), x)


def replacement303(a, b, c, e, m, x):
    return Dist(x**(-m)*(e*x)**m, Int(x**m/sqrt(a + b*x + c*x**S(2)), x), x)


def replacement304(a, b, c, d, e, m, x):
    return Dist(S(2)*sqrt(-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))*(S(2)*c*(d + e*x)/(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2))))**(-m)*(d + e*x)**m*Rt(-S(4)*a*c + b**S(2), S(2))/(c*sqrt(a + b*x + c*x**S(2))), Subst(Int((S(2)*e*x**S(2)*Rt(-S(4)*a*c + b**S(2), S(2))/(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**m/sqrt(S(1) - x**S(2)), x), x, sqrt(S(2))*sqrt((b + S(2)*c*x + Rt(-S(4)*a*c + b**S(2), S(2)))/Rt(-S(4)*a*c + b**S(2), S(2)))/S(2)), x)


def replacement305(a, c, d, e, m, x):
    return Dist(S(2)*a*(c*(d + e*x)/(-a*e*Rt(-c/a, S(2)) + c*d))**(-m)*sqrt(S(1) + c*x**S(2)/a)*(d + e*x)**m*Rt(-c/a, S(2))/(c*sqrt(a + c*x**S(2))), Subst(Int((S(2)*a*e*x**S(2)*Rt(-c/a, S(2))/(-a*e*Rt(-c/a, S(2)) + c*d) + S(1))**m/sqrt(S(1) - x**S(2)), x), x, sqrt(-x*Rt(-c/a, S(2))/S(2) + S(1)/2)), x)


def replacement306(a, b, c, d, e, m, p, x):
    return Dist(p*(-S(4)*a*c + b**S(2))/(S(2)*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x), x) - Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/(S(2)*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement307(a, c, d, e, m, p, x):
    return -Dist(S(2)*a*c*p/((m + S(1))*(a*e**S(2) + c*d**S(2))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2)), x), x) - Simp((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-S(2)*a*e + S(2)*c*d*x)/(S(2)*(m + S(1))*(a*e**S(2) + c*d**S(2))), x)


def replacement308(a, b, c, d, e, m, p, x):
    return -Dist(S(2)*(S(2)*p + S(3))*(a*e**S(2) - b*d*e + c*d**S(2))/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement309(a, c, d, e, m, p, x):
    return Dist((S(2)*p + S(3))*(a*e**S(2) + c*d**S(2))/(S(2)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2)), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))), x)


def replacement310(a, b, c, d, e, x):
    return Dist(S(-2), Subst(Int(S(1)/(S(4)*a*e**S(2) - S(4)*b*d*e + S(4)*c*d**S(2) - x**S(2)), x), x, (S(2)*a*e - b*d - x*(-b*e + S(2)*c*d))/sqrt(a + b*x + c*x**S(2))), x)


def replacement311(a, c, d, e, x):
    return -Subst(Int(S(1)/(a*e**S(2) + c*d**S(2) - x**S(2)), x), x, (a*e - c*d*x)/sqrt(a + c*x**S(2)))


def replacement312(a, b, c, d, e, m, p, x):
    return -Simp(((b + S(2)*c*x + Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d + e*Rt(-S(4)*a*c + b**S(2), S(2)))/((b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2)))))**(-p)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*Hypergeometric2F1(m + S(1), -p, m + S(2), -S(4)*c*(d + e*x)*Rt(-S(4)*a*c + b**S(2), S(2))/((b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2)))))/((m + S(1))*(-b*e + S(2)*c*d + e*Rt(-S(4)*a*c + b**S(2), S(2)))), x)


def replacement313(a, c, d, e, m, p, x):
    return Simp(((c*d + e*Rt(-a*c, S(2)))*(c*x + Rt(-a*c, S(2)))/((c*d - e*Rt(-a*c, S(2)))*(c*x - Rt(-a*c, S(2)))))**(-p)*(a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-c*x + Rt(-a*c, S(2)))*Hypergeometric2F1(m + S(1), -p, m + S(2), S(2)*c*(d + e*x)*Rt(-a*c, S(2))/((c*d - e*Rt(-a*c, S(2)))*(-c*x + Rt(-a*c, S(2)))))/((m + S(1))*(c*d + e*Rt(-a*c, S(2)))), x)


def replacement314(a, b, c, d, e, m, p, x):
    return Dist(m*(-b*e + S(2)*c*d)/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp((b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement315(a, c, d, e, m, p, x):
    return -Dist(d*m/(S(2)*a*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x), x) - Simp(x*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*(p + S(1))), x)


def replacement316(a, b, c, d, e, m, p, x):
    return Dist((-b*e + S(2)*c*d)/(S(2)*a*e**S(2) - S(2)*b*d*e + S(2)*c*d**S(2)), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x), x) + Simp(e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement317(a, c, d, e, m, p, x):
    return Dist(c*d/(a*e**S(2) + c*d**S(2)), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x), x) + Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))), x)


def replacement318(a, b, c, d, e, m, p, x):
    return -Dist(p/(e*(m + S(1))), Int((b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))), x)


def replacement319(a, c, d, e, m, p, x):
    return -Dist(S(2)*c*p/(e*(m + S(1))), Int(x*(a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1)), x), x) + Simp((a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(1))), x)


def replacement320(a, b, c, d, e, m, p, x):
    return -Dist(p/(e*(m + S(2)*p + S(1))), Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d), x), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))), x)


def replacement321(a, c, d, e, m, p, x):
    return Dist(S(2)*p/(e*(m + S(2)*p + S(1))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**m*Simp(a*e - c*d*x, x), x), x) + Simp((a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(2)*p + S(1))), x)


def replacement322(a, b, c, d, e, m, p, x):
    return -Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(b*e*m + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*e*x*(m + S(2)*p + S(3))), x), x) + Simp((b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement323(a, c, d, e, m, p, x):
    return Dist(S(1)/(S(2)*a*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(d*(S(2)*p + S(3)) + e*x*(m + S(2)*p + S(3))), x), x) - Simp(x*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*(p + S(1))), x)


def replacement324(a, b, c, d, e, m, p, x):
    return Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-S(2)*c*d**S(2)*(S(2)*p + S(3)) + e*x*(b*e - S(2)*c*d)*(m + S(2)*p + S(2)) + e*(S(2)*a*e*(m + S(-1)) + b*d*(-m + S(2)*p + S(4))), x), x), x) + Simp((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement325(a, c, d, e, m, p, x):
    return Dist(-S(1)/(S(2)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*Simp(a*e**S(2)*(m + S(-1)) - c*d**S(2)*(S(2)*p + S(3)) - c*d*e*x*(m + S(2)*p + S(2)), x), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))), x)


def replacement326(a, b, c, d, e, m, p, x):
    return Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-S(2)*a*c*e**S(2)*(m + S(2)*p + S(3)) + b**S(2)*e**S(2)*(m + p + S(2)) + b*c*d*e*(-m + S(2)*p + S(2)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(3)) - c*e*x*(-b*e + S(2)*c*d)*(m + S(2)*p + S(4)), x), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*e - b**S(2)*e + b*c*d + c*x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement327(a, c, d, e, m, p, x):
    return Dist(S(1)/(S(2)*a*(p + S(1))*(a*e**S(2) + c*d**S(2))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*Simp(a*e**S(2)*(m + S(2)*p + S(3)) + c*d**S(2)*(S(2)*p + S(3)) + c*d*e*x*(m + S(2)*p + S(4)), x), x), x) - Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(a*e + c*d*x)/(S(2)*a*(p + S(1))*(a*e**S(2) + c*d**S(2))), x)


def replacement328(a, b, c, d, e, m, p, x):
    return Dist(S(1)/(c*(m + S(2)*p + S(1))), Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**p*Simp(c*d**S(2)*(m + S(2)*p + S(1)) + e*x*(m + p)*(-b*e + S(2)*c*d) - e*(a*e*(m + S(-1)) + b*d*(p + S(1))), x), x), x) + Simp(e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))), x)


def replacement329(a, c, d, e, m, p, x):
    return Dist(S(1)/(c*(m + S(2)*p + S(1))), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-2))*Simp(-a*e**S(2)*(m + S(-1)) + c*d**S(2)*(m + S(2)*p + S(1)) + S(2)*c*d*e*x*(m + p), x), x), x) + Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))), x)


def replacement330(a, b, c, d, e, m, p, x):
    return Dist(S(1)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*e*(m + p + S(2)) + c*d*(m + S(1)) - c*e*x*(m + S(2)*p + S(3)), x), x), x) + Simp(e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement331(a, c, d, e, m, p, x):
    return Dist(c/((m + S(1))*(a*e**S(2) + c*d**S(2))), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(d*(m + S(1)) - e*x*(m + S(2)*p + S(3)), x), x), x) + Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))), x)


def With332(a, b, c, d, e, x):
    q = Rt(S(3)*c*e**S(2)*(-b*e + S(2)*c*d), S(3))
    return -Simp(S(3)*c*e*log(d + e*x)/(S(2)*q**S(2)), x) + Simp(S(3)*c*e*log(-b*e + c*d - c*e*x - q*(a + b*x + c*x**S(2))**(S(1)/3))/(S(2)*q**S(2)), x) - Simp(sqrt(S(3))*c*e*ArcTan(sqrt(S(3))/S(3) + S(2)*sqrt(S(3))*(-b*e + c*d - c*e*x)/(S(3)*q*(a + b*x + c*x**S(2))**(S(1)/3)))/q**S(2), x)


def With333(a, c, d, e, x):
    q = Rt(S(6)*c**S(2)*e**S(2)/d**S(2), S(3))
    return -Simp(S(3)*c*e*log(d + e*x)/(S(2)*d**S(2)*q**S(2)), x) + Simp(S(3)*c*e*log(c*d - c*e*x - d*q*(a + c*x**S(2))**(S(1)/3))/(S(2)*d**S(2)*q**S(2)), x) - Simp(sqrt(S(3))*c*e*ArcTan(S(2)*sqrt(S(3))*c*(d - e*x)/(S(3)*d*q*(a + c*x**S(2))**(S(1)/3)) + sqrt(S(3))/S(3))/(d**S(2)*q**S(2)), x)


def With334(a, b, c, d, e, x):
    q = Rt(-S(3)*c*e**S(2)*(-b*e + S(2)*c*d), S(3))
    return -Simp(S(3)*c*e*log(d + e*x)/(S(2)*q**S(2)), x) + Simp(S(3)*c*e*log(-b*e + c*d - c*e*x + q*(a + b*x + c*x**S(2))**(S(1)/3))/(S(2)*q**S(2)), x) - Simp(sqrt(S(3))*c*e*ArcTan(sqrt(S(3))/S(3) - S(2)*sqrt(S(3))*(-b*e + c*d - c*e*x)/(S(3)*q*(a + b*x + c*x**S(2))**(S(1)/3)))/q**S(2), x)


def With335(a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist((b + S(2)*c*x - q)**(S(1)/3)*(b + S(2)*c*x + q)**(S(1)/3)/(a + b*x + c*x**S(2))**(S(1)/3), Int(S(1)/((d + e*x)*(b + S(2)*c*x - q)**(S(1)/3)*(b + S(2)*c*x + q)**(S(1)/3)), x), x)


def replacement336(a, c, d, e, x):
    return Dist(d, Int(S(1)/((a + c*x**S(2))**(S(1)/4)*(d**S(2) - e**S(2)*x**S(2))), x), x) - Dist(e, Int(x/((a + c*x**S(2))**(S(1)/4)*(d**S(2) - e**S(2)*x**S(2))), x), x)


def replacement337(a, c, d, e, x):
    return Dist(d, Int(S(1)/((a + c*x**S(2))**(S(3)/4)*(d**S(2) - e**S(2)*x**S(2))), x), x) - Dist(e, Int(x/((a + c*x**S(2))**(S(3)/4)*(d**S(2) - e**S(2)*x**S(2))), x), x)


def replacement338(a, b, c, d, e, p, x):
    return Dist((-S(4)*c/(-S(4)*a*c + b**S(2)))**(-p), Subst(Int(Simp(-x**S(2)/(-S(4)*a*c + b**S(2)) + S(1), x)**p/Simp(-b*e + S(2)*c*d + e*x, x), x), x, b + S(2)*c*x), x)


def replacement339(a, b, c, d, e, p, x):
    return Dist((-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))**(-p)*(a + b*x + c*x**S(2))**p, Int((-a*c/(-S(4)*a*c + b**S(2)) - b*c*x/(-S(4)*a*c + b**S(2)) - c**S(2)*x**S(2)/(-S(4)*a*c + b**S(2)))**p/(d + e*x), x), x)


def replacement340(a, c, d, e, m, p, x):
    return Int((d + e*x)**m*(-x*Rt(-c, S(2)) + Rt(a, S(2)))**p*(x*Rt(-c, S(2)) + Rt(a, S(2)))**p, x)


def With341(a, b, c, d, e, m, p, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return -Dist((e*(b + S(2)*c*x - q)/(S(2)*c*(d + e*x)))**(-p)*(e*(b + S(2)*c*x + q)/(S(2)*c*(d + e*x)))**(-p)*(a + b*x + c*x**S(2))**p*(S(1)/(d + e*x))**(S(2)*p)/e, Subst(Int(x**(-m - S(2)*p + S(-2))*Simp(-x*(d - e*(b - q)/(S(2)*c)) + S(1), x)**p*Simp(-x*(d - e*(b + q)/(S(2)*c)) + S(1), x)**p, x), x, S(1)/(d + e*x)), x)


def With342(a, c, d, e, m, p, x):
    q = Rt(-a*c, S(2))
    return -Dist((e*(c*x + q)/(c*(d + e*x)))**(-p)*(-e*(-c*x + q)/(c*(d + e*x)))**(-p)*(a + c*x**S(2))**p*(S(1)/(d + e*x))**(S(2)*p)/e, Subst(Int(x**(-m - S(2)*p + S(-2))*Simp(-x*(d - e*q/c) + S(1), x)**p*Simp(-x*(d + e*q/c) + S(1), x)**p, x), x, S(1)/(d + e*x)), x)


def With343(a, b, c, d, e, m, p, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist((-(d + e*x)/(d - e*(b - q)/(S(2)*c)) + S(1))**(-p)*(-(d + e*x)/(d - e*(b + q)/(S(2)*c)) + S(1))**(-p)*(a + b*x + c*x**S(2))**p/e, Subst(Int(x**m*Simp(-x/(d - e*(b - q)/(S(2)*c)) + S(1), x)**p*Simp(-x/(d - e*(b + q)/(S(2)*c)) + S(1), x)**p, x), x, d + e*x), x)


def With344(a, c, d, e, m, p, x):
    q = Rt(-a*c, S(2))
    return Dist((a + c*x**S(2))**p*(-(d + e*x)/(d - e*q/c) + S(1))**(-p)*(-(d + e*x)/(d + e*q/c) + S(1))**(-p)/e, Subst(Int(x**m*Simp(-x/(d - e*q/c) + S(1), x)**p*Simp(-x/(d + e*q/c) + S(1), x)**p, x), x, d + e*x), x)


def replacement345(a, b, c, d, e, m, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x, u), x)


def replacement346(a, c, d, e, m, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + c*x**S(2))**p*(d + e*x)**m, x), x, u), x)


def replacement347(a, c, d, e, n, p, x):
    return Dist(d, Int(x**n*(a + c*x**S(2))**p, x), x) + Dist(e, Int(x**(n + S(1))*(a + c*x**S(2))**p, x), x)


def replacement348(a, b, c, d, e, f, g, m, x):
    return Dist((f + g*x)/sqrt(a + b*x + c*x**S(2)), Int((d + e*x)**m, x), x)


def replacement349(a, b, c, d, e, f, g, m, p, x):
    return -Simp(f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))*(-d*g + e*f)), x)


def replacement350(a, b, c, d, e, f, g, m, p, x):
    return -Dist(e*g*m/(S(2)*c*(p + S(1))), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp(g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))), x)


def replacement351(a, b, c, d, e, f, g, m, p, x):
    return Dist(e*f*g*(m + S(2)*p + S(3))/(b*(p + S(1))*(-d*g + e*f)), Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x), x) - Simp(f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))*(-d*g + e*f)), x)


def replacement352(a, b, c, d, e, f, g, m, p, x):
    return -Dist(g*(S(2)*p + S(1))/(e*(m + S(1))), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x), x) + Simp((d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p/(e*(m + S(1))), x)


def replacement353(a, b, c, d, e, f, g, m, p, x):
    return -Dist(g*(m + S(2)*p + S(3))/((m + S(1))*(-d*g + e*f)), Int((d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p, x), x) + Simp(S(2)*f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(m + S(1))*(-d*g + e*f)), x)


def replacement354(a, b, c, d, e, f, g, m, p, x):
    return -Dist(b*m*(-d*g + e*f)/(S(2)*c*f*(m + S(2)*p + S(2))), Int((d + e*x)**(m + S(-1))*(f + g*x)*(a + b*x + c*x**S(2))**p, x), x) + Simp(g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))), x)


def replacement355(a, b, c, d, e, f, g, m, p, x):
    return Dist((S(2)*p + S(1))*(-d*g + e*f)/(e*(m + S(2)*p + S(2))), Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x) + Simp((d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(2))), x)


def replacement356(a, b, c, d, e, f, g, p, x):
    return Dist((-b*g + S(2)*c*f)/(-b*e + S(2)*c*d), Int((a + b*x + c*x**S(2))**p, x), x) - Dist((-d*g + e*f)/(-b*e + S(2)*c*d), Int((b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(d + e*x), x), x)


def replacement357(a, b, c, d, e, f, g, m, p, x):
    return Dist(g/e, Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x), x) + Dist((-d*g + e*f)/e, Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x)


def replacement358(a, b, c, d, e, f, g, m, p, x):
    return Dist((-b*g + S(2)*c*f)/(-b*e + S(2)*c*d), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x), x) - Dist((-d*g + e*f)/(-b*e + S(2)*c*d), Int((b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x)


def replacement359(a, b, c, d, e, f, g, m, p, x):
    return Dist((S(2)*c*e*f*(m + S(2)*p + S(2)) - g*(b*e*(m + S(1)) + S(2)*c*d*(S(2)*p + S(1))))/(e*(m + S(1))*(-b*e + S(2)*c*d)), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x), x) - Simp((b + S(2)*c*x)*(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**p/(e*(m + S(1))*(-b*e + S(2)*c*d)), x)


def replacement360(a, b, c, d, e, f, g, m, p, x):
    return Dist((S(2)*c*e*f*(m + S(2)*p + S(2)) - g*(b*e*(m + S(1)) + S(2)*c*(S(2)*d*p + d)))/(S(2)*c*e*(m + S(2)*p + S(2))), Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x) + Simp(g*(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(S(2)*c*e*(m + S(2)*p + S(2))), x)


def replacement361(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(c**(-IntPart(p))*(b/S(2) + c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((b/S(2) + c*x)**(S(2)*p)*(d + e*x)**m*(f + g*x)**n, x), x)


def replacement362(a, b, c, d, e, f, g, m, n, p, x):
    return Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x)


def replacement363(a, c, d, e, f, g, m, n, p, x):
    return Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x)


def replacement364(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(d**m*e**m, Int((f + g*x)**n*(a*e + c*d*x)**(-m)*(a + b*x + c*x**S(2))**(m + p), x), x)


def replacement365(a, c, d, e, f, g, m, n, p, x):
    return Dist(d**m*e**m, Int((a + c*x**S(2))**(m + p)*(f + g*x)**n*(a*e + c*d*x)**(-m), x), x)


def replacement366(a, b, c, d, e, f, g, m, p, x):
    return Simp(g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))), x)


def replacement367(a, c, d, e, f, g, m, p, x):
    return Simp(g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))), x)


def replacement368(a, b, c, d, e, f, g, m, p, x):
    return -Dist(e*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))/(c*(p + S(1))*(-b*e + S(2)*c*d)), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp((d + e*x)**m*(c*e*f + g*(-b*e + c*d))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))*(-b*e + S(2)*c*d)), x)


def replacement369(a, c, d, e, f, g, m, p, x):
    return -Dist(e*(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))/(S(2)*c*d*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g + e*f)/(S(2)*c*d*(p + S(1))), x)


def replacement370(a, b, c, d, e, f, g, m, p, x):
    return -Dist(e*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))/(c*(p + S(1))*(-b*e + S(2)*c*d)), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp((d + e*x)**m*(c*e*f + g*(-b*e + c*d))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))*(-b*e + S(2)*c*d)), x)


def replacement371(a, c, d, e, f, g, m, p, x):
    return -Dist(e*(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))/(S(2)*c*d*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g + e*f)/(S(2)*c*d*(p + S(1))), x)


def replacement372(a, b, c, d, e, f, g, m, p, x):
    return Dist((e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))/(e*(-b*e + S(2)*c*d)*(m + p + S(1))), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x), x) + Simp((d + e*x)**m*(d*g - e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))), x)


def replacement373(a, c, d, e, f, g, m, p, x):
    return Dist((S(2)*c*e*f*(p + S(1)) + m*(c*d*g + c*e*f))/(S(2)*c*d*e*(m + p + S(1))), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g - e*f)/(S(2)*c*d*(m + p + S(1))), x)


def replacement374(a, b, c, d, e, f, g, m, p, x):
    return Dist((e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))/(c*e*(m + S(2)*p + S(2))), Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x) + Simp(g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))), x)


def replacement375(a, c, d, e, f, g, m, p, x):
    return Dist((S(2)*e*f*(p + S(1)) + m*(d*g + e*f))/(e*(m + S(2)*p + S(2))), Int((a + c*x**S(2))**p*(d + e*x)**m, x), x) + Simp(g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))), x)


def replacement376(a, c, f, g, p, x):
    return -Dist(S(1)/(S(2)*a*c*(p + S(1))), Int(x*(a + c*x**S(2))**(p + S(1))*Simp(S(2)*a*g - c*f*x*(S(2)*p + S(5)), x), x), x) + Simp(x**S(2)*(a + c*x**S(2))**(p + S(1))*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))), x)


def replacement377(a, c, f, g, p, x):
    return Dist(S(1)/c, Int((a + c*x**S(2))**(p + S(1))*(f + g*x), x), x) - Dist(f**S(2)/c, Int((a + c*x**S(2))**(p + S(1))/(f - g*x), x), x)


def replacement378(a, b, c, d, e, f, g, m, n, p, x):
    return Int((f + g*x)**n*(a/d + c*x/e)**(-m)*(a + b*x + c*x**S(2))**(m + p), x)


def replacement379(a, c, d, e, f, g, m, n, p, x):
    return Dist(a**(-m)*d**(S(2)*m), Int((a + c*x**S(2))**(m + p)*(d - e*x)**(-m)*(f + g*x)**n, x), x)


def replacement380(a, b, c, d, e, f, g, n, p, x):
    return -Dist(S(1)/(d*e*p*(-S(4)*a*c + b**S(2))), Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p*Simp(-S(2)*a*c*(d*g*n - e*f*(S(2)*p + S(1))) + b*(a*e*g*n - c*d*f*(S(2)*p + S(1))) - c*g*x*(-S(2)*a*e + b*d)*(n + S(2)*p + S(1)), x), x), x) - Simp((f + g*x)**n*(a*(-b*e + S(2)*c*d) + c*x*(-S(2)*a*e + b*d))*(a + b*x + c*x**S(2))**p/(d*e*p*(-S(4)*a*c + b**S(2))), x)


def replacement381(a, c, d, e, f, g, n, p, x):
    return -Dist(S(1)/(S(2)*d*e*p), Int((a + c*x**S(2))**p*(f + g*x)**(n + S(-1))*Simp(d*g*n - e*f*(S(2)*p + S(1)) - e*g*x*(n + S(2)*p + S(1)), x), x), x) + Simp((a + c*x**S(2))**p*(d - e*x)*(f + g*x)**n/(S(2)*d*e*p), x)


def replacement382(a, b, c, d, e, f, g, n, p, x):
    return -Dist(S(1)/(d*e*p*(-S(4)*a*c + b**S(2))*(a*g**S(2) - b*f*g + c*f**S(2))), Int((f + g*x)**n*(a + b*x + c*x**S(2))**p*Simp(S(2)*a*c*(a*e*g**S(2)*(n + S(2)*p + S(1)) + c*f*(-d*g*n + S(2)*e*f*p + e*f)) + b**S(2)*g*(-a*e*g*(n + p + S(1)) + c*d*f*p) + b*c*(a*g*(d*g*(n + S(1)) + e*f*(n - S(2)*p)) - c*d*f**S(2)*(S(2)*p + S(1))) + c*g*x*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f))*(n + S(2)*p + S(2)), x), x), x) - Simp((f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p*(a*c*d*(-b*g + S(2)*c*f) - a*e*(S(2)*a*c*g - b**S(2)*g + b*c*f) + c*x*(-a*e*(-b*g + S(2)*c*f) + c*d*(-S(2)*a*g + b*f)))/(d*e*p*(-S(4)*a*c + b**S(2))*(a*g**S(2) - b*f*g + c*f**S(2))), x)


def replacement383(a, c, d, e, f, g, n, p, x):
    return Dist(S(1)/(S(2)*d*e*p*(a*g**S(2) + c*f**S(2))), Int((a + c*x**S(2))**p*(f + g*x)**n*Simp(a*e*g**S(2)*(n + S(2)*p + S(1)) - c*f*(d*g*n - e*(S(2)*f*p + f)) + c*g*x*(d*g + e*f)*(n + S(2)*p + S(2)), x), x), x) + Simp((a + c*x**S(2))**p*(f + g*x)**(n + S(1))*(-a*e*g + c*d*f - c*x*(d*g + e*f))/(S(2)*d*e*p*(a*g**S(2) + c*f**S(2))), x)


def replacement384(a, b, c, d, e, f, g, m, n, p, x):
    return -Simp(e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m - n + S(-1))), x)


def replacement385(a, c, d, e, f, g, m, n, p, x):
    return -Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(m - n + S(-1))), x)


def replacement386(a, b, c, d, e, f, g, m, n, p, x):
    return -Simp(e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((n + S(1))*(-b*e*g + c*d*g + c*e*f)), x)


def replacement387(a, c, d, e, f, g, m, n, p, x):
    return -Simp(e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/(c*(n + S(1))*(d*g + e*f)), x)


def replacement388(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(c*m/(e*g*(n + S(1))), Int((d + e*x)**(m + S(1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x), x) + Simp((d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p/(g*(n + S(1))), x)


def replacement389(a, c, d, e, f, g, m, n, p, x):
    return Dist(c*m/(e*g*(n + S(1))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*(f + g*x)**(n + S(1)), x), x) + Simp((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1))/(g*(n + S(1))), x)


def replacement390(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(m*(-b*e*g + c*d*g + c*e*f)/(e**S(2)*g*(m - n + S(-1))), Int((d + e*x)**(m + S(1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(-1)), x), x) - Simp((d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p/(g*(m - n + S(-1))), x)


def replacement391(a, c, d, e, f, g, m, n, p, x):
    return -Dist(c*m*(d*g + e*f)/(e**S(2)*g*(m - n + S(-1))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*(f + g*x)**n, x), x) - Simp((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1))/(g*(m - n + S(-1))), x)


def replacement392(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(e*g*n/(c*(p + S(1))), Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp(e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))), x)


def replacement393(a, c, d, e, f, g, m, n, p, x):
    return -Dist(e*g*n/(c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x), x) + Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(p + S(1))), x)


def replacement394(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(e**S(2)*g*(m - n + S(-2))/((p + S(1))*(-b*e*g + c*d*g + c*e*f)), Int((d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp(e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-b*e*g + c*d*g + c*e*f)), x)


def replacement395(a, c, d, e, f, g, m, n, p, x):
    return Dist(e**S(2)*g*(m - n + S(-2))/(c*(p + S(1))*(d*g + e*f)), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n, x), x) + Simp(e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/(c*(p + S(1))*(d*g + e*f)), x)


def replacement396(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(n*(-b*e*g + c*d*g + c*e*f)/(c*e*(m - n + S(-1))), Int((d + e*x)**m*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x), x) - Simp(e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m - n + S(-1))), x)


def replacement397(a, c, d, e, f, g, m, n, p, x):
    return -Dist(n*(d*g + e*f)/(e*(m - n + S(-1))), Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(-1)), x), x) - Simp(e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(m - n + S(-1))), x)


def replacement398(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(c*e*(m - n + S(-2))/((n + S(1))*(-b*e*g + c*d*g + c*e*f)), Int((d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p, x), x) - Simp(e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((n + S(1))*(-b*e*g + c*d*g + c*e*f)), x)


def replacement399(a, c, d, e, f, g, m, n, p, x):
    return -Dist(e*(m - n + S(-2))/((n + S(1))*(d*g + e*f)), Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1)), x), x) - Simp(e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/((n + S(1))*(c*d*g + c*e*f)), x)


def replacement400(a, b, c, d, e, f, g, x):
    return Dist(S(2)*e**S(2), Subst(Int(S(1)/(-b*e*g + c*(d*g + e*f) + e**S(2)*g*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))/sqrt(d + e*x)), x)


def replacement401(a, c, d, e, f, g, x):
    return Dist(S(2)*e**S(2), Subst(Int(S(1)/(c*(d*g + e*f) + e**S(2)*g*x**S(2)), x), x, sqrt(a + c*x**S(2))/sqrt(d + e*x)), x)


def replacement402(a, b, c, d, e, f, g, m, n, p, x):
    return Simp(e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*g*(n + p + S(2))), x)


def replacement403(a, c, d, e, f, g, m, n, p, x):
    return Simp(e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))/(c*g*(n + p + S(2))), x)


def replacement404(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(e*(b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))/(g*(n + S(1))*(-b*e*g + c*d*g + c*e*f)), Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p, x), x) + Simp(e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(g*(n + S(1))*(-b*e*g + c*d*g + c*e*f)), x)


def replacement405(a, c, d, e, f, g, m, n, p, x):
    return -Dist(e*(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))/(g*(n + S(1))*(d*g + e*f)), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1)), x), x) + Simp(e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(-d*g + e*f)/(c*g*(n + S(1))*(d*g + e*f)), x)


def replacement406(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist((b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))/(c*g*(n + p + S(2))), Int((d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x) + Simp(e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*g*(n + p + S(2))), x)


def replacement407(a, c, d, e, f, g, m, n, p, x):
    return -Dist((-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))/(g*(n + p + S(2))), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*(f + g*x)**n, x), x) + Simp(e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))/(c*g*(n + p + S(2))), x)


def replacement408(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x)


def replacement409(a, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand(S(1)/sqrt(a + c*x**S(2)), (a + c*x**S(2))**(p + S(1)/2)*(d + e*x)**m*(f + g*x)**n, x), x)


def replacement410(a, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x)


def replacement411(a, b, c, d, e, p, x):
    return -Dist(e**(S(-2)), Int((d - e*x)*(a + b*x + c*x**S(2))**p, x), x) + Dist(d**S(2)/e**S(2), Int((a + b*x + c*x**S(2))**p/(d + e*x), x), x)


def replacement412(a, c, d, e, p, x):
    return -Dist(e**(S(-2)), Int((a + c*x**S(2))**p*(d - e*x), x), x) + Dist(d**S(2)/e**S(2), Int((a + c*x**S(2))**p/(d + e*x), x), x)


def replacement413(a, b, c, d, e, f, g, m, p, x):
    return -Dist(S(1)/(c*e**S(2)*(m + S(2)*p + S(3))), Int((d + e*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*e*g*(d*g + e*f*(m + p + S(1))) - c*(d**S(2)*g**S(2) + d*e*f*g*m + e**S(2)*f**S(2)*(m + S(2)*p + S(3))) + e*g*x*(b*e*g*(m + p + S(2)) - c*(d*g*m + e*f*(m + S(2)*p + S(4)))), x), x), x) + Simp(g*(d + e*x)**m*(f + g*x)*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(3))), x)


def replacement414(a, c, d, e, f, g, m, p, x):
    return -Dist(S(1)/(c*e**S(2)*(m + S(2)*p + S(3))), Int((a + c*x**S(2))**p*(d + e*x)**m*Simp(-c*e*g*x*(d*g*m + e*f*(m + S(2)*p + S(4))) - c*(d**S(2)*g**S(2) + d*e*f*g*m + e**S(2)*f**S(2)*(m + S(2)*p + S(3))), x), x), x) + Simp(g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(f + g*x)/(c*(m + S(2)*p + S(3))), x)


def replacement415(b, c, e, f, g, m, n, p, x):
    return Dist(x**(-m - p)*(e*x)**m*(b + c*x)**(-p)*(b*x + c*x**S(2))**p, Int(x**(m + p)*(b + c*x)**p*(f + g*x)**n, x), x)


def replacement416(a, c, d, e, f, g, m, n, p, x):
    return Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x)


def replacement417(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x), x)


def replacement418(a, c, d, e, f, g, m, n, p, x):
    return Dist((a + c*x**S(2))**FracPart(p)*(d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p)), Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x), x)


def replacement419(a, b, c, d, e, f, g, m, p, x):
    return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)*(a + b*x + c*x**S(2))**p, x), x)


def replacement420(a, c, d, e, f, g, m, p, x):
    return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x), x), x)


def replacement421(a, b, c, d, e, f, g, x):
    return Dist(e*(-d*g + e*f)/(a*e**S(2) - b*d*e + c*d**S(2)), Int(S(1)/(d + e*x), x), x) + Dist(S(1)/(a*e**S(2) - b*d*e + c*d**S(2)), Int(Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x), x)


def replacement422(a, c, d, e, f, g, x):
    return Dist(e*(-d*g + e*f)/(a*e**S(2) + c*d**S(2)), Int(S(1)/(d + e*x), x), x) + Dist(S(1)/(a*e**S(2) + c*d**S(2)), Int(Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x), x)


def replacement423(a, b, c, d, e, f, g, m, p, x):
    return -Simp((d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*(p + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement424(a, c, d, e, f, g, m, p, x):
    return -Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/(S(2)*(p + S(1))*(a*e**S(2) + c*d**S(2))), x)


def replacement425(a, b, c, d, e, f, g, m, p, x):
    return -Dist(m*(-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x), x) + Simp((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*g + b*f + x*(-b*g + S(2)*c*f))/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement426(a, c, d, e, f, g, m, p, x):
    return -Dist(m*(a*e*g + c*d*f)/(S(2)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))), x)


def replacement427(a, b, c, d, e, f, g, m, p, x):
    return -Dist((-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)/(S(2)*a*e**S(2) - S(2)*b*d*e + S(2)*c*d**S(2)), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x), x) - Simp((d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*(p + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement428(a, c, d, e, f, g, m, p, x):
    return Dist((a*e*g + c*d*f)/(a*e**S(2) + c*d**S(2)), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x), x) - Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/(S(2)*(p + S(1))*(a*e**S(2) + c*d**S(2))), x)


def replacement429(a, b, c, d, e, f, g, p, x):
    return -Simp((a + b*x + c*x**S(2))**(p + S(1))*(b*e*g*(p + S(2)) - S(2)*c*e*g*x*(p + S(1)) - c*(S(2)*p + S(3))*(d*g + e*f))/(S(2)*c**S(2)*(p + S(1))*(S(2)*p + S(3))), x)


def replacement430(a, c, d, e, f, g, p, x):
    return Simp((a + c*x**S(2))**(p + S(1))*(S(2)*e*g*x*(p + S(1)) + (S(2)*p + S(3))*(d*g + e*f))/(S(2)*c*(p + S(1))*(S(2)*p + S(3))), x)


def replacement431(a, b, c, d, e, f, g, p, x):
    return -Dist((-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x + c*x**S(2))**(p + S(1)), x), x) - Simp((a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f) - x*(b**S(2)*e*g - b*c*(d*g + e*f) + S(2)*c*(-a*e*g + c*d*f)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement432(a, c, d, e, f, g, p, x):
    return -Dist((a*e*g - c*d*f*(S(2)*p + S(3)))/(S(2)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1)), x), x) - Simp((a + c*x**S(2))**(p + S(1))*(-a*(d*g + e*(f + g*x)) + c*d*f*x)/(S(2)*a*c*(p + S(1))), x)


def replacement433(a, b, c, d, e, f, g, p, x):
    return Dist((-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))/(S(2)*c**S(2)*(S(2)*p + S(3))), Int((a + b*x + c*x**S(2))**p, x), x) - Simp((a + b*x + c*x**S(2))**(p + S(1))*(b*e*g*(p + S(2)) - S(2)*c*e*g*x*(p + S(1)) - c*(S(2)*p + S(3))*(d*g + e*f))/(S(2)*c**S(2)*(p + S(1))*(S(2)*p + S(3))), x)


def replacement434(a, c, d, e, f, g, p, x):
    return -Dist((a*e*g - c*d*f*(S(2)*p + S(3)))/(c*(S(2)*p + S(3))), Int((a + c*x**S(2))**p, x), x) + Simp((a + c*x**S(2))**(p + S(1))*(S(2)*e*g*x*(p + S(1)) + (S(2)*p + S(3))*(d*g + e*f))/(S(2)*c*(p + S(1))*(S(2)*p + S(3))), x)


def replacement435(a, c, e, f, g, m, p, x):
    return Dist(f, Int((e*x)**m*(a + c*x**S(2))**p, x), x) + Dist(g/e, Int((e*x)**(m + S(1))*(a + c*x**S(2))**p, x), x)


def replacement436(a, b, c, d, e, f, g, m, p, x):
    return Dist((d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((f + g*x)*(a*d + c*e*x**S(3))**p, x), x)


def replacement437(a, b, c, d, e, f, g, m, p, x):
    return -Dist(p/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(S(2)*a*c*e*(m + S(2))*(-d*g + e*f) + b**S(2)*e*(d*g*(p + S(1)) - e*f*(m + p + S(2))) + b*(a*e**S(2)*g*(m + S(1)) - c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2)))) - c*x*(S(2)*c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2))) - e*(S(2)*a*e*g*(m + S(1)) - b*(d*g*(m - S(2)*p) + e*f*(m + S(2)*p + S(2))))), x), x), x) - Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-d*p*(-b*e + S(2)*c*d)*(-d*g + e*f) - e*x*(g*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2)) + p*(-b*e + S(2)*c*d)*(-d*g + e*f)) + (d*g - e*f*(m + S(2)))*(a*e**S(2) - b*d*e + c*d**S(2)))/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement438(a, c, d, e, f, g, m, p, x):
    return -Dist(p/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) + c*d**S(2))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2))*Simp(S(2)*a*c*e*(m + S(2))*(-d*g + e*f) - c*x*(-S(2)*a*e**S(2)*g*(m + S(1)) + S(2)*c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2)))), x), x), x) - Simp((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-S(2)*c*d**S(2)*p*(-d*g + e*f) - e*x*(S(2)*c*d*p*(-d*g + e*f) + g*(m + S(1))*(a*e**S(2) + c*d**S(2))) + (a*e**S(2) + c*d**S(2))*(d*g - e*f*(m + S(2))))/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) + c*d**S(2))), x)


def replacement439(a, b, c, d, e, f, g, m, p, x):
    return Dist(p/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-b*e*f*(m + S(2)*p + S(2)) + g*(S(2)*a*e*m + S(2)*a*e + S(2)*b*d*p + b*d) + x*(-S(2)*c*e*f*(m + S(2)*p + S(2)) + g*(b*e*m + b*e + S(4)*c*d*p + S(2)*c*d)), x), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-d*g*(S(2)*p + S(1)) + e*f*(m + S(2)*p + S(2)) + e*g*x*(m + S(1)))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))), x)


def replacement440(a, c, d, e, f, g, m, p, x):
    return Dist(p/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*Simp(g*(S(2)*a*e*m + S(2)*a*e) + x*(-S(2)*c*e*f*(m + S(2)*p + S(2)) + g*(S(4)*c*d*p + S(2)*c*d)), x), x), x) + Simp((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-d*g*(S(2)*p + S(1)) + e*f*(m + S(2)*p + S(2)) + e*g*x*(m + S(1)))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))), x)


def replacement441(a, b, c, d, e, f, g, m, p, x):
    return -Dist(p/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))), Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(c*e*f*(-S(2)*a*e + b*d)*(m + S(2)*p + S(2)) + g*(a*e*(b*e*m + b*e - S(2)*c*d*m) + b*d*(b*e*p - S(2)*c*d*p - c*d)) + x*(c*e*f*(-b*e + S(2)*c*d)*(m + S(2)*p + S(2)) + g*(b**S(2)*e**S(2)*(m + p + S(1)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(1)) - c*e*(S(2)*a*e*(m + S(2)*p + S(1)) + b*d*(m - S(2)*p)))), x), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(c*e*f*(m + S(2)*p + S(2)) + c*e*g*x*(m + S(2)*p + S(1)) - g*(-b*e*p + S(2)*c*d*p + c*d))/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))), x)


def replacement442(a, c, d, e, f, g, m, p, x):
    return Dist(S(2)*p/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**m*Simp(a*c*d*e*g*m + a*c*e**S(2)*f*(m + S(2)*p + S(2)) - x*(c**S(2)*d*e*f*(m + S(2)*p + S(2)) - g*(a*c*e**S(2)*(m + S(2)*p + S(1)) + c**S(2)*d**S(2)*(S(2)*p + S(1)))), x), x), x) + Simp((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-c*d*g*(S(2)*p + S(1)) + c*e*f*(m + S(2)*p + S(2)) + c*e*g*x*(m + S(2)*p + S(1)))/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))), x)


def replacement443(a, b, c, d, e, f, g, m, p, x):
    return -Dist(S(1)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(b*e*g*(a*e*(m + S(-1)) + b*d*(p + S(2))) + S(2)*c**S(2)*d**S(2)*f*(S(2)*p + S(3)) - c*(S(2)*a*e*(d*g*m + e*f*(m + S(-1))) + b*d*(d*g*(S(2)*p + S(3)) - e*f*(m - S(2)*p + S(-4)))) + e*x*(b**S(2)*e*g*(m + p + S(1)) + S(2)*c**S(2)*d*f*(m + S(2)*p + S(2)) - c*(S(2)*a*e*g*m + b*(d*g + e*f)*(m + S(2)*p + S(2)))), x), x), x) - Simp((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f) - x*(b**S(2)*e*g + S(2)*c**S(2)*d*f - c*(S(2)*a*e*g + b*d*g + b*e*f)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement444(a, c, d, e, f, g, m, p, x):
    return -Dist(S(1)/(S(4)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*Simp(S(2)*a*e*(d*g*m + e*f*(m + S(-1))) - S(2)*c*d**S(2)*f*(S(2)*p + S(3)) + e*x*(S(2)*a*e*g*m - S(2)*c*d*f*(m + S(2)*p + S(2))), x), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(S(2)*a*(d*g + e*f) - x*(-S(2)*a*e*g + S(2)*c*d*f))/(S(4)*a*c*(p + S(1))), x)


def replacement445(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-e*x*(-b*g + S(2)*c*f)*(m + S(2)*p + S(3)) - f*(b*e*m + S(2)*c*d*(S(2)*p + S(3))) + g*(S(2)*a*e*m + b*d*(S(2)*p + S(3))), x), x), x) + Simp((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*g + b*f + x*(-b*g + S(2)*c*f))/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement446(a, c, d, e, f, g, m, p, x):
    return -Dist(S(1)/(S(2)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*Simp(a*e*g*m - c*d*f*(S(2)*p + S(3)) - c*e*f*x*(m + S(2)*p + S(3)), x), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))), x)


def replacement447(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(c*e*x*(-f*(-b*e + S(2)*c*d) + g*(-S(2)*a*e + b*d))*(m + S(2)*p + S(4)) + f*(-S(2)*a*c*e**S(2)*(m + S(2)*p + S(3)) + b**S(2)*e**S(2)*(m + p + S(2)) + b*c*d*e*(-m + S(2)*p + S(2)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(3))) - g*(a*e*(b*e*m + b*e - S(2)*c*d*m) - b*d*(-b*e*p - b*e + S(2)*c*d*p + S(3)*c*d)), x), x), x) + Simp((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-a*g*(-b*e + S(2)*c*d) + c*x*(f*(-b*e + S(2)*c*d) - g*(-S(2)*a*e + b*d)) + f*(S(2)*a*c*e - b**S(2)*e + b*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement448(a, c, d, e, f, g, m, p, x):
    return Dist(S(1)/(S(2)*a*c*(p + S(1))*(a*e**S(2) + c*d**S(2))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*Simp(-a*c*d*e*g*m + c*e*x*(a*e*g + c*d*f)*(m + S(2)*p + S(4)) + f*(a*c*e**S(2)*(m + S(2)*p + S(3)) + c**S(2)*d**S(2)*(S(2)*p + S(3))), x), x), x) - Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-a*c*d*g + a*c*e*f + c*x*(a*e*g + c*d*f))/(S(2)*a*c*(p + S(1))*(a*e**S(2) + c*d**S(2))), x)


def replacement449(a, b, c, d, e, f, g, m, x):
    return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)/(a + b*x + c*x**S(2)), x), x)


def replacement450(a, c, d, e, f, g, m, x):
    return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)/(a + c*x**S(2)), x), x)


def replacement451(a, b, c, d, e, f, g, m, x):
    return Dist(S(1)/c, Int((d + e*x)**(m + S(-1))*Simp(-a*e*g + c*d*f + x*(-b*e*g + c*d*g + c*e*f), x)/(a + b*x + c*x**S(2)), x), x) + Simp(g*(d + e*x)**m/(c*m), x)


def replacement452(a, c, d, e, f, g, m, x):
    return Dist(S(1)/c, Int((d + e*x)**(m + S(-1))*Simp(-a*e*g + c*d*f + x*(c*d*g + c*e*f), x)/(a + c*x**S(2)), x), x) + Simp(g*(d + e*x)**m/(c*m), x)


def replacement453(a, b, c, d, e, f, g, x):
    return Dist(S(2), Subst(Int((-d*g + e*f + g*x**S(2))/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)), x)


def replacement454(a, c, d, e, f, g, x):
    return Dist(S(2), Subst(Int((-d*g + e*f + g*x**S(2))/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)), x)


def replacement455(a, b, c, d, e, f, g, m, x):
    return Dist(S(1)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((d + e*x)**(m + S(1))*Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x), x) + Simp((d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement456(a, c, d, e, f, g, m, x):
    return Dist(S(1)/(a*e**S(2) + c*d**S(2)), Int((d + e*x)**(m + S(1))*Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x), x) + Simp((d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))), x)


def replacement457(a, b, c, d, e, f, g, m, x):
    return Int(ExpandIntegrand((d + e*x)**m, (f + g*x)/(a + b*x + c*x**S(2)), x), x)


def replacement458(a, c, d, e, f, g, m, x):
    return Int(ExpandIntegrand((d + e*x)**m, (f + g*x)/(a + c*x**S(2)), x), x)


def replacement459(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/(c*(m + S(2)*p + S(2))), Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p*Simp(d*(p + S(1))*(-b*g + S(2)*c*f) + m*(-a*e*g + c*d*f) + x*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(-b*e*g + c*d*g + c*e*f)), x), x), x) + Simp(g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))), x)


def replacement460(a, c, d, e, f, g, m, p, x):
    return Dist(S(1)/(c*(m + S(2)*p + S(2))), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*Simp(-a*e*g*m + c*d*f*(m + S(2)*p + S(2)) + c*x*(d*g*m + e*f*(m + S(2)*p + S(2))), x), x), x) + Simp(g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))), x)


def replacement461(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(b*(p + S(1))*(d*g - e*f) - c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g - b*e*f + c*d*f), x), x), x) + Simp((d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement462(a, c, d, e, f, g, m, p, x):
    return Dist(S(1)/((m + S(1))*(a*e**S(2) + c*d**S(2))), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(-c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g + c*d*f), x), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))), x)


def replacement463(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(b*(p + S(1))*(d*g - e*f) - c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g - b*e*f + c*d*f), x), x), x) + Simp((d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement464(a, c, d, e, f, g, m, p, x):
    return Dist(S(1)/((m + S(1))*(a*e**S(2) + c*d**S(2))), Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(-c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g + c*d*f), x), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))), x)


def replacement465(a, b, c, d, e, f, g, x):
    return Dist(S(4)*f*(a - d)/(-a*e + b*d), Subst(Int(S(1)/(S(4)*a - S(4)*d - x**S(2)), x), x, (S(2)*a - S(2)*d + x*(b - e))/sqrt(a + b*x + c*x**S(2))), x)


def replacement466(a, b, c, f, g, x):
    return Dist(S(2), Subst(Int((f + g*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x, sqrt(x)), x)


def replacement467(a, c, f, g, x):
    return Dist(S(2), Subst(Int((f + g*x**S(2))/sqrt(a + c*x**S(4)), x), x, sqrt(x)), x)


def replacement468(a, b, c, e, f, g, x):
    return Dist(sqrt(x)/sqrt(e*x), Int((f + g*x)/(sqrt(x)*sqrt(a + b*x + c*x**S(2))), x), x)


def replacement469(a, c, e, f, g, x):
    return Dist(sqrt(x)/sqrt(e*x), Int((f + g*x)/(sqrt(x)*sqrt(a + c*x**S(2))), x), x)


def replacement470(a, b, c, d, e, f, g, m, p, x):
    return Dist(g/e, Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x), x) + Dist((-d*g + e*f)/e, Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x)


def replacement471(a, c, d, e, f, g, m, p, x):
    return Dist(g/e, Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x), x) + Dist((-d*g + e*f)/e, Int((a + c*x**S(2))**p*(d + e*x)**m, x), x)


def replacement472(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x)


def replacement473(a, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x)


def replacement474(a, b, c, d, e, f, g, p, x):
    return -Dist(S(1)/(e*(-d*g + e*f)), Int((a + b*x + c*x**S(2))**(p + S(-1))*Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(f + g*x), x), x) + Dist((a*e**S(2) - b*d*e + c*d**S(2))/(e*(-d*g + e*f)), Int((a + b*x + c*x**S(2))**(p + S(-1))/(d + e*x), x), x)


def replacement475(a, c, d, e, f, g, p, x):
    return -Dist(S(1)/(e*(-d*g + e*f)), Int((a + c*x**S(2))**(p + S(-1))*Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(f + g*x), x), x) + Dist((a*e**S(2) + c*d**S(2))/(e*(-d*g + e*f)), Int((a + c*x**S(2))**(p + S(-1))/(d + e*x), x), x)


def With476(a, b, c, d, e, f, g, m, n, p, x):
    q = Denominator(m)
    return Dist(q/e, Subst(Int(x**(q*(m + S(1)) + S(-1))*(g*x**q/e + (-d*g + e*f)/e)**n*(c*x**(S(2)*q)/e**S(2) - x**q*(-b*e + S(2)*c*d)/e**S(2) + (a*e**S(2) - b*d*e + c*d**S(2))/e**S(2))**p, x), x, (d + e*x)**(S(1)/q)), x)


def With477(a, c, d, e, f, g, m, n, p, x):
    q = Denominator(m)
    return Dist(q/e, Subst(Int(x**(q*(m + S(1)) + S(-1))*(g*x**q/e + (-d*g + e*f)/e)**n*(-S(2)*c*d*x**q/e**S(2) + c*x**(S(2)*q)/e**S(2) + (a*e**S(2) + c*d**S(2))/e**S(2))**p, x), x, (d + e*x)**(S(1)/q)), x)


def replacement478(a, b, c, d, e, f, g, m, n, p, x):
    return Int((d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x)


def replacement479(a, c, d, e, f, g, m, n, p, x):
    return Int((a + c*x**S(2))**p*(d*f + e*g*x**S(2))**m, x)


def replacement480(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m)), Int((d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x), x)


def replacement481(a, c, d, e, f, g, m, n, p, x):
    return Dist((d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m)), Int((a + c*x**S(2))**p*(d*f + e*g*x**S(2))**m, x), x)


def replacement482(a, b, c, d, e, f, g, m, n, x):
    return Dist(c, Int(x**S(2)*(d + e*x)**m*(f + g*x)**n, x), x) + Int((a + b*x)*(d + e*x)**m*(f + g*x)**n, x)


def replacement483(a, c, d, e, f, g, m, n, x):
    return Dist(a, Int((d + e*x)**m*(f + g*x)**n, x), x) + Dist(c, Int(x**S(2)*(d + e*x)**m*(f + g*x)**n, x), x)


def replacement484(a, b, c, d, e, f, g, m, n, x):
    return Dist(c**(S(-2)), Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(a*b*e*g**S(2) - a*c*d*g**S(2) - S(2)*a*c*e*f*g + c**S(2)*d*f**S(2) + x*(-a*c*e*g**S(2) + b**S(2)*e*g**S(2) - b*c*d*g**S(2) - S(2)*b*c*e*f*g + S(2)*c**S(2)*d*f*g + c**S(2)*e*f**S(2)), x)/(a + b*x + c*x**S(2)), x), x) + Dist(g/c**S(2), Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(-b*e*g + c*d*g + S(2)*c*e*f + c*e*g*x, x), x), x)


def replacement485(a, c, d, e, f, g, m, n, x):
    return Dist(S(1)/c, Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(-a*d*g**S(2) - S(2)*a*e*f*g + c*d*f**S(2) + x*(-a*e*g**S(2) + S(2)*c*d*f*g + c*e*f**S(2)), x)/(a + c*x**S(2)), x), x) + Dist(g/c, Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(d*g + S(2)*e*f + e*g*x, x), x), x)


def replacement486(a, b, c, d, e, f, g, m, n, x):
    return Dist(S(1)/c, Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*Simp(-a*e*g + c*d*f + x*(-b*e*g + c*d*g + c*e*f), x)/(a + b*x + c*x**S(2)), x), x) + Dist(e*g/c, Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x), x)


def replacement487(a, c, d, e, f, g, m, n, x):
    return Dist(S(1)/c, Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*Simp(-a*e*g + c*d*f + x*(c*d*g + c*e*f), x)/(a + c*x**S(2)), x), x) + Dist(e*g/c, Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x), x)


def replacement488(a, b, c, d, e, f, g, m, n, x):
    return -Dist(g*(-d*g + e*f)/(a*g**S(2) - b*f*g + c*f**S(2)), Int((d + e*x)**(m + S(-1))*(f + g*x)**n, x), x) + Dist(S(1)/(a*g**S(2) - b*f*g + c*f**S(2)), Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*Simp(a*e*g - b*d*g + c*d*f + c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x), x)


def replacement489(a, c, d, e, f, g, m, n, x):
    return -Dist(g*(-d*g + e*f)/(a*g**S(2) + c*f**S(2)), Int((d + e*x)**(m + S(-1))*(f + g*x)**n, x), x) + Dist(S(1)/(a*g**S(2) + c*f**S(2)), Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*Simp(a*e*g + c*d*f + c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x), x)


def replacement490(a, b, c, d, e, f, g, m, x):
    return Int(ExpandIntegrand(S(1)/(sqrt(d + e*x)*sqrt(f + g*x)), (d + e*x)**(m + S(1)/2)/(a + b*x + c*x**S(2)), x), x)


def replacement491(a, c, d, e, f, g, m, x):
    return Int(ExpandIntegrand(S(1)/(sqrt(d + e*x)*sqrt(f + g*x)), (d + e*x)**(m + S(1)/2)/(a + c*x**S(2)), x), x)


def replacement492(a, b, c, d, e, f, g, m, n, x):
    return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n, S(1)/(a + b*x + c*x**S(2)), x), x)


def replacement493(a, c, d, e, f, g, m, n, x):
    return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n, S(1)/(a + c*x**S(2)), x), x)


def replacement494(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x)


def replacement495(a, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x)


def replacement496(a, b, c, d, e, g, m, n, p, x):
    return Dist((d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((g*x)**n*(a*d + c*e*x**S(3))**p, x), x)


def replacement497(a, b, c, d, e, f, g, n, p, x):
    return -Dist(S(1)/(e*(-d*g + e*f)), Int((f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(-1))*(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f)), x), x) + Dist((a*e**S(2) - b*d*e + c*d**S(2))/(e*(-d*g + e*f)), Int((f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(-1))/(d + e*x), x), x)


def replacement498(a, c, d, e, f, g, n, p, x):
    return -Dist(S(1)/(e*(-d*g + e*f)), Int((a + c*x**S(2))**(p + S(-1))*(f + g*x)**n*(a*e*g + c*d*f - c*x*(-d*g + e*f)), x), x) + Dist((a*e**S(2) + c*d**S(2))/(e*(-d*g + e*f)), Int((a + c*x**S(2))**(p + S(-1))*(f + g*x)**(n + S(1))/(d + e*x), x), x)


def replacement499(a, b, c, d, e, f, g, n, p, x):
    return Dist(e*(-d*g + e*f)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(d + e*x), x), x) + Dist(S(1)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p*(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f)), x), x)


def replacement500(a, c, d, e, f, g, n, p, x):
    return Dist(e*(-d*g + e*f)/(a*e**S(2) + c*d**S(2)), Int((a + c*x**S(2))**(p + S(1))*(f + g*x)**(n + S(-1))/(d + e*x), x), x) + Dist(S(1)/(a*e**S(2) + c*d**S(2)), Int((a + c*x**S(2))**p*(f + g*x)**(n + S(-1))*(a*e*g + c*d*f - c*x*(-d*g + e*f)), x), x)


def With501(a, b, c, d, e, f, g, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Simp(-sqrt(S(2))*sqrt(-g*(b + S(2)*c*x - q)/(-b*g + S(2)*c*f + g*q))*sqrt(-g*(b + S(2)*c*x + q)/(-b*g + S(2)*c*f - g*q))*EllipticPi(e*(-b*g + S(2)*c*f + g*q)/(S(2)*c*(-d*g + e*f)), asin(sqrt(S(2))*sqrt(c/(-b*g + S(2)*c*f + g*q))*sqrt(f + g*x)), (-b*g + S(2)*c*f + g*q)/(-b*g + S(2)*c*f - g*q))/(sqrt(c/(-b*g + S(2)*c*f + g*q))*(-d*g + e*f)*sqrt(a + b*x + c*x**S(2))), x)


def With502(a, c, d, e, f, g, x):
    q = Rt(-a*c, S(2))
    return Simp(-S(2)*sqrt(g*(-c*x + q)/(c*f + g*q))*sqrt(-g*(c*x + q)/(c*f - g*q))*EllipticPi(e*(c*f + g*q)/(c*(-d*g + e*f)), asin(sqrt(c/(c*f + g*q))*sqrt(f + g*x)), (c*f + g*q)/(c*f - g*q))/(sqrt(c/(c*f + g*q))*sqrt(a + c*x**S(2))*(-d*g + e*f)), x)


def replacement503(a, b, c, d, e, f, g, n, x):
    return Int(ExpandIntegrand(S(1)/(sqrt(f + g*x)*sqrt(a + b*x + c*x**S(2))), (f + g*x)**(n + S(1)/2)/(d + e*x), x), x)


def replacement504(a, c, d, e, f, g, n, x):
    return Int(ExpandIntegrand(S(1)/(sqrt(a + c*x**S(2))*sqrt(f + g*x)), (f + g*x)**(n + S(1)/2)/(d + e*x), x), x)


def replacement505(a, b, c, d, e, f, g, x):
    return Dist(-S(2)*sqrt((-d*g + e*f)**S(2)*(a + b*x + c*x**S(2))/((d + e*x)**S(2)*(a*g**S(2) - b*f*g + c*f**S(2))))*(d + e*x)/((-d*g + e*f)*sqrt(a + b*x + c*x**S(2))), Subst(Int(S(1)/sqrt(x**S(4)*(a*e**S(2) - b*d*e + c*d**S(2))/(a*g**S(2) - b*f*g + c*f**S(2)) - x**S(2)*(S(2)*a*e*g - b*d*g - b*e*f + S(2)*c*d*f)/(a*g**S(2) - b*f*g + c*f**S(2)) + S(1)), x), x, sqrt(f + g*x)/sqrt(d + e*x)), x)


def replacement506(a, c, d, e, f, g, x):
    return Dist(-S(2)*sqrt((a + c*x**S(2))*(-d*g + e*f)**S(2)/((d + e*x)**S(2)*(a*g**S(2) + c*f**S(2))))*(d + e*x)/(sqrt(a + c*x**S(2))*(-d*g + e*f)), Subst(Int(S(1)/sqrt(x**S(4)*(a*e**S(2) + c*d**S(2))/(a*g**S(2) + c*f**S(2)) - x**S(2)*(S(2)*a*e*g + S(2)*c*d*f)/(a*g**S(2) + c*f**S(2)) + S(1)), x), x, sqrt(f + g*x)/sqrt(d + e*x)), x)


def replacement507(a, c, e, f, g, m, p, x):
    return Dist(S(2)*f*g/e, Int((e*x)**(m + S(1))*(a + c*x**S(2))**p, x), x) + Int((e*x)**m*(a + c*x**S(2))**p*(f**S(2) + g**S(2)*x**S(2)), x)


def replacement508(a, c, e, f, g, m, p, x):
    return Dist(f, Int((e*x)**m*(a + c*x**S(2))**p*(f**S(2) + S(3)*g**S(2)*x**S(2)), x), x) + Dist(g/e, Int((e*x)**(m + S(1))*(a + c*x**S(2))**p*(S(3)*f**S(2) + g**S(2)*x**S(2)), x), x)


def replacement509(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g/e, Int((d + e*x)**(m + S(1))*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x), x) + Dist((-d*g + e*f)/e, Int((d + e*x)**m*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x), x)


def replacement510(a, c, d, e, f, g, m, n, p, x):
    return Dist(g/e, Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(f + g*x)**(n + S(-1)), x), x) + Dist((-d*g + e*f)/e, Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(-1)), x), x)


def replacement511(a, b, c, d, e, f, g, m, n, p, x):
    return Int((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x)


def replacement512(a, c, d, e, f, g, m, n, p, x):
    return Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x)


def replacement513(a, b, c, d, e, f, g, m, n, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x, u), x)


def replacement514(a, c, d, e, f, g, m, n, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x, u), x)


def replacement515(a, b, c, d, e, f, p, q, x):
    return Dist((c/f)**p, Int((d + e*x + f*x**S(2))**(p + q), x), x)


def replacement516(a, b, c, d, e, f, p, q, x):
    return Dist(a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p)), Int((d + e*x + f*x**S(2))**(p + q), x), x)


def replacement517(a, b, c, d, e, f, p, q, x):
    return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((b + S(2)*c*x)**(S(2)*p)*(d + e*x + f*x**S(2))**q, x), x)


def replacement518(a, b, c, d, f, p, q, x):
    return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q, x), x)


def replacement519(a, b, c, d, e, f, q, x):
    return Simp((d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) - c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))), x)


def replacement520(a, c, d, e, f, q, x):
    return Simp((-c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))), x)


def replacement521(a, b, c, d, f, q, x):
    return Simp((d + f*x**S(2))**(q + S(1))*(S(2)*a*f*x*(q + S(1)) + b*d)/(S(2)*d*f*(q + S(1))), x)


def replacement522(a, b, c, d, f, q, x):
    return Dist(b, Int(x*(d + f*x**S(2))**q, x), x) + Int((a + c*x**S(2))*(d + f*x**S(2))**q, x)


def replacement523(a, b, c, d, e, f, q, x):
    return Int(ExpandIntegrand((a + b*x + c*x**S(2))*(d + e*x + f*x**S(2))**q, x), x)


def replacement524(a, c, d, e, f, q, x):
    return Int(ExpandIntegrand((a + c*x**S(2))*(d + e*x + f*x**S(2))**q, x), x)


def replacement525(a, b, c, d, e, f, q, x):
    return -Dist((c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))/(f*(q + S(1))*(-S(4)*d*f + e**S(2))), Int((d + e*x + f*x**S(2))**(q + S(1)), x), x) + Simp((d + e*x + f*x**S(2))**(q + S(1))*(a*e*f - S(2)*b*d*f + c*d*e + x*(c*(-S(2)*d*f + e**S(2)) + f*(S(2)*a*f - b*e)))/(f*(q + S(1))*(-S(4)*d*f + e**S(2))), x)


def replacement526(a, c, d, e, f, q, x):
    return -Dist((S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))/(f*(q + S(1))*(-S(4)*d*f + e**S(2))), Int((d + e*x + f*x**S(2))**(q + S(1)), x), x) + Simp((d + e*x + f*x**S(2))**(q + S(1))*(a*e*f + c*d*e + x*(S(2)*a*f**S(2) + c*(-S(2)*d*f + e**S(2))))/(f*(q + S(1))*(-S(4)*d*f + e**S(2))), x)


def replacement527(a, b, c, d, f, q, x):
    return Dist((S(2)*a*f*q + S(3)*a*f - c*d)/(S(2)*d*f*(q + S(1))), Int((d + f*x**S(2))**(q + S(1)), x), x) + Simp((d + f*x**S(2))**(q + S(1))*(b*d - x*(a*f - c*d))/(S(2)*d*f*(q + S(1))), x)


def replacement528(a, b, c, d, e, f, q, x):
    return Dist((c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))/(S(2)*f**S(2)*(S(2)*q + S(3))), Int((d + e*x + f*x**S(2))**q, x), x) + Simp((d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) - c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))), x)


def replacement529(a, c, d, e, f, q, x):
    return Dist((S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))/(S(2)*f**S(2)*(S(2)*q + S(3))), Int((d + e*x + f*x**S(2))**q, x), x) + Simp((-c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))), x)


def replacement530(a, b, c, d, f, q, x):
    return Dist((S(2)*a*f*q + S(3)*a*f - c*d)/(f*(S(2)*q + S(3))), Int((d + f*x**S(2))**q, x), x) + Simp((d + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))), x)


def replacement531(a, b, c, d, e, f, p, q, x):
    return -Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(b*e*q + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)) + x*(S(2)*b*f*q + S(2)*c*e*(S(2)*p + q + S(3))), x), x), x) + Simp((b + S(2)*c*x)*(a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement532(a, b, c, d, f, p, q, x):
    return -Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(S(2)*b*f*q*x + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)), x), x), x) + Simp((b + S(2)*c*x)*(d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement533(a, c, d, e, f, p, q, x):
    return -Dist(-S(1)/(S(4)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*e*x*(S(2)*p + q + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)), x), x), x) + Simp(-x*(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q/(S(2)*a*(p + S(1))), x)


def replacement534(a, b, c, d, e, f, p, q, x):
    return -Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))), Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + S(2)*c*(p + S(1))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) - e*(p + q + S(2))*(-S(2)*a*c**S(2)*e - b**S(3)*f + b**S(2)*c*e - b*c*(-S(3)*a*f + c*d)) + x*(S(2)*f*(p + q + S(2))*(S(2)*a*c**S(2)*e + b**S(3)*f - b**S(2)*c*e + b*c*(-S(3)*a*f + c*d)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e))) - (a*f*(p + S(1)) - c*d*(p + S(2)))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)), x), x), x) + Simp((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(S(2)*a*c**S(2)*e + b**S(3)*f - b**S(2)*c*e + b*c*(-S(3)*a*f + c*d) + c*x*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))), x)


def replacement535(a, b, c, d, f, p, q, x):
    return -Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))), Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + S(2)*c*(p + S(1))*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + x*(-b*f*(p + S(1))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + S(2)*f*(b**S(3)*f + b*c*(-S(3)*a*f + c*d))*(p + q + S(2))) - (a*f*(p + S(1)) - c*d*(p + S(2)))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d), x), x), x) + Simp((d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(b**S(3)*f + b*c*(-S(3)*a*f + c*d) + c*x*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))), x)


def replacement536(a, c, d, e, f, p, q, x):
    return -Dist(-S(1)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(S(2)*a*c**S(2)*e**S(2)*(p + q + S(2)) + c*f*x**S(2)*(-S(2)*a*c*f + S(2)*c**S(2)*d)*(S(2)*p + S(2)*q + S(5)) + S(2)*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + x*(S(4)*a*c**S(2)*e*f*(p + q + S(2)) + c*e*(-S(2)*a*c*f + S(2)*c**S(2)*d)*(S(2)*p + q + S(4))) - (-S(2)*a*c*f + S(2)*c**S(2)*d)*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x), x) + Simp(-(a + c*x**S(2))**(p + S(1))*(S(2)*a*c**S(2)*e + c*x*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(d + e*x + f*x**S(2))**(q + S(1))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))), x)


def replacement537(a, b, c, d, e, f, p, q, x):
    return -Dist(S(1)/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))), Int((a + b*x + c*x**S(2))**(p + S(-2))*(d + e*x + f*x**S(2))**q*Simp(x**S(2)*(c*(p + q)*(-c*(S(2)*d*f*(S(1) - S(2)*p) + e**S(2)*(S(3)*p + q + S(-1))) + f*(-S(2)*a*f + b*e)*(S(4)*p + S(2)*q + S(-1))) + p*(S(1) - p)*(-b*f + c*e)**S(2)) + x*(S(2)*(S(1) - p)*(S(2)*p + q)*(-a*f + c*d)*(-b*f + c*e) - (p + q)*(b*(c*(S(2)*p + q)*(-S(4)*d*f + e**S(2)) + f*(S(2)*p + S(2)*q + S(1))*(S(2)*a*f - b*e + S(2)*c*d)) + e*f*(S(1) - p)*(-S(4)*a*c + b**S(2)))) + (S(1) - p)*(S(2)*p + q)*(-a*e + b*d)*(-b*f + c*e) - (p + q)*(-a*(c*(S(2)*d*f - e**S(2)*(S(2)*p + q)) + f*(-S(2)*a*f + b*e)*(S(2)*p + S(2)*q + S(1))) + b**S(2)*d*f*(S(1) - p)), x), x), x) + Simp((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(3)*p + S(2)*q) - c*e*(S(2)*p + q) + S(2)*c*f*x*(p + q))/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))), x)


def replacement538(a, b, c, d, f, p, q, x):
    return -Dist(S(1)/(S(2)*f*(p + q)*(S(2)*p + S(2)*q + S(1))), Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-2))*Simp(b**S(2)*d*(p + S(-1))*(S(2)*p + q) + x**S(2)*(b**S(2)*f*p*(S(1) - p) + S(2)*c*(p + q)*(-a*f*(S(4)*p + S(2)*q + S(-1)) + c*d*(S(2)*p + S(-1)))) - x*(S(2)*b*(S(1) - p)*(S(2)*p + q)*(-a*f + c*d) - S(2)*b*(p + q)*(S(2)*c*d*(S(2)*p + q) - (a*f + c*d)*(S(2)*p + S(2)*q + S(1)))) - (p + q)*(-S(2)*a*(-a*f*(S(2)*p + S(2)*q + S(1)) + c*d) + b**S(2)*d*(S(1) - p)), x), x), x) + Simp((d + f*x**S(2))**(q + S(1))*(b*(S(3)*p + S(2)*q) + S(2)*c*x*(p + q))*(a + b*x + c*x**S(2))**(p + S(-1))/(S(2)*f*(p + q)*(S(2)*p + S(2)*q + S(1))), x)


def replacement539(a, c, d, e, f, p, q, x):
    return -Dist(S(1)/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))), Int((a + c*x**S(2))**(p + S(-2))*(d + e*x + f*x**S(2))**q*Simp(-a*c*e**S(2)*(S(1) - p)*(S(2)*p + q) + a*(p + q)*(-S(2)*a*f**S(2)*(S(2)*p + S(2)*q + S(1)) + c*(S(2)*d*f - e**S(2)*(S(2)*p + q))) + x**S(2)*(c**S(2)*e**S(2)*p*(S(1) - p) - c*(p + q)*(S(2)*a*f**S(2)*(S(4)*p + S(2)*q + S(-1)) + c*(S(2)*d*f*(S(1) - S(2)*p) + e**S(2)*(S(3)*p + q + S(-1))))) + x*(S(4)*a*c*e*f*(S(1) - p)*(p + q) + S(2)*c*e*(S(1) - p)*(S(2)*p + q)*(-a*f + c*d)), x), x), x) - Simp(c*(a + c*x**S(2))**(p + S(-1))*(e*(S(2)*p + q) - S(2)*f*x*(p + q))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))), x)


def With540(a, b, c, d, e, f, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
    if NonzeroQ(q):
        return True
    return False


def replacement540(a, b, c, d, e, f, x):

    q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
    return Dist(S(1)/q, Int((-a*c*f + b**S(2)*f - b*c*e + c**S(2)*d - x*(-b*c*f + c**S(2)*e))/(a + b*x + c*x**S(2)), x), x) + Dist(S(1)/q, Int((a*f**S(2) - b*e*f - c*d*f + c*e**S(2) + x*(-b*f**S(2) + c*e*f))/(d + e*x + f*x**S(2)), x), x)


def With541(a, b, c, d, f, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
    if NonzeroQ(q):
        return True
    return False


def replacement541(a, b, c, d, f, x):

    q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
    return -Dist(S(1)/q, Int((-a*f**S(2) + b*f**S(2)*x + c*d*f)/(d + f*x**S(2)), x), x) + Dist(S(1)/q, Int((-a*c*f + b**S(2)*f + b*c*f*x + c**S(2)*d)/(a + b*x + c*x**S(2)), x), x)


def replacement542(a, b, c, d, e, f, x):
    return Dist(-S(2)*e, Subst(Int(S(1)/(e*(-S(4)*a*f + b*e) - x**S(2)*(-a*e + b*d)), x), x, (e + S(2)*f*x)/sqrt(d + e*x + f*x**S(2))), x)


def With543(a, b, c, d, e, f, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/q, Int(S(1)/((b + S(2)*c*x - q)*sqrt(d + e*x + f*x**S(2))), x), x) - Dist(S(2)*c/q, Int(S(1)/((b + S(2)*c*x + q)*sqrt(d + e*x + f*x**S(2))), x), x)


def replacement544(a, c, d, e, f, x):
    return Dist(S(1)/2, Int(S(1)/((a - x*Rt(-a*c, S(2)))*sqrt(d + e*x + f*x**S(2))), x), x) + Dist(S(1)/2, Int(S(1)/((a + x*Rt(-a*c, S(2)))*sqrt(d + e*x + f*x**S(2))), x), x)


def With545(a, b, c, d, f, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/q, Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x - q)), x), x) - Dist(S(2)*c/q, Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x + q)), x), x)


def With546(a, b, c, d, e, f, x):
    q = Rt(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2), S(2))
    return -Dist(S(1)/(S(2)*q), Int((-a*f + c*d - q + x*(-b*f + c*e))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x) + Dist(S(1)/(S(2)*q), Int((-a*f + c*d + q + x*(-b*f + c*e))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x)


def With547(a, c, d, e, f, x):
    q = Rt(a*c*e**S(2) + (-a*f + c*d)**S(2), S(2))
    return -Dist(S(1)/(S(2)*q), Int((-a*f + c*d + c*e*x - q)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x) + Dist(S(1)/(S(2)*q), Int((-a*f + c*d + c*e*x + q)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x)


def With548(a, b, c, d, f, x):
    q = Rt(b**S(2)*d*f + (-a*f + c*d)**S(2), S(2))
    return -Dist(S(1)/(S(2)*q), Int((-a*f - b*f*x + c*d - q)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x), x) + Dist(S(1)/(S(2)*q), Int((-a*f - b*f*x + c*d + q)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x), x)


def replacement549(a, b, c, d, e, f, x):
    return -Dist(S(1)/f, Int((-a*f + c*d + x*(-b*f + c*e))/(sqrt(a + b*x + c*x**S(2))*(d + e*x + f*x**S(2))), x), x) + Dist(c/f, Int(S(1)/sqrt(a + b*x + c*x**S(2)), x), x)


def replacement550(a, b, c, d, f, x):
    return -Dist(S(1)/f, Int((-a*f - b*f*x + c*d)/((d + f*x**S(2))*sqrt(a + b*x + c*x**S(2))), x), x) + Dist(c/f, Int(S(1)/sqrt(a + b*x + c*x**S(2)), x), x)


def replacement551(a, c, d, e, f, x):
    return -Dist(S(1)/f, Int((-a*f + c*d + c*e*x)/(sqrt(a + c*x**S(2))*(d + e*x + f*x**S(2))), x), x) + Dist(c/f, Int(S(1)/sqrt(a + c*x**S(2)), x), x)


def With552(a, b, c, d, e, f, x):
    r = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)/sqrt(a + b*x + c*x**S(2)), Int(S(1)/(sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*sqrt(d + e*x + f*x**S(2))), x), x)


def With553(a, b, c, d, f, x):
    r = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)/sqrt(a + b*x + c*x**S(2)), Int(S(1)/(sqrt(S(2)*a + x*(b + r))*sqrt(d + f*x**S(2))*sqrt(b + S(2)*c*x + r)), x), x)


def replacement554(a, b, c, d, e, f, p, q, x):
    return Int((a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)


def replacement555(a, c, d, e, f, p, q, x):
    return Int((a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)


def replacement556(a, b, c, d, e, f, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u), x)


def replacement557(a, c, d, e, f, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u), x)


def replacement558(a, b, c, d, e, f, g, h, m, p, q, x):
    return Dist((c/f)**p, Int((g + h*x)**m*(d + e*x + f*x**S(2))**(p + q), x), x)


def replacement559(a, b, c, d, e, f, g, h, m, p, q, x):
    return Dist(a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p)), Int((g + h*x)**m*(d + e*x + f*x**S(2))**(p + q), x), x)


def replacement560(a, b, c, d, e, f, g, h, m, p, q, x):
    return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((b + S(2)*c*x)**(S(2)*p)*(g + h*x)**m*(d + e*x + f*x**S(2))**q, x), x)


def replacement561(a, b, c, d, f, g, h, m, p, q, x):
    return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q*(g + h*x)**m, x), x)


def replacement562(a, b, c, d, e, f, g, h, m, p, x):
    return Int((f*h*x/c + d*g/a)**m*(a + b*x + c*x**S(2))**(m + p), x)


def replacement563(a, c, d, e, f, g, h, m, p, x):
    return Int((a + c*x**S(2))**(m + p)*(f*h*x/c + d*g/a)**m, x)


def replacement564(a, b, c, d, f, g, h, m, p, x):
    return Int((f*h*x/c + d*g/a)**m*(a + b*x + c*x**S(2))**(m + p), x)


def replacement565(a, c, d, f, g, h, m, p, x):
    return Int((a + c*x**S(2))**(m + p)*(f*h*x/c + d*g/a)**m, x)


def replacement566(a, b, c, e, f, p, q, x):
    return Int((a/e + c*x/f)**p*(e*x + f*x**S(2))**(p + q), x)


def replacement567(a, c, e, f, p, q, x):
    return Int((a/e + c*x/f)**p*(e*x + f*x**S(2))**(p + q), x)


def replacement568(a, b, c, d, e, f, g, h, m, p, x):
    return Simp(f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))), x)


def replacement569(a, c, d, e, f, g, h, m, p, x):
    return Simp(f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))), x)


def replacement570(a, b, c, d, f, g, h, m, p, x):
    return Simp(f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))), x)


def replacement571(a, b, c, d, e, f, g, h, m, p, x):
    return Int(ExpandIntegrand((g + h*x)**m*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2)), x), x)


def replacement572(a, c, d, e, f, g, h, m, p, x):
    return Int(ExpandIntegrand((a + c*x**S(2))**p*(g + h*x)**m*(d + e*x + f*x**S(2)), x), x)


def replacement573(a, b, c, d, f, g, h, m, p, x):
    return Int(ExpandIntegrand((d + f*x**S(2))*(g + h*x)**m*(a + b*x + c*x**S(2))**p, x), x)


def replacement574(a, c, d, f, g, h, m, p, x):
    return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + f*x**S(2))*(g + h*x)**m, x), x)


def replacement575(a, b, c, d, e, f, g, h, m, p, x):
    return Dist(S(1)/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))), Int((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*(f*g**S(2)*(p + S(1)) - h*(-d*h*(m + p + S(2)) + e*g*(p + S(1)))) + h*(m + S(1))*(a*e*h - a*f*g + c*d*g) - x*(c*(S(2)*f*g**S(2)*(p + S(1)) - h*(-d*h + e*g)*(m + S(2)*p + S(3))) + f*h*(m + S(1))*(-a*h + b*g)), x), x), x) + Simp((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(d*h**S(2) - e*g*h + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))), x)


def replacement576(a, c, d, e, f, g, h, m, p, x):
    return Dist(S(1)/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))), Int((a + c*x**S(2))**p*(g + h*x)**(m + S(1))*Simp(h*(m + S(1))*(a*e*h - a*f*g + c*d*g) + x*(a*f*h**S(2)*(m + S(1)) - c*(S(2)*f*g**S(2)*(p + S(1)) - h*(-d*h + e*g)*(m + S(2)*p + S(3)))), x), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(d*h**S(2) - e*g*h + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))), x)


def replacement577(a, b, c, d, f, g, h, m, p, x):
    return Dist(S(1)/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))), Int((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*(d*h**S(2)*(m + p + S(2)) + f*g**S(2)*(p + S(1))) + h*(m + S(1))*(-a*f*g + c*d*g) - x*(c*(d*h**S(2)*(m + S(2)*p + S(3)) + S(2)*f*g**S(2)*(p + S(1))) + f*h*(m + S(1))*(-a*h + b*g)), x), x), x) + Simp((g + h*x)**(m + S(1))*(d*h**S(2) + f*g**S(2))*(a + b*x + c*x**S(2))**(p + S(1))/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))), x)


def replacement578(a, c, d, f, g, h, m, p, x):
    return Dist(S(1)/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))), Int((a + c*x**S(2))**p*(g + h*x)**(m + S(1))*Simp(h*(m + S(1))*(-a*f*g + c*d*g) + x*(a*f*h**S(2)*(m + S(1)) - c*(d*h**S(2)*(m + S(2)*p + S(3)) + S(2)*f*g**S(2)*(p + S(1)))), x), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(d*h**S(2) + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))), x)


def replacement579(a, b, c, d, e, f, g, h, x):
    return Dist((d*h**S(2) - e*g*h + f*g**S(2))/(a*h**S(2) - b*g*h + c*g**S(2)), Int(S(1)/((g + h*x)*sqrt(a + b*x + c*x**S(2))), x), x) + Dist(S(1)/(a*h**S(2) - b*g*h + c*g**S(2)), Int((a*e*h - a*f*g - b*d*h + c*d*g + x*(a*f*h - b*f*g - c*d*h + c*e*g))/(a + b*x + c*x**S(2))**(S(3)/2), x), x)


def replacement580(a, c, d, e, f, g, h, x):
    return Dist((d*h**S(2) - e*g*h + f*g**S(2))/(a*h**S(2) + c*g**S(2)), Int(S(1)/(sqrt(a + c*x**S(2))*(g + h*x)), x), x) + Dist(S(1)/(a*h**S(2) + c*g**S(2)), Int((a*e*h - a*f*g + c*d*g + x*(a*f*h - c*d*h + c*e*g))/(a + c*x**S(2))**(S(3)/2), x), x)


def replacement581(a, b, c, d, f, g, h, x):
    return Dist((d*h**S(2) + f*g**S(2))/(a*h**S(2) - b*g*h + c*g**S(2)), Int(S(1)/((g + h*x)*sqrt(a + b*x + c*x**S(2))), x), x) + Dist(S(1)/(a*h**S(2) - b*g*h + c*g**S(2)), Int((-a*f*g - b*d*h + c*d*g - x*(-a*f*h + b*f*g + c*d*h))/(a + b*x + c*x**S(2))**(S(3)/2), x), x)


def replacement582(a, c, d, f, g, h, x):
    return Dist((d*h**S(2) + f*g**S(2))/(a*h**S(2) + c*g**S(2)), Int(S(1)/(sqrt(a + c*x**S(2))*(g + h*x)), x), x) + Dist(S(1)/(a*h**S(2) + c*g**S(2)), Int((-a*f*g + c*d*g - x*(-a*f*h + c*d*h))/(a + c*x**S(2))**(S(3)/2), x), x)


def replacement583(a, b, c, d, e, f, g, h, m, p, x):
    return -Dist(S(1)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((g + h*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(c*(S(2)*p + S(3))*(-b*e + S(2)*c*d) - f*(S(2)*a*c - b**S(2)*(p + S(2)))) + h*m*(a*b*f - S(2)*a*c*e + b*c*d) + h*x*(c*(-b*e + S(2)*c*d)*(m + S(2)*p + S(3)) - f*(S(2)*a*c*(m + S(1)) - b**S(2)*(m + p + S(2)))), x), x), x) + Simp((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(a*b*f - S(2)*a*c*e + b*c*d + x*(c*(-b*e + S(2)*c*d) + f*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement584(a, c, d, e, f, g, h, m, p, x):
    return -Dist(S(1)/(S(2)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(-1))*Simp(a*(e*h*m + f*g) - c*d*g*(S(2)*p + S(3)) + h*x*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(g + h*x)**m*(a*e - x*(-a*f + c*d))/(S(2)*a*c*(p + S(1))), x)


def replacement585(a, b, c, d, f, g, h, m, p, x):
    return -Dist(S(1)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((g + h*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(S(2)*c**S(2)*d*(S(2)*p + S(3)) - f*(S(2)*a*c - b**S(2)*(p + S(2)))) + h*m*(a*b*f + b*c*d) + h*x*(S(2)*c**S(2)*d*(m + S(2)*p + S(3)) - f*(S(2)*a*c*(m + S(1)) - b**S(2)*(m + p + S(2)))), x), x), x) + Simp((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(a*b*f + b*c*d + x*(S(2)*c**S(2)*d + f*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement586(a, c, d, f, g, h, m, p, x):
    return -Dist(S(1)/(S(2)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(-1))*Simp(a*f*g - c*d*g*(S(2)*p + S(3)) + h*x*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x), x) - Simp(x*(a + c*x**S(2))**(p + S(1))*(g + h*x)**m*(-a*f + c*d)/(S(2)*a*c*(p + S(1))), x)


def replacement587(a, b, c, d, e, f, g, h, m, p, x):
    return -Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))), Int((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(p + S(2))*(-S(2)*a*(-c*e*h + c*f*g) + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c**S(2)*d*g) + h*x*(m + S(2)*p + S(4))*(-S(2)*a*(-c*e*h + c*f*g) + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c**S(2)*d*g) - h*(-(-a*e + b*d)*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))*(m + p + S(2)) + (p + S(1))*(c*g**S(2) - h*(-a*h + b*g))*(S(2)*a*f - b*e + S(2)*c*d), x), x), x) - Simp((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-x*(b*f*(-a*h + b*g) + S(2)*c**S(2)*d*g - c*(-S(2)*a*e*h + S(2)*a*f*g + b*d*h + b*e*g)) - (-a*e + b*d)*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))), x)


def replacement588(a, c, d, e, f, g, h, m, p, x):
    return Dist(S(1)/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))), Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**m*Simp(g*(p + S(2))*(-a*(-c*e*h + c*f*g) + c**S(2)*d*g) + h*x*(-a*(-c*e*h + c*f*g) + c**S(2)*d*g)*(m + S(2)*p + S(4)) - h*(a*c*e*g - a*h*(-a*f + c*d))*(m + p + S(2)) + (p + S(1))*(a*f + c*d)*(a*h**S(2) + c*g**S(2)), x), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(a*c*e*g - a*h*(-a*f + c*d) - c*x*(a*e*h - a*f*g + c*d*g))/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))), x)


def replacement589(a, b, c, d, f, g, h, m, p, x):
    return -Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))), Int((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(p + S(2))*(-S(2)*a*c*f*g + b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c**S(2)*d*g) + h*x*(m + S(2)*p + S(4))*(-S(2)*a*c*f*g + b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c**S(2)*d*g) - h*(-b*d*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))*(m + p + S(2)) + S(2)*(p + S(1))*(a*f + c*d)*(c*g**S(2) - h*(-a*h + b*g)), x), x), x) - Simp((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*d*(-b*h + S(2)*c*g) - x*(b*f*(-a*h + b*g) + S(2)*c**S(2)*d*g - c*(S(2)*a*f*g + b*d*h)) + (-a*f + c*d)*(-S(2)*a*h + b*g))/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))), x)


def replacement590(a, c, d, f, g, h, m, p, x):
    return Dist(S(1)/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))), Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**m*Simp(a*h**S(2)*(-a*f + c*d)*(m + p + S(2)) + g*(p + S(2))*(-a*c*f*g + c**S(2)*d*g) + h*x*(-a*c*f*g + c**S(2)*d*g)*(m + S(2)*p + S(4)) + (p + S(1))*(a*f + c*d)*(a*h**S(2) + c*g**S(2)), x), x), x) - Simp((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(a*h*(-a*f + c*d) + c*x*(-a*f*g + c*d*g))/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))), x)


def replacement591(a, b, c, d, e, f, g, h, m, p, x):
    return -Dist(h**(S(-2)), Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(-d*h**S(2) + f*g**S(2) + h*x*(-e*h + S(2)*f*g)), x), x) + Dist(f/h**S(2), Int((g + h*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x), x)


def replacement592(a, c, d, e, f, g, h, m, p, x):
    return -Dist(h**(S(-2)), Int((a + c*x**S(2))**p*(g + h*x)**m*(-d*h**S(2) + f*g**S(2) + h*x*(-e*h + S(2)*f*g)), x), x) + Dist(f/h**S(2), Int((a + c*x**S(2))**p*(g + h*x)**(m + S(2)), x), x)


def replacement593(a, b, c, d, f, g, h, m, p, x):
    return -Dist(h**(S(-2)), Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(-d*h**S(2) + f*g**S(2) + S(2)*f*g*h*x), x), x) + Dist(f/h**S(2), Int((g + h*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x), x)


def replacement594(a, c, d, f, g, h, m, p, x):
    return -Dist(h**(S(-2)), Int((a + c*x**S(2))**p*(g + h*x)**m*(-d*h**S(2) + f*g**S(2) + S(2)*f*g*h*x), x), x) + Dist(f/h**S(2), Int((a + c*x**S(2))**p*(g + h*x)**(m + S(2)), x), x)


def replacement595(a, b, c, d, e, f, g, h, m, p, x):
    return -Dist(S(1)/(c*h*(m + S(2)*p + S(3))), Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*f*g*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))) + x*(b*f*h*(m + p + S(2)) + c*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1)))), x), x), x) + Simp(f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))), x)


def replacement596(a, c, d, e, f, g, h, m, p, x):
    return -Dist(S(1)/(c*h*(m + S(2)*p + S(3))), Int((a + c*x**S(2))**p*(g + h*x)**m*Simp(c*x*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1))) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x), x) + Simp(f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))), x)


def replacement597(a, b, c, d, f, g, h, m, p, x):
    return -Dist(S(1)/(c*h*(m + S(2)*p + S(3))), Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*f*g*(p + S(1)) + f*x*(b*h*(m + p + S(2)) + S(2)*c*g*(p + S(1))) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x), x) + Simp(f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))), x)


def replacement598(a, c, d, f, g, h, m, p, x):
    return -Dist(S(1)/(c*h*(m + S(2)*p + S(3))), Int((a + c*x**S(2))**p*(g + h*x)**m*Simp(S(2)*c*f*g*x*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x), x) + Simp(f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))), x)


def replacement599(a, b, c, d, e, f, g, h, p, q, x):
    return Int(ExpandIntegrand((g + h*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x)


def replacement600(a, c, d, e, f, g, h, p, q, x):
    return Int(ExpandIntegrand((a + c*x**S(2))**p*(g + h*x)*(d + e*x + f*x**S(2))**q, x), x)


def replacement601(a, b, c, d, e, f, g, h, p, q, x):
    return -Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-d*(S(2)*p + S(3))*(b*h - S(2)*c*g) + e*q*(-S(2)*a*h + b*g) - f*x**S(2)*(b*h - S(2)*c*g)*(S(2)*p + S(2)*q + S(3)) + x*(-e*(b*h - S(2)*c*g)*(S(2)*p + q + S(3)) + S(2)*f*q*(-S(2)*a*h + b*g)), x), x), x) + Simp((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g))/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement602(a, c, d, e, f, g, h, p, q, x):
    return Dist(S(1)/(S(2)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-a*e*h*q + c*d*g*(S(2)*p + S(3)) + c*f*g*x**S(2)*(S(2)*p + S(2)*q + S(3)) + x*(-S(2)*a*f*h*q + c*e*g*(S(2)*p + q + S(3))), x), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(a*h - c*g*x)*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))), x)


def replacement603(a, b, c, d, f, g, h, p, q, x):
    return -Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-d*(S(2)*p + S(3))*(b*h - S(2)*c*g) + S(2)*f*q*x*(-S(2)*a*h + b*g) - f*x**S(2)*(b*h - S(2)*c*g)*(S(2)*p + S(2)*q + S(3)), x), x), x) + Simp((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g))/((p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement604(a, b, c, d, e, f, g, h, p, q, x):
    return Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))), Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d)) - e*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))*(p + q + S(2)) - x*(S(2)*f*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d))) + (p + S(1))*(b*h - S(2)*c*g)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d)), x), x), x) + Simp((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + c*x*(g*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) - h*(a*b*f - S(2)*a*c*e + b*c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))), x)


def replacement605(a, c, d, e, f, g, h, p, q, x):
    return Dist(-S(1)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d))*(S(2)*p + S(2)*q + S(5)) - S(2)*c*g*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)) - e*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(p + q + S(2)) - x*(c*e*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d))*(S(2)*p + q + S(4)) + S(2)*f*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d)), x), x), x) + Simp(-(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d) + c*x*(S(2)*a*c*e*h + g*(-S(2)*a*c*f + S(2)*c**S(2)*d)))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))), x)


def replacement606(a, b, c, d, f, g, h, p, q, x):
    return Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))), Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)) - x*(-b*f*(p + S(1))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)) + S(2)*f*(-b*c*g*(a*f + c*d) + (-a*h + b*g)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (p + S(1))*(b*h - S(2)*c*g)*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)), x), x), x) + Simp((d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*c*g*(a*f + c*d) + c*x*(g*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - h*(a*b*f + b*c*d)) + (-a*h + b*g)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))), x)


def replacement607(a, b, c, d, e, f, g, h, p, q, x):
    return -Dist(S(1)/(S(2)*f*(p + q + S(1))), Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*(e*h - S(2)*f*g)*(p + q + S(1)) + h*p*(-a*e + b*d) + x**S(2)*(c*(e*h - S(2)*f*g)*(p + q + S(1)) + h*p*(-b*f + c*e)) + x*(b*(e*h - S(2)*f*g)*(p + q + S(1)) + S(2)*h*p*(-a*f + c*d)), x), x), x) + Simp(h*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f*(p + q + S(1))), x)


def replacement608(a, c, d, e, f, g, h, p, q, x):
    return Dist(S(1)/(S(2)*f*(p + q + S(1))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*e*h*p - a*(e*h - S(2)*f*g)*(p + q + S(1)) - S(2)*h*p*x*(-a*f + c*d) - x**S(2)*(c*e*h*p + c*(e*h - S(2)*f*g)*(p + q + S(1))), x), x), x) + Simp(h*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f*(p + q + S(1))), x)


def replacement609(a, b, c, d, f, g, h, p, q, x):
    return -Dist(S(1)/(S(2)*f*(p + q + S(1))), Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-S(2)*a*f*g*(p + q + S(1)) + b*d*h*p + x**S(2)*(-b*f*h*p - S(2)*c*f*g*(p + q + S(1))) + x*(-S(2)*b*f*g*(p + q + S(1)) + S(2)*h*p*(-a*f + c*d)), x), x), x) + Simp(h*(d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**p/(S(2)*f*(p + q + S(1))), x)


def With610(a, b, c, d, e, f, g, h, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
    if NonzeroQ(q):
        return True
    return False


def replacement610(a, b, c, d, e, f, g, h, x):

    q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
    return Dist(S(1)/q, Int(Simp(-a*b*f*h + a*c*e*h - a*c*f*g + b**S(2)*f*g - b*c*e*g + c**S(2)*d*g + c*x*(-a*f*h + b*f*g + c*d*h - c*e*g), x)/(a + b*x + c*x**S(2)), x), x) + Dist(S(1)/q, Int(Simp(a*f**S(2)*g + b*d*f*h - b*e*f*g - c*d*e*h - c*d*f*g + c*e**S(2)*g - f*x*(-a*f*h + b*f*g + c*d*h - c*e*g), x)/(d + e*x + f*x**S(2)), x), x)


def With611(a, b, c, d, f, g, h, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
    if NonzeroQ(q):
        return True
    return False


def replacement611(a, b, c, d, f, g, h, x):

    q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
    return Dist(S(1)/q, Int(Simp(a*f**S(2)*g + b*d*f*h - c*d*f*g - f*x*(-a*f*h + b*f*g + c*d*h), x)/(d + f*x**S(2)), x), x) + Dist(S(1)/q, Int(Simp(-a*b*f*h - a*c*f*g + b**S(2)*f*g + c**S(2)*d*g + c*x*(-a*f*h + b*f*g + c*d*h), x)/(a + b*x + c*x**S(2)), x), x)


def replacement612(a, c, d, f, g, h, x):
    return Dist(g, Int(S(1)/((a + c*x**S(2))*sqrt(d + f*x**S(2))), x), x) + Dist(h, Int(x/((a + c*x**S(2))*sqrt(d + f*x**S(2))), x), x)


def With613(a, c, d, f, g, h, x):
    q = Rt(-a*c, S(2))
    return -Dist((c*g - h*q)/(S(2)*q), Int(S(1)/(sqrt(d + f*x**S(2))*(c*x + q)), x), x) - Dist((c*g + h*q)/(S(2)*q), Int(S(1)/(sqrt(d + f*x**S(2))*(-c*x + q)), x), x)


def replacement614(a, b, c, d, e, f, g, h, x):
    return Dist(-S(2)*g, Subst(Int(S(1)/(-a*e + b*d - b*x**S(2)), x), x, sqrt(d + e*x + f*x**S(2))), x)


def replacement615(a, b, c, d, e, f, g, h, x):
    return Dist(h/(S(2)*f), Int((e + S(2)*f*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x) - Dist((e*h - S(2)*f*g)/(S(2)*f), Int(S(1)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x)


def replacement616(a, b, c, d, e, f, x):
    return Dist(-S(2)*e, Subst(Int((-d*x**S(2) + S(1))/(-b*f + c*e + d**S(2)*x**S(4)*(-b*f + c*e) - e*x**S(2)*(S(2)*a*f - b*e + S(2)*c*d)), x), x, (S(1) + x*(e + sqrt(-S(4)*d*f + e**S(2)))/(S(2)*d))/sqrt(d + e*x + f*x**S(2))), x)


def replacement617(a, b, c, d, e, f, g, h, x):
    return Dist(g, Subst(Int(S(1)/(a + x**S(2)*(-a*f + c*d)), x), x, x/sqrt(d + e*x + f*x**S(2))), x)


def replacement618(a, b, c, d, e, f, g, h, x):
    return Dist(h/e, Int((S(2)*d + e*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x) - Dist((S(2)*d*h - e*g)/e, Int(S(1)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x)


def replacement619(a, b, c, d, e, f, g, h, x):
    return Dist(-S(2)*g*(-S(2)*a*h + b*g), Subst(Int(S(1)/Simp(g*(-S(4)*a*c + b**S(2))*(-S(2)*a*h + b*g) - x**S(2)*(-a*e + b*d), x), x), x, Simp(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g), x)/sqrt(d + e*x + f*x**S(2))), x)


def replacement620(a, c, d, e, f, g, h, x):
    return Dist(-S(2)*a*g*h, Subst(Int(S(1)/Simp(S(2)*a**S(2)*c*g*h + a*e*x**S(2), x), x), x, Simp(a*h - c*g*x, x)/sqrt(d + e*x + f*x**S(2))), x)


def replacement621(a, b, c, d, f, g, h, x):
    return Dist(-S(2)*g*(-S(2)*a*h + b*g), Subst(Int(S(1)/Simp(-b*d*x**S(2) + g*(-S(4)*a*c + b**S(2))*(-S(2)*a*h + b*g), x), x), x, Simp(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g), x)/sqrt(d + f*x**S(2))), x)


def With622(a, b, c, d, e, f, g, h, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist((S(2)*c*g - h*(b - q))/q, Int(S(1)/((b + S(2)*c*x - q)*sqrt(d + e*x + f*x**S(2))), x), x) - Dist((S(2)*c*g - h*(b + q))/q, Int(S(1)/((b + S(2)*c*x + q)*sqrt(d + e*x + f*x**S(2))), x), x)


def With623(a, c, d, e, f, g, h, x):
    q = Rt(-a*c, S(2))
    return Dist(-c*g/(S(2)*q) + h/S(2), Int(S(1)/((c*x + q)*sqrt(d + e*x + f*x**S(2))), x), x) + Dist(c*g/(S(2)*q) + h/S(2), Int(S(1)/((c*x - q)*sqrt(d + e*x + f*x**S(2))), x), x)


def With624(a, b, c, d, f, g, h, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist((S(2)*c*g - h*(b - q))/q, Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x - q)), x), x) - Dist((S(2)*c*g - h*(b + q))/q, Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x + q)), x), x)


def With625(a, b, c, d, e, f, g, h, x):
    q = Rt(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2), S(2))
    return Dist(S(1)/(S(2)*q), Int(Simp(-g*(-a*f + c*d - q) + h*(-a*e + b*d) - x*(g*(-b*f + c*e) - h*(-a*f + c*d + q)), x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x) - Dist(S(1)/(S(2)*q), Int(Simp(-g*(-a*f + c*d + q) + h*(-a*e + b*d) - x*(g*(-b*f + c*e) - h*(-a*f + c*d - q)), x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x)


def With626(a, c, d, e, f, g, h, x):
    q = Rt(a*c*e**S(2) + (-a*f + c*d)**S(2), S(2))
    return Dist(S(1)/(S(2)*q), Int(Simp(-a*e*h - g*(-a*f + c*d - q) + x*(-c*e*g + h*(-a*f + c*d + q)), x)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x) - Dist(S(1)/(S(2)*q), Int(Simp(-a*e*h - g*(-a*f + c*d + q) + x*(-c*e*g + h*(-a*f + c*d - q)), x)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x)


def With627(a, b, c, d, f, g, h, x):
    q = Rt(b**S(2)*d*f + (-a*f + c*d)**S(2), S(2))
    return Dist(S(1)/(S(2)*q), Int(Simp(b*d*h - g*(-a*f + c*d - q) + x*(b*f*g + h*(-a*f + c*d + q)), x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x), x) - Dist(S(1)/(S(2)*q), Int(Simp(b*d*h - g*(-a*f + c*d + q) + x*(b*f*g + h*(-a*f + c*d - q)), x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x), x)


def With628(a, b, c, d, e, f, g, h, x):
    s = Rt(-S(4)*a*c + b**S(2), S(2))
    t = Rt(-S(4)*d*f + e**S(2), S(2))
    return Dist(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + x*(e + t))*sqrt(b + S(2)*c*x + s)*sqrt(e + S(2)*f*x + t)/(sqrt(a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), Int((g + h*x)/(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + x*(e + t))*sqrt(b + S(2)*c*x + s)*sqrt(e + S(2)*f*x + t)), x), x)


def With629(a, b, c, d, f, g, h, x):
    s = Rt(-S(4)*a*c + b**S(2), S(2))
    t = Rt(-S(4)*d*f, S(2))
    return Dist(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + t*x)*sqrt(S(2)*f*x + t)*sqrt(b + S(2)*c*x + s)/(sqrt(d + f*x**S(2))*sqrt(a + b*x + c*x**S(2))), Int((g + h*x)/(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + t*x)*sqrt(S(2)*f*x + t)*sqrt(b + S(2)*c*x + s)), x), x)


def With630(a, b, c, d, e, f, g, h, x):
    q = S(3)**(S(2)/3)*(-c*h**S(2)/(-b*h + S(2)*c*g)**S(2))**(S(1)/3)
    return -Simp(S(3)*h*q*log((-S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(2)/3) + S(2)**(S(1)/3)*(S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(1)/3))/(S(2)*f), x) + Simp(h*q*log(d + e*x + f*x**S(2))/(S(2)*f), x) + Simp(sqrt(S(3))*h*q*ArcTan(-S(2)**(S(2)/3)*sqrt(S(3))*(-S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(2)/3)/(S(3)*(S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(1)/3)) + sqrt(S(3))/S(3))/f, x)


def replacement631(a, c, d, f, g, h, x):
    return Simp(S(2)**(S(1)/3)*h*log(d + f*x**S(2))/(S(4)*a**(S(1)/3)*f), x) - Simp(S(3)*S(2)**(S(1)/3)*h*log((S(1) - S(3)*h*x/g)**(S(2)/3) + S(2)**(S(1)/3)*(S(1) + S(3)*h*x/g)**(S(1)/3))/(S(4)*a**(S(1)/3)*f), x) + Simp(S(2)**(S(1)/3)*sqrt(S(3))*h*ArcTan(-S(2)**(S(2)/3)*sqrt(S(3))*(S(1) - S(3)*h*x/g)**(S(2)/3)/(S(3)*(S(1) + S(3)*h*x/g)**(S(1)/3)) + sqrt(S(3))/S(3))/(S(2)*a**(S(1)/3)*f), x)


def With632(a, b, c, d, e, f, g, h, x):
    q = -c/(-S(4)*a*c + b**S(2))
    return Dist((q*(a + b*x + c*x**S(2)))**(S(1)/3)/(a + b*x + c*x**S(2))**(S(1)/3), Int((g + h*x)/((d + e*x + f*x**S(2))*(a*q + b*q*x + c*q*x**S(2))**(S(1)/3)), x), x)


def replacement633(a, c, d, f, g, h, x):
    return Dist((S(1) + c*x**S(2)/a)**(S(1)/3)/(a + c*x**S(2))**(S(1)/3), Int((g + h*x)/((S(1) + c*x**S(2)/a)**(S(1)/3)*(d + f*x**S(2))), x), x)


def replacement634(a, b, c, d, e, f, g, h, p, q, x):
    return Int((g + h*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)


def replacement635(a, c, d, e, f, g, h, p, q, x):
    return Int((a + c*x**S(2))**p*(g + h*x)*(d + e*x + f*x**S(2))**q, x)


def replacement636(a, b, c, d, e, f, g, h, m, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u), x)


def replacement637(a, c, d, e, f, g, h, m, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + c*x**S(2))**p*(g + h*x)**m*(d + e*x + f*x**S(2))**q, x), x, u), x)


def replacement638(m, p, q, u, v, x, z):
    return Int(ExpandToSum(u, x)**p*ExpandToSum(v, x)**q*ExpandToSum(z, x)**m, x)


def replacement639(a, b, c, d, e, f, g, h, i, m, n, p, q, x):
    return Int((h + i*x)**q*(d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x)


def replacement640(a, b, c, d, e, f, g, h, i, m, n, p, q, x):
    return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(h + i*x)**q*(a + b*x + c*x**S(2))**p, x), x)


def replacement641(a, b, c, d, e, f, g, h, i, m, n, p, q, x):
    return Dist((d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m)), Int((h + i*x)**q*(d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x), x)


def replacement642(A, B, C, a, b, c, d, e, f, p, q, x):
    return Dist((c/f)**p, Int((A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x), x)


def replacement643(A, C, a, b, c, d, e, f, p, q, x):
    return Dist((c/f)**p, Int((A + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x), x)


def replacement644(A, B, C, a, b, c, d, e, f, p, q, x):
    return Dist(a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p)), Int((A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x), x)


def replacement645(A, C, a, b, c, d, e, f, p, q, x):
    return Dist(a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p)), Int((A + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x), x)


def replacement646(A, B, C, a, b, c, d, e, f, p, q, x):
    return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((b + S(2)*c*x)**(S(2)*p)*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x), x)


def replacement647(A, C, a, b, c, d, e, f, p, q, x):
    return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((A + C*x**S(2))*(b + S(2)*c*x)**(S(2)*p)*(d + e*x + f*x**S(2))**q, x), x)


def replacement648(A, B, C, a, b, c, d, f, p, q, x):
    return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q*(A + B*x + C*x**S(2)), x), x)


def replacement649(A, C, a, b, c, d, f, p, q, x):
    return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p), Int((A + C*x**S(2))*(b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q, x), x)


def replacement650(A, B, C, a, b, c, d, e, f, p, q, x):
    return Int(ExpandIntegrand((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x)


def replacement651(A, C, a, b, c, d, e, f, p, q, x):
    return Int(ExpandIntegrand((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x)


def replacement652(A, B, C, a, c, d, e, f, p, q, x):
    return Int(ExpandIntegrand((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x), x)


def replacement653(A, C, a, c, d, e, f, p, q, x):
    return Int(ExpandIntegrand((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x)


def replacement654(A, B, C, a, b, c, d, e, f, p, q, x):
    return -Dist(S(1)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-d*(C*(S(2)*a*c - b**S(2)*(p + S(2))) + c*(S(2)*p + S(3))*(-S(2)*A*c + B*b)) + e*q*(A*b*c - S(2)*B*a*c + C*a*b) - f*x**S(2)*(C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + S(2)*q + S(3))) + x*(-e*(C*(S(2)*a*c*(q + S(1)) - b**S(2)*(p + q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + q + S(3))) + S(2)*f*q*(A*b*c - S(2)*B*a*c + C*a*b)), x), x), x) + Simp((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(A*b*c - S(2)*B*a*c + C*a*b - x*(-C*(-S(2)*a*c + b**S(2)) + c*(-S(2)*A*c + B*b)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement655(A, C, a, b, c, d, e, f, p, q, x):
    return -Dist(S(1)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*(b*e*q + S(2)*c*d*(S(2)*p + S(3))) - C*(-a*b*e*q + S(2)*a*c*d - b**S(2)*d*(p + S(2))) - f*x**S(2)*(-S(2)*A*c**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2)))) + x*(S(2)*A*c*(b*f*q + c*e*(S(2)*p + q + S(3))) + C*(S(2)*a*b*f*q - S(2)*a*c*e*(q + S(1)) + b**S(2)*e*(p + q + S(2)))), x), x), x) + Simp((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(A*b*c + C*a*b + x*(S(2)*A*c**S(2) + C*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement656(A, B, C, a, c, d, e, f, p, q, x):
    return -Dist(-S(1)/(S(2)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*d*(S(2)*p + S(3)) - a*(B*e*q + C*d) - f*x**S(2)*(-A*c*(S(2)*p + S(2)*q + S(3)) + C*a*(S(2)*q + S(1))) + x*(A*c*e*(S(2)*p + q + S(3)) - a*(S(2)*B*f*q + C*e*(q + S(1)))), x), x), x) + Simp((a + c*x**S(2))**(p + S(1))*(B*a - x*(A*c - C*a))*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))), x)


def replacement657(A, C, a, c, d, e, f, p, q, x):
    return Dist(S(1)/(S(2)*a*c*(p + S(1))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*d*(S(2)*p + S(3)) - C*a*d - f*x**S(2)*(-A*c*(S(2)*p + S(2)*q + S(3)) + C*a*(S(2)*q + S(1))) + x*(A*c*e*(S(2)*p + q + S(3)) - C*a*e*(q + S(1))), x), x), x) - Simp(x*(a + c*x**S(2))**(p + S(1))*(A*c - C*a)*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))), x)


def replacement658(A, B, C, a, b, c, d, f, p, q, x):
    return -Dist(S(1)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-d*(C*(S(2)*a*c - b**S(2)*(p + S(2))) + c*(S(2)*p + S(3))*(-S(2)*A*c + B*b)) + S(2)*f*q*x*(A*b*c - S(2)*B*a*c + C*a*b) - f*x**S(2)*(C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + S(2)*q + S(3))), x), x), x) + Simp((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*c - S(2)*B*a*c + C*a*b - x*(-C*(-S(2)*a*c + b**S(2)) + c*(-S(2)*A*c + B*b)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement659(A, C, a, b, c, d, f, p, q, x):
    return -Dist(S(1)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(S(2)*A*c**S(2)*d*(S(2)*p + S(3)) - C*(S(2)*a*c*d - b**S(2)*d*(p + S(2))) - f*x**S(2)*(-S(2)*A*c**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2)))) + x*(S(2)*A*b*c*f*q + S(2)*C*a*b*f*q), x), x), x) + Simp((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*c + C*a*b + x*(S(2)*A*c**S(2) + C*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement660(A, B, C, a, b, c, d, e, f, p, q, x):
    return Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))), Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e)) - e*((A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - x*(S(2)*f*((A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e))) + (p + S(1))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))*(-S(2)*A*c + B*b - S(2)*C*a) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e)), x), x), x) + Simp((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(c*x*(A*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) - B*(a*b*f - S(2)*a*c*e + b*c*d) + C*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))), x)


def replacement661(A, C, a, b, c, d, e, f, p, q, x):
    return Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))), Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e)) - e*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - x*(S(2)*f*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e)), x), x), x) + Simp((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + c*x*(A*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + C*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))), x)


def replacement662(A, B, C, a, c, d, e, f, p, q, x):
    return Dist(-S(1)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(S(2)*p + S(2)*q + S(5)) - e*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a))*(p + q + S(2)) - x*(c*e*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(S(2)*p + q + S(4)) + S(2)*f*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a))*(p + q + S(2))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + (S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x), x) + Simp(-(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a) + c*x*(A*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*B*a*c*e - S(2)*C*a*(-a*f + c*d)))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))), x)


def replacement663(A, C, a, c, d, e, f, p, q, x):
    return Dist(-S(1)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))), Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-S(2)*a*c*e**S(2)*(A*c - C*a)*(p + q + S(2)) - c*f*x**S(2)*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(S(2)*p + S(2)*q + S(5)) - x*(S(4)*a*c*e*f*(A*c - C*a)*(p + q + S(2)) + c*e*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(S(2)*p + q + S(4))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + (S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x), x) + Simp(-(a + c*x**S(2))**(p + S(1))*(S(2)*a*c*e*(A*c - C*a) + c*x*(A*(-S(2)*a*c*f + S(2)*c**S(2)*d) - S(2)*C*a*(-a*f + c*d)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))), x)


def replacement664(A, B, C, a, b, c, d, f, p, q, x):
    return Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))), Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)) - x*(-b*f*(p + S(1))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)) + S(2)*f*(-b*(A*c - C*a)*(a*f + c*d) + (A*b - B*a)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (p + S(1))*(b**S(2)*d*f + (-a*f + c*d)**S(2))*(-S(2)*A*c + B*b - S(2)*C*a) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)), x), x), x) + Simp((d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*(A*c - C*a)*(a*f + c*d) + c*x*(A*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - B*(a*b*f + b*c*d) + C*(-S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*b - B*a)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))), x)


def replacement665(A, C, a, b, c, d, f, p, q, x):
    return Dist(S(1)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))), Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)) - x*(-b*f*(p + S(1))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)) + S(2)*f*(A*b*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - b*(A*c - C*a)*(a*f + c*d))*(p + q + S(2))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)), x), x), x) + Simp((d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - b*(A*c - C*a)*(a*f + c*d) + c*x*(A*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + C*(-S(2)*a*(-a*f + c*d) + b**S(2)*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))), x)


def replacement666(A, B, C, a, b, c, d, e, f, p, q, x):
    return -Dist(S(1)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(p*(-a*e + b*d)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(p*(-b*f + c*e)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + x*(S(2)*p*(-a*f + c*d)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (C*e*f*p*(-S(4)*a*c + b**S(2)) - b*c*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))) + f*(-S(2)*A*f + B*e)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1)), x), x), x) + Simp((a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))*(B*c*f*(S(2)*p + S(2)*q + S(3)) + S(2)*C*c*f*x*(p + q + S(1)) + C*(b*f*p - c*e*(S(2)*p + q + S(2))))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), x)


def replacement667(A, C, a, b, c, d, e, f, p, q, x):
    return -Dist(S(1)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(p*(-a*e + b*d)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + x**S(2)*(p*(-b*f + c*e)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(S(2)*p*(-a*f + c*d)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + (C*e*f*p*(-S(4)*a*c + b**S(2)) - b*c*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2)))))*(p + q + S(1)), x), x), x) + Simp((S(2)*C*c*f*x*(p + q + S(1)) + C*(b*f*p - c*e*(S(2)*p + q + S(2))))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), x)


def replacement668(A, B, C, a, c, d, e, f, p, q, x):
    return -Dist(S(1)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*c*(C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))) + f*(-S(2)*A*f + B*e)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1)) - a*e*p*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(c*e*p*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (-S(4)*C*a*c*f**S(2)*p - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + x*(-S(4)*C*a*c*e*f*p*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3)))), x), x), x) + Simp((a + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*c*e*(S(2)*p + q + S(2)) + S(2)*C*c*f*x*(p + q + S(1)))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), x)


def replacement669(A, C, a, c, d, e, f, p, q, x):
    return -Dist(S(1)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))))*(p + q + S(1)) - a*e*p*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(c*e*p*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3))) + (-S(4)*C*a*c*f**S(2)*p - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-S(4)*C*a*c*e*f*p*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3)))), x), x), x) + Simp((a + c*x**S(2))**p*(-C*c*e*(S(2)*p + q + S(2)) + S(2)*C*c*f*x*(p + q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), x)


def replacement670(A, B, C, a, b, c, d, f, p, q, x):
    return -Dist(S(1)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(b*d*p*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1))) + x**S(2)*(-b*f*p*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1))) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-b*c*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1)))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + S(2)*C*d*f))*(p + q + S(1)), x), x), x) + Simp((d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**p*(B*c*f*(S(2)*p + S(2)*q + S(3)) + C*b*f*p + S(2)*C*c*f*x*(p + q + S(1)))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), x)


def replacement671(A, C, a, b, c, d, f, p, q, x):
    return -Dist(S(1)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-C*b**S(2)*d*f*p*(q + S(1)) + x**S(2)*(C*b**S(2)*f**S(2)*p*(q + S(1)) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-S(2)*C*b*f*p*(q + S(1))*(-a*f + c*d) - b*c*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + S(2)*C*d*f))*(p + q + S(1)), x), x), x) + Simp((d + f*x**S(2))**(q + S(1))*(C*b*f*p + S(2)*C*c*f*x*(p + q + S(1)))*(a + b*x + c*x**S(2))**p/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))), x)


def With672(A, B, C, a, b, c, d, e, f, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
    if NonzeroQ(q):
        return True
    return False


def replacement672(A, B, C, a, b, c, d, e, f, x):

    q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
    return Dist(S(1)/q, Int((-A*a*c*f + A*b**S(2)*f - A*b*c*e + A*c**S(2)*d - B*a*b*f + B*a*c*e + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - A*c*e - B*a*f + B*c*d + C*a*e - C*b*d))/(a + b*x + c*x**S(2)), x), x) + Dist(S(1)/q, Int((A*a*f**S(2) - A*b*e*f - A*c*d*f + A*c*e**S(2) + B*b*d*f - B*c*d*e - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - A*c*e - B*a*f + B*c*d + C*a*e - C*b*d))/(d + e*x + f*x**S(2)), x), x)


def With673(A, C, a, b, c, d, e, f, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
    if NonzeroQ(q):
        return True
    return False


def replacement673(A, C, a, b, c, d, e, f, x):

    q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
    return Dist(S(1)/q, Int((-A*a*c*f + A*b**S(2)*f - A*b*c*e + A*c**S(2)*d + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - A*c*e + C*a*e - C*b*d))/(a + b*x + c*x**S(2)), x), x) + Dist(S(1)/q, Int((A*a*f**S(2) - A*b*e*f - A*c*d*f + A*c*e**S(2) - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - A*c*e + C*a*e - C*b*d))/(d + e*x + f*x**S(2)), x), x)


def With674(A, B, C, a, b, c, d, f, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
    if NonzeroQ(q):
        return True
    return False


def replacement674(A, B, C, a, b, c, d, f, x):

    q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
    return Dist(S(1)/q, Int((A*a*f**S(2) - A*c*d*f + B*b*d*f - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - B*a*f + B*c*d - C*b*d))/(d + f*x**S(2)), x), x) + Dist(S(1)/q, Int((-A*a*c*f + A*b**S(2)*f + A*c**S(2)*d - B*a*b*f + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - B*a*f + B*c*d - C*b*d))/(a + b*x + c*x**S(2)), x), x)


def With675(A, C, a, b, c, d, f, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
    if NonzeroQ(q):
        return True
    return False


def replacement675(A, C, a, b, c, d, f, x):

    q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
    return Dist(S(1)/q, Int((A*a*f**S(2) - A*c*d*f - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - C*b*d))/(d + f*x**S(2)), x), x) + Dist(S(1)/q, Int((-A*a*c*f + A*b**S(2)*f + A*c**S(2)*d + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - C*b*d))/(a + b*x + c*x**S(2)), x), x)


def replacement676(A, B, C, a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int((A*c - C*a + x*(B*c - C*b))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x) + Dist(C/c, Int(S(1)/sqrt(d + e*x + f*x**S(2)), x), x)


def replacement677(A, C, a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int((A*c - C*a - C*b*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x) + Dist(C/c, Int(S(1)/sqrt(d + e*x + f*x**S(2)), x), x)


def replacement678(A, B, C, a, c, d, e, f, x):
    return Dist(S(1)/c, Int((A*c + B*c*x - C*a)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x) + Dist(C/c, Int(S(1)/sqrt(d + e*x + f*x**S(2)), x), x)


def replacement679(A, C, a, c, d, e, f, x):
    return Dist(C/c, Int(S(1)/sqrt(d + e*x + f*x**S(2)), x), x) + Dist((A*c - C*a)/c, Int(S(1)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x), x)


def replacement680(A, B, C, a, b, c, d, f, x):
    return Dist(S(1)/c, Int((A*c - C*a + x*(B*c - C*b))/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x), x) + Dist(C/c, Int(S(1)/sqrt(d + f*x**S(2)), x), x)


def replacement681(A, C, a, b, c, d, f, x):
    return Dist(S(1)/c, Int((A*c - C*a - C*b*x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x), x) + Dist(C/c, Int(S(1)/sqrt(d + f*x**S(2)), x), x)


def replacement682(A, B, C, a, b, c, d, e, f, p, q, x):
    return Int((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)


def replacement683(A, C, a, b, c, d, e, f, p, q, x):
    return Int((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)


def replacement684(A, B, C, a, c, d, e, f, p, q, x):
    return Int((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x)


def replacement685(A, C, a, c, d, e, f, p, q, x):
    return Int((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)


def replacement686(A, B, C, a, b, c, d, e, f, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u), x)


def replacement687(A, B, a, b, c, d, e, f, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((A + B*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u), x)


def replacement688(A, C, a, b, c, d, e, f, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u), x)


def replacement689(A, B, C, a, c, d, e, f, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x), x, u), x)


def replacement690(A, B, a, c, d, e, f, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((A + B*x)*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u), x)


def replacement691(A, C, a, c, d, e, f, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u), x)
