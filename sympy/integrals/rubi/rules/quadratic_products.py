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

def quadratic_products(rubi):
    from sympy.integrals.rubi.constraints import cons45, cons2, cons3, cons7, cons225, cons5, cons226, cons128, cons227, cons228, cons13, cons163, cons229, cons137, cons230, cons231, cons232, cons233, cons234, cons235, cons68, cons69, cons47, cons236, cons27, cons48, cons21, cons237, cons238, cons239, cons240, cons66, cons241, cons242, cons243, cons244, cons146, cons245, cons246, cons247, cons248, cons249, cons250, cons251, cons252, cons166, cons253, cons31, cons168, cons254, cons255, cons94, cons147, cons256, cons38, cons257, cons258, cons41, cons17, cons259, cons260, cons261, cons262, cons263, cons264, cons265, cons266, cons267, cons43, cons268, cons54, cons269, cons270, cons271, cons272, cons273, cons274, cons275, cons276, cons277, cons278, cons279, cons280, cons281, cons282, cons283, cons284, cons285, cons286, cons18, cons287, cons288, cons289, cons290, cons291, cons292, cons293, cons294, cons295, cons296, cons297, cons298, cons299, cons300, cons301, cons302, cons303, cons304, cons305, cons306, cons307, cons308, cons309, cons310, cons311, cons312, cons313, cons84, cons85, cons314, cons315, cons316, cons125, cons208, cons317, cons318, cons319, cons62, cons320, cons321, cons322, cons323, cons324, cons4, cons325, cons326, cons327, cons139, cons328, cons329, cons330, cons331, cons150, cons332, cons148, cons333, cons196, cons334, cons335, cons336, cons337, cons338, cons89, cons339, cons340, cons341, cons88, cons87, cons342, cons343, cons344, cons126, cons345, cons346, cons207, cons347, cons348, cons349, cons350, cons351, cons352, cons353, cons354, cons355, cons356, cons357, cons358, cons359, cons360, cons361, cons362, cons363, cons364, cons365, cons366, cons367, cons368, cons369, cons370, cons371, cons372, cons373, cons374, cons375, cons149, cons376, cons124, cons377, cons93, cons23, cons165, cons73, cons378, cons80, cons379, cons380, cons381, cons382, cons383, cons384, cons385, cons50, cons386, cons387, cons388, cons389, cons390, cons391, cons392, cons393, cons394, cons395, cons396, cons397, cons398, cons399, cons400, cons401, cons402, cons403, cons404, cons405, cons406, cons407, cons408, cons409, cons410, cons411, cons412, cons413, cons414, cons415, cons416, cons209, cons417, cons418, cons419, cons420, cons421, cons422, cons423, cons424, cons425, cons426, cons427, cons428, cons429, cons430, cons431, cons220, cons432, cons433, cons434, cons435, cons436, cons437, cons438, cons439, cons440, cons441, cons442, cons443, cons444, cons445, cons446, cons447, cons448, cons449, cons450, cons451, cons452, cons453, cons224, cons34, cons35, cons36, cons454, cons455, cons456, cons457, cons458

    pattern189 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons45)
    def replacement189(x, c, b, a):
        return (b/S(2) + c*x)*Int(S(1)/(b/S(2) + c*x), x)/sqrt(a + b*x + c*x**S(2))
    rule189 = ReplacementRule(pattern189, replacement189)
    rubi.add(rule189.pattern, label = 189)

    pattern190 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons45, cons225)
    def replacement190(p, x, b, c, a):
        return (b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1)))
    rule190 = ReplacementRule(pattern190, replacement190)
    rubi.add(rule190.pattern, label = 190)

    def With191(p, x, b, c, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return c**(-p)*Int(Simp(b/S(2) + c*x - q/S(2), x)**p*Simp(b/S(2) + c*x + q/S(2), x)**p, x)
    pattern191 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons226, cons128, cons227 )
    rule191 = ReplacementRule(pattern191, With191)
    rubi.add(rule191.pattern, label = 191)

    pattern192 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons226, cons128, cons228)
    def replacement192(p, x, b, c, a):
        return Int(ExpandIntegrand((a + b*x + c*x**S(2))**p, x), x)
    rule192 = ReplacementRule(pattern192, replacement192)
    rubi.add(rule192.pattern, label = 192)

    pattern193 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons226, cons13, cons163, cons229)
    def replacement193(p, x, b, c, a):
        return -p*(-S(4)*a*c + b**S(2))*Int((a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*c*(S(2)*p + S(1))) + (b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1)))
    rule193 = ReplacementRule(pattern193, replacement193)
    rubi.add(rule193.pattern, label = 193)

    pattern194 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(-3)/2), x_), cons2, cons3, cons7, cons226)
    def replacement194(x, c, b, a):
        return -S(2)*(b + S(2)*c*x)/((-S(4)*a*c + b**S(2))*sqrt(a + b*x + c*x**S(2)))
    rule194 = ReplacementRule(pattern194, replacement194)
    rubi.add(rule194.pattern, label = 194)

    pattern195 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons226, cons13, cons137, cons230, cons229)
    def replacement195(p, x, b, c, a):
        return -S(2)*c*(S(2)*p + S(3))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (b + S(2)*c*x)*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule195 = ReplacementRule(pattern195, replacement195)
    rubi.add(rule195.pattern, label = 195)

    def With196(x, c, b, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return c*Int(S(1)/Simp(b/S(2) + c*x - q/S(2), x), x)/q - c*Int(S(1)/Simp(b/S(2) + c*x + q/S(2), x), x)/q
    pattern196 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons226, cons231, cons227 )
    rule196 = ReplacementRule(pattern196, With196)
    rubi.add(rule196.pattern, label = 196)

    def With197(x, c, b, a):
        q = -S(4)*a*c/b**S(2) + S(1)
        if And(RationalQ(q), Or(EqQ(q**S(2), S(1)), Not(RationalQ(-S(4)*a*c + b**S(2))))):
            return -S(2)*Subst(Int(S(1)/(q - x**S(2)), x), x, S(1) + S(2)*c*x/b)/b
        return False
    pattern197 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons226, CustomConstraint(With197))
    rule197 = ReplacementRule(pattern197, With197)
    rubi.add(rule197.pattern, label = 197)

    pattern198 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons226)
    def replacement198(x, c, b, a):
        return -S(2)*Subst(Int(S(1)/Simp(-S(4)*a*c + b**S(2) - x**S(2), x), x), x, b + S(2)*c*x)
    rule198 = ReplacementRule(pattern198, replacement198)
    rubi.add(rule198.pattern, label = 198)

    pattern199 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons5, cons232)
    def replacement199(p, x, b, c, a):
        return (-S(4)*c/(-S(4)*a*c + b**S(2)))**(-p)*Subst(Int(Simp(-x**S(2)/(-S(4)*a*c + b**S(2)) + S(1), x)**p, x), x, b + S(2)*c*x)/(S(2)*c)
    rule199 = ReplacementRule(pattern199, replacement199)
    rubi.add(rule199.pattern, label = 199)

    pattern200 = Pattern(Integral(S(1)/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons3, cons7, cons233)
    def replacement200(x, c, b):
        return S(2)*Subst(Int(S(1)/(-c*x**S(2) + S(1)), x), x, x/sqrt(b*x + c*x**S(2)))
    rule200 = ReplacementRule(pattern200, replacement200)
    rubi.add(rule200.pattern, label = 200)

    pattern201 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons226)
    def replacement201(x, c, b, a):
        return S(2)*Subst(Int(S(1)/(S(4)*c - x**S(2)), x), x, (b + S(2)*c*x)/sqrt(a + b*x + c*x**S(2)))
    rule201 = ReplacementRule(pattern201, replacement201)
    rubi.add(rule201.pattern, label = 201)

    pattern202 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons3, cons7, cons13, cons234)
    def replacement202(p, x, c, b):
        return (-c*(b*x + c*x**S(2))/b**S(2))**(-p)*(b*x + c*x**S(2))**p*Int((-c*x/b - c**S(2)*x**S(2)/b**S(2))**p, x)
    rule202 = ReplacementRule(pattern202, replacement202)
    rubi.add(rule202.pattern, label = 202)

    def With203(p, x, b, c, a):
        d = Denominator(p)
        if LessEqual(S(3), d, S(4)):
            return d*sqrt((b + S(2)*c*x)**S(2))*Subst(Int(x**(d*(p + S(1)) + S(-1))/sqrt(-S(4)*a*c + b**S(2) + S(4)*c*x**d), x), x, (a + b*x + c*x**S(2))**(S(1)/d))/(b + S(2)*c*x)
        return False
    pattern203 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons226, cons13, CustomConstraint(With203))
    rule203 = ReplacementRule(pattern203, With203)
    rubi.add(rule203.pattern, label = 203)

    def With204(p, x, b, c, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -((-b - S(2)*c*x + q)/(S(2)*q))**(-p + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Hypergeometric2F1(-p, p + S(1), p + S(2), (b + S(2)*c*x + q)/(S(2)*q))/(q*(p + S(1)))
    pattern204 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons5, cons226, cons235 )
    rule204 = ReplacementRule(pattern204, With204)
    rubi.add(rule204.pattern, label = 204)

    pattern205 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons5, cons68, cons69)
    def replacement205(p, x, b, c, u, a):
        return Subst(Int((a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1))
    rule205 = ReplacementRule(pattern205, replacement205)
    rubi.add(rule205.pattern, label = 205)

    pattern206 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons45, cons47, cons236)
    def replacement206(m, p, e, x, b, d, c, a):
        return c**(-m/S(2) + S(-1)/2)*e**m*(a + b*x + c*x**S(2))**(m/S(2) + p + S(1)/2)/(m + S(2)*p + S(1))
    rule206 = ReplacementRule(pattern206, replacement206)
    rubi.add(rule206.pattern, label = 206)

    pattern207 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons45, cons47, cons237)
    def replacement207(m, p, x, e, b, d, c, a):
        return (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*log(RemoveContent(d + e*x, x))/e
    rule207 = ReplacementRule(pattern207, replacement207)
    rubi.add(rule207.pattern, label = 207)

    pattern208 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons45, cons47, cons238)
    def replacement208(m, p, e, x, b, d, c, a):
        return (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1)))
    rule208 = ReplacementRule(pattern208, replacement208)
    rubi.add(rule208.pattern, label = 208)

    pattern209 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons45, cons239, cons240, cons66)
    def replacement209(m, p, e, x, b, d, c, a):
        return -(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(1))*(-b*e + S(2)*c*d))
    rule209 = ReplacementRule(pattern209, replacement209)
    rubi.add(rule209.pattern, label = 209)

    pattern210 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0)))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons45, cons239)
    def replacement210(x, e, b, d, c, a):
        return sqrt(a + b*x + c*x**S(2))*Int((b + S(2)*c*x)/(d + e*x)**S(2), x)/(b + S(2)*c*x)
    rule210 = ReplacementRule(pattern210, replacement210)
    rubi.add(rule210.pattern, label = 210)

    pattern211 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons45, cons239, cons241)
    def replacement211(m, x, e, b, d, c, a):
        return (d + e*x)**(m + S(1))*sqrt(a + b*x + c*x**S(2))/(e*(m + S(2))) - (-b*e + S(2)*c*d)*sqrt(a + b*x + c*x**S(2))*Int((d + e*x)**m, x)/(e*(b + S(2)*c*x)*(m + S(2)))
    rule211 = ReplacementRule(pattern211, replacement211)
    rubi.add(rule211.pattern, label = 211)

    pattern212 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons45, cons239)
    def replacement212(x, e, b, d, c, a):
        return -S(4)*c*e*sqrt(a + b*x + c*x**S(2))/((d + e*x)*(-b*e + S(2)*c*d)**S(2)) + S(2)*c*Int(S(1)/((d + e*x)*sqrt(a + b*x + c*x**S(2))), x)/(-b*e + S(2)*c*d)
    rule212 = ReplacementRule(pattern212, replacement212)
    rubi.add(rule212.pattern, label = 212)

    pattern213 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons45, cons239, cons242, cons241)
    def replacement213(m, p, x, e, b, d, c, a):
        return -S(2)*c*e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)**S(2)*(m*p + S(-1))) - (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(2))*(-b*e + S(2)*c*d))
    rule213 = ReplacementRule(pattern213, replacement213)
    rubi.add(rule213.pattern, label = 213)

    pattern214 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons45, cons239, cons243)
    def replacement214(p, x, e, b, d, c, a):
        return e*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))) + (-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c)
    rule214 = ReplacementRule(pattern214, replacement214)
    rubi.add(rule214.pattern, label = 214)

    pattern215 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons45, cons239, cons244, cons146, cons245, cons246)
    def replacement215(m, p, x, e, b, d, c, a):
        return (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(1))) + p*(S(2)*p + S(-1))*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(1)))
    rule215 = ReplacementRule(pattern215, replacement215)
    rubi.add(rule215.pattern, label = 215)

    pattern216 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons45, cons239, cons244, cons146, cons247, cons246, cons248)
    def replacement216(m, p, x, e, b, d, c, a):
        return S(2)*c*p*(S(2)*p + S(-1))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(1))*(m + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))/(e**S(2)*(m + S(1))*(m + S(2)))
    rule216 = ReplacementRule(pattern216, replacement216)
    rubi.add(rule216.pattern, label = 216)

    pattern217 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons45, cons239, cons13, cons163, cons249, cons238, cons248, cons250, cons251, cons246)
    def replacement217(m, p, x, e, b, d, c, a):
        return (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(-b*e + S(2)*c*d)*(a + b*x + c*x**S(2))**(p + S(-1))/(S(2)*c*e**S(2)*(m + S(2)*p)*(m + S(2)*p + S(1))) + p*(S(2)*p + S(-1))*(-b*e + S(2)*c*d)**S(2)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*c*e**S(2)*(m + S(2)*p)*(m + S(2)*p + S(1)))
    rule217 = ReplacementRule(pattern217, replacement217)
    rubi.add(rule217.pattern, label = 217)

    pattern218 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons45, cons239, cons244, cons137, cons252, cons246)
    def replacement218(m, p, e, x, b, d, c, a):
        return e**S(2)*m*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(m + S(2)*p + S(2))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)) + (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((S(2)*p + S(1))*(-b*e + S(2)*c*d))
    rule218 = ReplacementRule(pattern218, replacement218)
    rubi.add(rule218.pattern, label = 218)

    pattern219 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons45, cons239, cons244, cons137, cons166, cons246)
    def replacement219(m, p, x, e, b, d, c, a):
        return e**S(2)*m*(m + S(-1))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(S(2)*c*(p + S(1))*(S(2)*p + S(1))) - e*m*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))*(S(2)*p + S(1))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1)))
    rule219 = ReplacementRule(pattern219, replacement219)
    rubi.add(rule219.pattern, label = 219)

    pattern220 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons45, cons239, cons244, cons137, cons253, cons246)
    def replacement220(m, p, x, e, b, d, c, a):
        return S(2)*c*e**S(2)*(m + S(2)*p + S(2))*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)**S(2)) - S(2)*c*e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(m + S(2)*p + S(2))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)**S(2)) + (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((S(2)*p + S(1))*(-b*e + S(2)*c*d))
    rule220 = ReplacementRule(pattern220, replacement220)
    rubi.add(rule220.pattern, label = 220)

    pattern221 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons45, cons239, cons31, cons168, cons238, cons254, cons255)
    def replacement221(m, p, e, x, b, d, c, a):
        return m*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*(m + S(2)*p + S(1))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p/(S(2)*c*(m + S(2)*p + S(1)))
    rule221 = ReplacementRule(pattern221, replacement221)
    rubi.add(rule221.pattern, label = 221)

    pattern222 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons45, cons239, cons31, cons94, cons246)
    def replacement222(m, p, x, e, b, d, c, a):
        return S(2)*c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((m + S(1))*(-b*e + S(2)*c*d)) - (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(1))*(-b*e + S(2)*c*d))
    rule222 = ReplacementRule(pattern222, replacement222)
    rubi.add(rule222.pattern, label = 222)

    pattern223 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons45, cons147, cons239)
    def replacement223(m, p, e, x, b, d, c, a):
        return c**(-IntPart(p))*(b/S(2) + c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b/S(2) + c*x)**(S(2)*p)*(d + e*x)**m, x)
    rule223 = ReplacementRule(pattern223, replacement223)
    rubi.add(rule223.pattern, label = 223)

    pattern224 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons226, cons256, cons38)
    def replacement224(m, p, e, x, b, d, c, a):
        return Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x)
    rule224 = ReplacementRule(pattern224, replacement224)
    rubi.add(rule224.pattern, label = 224)

    pattern225 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons5, cons257, cons258)
    def replacement225(m, p, e, x, d, c, a):
        return Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x)
    rule225 = ReplacementRule(pattern225, replacement225)
    rubi.add(rule225.pattern, label = 225)

    pattern226 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons256, cons147, cons41)
    def replacement226(m, p, e, x, b, d, c, a):
        return e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1)))
    rule226 = ReplacementRule(pattern226, replacement226)
    rubi.add(rule226.pattern, label = 226)

    pattern227 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons21, cons5, cons257, cons147, cons41)
    def replacement227(m, p, x, e, d, c, a):
        return e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(p + S(1)))
    rule227 = ReplacementRule(pattern227, replacement227)
    rubi.add(rule227.pattern, label = 227)

    pattern228 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons256, cons147, cons240)
    def replacement228(m, p, e, x, b, d, c, a):
        return e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-b*e + S(2)*c*d))
    rule228 = ReplacementRule(pattern228, replacement228)
    rubi.add(rule228.pattern, label = 228)

    pattern229 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons5, cons257, cons147, cons240)
    def replacement229(m, p, x, e, d, c, a):
        return e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(p + S(1)))
    rule229 = ReplacementRule(pattern229, replacement229)
    rubi.add(rule229.pattern, label = 229)

    pattern230 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons226, cons256, cons147, cons13, cons137)
    def replacement230(p, e, x, b, d, c, a):
        return -e**S(2)*(p + S(2))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1)))
    rule230 = ReplacementRule(pattern230, replacement230)
    rubi.add(rule230.pattern, label = 230)

    pattern231 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**S(2), x_), cons2, cons7, cons27, cons48, cons5, cons257, cons147, cons13, cons137)
    def replacement231(p, x, e, d, c, a):
        return -e**S(2)*(p + S(2))*Int((a + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)/(c*(p + S(1)))
    rule231 = ReplacementRule(pattern231, replacement231)
    rubi.add(rule231.pattern, label = 231)

    pattern232 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons226, cons256, cons147, cons17, cons13, cons259, cons260, cons261)
    def replacement232(m, p, e, x, b, d, c, a):
        return Int((a/d + c*x/e)**(-m)*(a + b*x + c*x**S(2))**(m + p), x)
    rule232 = ReplacementRule(pattern232, replacement232)
    rubi.add(rule232.pattern, label = 232)

    pattern233 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons21, cons5, cons257, cons147, cons17, cons13, cons259, cons260, cons261)
    def replacement233(m, p, x, e, d, c, a):
        return a**(-m)*d**(S(2)*m)*Int((a + c*x**S(2))**(m + p)*(d - e*x)**(-m), x)
    rule233 = ReplacementRule(pattern233, replacement233)
    rubi.add(rule233.pattern, label = 233)

    pattern234 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons256, cons147, cons262)
    def replacement234(m, p, e, x, b, d, c, a):
        return e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + (m + p)*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*(m + S(2)*p + S(1)))
    rule234 = ReplacementRule(pattern234, replacement234)
    rubi.add(rule234.pattern, label = 234)

    pattern235 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons5, cons257, cons147, cons262)
    def replacement235(m, p, x, e, d, c, a):
        return S(2)*d*(m + p)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1)), x)/(m + S(2)*p + S(1)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1)))
    rule235 = ReplacementRule(pattern235, replacement235)
    rubi.add(rule235.pattern, label = 235)

    pattern236 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons256, cons147, cons263)
    def replacement236(m, p, e, x, b, d, c, a):
        return c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((-b*e + S(2)*c*d)*(m + p + S(1))) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1)))
    rule236 = ReplacementRule(pattern236, replacement236)
    rubi.add(rule236.pattern, label = 236)

    pattern237 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons21, cons5, cons257, cons147, cons263)
    def replacement237(m, p, x, e, d, c, a):
        return (m + S(2)*p + S(2))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*d*(m + p + S(1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(m + p + S(1)))
    rule237 = ReplacementRule(pattern237, replacement237)
    rubi.add(rule237.pattern, label = 237)

    pattern238 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons256)
    def replacement238(e, x, b, d, c, a):
        return S(2)*e*Subst(Int(S(1)/(-b*e + S(2)*c*d + e**S(2)*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))/sqrt(d + e*x))
    rule238 = ReplacementRule(pattern238, replacement238)
    rubi.add(rule238.pattern, label = 238)

    pattern239 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons257)
    def replacement239(x, e, d, c, a):
        return S(2)*e*Subst(Int(S(1)/(S(2)*c*d + e**S(2)*x**S(2)), x), x, sqrt(a + c*x**S(2))/sqrt(d + e*x))
    rule239 = ReplacementRule(pattern239, replacement239)
    rubi.add(rule239.pattern, label = 239)

    pattern240 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons226, cons256, cons244, cons163, cons264, cons253, cons246)
    def replacement240(m, p, e, x, b, d, c, a):
        return -c*p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + p + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + p + S(1)))
    rule240 = ReplacementRule(pattern240, replacement240)
    rubi.add(rule240.pattern, label = 240)

    pattern241 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons257, cons244, cons163, cons264, cons253, cons246)
    def replacement241(m, p, x, e, d, c, a):
        return -c*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2)), x)/(e**S(2)*(m + p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + p + S(1)))
    rule241 = ReplacementRule(pattern241, replacement241)
    rubi.add(rule241.pattern, label = 241)

    pattern242 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons226, cons256, cons244, cons163, cons265, cons238, cons246)
    def replacement242(m, p, e, x, b, d, c, a):
        return (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - p*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(2)*p + S(1)))
    rule242 = ReplacementRule(pattern242, replacement242)
    rubi.add(rule242.pattern, label = 242)

    pattern243 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons257, cons244, cons163, cons265, cons238, cons246)
    def replacement243(m, p, x, e, d, c, a):
        return -S(2)*c*d*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1)), x)/(e**S(2)*(m + S(2)*p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(2)*p + S(1)))
    rule243 = ReplacementRule(pattern243, replacement243)
    rubi.add(rule243.pattern, label = 243)

    pattern244 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons226, cons256, cons244, cons137, cons252, cons246)
    def replacement244(m, p, e, x, b, d, c, a):
        return -(-b*e + S(2)*c*d)*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (d + e*x)**m*(-b*e + S(2)*c*d)*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule244 = ReplacementRule(pattern244, replacement244)
    rubi.add(rule244.pattern, label = 244)

    pattern245 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons257, cons244, cons137, cons252, cons246)
    def replacement245(m, p, x, e, d, c, a):
        return d*(m + S(2)*p + S(2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*(p + S(1))) - d*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*e*(p + S(1)))
    rule245 = ReplacementRule(pattern245, replacement245)
    rubi.add(rule245.pattern, label = 245)

    pattern246 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons226, cons256, cons244, cons137, cons166, cons246)
    def replacement246(m, p, e, x, b, d, c, a):
        return -e**S(2)*(m + p)*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1)))
    rule246 = ReplacementRule(pattern246, replacement246)
    rubi.add(rule246.pattern, label = 246)

    pattern247 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons257, cons244, cons137, cons166, cons246)
    def replacement247(m, p, x, e, d, c, a):
        return -e**S(2)*(m + p)*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(p + S(1)))
    rule247 = ReplacementRule(pattern247, replacement247)
    rubi.add(rule247.pattern, label = 247)

    pattern248 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons226, cons256, cons31, cons266, cons238, cons246)
    def replacement248(m, p, e, x, b, d, c, a):
        return e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + (m + p)*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*(m + S(2)*p + S(1)))
    rule248 = ReplacementRule(pattern248, replacement248)
    rubi.add(rule248.pattern, label = 248)

    pattern249 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons5, cons257, cons31, cons266, cons238, cons246)
    def replacement249(m, p, x, e, d, c, a):
        return S(2)*d*(m + p)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1)), x)/(m + S(2)*p + S(1)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1)))
    rule249 = ReplacementRule(pattern249, replacement249)
    rubi.add(rule249.pattern, label = 249)

    pattern250 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons226, cons256, cons31, cons267, cons253, cons246)
    def replacement250(m, p, e, x, b, d, c, a):
        return c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((-b*e + S(2)*c*d)*(m + p + S(1))) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1)))
    rule250 = ReplacementRule(pattern250, replacement250)
    rubi.add(rule250.pattern, label = 250)

    pattern251 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons5, cons257, cons31, cons267, cons253, cons246)
    def replacement251(m, p, x, e, d, c, a):
        return (m + S(2)*p + S(2))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*d*(m + p + S(1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(m + p + S(1)))
    rule251 = ReplacementRule(pattern251, replacement251)
    rubi.add(rule251.pattern, label = 251)

    pattern252 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons3, cons7, cons48, cons21, cons147)
    def replacement252(m, p, x, e, b, c):
        return x**(-m - p)*(e*x)**m*(b + c*x)**(-p)*(b*x + c*x**S(2))**p*Int(x**(m + p)*(b + c*x)**p, x)
    rule252 = ReplacementRule(pattern252, replacement252)
    rubi.add(rule252.pattern, label = 252)

    pattern253 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons5, cons257, cons147, cons43, cons268)
    def replacement253(m, p, x, e, d, c, a):
        return Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x)
    rule253 = ReplacementRule(pattern253, replacement253)
    rubi.add(rule253.pattern, label = 253)

    pattern254 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons226, cons256, cons147)
    def replacement254(m, p, e, x, b, d, c, a):
        return (d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x)
    rule254 = ReplacementRule(pattern254, replacement254)
    rubi.add(rule254.pattern, label = 254)

    pattern255 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons257, cons147)
    def replacement255(m, p, x, e, d, c, a):
        return (a + c*x**S(2))**FracPart(p)*(d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x)
    rule255 = ReplacementRule(pattern255, replacement255)
    rubi.add(rule255.pattern, label = 255)

    pattern256 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons47)
    def replacement256(e, x, b, d, c, a):
        return b**S(2)*Int((d + e*x)/(a + b*x + c*x**S(2)), x)/(d**S(2)*(-S(4)*a*c + b**S(2))) - S(4)*b*c*Int(S(1)/(b + S(2)*c*x), x)/(d*(-S(4)*a*c + b**S(2)))
    rule256 = ReplacementRule(pattern256, replacement256)
    rubi.add(rule256.pattern, label = 256)

    pattern257 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons47, cons242, cons54)
    def replacement257(m, p, e, x, b, d, c, a):
        return S(2)*c*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule257 = ReplacementRule(pattern257, replacement257)
    rubi.add(rule257.pattern, label = 257)

    pattern258 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons226, cons47, cons128, cons269)
    def replacement258(m, p, e, x, b, d, c, a):
        return Int(ExpandIntegrand((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x)
    rule258 = ReplacementRule(pattern258, replacement258)
    rubi.add(rule258.pattern, label = 258)

    pattern259 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons47, cons270, cons244, cons163, cons94, cons271, cons246)
    def replacement259(m, p, e, x, b, d, c, a):
        return -b*p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(d*e*(m + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1)))
    rule259 = ReplacementRule(pattern259, replacement259)
    rubi.add(rule259.pattern, label = 259)

    pattern260 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons226, cons47, cons270, cons13, cons163, cons272, cons273, cons31, cons246)
    def replacement260(m, p, e, x, b, d, c, a):
        return (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - d*p*(-S(4)*a*c + b**S(2))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(b*e*(m + S(2)*p + S(1)))
    rule260 = ReplacementRule(pattern260, replacement260)
    rubi.add(rule260.pattern, label = 260)

    pattern261 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons226, cons47, cons270, cons244, cons137, cons166, cons246)
    def replacement261(m, p, e, x, b, d, c, a):
        return -d*e*(m + S(-1))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(b*(p + S(1))) + d*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1)))
    rule261 = ReplacementRule(pattern261, replacement261)
    rubi.add(rule261.pattern, label = 261)

    pattern262 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons226, cons47, cons270, cons13, cons137, cons274, cons31, cons246)
    def replacement262(m, p, e, x, b, d, c, a):
        return -S(2)*c*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + S(2)*c*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule262 = ReplacementRule(pattern262, replacement262)
    rubi.add(rule262.pattern, label = 262)

    pattern263 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons47)
    def replacement263(e, x, b, d, c, a):
        return S(4)*c*Subst(Int(S(1)/(-S(4)*a*c*e + b**S(2)*e + S(4)*c*e*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2)))
    rule263 = ReplacementRule(pattern263, replacement263)
    rubi.add(rule263.pattern, label = 263)

    pattern264 = Pattern(Integral(S(1)/(sqrt(d_ + x_*WC('e', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons47, cons275)
    def replacement264(e, x, b, d, c, a):
        return S(4)*sqrt(-c/(-S(4)*a*c + b**S(2)))*Subst(Int(S(1)/sqrt(Simp(-b**S(2)*x**S(4)/(d**S(2)*(-S(4)*a*c + b**S(2))) + S(1), x)), x), x, sqrt(d + e*x))/e
    rule264 = ReplacementRule(pattern264, replacement264)
    rubi.add(rule264.pattern, label = 264)

    pattern265 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons47, cons275)
    def replacement265(e, x, b, d, c, a):
        return S(4)*sqrt(-c/(-S(4)*a*c + b**S(2)))*Subst(Int(x**S(2)/sqrt(Simp(-b**S(2)*x**S(4)/(d**S(2)*(-S(4)*a*c + b**S(2))) + S(1), x)), x), x, sqrt(d + e*x))/e
    rule265 = ReplacementRule(pattern265, replacement265)
    rubi.add(rule265.pattern, label = 265)

    pattern266 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons47, cons276)
    def replacement266(m, e, x, b, d, c, a):
        return sqrt(-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))*Int((d + e*x)**m/sqrt(-a*c/(-S(4)*a*c + b**S(2)) - b*c*x/(-S(4)*a*c + b**S(2)) - c**S(2)*x**S(2)/(-S(4)*a*c + b**S(2))), x)/sqrt(a + b*x + c*x**S(2))
    rule266 = ReplacementRule(pattern266, replacement266)
    rubi.add(rule266.pattern, label = 266)

    pattern267 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons226, cons47, cons270, cons31, cons166, cons238, cons277)
    def replacement267(m, p, e, x, b, d, c, a):
        return S(2)*d*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(m + S(2)*p + S(1))) + d**S(2)*(m + S(-1))*(-S(4)*a*c + b**S(2))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**p, x)/(b**S(2)*(m + S(2)*p + S(1)))
    rule267 = ReplacementRule(pattern267, replacement267)
    rubi.add(rule267.pattern, label = 267)

    pattern268 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons226, cons47, cons270, cons31, cons94, cons278)
    def replacement268(m, p, e, x, b, d, c, a):
        return b**S(2)*(m + S(2)*p + S(3))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/(d**S(2)*(m + S(1))*(-S(4)*a*c + b**S(2))) - S(2)*b*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(d*(m + S(1))*(-S(4)*a*c + b**S(2)))
    rule268 = ReplacementRule(pattern268, replacement268)
    rubi.add(rule268.pattern, label = 268)

    pattern269 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons47)
    def replacement269(m, p, e, x, b, d, c, a):
        return Subst(Int(x**m*(a - b**S(2)/(S(4)*c) + c*x**S(2)/e**S(2))**p, x), x, d + e*x)/e
    rule269 = ReplacementRule(pattern269, replacement269)
    rubi.add(rule269.pattern, label = 269)

    pattern270 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons226, cons279, cons239, cons128)
    def replacement270(m, p, e, x, b, d, c, a):
        return Int(ExpandIntegrand((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x)
    rule270 = ReplacementRule(pattern270, replacement270)
    rubi.add(rule270.pattern, label = 270)

    pattern271 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons280, cons128, cons281)
    def replacement271(m, p, e, x, d, c, a):
        return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m, x), x)
    rule271 = ReplacementRule(pattern271, replacement271)
    rubi.add(rule271.pattern, label = 271)

    def With272(x, e, b, d, c, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (c*d - e*(b/S(2) - q/S(2)))*Int(S(1)/(b/S(2) + c*x - q/S(2)), x)/q - (c*d - e*(b/S(2) + q/S(2)))*Int(S(1)/(b/S(2) + c*x + q/S(2)), x)/q
    pattern272 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239, cons282 )
    rule272 = ReplacementRule(pattern272, With272)
    rubi.add(rule272.pattern, label = 272)

    def With273(x, e, d, c, a):
        q = Rt(-a*c, S(2))
        return (-c*d/(S(2)*q) + e/S(2))*Int(S(1)/(c*x + q), x) + (c*d/(S(2)*q) + e/S(2))*Int(S(1)/(c*x - q), x)
    pattern273 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons280, cons283 )
    rule273 = ReplacementRule(pattern273, With273)
    rubi.add(rule273.pattern, label = 273)

    pattern274 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239, cons284)
    def replacement274(x, e, b, d, c, a):
        return e*Int((b + S(2)*c*x)/(a + b*x + c*x**S(2)), x)/(S(2)*c) + (-b*e + S(2)*c*d)*Int(S(1)/(a + b*x + c*x**S(2)), x)/(S(2)*c)
    rule274 = ReplacementRule(pattern274, replacement274)
    rubi.add(rule274.pattern, label = 274)

    pattern275 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons280, cons285)
    def replacement275(x, e, d, c, a):
        return d*Int(S(1)/(a + c*x**S(2)), x) + e*Int(x/(a + c*x**S(2)), x)
    rule275 = ReplacementRule(pattern275, replacement275)
    rubi.add(rule275.pattern, label = 275)

    pattern276 = Pattern(Integral(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239)
    def replacement276(e, x, b, d, c, a):
        return S(2)*e*Subst(Int(x**S(2)/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x))
    rule276 = ReplacementRule(pattern276, replacement276)
    rubi.add(rule276.pattern, label = 276)

    pattern277 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons280)
    def replacement277(x, e, d, c, a):
        return S(2)*e*Subst(Int(x**S(2)/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x))
    rule277 = ReplacementRule(pattern277, replacement277)
    rubi.add(rule277.pattern, label = 277)

    pattern278 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239, cons17, cons166, cons286)
    def replacement278(m, e, x, b, d, c, a):
        return Int(PolynomialDivide((d + e*x)**m, a + b*x + c*x**S(2), x), x)
    rule278 = ReplacementRule(pattern278, replacement278)
    rubi.add(rule278.pattern, label = 278)

    pattern279 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons280, cons17, cons166, cons286)
    def replacement279(m, x, e, d, c, a):
        return Int(PolynomialDivide((d + e*x)**m, a + c*x**S(2), x), x)
    rule279 = ReplacementRule(pattern279, replacement279)
    rubi.add(rule279.pattern, label = 279)

    pattern280 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239, cons31, cons166)
    def replacement280(m, e, x, b, d, c, a):
        return e*(d + e*x)**(m + S(-1))/(c*(m + S(-1))) + Int((d + e*x)**(m + S(-2))*Simp(-a*e**S(2) + c*d**S(2) + e*x*(-b*e + S(2)*c*d), x)/(a + b*x + c*x**S(2)), x)/c
    rule280 = ReplacementRule(pattern280, replacement280)
    rubi.add(rule280.pattern, label = 280)

    pattern281 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons280, cons31, cons166)
    def replacement281(m, x, e, d, c, a):
        return e*(d + e*x)**(m + S(-1))/(c*(m + S(-1))) + Int((d + e*x)**(m + S(-2))*Simp(-a*e**S(2) + c*d**S(2) + S(2)*c*d*e*x, x)/(a + c*x**S(2)), x)/c
    rule281 = ReplacementRule(pattern281, replacement281)
    rubi.add(rule281.pattern, label = 281)

    pattern282 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239)
    def replacement282(e, x, b, d, c, a):
        return e**S(2)*Int(S(1)/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int((-b*e + c*d - c*e*x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2))
    rule282 = ReplacementRule(pattern282, replacement282)
    rubi.add(rule282.pattern, label = 282)

    pattern283 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons280)
    def replacement283(x, e, d, c, a):
        return e**S(2)*Int(S(1)/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int((c*d - c*e*x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2))
    rule283 = ReplacementRule(pattern283, replacement283)
    rubi.add(rule283.pattern, label = 283)

    pattern284 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239)
    def replacement284(e, x, b, d, c, a):
        return S(2)*e*Subst(Int(S(1)/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x))
    rule284 = ReplacementRule(pattern284, replacement284)
    rubi.add(rule284.pattern, label = 284)

    pattern285 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons280)
    def replacement285(x, e, d, c, a):
        return S(2)*e*Subst(Int(S(1)/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x))
    rule285 = ReplacementRule(pattern285, replacement285)
    rubi.add(rule285.pattern, label = 285)

    pattern286 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons226, cons279, cons239, cons31, cons94)
    def replacement286(m, e, x, b, d, c, a):
        return e*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(-b*e + c*d - c*e*x, x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2))
    rule286 = ReplacementRule(pattern286, replacement286)
    rubi.add(rule286.pattern, label = 286)

    pattern287 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons21, cons280, cons31, cons94)
    def replacement287(m, x, e, d, c, a):
        return c*Int((d - e*x)*(d + e*x)**(m + S(1))/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)) + e*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2)))
    rule287 = ReplacementRule(pattern287, replacement287)
    rubi.add(rule287.pattern, label = 287)

    pattern288 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons226, cons279, cons239, cons18)
    def replacement288(m, e, x, b, d, c, a):
        return Int(ExpandIntegrand((d + e*x)**m, S(1)/(a + b*x + c*x**S(2)), x), x)
    rule288 = ReplacementRule(pattern288, replacement288)
    rubi.add(rule288.pattern, label = 288)

    pattern289 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons21, cons280, cons18)
    def replacement289(m, x, e, d, c, a):
        return Int(ExpandIntegrand((d + e*x)**m, S(1)/(a + c*x**S(2)), x), x)
    rule289 = ReplacementRule(pattern289, replacement289)
    rubi.add(rule289.pattern, label = 289)

    pattern290 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239)
    def replacement290(e, x, b, d, c, a):
        return -S(2)*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((-S(4)*a*c + b**S(2))*sqrt(a + b*x + c*x**S(2)))
    rule290 = ReplacementRule(pattern290, replacement290)
    rubi.add(rule290.pattern, label = 290)

    pattern291 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2), x_), cons2, cons7, cons27, cons48, cons280)
    def replacement291(x, e, d, c, a):
        return (-a*e + c*d*x)/(a*c*sqrt(a + c*x**S(2)))
    rule291 = ReplacementRule(pattern291, replacement291)
    rubi.add(rule291.pattern, label = 291)

    pattern292 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239, cons13, cons137, cons230)
    def replacement292(p, e, x, b, d, c, a):
        return -(S(2)*p + S(3))*(-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule292 = ReplacementRule(pattern292, replacement292)
    rubi.add(rule292.pattern, label = 292)

    pattern293 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons280, cons13, cons137, cons230)
    def replacement293(p, x, e, d, c, a):
        return d*(S(2)*p + S(3))*Int((a + c*x**S(2))**(p + S(1)), x)/(S(2)*a*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1)))
    rule293 = ReplacementRule(pattern293, replacement293)
    rubi.add(rule293.pattern, label = 293)

    pattern294 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons226, cons279, cons239, cons287)
    def replacement294(p, e, x, b, d, c, a):
        return e*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))) + (-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c)
    rule294 = ReplacementRule(pattern294, replacement294)
    rubi.add(rule294.pattern, label = 294)

    pattern295 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons5, cons280, cons287)
    def replacement295(p, x, e, d, c, a):
        return d*Int((a + c*x**S(2))**p, x) + e*(a + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1)))
    rule295 = ReplacementRule(pattern295, replacement295)
    rubi.add(rule295.pattern, label = 295)

    pattern296 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons288, cons289, cons290, cons147)
    def replacement296(m, p, x, e, b, d, c, a):
        return (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m - p)*(a*d + c*e*x**S(3))**p, x)
    rule296 = ReplacementRule(pattern296, replacement296)
    rubi.add(rule296.pattern, label = 296)

    pattern297 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons3, cons7, cons27, cons48, cons291, cons239, cons31, cons292, cons293, cons294)
    def replacement297(m, x, e, b, d, c):
        return Int((d + e*x)**m/(sqrt(b*x)*sqrt(S(1) + c*x/b)), x)
    rule297 = ReplacementRule(pattern297, replacement297)
    rubi.add(rule297.pattern, label = 297)

    pattern298 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons3, cons7, cons27, cons48, cons291, cons239, cons31, cons292)
    def replacement298(m, x, e, b, d, c):
        return sqrt(x)*sqrt(b + c*x)*Int((d + e*x)**m/(sqrt(x)*sqrt(b + c*x)), x)/sqrt(b*x + c*x**S(2))
    rule298 = ReplacementRule(pattern298, replacement298)
    rubi.add(rule298.pattern, label = 298)

    pattern299 = Pattern(Integral(x_**m_/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons226, cons295)
    def replacement299(m, x, b, c, a):
        return S(2)*Subst(Int(x**(S(2)*m + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x, sqrt(x))
    rule299 = ReplacementRule(pattern299, replacement299)
    rubi.add(rule299.pattern, label = 299)

    pattern300 = Pattern(Integral((e_*x_)**m_/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons48, cons226, cons295)
    def replacement300(m, x, e, b, c, a):
        return x**(-m)*(e*x)**m*Int(x**m/sqrt(a + b*x + c*x**S(2)), x)
    rule300 = ReplacementRule(pattern300, replacement300)
    rubi.add(rule300.pattern, label = 300)

    pattern301 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239, cons295)
    def replacement301(m, e, x, b, d, c, a):
        return S(2)*sqrt(-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))*(S(2)*c*(d + e*x)/(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2))))**(-m)*(d + e*x)**m*Rt(-S(4)*a*c + b**S(2), S(2))*Subst(Int((S(2)*e*x**S(2)*Rt(-S(4)*a*c + b**S(2), S(2))/(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**m/sqrt(-x**S(2) + S(1)), x), x, sqrt(S(2))*sqrt((b + S(2)*c*x + Rt(-S(4)*a*c + b**S(2), S(2)))/Rt(-S(4)*a*c + b**S(2), S(2)))/S(2))/(c*sqrt(a + b*x + c*x**S(2)))
    rule301 = ReplacementRule(pattern301, replacement301)
    rubi.add(rule301.pattern, label = 301)

    pattern302 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/sqrt(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons280, cons295)
    def replacement302(m, x, e, d, c, a):
        return S(2)*a*(c*(d + e*x)/(-a*e*Rt(-c/a, S(2)) + c*d))**(-m)*sqrt(S(1) + c*x**S(2)/a)*(d + e*x)**m*Rt(-c/a, S(2))*Subst(Int((S(2)*a*e*x**S(2)*Rt(-c/a, S(2))/(-a*e*Rt(-c/a, S(2)) + c*d) + S(1))**m/sqrt(-x**S(2) + S(1)), x), x, sqrt(-x*Rt(-c/a, S(2))/S(2) + S(1)/2))/(c*sqrt(a + c*x**S(2)))
    rule302 = ReplacementRule(pattern302, replacement302)
    rubi.add(rule302.pattern, label = 302)

    pattern303 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239, cons244, cons296, cons163)
    def replacement303(m, p, e, x, b, d, c, a):
        return p*(-S(4)*a*c + b**S(2))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) - (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/(S(2)*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2)))
    rule303 = ReplacementRule(pattern303, replacement303)
    rubi.add(rule303.pattern, label = 303)

    pattern304 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons280, cons244, cons296, cons163)
    def replacement304(m, p, x, e, d, c, a):
        return -S(2)*a*c*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2)), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))) - (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-S(2)*a*e + S(2)*c*d*x)/(S(2)*(m + S(1))*(a*e**S(2) + c*d**S(2)))
    rule304 = ReplacementRule(pattern304, replacement304)
    rubi.add(rule304.pattern, label = 304)

    pattern305 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239, cons244, cons296, cons137)
    def replacement305(m, p, e, x, b, d, c, a):
        return (d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))) - S(2)*(S(2)*p + S(3))*(a*e**S(2) - b*d*e + c*d**S(2))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule305 = ReplacementRule(pattern305, replacement305)
    rubi.add(rule305.pattern, label = 305)

    pattern306 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons280, cons244, cons296, cons137)
    def replacement306(m, p, x, e, d, c, a):
        return (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))) + (S(2)*p + S(3))*(a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2)), x)/(S(2)*a*c*(p + S(1)))
    rule306 = ReplacementRule(pattern306, replacement306)
    rubi.add(rule306.pattern, label = 306)

    pattern307 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons239)
    def replacement307(e, x, b, d, c, a):
        return -S(2)*Subst(Int(S(1)/(S(4)*a*e**S(2) - S(4)*b*d*e + S(4)*c*d**S(2) - x**S(2)), x), x, (S(2)*a*e - b*d - x*(-b*e + S(2)*c*d))/sqrt(a + b*x + c*x**S(2)))
    rule307 = ReplacementRule(pattern307, replacement307)
    rubi.add(rule307.pattern, label = 307)

    pattern308 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons297)
    def replacement308(x, e, d, c, a):
        return -Subst(Int(S(1)/(a*e**S(2) + c*d**S(2) - x**S(2)), x), x, (a*e - c*d*x)/sqrt(a + c*x**S(2)))
    rule308 = ReplacementRule(pattern308, replacement308)
    rubi.add(rule308.pattern, label = 308)

    pattern309 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons279, cons239, cons147, cons240)
    def replacement309(m, p, e, x, b, d, c, a):
        return -((b + S(2)*c*x + Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d + e*Rt(-S(4)*a*c + b**S(2), S(2)))/((b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2)))))**(-p)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*Hypergeometric2F1(m + S(1), -p, m + S(2), -S(4)*c*(d + e*x)*Rt(-S(4)*a*c + b**S(2), S(2))/((b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2)))))/((m + S(1))*(-b*e + S(2)*c*d + e*Rt(-S(4)*a*c + b**S(2), S(2))))
    rule309 = ReplacementRule(pattern309, replacement309)
    rubi.add(rule309.pattern, label = 309)

    pattern310 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons5, cons280, cons147, cons240)
    def replacement310(m, p, x, e, d, c, a):
        return ((c*d + e*Rt(-a*c, S(2)))*(c*x + Rt(-a*c, S(2)))/((c*d - e*Rt(-a*c, S(2)))*(c*x - Rt(-a*c, S(2)))))**(-p)*(a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-c*x + Rt(-a*c, S(2)))*Hypergeometric2F1(m + S(1), -p, m + S(2), S(2)*c*(d + e*x)*Rt(-a*c, S(2))/((c*d - e*Rt(-a*c, S(2)))*(-c*x + Rt(-a*c, S(2)))))/((m + S(1))*(c*d + e*Rt(-a*c, S(2))))
    rule310 = ReplacementRule(pattern310, replacement310)
    rubi.add(rule310.pattern, label = 310)

    pattern311 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons279, cons239, cons242, cons13, cons137)
    def replacement311(m, p, e, x, b, d, c, a):
        return m*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule311 = ReplacementRule(pattern311, replacement311)
    rubi.add(rule311.pattern, label = 311)

    pattern312 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons21, cons5, cons280, cons242, cons13, cons137)
    def replacement312(m, p, x, e, d, c, a):
        return -d*m*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*(p + S(1))) - x*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*(p + S(1)))
    rule312 = ReplacementRule(pattern312, replacement312)
    rubi.add(rule312.pattern, label = 312)

    pattern313 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons279, cons239, cons242)
    def replacement313(m, p, e, x, b, d, c, a):
        return e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + (-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*a*e**S(2) - S(2)*b*d*e + S(2)*c*d**S(2))
    rule313 = ReplacementRule(pattern313, replacement313)
    rubi.add(rule313.pattern, label = 313)

    pattern314 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons21, cons5, cons280, cons242)
    def replacement314(m, p, x, e, d, c, a):
        return c*d*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(a*e**S(2) + c*d**S(2)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2)))
    rule314 = ReplacementRule(pattern314, replacement314)
    rubi.add(rule314.pattern, label = 314)

    pattern315 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons226, cons279, cons239, cons13, cons163, cons298, cons66, cons299, cons300)
    def replacement315(m, p, e, x, b, d, c, a):
        return -p*Int((b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e*(m + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1)))
    rule315 = ReplacementRule(pattern315, replacement315)
    rubi.add(rule315.pattern, label = 315)

    pattern316 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons21, cons280, cons13, cons163, cons298, cons66, cons299, cons301)
    def replacement316(m, p, x, e, d, c, a):
        return -S(2)*c*p*Int(x*(a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1)), x)/(e*(m + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(1)))
    rule316 = ReplacementRule(pattern316, replacement316)
    rubi.add(rule316.pattern, label = 316)

    pattern317 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons226, cons279, cons239, cons13, cons163, cons238, cons302, cons303, cons300)
    def replacement317(m, p, e, x, b, d, c, a):
        return -p*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d), x), x)/(e*(m + S(2)*p + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1)))
    rule317 = ReplacementRule(pattern317, replacement317)
    rubi.add(rule317.pattern, label = 317)

    pattern318 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons21, cons280, cons13, cons163, cons238, cons302, cons303, cons301)
    def replacement318(m, p, x, e, d, c, a):
        return S(2)*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**m*Simp(a*e - c*d*x, x), x)/(e*(m + S(2)*p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(2)*p + S(1)))
    rule318 = ReplacementRule(pattern318, replacement318)
    rubi.add(rule318.pattern, label = 318)

    pattern319 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239, cons244, cons137, cons168, cons304, cons300)
    def replacement319(m, p, e, x, b, d, c, a):
        return (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(b*e*m + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*e*x*(m + S(2)*p + S(3))), x)/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule319 = ReplacementRule(pattern319, replacement319)
    rubi.add(rule319.pattern, label = 319)

    pattern320 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons280, cons244, cons137, cons168, cons304, cons301)
    def replacement320(m, p, x, e, d, c, a):
        return -x*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(d*(S(2)*p + S(3)) + e*x*(m + S(2)*p + S(3))), x)/(S(2)*a*(p + S(1)))
    rule320 = ReplacementRule(pattern320, replacement320)
    rubi.add(rule320.pattern, label = 320)

    pattern321 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons239, cons244, cons137, cons166, cons300)
    def replacement321(m, p, e, x, b, d, c, a):
        return (d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))) + Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-S(2)*c*d**S(2)*(S(2)*p + S(3)) + e*x*(b*e - S(2)*c*d)*(m + S(2)*p + S(2)) + e*(S(2)*a*e*(m + S(-1)) + b*d*(-m + S(2)*p + S(4))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule321 = ReplacementRule(pattern321, replacement321)
    rubi.add(rule321.pattern, label = 321)

    pattern322 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons280, cons244, cons137, cons166, cons301)
    def replacement322(m, p, x, e, d, c, a):
        return (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*Simp(a*e**S(2)*(m + S(-1)) - c*d**S(2)*(S(2)*p + S(3)) - c*d*e*x*(m + S(2)*p + S(2)), x), x)/(S(2)*a*c*(p + S(1)))
    rule322 = ReplacementRule(pattern322, replacement322)
    rubi.add(rule322.pattern, label = 322)

    pattern323 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons226, cons279, cons239, cons13, cons137, cons300)
    def replacement323(m, p, e, x, b, d, c, a):
        return (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*e - b**S(2)*e + b*c*d + c*x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-S(2)*a*c*e**S(2)*(m + S(2)*p + S(3)) + b**S(2)*e**S(2)*(m + p + S(2)) + b*c*d*e*(-m + S(2)*p + S(2)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(3)) - c*e*x*(-b*e + S(2)*c*d)*(m + S(2)*p + S(4)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2)))
    rule323 = ReplacementRule(pattern323, replacement323)
    rubi.add(rule323.pattern, label = 323)

    pattern324 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons21, cons280, cons13, cons137, cons301)
    def replacement324(m, p, x, e, d, c, a):
        return -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(a*e + c*d*x)/(S(2)*a*(p + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*Simp(a*e**S(2)*(m + S(2)*p + S(3)) + c*d**S(2)*(S(2)*p + S(3)) + c*d*e*x*(m + S(2)*p + S(4)), x), x)/(S(2)*a*(p + S(1))*(a*e**S(2) + c*d**S(2)))
    rule324 = ReplacementRule(pattern324, replacement324)
    rubi.add(rule324.pattern, label = 324)

    pattern325 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons279, cons239, cons305, cons238, cons300)
    def replacement325(m, p, e, x, b, d, c, a):
        return e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**p*Simp(c*d**S(2)*(m + S(2)*p + S(1)) + e*x*(m + p)*(-b*e + S(2)*c*d) - e*(a*e*(m + S(-1)) + b*d*(p + S(1))), x), x)/(c*(m + S(2)*p + S(1)))
    rule325 = ReplacementRule(pattern325, replacement325)
    rubi.add(rule325.pattern, label = 325)

    pattern326 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons21, cons5, cons280, cons305, cons238, cons301)
    def replacement326(m, p, x, e, d, c, a):
        return e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-2))*Simp(-a*e**S(2)*(m + S(-1)) + c*d**S(2)*(m + S(2)*p + S(1)) + S(2)*c*d*e*x*(m + p), x), x)/(c*(m + S(2)*p + S(1)))
    rule326 = ReplacementRule(pattern326, replacement326)
    rubi.add(rule326.pattern, label = 326)

    pattern327 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons279, cons239, cons306)
    def replacement327(m, p, e, x, b, d, c, a):
        return e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*e*(m + p + S(2)) + c*d*(m + S(1)) - c*e*x*(m + S(2)*p + S(3)), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2)))
    rule327 = ReplacementRule(pattern327, replacement327)
    rubi.add(rule327.pattern, label = 327)

    pattern328 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons7, cons27, cons48, cons21, cons5, cons280, cons307)
    def replacement328(m, p, x, e, d, c, a):
        return c*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(d*(m + S(1)) - e*x*(m + S(2)*p + S(3)), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2)))
    rule328 = ReplacementRule(pattern328, replacement328)
    rubi.add(rule328.pattern, label = 328)

    def With329(x, e, b, d, c, a):
        q = Rt(S(3)*c*e**S(2)*(-b*e + S(2)*c*d), S(3))
        return -sqrt(S(3))*c*e*ArcTan(sqrt(S(3))/S(3) + S(2)*sqrt(S(3))*(-b*e + c*d - c*e*x)/(S(3)*q*(a + b*x + c*x**S(2))**(S(1)/3)))/q**S(2) - S(3)*c*e*log(d + e*x)/(S(2)*q**S(2)) + S(3)*c*e*log(-b*e + c*d - c*e*x - q*(a + b*x + c*x**S(2))**(S(1)/3))/(S(2)*q**S(2))
    pattern329 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons2, cons3, cons7, cons27, cons48, cons239, cons308, cons309 )
    rule329 = ReplacementRule(pattern329, With329)
    rubi.add(rule329.pattern, label = 329)

    def With330(x, e, d, c, a):
        q = Rt(S(6)*c**S(2)*e**S(2)/d**S(2), S(3))
        return -sqrt(S(3))*c*e*ArcTan(S(2)*sqrt(S(3))*c*(d - e*x)/(S(3)*d*q*(a + c*x**S(2))**(S(1)/3)) + sqrt(S(3))/S(3))/(d**S(2)*q**S(2)) - S(3)*c*e*log(d + e*x)/(S(2)*d**S(2)*q**S(2)) + S(3)*c*e*log(c*d - c*e*x - d*q*(a + c*x**S(2))**(S(1)/3))/(S(2)*d**S(2)*q**S(2))
    pattern330 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons310 )
    rule330 = ReplacementRule(pattern330, With330)
    rubi.add(rule330.pattern, label = 330)

    def With331(x, e, b, d, c, a):
        q = Rt(-S(3)*c*e**S(2)*(-b*e + S(2)*c*d), S(3))
        return -sqrt(S(3))*c*e*ArcTan(sqrt(S(3))/S(3) - S(2)*sqrt(S(3))*(-b*e + c*d - c*e*x)/(S(3)*q*(a + b*x + c*x**S(2))**(S(1)/3)))/q**S(2) - S(3)*c*e*log(d + e*x)/(S(2)*q**S(2)) + S(3)*c*e*log(-b*e + c*d - c*e*x + q*(a + b*x + c*x**S(2))**(S(1)/3))/(S(2)*q**S(2))
    pattern331 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons2, cons3, cons7, cons27, cons48, cons239, cons308, cons311 )
    rule331 = ReplacementRule(pattern331, With331)
    rubi.add(rule331.pattern, label = 331)

    def With332(x, e, b, d, c, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (b + S(2)*c*x - q)**(S(1)/3)*(b + S(2)*c*x + q)**(S(1)/3)*Int(S(1)/((d + e*x)*(b + S(2)*c*x - q)**(S(1)/3)*(b + S(2)*c*x + q)**(S(1)/3)), x)/(a + b*x + c*x**S(2))**(S(1)/3)
    pattern332 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons312 )
    rule332 = ReplacementRule(pattern332, With332)
    rubi.add(rule332.pattern, label = 332)

    pattern333 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/4)*(d_ + x_*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons280)
    def replacement333(x, e, d, c, a):
        return d*Int(S(1)/((a + c*x**S(2))**(S(1)/4)*(d**S(2) - e**S(2)*x**S(2))), x) - e*Int(x/((a + c*x**S(2))**(S(1)/4)*(d**S(2) - e**S(2)*x**S(2))), x)
    rule333 = ReplacementRule(pattern333, replacement333)
    rubi.add(rule333.pattern, label = 333)

    pattern334 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/4)*(d_ + x_*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons280)
    def replacement334(x, e, d, c, a):
        return d*Int(S(1)/((a + c*x**S(2))**(S(3)/4)*(d**S(2) - e**S(2)*x**S(2))), x) - e*Int(x/((a + c*x**S(2))**(S(3)/4)*(d**S(2) - e**S(2)*x**S(2))), x)
    rule334 = ReplacementRule(pattern334, replacement334)
    rubi.add(rule334.pattern, label = 334)

    pattern335 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons232, cons229)
    def replacement335(p, e, x, b, d, c, a):
        return (-S(4)*c/(-S(4)*a*c + b**S(2)))**(-p)*Subst(Int(Simp(-x**S(2)/(-S(4)*a*c + b**S(2)) + S(1), x)**p/Simp(-b*e + S(2)*c*d + e*x, x), x), x, b + S(2)*c*x)
    rule335 = ReplacementRule(pattern335, replacement335)
    rubi.add(rule335.pattern, label = 335)

    pattern336 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons313, cons229)
    def replacement336(p, e, x, b, d, c, a):
        return (-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))**(-p)*(a + b*x + c*x**S(2))**p*Int((-a*c/(-S(4)*a*c + b**S(2)) - b*c*x/(-S(4)*a*c + b**S(2)) - c**S(2)*x**S(2)/(-S(4)*a*c + b**S(2)))**p/(d + e*x), x)
    rule336 = ReplacementRule(pattern336, replacement336)
    rubi.add(rule336.pattern, label = 336)

    pattern337 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons5, cons280, cons147, cons43, cons293)
    def replacement337(m, p, x, e, d, c, a):
        return Int((d + e*x)**m*(-x*Rt(-c, S(2)) + Rt(a, S(2)))**p*(x*Rt(-c, S(2)) + Rt(a, S(2)))**p, x)
    rule337 = ReplacementRule(pattern337, replacement337)
    rubi.add(rule337.pattern, label = 337)

    def With338(m, p, e, x, b, d, c, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -(e*(b + S(2)*c*x - q)/(S(2)*c*(d + e*x)))**(-p)*(e*(b + S(2)*c*x + q)/(S(2)*c*(d + e*x)))**(-p)*(a + b*x + c*x**S(2))**p*(S(1)/(d + e*x))**(S(2)*p)*Subst(Int(x**(-m - S(2)*p + S(-2))*Simp(-x*(d - e*(b - q)/(S(2)*c)) + S(1), x)**p*Simp(-x*(d - e*(b + q)/(S(2)*c)) + S(1), x)**p, x), x, S(1)/(d + e*x))/e
    pattern338 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons226, cons279, cons239, cons147, cons84 )
    rule338 = ReplacementRule(pattern338, With338)
    rubi.add(rule338.pattern, label = 338)

    def With339(m, p, x, e, d, c, a):
        q = Rt(-a*c, S(2))
        return -(e*(c*x + q)/(c*(d + e*x)))**(-p)*(-e*(-c*x + q)/(c*(d + e*x)))**(-p)*(a + c*x**S(2))**p*(S(1)/(d + e*x))**(S(2)*p)*Subst(Int(x**(-m - S(2)*p + S(-2))*Simp(-x*(d - e*q/c) + S(1), x)**p*Simp(-x*(d + e*q/c) + S(1), x)**p, x), x, S(1)/(d + e*x))/e
    pattern339 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons5, cons280, cons147, cons84 )
    rule339 = ReplacementRule(pattern339, With339)
    rubi.add(rule339.pattern, label = 339)

    def With340(m, p, e, x, b, d, c, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (-(d + e*x)/(d - e*(b - q)/(S(2)*c)) + S(1))**(-p)*(-(d + e*x)/(d - e*(b + q)/(S(2)*c)) + S(1))**(-p)*(a + b*x + c*x**S(2))**p*Subst(Int(x**m*Simp(-x/(d - e*(b - q)/(S(2)*c)) + S(1), x)**p*Simp(-x/(d - e*(b + q)/(S(2)*c)) + S(1), x)**p, x), x, d + e*x)/e
    pattern340 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons226, cons279, cons239, cons147 )
    rule340 = ReplacementRule(pattern340, With340)
    rubi.add(rule340.pattern, label = 340)

    def With341(m, p, x, e, d, c, a):
        q = Rt(-a*c, S(2))
        return (a + c*x**S(2))**p*(-(d + e*x)/(d - e*q/c) + S(1))**(-p)*(-(d + e*x)/(d + e*q/c) + S(1))**(-p)*Subst(Int(x**m*Simp(-x/(d - e*q/c) + S(1), x)**p*Simp(-x/(d + e*q/c) + S(1), x)**p, x), x, d + e*x)/e
    pattern341 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons5, cons280, cons147 )
    rule341 = ReplacementRule(pattern341, With341)
    rubi.add(rule341.pattern, label = 341)

    pattern342 = Pattern(Integral((u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons68, cons69)
    def replacement342(m, p, e, x, b, d, c, u, a):
        return Subst(Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1))
    rule342 = ReplacementRule(pattern342, replacement342)
    rubi.add(rule342.pattern, label = 342)

    pattern343 = Pattern(Integral((a_ + u_**S(2)*WC('c', S(1)))**WC('p', S(1))*(u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons5, cons68, cons69)
    def replacement343(m, p, e, x, d, c, u, a):
        return Subst(Int((a + c*x**S(2))**p*(d + e*x)**m, x), x, u)/Coefficient(u, x, S(1))
    rule343 = ReplacementRule(pattern343, replacement343)
    rubi.add(rule343.pattern, label = 343)

    pattern344 = Pattern(Integral(x_**WC('n', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons5, cons85, cons314)
    def replacement344(n, p, e, x, d, c, a):
        return d*Int(x**n*(a + c*x**S(2))**p, x) + e*Int(x**(n + S(1))*(a + c*x**S(2))**p, x)
    rule344 = ReplacementRule(pattern344, replacement344)
    rubi.add(rule344.pattern, label = 344)

    pattern345 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons315, cons45, cons316)
    def replacement345(m, g, e, x, b, d, f, c, a):
        return (f + g*x)*Int((d + e*x)**m, x)/sqrt(a + b*x + c*x**S(2))
    rule345 = ReplacementRule(pattern345, replacement345)
    rubi.add(rule345.pattern, label = 345)

    pattern346 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons45, cons316, cons147, cons242)
    def replacement346(m, p, g, e, x, b, d, f, c, a):
        return -f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))*(-d*g + e*f))
    rule346 = ReplacementRule(pattern346, replacement346)
    rubi.add(rule346.pattern, label = 346)

    pattern347 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons45, cons316, cons147, cons244, cons137, cons168)
    def replacement347(m, p, g, e, x, b, d, f, c, a):
        return -e*g*m*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(S(2)*c*(p + S(1))) + g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1)))
    rule347 = ReplacementRule(pattern347, replacement347)
    rubi.add(rule347.pattern, label = 347)

    pattern348 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons315, cons45, cons316, cons147, cons13, cons137, cons317)
    def replacement348(m, p, g, e, x, b, d, f, c, a):
        return e*f*g*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/(b*(p + S(1))*(-d*g + e*f)) - f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))*(-d*g + e*f))
    rule348 = ReplacementRule(pattern348, replacement348)
    rubi.add(rule348.pattern, label = 348)

    pattern349 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons45, cons316, cons147, cons31, cons94, cons225, cons318)
    def replacement349(m, p, g, e, x, b, d, f, c, a):
        return -g*(S(2)*p + S(1))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(1))) + (d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p/(e*(m + S(1)))
    rule349 = ReplacementRule(pattern349, replacement349)
    rubi.add(rule349.pattern, label = 349)

    pattern350 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons45, cons316, cons147, cons31, cons94, cons319)
    def replacement350(m, p, g, e, x, b, d, f, c, a):
        return -g*(m + S(2)*p + S(3))*Int((d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p, x)/((m + S(1))*(-d*g + e*f)) + S(2)*f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(m + S(1))*(-d*g + e*f))
    rule350 = ReplacementRule(pattern350, replacement350)
    rubi.add(rule350.pattern, label = 350)

    pattern351 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons45, cons316, cons147, cons62, cons319, cons320)
    def replacement351(m, p, g, e, x, b, d, f, c, a):
        return -b*m*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*f*(m + S(2)*p + S(2))) + g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2)))
    rule351 = ReplacementRule(pattern351, replacement351)
    rubi.add(rule351.pattern, label = 351)

    pattern352 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons45, cons316, cons147, cons319)
    def replacement352(m, p, g, e, x, b, d, f, c, a):
        return (d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(2))) + (S(2)*p + S(1))*(-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(2)*p + S(2)))
    rule352 = ReplacementRule(pattern352, replacement352)
    rubi.add(rule352.pattern, label = 352)

    pattern353 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons45, cons321, cons147, cons239, cons13, cons322)
    def replacement353(p, g, e, x, b, d, f, c, a):
        return (-b*g + S(2)*c*f)*Int((a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d) - (-d*g + e*f)*Int((b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(d + e*x), x)/(-b*e + S(2)*c*d)
    rule353 = ReplacementRule(pattern353, replacement353)
    rubi.add(rule353.pattern, label = 353)

    pattern354 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons45, cons321, cons147, cons323)
    def replacement354(m, p, g, e, x, b, d, f, c, a):
        return g*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/e
    rule354 = ReplacementRule(pattern354, replacement354)
    rubi.add(rule354.pattern, label = 354)

    pattern355 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons45, cons321, cons147, cons239, cons242)
    def replacement355(m, p, g, e, x, b, d, f, c, a):
        return (-b*g + S(2)*c*f)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d) - (-d*g + e*f)*Int((b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d)
    rule355 = ReplacementRule(pattern355, replacement355)
    rubi.add(rule355.pattern, label = 355)

    pattern356 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons45, cons321, cons147, cons239, cons319, cons270, cons31, cons94)
    def replacement356(m, p, g, e, x, b, d, f, c, a):
        return -(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**p/(e*(m + S(1))*(-b*e + S(2)*c*d)) + (S(2)*c*e*f*(m + S(2)*p + S(2)) - g*(b*e*(m + S(1)) + S(2)*c*d*(S(2)*p + S(1))))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(1))*(-b*e + S(2)*c*d))
    rule356 = ReplacementRule(pattern356, replacement356)
    rubi.add(rule356.pattern, label = 356)

    pattern357 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons45, cons321, cons147, cons239, cons319, cons270, cons272, cons324)
    def replacement357(m, p, g, e, x, b, d, f, c, a):
        return g*(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(S(2)*c*e*(m + S(2)*p + S(2))) + (S(2)*c*e*f*(m + S(2)*p + S(2)) - g*(b*e*(m + S(1)) + S(2)*c*(S(2)*d*p + d)))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*e*(m + S(2)*p + S(2)))
    rule357 = ReplacementRule(pattern357, replacement357)
    rubi.add(rule357.pattern, label = 357)

    pattern358 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons315, cons45, cons147)
    def replacement358(m, n, p, g, e, x, b, d, f, c, a):
        return c**(-IntPart(p))*(b/S(2) + c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b/S(2) + c*x)**(S(2)*p)*(d + e*x)**m*(f + g*x)**n, x)
    rule358 = ReplacementRule(pattern358, replacement358)
    rubi.add(rule358.pattern, label = 358)

    pattern359 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons315, cons226, cons256, cons38)
    def replacement359(m, n, p, d, g, e, x, b, f, c, a):
        return Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x)
    rule359 = ReplacementRule(pattern359, replacement359)
    rubi.add(rule359.pattern, label = 359)

    pattern360 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons315, cons257, cons325)
    def replacement360(m, n, p, d, g, e, x, f, c, a):
        return Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x)
    rule360 = ReplacementRule(pattern360, replacement360)
    rubi.add(rule360.pattern, label = 360)

    pattern361 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons256, cons84, cons314)
    def replacement361(m, n, p, d, g, e, x, b, f, c, a):
        return d**m*e**m*Int((f + g*x)**n*(a*e + c*d*x)**(-m)*(a + b*x + c*x**S(2))**(m + p), x)
    rule361 = ReplacementRule(pattern361, replacement361)
    rubi.add(rule361.pattern, label = 361)

    pattern362 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons257, cons84, cons314)
    def replacement362(m, n, p, d, g, e, x, f, c, a):
        return d**m*e**m*Int((a + c*x**S(2))**(m + p)*(f + g*x)**n*(a*e + c*d*x)**(-m), x)
    rule362 = ReplacementRule(pattern362, replacement362)
    rubi.add(rule362.pattern, label = 362)

    pattern363 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons256, cons326)
    def replacement363(m, p, g, e, x, b, d, f, c, a):
        return g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2)))
    rule363 = ReplacementRule(pattern363, replacement363)
    rubi.add(rule363.pattern, label = 363)

    pattern364 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons257, cons327)
    def replacement364(m, p, d, g, e, x, f, c, a):
        return g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2)))
    rule364 = ReplacementRule(pattern364, replacement364)
    rubi.add(rule364.pattern, label = 364)

    pattern365 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons256, cons244, cons137, cons168)
    def replacement365(m, p, g, e, x, b, d, f, c, a):
        return -e*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-b*e + S(2)*c*d)) + (d + e*x)**m*(c*e*f + g*(-b*e + c*d))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))*(-b*e + S(2)*c*d))
    rule365 = ReplacementRule(pattern365, replacement365)
    rubi.add(rule365.pattern, label = 365)

    pattern366 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons257, cons244, cons137, cons168)
    def replacement366(m, p, g, e, x, d, f, c, a):
        return -e*(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*c*d*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g + e*f)/(S(2)*c*d*(p + S(1)))
    rule366 = ReplacementRule(pattern366, replacement366)
    rubi.add(rule366.pattern, label = 366)

    pattern367 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons256, cons139, cons328, cons54)
    def replacement367(m, p, g, e, x, b, d, f, c, a):
        return -e*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-b*e + S(2)*c*d)) + (d + e*x)**m*(c*e*f + g*(-b*e + c*d))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))*(-b*e + S(2)*c*d))
    rule367 = ReplacementRule(pattern367, replacement367)
    rubi.add(rule367.pattern, label = 367)

    pattern368 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons257, cons139, cons328, cons54)
    def replacement368(m, p, d, g, e, x, f, c, a):
        return -e*(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*c*d*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g + e*f)/(S(2)*c*d*(p + S(1)))
    rule368 = ReplacementRule(pattern368, replacement368)
    rubi.add(rule368.pattern, label = 368)

    pattern369 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons256, cons329, cons253)
    def replacement369(m, p, g, e, x, b, d, f, c, a):
        return (d + e*x)**m*(d*g - e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))) + (e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(-b*e + S(2)*c*d)*(m + p + S(1)))
    rule369 = ReplacementRule(pattern369, replacement369)
    rubi.add(rule369.pattern, label = 369)

    pattern370 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons257, cons329, cons253)
    def replacement370(m, p, d, g, e, x, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g - e*f)/(S(2)*c*d*(m + p + S(1))) + (S(2)*c*e*f*(p + S(1)) + m*(c*d*g + c*e*f))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*c*d*e*(m + p + S(1)))
    rule370 = ReplacementRule(pattern370, replacement370)
    rubi.add(rule370.pattern, label = 370)

    pattern371 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons256, cons319)
    def replacement371(m, p, g, e, x, b, d, f, c, a):
        return g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))) + (e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(c*e*(m + S(2)*p + S(2)))
    rule371 = ReplacementRule(pattern371, replacement371)
    rubi.add(rule371.pattern, label = 371)

    pattern372 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons257, cons319)
    def replacement372(m, p, d, g, e, x, f, c, a):
        return (S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**p*(d + e*x)**m, x)/(e*(m + S(2)*p + S(2))) + g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2)))
    rule372 = ReplacementRule(pattern372, replacement372)
    rubi.add(rule372.pattern, label = 372)

    pattern373 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons2, cons7, cons125, cons208, cons330, cons13, cons331)
    def replacement373(p, g, x, f, c, a):
        return x**S(2)*(a + c*x**S(2))**(p + S(1))*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))) - Int(x*(a + c*x**S(2))**(p + S(1))*Simp(S(2)*a*g - c*f*x*(S(2)*p + S(5)), x), x)/(S(2)*a*c*(p + S(1)))
    rule373 = ReplacementRule(pattern373, replacement373)
    rubi.add(rule373.pattern, label = 373)

    pattern374 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons2, cons7, cons125, cons208, cons5, cons330)
    def replacement374(p, g, x, f, c, a):
        return -f**S(2)*Int((a + c*x**S(2))**(p + S(1))/(f - g*x), x)/c + Int((a + c*x**S(2))**(p + S(1))*(f + g*x), x)/c
    rule374 = ReplacementRule(pattern374, replacement374)
    rubi.add(rule374.pattern, label = 374)

    pattern375 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons256, cons147, cons150, cons13, cons332)
    def replacement375(m, n, p, d, g, e, x, b, f, c, a):
        return Int((f + g*x)**n*(a/d + c*x/e)**(-m)*(a + b*x + c*x**S(2))**(m + p), x)
    rule375 = ReplacementRule(pattern375, replacement375)
    rubi.add(rule375.pattern, label = 375)

    pattern376 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons257, cons147, cons150, cons13, cons332)
    def replacement376(m, n, p, d, g, e, x, f, c, a):
        return a**(-m)*d**(S(2)*m)*Int((a + c*x**S(2))**(m + p)*(d - e*x)**(-m)*(f + g*x)**n, x)
    rule376 = ReplacementRule(pattern376, replacement376)
    rubi.add(rule376.pattern, label = 376)

    pattern377 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons256, cons147, cons148, cons333)
    def replacement377(n, p, d, g, e, x, b, f, c, a):
        return -(f + g*x)**n*(a*(-b*e + S(2)*c*d) + c*x*(-S(2)*a*e + b*d))*(a + b*x + c*x**S(2))**p/(d*e*p*(-S(4)*a*c + b**S(2))) - Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p*Simp(-S(2)*a*c*(d*g*n - e*f*(S(2)*p + S(1))) + b*(a*e*g*n - c*d*f*(S(2)*p + S(1))) - c*g*x*(-S(2)*a*e + b*d)*(n + S(2)*p + S(1)), x), x)/(d*e*p*(-S(4)*a*c + b**S(2)))
    rule377 = ReplacementRule(pattern377, replacement377)
    rubi.add(rule377.pattern, label = 377)

    pattern378 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons257, cons147, cons148, cons333)
    def replacement378(n, p, d, g, e, x, f, c, a):
        return (a + c*x**S(2))**p*(d - e*x)*(f + g*x)**n/(S(2)*d*e*p) - Int((a + c*x**S(2))**p*(f + g*x)**(n + S(-1))*Simp(d*g*n - e*f*(S(2)*p + S(1)) - e*g*x*(n + S(2)*p + S(1)), x), x)/(S(2)*d*e*p)
    rule378 = ReplacementRule(pattern378, replacement378)
    rubi.add(rule378.pattern, label = 378)

    pattern379 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons256, cons147, cons196, cons333)
    def replacement379(n, p, d, g, e, x, b, f, c, a):
        return -(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p*(a*c*d*(-b*g + S(2)*c*f) - a*e*(S(2)*a*c*g - b**S(2)*g + b*c*f) + c*x*(-a*e*(-b*g + S(2)*c*f) + c*d*(-S(2)*a*g + b*f)))/(d*e*p*(-S(4)*a*c + b**S(2))*(a*g**S(2) - b*f*g + c*f**S(2))) - Int((f + g*x)**n*(a + b*x + c*x**S(2))**p*Simp(S(2)*a*c*(a*e*g**S(2)*(n + S(2)*p + S(1)) + c*f*(-d*g*n + S(2)*e*f*p + e*f)) + b**S(2)*g*(-a*e*g*(n + p + S(1)) + c*d*f*p) + b*c*(a*g*(d*g*(n + S(1)) + e*f*(n - S(2)*p)) - c*d*f**S(2)*(S(2)*p + S(1))) + c*g*x*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f))*(n + S(2)*p + S(2)), x), x)/(d*e*p*(-S(4)*a*c + b**S(2))*(a*g**S(2) - b*f*g + c*f**S(2)))
    rule379 = ReplacementRule(pattern379, replacement379)
    rubi.add(rule379.pattern, label = 379)

    pattern380 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons257, cons147, cons196, cons333)
    def replacement380(n, p, d, g, e, x, f, c, a):
        return (a + c*x**S(2))**p*(f + g*x)**(n + S(1))*(-a*e*g + c*d*f - c*x*(d*g + e*f))/(S(2)*d*e*p*(a*g**S(2) + c*f**S(2))) + Int((a + c*x**S(2))**p*(f + g*x)**n*Simp(a*e*g**S(2)*(n + S(2)*p + S(1)) - c*f*(d*g*n - e*(S(2)*f*p + f)) + c*g*x*(d*g + e*f)*(n + S(2)*p + S(2)), x), x)/(S(2)*d*e*p*(a*g**S(2) + c*f**S(2)))
    rule380 = ReplacementRule(pattern380, replacement380)
    rubi.add(rule380.pattern, label = 380)

    pattern381 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons315, cons226, cons256, cons147, cons41, cons334, cons335)
    def replacement381(m, n, p, d, g, e, x, b, f, c, a):
        return -e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m - n + S(-1)))
    rule381 = ReplacementRule(pattern381, replacement381)
    rubi.add(rule381.pattern, label = 381)

    pattern382 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons315, cons257, cons147, cons41, cons336, cons335)
    def replacement382(m, n, p, d, g, e, x, f, c, a):
        return -e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(m - n + S(-1)))
    rule382 = ReplacementRule(pattern382, replacement382)
    rubi.add(rule382.pattern, label = 382)

    pattern383 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons315, cons226, cons256, cons147, cons41, cons337)
    def replacement383(m, n, p, d, g, e, x, b, f, c, a):
        return -e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((n + S(1))*(-b*e*g + c*d*g + c*e*f))
    rule383 = ReplacementRule(pattern383, replacement383)
    rubi.add(rule383.pattern, label = 383)

    pattern384 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons315, cons257, cons147, cons41, cons337)
    def replacement384(m, n, p, d, g, e, x, f, c, a):
        return -e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/(c*(n + S(1))*(d*g + e*f))
    rule384 = ReplacementRule(pattern384, replacement384)
    rubi.add(rule384.pattern, label = 384)

    pattern385 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons256, cons147, cons41, cons338, cons163, cons89, cons339)
    def replacement385(m, n, p, d, g, e, x, b, f, c, a):
        return c*m*Int((d + e*x)**(m + S(1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e*g*(n + S(1))) + (d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p/(g*(n + S(1)))
    rule385 = ReplacementRule(pattern385, replacement385)
    rubi.add(rule385.pattern, label = 385)

    pattern386 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons257, cons147, cons41, cons338, cons163, cons89, cons339)
    def replacement386(m, n, p, d, g, e, x, f, c, a):
        return c*m*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*(f + g*x)**(n + S(1)), x)/(e*g*(n + S(1))) + (a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1))/(g*(n + S(1)))
    rule386 = ReplacementRule(pattern386, replacement386)
    rubi.add(rule386.pattern, label = 386)

    pattern387 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons315, cons226, cons256, cons147, cons41, cons338, cons163, cons335, cons340, cons341)
    def replacement387(m, n, p, d, g, e, x, b, f, c, a):
        return -(d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p/(g*(m - n + S(-1))) - m*(-b*e*g + c*d*g + c*e*f)*Int((d + e*x)**(m + S(1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*g*(m - n + S(-1)))
    rule387 = ReplacementRule(pattern387, replacement387)
    rubi.add(rule387.pattern, label = 387)

    pattern388 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons4, cons315, cons257, cons147, cons41, cons338, cons163, cons335, cons340, cons341)
    def replacement388(m, n, p, d, g, e, x, f, c, a):
        return -c*m*(d*g + e*f)*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*(f + g*x)**n, x)/(e**S(2)*g*(m - n + S(-1))) - (a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1))/(g*(m - n + S(-1)))
    rule388 = ReplacementRule(pattern388, replacement388)
    rubi.add(rule388.pattern, label = 388)

    pattern389 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons256, cons147, cons41, cons338, cons137, cons88)
    def replacement389(m, n, p, d, g, e, x, b, f, c, a):
        return -e*g*n*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1)))
    rule389 = ReplacementRule(pattern389, replacement389)
    rubi.add(rule389.pattern, label = 389)

    pattern390 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons257, cons147, cons41, cons338, cons137, cons88)
    def replacement390(m, n, p, d, g, e, x, f, c, a):
        return -e*g*n*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(p + S(1)))
    rule390 = ReplacementRule(pattern390, replacement390)
    rubi.add(rule390.pattern, label = 390)

    pattern391 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons315, cons226, cons256, cons147, cons41, cons338, cons137)
    def replacement391(m, n, p, d, g, e, x, b, f, c, a):
        return e**S(2)*g*(m - n + S(-2))*Int((d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-b*e*g + c*d*g + c*e*f)) + e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-b*e*g + c*d*g + c*e*f))
    rule391 = ReplacementRule(pattern391, replacement391)
    rubi.add(rule391.pattern, label = 391)

    pattern392 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons4, cons315, cons257, cons147, cons41, cons338, cons137)
    def replacement392(m, n, p, d, g, e, x, f, c, a):
        return e**S(2)*g*(m - n + S(-2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(c*(p + S(1))*(d*g + e*f)) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/(c*(p + S(1))*(d*g + e*f))
    rule392 = ReplacementRule(pattern392, replacement392)
    rubi.add(rule392.pattern, label = 392)

    pattern393 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons256, cons147, cons41, cons87, cons88, cons335, cons342)
    def replacement393(m, n, p, d, g, e, x, b, f, c, a):
        return -e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m - n + S(-1))) - n*(-b*e*g + c*d*g + c*e*f)*Int((d + e*x)**m*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*e*(m - n + S(-1)))
    rule393 = ReplacementRule(pattern393, replacement393)
    rubi.add(rule393.pattern, label = 393)

    pattern394 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons257, cons147, cons41, cons87, cons88, cons335, cons342)
    def replacement394(m, n, p, d, g, e, x, f, c, a):
        return -n*(d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(-1)), x)/(e*(m - n + S(-1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(m - n + S(-1)))
    rule394 = ReplacementRule(pattern394, replacement394)
    rubi.add(rule394.pattern, label = 394)

    pattern395 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons256, cons147, cons41, cons87, cons89, cons246)
    def replacement395(m, n, p, d, g, e, x, b, f, c, a):
        return -c*e*(m - n + S(-2))*Int((d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p, x)/((n + S(1))*(-b*e*g + c*d*g + c*e*f)) - e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((n + S(1))*(-b*e*g + c*d*g + c*e*f))
    rule395 = ReplacementRule(pattern395, replacement395)
    rubi.add(rule395.pattern, label = 395)

    pattern396 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons257, cons147, cons41, cons87, cons89, cons246)
    def replacement396(m, n, p, d, g, e, x, f, c, a):
        return -e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/((n + S(1))*(c*d*g + c*e*f)) - e*(m - n + S(-2))*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1)), x)/((n + S(1))*(d*g + e*f))
    rule396 = ReplacementRule(pattern396, replacement396)
    rubi.add(rule396.pattern, label = 396)

    pattern397 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/((x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons256)
    def replacement397(d, g, e, x, b, f, c, a):
        return S(2)*e**S(2)*Subst(Int(S(1)/(-b*e*g + c*(d*g + e*f) + e**S(2)*g*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))/sqrt(d + e*x))
    rule397 = ReplacementRule(pattern397, replacement397)
    rubi.add(rule397.pattern, label = 397)

    pattern398 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons257)
    def replacement398(d, g, e, x, f, c, a):
        return S(2)*e**S(2)*Subst(Int(S(1)/(c*(d*g + e*f) + e**S(2)*g*x**S(2)), x), x, sqrt(a + c*x**S(2))/sqrt(d + e*x))
    rule398 = ReplacementRule(pattern398, replacement398)
    rubi.add(rule398.pattern, label = 398)

    pattern399 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons315, cons226, cons256, cons147, cons343, cons344, cons126)
    def replacement399(m, n, p, d, g, e, x, b, f, c, a):
        return e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*g*(n + p + S(2)))
    rule399 = ReplacementRule(pattern399, replacement399)
    rubi.add(rule399.pattern, label = 399)

    pattern400 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons315, cons257, cons147, cons343, cons345, cons126)
    def replacement400(m, n, p, d, g, e, x, f, c, a):
        return e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))/(c*g*(n + p + S(2)))
    rule400 = ReplacementRule(pattern400, replacement400)
    rubi.add(rule400.pattern, label = 400)

    pattern401 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons256, cons147, cons343, cons87, cons89, cons246)
    def replacement401(m, n, p, d, g, e, x, b, f, c, a):
        return e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(g*(n + S(1))*(-b*e*g + c*d*g + c*e*f)) - e*(b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p, x)/(g*(n + S(1))*(-b*e*g + c*d*g + c*e*f))
    rule401 = ReplacementRule(pattern401, replacement401)
    rubi.add(rule401.pattern, label = 401)

    pattern402 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons257, cons147, cons343, cons87, cons89, cons246)
    def replacement402(m, n, p, d, g, e, x, f, c, a):
        return -e*(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1)), x)/(g*(n + S(1))*(d*g + e*f)) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(-d*g + e*f)/(c*g*(n + S(1))*(d*g + e*f))
    rule402 = ReplacementRule(pattern402, replacement402)
    rubi.add(rule402.pattern, label = 402)

    pattern403 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons315, cons226, cons256, cons147, cons343, cons346, cons246)
    def replacement403(m, n, p, d, g, e, x, b, f, c, a):
        return e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*g*(n + p + S(2))) - (b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))*Int((d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x)/(c*g*(n + p + S(2)))
    rule403 = ReplacementRule(pattern403, replacement403)
    rubi.add(rule403.pattern, label = 403)

    pattern404 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons315, cons257, cons147, cons343, cons346, cons246)
    def replacement404(m, n, p, d, g, e, x, f, c, a):
        return -(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(g*(n + p + S(2))) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))/(c*g*(n + p + S(2)))
    rule404 = ReplacementRule(pattern404, replacement404)
    rubi.add(rule404.pattern, label = 404)

    pattern405 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons315, cons226, cons256, cons147, cons207)
    def replacement405(m, n, p, d, g, e, x, b, f, c, a):
        return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x)
    rule405 = ReplacementRule(pattern405, replacement405)
    rubi.add(rule405.pattern, label = 405)

    pattern406 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons315, cons257, cons347, cons150, cons348, cons349)
    def replacement406(m, n, p, d, g, e, x, f, c, a):
        return Int(ExpandIntegrand(S(1)/sqrt(a + c*x**S(2)), (a + c*x**S(2))**(p + S(1)/2)*(d + e*x)**m*(f + g*x)**n, x), x)
    rule406 = ReplacementRule(pattern406, replacement406)
    rubi.add(rule406.pattern, label = 406)

    pattern407 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons315, cons257, cons147, cons207)
    def replacement407(m, n, p, d, g, e, x, f, c, a):
        return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x)
    rule407 = ReplacementRule(pattern407, replacement407)
    rubi.add(rule407.pattern, label = 407)

    pattern408 = Pattern(Integral(x_**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons226, cons256)
    def replacement408(p, x, e, b, d, c, a):
        return d**S(2)*Int((a + b*x + c*x**S(2))**p/(d + e*x), x)/e**S(2) - Int((d - e*x)*(a + b*x + c*x**S(2))**p, x)/e**S(2)
    rule408 = ReplacementRule(pattern408, replacement408)
    rubi.add(rule408.pattern, label = 408)

    pattern409 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons5, cons257)
    def replacement409(p, x, e, d, c, a):
        return d**S(2)*Int((a + c*x**S(2))**p/(d + e*x), x)/e**S(2) - Int((a + c*x**S(2))**p*(d - e*x), x)/e**S(2)
    rule409 = ReplacementRule(pattern409, replacement409)
    rubi.add(rule409.pattern, label = 409)

    pattern410 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons256, cons147, cons270)
    def replacement410(m, p, d, g, e, x, b, f, c, a):
        return g*(d + e*x)**m*(f + g*x)*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(3))) - Int((d + e*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*e*g*(d*g + e*f*(m + p + S(1))) - c*(d**S(2)*g**S(2) + d*e*f*g*m + e**S(2)*f**S(2)*(m + S(2)*p + S(3))) + e*g*x*(b*e*g*(m + p + S(2)) - c*(d*g*m + e*f*(m + S(2)*p + S(4)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(3)))
    rule410 = ReplacementRule(pattern410, replacement410)
    rubi.add(rule410.pattern, label = 410)

    pattern411 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**S(2), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons257, cons147, cons270)
    def replacement411(m, p, d, g, e, x, f, c, a):
        return g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(f + g*x)/(c*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(d + e*x)**m*Simp(-c*e*g*x*(d*g*m + e*f*(m + S(2)*p + S(4))) - c*(d**S(2)*g**S(2) + d*e*f*g*m + e**S(2)*f**S(2)*(m + S(2)*p + S(3))), x), x)/(c*e**S(2)*(m + S(2)*p + S(3)))
    rule411 = ReplacementRule(pattern411, replacement411)
    rubi.add(rule411.pattern, label = 411)

    pattern412 = Pattern(Integral((x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons3, cons7, cons48, cons125, cons208, cons21, cons4, cons147)
    def replacement412(m, n, p, g, e, x, b, f, c):
        return x**(-m - p)*(e*x)**m*(b + c*x)**(-p)*(b*x + c*x**S(2))**p*Int(x**(m + p)*(b + c*x)**p*(f + g*x)**n, x)
    rule412 = ReplacementRule(pattern412, replacement412)
    rubi.add(rule412.pattern, label = 412)

    pattern413 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons315, cons257, cons147, cons43, cons268)
    def replacement413(m, n, p, d, g, e, x, f, c, a):
        return Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x)
    rule413 = ReplacementRule(pattern413, replacement413)
    rubi.add(rule413.pattern, label = 413)

    pattern414 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons315, cons226, cons256, cons147)
    def replacement414(m, n, p, d, g, e, x, b, f, c, a):
        return (d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x)
    rule414 = ReplacementRule(pattern414, replacement414)
    rubi.add(rule414.pattern, label = 414)

    pattern415 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons315, cons257, cons147)
    def replacement415(m, n, p, d, g, e, x, f, c, a):
        return (a + c*x**S(2))**FracPart(p)*(d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x)
    rule415 = ReplacementRule(pattern415, replacement415)
    rubi.add(rule415.pattern, label = 415)

    pattern416 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons315, cons226, cons279, cons128)
    def replacement416(m, p, g, e, x, b, d, f, c, a):
        return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)*(a + b*x + c*x**S(2))**p, x), x)
    rule416 = ReplacementRule(pattern416, replacement416)
    rubi.add(rule416.pattern, label = 416)

    pattern417 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons280, cons128)
    def replacement417(m, p, g, e, x, d, f, c, a):
        return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x), x), x)
    rule417 = ReplacementRule(pattern417, replacement417)
    rubi.add(rule417.pattern, label = 417)

    pattern418 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279)
    def replacement418(g, e, x, b, d, f, c, a):
        return e*(-d*g + e*f)*Int(S(1)/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int(Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2))
    rule418 = ReplacementRule(pattern418, replacement418)
    rubi.add(rule418.pattern, label = 418)

    pattern419 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280)
    def replacement419(g, e, x, d, f, c, a):
        return e*(-d*g + e*f)*Int(S(1)/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int(Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2))
    rule419 = ReplacementRule(pattern419, replacement419)
    rubi.add(rule419.pattern, label = 419)

    pattern420 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons279, cons242, cons350)
    def replacement420(m, p, g, e, x, b, d, f, c, a):
        return -(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*(p + S(1))*(a*e**S(2) - b*d*e + c*d**S(2)))
    rule420 = ReplacementRule(pattern420, replacement420)
    rubi.add(rule420.pattern, label = 420)

    pattern421 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons280, cons242, cons351)
    def replacement421(m, p, g, e, x, d, f, c, a):
        return -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/(S(2)*(p + S(1))*(a*e**S(2) + c*d**S(2)))
    rule421 = ReplacementRule(pattern421, replacement421)
    rubi.add(rule421.pattern, label = 421)

    pattern422 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons242, cons13, cons137, cons352)
    def replacement422(m, p, g, e, x, b, d, f, c, a):
        return -m*(-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*g + b*f + x*(-b*g + S(2)*c*f))/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule422 = ReplacementRule(pattern422, replacement422)
    rubi.add(rule422.pattern, label = 422)

    pattern423 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons242, cons13, cons137, cons353)
    def replacement423(m, p, g, e, x, d, f, c, a):
        return -m*(a*e*g + c*d*f)*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*c*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(a*g - c*f*x)/(S(2)*a*c*(p + S(1)))
    rule423 = ReplacementRule(pattern423, replacement423)
    rubi.add(rule423.pattern, label = 423)

    pattern424 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons279, cons242)
    def replacement424(m, p, g, e, x, b, d, f, c, a):
        return -(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*(p + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) - (-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*a*e**S(2) - S(2)*b*d*e + S(2)*c*d**S(2))
    rule424 = ReplacementRule(pattern424, replacement424)
    rubi.add(rule424.pattern, label = 424)

    pattern425 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons280, cons242)
    def replacement425(m, p, g, e, x, d, f, c, a):
        return -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/(S(2)*(p + S(1))*(a*e**S(2) + c*d**S(2))) + (a*e*g + c*d*f)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(a*e**S(2) + c*d**S(2))
    rule425 = ReplacementRule(pattern425, replacement425)
    rubi.add(rule425.pattern, label = 425)

    pattern426 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons226, cons279, cons354)
    def replacement426(p, g, e, x, b, d, f, c, a):
        return -(a + b*x + c*x**S(2))**(p + S(1))*(b*e*g*(p + S(2)) - S(2)*c*e*g*x*(p + S(1)) - c*(S(2)*p + S(3))*(d*g + e*f))/(S(2)*c**S(2)*(p + S(1))*(S(2)*p + S(3)))
    rule426 = ReplacementRule(pattern426, replacement426)
    rubi.add(rule426.pattern, label = 426)

    pattern427 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons280, cons355)
    def replacement427(p, g, e, x, d, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(S(2)*e*g*x*(p + S(1)) + (S(2)*p + S(3))*(d*g + e*f))/(S(2)*c*(p + S(1))*(S(2)*p + S(3)))
    rule427 = ReplacementRule(pattern427, replacement427)
    rubi.add(rule427.pattern, label = 427)

    pattern428 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons13, cons137)
    def replacement428(p, g, e, x, b, d, f, c, a):
        return -(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f) - x*(b**S(2)*e*g - b*c*(d*g + e*f) + S(2)*c*(-a*e*g + c*d*f)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - (-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule428 = ReplacementRule(pattern428, replacement428)
    rubi.add(rule428.pattern, label = 428)

    pattern429 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons13, cons137)
    def replacement429(p, g, e, x, d, f, c, a):
        return -(a + c*x**S(2))**(p + S(1))*(-a*(d*g + e*(f + g*x)) + c*d*f*x)/(S(2)*a*c*(p + S(1))) - (a*e*g - c*d*f*(S(2)*p + S(3)))*Int((a + c*x**S(2))**(p + S(1)), x)/(S(2)*a*c*(p + S(1)))
    rule429 = ReplacementRule(pattern429, replacement429)
    rubi.add(rule429.pattern, label = 429)

    pattern430 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons226, cons279, cons287)
    def replacement430(p, g, e, x, b, d, f, c, a):
        return (-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c**S(2)*(S(2)*p + S(3))) - (a + b*x + c*x**S(2))**(p + S(1))*(b*e*g*(p + S(2)) - S(2)*c*e*g*x*(p + S(1)) - c*(S(2)*p + S(3))*(d*g + e*f))/(S(2)*c**S(2)*(p + S(1))*(S(2)*p + S(3)))
    rule430 = ReplacementRule(pattern430, replacement430)
    rubi.add(rule430.pattern, label = 430)

    pattern431 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons280, cons287)
    def replacement431(p, g, e, x, d, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(S(2)*e*g*x*(p + S(1)) + (S(2)*p + S(3))*(d*g + e*f))/(S(2)*c*(p + S(1))*(S(2)*p + S(3))) - (a*e*g - c*d*f*(S(2)*p + S(3)))*Int((a + c*x**S(2))**p, x)/(c*(S(2)*p + S(3)))
    rule431 = ReplacementRule(pattern431, replacement431)
    rubi.add(rule431.pattern, label = 431)

    pattern432 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons2, cons7, cons48, cons125, cons208, cons5, cons356, cons357)
    def replacement432(m, p, g, x, e, f, c, a):
        return f*Int((e*x)**m*(a + c*x**S(2))**p, x) + g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p, x)/e
    rule432 = ReplacementRule(pattern432, replacement432)
    rubi.add(rule432.pattern, label = 432)

    pattern433 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons358, cons288, cons289)
    def replacement433(m, p, g, e, x, b, d, f, c, a):
        return (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((f + g*x)*(a*d + c*e*x**S(3))**p, x)
    rule433 = ReplacementRule(pattern433, replacement433)
    rubi.add(rule433.pattern, label = 433)

    pattern434 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons244, cons163, cons247, cons359, cons360)
    def replacement434(m, p, g, e, x, b, d, f, c, a):
        return -p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(S(2)*a*c*e*(m + S(2))*(-d*g + e*f) + b**S(2)*e*(d*g*(p + S(1)) - e*f*(m + p + S(2))) + b*(a*e**S(2)*g*(m + S(1)) - c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2)))) - c*x*(S(2)*c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2))) - e*(S(2)*a*e*g*(m + S(1)) - b*(d*g*(m - S(2)*p) + e*f*(m + S(2)*p + S(2))))), x), x)/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) - (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-d*p*(-b*e + S(2)*c*d)*(-d*g + e*f) - e*x*(g*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2)) + p*(-b*e + S(2)*c*d)*(-d*g + e*f)) + (d*g - e*f*(m + S(2)))*(a*e**S(2) - b*d*e + c*d**S(2)))/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) - b*d*e + c*d**S(2)))
    rule434 = ReplacementRule(pattern434, replacement434)
    rubi.add(rule434.pattern, label = 434)

    pattern435 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons244, cons163, cons247, cons359, cons360)
    def replacement435(m, p, g, e, x, d, f, c, a):
        return -p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2))*Simp(S(2)*a*c*e*(m + S(2))*(-d*g + e*f) - c*x*(-S(2)*a*e**S(2)*g*(m + S(1)) + S(2)*c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2)))), x), x)/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) + c*d**S(2))) - (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-S(2)*c*d**S(2)*p*(-d*g + e*f) - e*x*(S(2)*c*d*p*(-d*g + e*f) + g*(m + S(1))*(a*e**S(2) + c*d**S(2))) + (a*e**S(2) + c*d**S(2))*(d*g - e*f*(m + S(2))))/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) + c*d**S(2)))
    rule435 = ReplacementRule(pattern435, replacement435)
    rubi.add(rule435.pattern, label = 435)

    pattern436 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons315, cons226, cons279, cons13, cons163, cons361, cons66, cons299, cons362)
    def replacement436(m, p, g, e, x, b, d, f, c, a):
        return p*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-b*e*f*(m + S(2)*p + S(2)) + g*(S(2)*a*e*m + S(2)*a*e + S(2)*b*d*p + b*d) + x*(-S(2)*c*e*f*(m + S(2)*p + S(2)) + g*(b*e*m + b*e + S(4)*c*d*p + S(2)*c*d)), x), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-d*g*(S(2)*p + S(1)) + e*f*(m + S(2)*p + S(2)) + e*g*x*(m + S(1)))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2)))
    rule436 = ReplacementRule(pattern436, replacement436)
    rubi.add(rule436.pattern, label = 436)

    pattern437 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons315, cons280, cons13, cons163, cons361, cons66, cons299, cons362)
    def replacement437(m, p, g, e, x, d, f, c, a):
        return p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*Simp(g*(S(2)*a*e*m + S(2)*a*e) + x*(-S(2)*c*e*f*(m + S(2)*p + S(2)) + g*(S(4)*c*d*p + S(2)*c*d)), x), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-d*g*(S(2)*p + S(1)) + e*f*(m + S(2)*p + S(2)) + e*g*x*(m + S(1)))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2)))
    rule437 = ReplacementRule(pattern437, replacement437)
    rubi.add(rule437.pattern, label = 437)

    pattern438 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons315, cons226, cons279, cons13, cons163, cons363, cons303, cons362)
    def replacement438(m, p, g, e, x, b, d, f, c, a):
        return -p*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(c*e*f*(-S(2)*a*e + b*d)*(m + S(2)*p + S(2)) + g*(a*e*(b*e*m + b*e - S(2)*c*d*m) + b*d*(b*e*p - S(2)*c*d*p - c*d)) + x*(c*e*f*(-b*e + S(2)*c*d)*(m + S(2)*p + S(2)) + g*(b**S(2)*e**S(2)*(m + p + S(1)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(1)) - c*e*(S(2)*a*e*(m + S(2)*p + S(1)) + b*d*(m - S(2)*p)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(c*e*f*(m + S(2)*p + S(2)) + c*e*g*x*(m + S(2)*p + S(1)) - g*(-b*e*p + S(2)*c*d*p + c*d))/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2)))
    rule438 = ReplacementRule(pattern438, replacement438)
    rubi.add(rule438.pattern, label = 438)

    pattern439 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons315, cons280, cons13, cons163, cons363, cons303, cons362)
    def replacement439(m, p, g, e, x, d, f, c, a):
        return S(2)*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**m*Simp(a*c*d*e*g*m + a*c*e**S(2)*f*(m + S(2)*p + S(2)) - x*(c**S(2)*d*e*f*(m + S(2)*p + S(2)) - g*(a*c*e**S(2)*(m + S(2)*p + S(1)) + c**S(2)*d**S(2)*(S(2)*p + S(1)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-c*d*g*(S(2)*p + S(1)) + c*e*f*(m + S(2)*p + S(2)) + c*e*g*x*(m + S(2)*p + S(1)))/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2)))
    rule439 = ReplacementRule(pattern439, replacement439)
    rubi.add(rule439.pattern, label = 439)

    pattern440 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons244, cons137, cons166, cons364)
    def replacement440(m, p, g, e, x, b, d, f, c, a):
        return -(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f) - x*(b**S(2)*e*g + S(2)*c**S(2)*d*f - c*(S(2)*a*e*g + b*d*g + b*e*f)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(b*e*g*(a*e*(m + S(-1)) + b*d*(p + S(2))) + S(2)*c**S(2)*d**S(2)*f*(S(2)*p + S(3)) - c*(S(2)*a*e*(d*g*m + e*f*(m + S(-1))) + b*d*(d*g*(S(2)*p + S(3)) - e*f*(m - S(2)*p + S(-4)))) + e*x*(b**S(2)*e*g*(m + p + S(1)) + S(2)*c**S(2)*d*f*(m + S(2)*p + S(2)) - c*(S(2)*a*e*g*m + b*(d*g + e*f)*(m + S(2)*p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule440 = ReplacementRule(pattern440, replacement440)
    rubi.add(rule440.pattern, label = 440)

    pattern441 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons244, cons137, cons166, cons365)
    def replacement441(m, p, g, e, x, d, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(S(2)*a*(d*g + e*f) - x*(-S(2)*a*e*g + S(2)*c*d*f))/(S(4)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*Simp(S(2)*a*e*(d*g*m + e*f*(m + S(-1))) - S(2)*c*d**S(2)*f*(S(2)*p + S(3)) + e*x*(S(2)*a*e*g*m - S(2)*c*d*f*(m + S(2)*p + S(2))), x), x)/(S(4)*a*c*(p + S(1)))
    rule441 = ReplacementRule(pattern441, replacement441)
    rubi.add(rule441.pattern, label = 441)

    pattern442 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons244, cons137, cons168, cons366, cons362)
    def replacement442(m, p, g, e, x, b, d, f, c, a):
        return (d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*g + b*f + x*(-b*g + S(2)*c*f))/((p + S(1))*(-S(4)*a*c + b**S(2))) + Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-e*x*(-b*g + S(2)*c*f)*(m + S(2)*p + S(3)) - f*(b*e*m + S(2)*c*d*(S(2)*p + S(3))) + g*(S(2)*a*e*m + b*d*(S(2)*p + S(3))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule442 = ReplacementRule(pattern442, replacement442)
    rubi.add(rule442.pattern, label = 442)

    pattern443 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons244, cons137, cons168, cons366, cons362)
    def replacement443(m, p, g, e, x, d, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*Simp(a*e*g*m - c*d*f*(S(2)*p + S(3)) - c*e*f*x*(m + S(2)*p + S(3)), x), x)/(S(2)*a*c*(p + S(1)))
    rule443 = ReplacementRule(pattern443, replacement443)
    rubi.add(rule443.pattern, label = 443)

    pattern444 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons315, cons226, cons279, cons13, cons137, cons362)
    def replacement444(m, p, g, e, x, b, d, f, c, a):
        return (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-a*g*(-b*e + S(2)*c*d) + c*x*(f*(-b*e + S(2)*c*d) - g*(-S(2)*a*e + b*d)) + f*(S(2)*a*c*e - b**S(2)*e + b*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(c*e*x*(-f*(-b*e + S(2)*c*d) + g*(-S(2)*a*e + b*d))*(m + S(2)*p + S(4)) + f*(-S(2)*a*c*e**S(2)*(m + S(2)*p + S(3)) + b**S(2)*e**S(2)*(m + p + S(2)) + b*c*d*e*(-m + S(2)*p + S(2)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(3))) - g*(a*e*(b*e*m + b*e - S(2)*c*d*m) - b*d*(-b*e*p - b*e + S(2)*c*d*p + S(3)*c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2)))
    rule444 = ReplacementRule(pattern444, replacement444)
    rubi.add(rule444.pattern, label = 444)

    pattern445 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons13, cons137, cons362)
    def replacement445(m, p, g, e, x, d, f, c, a):
        return -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-a*c*d*g + a*c*e*f + c*x*(a*e*g + c*d*f))/(S(2)*a*c*(p + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*Simp(-a*c*d*e*g*m + c*e*x*(a*e*g + c*d*f)*(m + S(2)*p + S(4)) + f*(a*c*e**S(2)*(m + S(2)*p + S(3)) + c**S(2)*d**S(2)*(S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1))*(a*e**S(2) + c*d**S(2)))
    rule445 = ReplacementRule(pattern445, replacement445)
    rubi.add(rule445.pattern, label = 445)

    pattern446 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons17)
    def replacement446(m, g, e, x, b, d, f, c, a):
        return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)/(a + b*x + c*x**S(2)), x), x)
    rule446 = ReplacementRule(pattern446, replacement446)
    rubi.add(rule446.pattern, label = 446)

    pattern447 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons17)
    def replacement447(m, g, e, x, d, f, c, a):
        return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)/(a + c*x**S(2)), x), x)
    rule447 = ReplacementRule(pattern447, replacement447)
    rubi.add(rule447.pattern, label = 447)

    pattern448 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons367, cons168)
    def replacement448(m, g, e, x, b, d, f, c, a):
        return g*(d + e*x)**m/(c*m) + Int((d + e*x)**(m + S(-1))*Simp(-a*e*g + c*d*f + x*(-b*e*g + c*d*g + c*e*f), x)/(a + b*x + c*x**S(2)), x)/c
    rule448 = ReplacementRule(pattern448, replacement448)
    rubi.add(rule448.pattern, label = 448)

    pattern449 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons367, cons168)
    def replacement449(m, g, e, x, d, f, c, a):
        return g*(d + e*x)**m/(c*m) + Int((d + e*x)**(m + S(-1))*Simp(-a*e*g + c*d*f + x*(c*d*g + c*e*f), x)/(a + c*x**S(2)), x)/c
    rule449 = ReplacementRule(pattern449, replacement449)
    rubi.add(rule449.pattern, label = 449)

    pattern450 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279)
    def replacement450(g, e, x, b, d, f, c, a):
        return S(2)*Subst(Int((-d*g + e*f + g*x**S(2))/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x))
    rule450 = ReplacementRule(pattern450, replacement450)
    rubi.add(rule450.pattern, label = 450)

    pattern451 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280)
    def replacement451(g, e, x, d, f, c, a):
        return S(2)*Subst(Int((-d*g + e*f + g*x**S(2))/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x))
    rule451 = ReplacementRule(pattern451, replacement451)
    rubi.add(rule451.pattern, label = 451)

    pattern452 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons315, cons226, cons279, cons367, cons94)
    def replacement452(m, g, e, x, b, d, f, c, a):
        return (d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2))
    rule452 = ReplacementRule(pattern452, replacement452)
    rubi.add(rule452.pattern, label = 452)

    pattern453 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons315, cons280, cons367, cons94)
    def replacement453(m, g, e, x, d, f, c, a):
        return (d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2))
    rule453 = ReplacementRule(pattern453, replacement453)
    rubi.add(rule453.pattern, label = 453)

    pattern454 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons356)
    def replacement454(m, g, e, x, b, d, f, c, a):
        return Int(ExpandIntegrand((d + e*x)**m, (f + g*x)/(a + b*x + c*x**S(2)), x), x)
    rule454 = ReplacementRule(pattern454, replacement454)
    rubi.add(rule454.pattern, label = 454)

    pattern455 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons356)
    def replacement455(m, g, e, x, d, f, c, a):
        return Int(ExpandIntegrand((d + e*x)**m, (f + g*x)/(a + c*x**S(2)), x), x)
    rule455 = ReplacementRule(pattern455, replacement455)
    rubi.add(rule455.pattern, label = 455)

    pattern456 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons226, cons279, cons31, cons168, cons319, cons368, cons362)
    def replacement456(m, p, g, e, x, b, d, f, c, a):
        return g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))) + Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p*Simp(d*(p + S(1))*(-b*g + S(2)*c*f) + m*(-a*e*g + c*d*f) + x*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(-b*e*g + c*d*g + c*e*f)), x), x)/(c*(m + S(2)*p + S(2)))
    rule456 = ReplacementRule(pattern456, replacement456)
    rubi.add(rule456.pattern, label = 456)

    pattern457 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons280, cons31, cons168, cons319, cons368, cons362)
    def replacement457(m, p, g, e, x, d, f, c, a):
        return g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*Simp(-a*e*g*m + c*d*f*(m + S(2)*p + S(2)) + c*x*(d*g*m + e*f*(m + S(2)*p + S(2))), x), x)/(c*(m + S(2)*p + S(2)))
    rule457 = ReplacementRule(pattern457, replacement457)
    rubi.add(rule457.pattern, label = 457)

    pattern458 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons226, cons279, cons31, cons94, cons362)
    def replacement458(m, p, g, e, x, b, d, f, c, a):
        return (d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(b*(p + S(1))*(d*g - e*f) - c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g - b*e*f + c*d*f), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2)))
    rule458 = ReplacementRule(pattern458, replacement458)
    rubi.add(rule458.pattern, label = 458)

    pattern459 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons5, cons315, cons280, cons31, cons94, cons362)
    def replacement459(m, p, g, e, x, d, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(-c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g + c*d*f), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2)))
    rule459 = ReplacementRule(pattern459, replacement459)
    rubi.add(rule459.pattern, label = 459)

    pattern460 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons279, cons369, cons66)
    def replacement460(m, p, g, e, x, b, d, f, c, a):
        return (d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(b*(p + S(1))*(d*g - e*f) - c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g - b*e*f + c*d*f), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2)))
    rule460 = ReplacementRule(pattern460, replacement460)
    rubi.add(rule460.pattern, label = 460)

    pattern461 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons280, cons369, cons66)
    def replacement461(m, p, g, e, x, d, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(-c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g + c*d*f), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2)))
    rule461 = ReplacementRule(pattern461, replacement461)
    rubi.add(rule461.pattern, label = 461)

    pattern462 = Pattern(Integral((f_ + x_*WC('g', S(1)))/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons370, cons371, cons372)
    def replacement462(g, e, x, b, d, f, c, a):
        return S(4)*f*(a - d)*Subst(Int(S(1)/(S(4)*a - S(4)*d - x**S(2)), x), x, (S(2)*a - S(2)*d + x*(b - e))/sqrt(a + b*x + c*x**S(2)))/(-a*e + b*d)
    rule462 = ReplacementRule(pattern462, replacement462)
    rubi.add(rule462.pattern, label = 462)

    pattern463 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons125, cons208, cons226)
    def replacement463(g, x, b, f, c, a):
        return S(2)*Subst(Int((f + g*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x, sqrt(x))
    rule463 = ReplacementRule(pattern463, replacement463)
    rubi.add(rule463.pattern, label = 463)

    pattern464 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)))), x_), cons2, cons7, cons125, cons208, cons373)
    def replacement464(g, x, f, c, a):
        return S(2)*Subst(Int((f + g*x**S(2))/sqrt(a + c*x**S(4)), x), x, sqrt(x))
    rule464 = ReplacementRule(pattern464, replacement464)
    rubi.add(rule464.pattern, label = 464)

    pattern465 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(e_*x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons48, cons125, cons208, cons226)
    def replacement465(g, x, e, b, f, c, a):
        return sqrt(x)*Int((f + g*x)/(sqrt(x)*sqrt(a + b*x + c*x**S(2))), x)/sqrt(e*x)
    rule465 = ReplacementRule(pattern465, replacement465)
    rubi.add(rule465.pattern, label = 465)

    pattern466 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(e_*x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)))), x_), cons2, cons7, cons48, cons125, cons208, cons374)
    def replacement466(g, x, e, f, c, a):
        return sqrt(x)*Int((f + g*x)/(sqrt(x)*sqrt(a + c*x**S(2))), x)/sqrt(e*x)
    rule466 = ReplacementRule(pattern466, replacement466)
    rubi.add(rule466.pattern, label = 466)

    pattern467 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons279)
    def replacement467(m, p, g, e, x, b, d, f, c, a):
        return g*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/e
    rule467 = ReplacementRule(pattern467, replacement467)
    rubi.add(rule467.pattern, label = 467)

    pattern468 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons280)
    def replacement468(m, p, g, e, x, d, f, c, a):
        return g*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/e + (-d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m, x)/e
    rule468 = ReplacementRule(pattern468, replacement468)
    rubi.add(rule468.pattern, label = 468)

    pattern469 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons226, cons279, cons375)
    def replacement469(m, n, p, g, e, x, b, d, f, c, a):
        return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x)
    rule469 = ReplacementRule(pattern469, replacement469)
    rubi.add(rule469.pattern, label = 469)

    pattern470 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons7, cons27, cons48, cons125, cons208, cons280, cons375)
    def replacement470(m, n, p, g, e, x, d, f, c, a):
        return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x)
    rule470 = ReplacementRule(pattern470, replacement470)
    rubi.add(rule470.pattern, label = 470)

    pattern471 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons149, cons163)
    def replacement471(p, g, e, x, b, d, f, c, a):
        return (a*e**S(2) - b*d*e + c*d**S(2))*Int((a + b*x + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + b*x + c*x**S(2))**(p + S(-1))*Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(f + g*x), x)/(e*(-d*g + e*f))
    rule471 = ReplacementRule(pattern471, replacement471)
    rubi.add(rule471.pattern, label = 471)

    pattern472 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons149, cons163)
    def replacement472(p, g, e, x, d, f, c, a):
        return (a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + c*x**S(2))**(p + S(-1))*Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(f + g*x), x)/(e*(-d*g + e*f))
    rule472 = ReplacementRule(pattern472, replacement472)
    rubi.add(rule472.pattern, label = 472)

    def With473(m, n, p, g, e, x, b, d, f, c, a):
        q = Denominator(m)
        return q*Subst(Int(x**(q*(m + S(1)) + S(-1))*(g*x**q/e + (-d*g + e*f)/e)**n*(c*x**(S(2)*q)/e**S(2) - x**q*(-b*e + S(2)*c*d)/e**S(2) + (a*e**S(2) - b*d*e + c*d**S(2))/e**S(2))**p, x), x, (d + e*x)**(S(1)/q))/e
    pattern473 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons376, cons367 )
    rule473 = ReplacementRule(pattern473, With473)
    rubi.add(rule473.pattern, label = 473)

    def With474(m, n, p, g, e, x, d, f, c, a):
        q = Denominator(m)
        return q*Subst(Int(x**(q*(m + S(1)) + S(-1))*(g*x**q/e + (-d*g + e*f)/e)**n*(-S(2)*c*d*x**q/e**S(2) + c*x**(S(2)*q)/e**S(2) + (a*e**S(2) + c*d**S(2))/e**S(2))**p, x), x, (d + e*x)**(S(1)/q))/e
    pattern474 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons376, cons367 )
    rule474 = ReplacementRule(pattern474, With474)
    rubi.add(rule474.pattern, label = 474)

    pattern475 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons124, cons336, cons377)
    def replacement475(m, n, p, g, e, x, b, d, f, c, a):
        return Int((d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x)
    rule475 = ReplacementRule(pattern475, replacement475)
    rubi.add(rule475.pattern, label = 475)

    pattern476 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons124, cons336, cons377)
    def replacement476(m, n, p, g, e, x, d, f, c, a):
        return Int((a + c*x**S(2))**p*(d*f + e*g*x**S(2))**m, x)
    rule476 = ReplacementRule(pattern476, replacement476)
    rubi.add(rule476.pattern, label = 476)

    pattern477 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons124, cons336)
    def replacement477(m, n, p, g, e, x, b, d, f, c, a):
        return (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x)
    rule477 = ReplacementRule(pattern477, replacement477)
    rubi.add(rule477.pattern, label = 477)

    pattern478 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons124, cons336)
    def replacement478(m, n, p, g, e, x, d, f, c, a):
        return (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((a + c*x**S(2))**p*(d*f + e*g*x**S(2))**m, x)
    rule478 = ReplacementRule(pattern478, replacement478)
    rubi.add(rule478.pattern, label = 478)

    pattern479 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons226, cons279, cons93)
    def replacement479(m, n, g, e, x, b, d, f, c, a):
        return c*Int(x**S(2)*(d + e*x)**m*(f + g*x)**n, x) + Int((a + b*x)*(d + e*x)**m*(f + g*x)**n, x)
    rule479 = ReplacementRule(pattern479, replacement479)
    rubi.add(rule479.pattern, label = 479)

    pattern480 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons7, cons27, cons48, cons125, cons208, cons280, cons93)
    def replacement480(m, n, g, e, x, d, f, c, a):
        return a*Int((d + e*x)**m*(f + g*x)**n, x) + c*Int(x**S(2)*(d + e*x)**m*(f + g*x)**n, x)
    rule480 = ReplacementRule(pattern480, replacement480)
    rubi.add(rule480.pattern, label = 480)

    pattern481 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons226, cons279, cons18, cons23, cons93, cons168, cons165)
    def replacement481(m, n, g, e, x, b, d, f, c, a):
        return g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(-b*e*g + c*d*g + S(2)*c*e*f + c*e*g*x, x), x)/c**S(2) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(a*b*e*g**S(2) - a*c*d*g**S(2) - S(2)*a*c*e*f*g + c**S(2)*d*f**S(2) + x*(-a*c*e*g**S(2) + b**S(2)*e*g**S(2) - b*c*d*g**S(2) - S(2)*b*c*e*f*g + S(2)*c**S(2)*d*f*g + c**S(2)*e*f**S(2)), x)/(a + b*x + c*x**S(2)), x)/c**S(2)
    rule481 = ReplacementRule(pattern481, replacement481)
    rubi.add(rule481.pattern, label = 481)

    pattern482 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons280, cons18, cons23, cons93, cons168, cons165)
    def replacement482(m, n, g, e, x, d, f, c, a):
        return g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(d*g + S(2)*e*f + e*g*x, x), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(-a*d*g**S(2) - S(2)*a*e*f*g + c*d*f**S(2) + x*(-a*e*g**S(2) + S(2)*c*d*f*g + c*e*f**S(2)), x)/(a + c*x**S(2)), x)/c
    rule482 = ReplacementRule(pattern482, replacement482)
    rubi.add(rule482.pattern, label = 482)

    pattern483 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons226, cons279, cons18, cons23, cons93, cons168, cons88)
    def replacement483(m, n, g, e, x, b, d, f, c, a):
        return e*g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*Simp(-a*e*g + c*d*f + x*(-b*e*g + c*d*g + c*e*f), x)/(a + b*x + c*x**S(2)), x)/c
    rule483 = ReplacementRule(pattern483, replacement483)
    rubi.add(rule483.pattern, label = 483)

    pattern484 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons280, cons18, cons23, cons93, cons168, cons88)
    def replacement484(m, n, g, e, x, d, f, c, a):
        return e*g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*Simp(-a*e*g + c*d*f + x*(c*d*g + c*e*f), x)/(a + c*x**S(2)), x)/c
    rule484 = ReplacementRule(pattern484, replacement484)
    rubi.add(rule484.pattern, label = 484)

    pattern485 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons226, cons279, cons18, cons23, cons93, cons168, cons89)
    def replacement485(m, n, g, e, x, b, d, f, c, a):
        return -g*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(a*g**S(2) - b*f*g + c*f**S(2)) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*Simp(a*e*g - b*d*g + c*d*f + c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*g**S(2) - b*f*g + c*f**S(2))
    rule485 = ReplacementRule(pattern485, replacement485)
    rubi.add(rule485.pattern, label = 485)

    pattern486 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons280, cons18, cons23, cons93, cons168, cons89)
    def replacement486(m, n, g, e, x, d, f, c, a):
        return -g*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(a*g**S(2) + c*f**S(2)) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*Simp(a*e*g + c*d*f + c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*g**S(2) + c*f**S(2))
    rule486 = ReplacementRule(pattern486, replacement486)
    rubi.add(rule486.pattern, label = 486)

    pattern487 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons226, cons279, cons73)
    def replacement487(m, g, e, x, b, d, f, c, a):
        return Int(ExpandIntegrand(S(1)/(sqrt(d + e*x)*sqrt(f + g*x)), (d + e*x)**(m + S(1)/2)/(a + b*x + c*x**S(2)), x), x)
    rule487 = ReplacementRule(pattern487, replacement487)
    rubi.add(rule487.pattern, label = 487)

    pattern488 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons280, cons73)
    def replacement488(m, g, e, x, d, f, c, a):
        return Int(ExpandIntegrand(S(1)/(sqrt(d + e*x)*sqrt(f + g*x)), (d + e*x)**(m + S(1)/2)/(a + c*x**S(2)), x), x)
    rule488 = ReplacementRule(pattern488, replacement488)
    rubi.add(rule488.pattern, label = 488)

    pattern489 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons226, cons279, cons18, cons23)
    def replacement489(m, n, g, e, x, b, d, f, c, a):
        return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n, S(1)/(a + b*x + c*x**S(2)), x), x)
    rule489 = ReplacementRule(pattern489, replacement489)
    rubi.add(rule489.pattern, label = 489)

    pattern490 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons280, cons18, cons23)
    def replacement490(m, n, g, e, x, d, f, c, a):
        return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n, S(1)/(a + c*x**S(2)), x), x)
    rule490 = ReplacementRule(pattern490, replacement490)
    rubi.add(rule490.pattern, label = 490)

    pattern491 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons378)
    def replacement491(m, n, p, g, e, x, b, d, f, c, a):
        return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x)
    rule491 = ReplacementRule(pattern491, replacement491)
    rubi.add(rule491.pattern, label = 491)

    pattern492 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons378)
    def replacement492(m, n, p, g, e, x, d, f, c, a):
        return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x)
    rule492 = ReplacementRule(pattern492, replacement492)
    rubi.add(rule492.pattern, label = 492)

    pattern493 = Pattern(Integral((x_*WC('g', S(1)))**WC('n', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons21, cons4, cons5, cons358, cons288, cons289)
    def replacement493(m, n, p, g, e, x, b, d, c, a):
        return (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((g*x)**n*(a*d + c*e*x**S(3))**p, x)
    rule493 = ReplacementRule(pattern493, replacement493)
    rubi.add(rule493.pattern, label = 493)

    pattern494 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons23, cons147, cons338, cons163, cons89)
    def replacement494(n, p, g, e, x, b, d, f, c, a):
        return (a*e**S(2) - b*d*e + c*d**S(2))*Int((f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(-1))*(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f)), x)/(e*(-d*g + e*f))
    rule494 = ReplacementRule(pattern494, replacement494)
    rubi.add(rule494.pattern, label = 494)

    pattern495 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons23, cons147, cons338, cons163, cons89)
    def replacement495(n, p, g, e, x, d, f, c, a):
        return (a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(-1))*(f + g*x)**(n + S(1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + c*x**S(2))**(p + S(-1))*(f + g*x)**n*(a*e*g + c*d*f - c*x*(-d*g + e*f)), x)/(e*(-d*g + e*f))
    rule495 = ReplacementRule(pattern495, replacement495)
    rubi.add(rule495.pattern, label = 495)

    pattern496 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons23, cons147, cons338, cons137, cons88)
    def replacement496(n, p, g, e, x, b, d, f, c, a):
        return e*(-d*g + e*f)*Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p*(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f)), x)/(a*e**S(2) - b*d*e + c*d**S(2))
    rule496 = ReplacementRule(pattern496, replacement496)
    rubi.add(rule496.pattern, label = 496)

    pattern497 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons23, cons147, cons338, cons137, cons88)
    def replacement497(n, p, g, e, x, d, f, c, a):
        return e*(-d*g + e*f)*Int((a + c*x**S(2))**(p + S(1))*(f + g*x)**(n + S(-1))/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int((a + c*x**S(2))**p*(f + g*x)**(n + S(-1))*(a*e*g + c*d*f - c*x*(-d*g + e*f)), x)/(a*e**S(2) + c*d**S(2))
    rule497 = ReplacementRule(pattern497, replacement497)
    rubi.add(rule497.pattern, label = 497)

    def With498(g, e, x, b, d, f, c, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -sqrt(S(2))*sqrt(-g*(b + S(2)*c*x - q)/(-b*g + S(2)*c*f + g*q))*sqrt(-g*(b + S(2)*c*x + q)/(-b*g + S(2)*c*f - g*q))*EllipticPi(e*(-b*g + S(2)*c*f + g*q)/(S(2)*c*(-d*g + e*f)), asin(sqrt(S(2))*sqrt(c/(-b*g + S(2)*c*f + g*q))*sqrt(f + g*x)), (-b*g + S(2)*c*f + g*q)/(-b*g + S(2)*c*f - g*q))/(sqrt(c/(-b*g + S(2)*c*f + g*q))*(-d*g + e*f)*sqrt(a + b*x + c*x**S(2)))
    pattern498 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279 )
    rule498 = ReplacementRule(pattern498, With498)
    rubi.add(rule498.pattern, label = 498)

    def With499(g, e, x, d, f, c, a):
        q = Rt(-a*c, S(2))
        return -S(2)*sqrt(g*(-c*x + q)/(c*f + g*q))*sqrt(-g*(c*x + q)/(c*f - g*q))*EllipticPi(e*(c*f + g*q)/(c*(-d*g + e*f)), asin(sqrt(c/(c*f + g*q))*sqrt(f + g*x)), (c*f + g*q)/(c*f - g*q))/(sqrt(c/(c*f + g*q))*sqrt(a + c*x**S(2))*(-d*g + e*f))
    pattern499 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280 )
    rule499 = ReplacementRule(pattern499, With499)
    rubi.add(rule499.pattern, label = 499)

    pattern500 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279, cons80)
    def replacement500(n, g, e, x, b, d, f, c, a):
        return Int(ExpandIntegrand(S(1)/(sqrt(f + g*x)*sqrt(a + b*x + c*x**S(2))), (f + g*x)**(n + S(1)/2)/(d + e*x), x), x)
    rule500 = ReplacementRule(pattern500, replacement500)
    rubi.add(rule500.pattern, label = 500)

    pattern501 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280, cons80)
    def replacement501(n, g, e, x, d, f, c, a):
        return Int(ExpandIntegrand(S(1)/(sqrt(a + c*x**S(2))*sqrt(f + g*x)), (f + g*x)**(n + S(1)/2)/(d + e*x), x), x)
    rule501 = ReplacementRule(pattern501, replacement501)
    rubi.add(rule501.pattern, label = 501)

    pattern502 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons315, cons226, cons279)
    def replacement502(g, e, x, b, d, f, c, a):
        return -S(2)*sqrt((-d*g + e*f)**S(2)*(a + b*x + c*x**S(2))/((d + e*x)**S(2)*(a*g**S(2) - b*f*g + c*f**S(2))))*(d + e*x)*Subst(Int(S(1)/sqrt(x**S(4)*(a*e**S(2) - b*d*e + c*d**S(2))/(a*g**S(2) - b*f*g + c*f**S(2)) - x**S(2)*(S(2)*a*e*g - b*d*g - b*e*f + S(2)*c*d*f)/(a*g**S(2) - b*f*g + c*f**S(2)) + S(1)), x), x, sqrt(f + g*x)/sqrt(d + e*x))/((-d*g + e*f)*sqrt(a + b*x + c*x**S(2)))
    rule502 = ReplacementRule(pattern502, replacement502)
    rubi.add(rule502.pattern, label = 502)

    pattern503 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons315, cons280)
    def replacement503(g, e, x, d, f, c, a):
        return -S(2)*sqrt((a + c*x**S(2))*(-d*g + e*f)**S(2)/((d + e*x)**S(2)*(a*g**S(2) + c*f**S(2))))*(d + e*x)*Subst(Int(S(1)/sqrt(x**S(4)*(a*e**S(2) + c*d**S(2))/(a*g**S(2) + c*f**S(2)) - x**S(2)*(S(2)*a*e*g + S(2)*c*d*f)/(a*g**S(2) + c*f**S(2)) + S(1)), x), x, sqrt(f + g*x)/sqrt(d + e*x))/(sqrt(a + c*x**S(2))*(-d*g + e*f))
    rule503 = ReplacementRule(pattern503, replacement503)
    rubi.add(rule503.pattern, label = 503)

    pattern504 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(f_ + x_*WC('g', S(1)))**S(2), x_), cons2, cons7, cons48, cons125, cons208, cons21, cons5, cons379)
    def replacement504(m, p, g, e, x, f, c, a):
        return Int((e*x)**m*(a + c*x**S(2))**p*(f**S(2) + g**S(2)*x**S(2)), x) + S(2)*f*g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p, x)/e
    rule504 = ReplacementRule(pattern504, replacement504)
    rubi.add(rule504.pattern, label = 504)

    pattern505 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(f_ + x_*WC('g', S(1)))**S(3), x_), cons2, cons7, cons48, cons125, cons208, cons21, cons5, cons379)
    def replacement505(m, p, g, e, x, f, c, a):
        return f*Int((e*x)**m*(a + c*x**S(2))**p*(f**S(2) + S(3)*g**S(2)*x**S(2)), x) + g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p*(S(3)*f**S(2) + g**S(2)*x**S(2)), x)/e
    rule505 = ReplacementRule(pattern505, replacement505)
    rubi.add(rule505.pattern, label = 505)

    pattern506 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons226, cons279, cons148)
    def replacement506(m, n, p, g, e, x, b, d, f, c, a):
        return g*Int((d + e*x)**(m + S(1))*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/e
    rule506 = ReplacementRule(pattern506, replacement506)
    rubi.add(rule506.pattern, label = 506)

    pattern507 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons315, cons280, cons148)
    def replacement507(m, n, p, g, e, x, d, f, c, a):
        return g*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(f + g*x)**(n + S(-1)), x)/e + (-d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(-1)), x)/e
    rule507 = ReplacementRule(pattern507, replacement507)
    rubi.add(rule507.pattern, label = 507)

    pattern508 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons380)
    def replacement508(m, n, p, g, e, x, b, d, f, c, a):
        return Int((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x)
    rule508 = ReplacementRule(pattern508, replacement508)
    rubi.add(rule508.pattern, label = 508)

    pattern509 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons381)
    def replacement509(m, n, p, g, e, x, d, f, c, a):
        return Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x)
    rule509 = ReplacementRule(pattern509, replacement509)
    rubi.add(rule509.pattern, label = 509)

    pattern510 = Pattern(Integral((u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(u_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(a_ + u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons68, cons69)
    def replacement510(m, n, p, g, e, x, b, d, f, c, u, a):
        return Subst(Int((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1))
    rule510 = ReplacementRule(pattern510, replacement510)
    rubi.add(rule510.pattern, label = 510)

    pattern511 = Pattern(Integral((a_ + u_**S(2)*WC('c', S(1)))**WC('p', S(1))*(u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(u_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons68, cons69)
    def replacement511(m, n, p, g, e, x, d, f, c, u, a):
        return Subst(Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x, u)/Coefficient(u, x, S(1))
    rule511 = ReplacementRule(pattern511, replacement511)
    rubi.add(rule511.pattern, label = 511)

    pattern512 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons382, cons383, cons384, cons385)
    def replacement512(p, q, d, e, x, b, f, c, a):
        return (c/f)**p*Int((d + e*x + f*x**S(2))**(p + q), x)
    rule512 = ReplacementRule(pattern512, replacement512)
    rubi.add(rule512.pattern, label = 512)

    pattern513 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons382, cons383, cons147, cons386, cons387)
    def replacement513(p, q, d, e, x, b, f, c, a):
        return a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((d + e*x + f*x**S(2))**(p + q), x)
    rule513 = ReplacementRule(pattern513, replacement513)
    rubi.add(rule513.pattern, label = 513)

    pattern514 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons45, cons147)
    def replacement514(q, p, d, e, x, b, f, c, a):
        return (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + e*x + f*x**S(2))**q, x)
    rule514 = ReplacementRule(pattern514, replacement514)
    rubi.add(rule514.pattern, label = 514)

    pattern515 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons125, cons5, cons50, cons45, cons147)
    def replacement515(q, p, d, x, b, f, c, a):
        return (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q, x)
    rule515 = ReplacementRule(pattern515, replacement515)
    rubi.add(rule515.pattern, label = 515)

    pattern516 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons50, cons388, cons389, cons390)
    def replacement516(q, e, x, b, d, f, c, a):
        return (d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) - c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3)))
    rule516 = ReplacementRule(pattern516, replacement516)
    rubi.add(rule516.pattern, label = 516)

    pattern517 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons7, cons27, cons48, cons125, cons50, cons391, cons389, cons390)
    def replacement517(q, e, x, d, f, c, a):
        return (-c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3)))
    rule517 = ReplacementRule(pattern517, replacement517)
    rubi.add(rule517.pattern, label = 517)

    pattern518 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons125, cons50, cons389, cons392)
    def replacement518(q, d, x, b, f, c, a):
        return (d + f*x**S(2))**(q + S(1))*(S(2)*a*f*x*(q + S(1)) + b*d)/(S(2)*d*f*(q + S(1)))
    rule518 = ReplacementRule(pattern518, replacement518)
    rubi.add(rule518.pattern, label = 518)

    pattern519 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons125, cons50, cons393)
    def replacement519(q, d, x, b, f, c, a):
        return b*Int(x*(d + f*x**S(2))**q, x) + Int((a + c*x**S(2))*(d + f*x**S(2))**q, x)
    rule519 = ReplacementRule(pattern519, replacement519)
    rubi.add(rule519.pattern, label = 519)

    pattern520 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons394, cons393)
    def replacement520(q, e, x, b, d, f, c, a):
        return Int(ExpandIntegrand((a + b*x + c*x**S(2))*(d + e*x + f*x**S(2))**q, x), x)
    rule520 = ReplacementRule(pattern520, replacement520)
    rubi.add(rule520.pattern, label = 520)

    pattern521 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons394, cons393)
    def replacement521(q, e, x, d, f, c, a):
        return Int(ExpandIntegrand((a + c*x**S(2))*(d + e*x + f*x**S(2))**q, x), x)
    rule521 = ReplacementRule(pattern521, replacement521)
    rubi.add(rule521.pattern, label = 521)

    pattern522 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons394, cons395, cons396, cons397)
    def replacement522(q, e, x, b, d, f, c, a):
        return -(c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))*Int((d + e*x + f*x**S(2))**(q + S(1)), x)/(f*(q + S(1))*(-S(4)*d*f + e**S(2))) + (d + e*x + f*x**S(2))**(q + S(1))*(a*e*f - S(2)*b*d*f + c*d*e + x*(c*(-S(2)*d*f + e**S(2)) + f*(S(2)*a*f - b*e)))/(f*(q + S(1))*(-S(4)*d*f + e**S(2)))
    rule522 = ReplacementRule(pattern522, replacement522)
    rubi.add(rule522.pattern, label = 522)

    pattern523 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons7, cons27, cons48, cons125, cons394, cons395, cons396, cons398)
    def replacement523(q, e, x, d, f, c, a):
        return -(S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))*Int((d + e*x + f*x**S(2))**(q + S(1)), x)/(f*(q + S(1))*(-S(4)*d*f + e**S(2))) + (d + e*x + f*x**S(2))**(q + S(1))*(a*e*f + c*d*e + x*(S(2)*a*f**S(2) + c*(-S(2)*d*f + e**S(2))))/(f*(q + S(1))*(-S(4)*d*f + e**S(2)))
    rule523 = ReplacementRule(pattern523, replacement523)
    rubi.add(rule523.pattern, label = 523)

    pattern524 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons125, cons395, cons396, cons399)
    def replacement524(q, d, x, b, f, c, a):
        return (d + f*x**S(2))**(q + S(1))*(b*d - x*(a*f - c*d))/(S(2)*d*f*(q + S(1))) + (S(2)*a*f*q + S(3)*a*f - c*d)*Int((d + f*x**S(2))**(q + S(1)), x)/(S(2)*d*f*(q + S(1)))
    rule524 = ReplacementRule(pattern524, replacement524)
    rubi.add(rule524.pattern, label = 524)

    pattern525 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons27, cons48, cons125, cons2, cons3, cons7, cons50, cons394, cons400, cons401, cons397)
    def replacement525(q, e, x, b, d, f, c, a):
        return (c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))*Int((d + e*x + f*x**S(2))**q, x)/(S(2)*f**S(2)*(S(2)*q + S(3))) + (d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) - c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3)))
    rule525 = ReplacementRule(pattern525, replacement525)
    rubi.add(rule525.pattern, label = 525)

    pattern526 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons27, cons48, cons125, cons2, cons7, cons50, cons394, cons400, cons401, cons398)
    def replacement526(q, e, x, d, f, c, a):
        return (S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))*Int((d + e*x + f*x**S(2))**q, x)/(S(2)*f**S(2)*(S(2)*q + S(3))) + (-c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3)))
    rule526 = ReplacementRule(pattern526, replacement526)
    rubi.add(rule526.pattern, label = 526)

    pattern527 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons27, cons125, cons2, cons3, cons7, cons50, cons400, cons401, cons399)
    def replacement527(q, d, x, b, f, c, a):
        return (S(2)*a*f*q + S(3)*a*f - c*d)*Int((d + f*x**S(2))**q, x)/(f*(S(2)*q + S(3))) + (d + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3)))
    rule527 = ReplacementRule(pattern527, replacement527)
    rubi.add(rule527.pattern, label = 527)

    pattern528 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons226, cons394, cons402, cons137, cons403)
    def replacement528(q, p, e, x, b, d, f, c, a):
        return (b + S(2)*c*x)*(a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(b*e*q + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)) + x*(S(2)*b*f*q + S(2)*c*e*(S(2)*p + q + S(3))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule528 = ReplacementRule(pattern528, replacement528)
    rubi.add(rule528.pattern, label = 528)

    pattern529 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons125, cons226, cons402, cons137, cons403)
    def replacement529(q, p, x, b, d, f, c, a):
        return (b + S(2)*c*x)*(d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(S(2)*b*f*q*x + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule529 = ReplacementRule(pattern529, replacement529)
    rubi.add(rule529.pattern, label = 529)

    pattern530 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons394, cons402, cons137, cons403)
    def replacement530(q, p, e, x, d, f, c, a):
        return -x*(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q/(S(2)*a*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*e*x*(S(2)*p + q + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)), x), x)/(S(4)*a*c*(p + S(1)))
    rule530 = ReplacementRule(pattern530, replacement530)
    rubi.add(rule530.pattern, label = 530)

    pattern531 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons50, cons226, cons394, cons13, cons137, cons404, cons405)
    def replacement531(q, p, e, x, b, d, f, c, a):
        return (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(S(2)*a*c**S(2)*e + b**S(3)*f - b**S(2)*c*e + b*c*(-S(3)*a*f + c*d) + c*x*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + S(2)*c*(p + S(1))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) - e*(p + q + S(2))*(-S(2)*a*c**S(2)*e - b**S(3)*f + b**S(2)*c*e - b*c*(-S(3)*a*f + c*d)) + x*(S(2)*f*(p + q + S(2))*(S(2)*a*c**S(2)*e + b**S(3)*f - b**S(2)*c*e + b*c*(-S(3)*a*f + c*d)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e))) - (a*f*(p + S(1)) - c*d*(p + S(2)))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)))
    rule531 = ReplacementRule(pattern531, replacement531)
    rubi.add(rule531.pattern, label = 531)

    pattern532 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons125, cons50, cons226, cons13, cons137, cons406, cons405)
    def replacement532(q, p, x, b, d, f, c, a):
        return (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(b**S(3)*f + b*c*(-S(3)*a*f + c*d) + c*x*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + S(2)*c*(p + S(1))*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + x*(-b*f*(p + S(1))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + S(2)*f*(b**S(3)*f + b*c*(-S(3)*a*f + c*d))*(p + q + S(2))) - (a*f*(p + S(1)) - c*d*(p + S(2)))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2)))
    rule532 = ReplacementRule(pattern532, replacement532)
    rubi.add(rule532.pattern, label = 532)

    pattern533 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons50, cons394, cons13, cons137, cons407, cons405)
    def replacement533(q, p, e, x, d, f, c, a):
        return -(a + c*x**S(2))**(p + S(1))*(S(2)*a*c**S(2)*e + c*x*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(d + e*x + f*x**S(2))**(q + S(1))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(S(2)*a*c**S(2)*e**S(2)*(p + q + S(2)) + c*f*x**S(2)*(-S(2)*a*c*f + S(2)*c**S(2)*d)*(S(2)*p + S(2)*q + S(5)) + S(2)*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + x*(S(4)*a*c**S(2)*e*f*(p + q + S(2)) + c*e*(-S(2)*a*c*f + S(2)*c**S(2)*d)*(S(2)*p + q + S(4))) - (-S(2)*a*c*f + S(2)*c**S(2)*d)*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)))
    rule533 = ReplacementRule(pattern533, replacement533)
    rubi.add(rule533.pattern, label = 533)

    pattern534 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons50, cons226, cons394, cons13, cons146, cons408, cons409)
    def replacement534(q, p, e, x, b, d, f, c, a):
        return (a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(3)*p + S(2)*q) - c*e*(S(2)*p + q) + S(2)*c*f*x*(p + q))/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((a + b*x + c*x**S(2))**(p + S(-2))*(d + e*x + f*x**S(2))**q*Simp(x**S(2)*(c*(p + q)*(-c*(S(2)*d*f*(-S(2)*p + S(1)) + e**S(2)*(S(3)*p + q + S(-1))) + f*(-S(2)*a*f + b*e)*(S(4)*p + S(2)*q + S(-1))) + p*(-p + S(1))*(-b*f + c*e)**S(2)) + x*(S(2)*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d)*(-b*f + c*e) - (p + q)*(b*(c*(S(2)*p + q)*(-S(4)*d*f + e**S(2)) + f*(S(2)*p + S(2)*q + S(1))*(S(2)*a*f - b*e + S(2)*c*d)) + e*f*(-p + S(1))*(-S(4)*a*c + b**S(2)))) + (-p + S(1))*(S(2)*p + q)*(-a*e + b*d)*(-b*f + c*e) - (p + q)*(-a*(c*(S(2)*d*f - e**S(2)*(S(2)*p + q)) + f*(-S(2)*a*f + b*e)*(S(2)*p + S(2)*q + S(1))) + b**S(2)*d*f*(-p + S(1))), x), x)/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1)))
    rule534 = ReplacementRule(pattern534, replacement534)
    rubi.add(rule534.pattern, label = 534)

    pattern535 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons125, cons50, cons226, cons13, cons146, cons408, cons409)
    def replacement535(q, p, x, b, d, f, c, a):
        return (d + f*x**S(2))**(q + S(1))*(b*(S(3)*p + S(2)*q) + S(2)*c*x*(p + q))*(a + b*x + c*x**S(2))**(p + S(-1))/(S(2)*f*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-2))*Simp(b**S(2)*d*(p + S(-1))*(S(2)*p + q) + x**S(2)*(b**S(2)*f*p*(-p + S(1)) + S(2)*c*(p + q)*(-a*f*(S(4)*p + S(2)*q + S(-1)) + c*d*(S(2)*p + S(-1)))) - x*(S(2)*b*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d) - S(2)*b*(p + q)*(S(2)*c*d*(S(2)*p + q) - (a*f + c*d)*(S(2)*p + S(2)*q + S(1)))) - (p + q)*(-S(2)*a*(-a*f*(S(2)*p + S(2)*q + S(1)) + c*d) + b**S(2)*d*(-p + S(1))), x), x)/(S(2)*f*(p + q)*(S(2)*p + S(2)*q + S(1)))
    rule535 = ReplacementRule(pattern535, replacement535)
    rubi.add(rule535.pattern, label = 535)

    pattern536 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons50, cons394, cons13, cons146, cons408, cons409)
    def replacement536(q, p, e, x, d, f, c, a):
        return -c*(a + c*x**S(2))**(p + S(-1))*(e*(S(2)*p + q) - S(2)*f*x*(p + q))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((a + c*x**S(2))**(p + S(-2))*(d + e*x + f*x**S(2))**q*Simp(-a*c*e**S(2)*(-p + S(1))*(S(2)*p + q) + a*(p + q)*(-S(2)*a*f**S(2)*(S(2)*p + S(2)*q + S(1)) + c*(S(2)*d*f - e**S(2)*(S(2)*p + q))) + x**S(2)*(c**S(2)*e**S(2)*p*(-p + S(1)) - c*(p + q)*(S(2)*a*f**S(2)*(S(4)*p + S(2)*q + S(-1)) + c*(S(2)*d*f*(-S(2)*p + S(1)) + e**S(2)*(S(3)*p + q + S(-1))))) + x*(S(4)*a*c*e*f*(-p + S(1))*(p + q) + S(2)*c*e*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d)), x), x)/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1)))
    rule536 = ReplacementRule(pattern536, replacement536)
    rubi.add(rule536.pattern, label = 536)

    def With537(d, e, x, b, f, c, a):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-a*c*f + b**S(2)*f - b*c*e + c**S(2)*d - x*(-b*c*f + c**S(2)*e))/(a + b*x + c*x**S(2)), x)/q + Int((a*f**S(2) - b*e*f - c*d*f + c*e**S(2) + x*(-b*f**S(2) + c*e*f))/(d + e*x + f*x**S(2)), x)/q
        return False
    pattern537 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons226, cons394, CustomConstraint(With537))
    rule537 = ReplacementRule(pattern537, With537)
    rubi.add(rule537.pattern, label = 537)

    def With538(d, x, b, f, c, a):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return -Int((-a*f**S(2) + b*f**S(2)*x + c*d*f)/(d + f*x**S(2)), x)/q + Int((-a*c*f + b**S(2)*f + b*c*f*x + c**S(2)*d)/(a + b*x + c*x**S(2)), x)/q
        return False
    pattern538 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons125, cons226, CustomConstraint(With538))
    rule538 = ReplacementRule(pattern538, With538)
    rubi.add(rule538.pattern, label = 538)

    pattern539 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons226, cons394, cons410)
    def replacement539(e, x, b, d, f, c, a):
        return -S(2)*e*Subst(Int(S(1)/(e*(-S(4)*a*f + b*e) - x**S(2)*(-a*e + b*d)), x), x, (e + S(2)*f*x)/sqrt(d + e*x + f*x**S(2)))
    rule539 = ReplacementRule(pattern539, replacement539)
    rubi.add(rule539.pattern, label = 539)

    def With540(e, x, b, d, f, c, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int(S(1)/((b + S(2)*c*x - q)*sqrt(d + e*x + f*x**S(2))), x)/q - S(2)*c*Int(S(1)/((b + S(2)*c*x + q)*sqrt(d + e*x + f*x**S(2))), x)/q
    pattern540 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons226, cons394, cons411, cons231 )
    rule540 = ReplacementRule(pattern540, With540)
    rubi.add(rule540.pattern, label = 540)

    pattern541 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons394, cons412)
    def replacement541(e, x, d, f, c, a):
        return Int(S(1)/((a - x*Rt(-a*c, S(2)))*sqrt(d + e*x + f*x**S(2))), x)/S(2) + Int(S(1)/((a + x*Rt(-a*c, S(2)))*sqrt(d + e*x + f*x**S(2))), x)/S(2)
    rule541 = ReplacementRule(pattern541, replacement541)
    rubi.add(rule541.pattern, label = 541)

    def With542(d, x, b, f, c, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x - q)), x)/q - S(2)*c*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x + q)), x)/q
    pattern542 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons125, cons226, cons231 )
    rule542 = ReplacementRule(pattern542, With542)
    rubi.add(rule542.pattern, label = 542)

    def With543(e, x, b, d, f, c, a):
        q = Rt(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f + c*d - q + x*(-b*f + c*e))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) + Int((-a*f + c*d + q + x*(-b*f + c*e))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    pattern543 = Pattern(Integral(S(1)/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons226, cons394, cons411, cons413 )
    rule543 = ReplacementRule(pattern543, With543)
    rubi.add(rule543.pattern, label = 543)

    def With544(e, x, d, f, c, a):
        q = Rt(a*c*e**S(2) + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f + c*d + c*e*x - q)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) + Int((-a*f + c*d + c*e*x + q)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    pattern544 = Pattern(Integral(S(1)/((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons394, cons414 )
    rule544 = ReplacementRule(pattern544, With544)
    rubi.add(rule544.pattern, label = 544)

    def With545(x, b, d, f, c, a):
        q = Rt(b**S(2)*d*f + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f - b*f*x + c*d - q)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q) + Int((-a*f - b*f*x + c*d + q)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q)
    pattern545 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons125, cons226, cons413 )
    rule545 = ReplacementRule(pattern545, With545)
    rubi.add(rule545.pattern, label = 545)

    pattern546 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons226, cons394)
    def replacement546(d, e, x, b, f, c, a):
        return c*Int(S(1)/sqrt(a + b*x + c*x**S(2)), x)/f - Int((-a*f + c*d + x*(-b*f + c*e))/(sqrt(a + b*x + c*x**S(2))*(d + e*x + f*x**S(2))), x)/f
    rule546 = ReplacementRule(pattern546, replacement546)
    rubi.add(rule546.pattern, label = 546)

    pattern547 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(d_ + x_**S(2)*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons125, cons226)
    def replacement547(d, x, b, f, c, a):
        return c*Int(S(1)/sqrt(a + b*x + c*x**S(2)), x)/f - Int((-a*f - b*f*x + c*d)/((d + f*x**S(2))*sqrt(a + b*x + c*x**S(2))), x)/f
    rule547 = ReplacementRule(pattern547, replacement547)
    rubi.add(rule547.pattern, label = 547)

    pattern548 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)))/(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons394)
    def replacement548(d, x, e, f, c, a):
        return c*Int(S(1)/sqrt(a + c*x**S(2)), x)/f - Int((-a*f + c*d + c*e*x)/(sqrt(a + c*x**S(2))*(d + e*x + f*x**S(2))), x)/f
    rule548 = ReplacementRule(pattern548, replacement548)
    rubi.add(rule548.pattern, label = 548)

    def With549(d, e, x, b, f, c, a):
        r = Rt(-S(4)*a*c + b**S(2), S(2))
        return sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*Int(S(1)/(sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*sqrt(d + e*x + f*x**S(2))), x)/sqrt(a + b*x + c*x**S(2))
    pattern549 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons226, cons394 )
    rule549 = ReplacementRule(pattern549, With549)
    rubi.add(rule549.pattern, label = 549)

    def With550(d, x, b, f, c, a):
        r = Rt(-S(4)*a*c + b**S(2), S(2))
        return sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*Int(S(1)/(sqrt(S(2)*a + x*(b + r))*sqrt(d + f*x**S(2))*sqrt(b + S(2)*c*x + r)), x)/sqrt(a + b*x + c*x**S(2))
    pattern550 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons125, cons226 )
    rule550 = ReplacementRule(pattern550, With550)
    rubi.add(rule550.pattern, label = 550)

    pattern551 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons415)
    def replacement551(p, q, e, x, b, d, f, c, a):
        return Int((a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)
    rule551 = ReplacementRule(pattern551, replacement551)
    rubi.add(rule551.pattern, label = 551)

    pattern552 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons7, cons27, cons48, cons125, cons5, cons50, cons416)
    def replacement552(p, q, e, x, d, f, c, a):
        return Int((a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)
    rule552 = ReplacementRule(pattern552, replacement552)
    rubi.add(rule552.pattern, label = 552)

    pattern553 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons68, cons69)
    def replacement553(p, q, e, x, b, d, f, c, u, a):
        return Subst(Int((a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1))
    rule553 = ReplacementRule(pattern553, replacement553)
    rubi.add(rule553.pattern, label = 553)

    pattern554 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons5, cons50, cons68, cons69)
    def replacement554(p, q, e, x, d, f, c, u, a):
        return Subst(Int((a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1))
    rule554 = ReplacementRule(pattern554, replacement554)
    rubi.add(rule554.pattern, label = 554)

    pattern555 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons382, cons383, cons384, cons385)
    def replacement555(m, h, p, q, d, g, e, x, b, f, c, a):
        return (c/f)**p*Int((g + h*x)**m*(d + e*x + f*x**S(2))**(p + q), x)
    rule555 = ReplacementRule(pattern555, replacement555)
    rubi.add(rule555.pattern, label = 555)

    pattern556 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons382, cons383, cons147, cons386, cons387)
    def replacement556(m, h, p, q, d, g, e, x, b, f, c, a):
        return a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((g + h*x)**m*(d + e*x + f*x**S(2))**(p + q), x)
    rule556 = ReplacementRule(pattern556, replacement556)
    rubi.add(rule556.pattern, label = 556)

    pattern557 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons50, cons45)
    def replacement557(m, h, q, p, d, g, e, x, b, f, c, a):
        return (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(g + h*x)**m*(d + e*x + f*x**S(2))**q, x)
    rule557 = ReplacementRule(pattern557, replacement557)
    rubi.add(rule557.pattern, label = 557)

    pattern558 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons21, cons5, cons50, cons45)
    def replacement558(m, h, q, p, d, g, x, b, f, c, a):
        return (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q*(g + h*x)**m, x)
    rule558 = ReplacementRule(pattern558, replacement558)
    rubi.add(rule558.pattern, label = 558)

    pattern559 = Pattern(Integral((g_ + x_*WC('h', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons417, cons418, cons17)
    def replacement559(m, h, p, g, e, x, b, d, f, c, a):
        return Int((f*h*x/c + d*g/a)**m*(a + b*x + c*x**S(2))**(m + p), x)
    rule559 = ReplacementRule(pattern559, replacement559)
    rubi.add(rule559.pattern, label = 559)

    pattern560 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons419, cons418, cons17)
    def replacement560(m, h, p, g, e, x, d, f, c, a):
        return Int((a + c*x**S(2))**(m + p)*(f*h*x/c + d*g/a)**m, x)
    rule560 = ReplacementRule(pattern560, replacement560)
    rubi.add(rule560.pattern, label = 560)

    pattern561 = Pattern(Integral((g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons5, cons417, cons420, cons17)
    def replacement561(m, h, p, g, x, b, d, f, c, a):
        return Int((f*h*x/c + d*g/a)**m*(a + b*x + c*x**S(2))**(m + p), x)
    rule561 = ReplacementRule(pattern561, replacement561)
    rubi.add(rule561.pattern, label = 561)

    pattern562 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons5, cons419, cons420, cons17)
    def replacement562(m, h, p, g, x, d, f, c, a):
        return Int((a + c*x**S(2))**(m + p)*(f*h*x/c + d*g/a)**m, x)
    rule562 = ReplacementRule(pattern562, replacement562)
    rubi.add(rule562.pattern, label = 562)

    pattern563 = Pattern(Integral(x_**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons48, cons125, cons50, cons226, cons421, cons38)
    def replacement563(p, q, e, x, b, f, c, a):
        return Int((a/e + c*x/f)**p*(e*x + f*x**S(2))**(p + q), x)
    rule563 = ReplacementRule(pattern563, replacement563)
    rubi.add(rule563.pattern, label = 563)

    pattern564 = Pattern(Integral(x_**WC('p', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons48, cons125, cons50, cons422, cons38)
    def replacement564(p, q, e, x, f, c, a):
        return Int((a/e + c*x/f)**p*(e*x + f*x**S(2))**(p + q), x)
    rule564 = ReplacementRule(pattern564, replacement564)
    rubi.add(rule564.pattern, label = 564)

    pattern565 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons423, cons424, cons270)
    def replacement565(m, h, p, g, e, x, b, d, f, c, a):
        return f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3)))
    rule565 = ReplacementRule(pattern565, replacement565)
    rubi.add(rule565.pattern, label = 565)

    pattern566 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons425, cons426, cons270)
    def replacement566(m, h, p, g, e, x, d, f, c, a):
        return f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3)))
    rule566 = ReplacementRule(pattern566, replacement566)
    rubi.add(rule566.pattern, label = 566)

    pattern567 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons21, cons5, cons427, cons424, cons270)
    def replacement567(m, h, p, d, g, x, b, f, c, a):
        return f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3)))
    rule567 = ReplacementRule(pattern567, replacement567)
    rubi.add(rule567.pattern, label = 567)

    pattern568 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons226, cons394, cons428)
    def replacement568(m, h, p, g, e, x, b, d, f, c, a):
        return Int(ExpandIntegrand((g + h*x)**m*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2)), x), x)
    rule568 = ReplacementRule(pattern568, replacement568)
    rubi.add(rule568.pattern, label = 568)

    pattern569 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons394, cons428)
    def replacement569(m, h, p, g, e, x, d, f, c, a):
        return Int(ExpandIntegrand((a + c*x**S(2))**p*(g + h*x)**m*(d + e*x + f*x**S(2)), x), x)
    rule569 = ReplacementRule(pattern569, replacement569)
    rubi.add(rule569.pattern, label = 569)

    pattern570 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons21, cons226, cons428)
    def replacement570(m, h, p, d, g, x, b, f, c, a):
        return Int(ExpandIntegrand((d + f*x**S(2))*(g + h*x)**m*(a + b*x + c*x**S(2))**p, x), x)
    rule570 = ReplacementRule(pattern570, replacement570)
    rubi.add(rule570.pattern, label = 570)

    pattern571 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons21, cons428)
    def replacement571(m, h, p, g, x, d, f, c, a):
        return Int(ExpandIntegrand((a + c*x**S(2))**p*(d + f*x**S(2))*(g + h*x)**m, x), x)
    rule571 = ReplacementRule(pattern571, replacement571)
    rubi.add(rule571.pattern, label = 571)

    pattern572 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons226, cons394, cons31, cons94, cons429)
    def replacement572(m, h, p, g, e, x, b, d, f, c, a):
        return (g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(d*h**S(2) - e*g*h + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))) + Int((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*(f*g**S(2)*(p + S(1)) - h*(-d*h*(m + p + S(2)) + e*g*(p + S(1)))) + h*(m + S(1))*(a*e*h - a*f*g + c*d*g) - x*(c*(S(2)*f*g**S(2)*(p + S(1)) - h*(-d*h + e*g)*(m + S(2)*p + S(3))) + f*h*(m + S(1))*(-a*h + b*g)), x), x)/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2)))
    rule572 = ReplacementRule(pattern572, replacement572)
    rubi.add(rule572.pattern, label = 572)

    pattern573 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons394, cons31, cons94, cons430)
    def replacement573(m, h, p, g, e, x, d, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(d*h**S(2) - e*g*h + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**p*(g + h*x)**(m + S(1))*Simp(h*(m + S(1))*(a*e*h - a*f*g + c*d*g) + x*(a*f*h**S(2)*(m + S(1)) - c*(S(2)*f*g**S(2)*(p + S(1)) - h*(-d*h + e*g)*(m + S(2)*p + S(3)))), x), x)/(h*(m + S(1))*(a*h**S(2) + c*g**S(2)))
    rule573 = ReplacementRule(pattern573, replacement573)
    rubi.add(rule573.pattern, label = 573)

    pattern574 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons5, cons226, cons31, cons94, cons429)
    def replacement574(m, h, p, g, x, b, d, f, c, a):
        return (g + h*x)**(m + S(1))*(d*h**S(2) + f*g**S(2))*(a + b*x + c*x**S(2))**(p + S(1))/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))) + Int((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*(d*h**S(2)*(m + p + S(2)) + f*g**S(2)*(p + S(1))) + h*(m + S(1))*(-a*f*g + c*d*g) - x*(c*(d*h**S(2)*(m + S(2)*p + S(3)) + S(2)*f*g**S(2)*(p + S(1))) + f*h*(m + S(1))*(-a*h + b*g)), x), x)/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2)))
    rule574 = ReplacementRule(pattern574, replacement574)
    rubi.add(rule574.pattern, label = 574)

    pattern575 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons5, cons31, cons94, cons430)
    def replacement575(m, h, p, d, g, x, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(d*h**S(2) + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**p*(g + h*x)**(m + S(1))*Simp(h*(m + S(1))*(-a*f*g + c*d*g) + x*(a*f*h**S(2)*(m + S(1)) - c*(d*h**S(2)*(m + S(2)*p + S(3)) + S(2)*f*g**S(2)*(p + S(1)))), x), x)/(h*(m + S(1))*(a*h**S(2) + c*g**S(2)))
    rule575 = ReplacementRule(pattern575, replacement575)
    rubi.add(rule575.pattern, label = 575)

    pattern576 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, cons429)
    def replacement576(h, g, e, x, b, d, f, c, a):
        return (d*h**S(2) - e*g*h + f*g**S(2))*Int(S(1)/((g + h*x)*sqrt(a + b*x + c*x**S(2))), x)/(a*h**S(2) - b*g*h + c*g**S(2)) + Int((a*e*h - a*f*g - b*d*h + c*d*g + x*(a*f*h - b*f*g - c*d*h + c*e*g))/(a + b*x + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) - b*g*h + c*g**S(2))
    rule576 = ReplacementRule(pattern576, replacement576)
    rubi.add(rule576.pattern, label = 576)

    pattern577 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2)*(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons394, cons430)
    def replacement577(h, g, e, x, d, f, c, a):
        return (d*h**S(2) - e*g*h + f*g**S(2))*Int(S(1)/(sqrt(a + c*x**S(2))*(g + h*x)), x)/(a*h**S(2) + c*g**S(2)) + Int((a*e*h - a*f*g + c*d*g + x*(a*f*h - c*d*h + c*e*g))/(a + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) + c*g**S(2))
    rule577 = ReplacementRule(pattern577, replacement577)
    rubi.add(rule577.pattern, label = 577)

    pattern578 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons226, cons429)
    def replacement578(h, g, x, b, d, f, c, a):
        return (d*h**S(2) + f*g**S(2))*Int(S(1)/((g + h*x)*sqrt(a + b*x + c*x**S(2))), x)/(a*h**S(2) - b*g*h + c*g**S(2)) + Int((-a*f*g - b*d*h + c*d*g - x*(-a*f*h + b*f*g + c*d*h))/(a + b*x + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) - b*g*h + c*g**S(2))
    rule578 = ReplacementRule(pattern578, replacement578)
    rubi.add(rule578.pattern, label = 578)

    pattern579 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2)*(g_ + x_*WC('h', S(1)))), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons430)
    def replacement579(h, d, g, x, f, c, a):
        return (d*h**S(2) + f*g**S(2))*Int(S(1)/(sqrt(a + c*x**S(2))*(g + h*x)), x)/(a*h**S(2) + c*g**S(2)) + Int((-a*f*g + c*d*g - x*(-a*f*h + c*d*h))/(a + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) + c*g**S(2))
    rule579 = ReplacementRule(pattern579, replacement579)
    rubi.add(rule579.pattern, label = 579)

    pattern580 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, cons244, cons137, cons166)
    def replacement580(m, h, p, g, e, x, b, d, f, c, a):
        return (g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(a*b*f - S(2)*a*c*e + b*c*d + x*(c*(-b*e + S(2)*c*d) + f*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((g + h*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(c*(S(2)*p + S(3))*(-b*e + S(2)*c*d) - f*(S(2)*a*c - b**S(2)*(p + S(2)))) + h*m*(a*b*f - S(2)*a*c*e + b*c*d) + h*x*(c*(-b*e + S(2)*c*d)*(m + S(2)*p + S(3)) - f*(S(2)*a*c*(m + S(1)) - b**S(2)*(m + p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule580 = ReplacementRule(pattern580, replacement580)
    rubi.add(rule580.pattern, label = 580)

    pattern581 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons394, cons244, cons137, cons166)
    def replacement581(m, h, p, g, e, x, d, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(g + h*x)**m*(a*e - x*(-a*f + c*d))/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(-1))*Simp(a*(e*h*m + f*g) - c*d*g*(S(2)*p + S(3)) + h*x*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1)))
    rule581 = ReplacementRule(pattern581, replacement581)
    rubi.add(rule581.pattern, label = 581)

    pattern582 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons226, cons244, cons137, cons166)
    def replacement582(m, h, p, g, x, b, d, f, c, a):
        return (g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(a*b*f + b*c*d + x*(S(2)*c**S(2)*d + f*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((g + h*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(S(2)*c**S(2)*d*(S(2)*p + S(3)) - f*(S(2)*a*c - b**S(2)*(p + S(2)))) + h*m*(a*b*f + b*c*d) + h*x*(S(2)*c**S(2)*d*(m + S(2)*p + S(3)) - f*(S(2)*a*c*(m + S(1)) - b**S(2)*(m + p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule582 = ReplacementRule(pattern582, replacement582)
    rubi.add(rule582.pattern, label = 582)

    pattern583 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(g_ + x_*WC('h', S(1)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons244, cons137, cons166)
    def replacement583(m, h, p, d, g, x, f, c, a):
        return -x*(a + c*x**S(2))**(p + S(1))*(g + h*x)**m*(-a*f + c*d)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(-1))*Simp(a*f*g - c*d*g*(S(2)*p + S(3)) + h*x*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1)))
    rule583 = ReplacementRule(pattern583, replacement583)
    rubi.add(rule583.pattern, label = 583)

    pattern584 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons226, cons394, cons13, cons137, cons431)
    def replacement584(m, h, p, g, e, x, b, d, f, c, a):
        return -(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-x*(b*f*(-a*h + b*g) + S(2)*c**S(2)*d*g - c*(-S(2)*a*e*h + S(2)*a*f*g + b*d*h + b*e*g)) - (-a*e + b*d)*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(p + S(2))*(-S(2)*a*(-c*e*h + c*f*g) + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c**S(2)*d*g) + h*x*(m + S(2)*p + S(4))*(-S(2)*a*(-c*e*h + c*f*g) + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c**S(2)*d*g) - h*(-(-a*e + b*d)*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))*(m + p + S(2)) + (p + S(1))*(c*g**S(2) - h*(-a*h + b*g))*(S(2)*a*f - b*e + S(2)*c*d), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g)))
    rule584 = ReplacementRule(pattern584, replacement584)
    rubi.add(rule584.pattern, label = 584)

    pattern585 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons394, cons13, cons137, cons430)
    def replacement585(m, h, p, g, e, x, d, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(a*c*e*g - a*h*(-a*f + c*d) - c*x*(a*e*h - a*f*g + c*d*g))/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**m*Simp(g*(p + S(2))*(-a*(-c*e*h + c*f*g) + c**S(2)*d*g) + h*x*(-a*(-c*e*h + c*f*g) + c**S(2)*d*g)*(m + S(2)*p + S(4)) - h*(a*c*e*g - a*h*(-a*f + c*d))*(m + p + S(2)) + (p + S(1))*(a*f + c*d)*(a*h**S(2) + c*g**S(2)), x), x)/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2)))
    rule585 = ReplacementRule(pattern585, replacement585)
    rubi.add(rule585.pattern, label = 585)

    pattern586 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons21, cons226, cons13, cons137, cons431)
    def replacement586(m, h, p, g, x, b, d, f, c, a):
        return -(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*d*(-b*h + S(2)*c*g) - x*(b*f*(-a*h + b*g) + S(2)*c**S(2)*d*g - c*(S(2)*a*f*g + b*d*h)) + (-a*f + c*d)*(-S(2)*a*h + b*g))/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(p + S(2))*(-S(2)*a*c*f*g + b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c**S(2)*d*g) + h*x*(m + S(2)*p + S(4))*(-S(2)*a*c*f*g + b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c**S(2)*d*g) - h*(-b*d*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))*(m + p + S(2)) + S(2)*(p + S(1))*(a*f + c*d)*(c*g**S(2) - h*(-a*h + b*g)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g)))
    rule586 = ReplacementRule(pattern586, replacement586)
    rubi.add(rule586.pattern, label = 586)

    pattern587 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons21, cons13, cons137, cons430)
    def replacement587(m, h, p, d, g, x, f, c, a):
        return -(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(a*h*(-a*f + c*d) + c*x*(-a*f*g + c*d*g))/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**m*Simp(a*h**S(2)*(-a*f + c*d)*(m + p + S(2)) + g*(p + S(2))*(-a*c*f*g + c**S(2)*d*g) + h*x*(-a*c*f*g + c**S(2)*d*g)*(m + S(2)*p + S(4)) + (p + S(1))*(a*f + c*d)*(a*h**S(2) + c*g**S(2)), x), x)/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2)))
    rule587 = ReplacementRule(pattern587, replacement587)
    rubi.add(rule587.pattern, label = 587)

    pattern588 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons226, cons394, cons242)
    def replacement588(m, h, p, g, e, x, b, d, f, c, a):
        return f*Int((g + h*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/h**S(2) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(-d*h**S(2) + f*g**S(2) + h*x*(-e*h + S(2)*f*g)), x)/h**S(2)
    rule588 = ReplacementRule(pattern588, replacement588)
    rubi.add(rule588.pattern, label = 588)

    pattern589 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons394, cons242)
    def replacement589(m, h, p, g, e, x, d, f, c, a):
        return f*Int((a + c*x**S(2))**p*(g + h*x)**(m + S(2)), x)/h**S(2) - Int((a + c*x**S(2))**p*(g + h*x)**m*(-d*h**S(2) + f*g**S(2) + h*x*(-e*h + S(2)*f*g)), x)/h**S(2)
    rule589 = ReplacementRule(pattern589, replacement589)
    rubi.add(rule589.pattern, label = 589)

    pattern590 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons21, cons5, cons226, cons242)
    def replacement590(m, h, p, g, x, b, d, f, c, a):
        return f*Int((g + h*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/h**S(2) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(-d*h**S(2) + f*g**S(2) + S(2)*f*g*h*x), x)/h**S(2)
    rule590 = ReplacementRule(pattern590, replacement590)
    rubi.add(rule590.pattern, label = 590)

    pattern591 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons21, cons5, cons242)
    def replacement591(m, h, p, g, x, d, f, c, a):
        return f*Int((a + c*x**S(2))**p*(g + h*x)**(m + S(2)), x)/h**S(2) - Int((a + c*x**S(2))**p*(g + h*x)**m*(-d*h**S(2) + f*g**S(2) + S(2)*f*g*h*x), x)/h**S(2)
    rule591 = ReplacementRule(pattern591, replacement591)
    rubi.add(rule591.pattern, label = 591)

    pattern592 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons226, cons394, cons270)
    def replacement592(m, h, p, g, e, x, b, d, f, c, a):
        return f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*f*g*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))) + x*(b*f*h*(m + p + S(2)) + c*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1)))), x), x)/(c*h*(m + S(2)*p + S(3)))
    rule592 = ReplacementRule(pattern592, replacement592)
    rubi.add(rule592.pattern, label = 592)

    pattern593 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons394, cons270)
    def replacement593(m, h, p, g, e, x, d, f, c, a):
        return f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(g + h*x)**m*Simp(c*x*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1))) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3)))
    rule593 = ReplacementRule(pattern593, replacement593)
    rubi.add(rule593.pattern, label = 593)

    pattern594 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons21, cons5, cons226, cons270)
    def replacement594(m, h, p, d, g, x, b, f, c, a):
        return f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*f*g*(p + S(1)) + f*x*(b*h*(m + p + S(2)) + S(2)*c*g*(p + S(1))) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3)))
    rule594 = ReplacementRule(pattern594, replacement594)
    rubi.add(rule594.pattern, label = 594)

    pattern595 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)))*(g_ + x_*WC('h', S(1)))**WC('m', S(1)), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons21, cons5, cons270)
    def replacement595(m, h, p, d, g, x, f, c, a):
        return f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(g + h*x)**m*Simp(S(2)*c*f*g*x*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3)))
    rule595 = ReplacementRule(pattern595, replacement595)
    rubi.add(rule595.pattern, label = 595)

    pattern596 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, cons220, cons163)
    def replacement596(h, p, q, d, g, e, x, b, f, c, a):
        return Int(ExpandIntegrand((g + h*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x)
    rule596 = ReplacementRule(pattern596, replacement596)
    rubi.add(rule596.pattern, label = 596)

    pattern597 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons394, cons220, cons432)
    def replacement597(h, p, q, d, g, e, x, f, c, a):
        return Int(ExpandIntegrand((a + c*x**S(2))**p*(g + h*x)*(d + e*x + f*x**S(2))**q, x), x)
    rule597 = ReplacementRule(pattern597, replacement597)
    rubi.add(rule597.pattern, label = 597)

    pattern598 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, cons402, cons137, cons403)
    def replacement598(h, p, q, d, g, e, x, b, f, c, a):
        return (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-d*(S(2)*p + S(3))*(b*h - S(2)*c*g) + e*q*(-S(2)*a*h + b*g) - f*x**S(2)*(b*h - S(2)*c*g)*(S(2)*p + S(2)*q + S(3)) + x*(-e*(b*h - S(2)*c*g)*(S(2)*p + q + S(3)) + S(2)*f*q*(-S(2)*a*h + b*g)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule598 = ReplacementRule(pattern598, replacement598)
    rubi.add(rule598.pattern, label = 598)

    pattern599 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons394, cons402, cons137, cons403)
    def replacement599(h, p, q, d, g, e, x, f, c, a):
        return (a + c*x**S(2))**(p + S(1))*(a*h - c*g*x)*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-a*e*h*q + c*d*g*(S(2)*p + S(3)) + c*f*g*x**S(2)*(S(2)*p + S(2)*q + S(3)) + x*(-S(2)*a*f*h*q + c*e*g*(S(2)*p + q + S(3))), x), x)/(S(2)*a*c*(p + S(1)))
    rule599 = ReplacementRule(pattern599, replacement599)
    rubi.add(rule599.pattern, label = 599)

    pattern600 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons226, cons402, cons137, cons403)
    def replacement600(h, p, q, d, g, x, b, f, c, a):
        return (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-d*(S(2)*p + S(3))*(b*h - S(2)*c*g) + S(2)*f*q*x*(-S(2)*a*h + b*g) - f*x**S(2)*(b*h - S(2)*c*g)*(S(2)*p + S(2)*q + S(3)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2)))
    rule600 = ReplacementRule(pattern600, replacement600)
    rubi.add(rule600.pattern, label = 600)

    pattern601 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons50, cons226, cons394, cons13, cons137, cons404, cons405)
    def replacement601(h, p, q, d, g, e, x, b, f, c, a):
        return (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + c*x*(g*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) - h*(a*b*f - S(2)*a*c*e + b*c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d)) - e*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))*(p + q + S(2)) - x*(S(2)*f*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d))) + (p + S(1))*(b*h - S(2)*c*g)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)))
    rule601 = ReplacementRule(pattern601, replacement601)
    rubi.add(rule601.pattern, label = 601)

    pattern602 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons50, cons394, cons13, cons137, cons407, cons405)
    def replacement602(h, p, q, d, g, e, x, f, c, a):
        return -(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d) + c*x*(S(2)*a*c*e*h + g*(-S(2)*a*c*f + S(2)*c**S(2)*d)))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d))*(S(2)*p + S(2)*q + S(5)) - S(2)*c*g*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)) - e*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(p + q + S(2)) - x*(c*e*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d))*(S(2)*p + q + S(4)) + S(2)*f*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d)), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)))
    rule602 = ReplacementRule(pattern602, replacement602)
    rubi.add(rule602.pattern, label = 602)

    pattern603 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons50, cons226, cons13, cons137, cons406, cons405)
    def replacement603(h, p, q, d, g, x, b, f, c, a):
        return (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*c*g*(a*f + c*d) + c*x*(g*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - h*(a*b*f + b*c*d)) + (-a*h + b*g)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)) - x*(-b*f*(p + S(1))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)) + S(2)*f*(-b*c*g*(a*f + c*d) + (-a*h + b*g)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (p + S(1))*(b*h - S(2)*c*g)*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2)))
    rule603 = ReplacementRule(pattern603, replacement603)
    rubi.add(rule603.pattern, label = 603)

    pattern604 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons50, cons226, cons394, cons13, cons163, cons433)
    def replacement604(h, p, q, d, g, e, x, b, f, c, a):
        return h*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f*(p + q + S(1))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*(e*h - S(2)*f*g)*(p + q + S(1)) + h*p*(-a*e + b*d) + x**S(2)*(c*(e*h - S(2)*f*g)*(p + q + S(1)) + h*p*(-b*f + c*e)) + x*(b*(e*h - S(2)*f*g)*(p + q + S(1)) + S(2)*h*p*(-a*f + c*d)), x), x)/(S(2)*f*(p + q + S(1)))
    rule604 = ReplacementRule(pattern604, replacement604)
    rubi.add(rule604.pattern, label = 604)

    pattern605 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons50, cons394, cons13, cons163, cons433)
    def replacement605(h, p, q, d, g, e, x, f, c, a):
        return h*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f*(p + q + S(1))) + Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*e*h*p - a*(e*h - S(2)*f*g)*(p + q + S(1)) - S(2)*h*p*x*(-a*f + c*d) - x**S(2)*(c*e*h*p + c*(e*h - S(2)*f*g)*(p + q + S(1))), x), x)/(S(2)*f*(p + q + S(1)))
    rule605 = ReplacementRule(pattern605, replacement605)
    rubi.add(rule605.pattern, label = 605)

    pattern606 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons50, cons226, cons13, cons163, cons433)
    def replacement606(h, p, q, d, g, x, b, f, c, a):
        return h*(d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**p/(S(2)*f*(p + q + S(1))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-S(2)*a*f*g*(p + q + S(1)) + b*d*h*p + x**S(2)*(-b*f*h*p - S(2)*c*f*g*(p + q + S(1))) + x*(-S(2)*b*f*g*(p + q + S(1)) + S(2)*h*p*(-a*f + c*d)), x), x)/(S(2)*f*(p + q + S(1)))
    rule606 = ReplacementRule(pattern606, replacement606)
    rubi.add(rule606.pattern, label = 606)

    def With607(h, d, g, e, x, b, f, c, a):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int(Simp(-a*b*f*h + a*c*e*h - a*c*f*g + b**S(2)*f*g - b*c*e*g + c**S(2)*d*g + c*x*(-a*f*h + b*f*g + c*d*h - c*e*g), x)/(a + b*x + c*x**S(2)), x)/q + Int(Simp(a*f**S(2)*g + b*d*f*h - b*e*f*g - c*d*e*h - c*d*f*g + c*e**S(2)*g - f*x*(-a*f*h + b*f*g + c*d*h - c*e*g), x)/(d + e*x + f*x**S(2)), x)/q
        return False
    pattern607 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, CustomConstraint(With607))
    rule607 = ReplacementRule(pattern607, With607)
    rubi.add(rule607.pattern, label = 607)

    def With608(h, d, g, x, b, f, c, a):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int(Simp(a*f**S(2)*g + b*d*f*h - c*d*f*g - f*x*(-a*f*h + b*f*g + c*d*h), x)/(d + f*x**S(2)), x)/q + Int(Simp(-a*b*f*h - a*c*f*g + b**S(2)*f*g + c**S(2)*d*g + c*x*(-a*f*h + b*f*g + c*d*h), x)/(a + b*x + c*x**S(2)), x)/q
        return False
    pattern608 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons226, CustomConstraint(With608))
    rule608 = ReplacementRule(pattern608, With608)
    rubi.add(rule608.pattern, label = 608)

    pattern609 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons434)
    def replacement609(h, d, g, x, f, c, a):
        return g*Int(S(1)/((a + c*x**S(2))*sqrt(d + f*x**S(2))), x) + h*Int(x/((a + c*x**S(2))*sqrt(d + f*x**S(2))), x)
    rule609 = ReplacementRule(pattern609, replacement609)
    rubi.add(rule609.pattern, label = 609)

    def With610(h, d, g, x, f, c, a):
        q = Rt(-a*c, S(2))
        return -(c*g - h*q)*Int(S(1)/(sqrt(d + f*x**S(2))*(c*x + q)), x)/(S(2)*q) - (c*g + h*q)*Int(S(1)/(sqrt(d + f*x**S(2))*(-c*x + q)), x)/(S(2)*q)
    pattern610 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons435 )
    rule610 = ReplacementRule(pattern610, With610)
    rubi.add(rule610.pattern, label = 610)

    pattern611 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, cons410, cons436)
    def replacement611(h, g, e, x, b, d, f, c, a):
        return -S(2)*g*Subst(Int(S(1)/(-a*e + b*d - b*x**S(2)), x), x, sqrt(d + e*x + f*x**S(2)))
    rule611 = ReplacementRule(pattern611, replacement611)
    rubi.add(rule611.pattern, label = 611)

    pattern612 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, cons410, cons437)
    def replacement612(h, g, e, x, b, d, f, c, a):
        return h*Int((e + S(2)*f*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*f) - (e*h - S(2)*f*g)*Int(S(1)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*f)
    rule612 = ReplacementRule(pattern612, replacement612)
    rubi.add(rule612.pattern, label = 612)

    pattern613 = Pattern(Integral(x_/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons226, cons394, cons383)
    def replacement613(d, x, e, b, f, c, a):
        return -S(2)*e*Subst(Int((-d*x**S(2) + S(1))/(-b*f + c*e + d**S(2)*x**S(4)*(-b*f + c*e) - e*x**S(2)*(S(2)*a*f - b*e + S(2)*c*d)), x), x, (S(1) + x*(e + sqrt(-S(4)*d*f + e**S(2)))/(S(2)*d))/sqrt(d + e*x + f*x**S(2)))
    rule613 = ReplacementRule(pattern613, replacement613)
    rubi.add(rule613.pattern, label = 613)

    pattern614 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, cons383, cons438)
    def replacement614(h, d, g, e, x, b, f, c, a):
        return g*Subst(Int(S(1)/(a + x**S(2)*(-a*f + c*d)), x), x, x/sqrt(d + e*x + f*x**S(2)))
    rule614 = ReplacementRule(pattern614, replacement614)
    rubi.add(rule614.pattern, label = 614)

    pattern615 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, cons383, cons439)
    def replacement615(h, d, g, e, x, b, f, c, a):
        return h*Int((S(2)*d + e*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/e - (S(2)*d*h - e*g)*Int(S(1)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/e
    rule615 = ReplacementRule(pattern615, replacement615)
    rubi.add(rule615.pattern, label = 615)

    pattern616 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, cons372, cons440)
    def replacement616(h, g, e, x, b, d, f, c, a):
        return -S(2)*g*(-S(2)*a*h + b*g)*Subst(Int(S(1)/Simp(g*(-S(4)*a*c + b**S(2))*(-S(2)*a*h + b*g) - x**S(2)*(-a*e + b*d), x), x), x, Simp(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g), x)/sqrt(d + e*x + f*x**S(2)))
    rule616 = ReplacementRule(pattern616, replacement616)
    rubi.add(rule616.pattern, label = 616)

    pattern617 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons441)
    def replacement617(h, d, g, e, x, f, c, a):
        return -S(2)*a*g*h*Subst(Int(S(1)/Simp(S(2)*a**S(2)*c*g*h + a*e*x**S(2), x), x), x, Simp(a*h - c*g*x, x)/sqrt(d + e*x + f*x**S(2)))
    rule617 = ReplacementRule(pattern617, replacement617)
    rubi.add(rule617.pattern, label = 617)

    pattern618 = Pattern(Integral((g_ + x_*WC('h', S(1)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons226, cons442)
    def replacement618(h, d, g, x, b, f, c, a):
        return -S(2)*g*(-S(2)*a*h + b*g)*Subst(Int(S(1)/Simp(-b*d*x**S(2) + g*(-S(4)*a*c + b**S(2))*(-S(2)*a*h + b*g), x), x), x, Simp(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g), x)/sqrt(d + f*x**S(2)))
    rule618 = ReplacementRule(pattern618, replacement618)
    rubi.add(rule618.pattern, label = 618)

    def With619(h, g, e, x, b, d, f, c, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (S(2)*c*g - h*(b - q))*Int(S(1)/((b + S(2)*c*x - q)*sqrt(d + e*x + f*x**S(2))), x)/q - (S(2)*c*g - h*(b + q))*Int(S(1)/((b + S(2)*c*x + q)*sqrt(d + e*x + f*x**S(2))), x)/q
    pattern619 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, cons231 )
    rule619 = ReplacementRule(pattern619, With619)
    rubi.add(rule619.pattern, label = 619)

    def With620(h, g, e, x, d, f, c, a):
        q = Rt(-a*c, S(2))
        return (-c*g/(S(2)*q) + h/S(2))*Int(S(1)/((c*x + q)*sqrt(d + e*x + f*x**S(2))), x) + (c*g/(S(2)*q) + h/S(2))*Int(S(1)/((c*x - q)*sqrt(d + e*x + f*x**S(2))), x)
    pattern620 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons394, cons412 )
    rule620 = ReplacementRule(pattern620, With620)
    rubi.add(rule620.pattern, label = 620)

    def With621(h, d, g, x, b, f, c, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (S(2)*c*g - h*(b - q))*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x - q)), x)/q - (S(2)*c*g - h*(b + q))*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x + q)), x)/q
    pattern621 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons226, cons231 )
    rule621 = ReplacementRule(pattern621, With621)
    rubi.add(rule621.pattern, label = 621)

    def With622(h, g, e, x, b, d, f, c, a):
        q = Rt(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(-g*(-a*f + c*d - q) + h*(-a*e + b*d) - x*(g*(-b*f + c*e) - h*(-a*f + c*d + q)), x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) - Int(Simp(-g*(-a*f + c*d + q) + h*(-a*e + b*d) - x*(g*(-b*f + c*e) - h*(-a*f + c*d - q)), x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    pattern622 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394, cons372, cons413 )
    rule622 = ReplacementRule(pattern622, With622)
    rubi.add(rule622.pattern, label = 622)

    def With623(h, g, e, x, d, f, c, a):
        q = Rt(a*c*e**S(2) + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(-a*e*h - g*(-a*f + c*d - q) + x*(-c*e*g + h*(-a*f + c*d + q)), x)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) - Int(Simp(-a*e*h - g*(-a*f + c*d + q) + x*(-c*e*g + h*(-a*f + c*d - q)), x)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    pattern623 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons394, cons414 )
    rule623 = ReplacementRule(pattern623, With623)
    rubi.add(rule623.pattern, label = 623)

    def With624(h, d, g, x, b, f, c, a):
        q = Rt(b**S(2)*d*f + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(b*d*h - g*(-a*f + c*d - q) + x*(b*f*g + h*(-a*f + c*d + q)), x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q) - Int(Simp(b*d*h - g*(-a*f + c*d + q) + x*(b*f*g + h*(-a*f + c*d - q)), x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q)
    pattern624 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons226, cons413 )
    rule624 = ReplacementRule(pattern624, With624)
    rubi.add(rule624.pattern, label = 624)

    def With625(h, d, g, e, x, b, f, c, a):
        s = Rt(-S(4)*a*c + b**S(2), S(2))
        t = Rt(-S(4)*d*f + e**S(2), S(2))
        return sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + x*(e + t))*sqrt(b + S(2)*c*x + s)*sqrt(e + S(2)*f*x + t)*Int((g + h*x)/(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + x*(e + t))*sqrt(b + S(2)*c*x + s)*sqrt(e + S(2)*f*x + t)), x)/(sqrt(a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2)))
    pattern625 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons394 )
    rule625 = ReplacementRule(pattern625, With625)
    rubi.add(rule625.pattern, label = 625)

    def With626(h, d, g, x, b, f, c, a):
        s = Rt(-S(4)*a*c + b**S(2), S(2))
        t = Rt(-S(4)*d*f, S(2))
        return sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + t*x)*sqrt(S(2)*f*x + t)*sqrt(b + S(2)*c*x + s)*Int((g + h*x)/(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + t*x)*sqrt(S(2)*f*x + t)*sqrt(b + S(2)*c*x + s)), x)/(sqrt(d + f*x**S(2))*sqrt(a + b*x + c*x**S(2)))
    pattern626 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons209, cons226 )
    rule626 = ReplacementRule(pattern626, With626)
    rubi.add(rule626.pattern, label = 626)

    def With627(h, g, e, x, b, d, f, c, a):
        q = S(3)**(S(2)/3)*(-c*h**S(2)/(-b*h + S(2)*c*g)**S(2))**(S(1)/3)
        return sqrt(S(3))*h*q*ArcTan(-S(2)**(S(2)/3)*sqrt(S(3))*(-S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(2)/3)/(S(3)*(S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(1)/3)) + sqrt(S(3))/S(3))/f - S(3)*h*q*log((-S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(2)/3) + S(2)**(S(1)/3)*(S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(1)/3))/(S(2)*f) + h*q*log(d + e*x + f*x**S(2))/(S(2)*f)
    pattern627 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons410, cons443, cons444, cons445 )
    rule627 = ReplacementRule(pattern627, With627)
    rubi.add(rule627.pattern, label = 627)

    pattern628 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons446, cons447, cons43)
    def replacement628(h, d, g, x, f, c, a):
        return S(2)**(S(1)/3)*sqrt(S(3))*h*ArcTan(-S(2)**(S(2)/3)*sqrt(S(3))*(S(1) - S(3)*h*x/g)**(S(2)/3)/(S(3)*(S(1) + S(3)*h*x/g)**(S(1)/3)) + sqrt(S(3))/S(3))/(S(2)*a**(S(1)/3)*f) + S(2)**(S(1)/3)*h*log(d + f*x**S(2))/(S(4)*a**(S(1)/3)*f) - S(3)*S(2)**(S(1)/3)*h*log((S(1) - S(3)*h*x/g)**(S(2)/3) + S(2)**(S(1)/3)*(S(1) + S(3)*h*x/g)**(S(1)/3))/(S(4)*a**(S(1)/3)*f)
    rule628 = ReplacementRule(pattern628, replacement628)
    rubi.add(rule628.pattern, label = 628)

    def With629(h, g, e, x, b, d, f, c, a):
        q = -c/(-S(4)*a*c + b**S(2))
        return (q*(a + b*x + c*x**S(2)))**(S(1)/3)*Int((g + h*x)/((d + e*x + f*x**S(2))*(a*q + b*q*x + c*q*x**S(2))**(S(1)/3)), x)/(a + b*x + c*x**S(2))**(S(1)/3)
    pattern629 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons410, cons443, cons444, cons313 )
    rule629 = ReplacementRule(pattern629, With629)
    rubi.add(rule629.pattern, label = 629)

    pattern630 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons7, cons27, cons125, cons208, cons209, cons446, cons447, cons448)
    def replacement630(h, d, g, x, f, c, a):
        return (S(1) + c*x**S(2)/a)**(S(1)/3)*Int((g + h*x)/((S(1) + c*x**S(2)/a)**(S(1)/3)*(d + f*x**S(2))), x)/(a + c*x**S(2))**(S(1)/3)
    rule630 = ReplacementRule(pattern630, replacement630)
    rubi.add(rule630.pattern, label = 630)

    pattern631 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons449)
    def replacement631(h, p, q, g, e, x, b, d, f, c, a):
        return Int((g + h*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)
    rule631 = ReplacementRule(pattern631, replacement631)
    rubi.add(rule631.pattern, label = 631)

    pattern632 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons450)
    def replacement632(h, p, q, g, e, x, d, f, c, a):
        return Int((a + c*x**S(2))**p*(g + h*x)*(d + e*x + f*x**S(2))**q, x)
    rule632 = ReplacementRule(pattern632, replacement632)
    rubi.add(rule632.pattern, label = 632)

    pattern633 = Pattern(Integral((u_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons50, cons68, cons69)
    def replacement633(m, h, p, q, g, e, x, b, d, f, c, u, a):
        return Subst(Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1))
    rule633 = ReplacementRule(pattern633, replacement633)
    rubi.add(rule633.pattern, label = 633)

    pattern634 = Pattern(Integral((u_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons50, cons68, cons69)
    def replacement634(m, h, p, q, g, e, x, d, f, c, u, a):
        return Subst(Int((a + c*x**S(2))**p*(g + h*x)**m*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1))
    rule634 = ReplacementRule(pattern634, replacement634)
    rubi.add(rule634.pattern, label = 634)

    pattern635 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*z_**WC('m', S(1)), x_), cons21, cons5, cons50, cons451, cons452, cons453)
    def replacement635(m, p, q, v, z, x, u):
        return Int(ExpandToSum(u, x)**p*ExpandToSum(v, x)**q*ExpandToSum(z, x)**m, x)
    rule635 = ReplacementRule(pattern635, replacement635)
    rubi.add(rule635.pattern, label = 635)

    pattern636 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons21, cons4, cons5, cons50, cons336, cons124, cons377)
    def replacement636(m, n, h, q, p, g, e, x, b, d, i, f, c, a):
        return Int((h + i*x)**q*(d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x)
    rule636 = ReplacementRule(pattern636, replacement636)
    rubi.add(rule636.pattern, label = 636)

    pattern637 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons21, cons4, cons5, cons50, cons128, cons84)
    def replacement637(m, n, h, q, p, i, g, e, x, b, d, f, c, a):
        return Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(h + i*x)**q*(a + b*x + c*x**S(2))**p, x), x)
    rule637 = ReplacementRule(pattern637, replacement637)
    rubi.add(rule637.pattern, label = 637)

    pattern638 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons21, cons4, cons5, cons50, cons336, cons124)
    def replacement638(m, n, h, q, p, g, e, x, b, d, i, f, c, a):
        return (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((h + i*x)**q*(d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x)
    rule638 = ReplacementRule(pattern638, replacement638)
    rubi.add(rule638.pattern, label = 638)

    pattern639 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons5, cons50, cons382, cons383, cons384, cons385)
    def replacement639(p, q, d, B, a, e, x, b, f, C, c, A):
        return (c/f)**p*Int((A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x)
    rule639 = ReplacementRule(pattern639, replacement639)
    rubi.add(rule639.pattern, label = 639)

    pattern640 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons36, cons5, cons50, cons382, cons383, cons384, cons385)
    def replacement640(p, q, d, a, e, x, b, f, C, c, A):
        return (c/f)**p*Int((A + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x)
    rule640 = ReplacementRule(pattern640, replacement640)
    rubi.add(rule640.pattern, label = 640)

    pattern641 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons5, cons50, cons382, cons383, cons147, cons386, cons387)
    def replacement641(p, q, d, B, a, e, x, b, f, C, c, A):
        return a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x)
    rule641 = ReplacementRule(pattern641, replacement641)
    rubi.add(rule641.pattern, label = 641)

    pattern642 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons36, cons5, cons50, cons382, cons383, cons147, cons386, cons387)
    def replacement642(p, q, d, a, e, x, b, f, C, c, A):
        return a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((A + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x)
    rule642 = ReplacementRule(pattern642, replacement642)
    rubi.add(rule642.pattern, label = 642)

    pattern643 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons5, cons50, cons45)
    def replacement643(p, q, B, a, e, x, b, d, C, f, c, A):
        return (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x)
    rule643 = ReplacementRule(pattern643, replacement643)
    rubi.add(rule643.pattern, label = 643)

    pattern644 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons36, cons5, cons50, cons45)
    def replacement644(p, q, a, e, x, b, d, C, f, c, A):
        return (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((A + C*x**S(2))*(b + S(2)*c*x)**(S(2)*p)*(d + e*x + f*x**S(2))**q, x)
    rule644 = ReplacementRule(pattern644, replacement644)
    rubi.add(rule644.pattern, label = 644)

    pattern645 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons35, cons36, cons5, cons50, cons45)
    def replacement645(p, q, B, a, x, b, d, C, f, c, A):
        return (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q*(A + B*x + C*x**S(2)), x)
    rule645 = ReplacementRule(pattern645, replacement645)
    rubi.add(rule645.pattern, label = 645)

    pattern646 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons36, cons5, cons50, cons45)
    def replacement646(p, q, a, x, b, d, C, f, c, A):
        return (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((A + C*x**S(2))*(b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q, x)
    rule646 = ReplacementRule(pattern646, replacement646)
    rubi.add(rule646.pattern, label = 646)

    pattern647 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons226, cons394, cons220, cons163)
    def replacement647(p, q, d, B, a, e, x, b, f, C, c, A):
        return Int(ExpandIntegrand((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x)
    rule647 = ReplacementRule(pattern647, replacement647)
    rubi.add(rule647.pattern, label = 647)

    pattern648 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons36, cons226, cons394, cons220, cons163)
    def replacement648(p, q, d, a, e, x, b, f, C, c, A):
        return Int(ExpandIntegrand((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x)
    rule648 = ReplacementRule(pattern648, replacement648)
    rubi.add(rule648.pattern, label = 648)

    pattern649 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons394, cons220, cons432)
    def replacement649(p, q, d, B, a, e, x, f, C, c, A):
        return Int(ExpandIntegrand((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x), x)
    rule649 = ReplacementRule(pattern649, replacement649)
    rubi.add(rule649.pattern, label = 649)

    pattern650 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons36, cons394, cons220, cons432)
    def replacement650(p, q, d, a, e, x, f, C, c, A):
        return Int(ExpandIntegrand((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x)
    rule650 = ReplacementRule(pattern650, replacement650)
    rubi.add(rule650.pattern, label = 650)

    pattern651 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons226, cons394, cons402, cons137, cons403)
    def replacement651(p, q, d, B, a, e, x, b, f, C, c, A):
        return (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(A*b*c - S(2)*B*a*c + C*a*b - x*(-C*(-S(2)*a*c + b**S(2)) + c*(-S(2)*A*c + B*b)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-d*(C*(S(2)*a*c - b**S(2)*(p + S(2))) + c*(S(2)*p + S(3))*(-S(2)*A*c + B*b)) + e*q*(A*b*c - S(2)*B*a*c + C*a*b) - f*x**S(2)*(C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + S(2)*q + S(3))) + x*(-e*(C*(S(2)*a*c*(q + S(1)) - b**S(2)*(p + q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + q + S(3))) + S(2)*f*q*(A*b*c - S(2)*B*a*c + C*a*b)), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule651 = ReplacementRule(pattern651, replacement651)
    rubi.add(rule651.pattern, label = 651)

    pattern652 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons36, cons226, cons394, cons402, cons137, cons403)
    def replacement652(p, q, d, a, e, x, b, f, C, c, A):
        return (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(A*b*c + C*a*b + x*(S(2)*A*c**S(2) + C*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*(b*e*q + S(2)*c*d*(S(2)*p + S(3))) - C*(-a*b*e*q + S(2)*a*c*d - b**S(2)*d*(p + S(2))) - f*x**S(2)*(-S(2)*A*c**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2)))) + x*(S(2)*A*c*(b*f*q + c*e*(S(2)*p + q + S(3))) + C*(S(2)*a*b*f*q - S(2)*a*c*e*(q + S(1)) + b**S(2)*e*(p + q + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule652 = ReplacementRule(pattern652, replacement652)
    rubi.add(rule652.pattern, label = 652)

    pattern653 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons394, cons402, cons137, cons403)
    def replacement653(p, q, d, B, a, e, x, f, C, c, A):
        return (a + c*x**S(2))**(p + S(1))*(B*a - x*(A*c - C*a))*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*d*(S(2)*p + S(3)) - a*(B*e*q + C*d) - f*x**S(2)*(-A*c*(S(2)*p + S(2)*q + S(3)) + C*a*(S(2)*q + S(1))) + x*(A*c*e*(S(2)*p + q + S(3)) - a*(S(2)*B*f*q + C*e*(q + S(1)))), x), x)/(S(2)*a*c*(p + S(1)))
    rule653 = ReplacementRule(pattern653, replacement653)
    rubi.add(rule653.pattern, label = 653)

    pattern654 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons36, cons394, cons402, cons137, cons403)
    def replacement654(p, q, d, a, e, x, f, C, c, A):
        return -x*(a + c*x**S(2))**(p + S(1))*(A*c - C*a)*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*d*(S(2)*p + S(3)) - C*a*d - f*x**S(2)*(-A*c*(S(2)*p + S(2)*q + S(3)) + C*a*(S(2)*q + S(1))) + x*(A*c*e*(S(2)*p + q + S(3)) - C*a*e*(q + S(1))), x), x)/(S(2)*a*c*(p + S(1)))
    rule654 = ReplacementRule(pattern654, replacement654)
    rubi.add(rule654.pattern, label = 654)

    pattern655 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons35, cons36, cons226, cons402, cons137, cons403)
    def replacement655(p, q, d, B, a, x, b, f, C, c, A):
        return (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*c - S(2)*B*a*c + C*a*b - x*(-C*(-S(2)*a*c + b**S(2)) + c*(-S(2)*A*c + B*b)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-d*(C*(S(2)*a*c - b**S(2)*(p + S(2))) + c*(S(2)*p + S(3))*(-S(2)*A*c + B*b)) + S(2)*f*q*x*(A*b*c - S(2)*B*a*c + C*a*b) - f*x**S(2)*(C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + S(2)*q + S(3))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule655 = ReplacementRule(pattern655, replacement655)
    rubi.add(rule655.pattern, label = 655)

    pattern656 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons36, cons226, cons402, cons137, cons403)
    def replacement656(p, q, d, a, x, b, f, C, c, A):
        return (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*c + C*a*b + x*(S(2)*A*c**S(2) + C*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(S(2)*A*c**S(2)*d*(S(2)*p + S(3)) - C*(S(2)*a*c*d - b**S(2)*d*(p + S(2))) - f*x**S(2)*(-S(2)*A*c**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2)))) + x*(S(2)*A*b*c*f*q + S(2)*C*a*b*f*q), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule656 = ReplacementRule(pattern656, replacement656)
    rubi.add(rule656.pattern, label = 656)

    pattern657 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons50, cons226, cons394, cons13, cons137, cons404, cons405)
    def replacement657(p, q, d, B, a, e, x, b, f, C, c, A):
        return (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(c*x*(A*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) - B*(a*b*f - S(2)*a*c*e + b*c*d) + C*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e)) - e*((A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - x*(S(2)*f*((A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e))) + (p + S(1))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))*(-S(2)*A*c + B*b - S(2)*C*a) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)))
    rule657 = ReplacementRule(pattern657, replacement657)
    rubi.add(rule657.pattern, label = 657)

    pattern658 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons36, cons50, cons226, cons394, cons13, cons137, cons404, cons405)
    def replacement658(p, q, d, a, e, x, b, f, C, c, A):
        return (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + c*x*(A*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + C*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e)) - e*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - x*(S(2)*f*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)))
    rule658 = ReplacementRule(pattern658, replacement658)
    rubi.add(rule658.pattern, label = 658)

    pattern659 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons50, cons394, cons13, cons137, cons407, cons405)
    def replacement659(p, q, d, B, a, e, x, f, C, c, A):
        return -(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a) + c*x*(A*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*B*a*c*e - S(2)*C*a*(-a*f + c*d)))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(S(2)*p + S(2)*q + S(5)) - e*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a))*(p + q + S(2)) - x*(c*e*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(S(2)*p + q + S(4)) + S(2)*f*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a))*(p + q + S(2))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + (S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)))
    rule659 = ReplacementRule(pattern659, replacement659)
    rubi.add(rule659.pattern, label = 659)

    pattern660 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons36, cons50, cons394, cons13, cons137, cons407, cons405)
    def replacement660(p, q, d, a, e, x, f, C, c, A):
        return -(a + c*x**S(2))**(p + S(1))*(S(2)*a*c*e*(A*c - C*a) + c*x*(A*(-S(2)*a*c*f + S(2)*c**S(2)*d) - S(2)*C*a*(-a*f + c*d)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-S(2)*a*c*e**S(2)*(A*c - C*a)*(p + q + S(2)) - c*f*x**S(2)*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(S(2)*p + S(2)*q + S(5)) - x*(S(4)*a*c*e*f*(A*c - C*a)*(p + q + S(2)) + c*e*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(S(2)*p + q + S(4))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + (S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)))
    rule660 = ReplacementRule(pattern660, replacement660)
    rubi.add(rule660.pattern, label = 660)

    pattern661 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons35, cons36, cons50, cons226, cons13, cons137, cons406, cons405)
    def replacement661(p, q, d, B, a, x, b, f, C, c, A):
        return (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*(A*c - C*a)*(a*f + c*d) + c*x*(A*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - B*(a*b*f + b*c*d) + C*(-S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*b - B*a)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)) - x*(-b*f*(p + S(1))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)) + S(2)*f*(-b*(A*c - C*a)*(a*f + c*d) + (A*b - B*a)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (p + S(1))*(b**S(2)*d*f + (-a*f + c*d)**S(2))*(-S(2)*A*c + B*b - S(2)*C*a) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2)))
    rule661 = ReplacementRule(pattern661, replacement661)
    rubi.add(rule661.pattern, label = 661)

    pattern662 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons36, cons50, cons226, cons13, cons137, cons406, cons405)
    def replacement662(p, q, d, a, x, b, f, C, c, A):
        return (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - b*(A*c - C*a)*(a*f + c*d) + c*x*(A*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + C*(-S(2)*a*(-a*f + c*d) + b**S(2)*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)) - x*(-b*f*(p + S(1))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)) + S(2)*f*(A*b*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - b*(A*c - C*a)*(a*f + c*d))*(p + q + S(2))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2)))
    rule662 = ReplacementRule(pattern662, replacement662)
    rubi.add(rule662.pattern, label = 662)

    pattern663 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons50, cons226, cons394, cons13, cons163, cons433, cons454)
    def replacement663(p, q, d, B, a, e, x, b, f, C, c, A):
        return (a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))*(B*c*f*(S(2)*p + S(2)*q + S(3)) + S(2)*C*c*f*x*(p + q + S(1)) + C*(b*f*p - c*e*(S(2)*p + q + S(2))))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(p*(-a*e + b*d)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(p*(-b*f + c*e)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + x*(S(2)*p*(-a*f + c*d)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (C*e*f*p*(-S(4)*a*c + b**S(2)) - b*c*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))) + f*(-S(2)*A*f + B*e)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3)))
    rule663 = ReplacementRule(pattern663, replacement663)
    rubi.add(rule663.pattern, label = 663)

    pattern664 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons36, cons50, cons226, cons394, cons13, cons163, cons433, cons454)
    def replacement664(p, q, d, a, e, x, b, f, C, c, A):
        return (S(2)*C*c*f*x*(p + q + S(1)) + C*(b*f*p - c*e*(S(2)*p + q + S(2))))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(p*(-a*e + b*d)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + x**S(2)*(p*(-b*f + c*e)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(S(2)*p*(-a*f + c*d)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + (C*e*f*p*(-S(4)*a*c + b**S(2)) - b*c*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2)))))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3)))
    rule664 = ReplacementRule(pattern664, replacement664)
    rubi.add(rule664.pattern, label = 664)

    pattern665 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons50, cons394, cons13, cons163, cons433, cons454)
    def replacement665(p, q, d, B, a, e, x, f, C, c, A):
        return (a + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*c*e*(S(2)*p + q + S(2)) + S(2)*C*c*f*x*(p + q + S(1)))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*c*(C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))) + f*(-S(2)*A*f + B*e)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1)) - a*e*p*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(c*e*p*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (-S(4)*C*a*c*f**S(2)*p - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + x*(-S(4)*C*a*c*e*f*p*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3)))), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3)))
    rule665 = ReplacementRule(pattern665, replacement665)
    rubi.add(rule665.pattern, label = 665)

    pattern666 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons36, cons50, cons394, cons13, cons163, cons433, cons454)
    def replacement666(p, q, d, a, e, x, f, C, c, A):
        return (a + c*x**S(2))**p*(-C*c*e*(S(2)*p + q + S(2)) + S(2)*C*c*f*x*(p + q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))))*(p + q + S(1)) - a*e*p*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(c*e*p*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3))) + (-S(4)*C*a*c*f**S(2)*p - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-S(4)*C*a*c*e*f*p*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3)))), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3)))
    rule666 = ReplacementRule(pattern666, replacement666)
    rubi.add(rule666.pattern, label = 666)

    pattern667 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons35, cons36, cons50, cons226, cons13, cons163, cons433, cons454)
    def replacement667(p, q, d, B, a, x, b, f, C, c, A):
        return (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**p*(B*c*f*(S(2)*p + S(2)*q + S(3)) + C*b*f*p + S(2)*C*c*f*x*(p + q + S(1)))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(b*d*p*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1))) + x**S(2)*(-b*f*p*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1))) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-b*c*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1)))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + S(2)*C*d*f))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3)))
    rule667 = ReplacementRule(pattern667, replacement667)
    rubi.add(rule667.pattern, label = 667)

    pattern668 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons36, cons50, cons226, cons13, cons163, cons433, cons454)
    def replacement668(p, q, d, a, x, b, f, C, c, A):
        return (d + f*x**S(2))**(q + S(1))*(C*b*f*p + S(2)*C*c*f*x*(p + q + S(1)))*(a + b*x + c*x**S(2))**p/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-C*b**S(2)*d*f*p*(q + S(1)) + x**S(2)*(C*b**S(2)*f**S(2)*p*(q + S(1)) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-S(2)*C*b*f*p*(q + S(1))*(-a*f + c*d) - b*c*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + S(2)*C*d*f))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3)))
    rule668 = ReplacementRule(pattern668, replacement668)
    rubi.add(rule668.pattern, label = 668)

    def With669(d, B, a, e, x, b, f, C, c, A):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-A*a*c*f + A*b**S(2)*f - A*b*c*e + A*c**S(2)*d - B*a*b*f + B*a*c*e + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - A*c*e - B*a*f + B*c*d + C*a*e - C*b*d))/(a + b*x + c*x**S(2)), x)/q + Int((A*a*f**S(2) - A*b*e*f - A*c*d*f + A*c*e**S(2) + B*b*d*f - B*c*d*e - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - A*c*e - B*a*f + B*c*d + C*a*e - C*b*d))/(d + e*x + f*x**S(2)), x)/q
        return False
    pattern669 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons226, cons394, CustomConstraint(With669))
    rule669 = ReplacementRule(pattern669, With669)
    rubi.add(rule669.pattern, label = 669)

    def With670(d, a, e, x, b, f, C, c, A):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-A*a*c*f + A*b**S(2)*f - A*b*c*e + A*c**S(2)*d + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - A*c*e + C*a*e - C*b*d))/(a + b*x + c*x**S(2)), x)/q + Int((A*a*f**S(2) - A*b*e*f - A*c*d*f + A*c*e**S(2) - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - A*c*e + C*a*e - C*b*d))/(d + e*x + f*x**S(2)), x)/q
        return False
    pattern670 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons36, cons226, cons394, CustomConstraint(With670))
    rule670 = ReplacementRule(pattern670, With670)
    rubi.add(rule670.pattern, label = 670)

    def With671(d, B, a, x, b, f, C, c, A):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((A*a*f**S(2) - A*c*d*f + B*b*d*f - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - B*a*f + B*c*d - C*b*d))/(d + f*x**S(2)), x)/q + Int((-A*a*c*f + A*b**S(2)*f + A*c**S(2)*d - B*a*b*f + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - B*a*f + B*c*d - C*b*d))/(a + b*x + c*x**S(2)), x)/q
        return False
    pattern671 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons35, cons36, cons226, CustomConstraint(With671))
    rule671 = ReplacementRule(pattern671, With671)
    rubi.add(rule671.pattern, label = 671)

    def With672(d, a, x, b, f, C, c, A):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((A*a*f**S(2) - A*c*d*f - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - C*b*d))/(d + f*x**S(2)), x)/q + Int((-A*a*c*f + A*b**S(2)*f + A*c**S(2)*d + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - C*b*d))/(a + b*x + c*x**S(2)), x)/q
        return False
    pattern672 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons36, cons226, CustomConstraint(With672))
    rule672 = ReplacementRule(pattern672, With672)
    rubi.add(rule672.pattern, label = 672)

    pattern673 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons226, cons394)
    def replacement673(B, a, e, x, b, d, C, f, c, A):
        return C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c - C*a + x*(B*c - C*b))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c
    rule673 = ReplacementRule(pattern673, replacement673)
    rubi.add(rule673.pattern, label = 673)

    pattern674 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons36, cons226, cons394)
    def replacement674(a, e, x, b, d, C, f, c, A):
        return C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c - C*a - C*b*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c
    rule674 = ReplacementRule(pattern674, replacement674)
    rubi.add(rule674.pattern, label = 674)

    pattern675 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons394)
    def replacement675(B, a, e, x, d, C, f, c, A):
        return C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c + B*c*x - C*a)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c
    rule675 = ReplacementRule(pattern675, replacement675)
    rubi.add(rule675.pattern, label = 675)

    pattern676 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons36, cons394)
    def replacement676(a, e, x, d, C, f, c, A):
        return C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + (A*c - C*a)*Int(S(1)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c
    rule676 = ReplacementRule(pattern676, replacement676)
    rubi.add(rule676.pattern, label = 676)

    pattern677 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons35, cons36, cons226)
    def replacement677(B, a, x, b, d, C, f, c, A):
        return C*Int(S(1)/sqrt(d + f*x**S(2)), x)/c + Int((A*c - C*a + x*(B*c - C*b))/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/c
    rule677 = ReplacementRule(pattern677, replacement677)
    rubi.add(rule677.pattern, label = 677)

    pattern678 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons36, cons226)
    def replacement678(a, x, b, d, C, f, c, A):
        return C*Int(S(1)/sqrt(d + f*x**S(2)), x)/c + Int((A*c - C*a - C*b*x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/c
    rule678 = ReplacementRule(pattern678, replacement678)
    rubi.add(rule678.pattern, label = 678)

    pattern679 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons5, cons50, cons455)
    def replacement679(p, q, B, a, e, x, b, d, C, f, c, A):
        return Int((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)
    rule679 = ReplacementRule(pattern679, replacement679)
    rubi.add(rule679.pattern, label = 679)

    pattern680 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons36, cons5, cons50, cons456)
    def replacement680(p, q, a, e, x, b, d, C, f, c, A):
        return Int((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)
    rule680 = ReplacementRule(pattern680, replacement680)
    rubi.add(rule680.pattern, label = 680)

    pattern681 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons5, cons50, cons457)
    def replacement681(p, q, B, a, e, x, d, C, f, c, A):
        return Int((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x)
    rule681 = ReplacementRule(pattern681, replacement681)
    rubi.add(rule681.pattern, label = 681)

    pattern682 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons7, cons27, cons48, cons125, cons34, cons36, cons5, cons50, cons458)
    def replacement682(p, q, a, e, x, d, C, f, c, A):
        return Int((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x)
    rule682 = ReplacementRule(pattern682, replacement682)
    rubi.add(rule682.pattern, label = 682)

    pattern683 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons5, cons50, cons68, cons69)
    def replacement683(p, q, a, B, e, x, b, d, C, f, c, u, A):
        return Subst(Int((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1))
    rule683 = ReplacementRule(pattern683, replacement683)
    rubi.add(rule683.pattern, label = 683)

    pattern684 = Pattern(Integral((u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons5, cons50, cons68, cons69)
    def replacement684(p, q, B, A, e, x, b, d, f, c, u, a):
        return Subst(Int((A + B*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1))
    rule684 = ReplacementRule(pattern684, replacement684)
    rubi.add(rule684.pattern, label = 684)

    pattern685 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons36, cons5, cons50, cons68, cons69)
    def replacement685(p, q, a, e, x, b, d, C, f, c, u, A):
        return Subst(Int((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1))
    rule685 = ReplacementRule(pattern685, replacement685)
    rubi.add(rule685.pattern, label = 685)

    pattern686 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('C', S(1)) + u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons5, cons50, cons68, cons69)
    def replacement686(p, q, B, a, e, x, d, C, f, c, u, A):
        return Subst(Int((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1))
    rule686 = ReplacementRule(pattern686, replacement686)
    rubi.add(rule686.pattern, label = 686)

    pattern687 = Pattern(Integral((u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons35, cons36, cons5, cons50, cons68, cons69)
    def replacement687(p, q, B, A, e, x, d, f, c, u, a):
        return Subst(Int((A + B*x)*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1))
    rule687 = ReplacementRule(pattern687, replacement687)
    rubi.add(rule687.pattern, label = 687)

    pattern688 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons36, cons5, cons50, cons68, cons69)
    def replacement688(p, q, a, e, x, d, C, f, c, u, A):
        return Subst(Int((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1))
    rule688 = ReplacementRule(pattern688, replacement688)
    rubi.add(rule688.pattern, label = 688)

    return (rubi, (rule189, rule190, rule191, rule192, rule193, rule194, rule195, rule196, rule197, rule198, rule199, rule200, rule201, rule202, rule203, rule204, rule205, rule206, rule207, rule208, rule209, rule210, rule211, rule212, rule213, rule214, rule215, rule216, rule217, rule218, rule219, rule220, rule221, rule222, rule223, rule224, rule225, rule226, rule227, rule228, rule229, rule230, rule231, rule232, rule233, rule234, rule235, rule236, rule237, rule238, rule239, rule240, rule241, rule242, rule243, rule244, rule245, rule246, rule247, rule248, rule249, rule250, rule251, rule252, rule253, rule254, rule255, rule256, rule257, rule258, rule259, rule260, rule261, rule262, rule263, rule264, rule265, rule266, rule267, rule268, rule269, rule270, rule271, rule272, rule273, rule274, rule275, rule276, rule277, rule278, rule279, rule280, rule281, rule282, rule283, rule284, rule285, rule286, rule287, rule288, rule289, rule290, rule291, rule292, rule293, rule294, rule295, rule296, rule297, rule298, rule299, rule300, rule301, rule302, rule303, rule304, rule305, rule306, rule307, rule308, rule309, rule310, rule311, rule312, rule313, rule314, rule315, rule316, rule317, rule318, rule319, rule320, rule321, rule322, rule323, rule324, rule325, rule326, rule327, rule328, rule329, rule330, rule331, rule332, rule333, rule334, rule335, rule336, rule337, rule338, rule339, rule340, rule341, rule342, rule343, rule344, rule345, rule346, rule347, rule348, rule349, rule350, rule351, rule352, rule353, rule354, rule355, rule356, rule357, rule358, rule359, rule360, rule361, rule362, rule363, rule364, rule365, rule366, rule367, rule368, rule369, rule370, rule371, rule372, rule373, rule374, rule375, rule376, rule377, rule378, rule379, rule380, rule381, rule382, rule383, rule384, rule385, rule386, rule387, rule388, rule389, rule390, rule391, rule392, rule393, rule394, rule395, rule396, rule397, rule398, rule399, rule400, rule401, rule402, rule403, rule404, rule405, rule406, rule407, rule408, rule409, rule410, rule411, rule412, rule413, rule414, rule415, rule416, rule417, rule418, rule419, rule420, rule421, rule422, rule423, rule424, rule425, rule426, rule427, rule428, rule429, rule430, rule431, rule432, rule433, rule434, rule435, rule436, rule437, rule438, rule439, rule440, rule441, rule442, rule443, rule444, rule445, rule446, rule447, rule448, rule449, rule450, rule451, rule452, rule453, rule454, rule455, rule456, rule457, rule458, rule459, rule460, rule461, rule462, rule463, rule464, rule465, rule466, rule467, rule468, rule469, rule470, rule471, rule472, rule473, rule474, rule475, rule476, rule477, rule478, rule479, rule480, rule481, rule482, rule483, rule484, rule485, rule486, rule487, rule488, rule489, rule490, rule491, rule492, rule493, rule494, rule495, rule496, rule497, rule498, rule499, rule500, rule501, rule502, rule503, rule504, rule505, rule506, rule507, rule508, rule509, rule510, rule511, rule512, rule513, rule514, rule515, rule516, rule517, rule518, rule519, rule520, rule521, rule522, rule523, rule524, rule525, rule526, rule527, rule528, rule529, rule530, rule531, rule532, rule533, rule534, rule535, rule536, rule537, rule538, rule539, rule540, rule541, rule542, rule543, rule544, rule545, rule546, rule547, rule548, rule549, rule550, rule551, rule552, rule553, rule554, rule555, rule556, rule557, rule558, rule559, rule560, rule561, rule562, rule563, rule564, rule565, rule566, rule567, rule568, rule569, rule570, rule571, rule572, rule573, rule574, rule575, rule576, rule577, rule578, rule579, rule580, rule581, rule582, rule583, rule584, rule585, rule586, rule587, rule588, rule589, rule590, rule591, rule592, rule593, rule594, rule595, rule596, rule597, rule598, rule599, rule600, rule601, rule602, rule603, rule604, rule605, rule606, rule607, rule608, rule609, rule610, rule611, rule612, rule613, rule614, rule615, rule616, rule617, rule618, rule619, rule620, rule621, rule622, rule623, rule624, rule625, rule626, rule627, rule628, rule629, rule630, rule631, rule632, rule633, rule634, rule635, rule636, rule637, rule638, rule639, rule640, rule641, rule642, rule643, rule644, rule645, rule646, rule647, rule648, rule649, rule650, rule651, rule652, rule653, rule654, rule655, rule656, rule657, rule658, rule659, rule660, rule661, rule662, rule663, rule664, rule665, rule666, rule667, rule668, rule669, rule670, rule671, rule672, rule673, rule674, rule675, rule676, rule677, rule678, rule679, rule680, rule681, rule682, rule683, rule684, rule685, rule686, rule687, rule688, ))
