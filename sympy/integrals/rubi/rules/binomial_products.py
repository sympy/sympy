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


def binomial_products():
    from sympy.integrals.rubi.constraints import cons461, cons3, cons4, cons5, cons462, cons2, cons463, cons56, cons464, cons89, cons465, cons40, cons466, cons150, cons13, cons165, cons467, cons468, cons45, cons450, cons69, cons139, cons469, cons470, cons471, cons472, cons473, cons474, cons475, cons476, cons477, cons478, cons479, cons480, cons481, cons482, cons483, cons484, cons485, cons486, cons107, cons487, cons488, cons489, cons490, cons198, cons491, cons130, cons359, cons492, cons493, cons494, cons495, cons70, cons71, cons57, cons496, cons59, cons60, cons61, cons62, cons497, cons498, cons499, cons500, cons149, cons8, cons19, cons501, cons502, cons503, cons21, cons504, cons505, cons68, cons506, cons507, cons508, cons509, cons20, cons246, cons96, cons510, cons511, cons512, cons513, cons514, cons515, cons516, cons517, cons518, cons519, cons520, cons521, cons522, cons523, cons64, cons524, cons525, cons526, cons527, cons528, cons529, cons530, cons531, cons33, cons532, cons533, cons534, cons535, cons536, cons537, cons538, cons369, cons539, cons540, cons541, cons542, cons358, cons543, cons25, cons544, cons545, cons546, cons547, cons548, cons549, cons550, cons551, cons552, cons553, cons554, cons555, cons556, cons73, cons557, cons29, cons222, cons52, cons558, cons87, cons559, cons397, cons405, cons65, cons560, cons561, cons562, cons563, cons564, cons565, cons566, cons567, cons568, cons569, cons570, cons571, cons72, cons572, cons573, cons574, cons575, cons404, cons576, cons577, cons578, cons407, cons579, cons580, cons581, cons582, cons583, cons179, cons584, cons585, cons119, cons586, cons587, cons588, cons589, cons388, cons590, cons591, cons592, cons593, cons50, cons55, cons594, cons595, cons596, cons597, cons598, cons95, cons599, cons600, cons601, cons602, cons603, cons604, cons605, cons606, cons90, cons607, cons608, cons609, cons610, cons611, cons612, cons613, cons614, cons615, cons616, cons617, cons618, cons619, cons620, cons621, cons622, cons623, cons624, cons625, cons626, cons627, cons628, cons629, cons48, cons630, cons127, cons631, cons632, cons633, cons155, cons634, cons635, cons178, cons636, cons637, cons638, cons639, cons640, cons180, cons641, cons642, cons398, cons643, cons54, cons644, cons645, cons646, cons647, cons648, cons649, cons650, cons651, cons652, cons653, cons654, cons655, cons656, cons657, cons658, cons210, cons659, cons660, cons661, cons662, cons663, cons382, cons664, cons665


    pattern692 = Pattern(Integral((x_**n_*WC('b', S(1)))**p_, x_), cons3, cons4, cons5, cons461)
    rule692 = ReplacementRule(pattern692, replacement692)

    pattern693 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons462)
    rule693 = ReplacementRule(pattern693, replacement693)

    pattern694 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons463, cons56)
    rule694 = ReplacementRule(pattern694, replacement694)

    pattern695 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons4, cons464)
    rule695 = ReplacementRule(pattern695, replacement695)

    pattern696 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons89, cons465, cons40)
    rule696 = ReplacementRule(pattern696, replacement696)

    pattern697 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons466)
    rule697 = ReplacementRule(pattern697, replacement697)

    pattern698 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons150, cons13, cons165, cons467)
    rule698 = ReplacementRule(pattern698, replacement698)

    pattern699 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-5)/4), x_), cons2, cons3, cons468, cons45)
    rule699 = ReplacementRule(pattern699, replacement699)

    pattern700 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-5)/4), x_), cons2, cons3, cons468, cons450)
    rule700 = ReplacementRule(pattern700, replacement700)

    pattern701 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-7)/6), x_), cons2, cons3, cons69)
    rule701 = ReplacementRule(pattern701, replacement701)

    pattern702 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons150, cons13, cons139, cons467)
    rule702 = ReplacementRule(pattern702, replacement702)

    pattern703 = Pattern(Integral(S(1)/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule703 = ReplacementRule(pattern703, replacement703)

    pattern704 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons469, cons470)
    rule704 = ReplacementRule(pattern704, With704)

    pattern705 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons469, cons471)
    rule705 = ReplacementRule(pattern705, With705)

    pattern706 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons470, cons472)
    rule706 = ReplacementRule(pattern706, replacement706)

    pattern707 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons470, cons473)
    rule707 = ReplacementRule(pattern707, replacement707)

    pattern708 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons470)
    rule708 = ReplacementRule(pattern708, replacement708)

    pattern709 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons471, cons474)
    rule709 = ReplacementRule(pattern709, replacement709)

    pattern710 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons471, cons475)
    rule710 = ReplacementRule(pattern710, replacement710)

    pattern711 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons471)
    rule711 = ReplacementRule(pattern711, replacement711)

    pattern712 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons476, cons470)
    rule712 = ReplacementRule(pattern712, With712)

    pattern713 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons476, cons471)
    rule713 = ReplacementRule(pattern713, With713)

    pattern714 = Pattern(Integral(S(1)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons477)
    rule714 = ReplacementRule(pattern714, With714)

    pattern715 = Pattern(Integral(S(1)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons478)
    rule715 = ReplacementRule(pattern715, With715)

    pattern716 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons479, cons480)
    rule716 = ReplacementRule(pattern716, With716)

    pattern717 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons479, cons478)
    rule717 = ReplacementRule(pattern717, With717)

    pattern718 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons45, cons481)
    rule718 = ReplacementRule(pattern718, replacement718)

    pattern719 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons45, cons482)
    rule719 = ReplacementRule(pattern719, replacement719)

    pattern720 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons450)
    rule720 = ReplacementRule(pattern720, replacement720)

    pattern721 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons483)
    rule721 = ReplacementRule(pattern721, With721)

    pattern722 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons484)
    rule722 = ReplacementRule(pattern722, With722)

    pattern723 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons468)
    rule723 = ReplacementRule(pattern723, With723)

    pattern724 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons485, cons45)
    rule724 = ReplacementRule(pattern724, replacement724)

    pattern725 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons486, cons107, CustomConstraint(With725))
    rule725 = ReplacementRule(pattern725, replacement725)

    pattern726 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons486, cons107)
    rule726 = ReplacementRule(pattern726, With726)

    pattern727 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons485, cons450)
    rule727 = ReplacementRule(pattern727, replacement727)

    pattern728 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule728 = ReplacementRule(pattern728, With728)

    pattern729 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule729 = ReplacementRule(pattern729, replacement729)

    pattern730 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), cons2, cons3, cons468)
    rule730 = ReplacementRule(pattern730, replacement730)

    pattern731 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), cons2, cons3, cons485, cons45)
    rule731 = ReplacementRule(pattern731, replacement731)

    pattern732 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), cons2, cons3, cons485, cons450)
    rule732 = ReplacementRule(pattern732, replacement732)

    pattern733 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), cons2, cons3, cons45, cons468)
    rule733 = ReplacementRule(pattern733, replacement733)

    pattern734 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), cons2, cons3, cons45, cons485)
    rule734 = ReplacementRule(pattern734, replacement734)

    pattern735 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), cons2, cons3, cons450)
    rule735 = ReplacementRule(pattern735, replacement735)

    pattern736 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/3), x_), cons2, cons3, cons69)
    rule736 = ReplacementRule(pattern736, replacement736)

    pattern737 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-2)/3), x_), cons2, cons3, cons69)
    rule737 = ReplacementRule(pattern737, replacement737)

    pattern738 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**(S(-3)/4), x_), cons2, cons3, cons69)
    rule738 = ReplacementRule(pattern738, replacement738)

    pattern739 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/6), x_), cons2, cons3, cons69)
    rule739 = ReplacementRule(pattern739, replacement739)

    pattern740 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons150, cons13, cons487, cons488, cons489)
    rule740 = ReplacementRule(pattern740, replacement740)

    pattern741 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons150, cons13, cons487, cons488, cons490)
    rule741 = ReplacementRule(pattern741, replacement741)

    pattern742 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons198)
    rule742 = ReplacementRule(pattern742, replacement742)

    pattern743 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons491)
    rule743 = ReplacementRule(pattern743, With743)

    pattern744 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons130)
    rule744 = ReplacementRule(pattern744, replacement744)

    pattern745 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons359, cons492, cons493, cons494)
    rule745 = ReplacementRule(pattern745, replacement745)

    pattern746 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons359, cons492, cons493, cons495)
    rule746 = ReplacementRule(pattern746, replacement746)

    pattern747 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons5, cons70, cons71)
    rule747 = ReplacementRule(pattern747, replacement747)

    pattern748 = Pattern(Integral((x_**n_*WC('b1', S(1)) + WC('a1', S(0)))**WC('p', S(1))*(x_**n_*WC('b2', S(1)) + WC('a2', S(0)))**WC('p', S(1)), x_), cons59, cons60, cons61, cons62, cons4, cons5, cons57, cons496)
    rule748 = ReplacementRule(pattern748, replacement748)

    pattern749 = Pattern(Integral((a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), cons59, cons60, cons61, cons62, cons57, cons497, cons13, cons165, cons498)
    rule749 = ReplacementRule(pattern749, replacement749)

    pattern750 = Pattern(Integral((a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons57, cons497, cons13, cons139, cons498)
    rule750 = ReplacementRule(pattern750, replacement750)

    pattern751 = Pattern(Integral((a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons5, cons57, cons499)
    rule751 = ReplacementRule(pattern751, replacement751)

    pattern752 = Pattern(Integral((a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons5, cons57, cons500)
    rule752 = ReplacementRule(pattern752, With752)

    pattern753 = Pattern(Integral((x_**n_*WC('b1', S(1)) + WC('a1', S(0)))**p_*(x_**n_*WC('b2', S(1)) + WC('a2', S(0)))**p_, x_), cons59, cons60, cons61, cons62, cons4, cons5, cons57, cons149)
    rule753 = ReplacementRule(pattern753, replacement753)

    pattern754 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons5, cons57, cons496)
    rule754 = ReplacementRule(pattern754, replacement754)

    pattern755 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_, x_), cons3, cons8, cons19, cons4, cons5, cons501, cons502)
    rule755 = ReplacementRule(pattern755, replacement755)

    pattern756 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons3, cons8, cons19, cons4, cons5, cons501, cons503)
    rule756 = ReplacementRule(pattern756, replacement756)

    pattern757 = Pattern(Integral((c_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons3, cons8, cons19, cons4, cons5, cons21)
    rule757 = ReplacementRule(pattern757, replacement757)

    pattern758 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons19, cons4, cons40, cons504)
    rule758 = ReplacementRule(pattern758, replacement758)

    pattern759 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons505, cons68)
    rule759 = ReplacementRule(pattern759, replacement759)

    pattern760 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons5, cons57, cons506, cons68)
    rule760 = ReplacementRule(pattern760, replacement760)

    pattern761 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons19, cons4, cons5, cons502)
    rule761 = ReplacementRule(pattern761, replacement761)

    pattern762 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons19, cons4, cons5, cons57, cons507)
    rule762 = ReplacementRule(pattern762, replacement762)

    pattern763 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons502)
    rule763 = ReplacementRule(pattern763, replacement763)

    pattern764 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons5, cons57, cons507)
    rule764 = ReplacementRule(pattern764, replacement764)

    pattern765 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons130)
    rule765 = ReplacementRule(pattern765, replacement765)

    pattern766 = Pattern(Integral(x_**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons19, cons4, cons5, cons508, cons68)
    rule766 = ReplacementRule(pattern766, replacement766)

    pattern767 = Pattern(Integral(x_**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons19, cons4, cons5, cons57, cons509, cons68)
    rule767 = ReplacementRule(pattern767, replacement767)

    pattern768 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons508, cons56)
    rule768 = ReplacementRule(pattern768, replacement768)

    pattern769 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons5, cons57, cons509, cons56)
    rule769 = ReplacementRule(pattern769, replacement769)

    pattern770 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons150, cons20, CustomConstraint(With770))
    rule770 = ReplacementRule(pattern770, replacement770)

    pattern771 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons5, cons57, cons497, cons20, CustomConstraint(With771))
    rule771 = ReplacementRule(pattern771, replacement771)

    pattern772 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons150, cons246, cons165, cons96, cons510, cons511)
    rule772 = ReplacementRule(pattern772, replacement772)

    pattern773 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons57, cons497, cons246, cons165, cons512, cons513)
    rule773 = ReplacementRule(pattern773, replacement773)

    pattern774 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons150, cons246, cons165, cons514, cons511)
    rule774 = ReplacementRule(pattern774, replacement774)

    pattern775 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons468)
    rule775 = ReplacementRule(pattern775, replacement775)

    pattern776 = Pattern(Integral(x_**m_/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons468, cons515)
    rule776 = ReplacementRule(pattern776, replacement776)

    pattern777 = Pattern(Integral(x_**m_/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons468, cons516)
    rule777 = ReplacementRule(pattern777, replacement777)

    pattern778 = Pattern(Integral(sqrt(x_*WC('c', S(1)))/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons8, cons468)
    rule778 = ReplacementRule(pattern778, replacement778)

    pattern779 = Pattern(Integral((x_*WC('c', S(1)))**m_/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons8, cons468, cons517, cons518)
    rule779 = ReplacementRule(pattern779, replacement779)

    pattern780 = Pattern(Integral((x_*WC('c', S(1)))**m_/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons8, cons468, cons517, cons96)
    rule780 = ReplacementRule(pattern780, replacement780)

    pattern781 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons485)
    rule781 = ReplacementRule(pattern781, replacement781)

    pattern782 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons150, cons246, cons139, cons519, cons520, cons511)
    rule782 = ReplacementRule(pattern782, replacement782)

    pattern783 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons57, cons497, cons246, cons139, cons521, cons522, cons513)
    rule783 = ReplacementRule(pattern783, replacement783)

    pattern784 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons150, cons246, cons139, cons511)
    rule784 = ReplacementRule(pattern784, replacement784)

    pattern785 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons57, cons497, cons246, cons139, cons513)
    rule785 = ReplacementRule(pattern785, replacement785)

    pattern786 = Pattern(Integral(x_/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule786 = ReplacementRule(pattern786, replacement786)

    pattern787 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons523, cons64, cons524, cons470)
    rule787 = ReplacementRule(pattern787, With787)

    pattern788 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons525, cons64, cons524, cons471)
    rule788 = ReplacementRule(pattern788, With788)

    pattern789 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons526, cons64, cons524, cons470)
    rule789 = ReplacementRule(pattern789, With789)

    pattern790 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons526, cons64, cons524, cons471)
    rule790 = ReplacementRule(pattern790, With790)

    pattern791 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons477)
    rule791 = ReplacementRule(pattern791, With791)

    pattern792 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons478)
    rule792 = ReplacementRule(pattern792, With792)

    pattern793 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons527, cons64, cons524, cons480)
    rule793 = ReplacementRule(pattern793, With793)

    pattern794 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons527, cons64, cons528, cons478)
    rule794 = ReplacementRule(pattern794, With794)

    pattern795 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons527, cons64, cons529, cons478)
    rule795 = ReplacementRule(pattern795, With795)

    pattern796 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons530, cons531)
    rule796 = ReplacementRule(pattern796, replacement796)

    pattern797 = Pattern(Integral(x_/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons483)
    rule797 = ReplacementRule(pattern797, With797)

    pattern798 = Pattern(Integral(x_/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons484)
    rule798 = ReplacementRule(pattern798, With798)

    pattern799 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons468)
    rule799 = ReplacementRule(pattern799, With799)

    pattern800 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons486, cons107)
    rule800 = ReplacementRule(pattern800, With800)

    pattern801 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons485)
    rule801 = ReplacementRule(pattern801, With801)

    pattern802 = Pattern(Integral(x_**S(4)/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule802 = ReplacementRule(pattern802, With802)

    pattern803 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule803 = ReplacementRule(pattern803, replacement803)

    pattern804 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4), x_), cons2, cons3, cons468)
    rule804 = ReplacementRule(pattern804, replacement804)

    pattern805 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4), x_), cons2, cons3, cons485)
    rule805 = ReplacementRule(pattern805, replacement805)

    pattern806 = Pattern(Integral(S(1)/(x_**S(2)*(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)), x_), cons2, cons3, cons468)
    rule806 = ReplacementRule(pattern806, replacement806)

    pattern807 = Pattern(Integral(S(1)/(x_**S(2)*(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)), x_), cons2, cons3, cons485)
    rule807 = ReplacementRule(pattern807, replacement807)

    pattern808 = Pattern(Integral(sqrt(c_*x_)/(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4), x_), cons2, cons3, cons8, cons468)
    rule808 = ReplacementRule(pattern808, replacement808)

    pattern809 = Pattern(Integral(sqrt(c_*x_)/(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4), x_), cons2, cons3, cons8, cons485)
    rule809 = ReplacementRule(pattern809, replacement809)

    pattern810 = Pattern(Integral(S(1)/((x_*WC('c', S(1)))**(S(3)/2)*(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)), x_), cons2, cons3, cons8, cons468)
    rule810 = ReplacementRule(pattern810, replacement810)

    pattern811 = Pattern(Integral(S(1)/((x_*WC('c', S(1)))**(S(3)/2)*(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)), x_), cons2, cons3, cons8, cons485)
    rule811 = ReplacementRule(pattern811, replacement811)

    pattern812 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons150, cons33, cons532, cons514, cons511)
    rule812 = ReplacementRule(pattern812, replacement812)

    pattern813 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons150, cons533, cons514, cons534)
    rule813 = ReplacementRule(pattern813, replacement813)

    pattern814 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons5, cons57, cons497, cons33, cons531, cons512, cons513)
    rule814 = ReplacementRule(pattern814, replacement814)

    pattern815 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons5, cons57, cons497, cons535, cons512, cons536)
    rule815 = ReplacementRule(pattern815, replacement815)

    pattern816 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons150, cons33, cons96, cons511)
    rule816 = ReplacementRule(pattern816, replacement816)

    pattern817 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons150, cons537, cons534)
    rule817 = ReplacementRule(pattern817, replacement817)

    pattern818 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons5, cons57, cons497, cons33, cons96, cons513)
    rule818 = ReplacementRule(pattern818, replacement818)

    pattern819 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons5, cons57, cons497, cons538, cons536)
    rule819 = ReplacementRule(pattern819, replacement819)

    pattern820 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons150, cons369, cons511)
    rule820 = ReplacementRule(pattern820, With820)

    pattern821 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons5, cons57, cons497, cons369, cons513)
    rule821 = ReplacementRule(pattern821, With821)

    pattern822 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons150, cons13, cons487, cons488, cons539)
    rule822 = ReplacementRule(pattern822, replacement822)

    pattern823 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons57, cons497, cons13, cons487, cons488, cons540)
    rule823 = ReplacementRule(pattern823, replacement823)

    pattern824 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons150, cons13, cons487, cons488, cons20, cons541)
    rule824 = ReplacementRule(pattern824, replacement824)

    pattern825 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons57, cons497, cons13, cons487, cons488, cons20, cons542)
    rule825 = ReplacementRule(pattern825, replacement825)

    pattern826 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons198, cons20)
    rule826 = ReplacementRule(pattern826, replacement826)

    pattern827 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons5, cons57, cons499, cons20)
    rule827 = ReplacementRule(pattern827, replacement827)

    pattern828 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons198, cons369)
    rule828 = ReplacementRule(pattern828, With828)

    pattern829 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons5, cons57, cons499, cons369)
    rule829 = ReplacementRule(pattern829, With829)

    pattern830 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons198, cons358)
    rule830 = ReplacementRule(pattern830, replacement830)

    pattern831 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons5, cons57, cons499, cons358)
    rule831 = ReplacementRule(pattern831, replacement831)

    pattern832 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons19, cons5, cons491)
    rule832 = ReplacementRule(pattern832, With832)

    pattern833 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons19, cons5, cons57, cons500)
    rule833 = ReplacementRule(pattern833, With833)

    pattern834 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons491)
    rule834 = ReplacementRule(pattern834, replacement834)

    pattern835 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons5, cons57, cons500)
    rule835 = ReplacementRule(pattern835, replacement835)

    pattern836 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons19, cons4, cons5, cons543, cons25)
    rule836 = ReplacementRule(pattern836, replacement836)

    pattern837 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons19, cons4, cons5, cons57, cons544, cons545)
    rule837 = ReplacementRule(pattern837, replacement837)

    pattern838 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons543, cons25)
    rule838 = ReplacementRule(pattern838, replacement838)

    pattern839 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons5, cons57, cons544, cons545)
    rule839 = ReplacementRule(pattern839, replacement839)

    pattern840 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons19, cons4, cons546, cons13, cons165)
    rule840 = ReplacementRule(pattern840, replacement840)

    pattern841 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons19, cons4, cons57, cons547, cons13, cons165)
    rule841 = ReplacementRule(pattern841, replacement841)

    pattern842 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons546, cons13, cons165)
    rule842 = ReplacementRule(pattern842, replacement842)

    pattern843 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons57, cons547, cons13, cons165)
    rule843 = ReplacementRule(pattern843, replacement843)

    pattern844 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons548, cons13, cons165, cons514)
    rule844 = ReplacementRule(pattern844, replacement844)

    pattern845 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons57, cons549, cons13, cons165, cons512)
    rule845 = ReplacementRule(pattern845, replacement845)

    pattern846 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons19, cons4, cons548, cons13, cons487)
    rule846 = ReplacementRule(pattern846, With846)

    pattern847 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons19, cons4, cons57, cons549, cons13, cons487)
    rule847 = ReplacementRule(pattern847, With847)

    pattern848 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons548, cons13, cons487)
    rule848 = ReplacementRule(pattern848, replacement848)

    pattern849 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons57, cons549, cons13, cons487)
    rule849 = ReplacementRule(pattern849, replacement849)

    pattern850 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons548, cons13, cons139)
    rule850 = ReplacementRule(pattern850, replacement850)

    pattern851 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons57, cons548, cons13, cons139)
    rule851 = ReplacementRule(pattern851, replacement851)

    pattern852 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons19, cons4, cons550, cons533)
    rule852 = ReplacementRule(pattern852, With852)

    pattern853 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons19, cons4, cons550, cons537)
    rule853 = ReplacementRule(pattern853, replacement853)

    pattern854 = Pattern(Integral((c_*x_)**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons19, cons4, cons550, cons551)
    rule854 = ReplacementRule(pattern854, replacement854)

    pattern855 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons359, cons552)
    rule855 = ReplacementRule(pattern855, replacement855)

    pattern856 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons359, cons553)
    rule856 = ReplacementRule(pattern856, replacement856)

    pattern857 = Pattern(Integral(x_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons5, cons554, cons20, cons555)
    rule857 = ReplacementRule(pattern857, replacement857)

    pattern858 = Pattern(Integral(u_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons19, cons4, cons5, cons556)
    rule858 = ReplacementRule(pattern858, replacement858)

    pattern859 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons5, cons57, cons149)
    rule859 = ReplacementRule(pattern859, replacement859)

    pattern860 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons73, cons557)
    rule860 = ReplacementRule(pattern860, replacement860)

    pattern861 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons73, cons222, cons504)
    rule861 = ReplacementRule(pattern861, replacement861)

    pattern862 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons52, cons73, cons198)
    rule862 = ReplacementRule(pattern862, replacement862)

    pattern863 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons52, cons73, cons491)
    rule863 = ReplacementRule(pattern863, With863)

    pattern864 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons73, cons558, cons87)
    rule864 = ReplacementRule(pattern864, replacement864)

    pattern865 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons73, cons559, cons397, cons405, cons56)
    rule865 = ReplacementRule(pattern865, replacement865)

    pattern866 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons52, cons73, cons559, cons65)
    rule866 = ReplacementRule(pattern866, replacement866)

    pattern867 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons5, cons52, cons73, cons559)
    rule867 = ReplacementRule(pattern867, replacement867)

    pattern868 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons5, cons52, cons73, cons560, cons561)
    rule868 = ReplacementRule(pattern868, replacement868)

    pattern869 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons52, cons73, cons560, cons562, cons56)
    rule869 = ReplacementRule(pattern869, replacement869)

    pattern870 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons73, cons563)
    rule870 = ReplacementRule(pattern870, replacement870)

    pattern871 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons73, cons564)
    rule871 = ReplacementRule(pattern871, replacement871)

    pattern872 = Pattern(Integral((c_ + x_**n_*WC('d', S(1)))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons4, cons73, cons89, cons465)
    rule872 = ReplacementRule(pattern872, replacement872)

    pattern873 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons4, cons73, cons565)
    rule873 = ReplacementRule(pattern873, replacement873)

    pattern874 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons73, cons466, cons566, cons567)
    rule874 = ReplacementRule(pattern874, replacement874)

    pattern875 = Pattern(Integral(S(1)/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons4, cons73)
    rule875 = ReplacementRule(pattern875, replacement875)

    pattern876 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/3)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73, cons568, cons468)
    rule876 = ReplacementRule(pattern876, replacement876)

    pattern877 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/3)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73, cons568, cons485)
    rule877 = ReplacementRule(pattern877, replacement877)

    pattern878 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(2)/3)/(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons73, cons568)
    rule878 = ReplacementRule(pattern878, replacement878)

    pattern879 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73)
    rule879 = ReplacementRule(pattern879, replacement879)

    pattern880 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73)
    rule880 = ReplacementRule(pattern880, replacement880)

    pattern881 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**WC('p', S(1))/(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons73, cons13, cons165, cons569)
    rule881 = ReplacementRule(pattern881, replacement881)

    pattern882 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**p_/(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons73, cons13, cons139, cons570, cons571)
    rule882 = ReplacementRule(pattern882, replacement882)

    pattern883 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons72, cons572)
    rule883 = ReplacementRule(pattern883, replacement883)

    pattern884 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons72, cons573)
    rule884 = ReplacementRule(pattern884, With884)

    pattern885 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons73)
    rule885 = ReplacementRule(pattern885, replacement885)

    pattern886 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)/(c_ + x_**S(4)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons73)
    rule886 = ReplacementRule(pattern886, replacement886)

    pattern887 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**p_/(c_ + x_**S(4)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons73, cons13, cons574)
    rule887 = ReplacementRule(pattern887, replacement887)

    pattern888 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('b', S(1)))*(c_ + x_**S(4)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73)
    rule888 = ReplacementRule(pattern888, replacement888)

    pattern889 = Pattern(Integral(S(1)/((a_ + x_**S(4)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(4)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73)
    rule889 = ReplacementRule(pattern889, replacement889)

    pattern890 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons468, cons575)
    rule890 = ReplacementRule(pattern890, replacement890)

    pattern891 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons73, cons404, cons139, cons576, cons577)
    rule891 = ReplacementRule(pattern891, replacement891)

    pattern892 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons73, cons404, cons139, cons578, cons577)
    rule892 = ReplacementRule(pattern892, replacement892)

    pattern893 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons52, cons73, cons13, cons139, cons407, cons577)
    rule893 = ReplacementRule(pattern893, replacement893)

    pattern894 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons73, cons150, cons222, cons579)
    rule894 = ReplacementRule(pattern894, replacement894)

    pattern895 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons5, cons73, cons397, cons578, cons580, cons581, cons577)
    rule895 = ReplacementRule(pattern895, replacement895)

    pattern896 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons73, cons404, cons405, cons165, cons577)
    rule896 = ReplacementRule(pattern896, replacement896)

    pattern897 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons575, cons468, cons582)
    rule897 = ReplacementRule(pattern897, replacement897)

    pattern898 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons583, cons179, cons45, cons584)
    rule898 = ReplacementRule(pattern898, replacement898)

    pattern899 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons583, cons179, cons585)
    rule899 = ReplacementRule(pattern899, replacement899)

    pattern900 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons119)
    rule900 = ReplacementRule(pattern900, replacement900)

    pattern901 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons575, cons468)
    rule901 = ReplacementRule(pattern901, replacement901)

    pattern902 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons575, cons485)
    rule902 = ReplacementRule(pattern902, replacement902)

    pattern903 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons583, cons179, cons45)
    rule903 = ReplacementRule(pattern903, replacement903)

    pattern904 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons583, cons179, cons585)
    rule904 = ReplacementRule(pattern904, replacement904)

    pattern905 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons583, cons179, cons450)
    rule905 = ReplacementRule(pattern905, replacement905)

    pattern906 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons583, cons119)
    rule906 = ReplacementRule(pattern906, replacement906)

    pattern907 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons52, cons73, cons130)
    rule907 = ReplacementRule(pattern907, replacement907)

    pattern908 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons5, cons52, cons73, cons586, cons45, cons179)
    rule908 = ReplacementRule(pattern908, replacement908)

    pattern909 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons5, cons52, cons73, cons586, cons450)
    rule909 = ReplacementRule(pattern909, replacement909)

    pattern910 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons587, cons588, cons589)
    rule910 = ReplacementRule(pattern910, replacement910)

    pattern911 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons4, cons5, cons52, cons587, cons388, cons149)
    rule911 = ReplacementRule(pattern911, replacement911)

    pattern912 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons52, cons70, cons71)
    rule912 = ReplacementRule(pattern912, replacement912)

    pattern913 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1)), x_), cons5, cons52, cons590)
    rule913 = ReplacementRule(pattern913, replacement913)

    pattern914 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*x_**WC('m', S(1)), x_), cons5, cons52, cons591, cons592)
    rule914 = ReplacementRule(pattern914, replacement914)

    pattern915 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons593, cons502)
    rule915 = ReplacementRule(pattern915, replacement915)

    pattern916 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons593, cons503)
    rule916 = ReplacementRule(pattern916, replacement916)

    pattern917 = Pattern(Integral((e_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons21)
    rule917 = ReplacementRule(pattern917, replacement917)

    pattern918 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons52, cons73, cons55)
    rule918 = ReplacementRule(pattern918, replacement918)

    pattern919 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons73, cons222, cons504)
    rule919 = ReplacementRule(pattern919, replacement919)

    pattern920 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons52, cons73, cons502)
    rule920 = ReplacementRule(pattern920, replacement920)

    pattern921 = Pattern(Integral((e_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons73, cons502)
    rule921 = ReplacementRule(pattern921, replacement921)

    pattern922 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons73, cons557)
    rule922 = ReplacementRule(pattern922, replacement922)

    pattern923 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons73, cons594, cons68)
    rule923 = ReplacementRule(pattern923, replacement923)

    pattern924 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons59, cons60, cons61, cons62, cons8, cons29, cons50, cons19, cons4, cons5, cons595, cons57, cons596, cons68)
    rule924 = ReplacementRule(pattern924, replacement924)

    pattern925 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons73, cons597, cons598, cons95, cons599)
    rule925 = ReplacementRule(pattern925, replacement925)

    pattern926 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons73, cons597, cons68)
    rule926 = ReplacementRule(pattern926, replacement926)

    pattern927 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons73, cons598, cons95, cons599, cons600)
    rule927 = ReplacementRule(pattern927, replacement927)

    pattern928 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons59, cons60, cons61, cons62, cons8, cons29, cons50, cons5, cons595, cons57, cons598, cons95, cons599, cons600)
    rule928 = ReplacementRule(pattern928, replacement928)

    pattern929 = Pattern(Integral(x_**m_*(a_ + x_**S(2)*WC('b', S(1)))**p_*(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons73, cons13, cons139, cons601, cons602)
    rule929 = ReplacementRule(pattern929, replacement929)

    pattern930 = Pattern(Integral(x_**m_*(a_ + x_**S(2)*WC('b', S(1)))**p_*(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons73, cons13, cons139, cons603, cons602)
    rule930 = ReplacementRule(pattern930, replacement930)

    pattern931 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons73, cons13, cons139, cons604)
    rule931 = ReplacementRule(pattern931, replacement931)

    pattern932 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons59, cons60, cons61, cons62, cons8, cons29, cons50, cons19, cons4, cons595, cons57, cons13, cons139, cons604)
    rule932 = ReplacementRule(pattern932, replacement932)

    pattern933 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons73, cons605)
    rule933 = ReplacementRule(pattern933, replacement933)

    pattern934 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons59, cons60, cons61, cons62, cons8, cons29, cons50, cons19, cons4, cons5, cons595, cons57, cons605)
    rule934 = ReplacementRule(pattern934, replacement934)

    pattern935 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons73, cons150, cons130, cons606)
    rule935 = ReplacementRule(pattern935, replacement935)

    pattern936 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons73, cons150, cons95, cons96, cons90)
    rule936 = ReplacementRule(pattern936, replacement936)

    pattern937 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons73, cons150, cons13, cons139)
    rule937 = ReplacementRule(pattern937, replacement937)

    pattern938 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons73, cons150, cons607)
    rule938 = ReplacementRule(pattern938, replacement938)

    pattern939 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons5, cons52, cons73, cons150, cons20, CustomConstraint(With939))
    rule939 = ReplacementRule(pattern939, replacement939)

    pattern940 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons52, cons73, cons150, cons369, cons40)
    rule940 = ReplacementRule(pattern940, With940)

    pattern941 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons73, cons150, cons608, cons139, cons405, cons609, cons610)
    rule941 = ReplacementRule(pattern941, replacement941)

    pattern942 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons73, cons150, cons404, cons139, cons578, cons610)
    rule942 = ReplacementRule(pattern942, replacement942)

    pattern943 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons73, cons150, cons404, cons139, cons576, cons610)
    rule943 = ReplacementRule(pattern943, replacement943)

    pattern944 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons52, cons73, cons150, cons246, cons139, cons611, cons610)
    rule944 = ReplacementRule(pattern944, replacement944)

    pattern945 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons52, cons73, cons150, cons246, cons139, cons612, cons610)
    rule945 = ReplacementRule(pattern945, replacement945)

    pattern946 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons52, cons73, cons150, cons13, cons139, cons610)
    rule946 = ReplacementRule(pattern946, replacement946)

    pattern947 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons73, cons150, cons608, cons405, cons96, cons165, cons610)
    rule947 = ReplacementRule(pattern947, replacement947)

    pattern948 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons73, cons150, cons613, cons578, cons96, cons610)
    rule948 = ReplacementRule(pattern948, replacement948)

    pattern949 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons73, cons150, cons613, cons576, cons96, cons610)
    rule949 = ReplacementRule(pattern949, replacement949)

    pattern950 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons73, cons150, cons404, cons405, cons165, cons610)
    rule950 = ReplacementRule(pattern950, replacement950)

    pattern951 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons73, cons150, cons397, cons578, cons610)
    rule951 = ReplacementRule(pattern951, replacement951)

    pattern952 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons73, cons150, cons613, cons405, cons609, cons610)
    rule952 = ReplacementRule(pattern952, replacement952)

    pattern953 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons52, cons73, cons150, cons33, cons611, cons610)
    rule953 = ReplacementRule(pattern953, replacement953)

    pattern954 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons52, cons73, cons150, cons33, cons96, cons610)
    rule954 = ReplacementRule(pattern954, replacement954)

    pattern955 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons73, cons150, cons33, cons614)
    rule955 = ReplacementRule(pattern955, replacement955)

    pattern956 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons73, cons150)
    rule956 = ReplacementRule(pattern956, replacement956)

    pattern957 = Pattern(Integral(x_**m_/((a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73, cons615, cons616, cons617)
    rule957 = ReplacementRule(pattern957, replacement957)

    pattern958 = Pattern(Integral(x_**S(2)/((a_ + x_**S(4)*WC('b', S(1)))*sqrt(c_ + x_**S(4)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73)
    rule958 = ReplacementRule(pattern958, With958)

    pattern959 = Pattern(Integral(x_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73, cons618)
    rule959 = ReplacementRule(pattern959, With959)

    pattern960 = Pattern(Integral(x_**m_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73, cons618, cons619)
    rule960 = ReplacementRule(pattern960, replacement960)

    pattern961 = Pattern(Integral(x_**m_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73, cons618, cons620)
    rule961 = ReplacementRule(pattern961, replacement961)

    pattern962 = Pattern(Integral(x_**S(2)*sqrt(c_ + x_**S(4)*WC('d', S(1)))/(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons73)
    rule962 = ReplacementRule(pattern962, replacement962)

    pattern963 = Pattern(Integral(x_**WC('m', S(1))*sqrt(c_ + x_**S(3)*WC('d', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons73, cons618, cons621)
    rule963 = ReplacementRule(pattern963, replacement963)

    pattern964 = Pattern(Integral(x_**S(2)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73, cons468, cons575, cons582)
    rule964 = ReplacementRule(pattern964, replacement964)

    pattern965 = Pattern(Integral(x_**n_/(sqrt(a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons73, cons622, cons623)
    rule965 = ReplacementRule(pattern965, replacement965)

    pattern966 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons150, cons246, cons624, cons487)
    rule966 = ReplacementRule(pattern966, With966)

    pattern967 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons5, cons52, cons73, cons198, cons20)
    rule967 = ReplacementRule(pattern967, replacement967)

    pattern968 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons52, cons198, cons369)
    rule968 = ReplacementRule(pattern968, With968)

    pattern969 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons52, cons73, cons198, cons358)
    rule969 = ReplacementRule(pattern969, replacement969)

    pattern970 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons19, cons5, cons52, cons73, cons491)
    rule970 = ReplacementRule(pattern970, With970)

    pattern971 = Pattern(Integral((e_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons52, cons73, cons491)
    rule971 = ReplacementRule(pattern971, replacement971)

    pattern972 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons52, cons73, cons543, cons25)
    rule972 = ReplacementRule(pattern972, replacement972)

    pattern973 = Pattern(Integral((e_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons73, cons543, cons25)
    rule973 = ReplacementRule(pattern973, replacement973)

    pattern974 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons73, cons404, cons139, cons578, cons610)
    rule974 = ReplacementRule(pattern974, replacement974)

    pattern975 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons73, cons404, cons139, cons576, cons610)
    rule975 = ReplacementRule(pattern975, replacement975)

    pattern976 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons52, cons73, cons13, cons139, cons610)
    rule976 = ReplacementRule(pattern976, replacement976)

    pattern977 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons73, cons404, cons405, cons165, cons610)
    rule977 = ReplacementRule(pattern977, replacement977)

    pattern978 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons73, cons397, cons578, cons610)
    rule978 = ReplacementRule(pattern978, replacement978)

    pattern979 = Pattern(Integral(x_**m_/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons73, cons625)
    rule979 = ReplacementRule(pattern979, replacement979)

    pattern980 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons19, cons73)
    rule980 = ReplacementRule(pattern980, replacement980)

    pattern981 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons73, cons626, cons627, cons628)
    rule981 = ReplacementRule(pattern981, replacement981)

    pattern982 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons587, cons588, cons589)
    rule982 = ReplacementRule(pattern982, replacement982)

    pattern983 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons52, cons587, cons388, cons149)
    rule983 = ReplacementRule(pattern983, replacement983)

    pattern984 = Pattern(Integral((e_*x_)**m_*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons587)
    rule984 = ReplacementRule(pattern984, replacement984)

    pattern985 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons73, cons68, cons629, cons45, cons179)
    rule985 = ReplacementRule(pattern985, replacement985)

    pattern986 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons73, cons68, cons629, cons450)
    rule986 = ReplacementRule(pattern986, replacement986)

    pattern987 = Pattern(Integral(x_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons52, cons554, cons20, cons555)
    rule987 = ReplacementRule(pattern987, replacement987)

    pattern988 = Pattern(Integral(u_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons52, cons556)
    rule988 = ReplacementRule(pattern988, replacement988)

    pattern989 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons59, cons60, cons61, cons62, cons8, cons29, cons4, cons5, cons52, cons595, cons57, cons496)
    rule989 = ReplacementRule(pattern989, replacement989)

    pattern990 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)) + x_**WC('n2', S(1))*WC('e', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons59, cons60, cons61, cons62, cons8, cons29, cons50, cons4, cons5, cons52, cons595, cons48, cons57, cons496)
    rule990 = ReplacementRule(pattern990, replacement990)

    pattern991 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**p_*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons59, cons60, cons61, cons62, cons8, cons29, cons4, cons5, cons52, cons595, cons57)
    rule991 = ReplacementRule(pattern991, replacement991)

    pattern992 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)) + x_**WC('n2', S(1))*WC('e', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons59, cons60, cons61, cons62, cons8, cons29, cons50, cons4, cons5, cons52, cons595, cons48, cons57)
    rule992 = ReplacementRule(pattern992, replacement992)

    pattern993 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons630)
    rule993 = ReplacementRule(pattern993, replacement993)

    pattern994 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons631)
    rule994 = ReplacementRule(pattern994, replacement994)

    pattern995 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/((a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons631)
    rule995 = ReplacementRule(pattern995, replacement995)

    pattern996 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/(sqrt(a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons632)
    rule996 = ReplacementRule(pattern996, replacement996)

    pattern997 = Pattern(Integral((e_ + x_**S(2)*WC('f', S(1)))/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons468, cons575)
    rule997 = ReplacementRule(pattern997, replacement997)

    pattern998 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons404, cons139, cons405)
    rule998 = ReplacementRule(pattern998, replacement998)

    pattern999 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons52, cons13, cons139)
    rule999 = ReplacementRule(pattern999, replacement999)

    pattern1000 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons397, cons405, cons633)
    rule1000 = ReplacementRule(pattern1000, replacement1000)

    pattern1001 = Pattern(Integral((e_ + x_**S(4)*WC('f', S(1)))/((a_ + x_**S(4)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(4)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons155)
    rule1001 = ReplacementRule(pattern1001, replacement1001)

    pattern1002 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons4, cons634)
    rule1002 = ReplacementRule(pattern1002, replacement1002)

    pattern1003 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons52, cons635)
    rule1003 = ReplacementRule(pattern1003, replacement1003)

    pattern1004 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons155)
    rule1004 = ReplacementRule(pattern1004, replacement1004)

    pattern1005 = Pattern(Integral(S(1)/(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons8, cons29, cons50, cons127, cons178)
    rule1005 = ReplacementRule(pattern1005, replacement1005)

    pattern1006 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons636, cons637, cons638)
    rule1006 = ReplacementRule(pattern1006, replacement1006)

    pattern1007 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons639)
    rule1007 = ReplacementRule(pattern1007, replacement1007)

    pattern1008 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons575, cons640, cons638)
    rule1008 = ReplacementRule(pattern1008, replacement1008)

    pattern1009 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons583, cons179, cons180, cons641)
    rule1009 = ReplacementRule(pattern1009, replacement1009)

    pattern1010 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons119)
    rule1010 = ReplacementRule(pattern1010, replacement1010)

    pattern1011 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons575)
    rule1011 = ReplacementRule(pattern1011, replacement1011)

    pattern1012 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons583)
    rule1012 = ReplacementRule(pattern1012, replacement1012)

    pattern1013 = Pattern(Integral(sqrt(e_ + x_**S(2)*WC('f', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons575, cons640)
    rule1013 = ReplacementRule(pattern1013, replacement1013)

    pattern1014 = Pattern(Integral((e_ + x_**S(2)*WC('f', S(1)))**(S(3)/2)/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons575, cons640)
    rule1014 = ReplacementRule(pattern1014, replacement1014)

    pattern1015 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons575, cons640)
    rule1015 = ReplacementRule(pattern1015, replacement1015)

    pattern1016 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons642, cons398, cons643)
    rule1016 = ReplacementRule(pattern1016, replacement1016)

    pattern1017 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons54, cons397, cons578)
    rule1017 = ReplacementRule(pattern1017, replacement1017)

    pattern1018 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons54, cons397, cons398)
    rule1018 = ReplacementRule(pattern1018, replacement1018)

    pattern1019 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons54, cons397, cons644)
    rule1019 = ReplacementRule(pattern1019, replacement1019)

    pattern1020 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons155)
    rule1020 = ReplacementRule(pattern1020, replacement1020)

    pattern1021 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**S(2)*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons155)
    rule1021 = ReplacementRule(pattern1021, replacement1021)

    pattern1022 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons54, cons65, cons397, cons405)
    rule1022 = ReplacementRule(pattern1022, replacement1022)

    pattern1023 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons52, cons65, cons397, cons644)
    rule1023 = ReplacementRule(pattern1023, replacement1023)

    pattern1024 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons155)
    rule1024 = ReplacementRule(pattern1024, replacement1024)

    pattern1025 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons155)
    rule1025 = ReplacementRule(pattern1025, replacement1025)

    pattern1026 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))**(S(3)/2)*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons155)
    rule1026 = ReplacementRule(pattern1026, replacement1026)

    pattern1027 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/sqrt(e_ + x_**S(2)*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons645)
    rule1027 = ReplacementRule(pattern1027, replacement1027)

    pattern1028 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/sqrt(e_ + x_**S(2)*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons646)
    rule1028 = ReplacementRule(pattern1028, replacement1028)

    pattern1029 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/(e_ + x_**S(2)*WC('f', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons155)
    rule1029 = ReplacementRule(pattern1029, replacement1029)

    pattern1030 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons54, cons150, CustomConstraint(With1030))
    rule1030 = ReplacementRule(pattern1030, replacement1030)

    pattern1031 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons54, cons198)
    rule1031 = ReplacementRule(pattern1031, replacement1031)

    pattern1032 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons52, cons54, cons647)
    rule1032 = ReplacementRule(pattern1032, replacement1032)

    pattern1033 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons4, cons52, cons54, cons648, cons649, cons70, cons71)
    rule1033 = ReplacementRule(pattern1033, replacement1033)

    pattern1034 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons54, cons587, cons588)
    rule1034 = ReplacementRule(pattern1034, replacement1034)

    pattern1035 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons52, cons587, cons40, cons650)
    rule1035 = ReplacementRule(pattern1035, replacement1035)

    pattern1036 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons52, cons54, cons587, cons388)
    rule1036 = ReplacementRule(pattern1036, replacement1036)

    pattern1037 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons654, cons655, cons656, cons657, cons4, cons5, cons52, cons54, cons651, cons652, cons653)
    rule1037 = ReplacementRule(pattern1037, replacement1037)

    pattern1038 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons654, cons655, cons656, cons657, cons4, cons5, cons52, cons54, cons651, cons652)
    rule1038 = ReplacementRule(pattern1038, replacement1038)

    pattern1039 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons52, cons54, cons658, cons502)
    rule1039 = ReplacementRule(pattern1039, replacement1039)

    pattern1040 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons52, cons54, cons658, cons503)
    rule1040 = ReplacementRule(pattern1040, replacement1040)

    pattern1041 = Pattern(Integral((g_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons52, cons54, cons21)
    rule1041 = ReplacementRule(pattern1041, replacement1041)

    pattern1042 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons659)
    rule1042 = ReplacementRule(pattern1042, replacement1042)

    pattern1043 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons54, cons55)
    rule1043 = ReplacementRule(pattern1043, replacement1043)

    pattern1044 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons660, cons504)
    rule1044 = ReplacementRule(pattern1044, replacement1044)

    pattern1045 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons54, cons502)
    rule1045 = ReplacementRule(pattern1045, replacement1045)

    pattern1046 = Pattern(Integral((g_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons52, cons54, cons502)
    rule1046 = ReplacementRule(pattern1046, replacement1046)

    pattern1047 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons54, cons150, cons20, CustomConstraint(With1047))
    rule1047 = ReplacementRule(pattern1047, replacement1047)

    pattern1048 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons52, cons54, cons150, cons369)
    rule1048 = ReplacementRule(pattern1048, With1048)

    pattern1049 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons150, cons404, cons139, cons405, cons661)
    rule1049 = ReplacementRule(pattern1049, replacement1049)

    pattern1050 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons52, cons150, cons246, cons139, cons609)
    rule1050 = ReplacementRule(pattern1050, replacement1050)

    pattern1051 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons52, cons150, cons13, cons139)
    rule1051 = ReplacementRule(pattern1051, replacement1051)

    pattern1052 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons150, cons613, cons405, cons96, cons662)
    rule1052 = ReplacementRule(pattern1052, replacement1052)

    pattern1053 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons150, cons397, cons405, cons662)
    rule1053 = ReplacementRule(pattern1053, replacement1053)

    pattern1054 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons52, cons150, cons33, cons532)
    rule1054 = ReplacementRule(pattern1054, replacement1054)

    pattern1055 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons52, cons150, cons33, cons96)
    rule1055 = ReplacementRule(pattern1055, replacement1055)

    pattern1056 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons150)
    rule1056 = ReplacementRule(pattern1056, replacement1056)

    pattern1057 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons52, cons150)
    rule1057 = ReplacementRule(pattern1057, replacement1057)

    pattern1058 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons52, cons150, cons663)
    rule1058 = ReplacementRule(pattern1058, replacement1058)

    pattern1059 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons54, cons198, cons20)
    rule1059 = ReplacementRule(pattern1059, replacement1059)

    pattern1060 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons52, cons54, cons198, cons369)
    rule1060 = ReplacementRule(pattern1060, With1060)

    pattern1061 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons52, cons54, cons198, cons358)
    rule1061 = ReplacementRule(pattern1061, replacement1061)

    pattern1062 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons52, cons54, cons491)
    rule1062 = ReplacementRule(pattern1062, With1062)

    pattern1063 = Pattern(Integral((g_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons52, cons54, cons491)
    rule1063 = ReplacementRule(pattern1063, replacement1063)

    pattern1064 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons54, cons543)
    rule1064 = ReplacementRule(pattern1064, replacement1064)

    pattern1065 = Pattern(Integral((g_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons52, cons54, cons543)
    rule1065 = ReplacementRule(pattern1065, replacement1065)

    pattern1066 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons404, cons139, cons405, cons661)
    rule1066 = ReplacementRule(pattern1066, replacement1066)

    pattern1067 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons52, cons13, cons139)
    rule1067 = ReplacementRule(pattern1067, replacement1067)

    pattern1068 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons397, cons405, cons662)
    rule1068 = ReplacementRule(pattern1068, replacement1068)

    pattern1069 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons382)
    rule1069 = ReplacementRule(pattern1069, replacement1069)

    pattern1070 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons52, cons664)
    rule1070 = ReplacementRule(pattern1070, replacement1070)

    pattern1071 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons54, cons587, cons588)
    rule1071 = ReplacementRule(pattern1071, replacement1071)

    pattern1072 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons52, cons587, cons40, cons650)
    rule1072 = ReplacementRule(pattern1072, replacement1072)

    pattern1073 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons54, cons587, cons388)
    rule1073 = ReplacementRule(pattern1073, replacement1073)

    pattern1074 = Pattern(Integral((g_*x_)**m_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons52, cons54, cons587)
    rule1074 = ReplacementRule(pattern1074, replacement1074)

    pattern1075 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons52, cons54, cons665)
    rule1075 = ReplacementRule(pattern1075, replacement1075)

    pattern1076 = Pattern(Integral(u_**WC('m', S(1))*(e_ + v_**n_*WC('f', S(1)))**WC('r', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons54, cons556)
    rule1076 = ReplacementRule(pattern1076, replacement1076)

    pattern1077 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons654, cons655, cons656, cons657, cons210, cons19, cons4, cons5, cons52, cons54, cons651, cons652, cons653)
    rule1077 = ReplacementRule(pattern1077, replacement1077)

    pattern1078 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons29, cons654, cons655, cons656, cons657, cons210, cons19, cons4, cons5, cons52, cons54, cons651, cons652)
    rule1078 = ReplacementRule(pattern1078, replacement1078)
    return [rule692, rule693, rule694, rule695, rule696, rule697, rule698, rule699, rule700, rule701, rule702, rule703, rule704, rule705, rule706, rule707, rule708, rule709, rule710, rule711, rule712, rule713, rule714, rule715, rule716, rule717, rule718, rule719, rule720, rule721, rule722, rule723, rule724, rule725, rule726, rule727, rule728, rule729, rule730, rule731, rule732, rule733, rule734, rule735, rule736, rule737, rule738, rule739, rule740, rule741, rule742, rule743, rule744, rule745, rule746, rule747, rule748, rule749, rule750, rule751, rule752, rule753, rule754, rule755, rule756, rule757, rule758, rule759, rule760, rule761, rule762, rule763, rule764, rule765, rule766, rule767, rule768, rule769, rule770, rule771, rule772, rule773, rule774, rule775, rule776, rule777, rule778, rule779, rule780, rule781, rule782, rule783, rule784, rule785, rule786, rule787, rule788, rule789, rule790, rule791, rule792, rule793, rule794, rule795, rule796, rule797, rule798, rule799, rule800, rule801, rule802, rule803, rule804, rule805, rule806, rule807, rule808, rule809, rule810, rule811, rule812, rule813, rule814, rule815, rule816, rule817, rule818, rule819, rule820, rule821, rule822, rule823, rule824, rule825, rule826, rule827, rule828, rule829, rule830, rule831, rule832, rule833, rule834, rule835, rule836, rule837, rule838, rule839, rule840, rule841, rule842, rule843, rule844, rule845, rule846, rule847, rule848, rule849, rule850, rule851, rule852, rule853, rule854, rule855, rule856, rule857, rule858, rule859, rule860, rule861, rule862, rule863, rule864, rule865, rule866, rule867, rule868, rule869, rule870, rule871, rule872, rule873, rule874, rule875, rule876, rule877, rule878, rule879, rule880, rule881, rule882, rule883, rule884, rule885, rule886, rule887, rule888, rule889, rule890, rule891, rule892, rule893, rule894, rule895, rule896, rule897, rule898, rule899, rule900, rule901, rule902, rule903, rule904, rule905, rule906, rule907, rule908, rule909, rule910, rule911, rule912, rule913, rule914, rule915, rule916, rule917, rule918, rule919, rule920, rule921, rule922, rule923, rule924, rule925, rule926, rule927, rule928, rule929, rule930, rule931, rule932, rule933, rule934, rule935, rule936, rule937, rule938, rule939, rule940, rule941, rule942, rule943, rule944, rule945, rule946, rule947, rule948, rule949, rule950, rule951, rule952, rule953, rule954, rule955, rule956, rule957, rule958, rule959, rule960, rule961, rule962, rule963, rule964, rule965, rule966, rule967, rule968, rule969, rule970, rule971, rule972, rule973, rule974, rule975, rule976, rule977, rule978, rule979, rule980, rule981, rule982, rule983, rule984, rule985, rule986, rule987, rule988, rule989, rule990, rule991, rule992, rule993, rule994, rule995, rule996, rule997, rule998, rule999, rule1000, rule1001, rule1002, rule1003, rule1004, rule1005, rule1006, rule1007, rule1008, rule1009, rule1010, rule1011, rule1012, rule1013, rule1014, rule1015, rule1016, rule1017, rule1018, rule1019, rule1020, rule1021, rule1022, rule1023, rule1024, rule1025, rule1026, rule1027, rule1028, rule1029, rule1030, rule1031, rule1032, rule1033, rule1034, rule1035, rule1036, rule1037, rule1038, rule1039, rule1040, rule1041, rule1042, rule1043, rule1044, rule1045, rule1046, rule1047, rule1048, rule1049, rule1050, rule1051, rule1052, rule1053, rule1054, rule1055, rule1056, rule1057, rule1058, rule1059, rule1060, rule1061, rule1062, rule1063, rule1064, rule1065, rule1066, rule1067, rule1068, rule1069, rule1070, rule1071, rule1072, rule1073, rule1074, rule1075, rule1076, rule1077, rule1078, ]





def replacement692(b, n, p, x):
    return Dist(b**IntPart(p)*x**(-n*FracPart(p))*(b*x**n)**FracPart(p), Int(x**(n*p), x), x)


def replacement693(a, b, n, p, x):
    return Simp(x*(a + b*x**n)**(p + S(1))/a, x)


def replacement694(a, b, n, p, x):
    return Dist((n*(p + S(1)) + S(1))/(a*n*(p + S(1))), Int((a + b*x**n)**(p + S(1)), x), x) - Simp(x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))), x)


def replacement695(a, b, n, x):
    return Int(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n), x)


def replacement696(a, b, n, p, x):
    return Int(x**(n*p)*(a*x**(-n) + b)**p, x)


def replacement697(a, b, n, p, x):
    return Int(ExpandIntegrand((a + b*x**n)**p, x), x)


def replacement698(a, b, n, p, x):
    return Dist(a*n*p/(n*p + S(1)), Int((a + b*x**n)**(p + S(-1)), x), x) + Simp(x*(a + b*x**n)**p/(n*p + S(1)), x)


def replacement699(a, b, x):
    return Simp(S(2)*EllipticE(ArcTan(x*Rt(b/a, S(2)))/S(2), S(2))/(a**(S(5)/4)*Rt(b/a, S(2))), x)


def replacement700(a, b, x):
    return Dist((S(1) + b*x**S(2)/a)**(S(1)/4)/(a*(a + b*x**S(2))**(S(1)/4)), Int((S(1) + b*x**S(2)/a)**(S(-5)/4), x), x)


def replacement701(a, b, x):
    return Dist(S(1)/((a/(a + b*x**S(2)))**(S(2)/3)*(a + b*x**S(2))**(S(2)/3)), Subst(Int((-b*x**S(2) + S(1))**(S(-1)/3), x), x, x/sqrt(a + b*x**S(2))), x)


def replacement702(a, b, n, p, x):
    return Dist((n*(p + S(1)) + S(1))/(a*n*(p + S(1))), Int((a + b*x**n)**(p + S(1)), x), x) - Simp(x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))), x)


def replacement703(a, b, x):
    return Dist(S(1)/(S(3)*Rt(a, S(3))**S(2)), Int((-x*Rt(b, S(3)) + S(2)*Rt(a, S(3)))/(x**S(2)*Rt(b, S(3))**S(2) - x*Rt(a, S(3))*Rt(b, S(3)) + Rt(a, S(3))**S(2)), x), x) + Dist(S(1)/(S(3)*Rt(a, S(3))**S(2)), Int(S(1)/(x*Rt(b, S(3)) + Rt(a, S(3))), x), x)


def With704(a, b, n, x):
    r = Numerator(Rt(a/b, n))
    s = Denominator(Rt(a/b, n))
    k = Symbol('k')
    u = Symbol('u')
    u = Int((r - s*x*cos(Pi*(S(2)*k + S(-1))/n))/(r**S(2) - S(2)*r*s*x*cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)
    u = Int((r - s*x*cos(Pi*(2*k - 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
    return Simp(Dist(2*r/(a*n), Sum_doit(u, List(k, 1, n/2 - 1/2)), x) + r*Int(1/(r + s*x), x)/(a*n), x)


def With705(a, b, n, x):
    r = Numerator(Rt(-a/b, n))
    s = Denominator(Rt(-a/b, n))
    k = Symbol('k')
    u = Symbol('u')
    u = Int((r + s*x*cos(Pi*(S(2)*k + S(-1))/n))/(r**S(2) + S(2)*r*s*x*cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)
    u = Int((r + s*x*cos(Pi*(2*k - 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
    return Simp(Dist(2*r/(a*n), Sum_doit(u, List(k, 1, n/2 - 1/2)), x) + r*Int(1/(r - s*x), x)/(a*n), x)


def replacement706(a, b, x):
    return Simp(ArcTan(x*Rt(b, S(2))/Rt(a, S(2)))/(Rt(a, S(2))*Rt(b, S(2))), x)


def replacement707(a, b, x):
    return -Simp(ArcTan(x*Rt(-b, S(2))/Rt(-a, S(2)))/(Rt(-a, S(2))*Rt(-b, S(2))), x)


def replacement708(a, b, x):
    return Simp(ArcTan(x/Rt(a/b, S(2)))*Rt(a/b, S(2))/a, x)


def replacement709(a, b, x):
    return Simp(atanh(x*Rt(-b, S(2))/Rt(a, S(2)))/(Rt(a, S(2))*Rt(-b, S(2))), x)


def replacement710(a, b, x):
    return -Simp(atanh(x*Rt(b, S(2))/Rt(-a, S(2)))/(Rt(-a, S(2))*Rt(b, S(2))), x)


def replacement711(a, b, x):
    return Simp(Rt(-a/b, S(2))*atanh(x/Rt(-a/b, S(2)))/a, x)


def With712(a, b, n, x):
    r = Numerator(Rt(a/b, n))
    s = Denominator(Rt(a/b, n))
    k = Symbol('k')
    u = Symbol('u')
    v = Symbol('v')
    u = Int((r - s*x*cos(Pi*(S(2)*k + S(-1))/n))/(r**S(2) - S(2)*r*s*x*cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x) + Int((r + s*x*cos(Pi*(S(2)*k + S(-1))/n))/(r**S(2) + S(2)*r*s*x*cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)
    u = Int((r - s*x*cos(Pi*(2*k - 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x) + Int((r + s*x*cos(Pi*(2*k - 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
    return Simp(Dist(2*r/(a*n), Sum_doit(u, List(k, 1, n/4 - 1/2)), x) + 2*r**2*Int(1/(r**2 + s**2*x**2), x)/(a*n), x)


def With713(a, b, n, x):
    r = Numerator(Rt(-a/b, n))
    s = Denominator(Rt(-a/b, n))
    k = Symbol('k')
    u = Symbol('u')
    u = Int((r - s*x*cos(S(2)*Pi*k/n))/(r**S(2) - S(2)*r*s*x*cos(S(2)*Pi*k/n) + s**S(2)*x**S(2)), x) + Int((r + s*x*cos(S(2)*Pi*k/n))/(r**S(2) + S(2)*r*s*x*cos(S(2)*Pi*k/n) + s**S(2)*x**S(2)), x)
    u = Int((r - s*x*cos(2*Pi*k/n))/(r**2 - 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x) + Int((r + s*x*cos(2*Pi*k/n))/(r**2 + 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x)
    return Simp(Dist(2*r/(a*n), Sum_doit(u, List(k, 1, n/4 - 1/2)), x) + 2*r**2*Int(1/(r**2 - s**2*x**2), x)/(a*n), x)


def With714(a, b, x):
    r = Numerator(Rt(a/b, S(2)))
    s = Denominator(Rt(a/b, S(2)))
    return Dist(S(1)/(S(2)*r), Int((r - s*x**S(2))/(a + b*x**S(4)), x), x) + Dist(S(1)/(S(2)*r), Int((r + s*x**S(2))/(a + b*x**S(4)), x), x)


def With715(a, b, x):
    r = Numerator(Rt(-a/b, S(2)))
    s = Denominator(Rt(-a/b, S(2)))
    return Dist(r/(S(2)*a), Int(S(1)/(r - s*x**S(2)), x), x) + Dist(r/(S(2)*a), Int(S(1)/(r + s*x**S(2)), x), x)


def With716(a, b, n, x):
    r = Numerator(Rt(a/b, S(4)))
    s = Denominator(Rt(a/b, S(4)))
    return Dist(sqrt(S(2))*r/(S(4)*a), Int((sqrt(S(2))*r - s*x**(n/S(4)))/(r**S(2) - sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x), x) + Dist(sqrt(S(2))*r/(S(4)*a), Int((sqrt(S(2))*r + s*x**(n/S(4)))/(r**S(2) + sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x), x)


def With717(a, b, n, x):
    r = Numerator(Rt(-a/b, S(2)))
    s = Denominator(Rt(-a/b, S(2)))
    return Dist(r/(S(2)*a), Int(S(1)/(r - s*x**(n/S(2))), x), x) + Dist(r/(S(2)*a), Int(S(1)/(r + s*x**(n/S(2))), x), x)


def replacement718(a, b, x):
    return Simp(asinh(x*Rt(b, S(2))/sqrt(a))/Rt(b, S(2)), x)


def replacement719(a, b, x):
    return Simp(asin(x*Rt(-b, S(2))/sqrt(a))/Rt(-b, S(2)), x)


def replacement720(a, b, x):
    return Subst(Int(S(1)/(-b*x**S(2) + S(1)), x), x, x/sqrt(a + b*x**S(2)))


def With721(a, b, x):
    r = Numer(Rt(b/a, S(3)))
    s = Denom(Rt(b/a, S(3)))
    return Simp(S(2)*S(3)**(S(3)/4)*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(sqrt(S(3)) + S(2))*(r*x + s)*EllipticF(asin((r*x + s*(S(1) - sqrt(S(3))))/(r*x + s*(S(1) + sqrt(S(3))))), S(-7) - S(4)*sqrt(S(3)))/(S(3)*r*sqrt(s*(r*x + s)/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(a + b*x**S(3))), x)


def With722(a, b, x):
    r = Numer(Rt(b/a, S(3)))
    s = Denom(Rt(b/a, S(3)))
    return Simp(S(2)*S(3)**(S(3)/4)*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(S(1) - sqrt(S(3))))**S(2))*sqrt(S(2) - sqrt(S(3)))*(r*x + s)*EllipticF(asin((r*x + s*(S(1) + sqrt(S(3))))/(r*x + s*(S(1) - sqrt(S(3))))), S(-7) + S(4)*sqrt(S(3)))/(S(3)*r*sqrt(-s*(r*x + s)/(r*x + s*(S(1) - sqrt(S(3))))**S(2))*sqrt(a + b*x**S(3))), x)


def With723(a, b, x):
    q = Rt(b/a, S(4))
    return Simp(sqrt((a + b*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticF(S(2)*ArcTan(q*x), S(1)/2)/(S(2)*q*sqrt(a + b*x**S(4))), x)


def replacement724(a, b, x):
    return Simp(EllipticF(asin(x*Rt(-b, S(4))/Rt(a, S(4))), S(-1))/(Rt(a, S(4))*Rt(-b, S(4))), x)


def With725(a, b, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-a*b, S(2))
    if IntegerQ(q):
        return True
    return False


def replacement725(a, b, x):

    q = Rt(-a*b, S(2))
    return Simp(sqrt(S(2))*sqrt((a + q*x**S(2))/q)*sqrt(-a + q*x**S(2))*EllipticF(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(S(2)*sqrt(-a)*sqrt(a + b*x**S(4))), x)


def With726(a, b, x):
    q = Rt(-a*b, S(2))
    return Simp(sqrt(S(2))*sqrt((a + q*x**S(2))/q)*sqrt((a - q*x**S(2))/(a + q*x**S(2)))*EllipticF(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(S(2)*sqrt(a/(a + q*x**S(2)))*sqrt(a + b*x**S(4))), x)


def replacement727(a, b, x):
    return Dist(sqrt(S(1) + b*x**S(4)/a)/sqrt(a + b*x**S(4)), Int(S(1)/sqrt(S(1) + b*x**S(4)/a), x), x)


def With728(a, b, x):
    r = Numer(Rt(b/a, S(3)))
    s = Denom(Rt(b/a, S(3)))
    return Simp(S(3)**(S(3)/4)*x*sqrt((r**S(2)*x**S(4) - r*s*x**S(2) + s**S(2))/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*(r*x**S(2) + s)*EllipticF(acos((r*x**S(2)*(S(1) - sqrt(S(3))) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)), sqrt(S(3))/S(4) + S(1)/2)/(S(6)*s*sqrt(r*x**S(2)*(r*x**S(2) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*sqrt(a + b*x**S(6))), x)


def replacement729(a, b, x):
    return Dist(S(1)/2, Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x), x) + Dist(S(1)/2, Int((x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x), x)


def replacement730(a, b, x):
    return -Dist(a, Int((a + b*x**S(2))**(S(-5)/4), x), x) + Simp(S(2)*x/(a + b*x**S(2))**(S(1)/4), x)


def replacement731(a, b, x):
    return Simp(S(2)*EllipticE(asin(x*Rt(-b/a, S(2)))/S(2), S(2))/(a**(S(1)/4)*Rt(-b/a, S(2))), x)


def replacement732(a, b, x):
    return Dist((S(1) + b*x**S(2)/a)**(S(1)/4)/(a + b*x**S(2))**(S(1)/4), Int((S(1) + b*x**S(2)/a)**(S(-1)/4), x), x)


def replacement733(a, b, x):
    return Simp(S(2)*EllipticF(ArcTan(x*Rt(b/a, S(2)))/S(2), S(2))/(a**(S(3)/4)*Rt(b/a, S(2))), x)


def replacement734(a, b, x):
    return Simp(S(2)*EllipticF(asin(x*Rt(-b/a, S(2)))/S(2), S(2))/(a**(S(3)/4)*Rt(-b/a, S(2))), x)


def replacement735(a, b, x):
    return Dist((S(1) + b*x**S(2)/a)**(S(3)/4)/(a + b*x**S(2))**(S(3)/4), Int((S(1) + b*x**S(2)/a)**(S(-3)/4), x), x)


def replacement736(a, b, x):
    return Dist(S(3)*sqrt(b*x**S(2))/(S(2)*b*x), Subst(Int(x/sqrt(-a + x**S(3)), x), x, (a + b*x**S(2))**(S(1)/3)), x)


def replacement737(a, b, x):
    return Dist(S(3)*sqrt(b*x**S(2))/(S(2)*b*x), Subst(Int(S(1)/sqrt(-a + x**S(3)), x), x, (a + b*x**S(2))**(S(1)/3)), x)


def replacement738(a, b, x):
    return Dist(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(3)/4)/(a + b*x**S(4))**(S(3)/4), Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(3)/4)), x), x)


def replacement739(a, b, x):
    return -Dist(a/S(2), Int((a + b*x**S(2))**(S(-7)/6), x), x) + Simp(S(3)*x/(S(2)*(a + b*x**S(2))**(S(1)/6)), x)


def replacement740(a, b, n, p, x):
    return Dist(a**(p + S(1)/n), Subst(Int((-b*x**n + S(1))**(-p + S(-1) - S(1)/n), x), x, x*(a + b*x**n)**(-S(1)/n)), x)


def replacement741(a, b, n, p, x):
    return Dist((a/(a + b*x**n))**(p + S(1)/n)*(a + b*x**n)**(p + S(1)/n), Subst(Int((-b*x**n + S(1))**(-p + S(-1) - S(1)/n), x), x, x*(a + b*x**n)**(-S(1)/n)), x)


def replacement742(a, b, n, p, x):
    return -Subst(Int((a + b*x**(-n))**p/x**S(2), x), x, S(1)/x)


def With743(a, b, n, p, x):
    k = Denominator(n)
    return Dist(k, Subst(Int(x**(k + S(-1))*(a + b*x**(k*n))**p, x), x, x**(S(1)/k)), x)


def replacement744(a, b, n, p, x):
    return Int(ExpandIntegrand((a + b*x**n)**p, x), x)


def replacement745(a, b, n, p, x):
    return Simp(a**p*x*Hypergeometric2F1(-p, S(1)/n, S(1) + S(1)/n, -b*x**n/a), x)


def replacement746(a, b, n, p, x):
    return Dist(a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p), Int((S(1) + b*x**n/a)**p, x), x)


def replacement747(a, b, n, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*x**n)**p, x), x, u), x)


def replacement748(a1, a2, b1, b2, n, p, x):
    return Int((a1*a2 + b1*b2*x**(S(2)*n))**p, x)


def replacement749(a1, a2, b1, b2, n, p, x):
    return Dist(S(2)*a1*a2*n*p/(S(2)*n*p + S(1)), Int((a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x), x) + Simp(x*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(S(2)*n*p + S(1)), x)


def replacement750(a1, a2, b1, b2, n, p, x):
    return Dist((S(2)*n*(p + S(1)) + S(1))/(S(2)*a1*a2*n*(p + S(1))), Int((a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x), x) - Simp(x*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*n*(p + S(1))), x)


def replacement751(a1, a2, b1, b2, n, p, x):
    return -Subst(Int((a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p/x**S(2), x), x, S(1)/x)


def With752(a1, a2, b1, b2, n, p, x):
    k = Denominator(S(2)*n)
    return Dist(k, Subst(Int(x**(k + S(-1))*(a1 + b1*x**(k*n))**p*(a2 + b2*x**(k*n))**p, x), x, x**(S(1)/k)), x)


def replacement753(a1, a2, b1, b2, n, p, x):
    return Dist((a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p)), Int((a1*a2 + b1*b2*x**(S(2)*n))**p, x), x)


def replacement754(a1, a2, b1, b2, c, m, n, p, x):
    return Int((c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x)


def replacement755(b, c, m, n, p, x):
    return Dist(b**(S(1) - (m + S(1))/n)*c**m/n, Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n), x), x, x**n), x)


def replacement756(b, c, m, n, p, x):
    return Dist(b**IntPart(p)*c**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p), Int(x**(m + n*p), x), x)


def replacement757(b, c, m, n, p, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(b*x**n)**p, x), x)


def replacement758(a, b, m, n, p, x):
    return Int(x**(m + n*p)*(a*x**(-n) + b)**p, x)


def replacement759(a, b, c, m, n, p, x):
    return Simp((c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))), x)


def replacement760(a1, a2, b1, b2, c, m, n, p, x):
    return Simp((c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))), x)


def replacement761(a, b, m, n, p, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p, x), x, x**n), x)


def replacement762(a1, a2, b1, b2, m, n, p, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a1 + b1*x)**p*(a2 + b2*x)**p, x), x, x**n), x)


def replacement763(a, b, c, m, n, p, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a + b*x**n)**p, x), x)


def replacement764(a1, a2, b1, b2, c, m, n, p, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x), x)


def replacement765(a, b, c, m, n, p, x):
    return Int(ExpandIntegrand((c*x)**m*(a + b*x**n)**p, x), x)


def replacement766(a, b, m, n, p, x):
    return -Dist(b*(m + n*(p + S(1)) + S(1))/(a*(m + S(1))), Int(x**(m + n)*(a + b*x**n)**p, x), x) + Simp(x**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*(m + S(1))), x)


def replacement767(a1, a2, b1, b2, m, n, p, x):
    return -Dist(b1*b2*(m + S(2)*n*(p + S(1)) + S(1))/(a1*a2*(m + S(1))), Int(x**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x), x) + Simp(x**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*(m + S(1))), x)


def replacement768(a, b, c, m, n, p, x):
    return Dist((m + n*(p + S(1)) + S(1))/(a*n*(p + S(1))), Int((c*x)**m*(a + b*x**n)**(p + S(1)), x), x) - Simp((c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))), x)


def replacement769(a1, a2, b1, b2, c, m, n, p, x):
    return Dist((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*a1*a2*n*(p + S(1))), Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x), x) - Simp((c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))), x)


def With770(a, b, m, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    k = GCD(m + S(1), n)
    if Unequal(k, S(1)):
        return True
    return False


def replacement770(a, b, m, n, p, x):

    k = GCD(m + S(1), n)
    return Dist(S(1)/k, Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p, x), x, x**k), x)


def With771(a1, a2, b1, b2, m, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    k = GCD(m + S(1), S(2)*n)
    if Unequal(k, S(1)):
        return True
    return False


def replacement771(a1, a2, b1, b2, m, n, p, x):

    k = GCD(m + S(1), S(2)*n)
    return Dist(S(1)/k, Subst(Int(x**(S(-1) + (m + S(1))/k)*(a1 + b1*x**(n/k))**p*(a2 + b2*x**(n/k))**p, x), x, x**k), x)


def replacement772(a, b, c, m, n, p, x):
    return -Dist(b*c**(-n)*n*p/(m + S(1)), Int((c*x)**(m + n)*(a + b*x**n)**(p + S(-1)), x), x) + Simp((c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + S(1))), x)


def replacement773(a1, a2, b1, b2, c, m, n, p, x):
    return Dist(S(2)*a1*a2*n*p/(m + S(2)*n*p + S(1)), Int((c*x)**m*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x), x) + Simp((c*x)**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(c*(m + S(2)*n*p + S(1))), x)


def replacement774(a, b, c, m, n, p, x):
    return Dist(a*n*p/(m + n*p + S(1)), Int((c*x)**m*(a + b*x**n)**(p + S(-1)), x), x) + Simp((c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + n*p + S(1))), x)


def replacement775(a, b, x):
    return Dist(x*(a/(b*x**S(4)) + S(1))**(S(1)/4)/(b*(a + b*x**S(4))**(S(1)/4)), Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(5)/4)), x), x)


def replacement776(a, b, m, x):
    return -Dist(a*(m + S(-3))/(b*(m + S(-4))), Int(x**(m + S(-4))/(a + b*x**S(4))**(S(5)/4), x), x) + Simp(x**(m + S(-3))/(b*(a + b*x**S(4))**(S(1)/4)*(m + S(-4))), x)


def replacement777(a, b, m, x):
    return -Dist(b*m/(a*(m + S(1))), Int(x**(m + S(4))/(a + b*x**S(4))**(S(5)/4), x), x) + Simp(x**(m + S(1))/(a*(a + b*x**S(4))**(S(1)/4)*(m + S(1))), x)


def replacement778(a, b, c, x):
    return Dist(sqrt(c*x)*(a/(b*x**S(2)) + S(1))**(S(1)/4)/(b*(a + b*x**S(2))**(S(1)/4)), Int(S(1)/(x**S(2)*(a/(b*x**S(2)) + S(1))**(S(5)/4)), x), x)


def replacement779(a, b, c, m, x):
    return -Dist(S(2)*a*c**S(2)*(m + S(-1))/(b*(S(2)*m + S(-3))), Int((c*x)**(m + S(-2))/(a + b*x**S(2))**(S(5)/4), x), x) + Simp(S(2)*c*(c*x)**(m + S(-1))/(b*(a + b*x**S(2))**(S(1)/4)*(S(2)*m + S(-3))), x)


def replacement780(a, b, c, m, x):
    return -Dist(b*(S(2)*m + S(1))/(S(2)*a*c**S(2)*(m + S(1))), Int((c*x)**(m + S(2))/(a + b*x**S(2))**(S(5)/4), x), x) + Simp((c*x)**(m + S(1))/(a*c*(a + b*x**S(2))**(S(1)/4)*(m + S(1))), x)


def replacement781(a, b, x):
    return -Dist(S(1)/b, Int(S(1)/(x**S(2)*(a + b*x**S(4))**(S(1)/4)), x), x) - Simp(S(1)/(b*x*(a + b*x**S(4))**(S(1)/4)), x)


def replacement782(a, b, c, m, n, p, x):
    return -Dist(c**n*(m - n + S(1))/(b*n*(p + S(1))), Int((c*x)**(m - n)*(a + b*x**n)**(p + S(1)), x), x) + Simp(c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)


def replacement783(a1, a2, b1, b2, c, m, n, p, x):
    return -Dist(c**(S(2)*n)*(m - S(2)*n + S(1))/(S(2)*b1*b2*n*(p + S(1))), Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x), x) + Simp(c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*b1*b2*n*(p + S(1))), x)


def replacement784(a, b, c, m, n, p, x):
    return Dist((m + n*(p + S(1)) + S(1))/(a*n*(p + S(1))), Int((c*x)**m*(a + b*x**n)**(p + S(1)), x), x) - Simp((c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))), x)


def replacement785(a1, a2, b1, b2, c, m, n, p, x):
    return Dist((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*a1*a2*n*(p + S(1))), Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x), x) - Simp((c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))), x)


def replacement786(a, b, x):
    return Dist(S(1)/(S(3)*Rt(a, S(3))*Rt(b, S(3))), Int((x*Rt(b, S(3)) + Rt(a, S(3)))/(x**S(2)*Rt(b, S(3))**S(2) - x*Rt(a, S(3))*Rt(b, S(3)) + Rt(a, S(3))**S(2)), x), x) - Dist(S(1)/(S(3)*Rt(a, S(3))*Rt(b, S(3))), Int(S(1)/(x*Rt(b, S(3)) + Rt(a, S(3))), x), x)


def With787(a, b, m, n, x):
    r = Numerator(Rt(a/b, n))
    s = Denominator(Rt(a/b, n))
    k = Symbol('k')
    u = Symbol('u')
    u = Int((r*cos(Pi*m*(S(2)*k + S(-1))/n) - s*x*cos(Pi*(S(2)*k + S(-1))*(m + S(1))/n))/(r**S(2) - S(2)*r*s*x*cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)
    u = Int((r*cos(Pi*m*(2*k - 1)/n) - s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
    return Simp(Dist(2*r**(m + 1)*s**(-m)/(a*n), Sum_doit(u, List(k, 1, n/2 - 1/2)), x) - s**(-m)*(-r)**(m + 1)*Int(1/(r + s*x), x)/(a*n), x)


def With788(a, b, m, n, x):
    r = Numerator(Rt(-a/b, n))
    s = Denominator(Rt(-a/b, n))
    k = Symbol('k')
    u = Symbol('u')
    u = Int((r*cos(Pi*m*(S(2)*k + S(-1))/n) + s*x*cos(Pi*(S(2)*k + S(-1))*(m + S(1))/n))/(r**S(2) + S(2)*r*s*x*cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)
    u = Int((r*cos(Pi*m*(2*k - 1)/n) + s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
    return Simp(-Dist(2*s**(-m)*(-r)**(m + 1)/(a*n), Sum_doit(u, List(k, 1, n/2 - 1/2)), x) + r**(m + 1)*s**(-m)*Int(1/(r - s*x), x)/(a*n), x)


def With789(a, b, m, n, x):
    r = Numerator(Rt(a/b, n))
    s = Denominator(Rt(a/b, n))
    k = Symbol('k')
    u = Symbol('u')
    u = Int((r*cos(Pi*m*(S(2)*k + S(-1))/n) - s*x*cos(Pi*(S(2)*k + S(-1))*(m + S(1))/n))/(r**S(2) - S(2)*r*s*x*cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x) + Int((r*cos(Pi*m*(S(2)*k + S(-1))/n) + s*x*cos(Pi*(S(2)*k + S(-1))*(m + S(1))/n))/(r**S(2) + S(2)*r*s*x*cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)
    u = Int((r*cos(Pi*m*(2*k - 1)/n) - s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x) + Int((r*cos(Pi*m*(2*k - 1)/n) + s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
    return Simp(2*(-1)**(m/2)*r**(m + 2)*s**(-m)*Int(1/(r**2 + s**2*x**2), x)/(a*n) + Dist(2*r**(m + 1)*s**(-m)/(a*n), Sum_doit(u, List(k, 1, n/4 - 1/2)), x), x)


def With790(a, b, m, n, x):
    r = Numerator(Rt(-a/b, n))
    s = Denominator(Rt(-a/b, n))
    k = Symbol('k')
    u = Symbol('u')
    u = Int((r*cos(S(2)*Pi*k*m/n) - s*x*cos(S(2)*Pi*k*(m + S(1))/n))/(r**S(2) - S(2)*r*s*x*cos(S(2)*Pi*k/n) + s**S(2)*x**S(2)), x) + Int((r*cos(S(2)*Pi*k*m/n) + s*x*cos(S(2)*Pi*k*(m + S(1))/n))/(r**S(2) + S(2)*r*s*x*cos(S(2)*Pi*k/n) + s**S(2)*x**S(2)), x)
    u = Int((r*cos(2*Pi*k*m/n) - s*x*cos(2*Pi*k*(m + 1)/n))/(r**2 - 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x) + Int((r*cos(2*Pi*k*m/n) + s*x*cos(2*Pi*k*(m + 1)/n))/(r**2 + 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x)
    return Simp(Dist(2*r**(m + 1)*s**(-m)/(a*n), Sum_doit(u, List(k, 1, n/4 - 1/2)), x) + 2*r**(m + 2)*s**(-m)*Int(1/(r**2 - s**2*x**2), x)/(a*n), x)


def With791(a, b, x):
    r = Numerator(Rt(a/b, S(2)))
    s = Denominator(Rt(a/b, S(2)))
    return -Dist(S(1)/(S(2)*s), Int((r - s*x**S(2))/(a + b*x**S(4)), x), x) + Dist(S(1)/(S(2)*s), Int((r + s*x**S(2))/(a + b*x**S(4)), x), x)


def With792(a, b, x):
    r = Numerator(Rt(-a/b, S(2)))
    s = Denominator(Rt(-a/b, S(2)))
    return -Dist(s/(S(2)*b), Int(S(1)/(r - s*x**S(2)), x), x) + Dist(s/(S(2)*b), Int(S(1)/(r + s*x**S(2)), x), x)


def With793(a, b, m, n, x):
    r = Numerator(Rt(a/b, S(4)))
    s = Denominator(Rt(a/b, S(4)))
    return Dist(sqrt(S(2))*s**S(3)/(S(4)*b*r), Int(x**(m - n/S(4))/(r**S(2) - sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x), x) - Dist(sqrt(S(2))*s**S(3)/(S(4)*b*r), Int(x**(m - n/S(4))/(r**S(2) + sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x), x)


def With794(a, b, m, n, x):
    r = Numerator(Rt(-a/b, S(2)))
    s = Denominator(Rt(-a/b, S(2)))
    return Dist(r/(S(2)*a), Int(x**m/(r - s*x**(n/S(2))), x), x) + Dist(r/(S(2)*a), Int(x**m/(r + s*x**(n/S(2))), x), x)


def With795(a, b, m, n, x):
    r = Numerator(Rt(-a/b, S(2)))
    s = Denominator(Rt(-a/b, S(2)))
    return -Dist(s/(S(2)*b), Int(x**(m - n/S(2))/(r - s*x**(n/S(2))), x), x) + Dist(s/(S(2)*b), Int(x**(m - n/S(2))/(r + s*x**(n/S(2))), x), x)


def replacement796(a, b, m, n, x):
    return Int(PolynomialDivide(x**m, a + b*x**n, x), x)


def With797(a, b, x):
    r = Numer(Rt(b/a, S(3)))
    s = Denom(Rt(b/a, S(3)))
    return Dist(S(1)/r, Int((r*x + s*(S(1) - sqrt(S(3))))/sqrt(a + b*x**S(3)), x), x) + Dist(sqrt(S(2))*s/(r*sqrt(sqrt(S(3)) + S(2))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)


def With798(a, b, x):
    r = Numer(Rt(b/a, S(3)))
    s = Denom(Rt(b/a, S(3)))
    return Dist(S(1)/r, Int((r*x + s*(S(1) + sqrt(S(3))))/sqrt(a + b*x**S(3)), x), x) - Dist(sqrt(S(2))*s/(r*sqrt(S(2) - sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)


def With799(a, b, x):
    q = Rt(b/a, S(2))
    return -Dist(S(1)/q, Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x), x) + Dist(S(1)/q, Int(S(1)/sqrt(a + b*x**S(4)), x), x)


def With800(a, b, x):
    q = Rt(-b/a, S(2))
    return -Dist(S(1)/q, Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x), x) + Dist(S(1)/q, Int(S(1)/sqrt(a + b*x**S(4)), x), x)


def With801(a, b, x):
    q = Rt(-b/a, S(2))
    return Dist(S(1)/q, Int((q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x), x) - Dist(S(1)/q, Int(S(1)/sqrt(a + b*x**S(4)), x), x)


def With802(a, b, x):
    r = Numer(Rt(b/a, S(3)))
    s = Denom(Rt(b/a, S(3)))
    return -Dist(S(1)/(S(2)*r**S(2)), Int((-S(2)*r**S(2)*x**S(4) + s**S(2)*(S(-1) + sqrt(S(3))))/sqrt(a + b*x**S(6)), x), x) + Dist(s**S(2)*(S(-1) + sqrt(S(3)))/(S(2)*r**S(2)), Int(S(1)/sqrt(a + b*x**S(6)), x), x)


def replacement803(a, b, x):
    return -Dist(S(1)/(S(2)*Rt(b/a, S(4))), Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x), x) + Dist(S(1)/(S(2)*Rt(b/a, S(4))), Int((x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x), x)


def replacement804(a, b, x):
    return -Dist(a/S(2), Int(x**S(2)/(a + b*x**S(4))**(S(5)/4), x), x) + Simp(x**S(3)/(S(2)*(a + b*x**S(4))**(S(1)/4)), x)


def replacement805(a, b, x):
    return Dist(a/(S(2)*b), Int(S(1)/(x**S(2)*(a + b*x**S(4))**(S(1)/4)), x), x) + Simp((a + b*x**S(4))**(S(3)/4)/(S(2)*b*x), x)


def replacement806(a, b, x):
    return -Dist(b, Int(x**S(2)/(a + b*x**S(4))**(S(5)/4), x), x) - Simp(S(1)/(x*(a + b*x**S(4))**(S(1)/4)), x)


def replacement807(a, b, x):
    return Dist(x*(a/(b*x**S(4)) + S(1))**(S(1)/4)/(a + b*x**S(4))**(S(1)/4), Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(1)/4)), x), x)


def replacement808(a, b, c, x):
    return -Dist(a/S(2), Int(sqrt(c*x)/(a + b*x**S(2))**(S(5)/4), x), x) + Simp(x*sqrt(c*x)/(a + b*x**S(2))**(S(1)/4), x)


def replacement809(a, b, c, x):
    return Dist(a*c**S(2)/(S(2)*b), Int(S(1)/((c*x)**(S(3)/2)*(a + b*x**S(2))**(S(1)/4)), x), x) + Simp(c*(a + b*x**S(2))**(S(3)/4)/(b*sqrt(c*x)), x)


def replacement810(a, b, c, x):
    return -Dist(b/c**S(2), Int(sqrt(c*x)/(a + b*x**S(2))**(S(5)/4), x), x) + Simp(-S(2)/(c*sqrt(c*x)*(a + b*x**S(2))**(S(1)/4)), x)


def replacement811(a, b, c, x):
    return Dist(sqrt(c*x)*(a/(b*x**S(2)) + S(1))**(S(1)/4)/(c**S(2)*(a + b*x**S(2))**(S(1)/4)), Int(S(1)/(x**S(2)*(a/(b*x**S(2)) + S(1))**(S(1)/4)), x), x)


def replacement812(a, b, c, m, n, p, x):
    return -Dist(a*c**n*(m - n + S(1))/(b*(m + n*p + S(1))), Int((c*x)**(m - n)*(a + b*x**n)**p, x), x) + Simp(c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*p + S(1))), x)


def replacement813(a, b, c, m, n, p, x):
    return -Dist(a*c**n*(m - n + S(1))/(b*(m + n*p + S(1))), Int((c*x)**(m - n)*(a + b*x**n)**p, x), x) + Simp(c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*p + S(1))), x)


def replacement814(a1, a2, b1, b2, c, m, n, p, x):
    return -Dist(a1*a2*c**(S(2)*n)*(m - S(2)*n + S(1))/(b1*b2*(m + S(2)*n*p + S(1))), Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x), x) + Simp(c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(b1*b2*(m + S(2)*n*p + S(1))), x)


def replacement815(a1, a2, b1, b2, c, m, n, p, x):
    return -Dist(a1*a2*c**(S(2)*n)*(m - S(2)*n + S(1))/(b1*b2*(m + S(2)*n*p + S(1))), Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x), x) + Simp(c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(b1*b2*(m + S(2)*n*p + S(1))), x)


def replacement816(a, b, c, m, n, p, x):
    return -Dist(b*c**(-n)*(m + n*(p + S(1)) + S(1))/(a*(m + S(1))), Int((c*x)**(m + n)*(a + b*x**n)**p, x), x) + Simp((c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))), x)


def replacement817(a, b, c, m, n, p, x):
    return -Dist(b*c**(-n)*(m + n*(p + S(1)) + S(1))/(a*(m + S(1))), Int((c*x)**(m + n)*(a + b*x**n)**p, x), x) + Simp((c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))), x)


def replacement818(a1, a2, b1, b2, c, m, n, p, x):
    return -Dist(b1*b2*c**(-S(2)*n)*(m + S(2)*n*(p + S(1)) + S(1))/(a1*a2*(m + S(1))), Int((c*x)**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x), x) + Simp((c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))), x)


def replacement819(a1, a2, b1, b2, c, m, n, p, x):
    return -Dist(b1*b2*c**(-S(2)*n)*(m + S(2)*n*(p + S(1)) + S(1))/(a1*a2*(m + S(1))), Int((c*x)**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x), x) + Simp((c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))), x)


def With820(a, b, c, m, n, p, x):
    k = Denominator(m)
    return Dist(k/c, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*c**(-n)*x**(k*n))**p, x), x, (c*x)**(S(1)/k)), x)


def With821(a1, a2, b1, b2, c, m, n, p, x):
    k = Denominator(m)
    return Dist(k/c, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a1 + b1*c**(-n)*x**(k*n))**p*(a2 + b2*c**(-n)*x**(k*n))**p, x), x, (c*x)**(S(1)/k)), x)


def replacement822(a, b, m, n, p, x):
    return Dist(a**(p + (m + S(1))/n), Subst(Int(x**m*(-b*x**n + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x*(a + b*x**n)**(-S(1)/n)), x)


def replacement823(a1, a2, b1, b2, m, n, p, x):
    return Dist((a1*a2)**(p + (m + S(1))/(S(2)*n)), Subst(Int(x**m*(-b1*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x*(a1 + b1*x**n)**(-S(1)/(S(2)*n))*(a2 + b2*x**n)**(-S(1)/(S(2)*n))), x)


def replacement824(a, b, m, n, p, x):
    return Dist((a/(a + b*x**n))**(p + (m + S(1))/n)*(a + b*x**n)**(p + (m + S(1))/n), Subst(Int(x**m*(-b*x**n + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x*(a + b*x**n)**(-S(1)/n)), x)


def replacement825(a1, a2, b1, b2, m, n, p, x):
    return Dist((a1/(a1 + b1*x**n))**(p + (m + S(1))/(S(2)*n))*(a2/(a2 + b2*x**n))**(p + (m + S(1))/(S(2)*n))*(a1 + b1*x**n)**(p + (m + S(1))/(S(2)*n))*(a2 + b2*x**n)**(p + (m + S(1))/(S(2)*n)), Subst(Int(x**m*(-b1*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x*(a1 + b1*x**n)**(-S(1)/(S(2)*n))*(a2 + b2*x**n)**(-S(1)/(S(2)*n))), x)


def replacement826(a, b, m, n, p, x):
    return -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p, x), x, S(1)/x)


def replacement827(a1, a2, b1, b2, m, n, p, x):
    return -Subst(Int(x**(-m + S(-2))*(a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p, x), x, S(1)/x)


def With828(a, b, c, m, n, p, x):
    k = Denominator(m)
    return -Dist(k/c, Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*c**(-n)*x**(-k*n))**p, x), x, (c*x)**(-S(1)/k)), x)


def With829(a1, a2, b1, b2, c, m, n, p, x):
    k = Denominator(m)
    return -Dist(k/c, Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a1 + b1*c**(-n)*x**(-k*n))**p*(a2 + b2*c**(-n)*x**(-k*n))**p, x), x, (c*x)**(-S(1)/k)), x)


def replacement830(a, b, c, m, n, p, x):
    return -Dist((c*x)**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p, x), x, S(1)/x), x)


def replacement831(a1, a2, b1, b2, c, m, n, p, x):
    return -Dist((c*x)**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p, x), x, S(1)/x), x)


def With832(a, b, m, n, p, x):
    k = Denominator(n)
    return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n))**p, x), x, x**(S(1)/k)), x)


def With833(a1, a2, b1, b2, m, n, p, x):
    k = Denominator(S(2)*n)
    return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a1 + b1*x**(k*n))**p*(a2 + b2*x**(k*n))**p, x), x, x**(S(1)/k)), x)


def replacement834(a, b, c, m, n, p, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a + b*x**n)**p, x), x)


def replacement835(a1, a2, b1, b2, c, m, n, p, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x), x)


def replacement836(a, b, m, n, p, x):
    return Dist(S(1)/(m + S(1)), Subst(Int((a + b*x**(n/(m + S(1))))**p, x), x, x**(m + S(1))), x)


def replacement837(a1, a2, b1, b2, m, n, p, x):
    return Dist(S(1)/(m + S(1)), Subst(Int((a1 + b1*x**(n/(m + S(1))))**p*(a2 + b2*x**(n/(m + S(1))))**p, x), x, x**(m + S(1))), x)


def replacement838(a, b, c, m, n, p, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a + b*x**n)**p, x), x)


def replacement839(a1, a2, b1, b2, c, m, n, p, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x), x)


def replacement840(a, b, m, n, p, x):
    return -Dist(b*n*p/(m + S(1)), Int(x**(m + n)*(a + b*x**n)**(p + S(-1)), x), x) + Simp(x**(m + S(1))*(a + b*x**n)**p/(m + S(1)), x)


def replacement841(a1, a2, b1, b2, m, n, p, x):
    return -Dist(S(2)*b1*b2*n*p/(m + S(1)), Int(x**(m + n)*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x), x) + Simp(x**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(m + S(1)), x)


def replacement842(a, b, c, m, n, p, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a + b*x**n)**p, x), x)


def replacement843(a1, a2, b1, b2, c, m, n, p, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x), x)


def replacement844(a, b, c, m, n, p, x):
    return Dist(a*n*p/(m + n*p + S(1)), Int((c*x)**m*(a + b*x**n)**(p + S(-1)), x), x) + Simp((c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + n*p + S(1))), x)


def replacement845(a1, a2, b1, b2, c, m, n, p, x):
    return Dist(S(2)*a1*a2*n*p/(m + S(2)*n*p + S(1)), Int((c*x)**m*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x), x) + Simp((c*x)**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(c*(m + S(2)*n*p + S(1))), x)


def With846(a, b, m, n, p, x):
    k = Denominator(p)
    return Dist(a**(p + (m + S(1))/n)*k/n, Subst(Int(x**(k*(m + S(1))/n + S(-1))*(-b*x**k + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x**(n/k)*(a + b*x**n)**(-S(1)/k)), x)


def With847(a1, a2, b1, b2, m, n, p, x):
    k = Denominator(p)
    return Dist(k*(a1*a2)**(p + (m + S(1))/(S(2)*n))/(S(2)*n), Subst(Int(x**(k*(m + S(1))/(S(2)*n) + S(-1))*(-b1*x**k + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**k + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x**(S(2)*n/k)*(a1 + b1*x**n)**(-S(1)/k)*(a2 + b2*x**n)**(-S(1)/k)), x)


def replacement848(a, b, c, m, n, p, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a + b*x**n)**p, x), x)


def replacement849(a1, a2, b1, b2, c, m, n, p, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x), x)


def replacement850(a, b, c, m, n, p, x):
    return Dist((m + n*(p + S(1)) + S(1))/(a*n*(p + S(1))), Int((c*x)**m*(a + b*x**n)**(p + S(1)), x), x) - Simp((c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))), x)


def replacement851(a1, a2, b1, b2, c, m, n, p, x):
    return Dist((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*a1*a2*n*(p + S(1))), Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x), x) - Simp((c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))), x)


def With852(a, b, m, n, x):
    mn = m - n
    return -Dist(a/b, Int(x**mn/(a + b*x**n), x), x) + Simp(x**(mn + S(1))/(b*(mn + S(1))), x)


def replacement853(a, b, m, n, x):
    return -Dist(b/a, Int(x**(m + n)/(a + b*x**n), x), x) + Simp(x**(m + S(1))/(a*(m + S(1))), x)


def replacement854(a, b, c, m, n, x):
    return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m/(a + b*x**n), x), x)


def replacement855(a, b, c, m, n, p, x):
    return Simp(a**p*(c*x)**(m + S(1))*Hypergeometric2F1(-p, (m + S(1))/n, S(1) + (m + S(1))/n, -b*x**n/a)/(c*(m + S(1))), x)


def replacement856(a, b, c, m, n, p, x):
    return Dist(a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p), Int((c*x)**m*(S(1) + b*x**n/a)**p, x), x)


def replacement857(a, b, m, n, p, v, x):
    return Dist(Coefficient(v, x, S(1))**(-m + S(-1)), Subst(Int(SimplifyIntegrand((a + b*x**n)**p*(x - Coefficient(v, x, S(0)))**m, x), x), x, v), x)


def replacement858(a, b, m, n, p, u, v, x):
    return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(a + b*x**n)**p, x), x, v), x)


def replacement859(a1, a2, b1, b2, c, m, n, p, x):
    return Dist((a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p)), Int((c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x), x)


def replacement860(a, b, c, d, n, p, q, x):
    return Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement861(a, b, c, d, n, p, q, x):
    return Int(x**(n*(p + q))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q, x)


def replacement862(a, b, c, d, n, p, q, x):
    return -Subst(Int((a + b*x**(-n))**p*(c + d*x**(-n))**q/x**S(2), x), x, S(1)/x)


def With863(a, b, c, d, n, p, q, x):
    g = Denominator(n)
    return Dist(g, Subst(Int(x**(g + S(-1))*(a + b*x**(g*n))**p*(c + d*x**(g*n))**q, x), x, x**(S(1)/g)), x)


def replacement864(a, b, c, d, n, p, x):
    return Subst(Int(S(1)/(c - x**n*(-a*d + b*c)), x), x, x*(a + b*x**n)**(-S(1)/n))


def replacement865(a, b, c, d, n, p, q, x):
    return -Dist(c*q/(a*(p + S(1))), Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1)), x), x) - Simp(x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*n*(p + S(1))), x)


def replacement866(a, b, c, d, n, p, q, x):
    return Simp(a**p*c**(-p + S(-1))*x*(c + d*x**n)**(-S(1)/n)*Hypergeometric2F1(S(1)/n, -p, S(1) + S(1)/n, -x**n*(-a*d + b*c)/(a*(c + d*x**n))), x)


def replacement867(a, b, c, d, n, p, q, x):
    return Simp(x*(c*(a + b*x**n)/(a*(c + d*x**n)))**(-p)*(a + b*x**n)**p*(c + d*x**n)**(-p - S(1)/n)*Hypergeometric2F1(S(1)/n, -p, S(1) + S(1)/n, -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/c, x)


def replacement868(a, b, c, d, n, p, q, x):
    return Simp(x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c), x)


def replacement869(a, b, c, d, n, p, q, x):
    return Dist((b*c + n*(p + S(1))*(-a*d + b*c))/(a*n*(p + S(1))*(-a*d + b*c)), Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q, x), x) - Simp(b*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*n*(p + S(1))*(-a*d + b*c)), x)


def replacement870(a, b, c, d, n, p, x):
    return Simp(c*x*(a + b*x**n)**(p + S(1))/a, x)


def replacement871(a, b, c, d, n, p, x):
    return -Dist((a*d - b*c*(n*(p + S(1)) + S(1)))/(a*b*n*(p + S(1))), Int((a + b*x**n)**(p + S(1)), x), x) - Simp(x*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*n*(p + S(1))), x)


def replacement872(a, b, c, d, n, x):
    return -Dist((-a*d + b*c)/a, Int(S(1)/(a*x**(-n) + b), x), x) + Simp(c*x/a, x)


def replacement873(a, b, c, d, n, p, x):
    return -Dist((a*d - b*c*(n*(p + S(1)) + S(1)))/(b*(n*(p + S(1)) + S(1))), Int((a + b*x**n)**p, x), x) + Simp(d*x*(a + b*x**n)**(p + S(1))/(b*(n*(p + S(1)) + S(1))), x)


def replacement874(a, b, c, d, n, p, q, x):
    return Int(PolynomialDivide((a + b*x**n)**p, (c + d*x**n)**(-q), x), x)


def replacement875(a, b, c, d, n, x):
    return Dist(b/(-a*d + b*c), Int(S(1)/(a + b*x**n), x), x) - Dist(d/(-a*d + b*c), Int(S(1)/(c + d*x**n), x), x)


def replacement876(a, b, c, d, x):
    return Dist(sqrt(S(3))/(S(2)*c), Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(-x*Rt(b/a, S(2)) + sqrt(S(3)))), x), x) + Dist(sqrt(S(3))/(S(2)*c), Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(x*Rt(b/a, S(2)) + sqrt(S(3)))), x), x)


def replacement877(a, b, c, d, x):
    return Dist(S(1)/6, Int((-x*Rt(-b/a, S(2)) + S(3))/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x), x) + Dist(S(1)/6, Int((x*Rt(-b/a, S(2)) + S(3))/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x), x)


def replacement878(a, b, c, d, x):
    return Dist(b/d, Int((a + b*x**S(2))**(S(-1)/3), x), x) - Dist((-a*d + b*c)/d, Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x), x)


def replacement879(a, b, c, d, x):
    return Dist(sqrt(-b*x**S(2)/a)/(S(2)*x), Subst(Int(S(1)/(sqrt(-b*x/a)*(a + b*x)**(S(1)/4)*(c + d*x)), x), x, x**S(2)), x)


def replacement880(a, b, c, d, x):
    return Dist(sqrt(-b*x**S(2)/a)/(S(2)*x), Subst(Int(S(1)/(sqrt(-b*x/a)*(a + b*x)**(S(3)/4)*(c + d*x)), x), x, x**S(2)), x)


def replacement881(a, b, c, d, p, x):
    return Dist(b/d, Int((a + b*x**S(2))**(p + S(-1)), x), x) - Dist((-a*d + b*c)/d, Int((a + b*x**S(2))**(p + S(-1))/(c + d*x**S(2)), x), x)


def replacement882(a, b, c, d, p, x):
    return Dist(b/(-a*d + b*c), Int((a + b*x**S(2))**p, x), x) - Dist(d/(-a*d + b*c), Int((a + b*x**S(2))**(p + S(1))/(c + d*x**S(2)), x), x)


def replacement883(a, b, c, d, x):
    return Dist(a/c, Subst(Int(S(1)/(-S(4)*a*b*x**S(4) + S(1)), x), x, x/sqrt(a + b*x**S(4))), x)


def With884(a, b, c, d, x):
    q = Rt(-a*b, S(4))
    return Simp(a*ArcTan(q*x*(a + q**S(2)*x**S(2))/(a*sqrt(a + b*x**S(4))))/(S(2)*c*q), x) + Simp(a*atanh(q*x*(a - q**S(2)*x**S(2))/(a*sqrt(a + b*x**S(4))))/(S(2)*c*q), x)


def replacement885(a, b, c, d, x):
    return Dist(b/d, Int(S(1)/sqrt(a + b*x**S(4)), x), x) - Dist((-a*d + b*c)/d, Int(S(1)/(sqrt(a + b*x**S(4))*(c + d*x**S(4))), x), x)


def replacement886(a, b, c, d, x):
    return Dist(sqrt(a/(a + b*x**S(4)))*sqrt(a + b*x**S(4)), Subst(Int(S(1)/((c - x**S(4)*(-a*d + b*c))*sqrt(-b*x**S(4) + S(1))), x), x, x/(a + b*x**S(4))**(S(1)/4)), x)


def replacement887(a, b, c, d, p, x):
    return Dist(b/d, Int((a + b*x**S(4))**(p + S(-1)), x), x) - Dist((-a*d + b*c)/d, Int((a + b*x**S(4))**(p + S(-1))/(c + d*x**S(4)), x), x)


def replacement888(a, b, c, d, x):
    return Dist(S(1)/(S(2)*c), Int(S(1)/(sqrt(a + b*x**S(4))*(-x**S(2)*Rt(-d/c, S(2)) + S(1))), x), x) + Dist(S(1)/(S(2)*c), Int(S(1)/(sqrt(a + b*x**S(4))*(x**S(2)*Rt(-d/c, S(2)) + S(1))), x), x)


def replacement889(a, b, c, d, x):
    return Dist(b/(-a*d + b*c), Int((a + b*x**S(4))**(S(-3)/4), x), x) - Dist(d/(-a*d + b*c), Int((a + b*x**S(4))**(S(1)/4)/(c + d*x**S(4)), x), x)


def replacement890(a, b, c, d, x):
    return Simp(sqrt(a + b*x**S(2))*EllipticE(ArcTan(x*Rt(d/c, S(2))), S(1) - b*c/(a*d))/(c*sqrt(c*(a + b*x**S(2))/(a*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2))), x)


def replacement891(a, b, c, d, n, p, q, x):
    return Dist(S(1)/(a*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(n*(p + S(1)) + S(1)) + d*x**n*(n*(p + q + S(1)) + S(1)), x), x), x) - Simp(x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*n*(p + S(1))), x)


def replacement892(a, b, c, d, n, p, q, x):
    return -Dist(S(1)/(a*b*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(a*d - b*c*(n*(p + S(1)) + S(1))) + d*x**n*(a*d*(n*(q + S(-1)) + S(1)) - b*c*(n*(p + q) + S(1))), x), x), x) + Simp(x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(a*d - b*c)/(a*b*n*(p + S(1))), x)


def replacement893(a, b, c, d, n, p, q, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-a*d + b*c)), Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c + b*d*x**n*(n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x), x) - Simp(b*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*n*(p + S(1))*(-a*d + b*c)), x)


def replacement894(a, b, c, d, n, p, q, x):
    return Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement895(a, b, c, d, n, p, q, x):
    return Dist(S(1)/(b*(n*(p + q) + S(1))), Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(-a*d + b*c*(n*(p + q) + S(1))) + d*x**n*(-a*d*(n*(q + S(-1)) + S(1)) + b*c*(n*(p + S(2)*q + S(-1)) + S(1))), x), x), x) + Simp(d*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*(n*(p + q) + S(1))), x)


def replacement896(a, b, c, d, n, p, q, x):
    return Dist(n/(n*(p + q) + S(1)), Int((a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x), x) + Simp(x*(a + b*x**n)**p*(c + d*x**n)**q/(n*(p + q) + S(1)), x)


def replacement897(a, b, c, d, x):
    return Simp(sqrt(a + b*x**S(2))*EllipticF(ArcTan(x*Rt(d/c, S(2))), S(1) - b*c/(a*d))/(a*sqrt(c*(a + b*x**S(2))/(a*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2))), x)


def replacement898(a, b, c, d, x):
    return Simp(EllipticF(asin(x*Rt(-d/c, S(2))), b*c/(a*d))/(sqrt(a)*sqrt(c)*Rt(-d/c, S(2))), x)


def replacement899(a, b, c, d, x):
    return -Simp(EllipticF(acos(x*Rt(-d/c, S(2))), b*c/(-a*d + b*c))/(sqrt(c)*sqrt(a - b*c/d)*Rt(-d/c, S(2))), x)


def replacement900(a, b, c, d, x):
    return Dist(sqrt(S(1) + d*x**S(2)/c)/sqrt(c + d*x**S(2)), Int(S(1)/(sqrt(S(1) + d*x**S(2)/c)*sqrt(a + b*x**S(2))), x), x)


def replacement901(a, b, c, d, x):
    return Dist(a, Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x), x) + Dist(b, Int(x**S(2)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x), x)


def replacement902(a, b, c, d, x):
    return Dist(b/d, Int(sqrt(c + d*x**S(2))/sqrt(a + b*x**S(2)), x), x) - Dist((-a*d + b*c)/d, Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x), x)


def replacement903(a, b, c, d, x):
    return Simp(sqrt(a)*EllipticE(asin(x*Rt(-d/c, S(2))), b*c/(a*d))/(sqrt(c)*Rt(-d/c, S(2))), x)


def replacement904(a, b, c, d, x):
    return -Simp(sqrt(a - b*c/d)*EllipticE(acos(x*Rt(-d/c, S(2))), b*c/(-a*d + b*c))/(sqrt(c)*Rt(-d/c, S(2))), x)


def replacement905(a, b, c, d, x):
    return Dist(sqrt(a + b*x**S(2))/sqrt(S(1) + b*x**S(2)/a), Int(sqrt(S(1) + b*x**S(2)/a)/sqrt(c + d*x**S(2)), x), x)


def replacement906(a, b, c, d, x):
    return Dist(sqrt(S(1) + d*x**S(2)/c)/sqrt(c + d*x**S(2)), Int(sqrt(a + b*x**S(2))/sqrt(S(1) + d*x**S(2)/c), x), x)


def replacement907(a, b, c, d, n, p, q, x):
    return Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement908(a, b, c, d, n, p, q, x):
    return Simp(a**p*c**q*x*AppellF1(S(1)/n, -p, -q, S(1) + S(1)/n, -b*x**n/a, -d*x**n/c), x)


def replacement909(a, b, c, d, n, p, q, x):
    return Dist(a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p), Int((S(1) + b*x**n/a)**p*(c + d*x**n)**q, x), x)


def replacement910(a, b, c, d, mn, n, p, q, x):
    return Int(x**(-n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x)


def replacement911(a, b, c, d, mn, n, p, q, x):
    return Dist(x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q)), Int(x**(-n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x), x)


def replacement912(a, b, c, d, n, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*x**n)**p*(c + d*x**n)**q, x), x, u), x)


def replacement913(p, q, u, v, x):
    return Int(NormalizePseudoBinomial(u, x)**p*NormalizePseudoBinomial(v, x)**q, x)


def replacement914(m, p, q, u, v, x):
    return Int(NormalizePseudoBinomial(v, x)**q*NormalizePseudoBinomial(u*x**(m/p), x)**p, x)


def replacement915(b, c, d, e, m, n, p, q, x):
    return Dist(b**(S(1) - (m + S(1))/n)*e**m/n, Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n)*(c + d*x)**q, x), x, x**n), x)


def replacement916(b, c, d, e, m, n, p, q, x):
    return Dist(b**IntPart(p)*e**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p), Int(x**(m + n*p)*(c + d*x**n)**q, x), x)


def replacement917(b, c, d, e, m, n, p, q, x):
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement918(a, b, c, d, m, n, p, q, x):
    return Dist(S(1)/n, Subst(Int((a + b*x)**p*(c + d*x)**q, x), x, x**n), x)


def replacement919(a, b, c, d, m, n, p, q, x):
    return Int(x**(m + n*(p + q))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q, x)


def replacement920(a, b, c, d, m, n, p, q, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p*(c + d*x)**q, x), x, x**n), x)


def replacement921(a, b, c, d, e, m, n, p, q, x):
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement922(a, b, c, d, e, m, n, p, q, x):
    return Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement923(a, b, c, d, e, m, n, p, x):
    return Simp(c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))), x)


def replacement924(a1, a2, b1, b2, c, d, e, m, n, non2, p, x):
    return Simp(c*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(a1*a2*e*(m + S(1))), x)


def replacement925(a, b, c, d, e, m, n, p, x):
    return Dist(d*e**(-n), Int((e*x)**(m + n)*(a + b*x**n)**p, x), x) + Simp(c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))), x)


def replacement926(a, b, c, d, e, m, n, p, x):
    return Dist(d/b, Int((e*x)**m*(a + b*x**n)**(p + S(1)), x), x) + Simp((e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*e*(m + S(1))), x)


def replacement927(a, b, c, d, e, m, n, p, x):
    return Dist(e**(-n)*(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))/(a*(m + S(1))), Int((e*x)**(m + n)*(a + b*x**n)**p, x), x) + Simp(c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))), x)


def replacement928(a1, a2, b1, b2, c, d, e, m, n, non2, p, x):
    return Dist(e**(-n)*(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))/(a1*a2*(m + S(1))), Int((e*x)**(m + n)*(a1 + b1*x**(n/S(2)))**p*(a2 + b2*x**(n/S(2)))**p, x), x) + Simp(c*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(a1*a2*e*(m + S(1))), x)


def replacement929(a, b, c, d, m, p, x):
    return Dist(b**(-m/S(2) + S(-1))/(S(2)*(p + S(1))), Int((a + b*x**S(2))**(p + S(1))*ExpandToSum(S(2)*b*x**S(2)*(p + S(1))*Together((b**(m/S(2))*x**(m + S(-2))*(c + d*x**S(2)) - (-a)**(m/S(2) + S(-1))*(-a*d + b*c))/(a + b*x**S(2))) - (-a)**(m/S(2) + S(-1))*(-a*d + b*c), x), x), x) + Simp(b**(-m/S(2) + S(-1))*x*(-a)**(m/S(2) + S(-1))*(a + b*x**S(2))**(p + S(1))*(-a*d + b*c)/(S(2)*(p + S(1))), x)


def replacement930(a, b, c, d, m, p, x):
    return Dist(b**(-m/S(2) + S(-1))/(S(2)*(p + S(1))), Int(x**m*(a + b*x**S(2))**(p + S(1))*ExpandToSum(S(2)*b*(p + S(1))*Together((b**(m/S(2))*(c + d*x**S(2)) - x**(S(2) - m)*(-a)**(m/S(2) + S(-1))*(-a*d + b*c))/(a + b*x**S(2))) - x**(-m)*(-a)**(m/S(2) + S(-1))*(-a*d + b*c), x), x), x) + Simp(b**(-m/S(2) + S(-1))*x*(-a)**(m/S(2) + S(-1))*(a + b*x**S(2))**(p + S(1))*(-a*d + b*c)/(S(2)*(p + S(1))), x)


def replacement931(a, b, c, d, e, m, n, p, x):
    return -Dist((a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))/(a*b*n*(p + S(1))), Int((e*x)**m*(a + b*x**n)**(p + S(1)), x), x) - Simp((e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*e*n*(p + S(1))), x)


def replacement932(a1, a2, b1, b2, c, d, e, m, n, non2, p, x):
    return -Dist((a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))/(a1*a2*b1*b2*n*(p + S(1))), Int((e*x)**m*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1)), x), x) - Simp((e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))*(-a1*a2*d + b1*b2*c)/(a1*a2*b1*b2*e*n*(p + S(1))), x)


def replacement933(a, b, c, d, e, m, n, p, x):
    return -Dist((a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))/(b*(m + n*(p + S(1)) + S(1))), Int((e*x)**m*(a + b*x**n)**p, x), x) + Simp(d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(b*e*(m + n*(p + S(1)) + S(1))), x)


def replacement934(a1, a2, b1, b2, c, d, e, m, n, non2, p, x):
    return -Dist((a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))/(b1*b2*(m + n*(p + S(1)) + S(1))), Int((e*x)**m*(a1 + b1*x**(n/S(2)))**p*(a2 + b2*x**(n/S(2)))**p, x), x) + Simp(d*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(b1*b2*e*(m + n*(p + S(1)) + S(1))), x)


def replacement935(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p/(c + d*x**n), x), x)


def replacement936(a, b, c, d, e, m, n, p, x):
    return -Dist(e**(-n)/(a*(m + S(1))), Int((e*x)**(m + n)*(a + b*x**n)**p*Simp(-a*d**S(2)*x**n*(m + S(1)) + b*c**S(2)*n*(p + S(1)) + c*(m + S(1))*(-S(2)*a*d + b*c), x), x), x) + Simp(c**S(2)*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))), x)


def replacement937(a, b, c, d, e, m, n, p, x):
    return Dist(S(1)/(a*b**S(2)*n*(p + S(1))), Int((e*x)**m*(a + b*x**n)**(p + S(1))*Simp(a*b*d**S(2)*n*x**n*(p + S(1)) + b**S(2)*c**S(2)*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)**S(2), x), x), x) - Simp((e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)**S(2)/(a*b**S(2)*e*n*(p + S(1))), x)


def replacement938(a, b, c, d, e, m, n, p, x):
    return Dist(S(1)/(b*(m + n*(p + S(2)) + S(1))), Int((e*x)**m*(a + b*x**n)**p*Simp(b*c**S(2)*(m + n*(p + S(2)) + S(1)) + d*x**n*(S(2)*b*c*n*(p + S(1)) + (-a*d + S(2)*b*c)*(m + n + S(1))), x), x), x) + Simp(d**S(2)*e**(-n + S(-1))*(e*x)**(m + n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*(p + S(2)) + S(1))), x)


def With939(a, b, c, d, m, n, p, q, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    k = GCD(m + S(1), n)
    if Unequal(k, S(1)):
        return True
    return False


def replacement939(a, b, c, d, m, n, p, q, x):

    k = GCD(m + S(1), n)
    return Dist(S(1)/k, Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p*(c + d*x**(n/k))**q, x), x, x**k), x)


def With940(a, b, c, d, e, m, n, p, q, x):
    k = Denominator(m)
    return Dist(k/e, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*e**(-n)*x**(k*n))**p*(c + d*e**(-n)*x**(k*n))**q, x), x, (e*x)**(S(1)/k)), x)


def replacement941(a, b, c, d, e, m, n, p, q, x):
    return -Dist(e**n/(b*n*(p + S(1))), Int((e*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m - n + S(1)) + d*x**n*(m + n*(q + S(-1)) + S(1)), x), x), x) + Simp(e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*n*(p + S(1))), x)


def replacement942(a, b, c, d, e, m, n, p, q, x):
    return Dist(S(1)/(a*b*n*(p + S(1))), Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)) + d*x**n*(b*c*n*(p + S(1)) + (-a*d + b*c)*(m + n*(q + S(-1)) + S(1))), x), x), x) - Simp((e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(-a*d + b*c)/(a*b*e*n*(p + S(1))), x)


def replacement943(a, b, c, d, e, m, n, p, q, x):
    return Dist(S(1)/(a*n*(p + S(1))), Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m + n*(p + S(1)) + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x), x) - Simp((e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*n*(p + S(1))), x)


def replacement944(a, b, c, d, e, m, n, p, q, x):
    return Dist(e**(S(2)*n)/(b*n*(p + S(1))*(-a*d + b*c)), Int((e*x)**(m - S(2)*n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(a*c*(m - S(2)*n + S(1)) + x**n*(a*d*(m + n*q - n + S(1)) + b*c*n*(p + S(1))), x), x), x) - Simp(a*e**(S(2)*n + S(-1))*(e*x)**(m - S(2)*n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*n*(p + S(1))*(-a*d + b*c)), x)


def replacement945(a, b, c, d, e, m, n, p, q, x):
    return -Dist(e**n/(n*(p + S(1))*(-a*d + b*c)), Int((e*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m - n + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x), x) + Simp(e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(n*(p + S(1))*(-a*d + b*c)), x)


def replacement946(a, b, c, d, e, m, n, p, q, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-a*d + b*c)), Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c*(m + S(1)) + b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x), x) - Simp(b*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*e*n*(p + S(1))*(-a*d + b*c)), x)


def replacement947(a, b, c, d, e, m, n, p, q, x):
    return -Dist(e**(-n)*n/(m + S(1)), Int((e*x)**(m + n)*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*d*q + b*c*p + b*d*x**n*(p + q), x), x), x) + Simp((e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + S(1))), x)


def replacement948(a, b, c, d, e, m, n, p, q, x):
    return -Dist(e**(-n)/(a*(m + S(1))), Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*n*(a*d*(q + S(-1)) + b*c*(p + S(1))) + c*(m + S(1))*(-a*d + b*c) + d*x**n*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)), x), x), x) + Simp(c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(a*e*(m + S(1))), x)


def replacement949(a, b, c, d, e, m, n, p, q, x):
    return -Dist(e**(-n)/(a*(m + S(1))), Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(b*c*(m + S(1)) + d*x**n*(b*n*(p + q + S(1)) + b*(m + S(1))) + n*(a*d*q + b*c*(p + S(1))), x), x), x) + Simp((e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*(m + S(1))), x)


def replacement950(a, b, c, d, e, m, n, p, q, x):
    return Dist(n/(m + n*(p + q) + S(1)), Int((e*x)**m*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x), x) + Simp((e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + n*(p + q) + S(1))), x)


def replacement951(a, b, c, d, e, m, n, p, q, x):
    return Dist(S(1)/(b*(m + n*(p + q) + S(1))), Int((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)) + x**n*(b*c*d*n*(p + q) + d*n*(q + S(-1))*(-a*d + b*c) + d*(m + S(1))*(-a*d + b*c)), x), x), x) + Simp(d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*e*(m + n*(p + q) + S(1))), x)


def replacement952(a, b, c, d, e, m, n, p, q, x):
    return -Dist(e**n/(b*(m + n*(p + q) + S(1))), Int((e*x)**(m - n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(a*c*(m - n + S(1)) + x**n*(a*d*(m - n + S(1)) - n*q*(-a*d + b*c)), x), x), x) + Simp(e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*(m + n*(p + q) + S(1))), x)


def replacement953(a, b, c, d, e, m, n, p, q, x):
    return -Dist(e**(S(2)*n)/(b*d*(m + n*(p + q) + S(1))), Int((e*x)**(m - S(2)*n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*(m - S(2)*n + S(1)) + x**n*(a*d*(m + n*(q + S(-1)) + S(1)) + b*c*(m + n*(p + S(-1)) + S(1))), x), x), x) + Simp(e**(S(2)*n + S(-1))*(e*x)**(m - S(2)*n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*d*(m + n*(p + q) + S(1))), x)


def replacement954(a, b, c, d, e, m, n, p, q, x):
    return -Dist(e**(-n)/(a*c*(m + S(1))), Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(a*d*q + b*c*p) + (a*d + b*c)*(m + n + S(1)), x), x), x) + Simp((e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c*e*(m + S(1))), x)


def replacement955(a, b, c, d, e, m, n, x):
    return -Dist(a*e**n/(-a*d + b*c), Int((e*x)**(m - n)/(a + b*x**n), x), x) + Dist(c*e**n/(-a*d + b*c), Int((e*x)**(m - n)/(c + d*x**n), x), x)


def replacement956(a, b, c, d, e, m, n, x):
    return Dist(b/(-a*d + b*c), Int((e*x)**m/(a + b*x**n), x), x) - Dist(d/(-a*d + b*c), Int((e*x)**m/(c + d*x**n), x), x)


def replacement957(a, b, c, d, m, n, x):
    return Dist(S(1)/b, Int(x**(m - n)/sqrt(c + d*x**n), x), x) - Dist(a/b, Int(x**(m - n)/((a + b*x**n)*sqrt(c + d*x**n)), x), x)


def With958(a, b, c, d, x):
    r = Numerator(Rt(-a/b, S(2)))
    s = Denominator(Rt(-a/b, S(2)))
    return -Dist(s/(S(2)*b), Int(S(1)/(sqrt(c + d*x**S(4))*(r - s*x**S(2))), x), x) + Dist(s/(S(2)*b), Int(S(1)/(sqrt(c + d*x**S(4))*(r + s*x**S(2))), x), x)


def With959(a, b, c, d, x):
    q = Rt(d/c, S(3))
    return Simp(S(2)**(S(1)/3)*q*log(-S(2)**(S(1)/3)*q*x + S(1) - sqrt(c + d*x**S(3))/sqrt(c))/(S(12)*b*sqrt(c)), x) - Simp(S(2)**(S(1)/3)*q*log(-S(2)**(S(1)/3)*q*x + S(1) + sqrt(c + d*x**S(3))/sqrt(c))/(S(12)*b*sqrt(c)), x) + Simp(S(2)**(S(1)/3)*q*atanh(sqrt(c + d*x**S(3))/sqrt(c))/(S(18)*b*sqrt(c)), x) - Simp(S(2)**(S(1)/3)*sqrt(S(3))*q*ArcTan(sqrt(S(3))/S(3) + S(2)**(S(2)/3)*sqrt(S(3))*(sqrt(c) - sqrt(c + d*x**S(3)))/(S(3)*sqrt(c)*q*x))/(S(18)*b*sqrt(c)), x) + Simp(S(2)**(S(1)/3)*sqrt(S(3))*q*ArcTan(sqrt(S(3))/S(3) + S(2)**(S(2)/3)*sqrt(S(3))*(sqrt(c) + sqrt(c + d*x**S(3)))/(S(3)*sqrt(c)*q*x))/(S(18)*b*sqrt(c)), x)


def replacement960(a, b, c, d, m, x):
    return Dist(S(1)/b, Int(x**(m + S(-3))/sqrt(c + d*x**S(3)), x), x) - Dist(a/b, Int(x**(m + S(-3))/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x), x)


def replacement961(a, b, c, d, m, x):
    return Dist(S(1)/a, Int(x**m/sqrt(c + d*x**S(3)), x), x) - Dist(b/a, Int(x**(m + S(3))/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x), x)


def replacement962(a, b, c, d, x):
    return Dist(d/b, Int(x**S(2)/sqrt(c + d*x**S(4)), x), x) + Dist((-a*d + b*c)/b, Int(x**S(2)/((a + b*x**S(4))*sqrt(c + d*x**S(4))), x), x)


def replacement963(a, b, c, d, m, x):
    return Dist(d/b, Int(x**m/sqrt(c + d*x**S(3)), x), x) + Dist((-a*d + b*c)/b, Int(x**m/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x), x)


def replacement964(a, b, c, d, x):
    return -Dist(c/b, Int(sqrt(a + b*x**S(2))/(c + d*x**S(2))**(S(3)/2), x), x) + Simp(x*sqrt(a + b*x**S(2))/(b*sqrt(c + d*x**S(2))), x)


def replacement965(a, b, c, d, n, x):
    return Dist(S(1)/b, Int(sqrt(a + b*x**n)/sqrt(c + d*x**n), x), x) - Dist(a/b, Int(S(1)/(sqrt(a + b*x**n)*sqrt(c + d*x**n)), x), x)


def With966(a, b, c, d, m, n, p, q, x):
    k = Denominator(p)
    return Dist(a**(p + (m + S(1))/n)*k/n, Subst(Int(x**(k*(m + S(1))/n + S(-1))*(c - x**k*(-a*d + b*c))**q*(-b*x**k + S(1))**(-p - q + S(-1) - (m + S(1))/n), x), x, x**(n/k)*(a + b*x**n)**(-S(1)/k)), x)


def replacement967(a, b, c, d, m, n, p, q, x):
    return -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q, x), x, S(1)/x)


def With968(a, b, c, d, e, m, n, p, q, x):
    g = Denominator(m)
    return -Dist(g/e, Subst(Int(x**(-g*(m + S(1)) + S(-1))*(a + b*e**(-n)*x**(-g*n))**p*(c + d*e**(-n)*x**(-g*n))**q, x), x, (e*x)**(-S(1)/g)), x)


def replacement969(a, b, c, d, e, m, n, p, q, x):
    return -Dist((e*x)**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q, x), x, S(1)/x), x)


def With970(a, b, c, d, m, n, p, q, x):
    g = Denominator(n)
    return Dist(g, Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + b*x**(g*n))**p*(c + d*x**(g*n))**q, x), x, x**(S(1)/g)), x)


def replacement971(a, b, c, d, e, m, n, p, q, x):
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement972(a, b, c, d, m, n, p, q, x):
    return Dist(S(1)/(m + S(1)), Subst(Int((a + b*x**(n/(m + S(1))))**p*(c + d*x**(n/(m + S(1))))**q, x), x, x**(m + S(1))), x)


def replacement973(a, b, c, d, e, m, n, p, q, x):
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement974(a, b, c, d, e, m, n, p, q, x):
    return Dist(S(1)/(a*b*n*(p + S(1))), Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)) + d*x**n*(b*c*n*(p + S(1)) + (-a*d + b*c)*(m + n*(q + S(-1)) + S(1))), x), x), x) - Simp((e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(-a*d + b*c)/(a*b*e*n*(p + S(1))), x)


def replacement975(a, b, c, d, e, m, n, p, q, x):
    return Dist(S(1)/(a*n*(p + S(1))), Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m + n*(p + S(1)) + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x), x) - Simp((e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*n*(p + S(1))), x)


def replacement976(a, b, c, d, e, m, n, p, q, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-a*d + b*c)), Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c*(m + S(1)) + b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x), x) - Simp(b*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*e*n*(p + S(1))*(-a*d + b*c)), x)


def replacement977(a, b, c, d, e, m, n, p, q, x):
    return Dist(n/(m + n*(p + q) + S(1)), Int((e*x)**m*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x), x) + Simp((e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + n*(p + q) + S(1))), x)


def replacement978(a, b, c, d, e, m, n, p, q, x):
    return Dist(S(1)/(b*(m + n*(p + q) + S(1))), Int((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)) + x**n*(b*c*d*n*(p + q) + d*n*(q + S(-1))*(-a*d + b*c) + d*(m + S(1))*(-a*d + b*c)), x), x), x) + Simp(d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*e*(m + n*(p + q) + S(1))), x)


def replacement979(a, b, c, d, m, n, x):
    return -Dist(a/(-a*d + b*c), Int(x**(m - n)/(a + b*x**n), x), x) + Dist(c/(-a*d + b*c), Int(x**(m - n)/(c + d*x**n), x), x)


def replacement980(a, b, c, d, e, m, n, x):
    return Dist(b/(-a*d + b*c), Int((e*x)**m/(a + b*x**n), x), x) - Dist(d/(-a*d + b*c), Int((e*x)**m/(c + d*x**n), x), x)


def replacement981(a, b, c, d, e, m, n, p, q, x):
    return Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement982(a, b, c, d, m, mn, n, p, q, x):
    return Int(x**(m - n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x)


def replacement983(a, b, c, d, m, mn, n, p, q, x):
    return Dist(x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q)), Int(x**(m - n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x), x)


def replacement984(a, b, c, d, e, m, mn, n, p, q, x):
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*x**n)**p*(c + d*x**(-n))**q, x), x)


def replacement985(a, b, c, d, e, m, n, p, q, x):
    return Simp(a**p*c**q*(e*x)**(m + S(1))*AppellF1((m + S(1))/n, -p, -q, S(1) + (m + S(1))/n, -b*x**n/a, -d*x**n/c)/(e*(m + S(1))), x)


def replacement986(a, b, c, d, e, m, n, p, q, x):
    return Dist(a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p), Int((e*x)**m*(S(1) + b*x**n/a)**p*(c + d*x**n)**q, x), x)


def replacement987(a, b, c, d, m, n, p, q, v, x):
    return Dist(Coefficient(v, x, S(1))**(-m + S(-1)), Subst(Int(SimplifyIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(x - Coefficient(v, x, S(0)))**m, x), x), x, v), x)


def replacement988(a, b, c, d, m, n, p, q, u, v, x):
    return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x, v), x)


def replacement989(a1, a2, b1, b2, c, d, n, non2, p, q, u, x):
    return Int(u*(c + d*x**n)**q*(a1*a2 + b1*b2*x**n)**p, x)


def replacement990(a1, a2, b1, b2, c, d, e, n, n2, non2, p, q, u, x):
    return Int(u*(a1*a2 + b1*b2*x**n)**p*(c + d*x**n + e*x**(S(2)*n))**q, x)


def replacement991(a1, a2, b1, b2, c, d, n, non2, p, q, u, x):
    return Dist((a1 + b1*x**(n/S(2)))**FracPart(p)*(a2 + b2*x**(n/S(2)))**FracPart(p)*(a1*a2 + b1*b2*x**n)**(-FracPart(p)), Int(u*(c + d*x**n)**q*(a1*a2 + b1*b2*x**n)**p, x), x)


def replacement992(a1, a2, b1, b2, c, d, e, n, n2, non2, p, q, u, x):
    return Dist((a1 + b1*x**(n/S(2)))**FracPart(p)*(a2 + b2*x**(n/S(2)))**FracPart(p)*(a1*a2 + b1*b2*x**n)**(-FracPart(p)), Int(u*(a1*a2 + b1*b2*x**n)**p*(c + d*x**n + e*x**(S(2)*n))**q, x), x)


def replacement993(a, b, c, d, e, f, n, p, q, r, x):
    return Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x)


def replacement994(a, b, c, d, e, f, n, x):
    return Dist((-a*f + b*e)/(-a*d + b*c), Int(S(1)/(a + b*x**n), x), x) - Dist((-c*f + d*e)/(-a*d + b*c), Int(S(1)/(c + d*x**n), x), x)


def replacement995(a, b, c, d, e, f, n, x):
    return Dist(f/b, Int(S(1)/sqrt(c + d*x**n), x), x) + Dist((-a*f + b*e)/b, Int(S(1)/((a + b*x**n)*sqrt(c + d*x**n)), x), x)


def replacement996(a, b, c, d, e, f, n, x):
    return Dist(f/b, Int(sqrt(a + b*x**n)/sqrt(c + d*x**n), x), x) + Dist((-a*f + b*e)/b, Int(S(1)/(sqrt(a + b*x**n)*sqrt(c + d*x**n)), x), x)


def replacement997(a, b, c, d, e, f, x):
    return Dist((-a*f + b*e)/(-a*d + b*c), Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x), x) - Dist((-c*f + d*e)/(-a*d + b*c), Int(sqrt(a + b*x**S(2))/(c + d*x**S(2))**(S(3)/2), x), x)


def replacement998(a, b, c, d, e, f, n, p, q, x):
    return Dist(S(1)/(a*b*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(-a*f + b*e*n*(p + S(1)) + b*e) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(n*q + S(1))), x), x), x) - Simp(x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*n*(p + S(1))), x)


def replacement999(a, b, c, d, e, f, n, p, q, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-a*d + b*c)), Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x), x) - Simp(x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*n*(p + S(1))*(-a*d + b*c)), x)


def replacement1000(a, b, c, d, e, f, n, p, q, x):
    return Dist(S(1)/(b*(n*(p + q + S(1)) + S(1))), Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(-a*f + b*e*n*(p + q + S(1)) + b*e) + x**n*(b*d*e*n*(p + q + S(1)) + d*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x), x) + Simp(f*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*(n*(p + q + S(1)) + S(1))), x)


def replacement1001(a, b, c, d, e, f, x):
    return Dist((-a*f + b*e)/(-a*d + b*c), Int((a + b*x**S(4))**(S(-3)/4), x), x) - Dist((-c*f + d*e)/(-a*d + b*c), Int((a + b*x**S(4))**(S(1)/4)/(c + d*x**S(4)), x), x)


def replacement1002(a, b, c, d, e, f, n, p, x):
    return Dist(f/d, Int((a + b*x**n)**p, x), x) + Dist((-c*f + d*e)/d, Int((a + b*x**n)**p/(c + d*x**n), x), x)


def replacement1003(a, b, c, d, e, f, n, p, q, x):
    return Dist(e, Int((a + b*x**n)**p*(c + d*x**n)**q, x), x) + Dist(f, Int(x**n*(a + b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement1004(a, b, c, d, e, f, x):
    return Dist(b/(-a*d + b*c), Int(S(1)/((a + b*x**S(2))*sqrt(e + f*x**S(2))), x), x) - Dist(d/(-a*d + b*c), Int(S(1)/((c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x)


def replacement1005(c, d, e, f, x):
    return Dist(S(1)/c, Int(S(1)/(x**S(2)*sqrt(e + f*x**S(2))), x), x) - Dist(d/c, Int(S(1)/((c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x)


def replacement1006(a, b, c, d, e, f, x):
    return Dist(d/b, Int(sqrt(e + f*x**S(2))/sqrt(c + d*x**S(2)), x), x) + Dist((-a*d + b*c)/b, Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x), x)


def replacement1007(a, b, c, d, e, f, x):
    return Dist(d/b, Int(sqrt(e + f*x**S(2))/sqrt(c + d*x**S(2)), x), x) + Dist((-a*d + b*c)/b, Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x), x)


def replacement1008(a, b, c, d, e, f, x):
    return Dist(b/(-a*f + b*e), Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x), x) - Dist(f/(-a*f + b*e), Int(S(1)/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x)


def replacement1009(a, b, c, d, e, f, x):
    return Simp(EllipticPi(b*c/(a*d), asin(x*Rt(-d/c, S(2))), c*f/(d*e))/(a*sqrt(c)*sqrt(e)*Rt(-d/c, S(2))), x)


def replacement1010(a, b, c, d, e, f, x):
    return Dist(sqrt(S(1) + d*x**S(2)/c)/sqrt(c + d*x**S(2)), Int(S(1)/(sqrt(S(1) + d*x**S(2)/c)*(a + b*x**S(2))*sqrt(e + f*x**S(2))), x), x)


def replacement1011(a, b, c, d, e, f, x):
    return Simp(c*sqrt(e + f*x**S(2))*EllipticPi(S(1) - b*c/(a*d), ArcTan(x*Rt(d/c, S(2))), -c*f/(d*e) + S(1))/(a*e*sqrt(c*(e + f*x**S(2))/(e*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2))), x)


def replacement1012(a, b, c, d, e, f, x):
    return Dist(d/b, Int(S(1)/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x) + Dist((-a*d + b*c)/b, Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x)


def replacement1013(a, b, c, d, e, f, x):
    return Dist(b/(-a*d + b*c), Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x), x) - Dist(d/(-a*d + b*c), Int(sqrt(e + f*x**S(2))/(c + d*x**S(2))**(S(3)/2), x), x)


def replacement1014(a, b, c, d, e, f, x):
    return Dist((-a*f + b*e)/(-a*d + b*c), Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x), x) - Dist((-c*f + d*e)/(-a*d + b*c), Int(sqrt(e + f*x**S(2))/(c + d*x**S(2))**(S(3)/2), x), x)


def replacement1015(a, b, c, d, e, f, x):
    return Dist(d/b**S(2), Int(sqrt(e + f*x**S(2))*(-a*d + S(2)*b*c + b*d*x**S(2))/sqrt(c + d*x**S(2)), x), x) + Dist((-a*d + b*c)**S(2)/b**S(2), Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x), x)


def replacement1016(a, b, c, d, e, f, q, r, x):
    return Dist(b*(-a*f + b*e)/(-a*d + b*c)**S(2), Int((c + d*x**S(2))**(q + S(2))*(e + f*x**S(2))**(r + S(-1))/(a + b*x**S(2)), x), x) - Dist((-a*d + b*c)**(S(-2)), Int((c + d*x**S(2))**q*(e + f*x**S(2))**(r + S(-1))*(-a*d**S(2)*e - b*c**S(2)*f + S(2)*b*c*d*e + d**S(2)*x**S(2)*(-a*f + b*e)), x), x)


def replacement1017(a, b, c, d, e, f, q, r, x):
    return Dist(d/b, Int((c + d*x**S(2))**(q + S(-1))*(e + f*x**S(2))**r, x), x) + Dist((-a*d + b*c)/b, Int((c + d*x**S(2))**(q + S(-1))*(e + f*x**S(2))**r/(a + b*x**S(2)), x), x)


def replacement1018(a, b, c, d, e, f, q, r, x):
    return Dist(b**S(2)/(-a*d + b*c)**S(2), Int((c + d*x**S(2))**(q + S(2))*(e + f*x**S(2))**r/(a + b*x**S(2)), x), x) - Dist(d/(-a*d + b*c)**S(2), Int((c + d*x**S(2))**q*(e + f*x**S(2))**r*(-a*d + S(2)*b*c + b*d*x**S(2)), x), x)


def replacement1019(a, b, c, d, e, f, q, r, x):
    return Dist(b/(-a*d + b*c), Int((c + d*x**S(2))**(q + S(1))*(e + f*x**S(2))**r/(a + b*x**S(2)), x), x) - Dist(d/(-a*d + b*c), Int((c + d*x**S(2))**q*(e + f*x**S(2))**r, x), x)


def replacement1020(a, b, c, d, e, f, x):
    return Dist((-a**S(2)*d*f + b**S(2)*c*e)/(S(2)*a*b**S(2)), Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x) + Dist(d*f/(S(2)*a*b**S(2)), Int((a - b*x**S(2))/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x) + Simp(x*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*a*(a + b*x**S(2))), x)


def replacement1021(a, b, c, d, e, f, x):
    return Dist((S(3)*a**S(2)*d*f - S(2)*a*b*(c*f + d*e) + b**S(2)*c*e)/(S(2)*a*(-a*d + b*c)*(-a*f + b*e)), Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x) - Dist(d*f/(S(2)*a*(-a*d + b*c)*(-a*f + b*e)), Int((a + b*x**S(2))/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x) + Simp(b**S(2)*x*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*a*(a + b*x**S(2))*(-a*d + b*c)*(-a*f + b*e)), x)


def replacement1022(a, b, c, d, e, f, n, p, q, r, x):
    return Dist(d/b, Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(e + f*x**n)**r, x), x) + Dist((-a*d + b*c)/b, Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*(e + f*x**n)**r, x), x)


def replacement1023(a, b, c, d, e, f, n, p, q, r, x):
    return Dist(b/(-a*d + b*c), Int((a + b*x**n)**p*(c + d*x**n)**(q + S(1))*(e + f*x**n)**r, x), x) - Dist(d/(-a*d + b*c), Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(e + f*x**n)**r, x), x)


def replacement1024(a, b, c, d, e, f, x):
    return Dist(sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))/(c*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2))), Subst(Int(S(1)/(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)*sqrt(S(1) - x**S(2)*(-a*f + b*e)/e)), x), x, x/sqrt(a + b*x**S(2))), x)


def replacement1025(a, b, c, d, e, f, x):
    return Dist(a*sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))/(c*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2))), Subst(Int(S(1)/(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)*sqrt(S(1) - x**S(2)*(-a*f + b*e)/e)*(-b*x**S(2) + S(1))), x), x, x/sqrt(a + b*x**S(2))), x)


def replacement1026(a, b, c, d, e, f, x):
    return Dist(sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))/(a*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2))), Subst(Int(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)/sqrt(S(1) - x**S(2)*(-a*f + b*e)/e), x), x, x/sqrt(a + b*x**S(2))), x)


def replacement1027(a, b, c, d, e, f, x):
    return -Dist(c*(-c*f + d*e)/(S(2)*f), Int(sqrt(a + b*x**S(2))/((c + d*x**S(2))**(S(3)/2)*sqrt(e + f*x**S(2))), x), x) - Dist((-a*d*f - b*c*f + b*d*e)/(S(2)*d*f), Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))), x), x) + Dist(b*c*(-c*f + d*e)/(S(2)*d*f), Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x) + Simp(d*x*sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*f*sqrt(c + d*x**S(2))), x)


def replacement1028(a, b, c, d, e, f, x):
    return -Dist((-a*d*f - b*c*f + b*d*e)/(S(2)*f**S(2)), Int(sqrt(e + f*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x), x) + Dist(e*(-a*f + b*e)/(S(2)*f), Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*(e + f*x**S(2))**(S(3)/2)), x), x) + Dist((-a*f + b*e)*(-S(2)*c*f + d*e)/(S(2)*f**S(2)), Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x) + Simp(x*sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))/(S(2)*sqrt(e + f*x**S(2))), x)


def replacement1029(a, b, c, d, e, f, x):
    return Dist(b/f, Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))), x), x) - Dist((-a*f + b*e)/f, Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*(e + f*x**S(2))**(S(3)/2)), x), x)


def With1030(a, b, c, d, e, f, n, p, q, r, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)
    if SumQ(u):
        return True
    return False


def replacement1030(a, b, c, d, e, f, n, p, q, r, x):

    u = ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)
    return Int(u, x)


def replacement1031(a, b, c, d, e, f, n, p, q, r, x):
    return -Subst(Int((a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r/x**S(2), x), x, S(1)/x)


def replacement1032(a, b, c, d, e, f, n, p, q, r, x):
    return Int((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)


def replacement1033(a, b, c, d, e, f, n, p, q, r, u, v, w, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x, u), x)


def replacement1034(a, b, c, d, e, f, mn, n, p, q, r, x):
    return Int(x**(-n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x)


def replacement1035(a, b, c, d, e, f, mn, n, p, q, r, x):
    return Int(x**(n*(p + r))*(c + d*x**(-n))**q*(a*x**(-n) + b)**p*(e*x**(-n) + f)**r, x)


def replacement1036(a, b, c, d, e, f, mn, n, p, q, r, x):
    return Dist(x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q)), Int(x**(-n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x), x)


def replacement1037(a, b, c, d, e1, e2, f1, f2, n, n2, p, q, r, x):
    return Int((a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x)


def replacement1038(a, b, c, d, e1, e2, f1, f2, n, n2, p, q, r, x):
    return Dist((e1 + f1*x**(n/S(2)))**FracPart(r)*(e2 + f2*x**(n/S(2)))**FracPart(r)*(e1*e2 + f1*f2*x**n)**(-FracPart(r)), Int((a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x), x)


def replacement1039(b, c, d, e, f, g, m, n, p, q, r, x):
    return Dist(b**(S(1) - (m + S(1))/n)*g**m/n, Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n)*(c + d*x)**q*(e + f*x)**r, x), x, x**n), x)


def replacement1040(b, c, d, e, f, g, m, n, p, q, r, x):
    return Dist(b**IntPart(p)*g**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p), Int(x**(m + n*p)*(c + d*x**n)**q*(e + f*x**n)**r, x), x)


def replacement1041(b, c, d, e, f, g, m, n, p, q, r, x):
    return Dist(g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m), Int(x**m*(b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x)


def replacement1042(a, b, c, d, e, f, g, m, n, p, q, r, x):
    return Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x)


def replacement1043(a, b, c, d, e, f, m, n, p, q, r, x):
    return Dist(S(1)/n, Subst(Int((a + b*x)**p*(c + d*x)**q*(e + f*x)**r, x), x, x**n), x)


def replacement1044(a, b, c, d, e, f, m, n, p, q, r, x):
    return Int(x**(m + n*(p + q + r))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q*(e*x**(-n) + f)**r, x)


def replacement1045(a, b, c, d, e, f, m, n, p, q, r, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p*(c + d*x)**q*(e + f*x)**r, x), x, x**n), x)


def replacement1046(a, b, c, d, e, f, g, m, n, p, q, r, x):
    return Dist(g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m), Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x)


def With1047(a, b, c, d, e, f, m, n, p, q, r, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    k = GCD(m + S(1), n)
    if Unequal(k, S(1)):
        return True
    return False


def replacement1047(a, b, c, d, e, f, m, n, p, q, r, x):

    k = GCD(m + S(1), n)
    return Dist(S(1)/k, Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p*(c + d*x**(n/k))**q*(e + f*x**(n/k))**r, x), x, x**k), x)


def With1048(a, b, c, d, e, f, g, m, n, p, q, r, x):
    k = Denominator(m)
    return Dist(k/g, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*g**(-n)*x**(k*n))**p*(c + d*g**(-n)*x**(k*n))**q*(e + f*g**(-n)*x**(k*n))**r, x), x, (g*x)**(S(1)/k)), x)


def replacement1049(a, b, c, d, e, f, g, m, n, p, q, x):
    return Dist(S(1)/(a*b*n*(p + S(1))), Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + S(1)) + (m + S(1))*(-a*f + b*e)) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(m + n*q + S(1))), x), x), x) - Simp((g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*g*n*(p + S(1))), x)


def replacement1050(a, b, c, d, e, f, g, m, n, p, q, x):
    return -Dist(g**n/(b*n*(p + S(1))*(-a*d + b*c)), Int((g*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(-a*f + b*e)*(m - n + S(1)) + x**n*(-b*n*(p + S(1))*(c*f - d*e) + d*(-a*f + b*e)*(m + n*q + S(1))), x), x), x) + Simp(g**(n + S(-1))*(g*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(b*n*(p + S(1))*(-a*d + b*c)), x)


def replacement1051(a, b, c, d, e, f, g, m, n, p, q, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-a*d + b*c)), Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(m + n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x), x) - Simp((g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*g*n*(p + S(1))*(-a*d + b*c)), x)


def replacement1052(a, b, c, d, e, f, g, m, n, p, q, x):
    return -Dist(g**(-n)/(a*(m + S(1))), Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + e*n*(a*d*q + b*c*(p + S(1))), x), x), x) + Simp(e*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*g*(m + S(1))), x)


def replacement1053(a, b, c, d, e, f, g, m, n, p, q, x):
    return Dist(S(1)/(b*(m + n*(p + q + S(1)) + S(1))), Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + x**n*(b*d*e*n*(p + q + S(1)) + d*(m + S(1))*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x), x) + Simp(f*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*g*(m + n*(p + q + S(1)) + S(1))), x)


def replacement1054(a, b, c, d, e, f, g, m, n, p, q, x):
    return -Dist(g**n/(b*d*(m + n*(p + q + S(1)) + S(1))), Int((g*x)**(m - n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*f*(m - n + S(1)) + x**n*(a*d*f*(m + n*q + S(1)) + b*(c*f*(m + n*p + S(1)) - d*e*(m + n*(p + q + S(1)) + S(1)))), x), x), x) + Simp(f*g**(n + S(-1))*(g*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*d*(m + n*(p + q + S(1)) + S(1))), x)


def replacement1055(a, b, c, d, e, f, g, m, n, p, q, x):
    return Dist(g**(-n)/(a*c*(m + S(1))), Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*f*(m + S(1)) - b*d*e*x**n*(m + n*(p + q + S(2)) + S(1)) - e*n*(a*d*q + b*c*p) - e*(a*d + b*c)*(m + n + S(1)), x), x), x) + Simp(e*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c*g*(m + S(1))), x)


def replacement1056(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(e + f*x**n)/(c + d*x**n), x), x)


def replacement1057(a, b, c, d, e, f, g, m, n, p, q, x):
    return Dist(e, Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x) + Dist(e**(-n)*f, Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement1058(a, b, c, d, e, f, g, m, n, p, q, r, x):
    return Dist(e, Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**(r + S(-1)), x), x) + Dist(e**(-n)*f, Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**(r + S(-1)), x), x)


def replacement1059(a, b, c, d, e, f, m, n, p, q, r, x):
    return -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r, x), x, S(1)/x)


def With1060(a, b, c, d, e, f, g, m, n, p, q, r, x):
    k = Denominator(m)
    return -Dist(k/g, Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*g**(-n)*x**(-k*n))**p*(c + d*g**(-n)*x**(-k*n))**q*(e + f*g**(-n)*x**(-k*n))**r, x), x, (g*x)**(-S(1)/k)), x)


def replacement1061(a, b, c, d, e, f, g, m, n, p, q, r, x):
    return -Dist((g*x)**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r, x), x, S(1)/x), x)


def With1062(a, b, c, d, e, f, m, n, p, q, r, x):
    k = Denominator(n)
    return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n))**p*(c + d*x**(k*n))**q*(e + f*x**(k*n))**r, x), x, x**(S(1)/k)), x)


def replacement1063(a, b, c, d, e, f, g, m, n, p, q, r, x):
    return Dist(g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m), Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x)


def replacement1064(a, b, c, d, e, f, m, n, p, q, r, x):
    return Dist(S(1)/(m + S(1)), Subst(Int((a + b*x**(n/(m + S(1))))**p*(c + d*x**(n/(m + S(1))))**q*(e + f*x**(n/(m + S(1))))**r, x), x, x**(m + S(1))), x)


def replacement1065(a, b, c, d, e, f, g, m, n, p, q, r, x):
    return Dist(g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m), Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x)


def replacement1066(a, b, c, d, e, f, g, m, n, p, q, x):
    return Dist(S(1)/(a*b*n*(p + S(1))), Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + S(1)) + (m + S(1))*(-a*f + b*e)) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(m + n*q + S(1))), x), x), x) - Simp((g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*g*n*(p + S(1))), x)


def replacement1067(a, b, c, d, e, f, g, m, n, p, q, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-a*d + b*c)), Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(m + n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x), x) - Simp((g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*g*n*(p + S(1))*(-a*d + b*c)), x)


def replacement1068(a, b, c, d, e, f, g, m, n, p, q, x):
    return Dist(S(1)/(b*(m + n*(p + q + S(1)) + S(1))), Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + x**n*(b*d*e*n*(p + q + S(1)) + d*(m + S(1))*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x), x) + Simp(f*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*g*(m + n*(p + q + S(1)) + S(1))), x)


def replacement1069(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(e + f*x**n)/(c + d*x**n), x), x)


def replacement1070(a, b, c, d, e, f, g, m, n, p, q, x):
    return Dist(e, Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x) + Dist(f*x**(-m)*(g*x)**m, Int(x**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q, x), x)


def replacement1071(a, b, c, d, e, f, m, mn, n, p, q, r, x):
    return Int(x**(m - n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x)


def replacement1072(a, b, c, d, e, f, m, mn, n, p, q, r, x):
    return Int(x**(m + n*(p + r))*(c + d*x**(-n))**q*(a*x**(-n) + b)**p*(e*x**(-n) + f)**r, x)


def replacement1073(a, b, c, d, e, f, m, mn, n, p, q, r, x):
    return Dist(x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q)), Int(x**(m - n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x), x)


def replacement1074(a, b, c, d, e, f, g, m, mn, n, p, q, r, x):
    return Dist(g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m), Int(x**m*(a + b*x**n)**p*(c + d*x**(-n))**q*(e + f*x**n)**r, x), x)


def replacement1075(a, b, c, d, e, f, g, m, n, p, q, r, x):
    return Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)


def replacement1076(a, b, c, d, e, f, m, n, p, q, r, u, v, x):
    return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x, v), x)


def replacement1077(a, b, c, d, e1, e2, f1, f2, g, m, n, n2, p, q, r, x):
    return Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x)


def replacement1078(a, b, c, d, e1, e2, f1, f2, g, m, n, n2, p, q, r, x):
    return Dist((e1 + f1*x**(n/S(2)))**FracPart(r)*(e2 + f2*x**(n/S(2)))**FracPart(r)*(e1*e2 + f1*f2*x**n)**(-FracPart(r)), Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x), x)
