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

def binomial_products(rubi):
    from sympy.integrals.rubi.constraints import cons459, cons3, cons4, cons5, cons460, cons2, cons461, cons54, cons462, cons87, cons463, cons38, cons464, cons148, cons13, cons163, cons465, cons466, cons43, cons448, cons67, cons137, cons467, cons468, cons469, cons470, cons471, cons472, cons473, cons474, cons475, cons476, cons477, cons478, cons479, cons480, cons481, cons482, cons483, cons484, cons105, cons485, cons486, cons487, cons488, cons196, cons489, cons128, cons357, cons490, cons491, cons492, cons493, cons68, cons69, cons55, cons494, cons57, cons58, cons59, cons60, cons495, cons496, cons497, cons498, cons147, cons7, cons21, cons499, cons500, cons501, cons18, cons502, cons503, cons66, cons504, cons505, cons506, cons507, cons17, cons244, cons94, cons508, cons509, cons510, cons511, cons512, cons513, cons514, cons515, cons516, cons517, cons518, cons519, cons520, cons521, cons62, cons522, cons523, cons524, cons525, cons526, cons527, cons528, cons529, cons31, cons530, cons531, cons532, cons533, cons534, cons535, cons536, cons367, cons537, cons538, cons539, cons540, cons356, cons541, cons23, cons542, cons543, cons544, cons545, cons546, cons547, cons548, cons549, cons550, cons551, cons552, cons553, cons554, cons71, cons555, cons27, cons220, cons50, cons556, cons85, cons557, cons395, cons403, cons63, cons558, cons559, cons560, cons561, cons562, cons563, cons564, cons565, cons566, cons567, cons568, cons569, cons70, cons570, cons571, cons572, cons573, cons402, cons574, cons575, cons576, cons405, cons577, cons578, cons579, cons580, cons581, cons177, cons582, cons583, cons117, cons584, cons585, cons586, cons587, cons386, cons588, cons589, cons590, cons591, cons48, cons53, cons592, cons593, cons594, cons595, cons596, cons93, cons597, cons598, cons599, cons600, cons601, cons602, cons603, cons604, cons88, cons605, cons606, cons607, cons608, cons609, cons610, cons611, cons612, cons613, cons614, cons615, cons616, cons617, cons618, cons619, cons620, cons621, cons622, cons623, cons624, cons625, cons626, cons627, cons46, cons628, cons125, cons629, cons630, cons631, cons153, cons632, cons633, cons176, cons634, cons635, cons636, cons637, cons638, cons178, cons639, cons640, cons396, cons641, cons52, cons642, cons643, cons644, cons645, cons646, cons647, cons648, cons649, cons650, cons651, cons652, cons653, cons654, cons655, cons656, cons208, cons657, cons658, cons659, cons660, cons661, cons380, cons662, cons663

    pattern689 = Pattern(Integral((x_**n_*WC('b', S(1)))**p_, x_), cons3, cons4, cons5, cons459)
    def replacement689(n, x, p, b):
        return b**IntPart(p)*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(n*p), x)
    rule689 = ReplacementRule(pattern689, replacement689)
    rubi.add(rule689.pattern, label = 689)

    pattern690 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons460)
    def replacement690(n, p, x, b, a):
        return x*(a + b*x**n)**(p + S(1))/a
    rule690 = ReplacementRule(pattern690, replacement690)
    rubi.add(rule690.pattern, label = 690)

    pattern691 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons461, cons54)
    def replacement691(n, p, x, b, a):
        return -x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))) + (n*(p + S(1)) + S(1))*Int((a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1)))
    rule691 = ReplacementRule(pattern691, replacement691)
    rubi.add(rule691.pattern, label = 691)

    pattern692 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons4, cons462)
    def replacement692(n, x, b, a):
        return Int(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n), x)
    rule692 = ReplacementRule(pattern692, replacement692)
    rubi.add(rule692.pattern, label = 692)

    pattern693 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons87, cons463, cons38)
    def replacement693(n, p, x, b, a):
        return Int(x**(n*p)*(a*x**(-n) + b)**p, x)
    rule693 = ReplacementRule(pattern693, replacement693)
    rubi.add(rule693.pattern, label = 693)

    pattern694 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons464)
    def replacement694(n, p, x, b, a):
        return Int(ExpandIntegrand((a + b*x**n)**p, x), x)
    rule694 = ReplacementRule(pattern694, replacement694)
    rubi.add(rule694.pattern, label = 694)

    pattern695 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons148, cons13, cons163, cons465)
    def replacement695(n, p, x, b, a):
        return a*n*p*Int((a + b*x**n)**(p + S(-1)), x)/(n*p + S(1)) + x*(a + b*x**n)**p/(n*p + S(1))
    rule695 = ReplacementRule(pattern695, replacement695)
    rubi.add(rule695.pattern, label = 695)

    pattern696 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-5)/4), x_), cons2, cons3, cons466, cons43)
    def replacement696(x, b, a):
        return S(2)*EllipticE(ArcTan(x*Rt(b/a, S(2)))/S(2), S(2))/(a**(S(5)/4)*Rt(b/a, S(2)))
    rule696 = ReplacementRule(pattern696, replacement696)
    rubi.add(rule696.pattern, label = 696)

    pattern697 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-5)/4), x_), cons2, cons3, cons466, cons448)
    def replacement697(x, b, a):
        return (S(1) + b*x**S(2)/a)**(S(1)/4)*Int((S(1) + b*x**S(2)/a)**(S(-5)/4), x)/(a*(a + b*x**S(2))**(S(1)/4))
    rule697 = ReplacementRule(pattern697, replacement697)
    rubi.add(rule697.pattern, label = 697)

    pattern698 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-7)/6), x_), cons2, cons3, cons67)
    def replacement698(x, b, a):
        return Subst(Int((-b*x**S(2) + S(1))**(S(-1)/3), x), x, x/sqrt(a + b*x**S(2)))/((a/(a + b*x**S(2)))**(S(2)/3)*(a + b*x**S(2))**(S(2)/3))
    rule698 = ReplacementRule(pattern698, replacement698)
    rubi.add(rule698.pattern, label = 698)

    pattern699 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons148, cons13, cons137, cons465)
    def replacement699(n, p, x, b, a):
        return -x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))) + (n*(p + S(1)) + S(1))*Int((a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1)))
    rule699 = ReplacementRule(pattern699, replacement699)
    rubi.add(rule699.pattern, label = 699)

    pattern700 = Pattern(Integral(S(1)/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons67)
    def replacement700(x, b, a):
        return Int((-x*Rt(b, S(3)) + S(2)*Rt(a, S(3)))/(x**S(2)*Rt(b, S(3))**S(2) - x*Rt(a, S(3))*Rt(b, S(3)) + Rt(a, S(3))**S(2)), x)/(S(3)*Rt(a, S(3))**S(2)) + Int(S(1)/(x*Rt(b, S(3)) + Rt(a, S(3))), x)/(S(3)*Rt(a, S(3))**S(2))
    rule700 = ReplacementRule(pattern700, replacement700)
    rubi.add(rule700.pattern, label = 700)

    def With701(n, x, b, a):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r - s*x*cos(Pi*(2*k - 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r*Int(S(1)/(r + s*x), x)/(a*n)
    pattern701 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons467, cons468 )
    rule701 = ReplacementRule(pattern701, With701)
    rubi.add(rule701.pattern, label = 701)

    def With702(n, x, b, a):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r + s*x*cos(Pi*(2*k - 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r*Int(S(1)/(r - s*x), x)/(a*n)
    pattern702 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons467, cons469 )
    rule702 = ReplacementRule(pattern702, With702)
    rubi.add(rule702.pattern, label = 702)

    pattern703 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons468, cons470)
    def replacement703(x, b, a):
        return ArcTan(x*Rt(b, S(2))/Rt(a, S(2)))/(Rt(a, S(2))*Rt(b, S(2)))
    rule703 = ReplacementRule(pattern703, replacement703)
    rubi.add(rule703.pattern, label = 703)

    pattern704 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons468, cons471)
    def replacement704(x, b, a):
        return -ArcTan(x*Rt(-b, S(2))/Rt(-a, S(2)))/(Rt(-a, S(2))*Rt(-b, S(2)))
    rule704 = ReplacementRule(pattern704, replacement704)
    rubi.add(rule704.pattern, label = 704)

    pattern705 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons468)
    def replacement705(x, b, a):
        return ArcTan(x/Rt(a/b, S(2)))*Rt(a/b, S(2))/a
    rule705 = ReplacementRule(pattern705, replacement705)
    rubi.add(rule705.pattern, label = 705)

    pattern706 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons469, cons472)
    def replacement706(x, b, a):
        return atanh(x*Rt(-b, S(2))/Rt(a, S(2)))/(Rt(a, S(2))*Rt(-b, S(2)))
    rule706 = ReplacementRule(pattern706, replacement706)
    rubi.add(rule706.pattern, label = 706)

    pattern707 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons469, cons473)
    def replacement707(x, b, a):
        return -atanh(x*Rt(b, S(2))/Rt(-a, S(2)))/(Rt(-a, S(2))*Rt(b, S(2)))
    rule707 = ReplacementRule(pattern707, replacement707)
    rubi.add(rule707.pattern, label = 707)

    pattern708 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons469)
    def replacement708(x, b, a):
        return Rt(-a/b, S(2))*atanh(x/Rt(-a/b, S(2)))/a
    rule708 = ReplacementRule(pattern708, replacement708)
    rubi.add(rule708.pattern, label = 708)

    def With709(n, x, b, a):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        v = Symbol('v')
        u = Int((r - s*x*cos(Pi*(2*k - 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x) + Int((r + s*x*cos(Pi*(2*k - 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**S(2)*Int(S(1)/(r**S(2) + s**S(2)*x**S(2)), x)/(a*n)
    pattern709 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons474, cons468 )
    rule709 = ReplacementRule(pattern709, With709)
    rubi.add(rule709.pattern, label = 709)

    def With710(n, x, b, a):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r - s*x*cos(2*Pi*k/n))/(r**2 - 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x) + Int((r + s*x*cos(2*Pi*k/n))/(r**2 + 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**S(2)*Int(S(1)/(r**S(2) - s**S(2)*x**S(2)), x)/(a*n)
    pattern710 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons474, cons469 )
    rule710 = ReplacementRule(pattern710, With710)
    rubi.add(rule710.pattern, label = 710)

    def With711(x, b, a):
        r = Numerator(Rt(a/b, S(2)))
        s = Denominator(Rt(a/b, S(2)))
        return Int((r - s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*r) + Int((r + s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*r)
    pattern711 = Pattern(Integral(S(1)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons475 )
    rule711 = ReplacementRule(pattern711, With711)
    rubi.add(rule711.pattern, label = 711)

    def With712(x, b, a):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return r*Int(S(1)/(r - s*x**S(2)), x)/(S(2)*a) + r*Int(S(1)/(r + s*x**S(2)), x)/(S(2)*a)
    pattern712 = Pattern(Integral(S(1)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons476 )
    rule712 = ReplacementRule(pattern712, With712)
    rubi.add(rule712.pattern, label = 712)

    def With713(n, x, b, a):
        r = Numerator(Rt(a/b, S(4)))
        s = Denominator(Rt(a/b, S(4)))
        return sqrt(S(2))*r*Int((sqrt(S(2))*r - s*x**(n/S(4)))/(r**S(2) - sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*a) + sqrt(S(2))*r*Int((sqrt(S(2))*r + s*x**(n/S(4)))/(r**S(2) + sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*a)
    pattern713 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons477, cons478 )
    rule713 = ReplacementRule(pattern713, With713)
    rubi.add(rule713.pattern, label = 713)

    def With714(n, x, b, a):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return r*Int(S(1)/(r - s*x**(n/S(2))), x)/(S(2)*a) + r*Int(S(1)/(r + s*x**(n/S(2))), x)/(S(2)*a)
    pattern714 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons477, cons476 )
    rule714 = ReplacementRule(pattern714, With714)
    rubi.add(rule714.pattern, label = 714)

    pattern715 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons43, cons479)
    def replacement715(x, b, a):
        return asinh(x*Rt(b, S(2))/sqrt(a))/Rt(b, S(2))
    rule715 = ReplacementRule(pattern715, replacement715)
    rubi.add(rule715.pattern, label = 715)

    pattern716 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons43, cons480)
    def replacement716(x, b, a):
        return asin(x*Rt(-b, S(2))/sqrt(a))/Rt(-b, S(2))
    rule716 = ReplacementRule(pattern716, replacement716)
    rubi.add(rule716.pattern, label = 716)

    pattern717 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons448)
    def replacement717(x, b, a):
        return Subst(Int(S(1)/(-b*x**S(2) + S(1)), x), x, x/sqrt(a + b*x**S(2)))
    rule717 = ReplacementRule(pattern717, replacement717)
    rubi.add(rule717.pattern, label = 717)

    def With718(x, b, a):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return S(2)*S(3)**(S(3)/4)*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(sqrt(S(3)) + S(2))*(r*x + s)*EllipticF(asin((r*x + s*(-sqrt(S(3)) + S(1)))/(r*x + s*(S(1) + sqrt(S(3))))), S(-7) - S(4)*sqrt(S(3)))/(S(3)*r*sqrt(s*(r*x + s)/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(a + b*x**S(3)))
    pattern718 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons481 )
    rule718 = ReplacementRule(pattern718, With718)
    rubi.add(rule718.pattern, label = 718)

    def With719(x, b, a):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return S(2)*S(3)**(S(3)/4)*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(-sqrt(S(3)) + S(1)))**S(2))*sqrt(-sqrt(S(3)) + S(2))*(r*x + s)*EllipticF(asin((r*x + s*(S(1) + sqrt(S(3))))/(r*x + s*(-sqrt(S(3)) + S(1)))), S(-7) + S(4)*sqrt(S(3)))/(S(3)*r*sqrt(-s*(r*x + s)/(r*x + s*(-sqrt(S(3)) + S(1)))**S(2))*sqrt(a + b*x**S(3)))
    pattern719 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons482 )
    rule719 = ReplacementRule(pattern719, With719)
    rubi.add(rule719.pattern, label = 719)

    def With720(x, b, a):
        q = Rt(b/a, S(4))
        return sqrt((a + b*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticF(S(2)*ArcTan(q*x), S(1)/2)/(S(2)*q*sqrt(a + b*x**S(4)))
    pattern720 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons466 )
    rule720 = ReplacementRule(pattern720, With720)
    rubi.add(rule720.pattern, label = 720)

    pattern721 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons483, cons43)
    def replacement721(x, b, a):
        return EllipticF(asin(x*Rt(-b, S(4))/Rt(a, S(4))), S(-1))/(Rt(a, S(4))*Rt(-b, S(4)))
    rule721 = ReplacementRule(pattern721, replacement721)
    rubi.add(rule721.pattern, label = 721)

    def With722(x, b, a):
        q = Rt(-a*b, S(2))
        if IntegerQ(q):
            return sqrt(S(2))*sqrt((a + q*x**S(2))/q)*sqrt(-a + q*x**S(2))*EllipticF(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(S(2)*sqrt(-a)*sqrt(a + b*x**S(4)))
        return False
    pattern722 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons484, cons105, CustomConstraint(With722))
    rule722 = ReplacementRule(pattern722, With722)
    rubi.add(rule722.pattern, label = 722)

    def With723(x, b, a):
        q = Rt(-a*b, S(2))
        return sqrt(S(2))*sqrt((a + q*x**S(2))/q)*sqrt((a - q*x**S(2))/(a + q*x**S(2)))*EllipticF(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(S(2)*sqrt(a/(a + q*x**S(2)))*sqrt(a + b*x**S(4)))
    pattern723 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons484, cons105 )
    rule723 = ReplacementRule(pattern723, With723)
    rubi.add(rule723.pattern, label = 723)

    pattern724 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons483, cons448)
    def replacement724(x, b, a):
        return sqrt(S(1) + b*x**S(4)/a)*Int(S(1)/sqrt(S(1) + b*x**S(4)/a), x)/sqrt(a + b*x**S(4))
    rule724 = ReplacementRule(pattern724, replacement724)
    rubi.add(rule724.pattern, label = 724)

    def With725(x, b, a):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return S(3)**(S(3)/4)*x*sqrt((r**S(2)*x**S(4) - r*s*x**S(2) + s**S(2))/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*(r*x**S(2) + s)*EllipticF(acos((r*x**S(2)*(-sqrt(S(3)) + S(1)) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)), sqrt(S(3))/S(4) + S(1)/2)/(S(6)*s*sqrt(r*x**S(2)*(r*x**S(2) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*sqrt(a + b*x**S(6)))
    pattern725 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), cons2, cons3, cons67 )
    rule725 = ReplacementRule(pattern725, With725)
    rubi.add(rule725.pattern, label = 725)

    pattern726 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), cons2, cons3, cons67)
    def replacement726(x, b, a):
        return Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/S(2) + Int((x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/S(2)
    rule726 = ReplacementRule(pattern726, replacement726)
    rubi.add(rule726.pattern, label = 726)

    pattern727 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), cons2, cons3, cons466)
    def replacement727(x, b, a):
        return -a*Int((a + b*x**S(2))**(S(-5)/4), x) + S(2)*x/(a + b*x**S(2))**(S(1)/4)
    rule727 = ReplacementRule(pattern727, replacement727)
    rubi.add(rule727.pattern, label = 727)

    pattern728 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), cons2, cons3, cons483, cons43)
    def replacement728(x, b, a):
        return S(2)*EllipticE(asin(x*Rt(-b/a, S(2)))/S(2), S(2))/(a**(S(1)/4)*Rt(-b/a, S(2)))
    rule728 = ReplacementRule(pattern728, replacement728)
    rubi.add(rule728.pattern, label = 728)

    pattern729 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), cons2, cons3, cons483, cons448)
    def replacement729(x, b, a):
        return (S(1) + b*x**S(2)/a)**(S(1)/4)*Int((S(1) + b*x**S(2)/a)**(S(-1)/4), x)/(a + b*x**S(2))**(S(1)/4)
    rule729 = ReplacementRule(pattern729, replacement729)
    rubi.add(rule729.pattern, label = 729)

    pattern730 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), cons2, cons3, cons43, cons466)
    def replacement730(x, b, a):
        return S(2)*EllipticF(ArcTan(x*Rt(b/a, S(2)))/S(2), S(2))/(a**(S(3)/4)*Rt(b/a, S(2)))
    rule730 = ReplacementRule(pattern730, replacement730)
    rubi.add(rule730.pattern, label = 730)

    pattern731 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), cons2, cons3, cons43, cons483)
    def replacement731(x, b, a):
        return S(2)*EllipticF(asin(x*Rt(-b/a, S(2)))/S(2), S(2))/(a**(S(3)/4)*Rt(-b/a, S(2)))
    rule731 = ReplacementRule(pattern731, replacement731)
    rubi.add(rule731.pattern, label = 731)

    pattern732 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), cons2, cons3, cons448)
    def replacement732(x, b, a):
        return (S(1) + b*x**S(2)/a)**(S(3)/4)*Int((S(1) + b*x**S(2)/a)**(S(-3)/4), x)/(a + b*x**S(2))**(S(3)/4)
    rule732 = ReplacementRule(pattern732, replacement732)
    rubi.add(rule732.pattern, label = 732)

    pattern733 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/3), x_), cons2, cons3, cons67)
    def replacement733(x, b, a):
        return S(3)*sqrt(b*x**S(2))*Subst(Int(x/sqrt(-a + x**S(3)), x), x, (a + b*x**S(2))**(S(1)/3))/(S(2)*b*x)
    rule733 = ReplacementRule(pattern733, replacement733)
    rubi.add(rule733.pattern, label = 733)

    pattern734 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-2)/3), x_), cons2, cons3, cons67)
    def replacement734(x, b, a):
        return S(3)*sqrt(b*x**S(2))*Subst(Int(S(1)/sqrt(-a + x**S(3)), x), x, (a + b*x**S(2))**(S(1)/3))/(S(2)*b*x)
    rule734 = ReplacementRule(pattern734, replacement734)
    rubi.add(rule734.pattern, label = 734)

    pattern735 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**(S(-3)/4), x_), cons2, cons3, cons67)
    def replacement735(x, b, a):
        return x**S(3)*(a/(b*x**S(4)) + S(1))**(S(3)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(3)/4)), x)/(a + b*x**S(4))**(S(3)/4)
    rule735 = ReplacementRule(pattern735, replacement735)
    rubi.add(rule735.pattern, label = 735)

    pattern736 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/6), x_), cons2, cons3, cons67)
    def replacement736(x, b, a):
        return -a*Int((a + b*x**S(2))**(S(-7)/6), x)/S(2) + S(3)*x/(S(2)*(a + b*x**S(2))**(S(1)/6))
    rule736 = ReplacementRule(pattern736, replacement736)
    rubi.add(rule736.pattern, label = 736)

    pattern737 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons148, cons13, cons485, cons486, cons487)
    def replacement737(n, p, x, b, a):
        return a**(p + S(1)/n)*Subst(Int((-b*x**n + S(1))**(-p + S(-1) - S(1)/n), x), x, x*(a + b*x**n)**(-S(1)/n))
    rule737 = ReplacementRule(pattern737, replacement737)
    rubi.add(rule737.pattern, label = 737)

    pattern738 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons148, cons13, cons485, cons486, cons488)
    def replacement738(n, p, x, b, a):
        return (a/(a + b*x**n))**(p + S(1)/n)*(a + b*x**n)**(p + S(1)/n)*Subst(Int((-b*x**n + S(1))**(-p + S(-1) - S(1)/n), x), x, x*(a + b*x**n)**(-S(1)/n))
    rule738 = ReplacementRule(pattern738, replacement738)
    rubi.add(rule738.pattern, label = 738)

    pattern739 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons196)
    def replacement739(n, p, x, b, a):
        return -Subst(Int((a + b*x**(-n))**p/x**S(2), x), x, S(1)/x)
    rule739 = ReplacementRule(pattern739, replacement739)
    rubi.add(rule739.pattern, label = 739)

    def With740(n, p, x, b, a):
        k = Denominator(n)
        return k*Subst(Int(x**(k + S(-1))*(a + b*x**(k*n))**p, x), x, x**(S(1)/k))
    pattern740 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons489 )
    rule740 = ReplacementRule(pattern740, With740)
    rubi.add(rule740.pattern, label = 740)

    pattern741 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons128)
    def replacement741(n, p, x, b, a):
        return Int(ExpandIntegrand((a + b*x**n)**p, x), x)
    rule741 = ReplacementRule(pattern741, replacement741)
    rubi.add(rule741.pattern, label = 741)

    pattern742 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons357, cons490, cons491, cons492)
    def replacement742(n, p, x, b, a):
        return a**p*x*Hypergeometric2F1(-p, S(1)/n, S(1) + S(1)/n, -b*x**n/a)
    rule742 = ReplacementRule(pattern742, replacement742)
    rubi.add(rule742.pattern, label = 742)

    pattern743 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons357, cons490, cons491, cons493)
    def replacement743(n, p, x, b, a):
        return a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((S(1) + b*x**n/a)**p, x)
    rule743 = ReplacementRule(pattern743, replacement743)
    rubi.add(rule743.pattern, label = 743)

    pattern744 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons5, cons68, cons69)
    def replacement744(n, p, x, b, u, a):
        return Subst(Int((a + b*x**n)**p, x), x, u)/Coefficient(u, x, S(1))
    rule744 = ReplacementRule(pattern744, replacement744)
    rubi.add(rule744.pattern, label = 744)

    pattern745 = Pattern(Integral((x_**n_*WC('b1', S(1)) + WC('a1', S(0)))**WC('p', S(1))*(x_**n_*WC('b2', S(1)) + WC('a2', S(0)))**WC('p', S(1)), x_), cons57, cons58, cons59, cons60, cons4, cons5, cons55, cons494)
    def replacement745(n, a1, p, b2, a2, b1, x):
        return Int((a1*a2 + b1*b2*x**(S(2)*n))**p, x)
    rule745 = ReplacementRule(pattern745, replacement745)
    rubi.add(rule745.pattern, label = 745)

    pattern746 = Pattern(Integral((a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), cons57, cons58, cons59, cons60, cons55, cons495, cons13, cons163, cons496)
    def replacement746(n, a1, p, b2, a2, b1, x):
        return S(2)*a1*a2*n*p*Int((a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(S(2)*n*p + S(1)) + x*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(S(2)*n*p + S(1))
    rule746 = ReplacementRule(pattern746, replacement746)
    rubi.add(rule746.pattern, label = 746)

    pattern747 = Pattern(Integral((a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons55, cons495, cons13, cons137, cons496)
    def replacement747(n, a1, p, b2, a2, b1, x):
        return -x*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*n*(p + S(1))) + (S(2)*n*(p + S(1)) + S(1))*Int((a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1)))
    rule747 = ReplacementRule(pattern747, replacement747)
    rubi.add(rule747.pattern, label = 747)

    pattern748 = Pattern(Integral((a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons5, cons55, cons497)
    def replacement748(n, a1, p, b2, a2, x, b1):
        return -Subst(Int((a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p/x**S(2), x), x, S(1)/x)
    rule748 = ReplacementRule(pattern748, replacement748)
    rubi.add(rule748.pattern, label = 748)

    def With749(n, a1, p, b2, a2, x, b1):
        k = Denominator(S(2)*n)
        return k*Subst(Int(x**(k + S(-1))*(a1 + b1*x**(k*n))**p*(a2 + b2*x**(k*n))**p, x), x, x**(S(1)/k))
    pattern749 = Pattern(Integral((a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons5, cons55, cons498 )
    rule749 = ReplacementRule(pattern749, With749)
    rubi.add(rule749.pattern, label = 749)

    pattern750 = Pattern(Integral((x_**n_*WC('b1', S(1)) + WC('a1', S(0)))**p_*(x_**n_*WC('b2', S(1)) + WC('a2', S(0)))**p_, x_), cons57, cons58, cons59, cons60, cons4, cons5, cons55, cons147)
    def replacement750(n, a1, p, b2, a2, b1, x):
        return (a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p))*Int((a1*a2 + b1*b2*x**(S(2)*n))**p, x)
    rule750 = ReplacementRule(pattern750, replacement750)
    rubi.add(rule750.pattern, label = 750)

    pattern751 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons5, cons55, cons494)
    def replacement751(m, n, a1, p, b2, a2, x, b1, c):
        return Int((c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x)
    rule751 = ReplacementRule(pattern751, replacement751)
    rubi.add(rule751.pattern, label = 751)

    pattern752 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_, x_), cons3, cons7, cons21, cons4, cons5, cons499, cons500)
    def replacement752(m, n, p, x, b, c):
        return b**(S(1) - (m + S(1))/n)*c**m*Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n), x), x, x**n)/n
    rule752 = ReplacementRule(pattern752, replacement752)
    rubi.add(rule752.pattern, label = 752)

    pattern753 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons3, cons7, cons21, cons4, cons5, cons499, cons501)
    def replacement753(m, n, p, x, b, c):
        return b**IntPart(p)*c**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p), x)
    rule753 = ReplacementRule(pattern753, replacement753)
    rubi.add(rule753.pattern, label = 753)

    pattern754 = Pattern(Integral((c_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons3, cons7, cons21, cons4, cons5, cons18)
    def replacement754(m, n, p, x, b, c):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(b*x**n)**p, x)
    rule754 = ReplacementRule(pattern754, replacement754)
    rubi.add(rule754.pattern, label = 754)

    pattern755 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons4, cons38, cons502)
    def replacement755(m, n, p, x, b, a):
        return Int(x**(m + n*p)*(a*x**(-n) + b)**p, x)
    rule755 = ReplacementRule(pattern755, replacement755)
    rubi.add(rule755.pattern, label = 755)

    pattern756 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons503, cons66)
    def replacement756(m, n, p, x, b, c, a):
        return (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1)))
    rule756 = ReplacementRule(pattern756, replacement756)
    rubi.add(rule756.pattern, label = 756)

    pattern757 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons5, cons55, cons504, cons66)
    def replacement757(m, n, a1, p, b2, a2, x, b1, c):
        return (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1)))
    rule757 = ReplacementRule(pattern757, replacement757)
    rubi.add(rule757.pattern, label = 757)

    pattern758 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons21, cons4, cons5, cons500)
    def replacement758(m, n, p, x, b, a):
        return Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p, x), x, x**n)/n
    rule758 = ReplacementRule(pattern758, replacement758)
    rubi.add(rule758.pattern, label = 758)

    pattern759 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons21, cons4, cons5, cons55, cons505)
    def replacement759(m, n, a1, p, b2, a2, x, b1):
        return Subst(Int(x**(S(-1) + (m + S(1))/n)*(a1 + b1*x)**p*(a2 + b2*x)**p, x), x, x**n)/n
    rule759 = ReplacementRule(pattern759, replacement759)
    rubi.add(rule759.pattern, label = 759)

    pattern760 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons500)
    def replacement760(m, n, p, x, b, c, a):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x)
    rule760 = ReplacementRule(pattern760, replacement760)
    rubi.add(rule760.pattern, label = 760)

    pattern761 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons5, cons55, cons505)
    def replacement761(m, n, a1, p, b2, a2, x, b1, c):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)
    rule761 = ReplacementRule(pattern761, replacement761)
    rubi.add(rule761.pattern, label = 761)

    pattern762 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons128)
    def replacement762(m, n, p, x, b, c, a):
        return Int(ExpandIntegrand((c*x)**m*(a + b*x**n)**p, x), x)
    rule762 = ReplacementRule(pattern762, replacement762)
    rubi.add(rule762.pattern, label = 762)

    pattern763 = Pattern(Integral(x_**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons4, cons5, cons506, cons66)
    def replacement763(m, n, p, x, b, a):
        return -b*(m + n*(p + S(1)) + S(1))*Int(x**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + x**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*(m + S(1)))
    rule763 = ReplacementRule(pattern763, replacement763)
    rubi.add(rule763.pattern, label = 763)

    pattern764 = Pattern(Integral(x_**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons21, cons4, cons5, cons55, cons507, cons66)
    def replacement764(m, n, a1, p, b2, a2, x, b1):
        return -b1*b2*(m + S(2)*n*(p + S(1)) + S(1))*Int(x**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + x**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*(m + S(1)))
    rule764 = ReplacementRule(pattern764, replacement764)
    rubi.add(rule764.pattern, label = 764)

    pattern765 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons506, cons54)
    def replacement765(m, n, p, x, b, c, a):
        return (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1)))
    rule765 = ReplacementRule(pattern765, replacement765)
    rubi.add(rule765.pattern, label = 765)

    pattern766 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons5, cons55, cons507, cons54)
    def replacement766(m, n, a1, p, b2, a2, x, b1, c):
        return (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1)))
    rule766 = ReplacementRule(pattern766, replacement766)
    rubi.add(rule766.pattern, label = 766)

    def With767(m, n, p, x, b, a):
        k = GCD(m + S(1), n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p, x), x, x**k)/k
        return False
    pattern767 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons148, cons17, CustomConstraint(With767))
    rule767 = ReplacementRule(pattern767, With767)
    rubi.add(rule767.pattern, label = 767)

    def With768(m, n, a1, p, b2, a2, x, b1):
        k = GCD(m + S(1), S(2)*n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a1 + b1*x**(n/k))**p*(a2 + b2*x**(n/k))**p, x), x, x**k)/k
        return False
    pattern768 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons5, cons55, cons495, cons17, CustomConstraint(With768))
    rule768 = ReplacementRule(pattern768, With768)
    rubi.add(rule768.pattern, label = 768)

    pattern769 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons148, cons244, cons163, cons94, cons508, cons509)
    def replacement769(m, n, p, x, b, c, a):
        return -b*c**(-n)*n*p*Int((c*x)**(m + n)*(a + b*x**n)**(p + S(-1)), x)/(m + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + S(1)))
    rule769 = ReplacementRule(pattern769, replacement769)
    rubi.add(rule769.pattern, label = 769)

    pattern770 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons55, cons495, cons244, cons163, cons510, cons511)
    def replacement770(m, n, a1, p, b2, a2, x, b1, c):
        return S(2)*a1*a2*n*p*Int((c*x)**m*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(2)*n*p + S(1)) + (c*x)**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(c*(m + S(2)*n*p + S(1)))
    rule770 = ReplacementRule(pattern770, replacement770)
    rubi.add(rule770.pattern, label = 770)

    pattern771 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons148, cons244, cons163, cons512, cons509)
    def replacement771(m, n, p, x, b, c, a):
        return a*n*p*Int((c*x)**m*(a + b*x**n)**(p + S(-1)), x)/(m + n*p + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + n*p + S(1)))
    rule771 = ReplacementRule(pattern771, replacement771)
    rubi.add(rule771.pattern, label = 771)

    pattern772 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons466)
    def replacement772(x, b, a):
        return x*(a/(b*x**S(4)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(5)/4)), x)/(b*(a + b*x**S(4))**(S(1)/4))
    rule772 = ReplacementRule(pattern772, replacement772)
    rubi.add(rule772.pattern, label = 772)

    pattern773 = Pattern(Integral(x_**m_/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons466, cons513)
    def replacement773(m, x, b, a):
        return -a*(m + S(-3))*Int(x**(m + S(-4))/(a + b*x**S(4))**(S(5)/4), x)/(b*(m + S(-4))) + x**(m + S(-3))/(b*(a + b*x**S(4))**(S(1)/4)*(m + S(-4)))
    rule773 = ReplacementRule(pattern773, replacement773)
    rubi.add(rule773.pattern, label = 773)

    pattern774 = Pattern(Integral(x_**m_/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons466, cons514)
    def replacement774(m, x, b, a):
        return -b*m*Int(x**(m + S(4))/(a + b*x**S(4))**(S(5)/4), x)/(a*(m + S(1))) + x**(m + S(1))/(a*(a + b*x**S(4))**(S(1)/4)*(m + S(1)))
    rule774 = ReplacementRule(pattern774, replacement774)
    rubi.add(rule774.pattern, label = 774)

    pattern775 = Pattern(Integral(sqrt(x_*WC('c', S(1)))/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons7, cons466)
    def replacement775(x, c, b, a):
        return sqrt(c*x)*(a/(b*x**S(2)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(2)*(a/(b*x**S(2)) + S(1))**(S(5)/4)), x)/(b*(a + b*x**S(2))**(S(1)/4))
    rule775 = ReplacementRule(pattern775, replacement775)
    rubi.add(rule775.pattern, label = 775)

    pattern776 = Pattern(Integral((x_*WC('c', S(1)))**m_/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons7, cons466, cons515, cons516)
    def replacement776(m, x, b, c, a):
        return -S(2)*a*c**S(2)*(m + S(-1))*Int((c*x)**(m + S(-2))/(a + b*x**S(2))**(S(5)/4), x)/(b*(S(2)*m + S(-3))) + S(2)*c*(c*x)**(m + S(-1))/(b*(a + b*x**S(2))**(S(1)/4)*(S(2)*m + S(-3)))
    rule776 = ReplacementRule(pattern776, replacement776)
    rubi.add(rule776.pattern, label = 776)

    pattern777 = Pattern(Integral((x_*WC('c', S(1)))**m_/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons7, cons466, cons515, cons94)
    def replacement777(m, x, b, c, a):
        return -b*(S(2)*m + S(1))*Int((c*x)**(m + S(2))/(a + b*x**S(2))**(S(5)/4), x)/(S(2)*a*c**S(2)*(m + S(1))) + (c*x)**(m + S(1))/(a*c*(a + b*x**S(2))**(S(1)/4)*(m + S(1)))
    rule777 = ReplacementRule(pattern777, replacement777)
    rubi.add(rule777.pattern, label = 777)

    pattern778 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons2, cons3, cons483)
    def replacement778(x, b, a):
        return -Int(S(1)/(x**S(2)*(a + b*x**S(4))**(S(1)/4)), x)/b - S(1)/(b*x*(a + b*x**S(4))**(S(1)/4))
    rule778 = ReplacementRule(pattern778, replacement778)
    rubi.add(rule778.pattern, label = 778)

    pattern779 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons148, cons244, cons137, cons517, cons518, cons509)
    def replacement779(m, n, p, x, b, c, a):
        return -c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**(p + S(1)), x)/(b*n*(p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*n*(p + S(1)))
    rule779 = ReplacementRule(pattern779, replacement779)
    rubi.add(rule779.pattern, label = 779)

    pattern780 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons55, cons495, cons244, cons137, cons519, cons520, cons511)
    def replacement780(m, n, a1, p, b2, a2, x, b1, c):
        return -c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*b1*b2*n*(p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*b1*b2*n*(p + S(1)))
    rule780 = ReplacementRule(pattern780, replacement780)
    rubi.add(rule780.pattern, label = 780)

    pattern781 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons148, cons244, cons137, cons509)
    def replacement781(m, n, p, x, b, c, a):
        return (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1)))
    rule781 = ReplacementRule(pattern781, replacement781)
    rubi.add(rule781.pattern, label = 781)

    pattern782 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons55, cons495, cons244, cons137, cons511)
    def replacement782(m, n, a1, p, b2, a2, x, b1, c):
        return (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1)))
    rule782 = ReplacementRule(pattern782, replacement782)
    rubi.add(rule782.pattern, label = 782)

    pattern783 = Pattern(Integral(x_/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons67)
    def replacement783(x, b, a):
        return Int((x*Rt(b, S(3)) + Rt(a, S(3)))/(x**S(2)*Rt(b, S(3))**S(2) - x*Rt(a, S(3))*Rt(b, S(3)) + Rt(a, S(3))**S(2)), x)/(S(3)*Rt(a, S(3))*Rt(b, S(3))) - Int(S(1)/(x*Rt(b, S(3)) + Rt(a, S(3))), x)/(S(3)*Rt(a, S(3))*Rt(b, S(3)))
    rule783 = ReplacementRule(pattern783, replacement783)
    rubi.add(rule783.pattern, label = 783)

    def With784(m, n, x, b, a):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(Pi*m*(2*k - 1)/n) - s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) - s**(-m)*(-r)**(m + S(1))*Int(S(1)/(r + s*x), x)/(a*n)
    pattern784 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons521, cons62, cons522, cons468 )
    rule784 = ReplacementRule(pattern784, With784)
    rubi.add(rule784.pattern, label = 784)

    def With785(m, n, x, b, a):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(Pi*m*(2*k - 1)/n) + s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return -Dist(S(2)*s**(-m)*(-r)**(m + S(1))/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r**(m + S(1))*s**(-m)*Int(S(1)/(r - s*x), x)/(a*n)
    pattern785 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons523, cons62, cons522, cons469 )
    rule785 = ReplacementRule(pattern785, With785)
    rubi.add(rule785.pattern, label = 785)

    def With786(m, n, x, b, a):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(Pi*m*(2*k - 1)/n) - s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x) + Int((r*cos(Pi*m*(2*k - 1)/n) + s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return S(2)*(S(-1))**(m/S(2))*r**(m + S(2))*s**(-m)*Int(S(1)/(r**S(2) + s**S(2)*x**S(2)), x)/(a*n) + Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x)
    pattern786 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons524, cons62, cons522, cons468 )
    rule786 = ReplacementRule(pattern786, With786)
    rubi.add(rule786.pattern, label = 786)

    def With787(m, n, x, b, a):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(2*Pi*k*m/n) - s*x*cos(2*Pi*k*(m + 1)/n))/(r**2 - 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x) + Int((r*cos(2*Pi*k*m/n) + s*x*cos(2*Pi*k*(m + 1)/n))/(r**2 + 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x)
        return Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**(m + S(2))*s**(-m)*Int(S(1)/(r**S(2) - s**S(2)*x**S(2)), x)/(a*n)
    pattern787 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons524, cons62, cons522, cons469 )
    rule787 = ReplacementRule(pattern787, With787)
    rubi.add(rule787.pattern, label = 787)

    def With788(x, b, a):
        r = Numerator(Rt(a/b, S(2)))
        s = Denominator(Rt(a/b, S(2)))
        return -Int((r - s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*s) + Int((r + s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*s)
    pattern788 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons475 )
    rule788 = ReplacementRule(pattern788, With788)
    rubi.add(rule788.pattern, label = 788)

    def With789(x, b, a):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return -s*Int(S(1)/(r - s*x**S(2)), x)/(S(2)*b) + s*Int(S(1)/(r + s*x**S(2)), x)/(S(2)*b)
    pattern789 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons476 )
    rule789 = ReplacementRule(pattern789, With789)
    rubi.add(rule789.pattern, label = 789)

    def With790(m, n, x, b, a):
        r = Numerator(Rt(a/b, S(4)))
        s = Denominator(Rt(a/b, S(4)))
        return sqrt(S(2))*s**S(3)*Int(x**(m - n/S(4))/(r**S(2) - sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*b*r) - sqrt(S(2))*s**S(3)*Int(x**(m - n/S(4))/(r**S(2) + sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*b*r)
    pattern790 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons525, cons62, cons522, cons478 )
    rule790 = ReplacementRule(pattern790, With790)
    rubi.add(rule790.pattern, label = 790)

    def With791(m, n, x, b, a):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return r*Int(x**m/(r - s*x**(n/S(2))), x)/(S(2)*a) + r*Int(x**m/(r + s*x**(n/S(2))), x)/(S(2)*a)
    pattern791 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons525, cons62, cons526, cons476 )
    rule791 = ReplacementRule(pattern791, With791)
    rubi.add(rule791.pattern, label = 791)

    def With792(m, n, x, b, a):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return -s*Int(x**(m - n/S(2))/(r - s*x**(n/S(2))), x)/(S(2)*b) + s*Int(x**(m - n/S(2))/(r + s*x**(n/S(2))), x)/(S(2)*b)
    pattern792 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons525, cons62, cons527, cons476 )
    rule792 = ReplacementRule(pattern792, With792)
    rubi.add(rule792.pattern, label = 792)

    pattern793 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons528, cons529)
    def replacement793(m, n, x, b, a):
        return Int(PolynomialDivide(x**m, a + b*x**n, x), x)
    rule793 = ReplacementRule(pattern793, replacement793)
    rubi.add(rule793.pattern, label = 793)

    def With794(x, b, a):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return sqrt(S(2))*s*Int(S(1)/sqrt(a + b*x**S(3)), x)/(r*sqrt(sqrt(S(3)) + S(2))) + Int((r*x + s*(-sqrt(S(3)) + S(1)))/sqrt(a + b*x**S(3)), x)/r
    pattern794 = Pattern(Integral(x_/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons481 )
    rule794 = ReplacementRule(pattern794, With794)
    rubi.add(rule794.pattern, label = 794)

    def With795(x, b, a):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return -sqrt(S(2))*s*Int(S(1)/sqrt(a + b*x**S(3)), x)/(r*sqrt(-sqrt(S(3)) + S(2))) + Int((r*x + s*(S(1) + sqrt(S(3))))/sqrt(a + b*x**S(3)), x)/r
    pattern795 = Pattern(Integral(x_/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons482 )
    rule795 = ReplacementRule(pattern795, With795)
    rubi.add(rule795.pattern, label = 795)

    def With796(x, b, a):
        q = Rt(b/a, S(2))
        return -Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x)/q + Int(S(1)/sqrt(a + b*x**S(4)), x)/q
    pattern796 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons466 )
    rule796 = ReplacementRule(pattern796, With796)
    rubi.add(rule796.pattern, label = 796)

    def With797(x, b, a):
        q = Rt(-b/a, S(2))
        return -Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x)/q + Int(S(1)/sqrt(a + b*x**S(4)), x)/q
    pattern797 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons484, cons105 )
    rule797 = ReplacementRule(pattern797, With797)
    rubi.add(rule797.pattern, label = 797)

    def With798(x, b, a):
        q = Rt(-b/a, S(2))
        return Int((q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x)/q - Int(S(1)/sqrt(a + b*x**S(4)), x)/q
    pattern798 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons483 )
    rule798 = ReplacementRule(pattern798, With798)
    rubi.add(rule798.pattern, label = 798)

    def With799(x, b, a):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return s**S(2)*(S(-1) + sqrt(S(3)))*Int(S(1)/sqrt(a + b*x**S(6)), x)/(S(2)*r**S(2)) - Int((-S(2)*r**S(2)*x**S(4) + s**S(2)*(S(-1) + sqrt(S(3))))/sqrt(a + b*x**S(6)), x)/(S(2)*r**S(2))
    pattern799 = Pattern(Integral(x_**S(4)/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), cons2, cons3, cons67 )
    rule799 = ReplacementRule(pattern799, With799)
    rubi.add(rule799.pattern, label = 799)

    pattern800 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), cons2, cons3, cons67)
    def replacement800(x, b, a):
        return -Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/(S(2)*Rt(b/a, S(4))) + Int((x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/(S(2)*Rt(b/a, S(4)))
    rule800 = ReplacementRule(pattern800, replacement800)
    rubi.add(rule800.pattern, label = 800)

    pattern801 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4), x_), cons2, cons3, cons466)
    def replacement801(x, b, a):
        return -a*Int(x**S(2)/(a + b*x**S(4))**(S(5)/4), x)/S(2) + x**S(3)/(S(2)*(a + b*x**S(4))**(S(1)/4))
    rule801 = ReplacementRule(pattern801, replacement801)
    rubi.add(rule801.pattern, label = 801)

    pattern802 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4), x_), cons2, cons3, cons483)
    def replacement802(x, b, a):
        return a*Int(S(1)/(x**S(2)*(a + b*x**S(4))**(S(1)/4)), x)/(S(2)*b) + (a + b*x**S(4))**(S(3)/4)/(S(2)*b*x)
    rule802 = ReplacementRule(pattern802, replacement802)
    rubi.add(rule802.pattern, label = 802)

    pattern803 = Pattern(Integral(S(1)/(x_**S(2)*(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)), x_), cons2, cons3, cons466)
    def replacement803(x, b, a):
        return -b*Int(x**S(2)/(a + b*x**S(4))**(S(5)/4), x) - S(1)/(x*(a + b*x**S(4))**(S(1)/4))
    rule803 = ReplacementRule(pattern803, replacement803)
    rubi.add(rule803.pattern, label = 803)

    pattern804 = Pattern(Integral(S(1)/(x_**S(2)*(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)), x_), cons2, cons3, cons483)
    def replacement804(x, b, a):
        return x*(a/(b*x**S(4)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(1)/4)), x)/(a + b*x**S(4))**(S(1)/4)
    rule804 = ReplacementRule(pattern804, replacement804)
    rubi.add(rule804.pattern, label = 804)

    pattern805 = Pattern(Integral(sqrt(c_*x_)/(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4), x_), cons2, cons3, cons7, cons466)
    def replacement805(x, c, b, a):
        return -a*Int(sqrt(c*x)/(a + b*x**S(2))**(S(5)/4), x)/S(2) + x*sqrt(c*x)/(a + b*x**S(2))**(S(1)/4)
    rule805 = ReplacementRule(pattern805, replacement805)
    rubi.add(rule805.pattern, label = 805)

    pattern806 = Pattern(Integral(sqrt(c_*x_)/(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4), x_), cons2, cons3, cons7, cons483)
    def replacement806(x, c, b, a):
        return a*c**S(2)*Int(S(1)/((c*x)**(S(3)/2)*(a + b*x**S(2))**(S(1)/4)), x)/(S(2)*b) + c*(a + b*x**S(2))**(S(3)/4)/(b*sqrt(c*x))
    rule806 = ReplacementRule(pattern806, replacement806)
    rubi.add(rule806.pattern, label = 806)

    pattern807 = Pattern(Integral(S(1)/((x_*WC('c', S(1)))**(S(3)/2)*(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)), x_), cons2, cons3, cons7, cons466)
    def replacement807(x, c, b, a):
        return -b*Int(sqrt(c*x)/(a + b*x**S(2))**(S(5)/4), x)/c**S(2) - S(2)/(c*sqrt(c*x)*(a + b*x**S(2))**(S(1)/4))
    rule807 = ReplacementRule(pattern807, replacement807)
    rubi.add(rule807.pattern, label = 807)

    pattern808 = Pattern(Integral(S(1)/((x_*WC('c', S(1)))**(S(3)/2)*(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)), x_), cons2, cons3, cons7, cons483)
    def replacement808(x, c, b, a):
        return sqrt(c*x)*(a/(b*x**S(2)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(2)*(a/(b*x**S(2)) + S(1))**(S(1)/4)), x)/(c**S(2)*(a + b*x**S(2))**(S(1)/4))
    rule808 = ReplacementRule(pattern808, replacement808)
    rubi.add(rule808.pattern, label = 808)

    pattern809 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons148, cons31, cons530, cons512, cons509)
    def replacement809(m, n, p, x, b, c, a):
        return -a*c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**p, x)/(b*(m + n*p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*p + S(1)))
    rule809 = ReplacementRule(pattern809, replacement809)
    rubi.add(rule809.pattern, label = 809)

    pattern810 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons148, cons531, cons512, cons532)
    def replacement810(m, n, p, x, b, c, a):
        return -a*c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**p, x)/(b*(m + n*p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*p + S(1)))
    rule810 = ReplacementRule(pattern810, replacement810)
    rubi.add(rule810.pattern, label = 810)

    pattern811 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons5, cons55, cons495, cons31, cons529, cons510, cons511)
    def replacement811(m, n, a1, p, b2, a2, x, b1, c):
        return -a1*a2*c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(b1*b2*(m + S(2)*n*p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(b1*b2*(m + S(2)*n*p + S(1)))
    rule811 = ReplacementRule(pattern811, replacement811)
    rubi.add(rule811.pattern, label = 811)

    pattern812 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons5, cons55, cons495, cons533, cons510, cons534)
    def replacement812(m, n, a1, p, b2, a2, x, b1, c):
        return -a1*a2*c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(b1*b2*(m + S(2)*n*p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(b1*b2*(m + S(2)*n*p + S(1)))
    rule812 = ReplacementRule(pattern812, replacement812)
    rubi.add(rule812.pattern, label = 812)

    pattern813 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons148, cons31, cons94, cons509)
    def replacement813(m, n, p, x, b, c, a):
        return -b*c**(-n)*(m + n*(p + S(1)) + S(1))*Int((c*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1)))
    rule813 = ReplacementRule(pattern813, replacement813)
    rubi.add(rule813.pattern, label = 813)

    pattern814 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons148, cons535, cons532)
    def replacement814(m, n, p, x, b, c, a):
        return -b*c**(-n)*(m + n*(p + S(1)) + S(1))*Int((c*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1)))
    rule814 = ReplacementRule(pattern814, replacement814)
    rubi.add(rule814.pattern, label = 814)

    pattern815 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons5, cons55, cons495, cons31, cons94, cons511)
    def replacement815(m, n, a1, p, b2, a2, x, b1, c):
        return -b1*b2*c**(-S(2)*n)*(m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1)))
    rule815 = ReplacementRule(pattern815, replacement815)
    rubi.add(rule815.pattern, label = 815)

    pattern816 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons5, cons55, cons495, cons536, cons534)
    def replacement816(m, n, a1, p, b2, a2, x, b1, c):
        return -b1*b2*c**(-S(2)*n)*(m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1)))
    rule816 = ReplacementRule(pattern816, replacement816)
    rubi.add(rule816.pattern, label = 816)

    def With817(m, n, p, x, b, c, a):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*c**(-n)*x**(k*n))**p, x), x, (c*x)**(S(1)/k))/c
    pattern817 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons148, cons367, cons509 )
    rule817 = ReplacementRule(pattern817, With817)
    rubi.add(rule817.pattern, label = 817)

    def With818(m, n, a1, p, b2, a2, x, b1, c):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a1 + b1*c**(-n)*x**(k*n))**p*(a2 + b2*c**(-n)*x**(k*n))**p, x), x, (c*x)**(S(1)/k))/c
    pattern818 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons5, cons55, cons495, cons367, cons511 )
    rule818 = ReplacementRule(pattern818, With818)
    rubi.add(rule818.pattern, label = 818)

    pattern819 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons148, cons13, cons485, cons486, cons537)
    def replacement819(m, n, p, x, b, a):
        return a**(p + (m + S(1))/n)*Subst(Int(x**m*(-b*x**n + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x*(a + b*x**n)**(-S(1)/n))
    rule819 = ReplacementRule(pattern819, replacement819)
    rubi.add(rule819.pattern, label = 819)

    pattern820 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons55, cons495, cons13, cons485, cons486, cons538)
    def replacement820(m, n, a1, p, b2, a2, x, b1):
        return (a1*a2)**(p + (m + S(1))/(S(2)*n))*Subst(Int(x**m*(-b1*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x*(a1 + b1*x**n)**(-S(1)/(S(2)*n))*(a2 + b2*x**n)**(-S(1)/(S(2)*n)))
    rule820 = ReplacementRule(pattern820, replacement820)
    rubi.add(rule820.pattern, label = 820)

    pattern821 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons148, cons13, cons485, cons486, cons17, cons539)
    def replacement821(m, n, p, x, b, a):
        return (a/(a + b*x**n))**(p + (m + S(1))/n)*(a + b*x**n)**(p + (m + S(1))/n)*Subst(Int(x**m*(-b*x**n + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x*(a + b*x**n)**(-S(1)/n))
    rule821 = ReplacementRule(pattern821, replacement821)
    rubi.add(rule821.pattern, label = 821)

    pattern822 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons55, cons495, cons13, cons485, cons486, cons17, cons540)
    def replacement822(m, n, a1, p, b2, a2, x, b1):
        return (a1/(a1 + b1*x**n))**(p + (m + S(1))/(S(2)*n))*(a2/(a2 + b2*x**n))**(p + (m + S(1))/(S(2)*n))*(a1 + b1*x**n)**(p + (m + S(1))/(S(2)*n))*(a2 + b2*x**n)**(p + (m + S(1))/(S(2)*n))*Subst(Int(x**m*(-b1*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x*(a1 + b1*x**n)**(-S(1)/(S(2)*n))*(a2 + b2*x**n)**(-S(1)/(S(2)*n)))
    rule822 = ReplacementRule(pattern822, replacement822)
    rubi.add(rule822.pattern, label = 822)

    pattern823 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons196, cons17)
    def replacement823(m, n, p, x, b, a):
        return -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p, x), x, S(1)/x)
    rule823 = ReplacementRule(pattern823, replacement823)
    rubi.add(rule823.pattern, label = 823)

    pattern824 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons5, cons55, cons497, cons17)
    def replacement824(m, n, a1, p, b2, a2, x, b1):
        return -Subst(Int(x**(-m + S(-2))*(a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p, x), x, S(1)/x)
    rule824 = ReplacementRule(pattern824, replacement824)
    rubi.add(rule824.pattern, label = 824)

    def With825(m, n, p, x, b, c, a):
        k = Denominator(m)
        return -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*c**(-n)*x**(-k*n))**p, x), x, (c*x)**(-S(1)/k))/c
    pattern825 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons196, cons367 )
    rule825 = ReplacementRule(pattern825, With825)
    rubi.add(rule825.pattern, label = 825)

    def With826(m, n, a1, p, b2, a2, x, b1, c):
        k = Denominator(m)
        return -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a1 + b1*c**(-n)*x**(-k*n))**p*(a2 + b2*c**(-n)*x**(-k*n))**p, x), x, (c*x)**(-S(1)/k))/c
    pattern826 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons5, cons55, cons497, cons367 )
    rule826 = ReplacementRule(pattern826, With826)
    rubi.add(rule826.pattern, label = 826)

    pattern827 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons196, cons356)
    def replacement827(m, n, p, x, b, c, a):
        return -(c*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p, x), x, S(1)/x)
    rule827 = ReplacementRule(pattern827, replacement827)
    rubi.add(rule827.pattern, label = 827)

    pattern828 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons5, cons55, cons497, cons356)
    def replacement828(m, n, a1, p, b2, a2, x, b1, c):
        return -(c*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p, x), x, S(1)/x)
    rule828 = ReplacementRule(pattern828, replacement828)
    rubi.add(rule828.pattern, label = 828)

    def With829(m, n, p, x, b, a):
        k = Denominator(n)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n))**p, x), x, x**(S(1)/k))
    pattern829 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons5, cons489 )
    rule829 = ReplacementRule(pattern829, With829)
    rubi.add(rule829.pattern, label = 829)

    def With830(m, n, a1, p, b2, a2, x, b1):
        k = Denominator(S(2)*n)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a1 + b1*x**(k*n))**p*(a2 + b2*x**(k*n))**p, x), x, x**(S(1)/k))
    pattern830 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons21, cons5, cons55, cons498 )
    rule830 = ReplacementRule(pattern830, With830)
    rubi.add(rule830.pattern, label = 830)

    pattern831 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons489)
    def replacement831(m, n, p, x, b, c, a):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x)
    rule831 = ReplacementRule(pattern831, replacement831)
    rubi.add(rule831.pattern, label = 831)

    pattern832 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons5, cons55, cons498)
    def replacement832(m, n, a1, p, b2, a2, x, b1, c):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)
    rule832 = ReplacementRule(pattern832, replacement832)
    rubi.add(rule832.pattern, label = 832)

    pattern833 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons4, cons5, cons541, cons23)
    def replacement833(m, n, p, x, b, a):
        return Subst(Int((a + b*x**(n/(m + S(1))))**p, x), x, x**(m + S(1)))/(m + S(1))
    rule833 = ReplacementRule(pattern833, replacement833)
    rubi.add(rule833.pattern, label = 833)

    pattern834 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons21, cons4, cons5, cons55, cons542, cons543)
    def replacement834(m, n, a1, p, b2, a2, x, b1):
        return Subst(Int((a1 + b1*x**(n/(m + S(1))))**p*(a2 + b2*x**(n/(m + S(1))))**p, x), x, x**(m + S(1)))/(m + S(1))
    rule834 = ReplacementRule(pattern834, replacement834)
    rubi.add(rule834.pattern, label = 834)

    pattern835 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons541, cons23)
    def replacement835(m, n, p, x, b, c, a):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x)
    rule835 = ReplacementRule(pattern835, replacement835)
    rubi.add(rule835.pattern, label = 835)

    pattern836 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons5, cons55, cons542, cons543)
    def replacement836(m, n, a1, p, b2, a2, x, b1, c):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)
    rule836 = ReplacementRule(pattern836, replacement836)
    rubi.add(rule836.pattern, label = 836)

    pattern837 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons4, cons544, cons13, cons163)
    def replacement837(m, n, p, x, b, a):
        return -b*n*p*Int(x**(m + n)*(a + b*x**n)**(p + S(-1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*x**n)**p/(m + S(1))
    rule837 = ReplacementRule(pattern837, replacement837)
    rubi.add(rule837.pattern, label = 837)

    pattern838 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons21, cons4, cons55, cons545, cons13, cons163)
    def replacement838(m, n, a1, p, b2, a2, x, b1):
        return -S(2)*b1*b2*n*p*Int(x**(m + n)*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(1)) + x**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(m + S(1))
    rule838 = ReplacementRule(pattern838, replacement838)
    rubi.add(rule838.pattern, label = 838)

    pattern839 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons544, cons13, cons163)
    def replacement839(m, n, p, x, b, c, a):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x)
    rule839 = ReplacementRule(pattern839, replacement839)
    rubi.add(rule839.pattern, label = 839)

    pattern840 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons55, cons545, cons13, cons163)
    def replacement840(m, n, a1, p, b2, a2, x, b1, c):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)
    rule840 = ReplacementRule(pattern840, replacement840)
    rubi.add(rule840.pattern, label = 840)

    pattern841 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons546, cons13, cons163, cons512)
    def replacement841(m, n, p, x, b, c, a):
        return a*n*p*Int((c*x)**m*(a + b*x**n)**(p + S(-1)), x)/(m + n*p + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + n*p + S(1)))
    rule841 = ReplacementRule(pattern841, replacement841)
    rubi.add(rule841.pattern, label = 841)

    pattern842 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons55, cons547, cons13, cons163, cons510)
    def replacement842(m, n, a1, p, b2, a2, x, b1, c):
        return S(2)*a1*a2*n*p*Int((c*x)**m*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(2)*n*p + S(1)) + (c*x)**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(c*(m + S(2)*n*p + S(1)))
    rule842 = ReplacementRule(pattern842, replacement842)
    rubi.add(rule842.pattern, label = 842)

    def With843(m, n, p, x, b, a):
        k = Denominator(p)
        return a**(p + (m + S(1))/n)*k*Subst(Int(x**(k*(m + S(1))/n + S(-1))*(-b*x**k + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x**(n/k)*(a + b*x**n)**(-S(1)/k))/n
    pattern843 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons4, cons546, cons13, cons485 )
    rule843 = ReplacementRule(pattern843, With843)
    rubi.add(rule843.pattern, label = 843)

    def With844(m, n, a1, p, b2, a2, x, b1):
        k = Denominator(p)
        return k*(a1*a2)**(p + (m + S(1))/(S(2)*n))*Subst(Int(x**(k*(m + S(1))/(S(2)*n) + S(-1))*(-b1*x**k + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**k + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x**(S(2)*n/k)*(a1 + b1*x**n)**(-S(1)/k)*(a2 + b2*x**n)**(-S(1)/k))/(S(2)*n)
    pattern844 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons21, cons4, cons55, cons547, cons13, cons485 )
    rule844 = ReplacementRule(pattern844, With844)
    rubi.add(rule844.pattern, label = 844)

    pattern845 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons546, cons13, cons485)
    def replacement845(m, n, p, x, b, c, a):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x)
    rule845 = ReplacementRule(pattern845, replacement845)
    rubi.add(rule845.pattern, label = 845)

    pattern846 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons55, cons547, cons13, cons485)
    def replacement846(m, n, a1, p, b2, a2, x, b1, c):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)
    rule846 = ReplacementRule(pattern846, replacement846)
    rubi.add(rule846.pattern, label = 846)

    pattern847 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons546, cons13, cons137)
    def replacement847(m, n, p, x, b, c, a):
        return (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1)))
    rule847 = ReplacementRule(pattern847, replacement847)
    rubi.add(rule847.pattern, label = 847)

    pattern848 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons55, cons546, cons13, cons137)
    def replacement848(m, n, a1, p, b2, a2, x, b1, c):
        return (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1)))
    rule848 = ReplacementRule(pattern848, replacement848)
    rubi.add(rule848.pattern, label = 848)

    def With849(m, n, x, b, a):
        mn = m - n
        return -a*Int(x**mn/(a + b*x**n), x)/b + x**(mn + S(1))/(b*(mn + S(1)))
    pattern849 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons21, cons4, cons548, cons531 )
    rule849 = ReplacementRule(pattern849, With849)
    rubi.add(rule849.pattern, label = 849)

    pattern850 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons21, cons4, cons548, cons535)
    def replacement850(m, n, x, b, a):
        return -b*Int(x**(m + n)/(a + b*x**n), x)/a + x**(m + S(1))/(a*(m + S(1)))
    rule850 = ReplacementRule(pattern850, replacement850)
    rubi.add(rule850.pattern, label = 850)

    pattern851 = Pattern(Integral((c_*x_)**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons21, cons4, cons548, cons549)
    def replacement851(m, n, x, b, c, a):
        return c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m/(a + b*x**n), x)
    rule851 = ReplacementRule(pattern851, replacement851)
    rubi.add(rule851.pattern, label = 851)

    pattern852 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons357, cons550)
    def replacement852(m, n, p, x, b, c, a):
        return a**p*(c*x)**(m + S(1))*Hypergeometric2F1(-p, (m + S(1))/n, S(1) + (m + S(1))/n, -b*x**n/a)/(c*(m + S(1)))
    rule852 = ReplacementRule(pattern852, replacement852)
    rubi.add(rule852.pattern, label = 852)

    pattern853 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons357, cons551)
    def replacement853(m, n, p, x, b, c, a):
        return a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((c*x)**m*(S(1) + b*x**n/a)**p, x)
    rule853 = ReplacementRule(pattern853, replacement853)
    rubi.add(rule853.pattern, label = 853)

    pattern854 = Pattern(Integral(x_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons5, cons552, cons17, cons553)
    def replacement854(m, n, p, v, x, b, a):
        return Coefficient(v, x, S(1))**(-m + S(-1))*Subst(Int(SimplifyIntegrand((a + b*x**n)**p*(x - Coefficient(v, x, S(0)))**m, x), x), x, v)
    rule854 = ReplacementRule(pattern854, replacement854)
    rubi.add(rule854.pattern, label = 854)

    pattern855 = Pattern(Integral(u_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons21, cons4, cons5, cons554)
    def replacement855(m, n, p, v, x, b, u, a):
        return u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p, x), x, v)/Coefficient(v, x, S(1))
    rule855 = ReplacementRule(pattern855, replacement855)
    rubi.add(rule855.pattern, label = 855)

    pattern856 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons5, cons55, cons147)
    def replacement856(m, n, a1, p, b2, a2, x, b1, c):
        return (a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p))*Int((c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x)
    rule856 = ReplacementRule(pattern856, replacement856)
    rubi.add(rule856.pattern, label = 856)

    pattern857 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons71, cons555)
    def replacement857(n, p, q, x, b, d, c, a):
        return Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x)
    rule857 = ReplacementRule(pattern857, replacement857)
    rubi.add(rule857.pattern, label = 857)

    pattern858 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons71, cons220, cons502)
    def replacement858(n, p, q, x, b, d, c, a):
        return Int(x**(n*(p + q))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q, x)
    rule858 = ReplacementRule(pattern858, replacement858)
    rubi.add(rule858.pattern, label = 858)

    pattern859 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons5, cons50, cons71, cons196)
    def replacement859(n, p, q, x, b, d, c, a):
        return -Subst(Int((a + b*x**(-n))**p*(c + d*x**(-n))**q/x**S(2), x), x, S(1)/x)
    rule859 = ReplacementRule(pattern859, replacement859)
    rubi.add(rule859.pattern, label = 859)

    def With860(n, p, q, x, b, d, c, a):
        g = Denominator(n)
        return g*Subst(Int(x**(g + S(-1))*(a + b*x**(g*n))**p*(c + d*x**(g*n))**q, x), x, x**(S(1)/g))
    pattern860 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons5, cons50, cons71, cons489 )
    rule860 = ReplacementRule(pattern860, With860)
    rubi.add(rule860.pattern, label = 860)

    pattern861 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons71, cons556, cons85)
    def replacement861(n, p, x, b, d, c, a):
        return Subst(Int(S(1)/(c - x**n*(-a*d + b*c)), x), x, x*(a + b*x**n)**(-S(1)/n))
    rule861 = ReplacementRule(pattern861, replacement861)
    rubi.add(rule861.pattern, label = 861)

    pattern862 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons71, cons557, cons395, cons403, cons54)
    def replacement862(n, q, p, x, b, d, c, a):
        return -c*q*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1)), x)/(a*(p + S(1))) - x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*n*(p + S(1)))
    rule862 = ReplacementRule(pattern862, replacement862)
    rubi.add(rule862.pattern, label = 862)

    pattern863 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons50, cons71, cons557, cons63)
    def replacement863(n, p, q, x, b, d, c, a):
        return a**p*c**(-p + S(-1))*x*(c + d*x**n)**(-S(1)/n)*Hypergeometric2F1(S(1)/n, -p, S(1) + S(1)/n, -x**n*(-a*d + b*c)/(a*(c + d*x**n)))
    rule863 = ReplacementRule(pattern863, replacement863)
    rubi.add(rule863.pattern, label = 863)

    pattern864 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons71, cons557)
    def replacement864(n, p, q, x, b, d, c, a):
        return x*(c*(a + b*x**n)/(a*(c + d*x**n)))**(-p)*(a + b*x**n)**p*(c + d*x**n)**(-p - S(1)/n)*Hypergeometric2F1(S(1)/n, -p, S(1) + S(1)/n, -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/c
    rule864 = ReplacementRule(pattern864, replacement864)
    rubi.add(rule864.pattern, label = 864)

    pattern865 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons71, cons558, cons559)
    def replacement865(n, p, q, x, b, d, c, a):
        return x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c)
    rule865 = ReplacementRule(pattern865, replacement865)
    rubi.add(rule865.pattern, label = 865)

    pattern866 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons50, cons71, cons558, cons560, cons54)
    def replacement866(n, p, q, x, b, d, c, a):
        return -b*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*n*(p + S(1))*(-a*d + b*c)) + (b*c + n*(p + S(1))*(-a*d + b*c))*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q, x)/(a*n*(p + S(1))*(-a*d + b*c))
    rule866 = ReplacementRule(pattern866, replacement866)
    rubi.add(rule866.pattern, label = 866)

    pattern867 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons71, cons561)
    def replacement867(n, p, x, b, d, c, a):
        return c*x*(a + b*x**n)**(p + S(1))/a
    rule867 = ReplacementRule(pattern867, replacement867)
    rubi.add(rule867.pattern, label = 867)

    pattern868 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons71, cons562)
    def replacement868(n, p, x, b, d, c, a):
        return -x*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*n*(p + S(1))) - (a*d - b*c*(n*(p + S(1)) + S(1)))*Int((a + b*x**n)**(p + S(1)), x)/(a*b*n*(p + S(1)))
    rule868 = ReplacementRule(pattern868, replacement868)
    rubi.add(rule868.pattern, label = 868)

    pattern869 = Pattern(Integral((c_ + x_**n_*WC('d', S(1)))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons4, cons71, cons87, cons463)
    def replacement869(n, x, b, d, c, a):
        return c*x/a - (-a*d + b*c)*Int(S(1)/(a*x**(-n) + b), x)/a
    rule869 = ReplacementRule(pattern869, replacement869)
    rubi.add(rule869.pattern, label = 869)

    pattern870 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons4, cons71, cons563)
    def replacement870(n, p, x, b, d, c, a):
        return d*x*(a + b*x**n)**(p + S(1))/(b*(n*(p + S(1)) + S(1))) - (a*d - b*c*(n*(p + S(1)) + S(1)))*Int((a + b*x**n)**p, x)/(b*(n*(p + S(1)) + S(1)))
    rule870 = ReplacementRule(pattern870, replacement870)
    rubi.add(rule870.pattern, label = 870)

    pattern871 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons71, cons464, cons564, cons565)
    def replacement871(n, p, q, x, b, d, c, a):
        return Int(PolynomialDivide((a + b*x**n)**p, (c + d*x**n)**(-q), x), x)
    rule871 = ReplacementRule(pattern871, replacement871)
    rubi.add(rule871.pattern, label = 871)

    pattern872 = Pattern(Integral(S(1)/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons4, cons71)
    def replacement872(n, x, b, d, c, a):
        return b*Int(S(1)/(a + b*x**n), x)/(-a*d + b*c) - d*Int(S(1)/(c + d*x**n), x)/(-a*d + b*c)
    rule872 = ReplacementRule(pattern872, replacement872)
    rubi.add(rule872.pattern, label = 872)

    pattern873 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/3)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71, cons566, cons466)
    def replacement873(x, b, d, c, a):
        return sqrt(S(3))*Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(-x*Rt(b/a, S(2)) + sqrt(S(3)))), x)/(S(2)*c) + sqrt(S(3))*Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(x*Rt(b/a, S(2)) + sqrt(S(3)))), x)/(S(2)*c)
    rule873 = ReplacementRule(pattern873, replacement873)
    rubi.add(rule873.pattern, label = 873)

    pattern874 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/3)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71, cons566, cons483)
    def replacement874(x, b, d, c, a):
        return Int((-x*Rt(-b/a, S(2)) + S(3))/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/S(6) + Int((x*Rt(-b/a, S(2)) + S(3))/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/S(6)
    rule874 = ReplacementRule(pattern874, replacement874)
    rubi.add(rule874.pattern, label = 874)

    pattern875 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(2)/3)/(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons71, cons566)
    def replacement875(x, b, d, c, a):
        return b*Int((a + b*x**S(2))**(S(-1)/3), x)/d - (-a*d + b*c)*Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/d
    rule875 = ReplacementRule(pattern875, replacement875)
    rubi.add(rule875.pattern, label = 875)

    pattern876 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71)
    def replacement876(x, b, d, c, a):
        return sqrt(-b*x**S(2)/a)*Subst(Int(S(1)/(sqrt(-b*x/a)*(a + b*x)**(S(1)/4)*(c + d*x)), x), x, x**S(2))/(S(2)*x)
    rule876 = ReplacementRule(pattern876, replacement876)
    rubi.add(rule876.pattern, label = 876)

    pattern877 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71)
    def replacement877(x, b, d, c, a):
        return sqrt(-b*x**S(2)/a)*Subst(Int(S(1)/(sqrt(-b*x/a)*(a + b*x)**(S(3)/4)*(c + d*x)), x), x, x**S(2))/(S(2)*x)
    rule877 = ReplacementRule(pattern877, replacement877)
    rubi.add(rule877.pattern, label = 877)

    pattern878 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**WC('p', S(1))/(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons71, cons13, cons163, cons567)
    def replacement878(p, x, b, d, c, a):
        return b*Int((a + b*x**S(2))**(p + S(-1)), x)/d - (-a*d + b*c)*Int((a + b*x**S(2))**(p + S(-1))/(c + d*x**S(2)), x)/d
    rule878 = ReplacementRule(pattern878, replacement878)
    rubi.add(rule878.pattern, label = 878)

    pattern879 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**p_/(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons71, cons13, cons137, cons568, cons569)
    def replacement879(p, x, b, d, c, a):
        return b*Int((a + b*x**S(2))**p, x)/(-a*d + b*c) - d*Int((a + b*x**S(2))**(p + S(1))/(c + d*x**S(2)), x)/(-a*d + b*c)
    rule879 = ReplacementRule(pattern879, replacement879)
    rubi.add(rule879.pattern, label = 879)

    pattern880 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons70, cons570)
    def replacement880(x, b, d, c, a):
        return a*Subst(Int(S(1)/(-S(4)*a*b*x**S(4) + S(1)), x), x, x/sqrt(a + b*x**S(4)))/c
    rule880 = ReplacementRule(pattern880, replacement880)
    rubi.add(rule880.pattern, label = 880)

    def With881(x, b, d, c, a):
        q = Rt(-a*b, S(4))
        return a*ArcTan(q*x*(a + q**S(2)*x**S(2))/(a*sqrt(a + b*x**S(4))))/(S(2)*c*q) + a*atanh(q*x*(a - q**S(2)*x**S(2))/(a*sqrt(a + b*x**S(4))))/(S(2)*c*q)
    pattern881 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons70, cons571 )
    rule881 = ReplacementRule(pattern881, With881)
    rubi.add(rule881.pattern, label = 881)

    pattern882 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons71)
    def replacement882(x, b, d, c, a):
        return b*Int(S(1)/sqrt(a + b*x**S(4)), x)/d - (-a*d + b*c)*Int(S(1)/(sqrt(a + b*x**S(4))*(c + d*x**S(4))), x)/d
    rule882 = ReplacementRule(pattern882, replacement882)
    rubi.add(rule882.pattern, label = 882)

    pattern883 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)/(c_ + x_**S(4)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons71)
    def replacement883(x, b, d, c, a):
        return sqrt(a/(a + b*x**S(4)))*sqrt(a + b*x**S(4))*Subst(Int(S(1)/((c - x**S(4)*(-a*d + b*c))*sqrt(-b*x**S(4) + S(1))), x), x, x/(a + b*x**S(4))**(S(1)/4))
    rule883 = ReplacementRule(pattern883, replacement883)
    rubi.add(rule883.pattern, label = 883)

    pattern884 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**p_/(c_ + x_**S(4)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons71, cons13, cons572)
    def replacement884(p, x, b, d, c, a):
        return b*Int((a + b*x**S(4))**(p + S(-1)), x)/d - (-a*d + b*c)*Int((a + b*x**S(4))**(p + S(-1))/(c + d*x**S(4)), x)/d
    rule884 = ReplacementRule(pattern884, replacement884)
    rubi.add(rule884.pattern, label = 884)

    pattern885 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('b', S(1)))*(c_ + x_**S(4)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71)
    def replacement885(x, b, d, c, a):
        return Int(S(1)/(sqrt(a + b*x**S(4))*(-x**S(2)*Rt(-d/c, S(2)) + S(1))), x)/(S(2)*c) + Int(S(1)/(sqrt(a + b*x**S(4))*(x**S(2)*Rt(-d/c, S(2)) + S(1))), x)/(S(2)*c)
    rule885 = ReplacementRule(pattern885, replacement885)
    rubi.add(rule885.pattern, label = 885)

    pattern886 = Pattern(Integral(S(1)/((a_ + x_**S(4)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(4)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71)
    def replacement886(x, b, d, c, a):
        return b*Int((a + b*x**S(4))**(S(-3)/4), x)/(-a*d + b*c) - d*Int((a + b*x**S(4))**(S(1)/4)/(c + d*x**S(4)), x)/(-a*d + b*c)
    rule886 = ReplacementRule(pattern886, replacement886)
    rubi.add(rule886.pattern, label = 886)

    pattern887 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons27, cons466, cons573)
    def replacement887(x, b, d, c, a):
        return sqrt(a + b*x**S(2))*EllipticE(ArcTan(x*Rt(d/c, S(2))), S(1) - b*c/(a*d))/(c*sqrt(c*(a + b*x**S(2))/(a*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2)))
    rule887 = ReplacementRule(pattern887, replacement887)
    rubi.add(rule887.pattern, label = 887)

    pattern888 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons71, cons402, cons137, cons574, cons575)
    def replacement888(n, p, q, x, b, d, c, a):
        return -x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*n*(p + S(1))) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(n*(p + S(1)) + S(1)) + d*x**n*(n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1)))
    rule888 = ReplacementRule(pattern888, replacement888)
    rubi.add(rule888.pattern, label = 888)

    pattern889 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons71, cons402, cons137, cons576, cons575)
    def replacement889(n, p, q, x, b, d, c, a):
        return x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(a*d - b*c)/(a*b*n*(p + S(1))) - Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(a*d - b*c*(n*(p + S(1)) + S(1))) + d*x**n*(a*d*(n*(q + S(-1)) + S(1)) - b*c*(n*(p + q) + S(1))), x), x)/(a*b*n*(p + S(1)))
    rule889 = ReplacementRule(pattern889, replacement889)
    rubi.add(rule889.pattern, label = 889)

    pattern890 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons50, cons71, cons13, cons137, cons405, cons575)
    def replacement890(n, p, q, x, b, d, c, a):
        return -b*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*n*(p + S(1))*(-a*d + b*c)) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c + b*d*x**n*(n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c))
    rule890 = ReplacementRule(pattern890, replacement890)
    rubi.add(rule890.pattern, label = 890)

    pattern891 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons71, cons148, cons220, cons577)
    def replacement891(n, p, q, x, b, d, c, a):
        return Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x)
    rule891 = ReplacementRule(pattern891, replacement891)
    rubi.add(rule891.pattern, label = 891)

    pattern892 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons5, cons71, cons395, cons576, cons578, cons579, cons575)
    def replacement892(n, p, q, x, b, d, c, a):
        return d*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*(n*(p + q) + S(1))) + Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(-a*d + b*c*(n*(p + q) + S(1))) + d*x**n*(-a*d*(n*(q + S(-1)) + S(1)) + b*c*(n*(p + S(2)*q + S(-1)) + S(1))), x), x)/(b*(n*(p + q) + S(1)))
    rule892 = ReplacementRule(pattern892, replacement892)
    rubi.add(rule892.pattern, label = 892)

    pattern893 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons71, cons402, cons403, cons163, cons575)
    def replacement893(n, p, q, x, b, d, c, a):
        return n*Int((a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(n*(p + q) + S(1)) + x*(a + b*x**n)**p*(c + d*x**n)**q/(n*(p + q) + S(1))
    rule893 = ReplacementRule(pattern893, replacement893)
    rubi.add(rule893.pattern, label = 893)

    pattern894 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons573, cons466, cons580)
    def replacement894(x, b, d, c, a):
        return sqrt(a + b*x**S(2))*EllipticF(ArcTan(x*Rt(d/c, S(2))), S(1) - b*c/(a*d))/(a*sqrt(c*(a + b*x**S(2))/(a*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2)))
    rule894 = ReplacementRule(pattern894, replacement894)
    rubi.add(rule894.pattern, label = 894)

    pattern895 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons581, cons177, cons43, cons582)
    def replacement895(x, b, d, c, a):
        return EllipticF(asin(x*Rt(-d/c, S(2))), b*c/(a*d))/(sqrt(a)*sqrt(c)*Rt(-d/c, S(2)))
    rule895 = ReplacementRule(pattern895, replacement895)
    rubi.add(rule895.pattern, label = 895)

    pattern896 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons581, cons177, cons583)
    def replacement896(x, b, d, c, a):
        return -EllipticF(acos(x*Rt(-d/c, S(2))), b*c/(-a*d + b*c))/(sqrt(c)*sqrt(a - b*c/d)*Rt(-d/c, S(2)))
    rule896 = ReplacementRule(pattern896, replacement896)
    rubi.add(rule896.pattern, label = 896)

    pattern897 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons117)
    def replacement897(x, b, d, c, a):
        return sqrt(S(1) + d*x**S(2)/c)*Int(S(1)/(sqrt(S(1) + d*x**S(2)/c)*sqrt(a + b*x**S(2))), x)/sqrt(c + d*x**S(2))
    rule897 = ReplacementRule(pattern897, replacement897)
    rubi.add(rule897.pattern, label = 897)

    pattern898 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons573, cons466)
    def replacement898(x, b, d, c, a):
        return a*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x) + b*Int(x**S(2)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x)
    rule898 = ReplacementRule(pattern898, replacement898)
    rubi.add(rule898.pattern, label = 898)

    pattern899 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons573, cons483)
    def replacement899(x, b, d, c, a):
        return b*Int(sqrt(c + d*x**S(2))/sqrt(a + b*x**S(2)), x)/d - (-a*d + b*c)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/d
    rule899 = ReplacementRule(pattern899, replacement899)
    rubi.add(rule899.pattern, label = 899)

    pattern900 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons581, cons177, cons43)
    def replacement900(x, b, d, c, a):
        return sqrt(a)*EllipticE(asin(x*Rt(-d/c, S(2))), b*c/(a*d))/(sqrt(c)*Rt(-d/c, S(2)))
    rule900 = ReplacementRule(pattern900, replacement900)
    rubi.add(rule900.pattern, label = 900)

    pattern901 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons581, cons177, cons583)
    def replacement901(x, b, d, c, a):
        return -sqrt(a - b*c/d)*EllipticE(acos(x*Rt(-d/c, S(2))), b*c/(-a*d + b*c))/(sqrt(c)*Rt(-d/c, S(2)))
    rule901 = ReplacementRule(pattern901, replacement901)
    rubi.add(rule901.pattern, label = 901)

    pattern902 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons581, cons177, cons448)
    def replacement902(x, b, d, c, a):
        return sqrt(a + b*x**S(2))*Int(sqrt(S(1) + b*x**S(2)/a)/sqrt(c + d*x**S(2)), x)/sqrt(S(1) + b*x**S(2)/a)
    rule902 = ReplacementRule(pattern902, replacement902)
    rubi.add(rule902.pattern, label = 902)

    pattern903 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons581, cons117)
    def replacement903(x, b, d, c, a):
        return sqrt(S(1) + d*x**S(2)/c)*Int(sqrt(a + b*x**S(2))/sqrt(S(1) + d*x**S(2)/c), x)/sqrt(c + d*x**S(2))
    rule903 = ReplacementRule(pattern903, replacement903)
    rubi.add(rule903.pattern, label = 903)

    pattern904 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons50, cons71, cons128)
    def replacement904(n, p, q, x, b, d, c, a):
        return Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x)
    rule904 = ReplacementRule(pattern904, replacement904)
    rubi.add(rule904.pattern, label = 904)

    pattern905 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons71, cons584, cons43, cons177)
    def replacement905(n, p, q, x, b, d, c, a):
        return a**p*c**q*x*AppellF1(S(1)/n, -p, -q, S(1) + S(1)/n, -b*x**n/a, -d*x**n/c)
    rule905 = ReplacementRule(pattern905, replacement905)
    rubi.add(rule905.pattern, label = 905)

    pattern906 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons71, cons584, cons448)
    def replacement906(n, p, q, x, b, d, c, a):
        return a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((S(1) + b*x**n/a)**p*(c + d*x**n)**q, x)
    rule906 = ReplacementRule(pattern906, replacement906)
    rubi.add(rule906.pattern, label = 906)

    pattern907 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons585, cons586, cons587)
    def replacement907(n, q, p, x, b, d, c, mn, a):
        return Int(x**(-n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x)
    rule907 = ReplacementRule(pattern907, replacement907)
    rubi.add(rule907.pattern, label = 907)

    pattern908 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons585, cons386, cons147)
    def replacement908(n, q, p, x, b, d, c, mn, a):
        return x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(-n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x)
    rule908 = ReplacementRule(pattern908, replacement908)
    rubi.add(rule908.pattern, label = 908)

    pattern909 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons68, cons69)
    def replacement909(n, p, q, x, b, d, c, u, a):
        return Subst(Int((a + b*x**n)**p*(c + d*x**n)**q, x), x, u)/Coefficient(u, x, S(1))
    rule909 = ReplacementRule(pattern909, replacement909)
    rubi.add(rule909.pattern, label = 909)

    pattern910 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1)), x_), cons5, cons50, cons588)
    def replacement910(p, q, v, x, u):
        return Int(NormalizePseudoBinomial(u, x)**p*NormalizePseudoBinomial(v, x)**q, x)
    rule910 = ReplacementRule(pattern910, replacement910)
    rubi.add(rule910.pattern, label = 910)

    pattern911 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*x_**WC('m', S(1)), x_), cons5, cons50, cons589, cons590)
    def replacement911(m, p, q, v, x, u):
        return Int(NormalizePseudoBinomial(v, x)**q*NormalizePseudoBinomial(u*x**(m/p), x)**p, x)
    rule911 = ReplacementRule(pattern911, replacement911)
    rubi.add(rule911.pattern, label = 911)

    pattern912 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons591, cons500)
    def replacement912(m, n, q, p, e, x, b, d, c):
        return b**(S(1) - (m + S(1))/n)*e**m*Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n)*(c + d*x)**q, x), x, x**n)/n
    rule912 = ReplacementRule(pattern912, replacement912)
    rubi.add(rule912.pattern, label = 912)

    pattern913 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons591, cons501)
    def replacement913(m, n, q, p, e, x, b, d, c):
        return b**IntPart(p)*e**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p)*(c + d*x**n)**q, x)
    rule913 = ReplacementRule(pattern913, replacement913)
    rubi.add(rule913.pattern, label = 913)

    pattern914 = Pattern(Integral((e_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons18)
    def replacement914(m, n, q, p, e, x, b, d, c):
        return e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(b*x**n)**p*(c + d*x**n)**q, x)
    rule914 = ReplacementRule(pattern914, replacement914)
    rubi.add(rule914.pattern, label = 914)

    pattern915 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons50, cons71, cons53)
    def replacement915(m, n, p, q, x, b, d, c, a):
        return Subst(Int((a + b*x)**p*(c + d*x)**q, x), x, x**n)/n
    rule915 = ReplacementRule(pattern915, replacement915)
    rubi.add(rule915.pattern, label = 915)

    pattern916 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons71, cons220, cons502)
    def replacement916(m, n, p, q, x, b, d, c, a):
        return Int(x**(m + n*(p + q))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q, x)
    rule916 = ReplacementRule(pattern916, replacement916)
    rubi.add(rule916.pattern, label = 916)

    pattern917 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons50, cons71, cons500)
    def replacement917(m, n, p, q, x, b, d, c, a):
        return Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p*(c + d*x)**q, x), x, x**n)/n
    rule917 = ReplacementRule(pattern917, replacement917)
    rubi.add(rule917.pattern, label = 917)

    pattern918 = Pattern(Integral((e_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons71, cons500)
    def replacement918(m, n, p, q, e, x, b, d, c, a):
        return e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x)
    rule918 = ReplacementRule(pattern918, replacement918)
    rubi.add(rule918.pattern, label = 918)

    pattern919 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons71, cons555)
    def replacement919(m, n, p, q, e, x, b, d, c, a):
        return Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x)
    rule919 = ReplacementRule(pattern919, replacement919)
    rubi.add(rule919.pattern, label = 919)

    pattern920 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons71, cons592, cons66)
    def replacement920(m, n, p, e, x, b, d, c, a):
        return c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1)))
    rule920 = ReplacementRule(pattern920, replacement920)
    rubi.add(rule920.pattern, label = 920)

    pattern921 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons57, cons58, cons59, cons60, cons7, cons27, cons48, cons21, cons4, cons5, cons593, cons55, cons594, cons66)
    def replacement921(m, n, a1, p, non2, b2, a2, e, b1, x, d, c):
        return c*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(a1*a2*e*(m + S(1)))
    rule921 = ReplacementRule(pattern921, replacement921)
    rubi.add(rule921.pattern, label = 921)

    pattern922 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons71, cons595, cons596, cons93, cons597)
    def replacement922(m, n, p, e, x, b, d, c, a):
        return d*e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p, x) + c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1)))
    rule922 = ReplacementRule(pattern922, replacement922)
    rubi.add(rule922.pattern, label = 922)

    pattern923 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons71, cons595, cons66)
    def replacement923(m, n, p, e, x, b, d, c, a):
        return d*Int((e*x)**m*(a + b*x**n)**(p + S(1)), x)/b + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*e*(m + S(1)))
    rule923 = ReplacementRule(pattern923, replacement923)
    rubi.add(rule923.pattern, label = 923)

    pattern924 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons71, cons596, cons93, cons597, cons598)
    def replacement924(m, n, p, e, x, b, d, c, a):
        return c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))) + e**(-n)*(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1)))
    rule924 = ReplacementRule(pattern924, replacement924)
    rubi.add(rule924.pattern, label = 924)

    pattern925 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons57, cons58, cons59, cons60, cons7, cons27, cons48, cons5, cons593, cons55, cons596, cons93, cons597, cons598)
    def replacement925(m, n, a1, p, non2, b2, a2, e, b1, x, d, c):
        return c*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(a1*a2*e*(m + S(1))) + e**(-n)*(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**(m + n)*(a1 + b1*x**(n/S(2)))**p*(a2 + b2*x**(n/S(2)))**p, x)/(a1*a2*(m + S(1)))
    rule925 = ReplacementRule(pattern925, replacement925)
    rubi.add(rule925.pattern, label = 925)

    pattern926 = Pattern(Integral(x_**m_*(a_ + x_**S(2)*WC('b', S(1)))**p_*(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons71, cons13, cons137, cons599, cons600)
    def replacement926(m, p, x, b, d, c, a):
        return b**(-m/S(2) + S(-1))*x*(-a)**(m/S(2) + S(-1))*(a + b*x**S(2))**(p + S(1))*(-a*d + b*c)/(S(2)*(p + S(1))) + b**(-m/S(2) + S(-1))*Int((a + b*x**S(2))**(p + S(1))*ExpandToSum(S(2)*b*x**S(2)*(p + S(1))*(b**(m/S(2))*x**(m + S(-2))*(c + d*x**S(2)) - (-a)**(m/S(2) + S(-1))*(-a*d + b*c))/(a + b*x**S(2)) - (-a)**(m/S(2) + S(-1))*(-a*d + b*c), x), x)/(S(2)*(p + S(1)))
    rule926 = ReplacementRule(pattern926, replacement926)
    rubi.add(rule926.pattern, label = 926)

    pattern927 = Pattern(Integral(x_**m_*(a_ + x_**S(2)*WC('b', S(1)))**p_*(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons71, cons13, cons137, cons601, cons600)
    def replacement927(m, p, x, b, d, c, a):
        return b**(-m/S(2) + S(-1))*x*(-a)**(m/S(2) + S(-1))*(a + b*x**S(2))**(p + S(1))*(-a*d + b*c)/(S(2)*(p + S(1))) + b**(-m/S(2) + S(-1))*Int(x**m*(a + b*x**S(2))**(p + S(1))*ExpandToSum(S(2)*b*(p + S(1))*(b**(m/S(2))*(c + d*x**S(2)) - x**(-m + S(2))*(-a)**(m/S(2) + S(-1))*(-a*d + b*c))/(a + b*x**S(2)) - x**(-m)*(-a)**(m/S(2) + S(-1))*(-a*d + b*c), x), x)/(S(2)*(p + S(1)))
    rule927 = ReplacementRule(pattern927, replacement927)
    rubi.add(rule927.pattern, label = 927)

    pattern928 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons71, cons13, cons137, cons602)
    def replacement928(m, n, p, e, x, b, d, c, a):
        return -(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*b*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*e*n*(p + S(1)))
    rule928 = ReplacementRule(pattern928, replacement928)
    rubi.add(rule928.pattern, label = 928)

    pattern929 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons57, cons58, cons59, cons60, cons7, cons27, cons48, cons21, cons4, cons593, cons55, cons13, cons137, cons602)
    def replacement929(m, n, a1, p, non2, b2, a2, e, b1, x, d, c):
        return -(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1)), x)/(a1*a2*b1*b2*n*(p + S(1))) - (e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))*(-a1*a2*d + b1*b2*c)/(a1*a2*b1*b2*e*n*(p + S(1)))
    rule929 = ReplacementRule(pattern929, replacement929)
    rubi.add(rule929.pattern, label = 929)

    pattern930 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons71, cons603)
    def replacement930(m, n, p, e, x, b, d, c, a):
        return d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(b*e*(m + n*(p + S(1)) + S(1))) - (a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a + b*x**n)**p, x)/(b*(m + n*(p + S(1)) + S(1)))
    rule930 = ReplacementRule(pattern930, replacement930)
    rubi.add(rule930.pattern, label = 930)

    pattern931 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons57, cons58, cons59, cons60, cons7, cons27, cons48, cons21, cons4, cons5, cons593, cons55, cons603)
    def replacement931(m, n, a1, p, non2, b2, a2, e, b1, x, d, c):
        return d*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(b1*b2*e*(m + n*(p + S(1)) + S(1))) - (a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a1 + b1*x**(n/S(2)))**p*(a2 + b2*x**(n/S(2)))**p, x)/(b1*b2*(m + n*(p + S(1)) + S(1)))
    rule931 = ReplacementRule(pattern931, replacement931)
    rubi.add(rule931.pattern, label = 931)

    pattern932 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons71, cons148, cons128, cons604)
    def replacement932(m, n, p, x, e, b, d, c, a):
        return Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p/(c + d*x**n), x), x)
    rule932 = ReplacementRule(pattern932, replacement932)
    rubi.add(rule932.pattern, label = 932)

    pattern933 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons71, cons148, cons93, cons94, cons88)
    def replacement933(m, n, p, x, e, b, d, c, a):
        return c**S(2)*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))) - e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*Simp(-a*d**S(2)*x**n*(m + S(1)) + b*c**S(2)*n*(p + S(1)) + c*(m + S(1))*(-S(2)*a*d + b*c), x), x)/(a*(m + S(1)))
    rule933 = ReplacementRule(pattern933, replacement933)
    rubi.add(rule933.pattern, label = 933)

    pattern934 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons71, cons148, cons13, cons137)
    def replacement934(m, n, p, x, e, b, d, c, a):
        return Int((e*x)**m*(a + b*x**n)**(p + S(1))*Simp(a*b*d**S(2)*n*x**n*(p + S(1)) + b**S(2)*c**S(2)*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)**S(2), x), x)/(a*b**S(2)*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)**S(2)/(a*b**S(2)*e*n*(p + S(1)))
    rule934 = ReplacementRule(pattern934, replacement934)
    rubi.add(rule934.pattern, label = 934)

    pattern935 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons71, cons148, cons605)
    def replacement935(m, n, p, x, e, b, d, c, a):
        return d**S(2)*e**(-n + S(-1))*(e*x)**(m + n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*(p + S(2)) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*Simp(b*c**S(2)*(m + n*(p + S(2)) + S(1)) + d*x**n*(S(2)*b*c*n*(p + S(1)) + (-a*d + S(2)*b*c)*(m + n + S(1))), x), x)/(b*(m + n*(p + S(2)) + S(1)))
    rule935 = ReplacementRule(pattern935, replacement935)
    rubi.add(rule935.pattern, label = 935)

    def With936(m, n, p, q, x, b, d, c, a):
        k = GCD(m + S(1), n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p*(c + d*x**(n/k))**q, x), x, x**k)/k
        return False
    pattern936 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons5, cons50, cons71, cons148, cons17, CustomConstraint(With936))
    rule936 = ReplacementRule(pattern936, With936)
    rubi.add(rule936.pattern, label = 936)

    def With937(m, n, p, q, x, e, b, d, c, a):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*e**(-n)*x**(k*n))**p*(c + d*e**(-n)*x**(k*n))**q, x), x, (e*x)**(S(1)/k))/e
    pattern937 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons50, cons71, cons148, cons367, cons38 )
    rule937 = ReplacementRule(pattern937, With937)
    rubi.add(rule937.pattern, label = 937)

    pattern938 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons71, cons148, cons606, cons137, cons403, cons607, cons608)
    def replacement938(m, n, p, q, x, e, b, d, c, a):
        return -e**n*Int((e*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m - n + S(1)) + d*x**n*(m + n*(q + S(-1)) + S(1)), x), x)/(b*n*(p + S(1))) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*n*(p + S(1)))
    rule938 = ReplacementRule(pattern938, replacement938)
    rubi.add(rule938.pattern, label = 938)

    pattern939 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons71, cons148, cons402, cons137, cons576, cons608)
    def replacement939(m, n, p, q, x, e, b, d, c, a):
        return Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)) + d*x**n*(b*c*n*(p + S(1)) + (-a*d + b*c)*(m + n*(q + S(-1)) + S(1))), x), x)/(a*b*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(-a*d + b*c)/(a*b*e*n*(p + S(1)))
    rule939 = ReplacementRule(pattern939, replacement939)
    rubi.add(rule939.pattern, label = 939)

    pattern940 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons71, cons148, cons402, cons137, cons574, cons608)
    def replacement940(m, n, p, q, x, e, b, d, c, a):
        return Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m + n*(p + S(1)) + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*n*(p + S(1)))
    rule940 = ReplacementRule(pattern940, replacement940)
    rubi.add(rule940.pattern, label = 940)

    pattern941 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons50, cons71, cons148, cons244, cons137, cons609, cons608)
    def replacement941(m, n, p, q, x, e, b, d, c, a):
        return -a*e**(S(2)*n + S(-1))*(e*x)**(m - S(2)*n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*n*(p + S(1))*(-a*d + b*c)) + e**(S(2)*n)*Int((e*x)**(m - S(2)*n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(a*c*(m - S(2)*n + S(1)) + x**n*(a*d*(m + n*q - n + S(1)) + b*c*n*(p + S(1))), x), x)/(b*n*(p + S(1))*(-a*d + b*c))
    rule941 = ReplacementRule(pattern941, replacement941)
    rubi.add(rule941.pattern, label = 941)

    pattern942 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons50, cons71, cons148, cons244, cons137, cons610, cons608)
    def replacement942(m, n, p, q, x, e, b, d, c, a):
        return -e**n*Int((e*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m - n + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(n*(p + S(1))*(-a*d + b*c)) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(n*(p + S(1))*(-a*d + b*c))
    rule942 = ReplacementRule(pattern942, replacement942)
    rubi.add(rule942.pattern, label = 942)

    pattern943 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons50, cons71, cons148, cons13, cons137, cons608)
    def replacement943(m, n, p, q, x, e, b, d, c, a):
        return -b*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*e*n*(p + S(1))*(-a*d + b*c)) + Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c*(m + S(1)) + b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c))
    rule943 = ReplacementRule(pattern943, replacement943)
    rubi.add(rule943.pattern, label = 943)

    pattern944 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons71, cons148, cons606, cons403, cons94, cons163, cons608)
    def replacement944(m, n, p, q, x, e, b, d, c, a):
        return -e**(-n)*n*Int((e*x)**(m + n)*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*d*q + b*c*p + b*d*x**n*(p + q), x), x)/(m + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + S(1)))
    rule944 = ReplacementRule(pattern944, replacement944)
    rubi.add(rule944.pattern, label = 944)

    pattern945 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons71, cons148, cons611, cons576, cons94, cons608)
    def replacement945(m, n, p, q, x, e, b, d, c, a):
        return c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(a*e*(m + S(1))) - e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*n*(a*d*(q + S(-1)) + b*c*(p + S(1))) + c*(m + S(1))*(-a*d + b*c) + d*x**n*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)), x), x)/(a*(m + S(1)))
    rule945 = ReplacementRule(pattern945, replacement945)
    rubi.add(rule945.pattern, label = 945)

    pattern946 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons71, cons148, cons611, cons574, cons94, cons608)
    def replacement946(m, n, p, q, x, e, b, d, c, a):
        return -e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(b*c*(m + S(1)) + d*x**n*(b*n*(p + q + S(1)) + b*(m + S(1))) + n*(a*d*q + b*c*(p + S(1))), x), x)/(a*(m + S(1))) + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*(m + S(1)))
    rule946 = ReplacementRule(pattern946, replacement946)
    rubi.add(rule946.pattern, label = 946)

    pattern947 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons71, cons148, cons402, cons403, cons163, cons608)
    def replacement947(m, n, p, q, x, e, b, d, c, a):
        return n*Int((e*x)**m*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(m + n*(p + q) + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + n*(p + q) + S(1)))
    rule947 = ReplacementRule(pattern947, replacement947)
    rubi.add(rule947.pattern, label = 947)

    pattern948 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons71, cons148, cons395, cons576, cons608)
    def replacement948(m, n, p, q, x, e, b, d, c, a):
        return d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*e*(m + n*(p + q) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)) + x**n*(b*c*d*n*(p + q) + d*n*(q + S(-1))*(-a*d + b*c) + d*(m + S(1))*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1)))
    rule948 = ReplacementRule(pattern948, replacement948)
    rubi.add(rule948.pattern, label = 948)

    pattern949 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons71, cons148, cons611, cons403, cons607, cons608)
    def replacement949(m, n, p, q, x, e, b, d, c, a):
        return -e**n*Int((e*x)**(m - n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(a*c*(m - n + S(1)) + x**n*(a*d*(m - n + S(1)) - n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1))) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*(m + n*(p + q) + S(1)))
    rule949 = ReplacementRule(pattern949, replacement949)
    rubi.add(rule949.pattern, label = 949)

    pattern950 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons50, cons71, cons148, cons31, cons609, cons608)
    def replacement950(m, n, p, q, x, e, b, d, c, a):
        return -e**(S(2)*n)*Int((e*x)**(m - S(2)*n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*(m - S(2)*n + S(1)) + x**n*(a*d*(m + n*(q + S(-1)) + S(1)) + b*c*(m + n*(p + S(-1)) + S(1))), x), x)/(b*d*(m + n*(p + q) + S(1))) + e**(S(2)*n + S(-1))*(e*x)**(m - S(2)*n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*d*(m + n*(p + q) + S(1)))
    rule950 = ReplacementRule(pattern950, replacement950)
    rubi.add(rule950.pattern, label = 950)

    pattern951 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons50, cons71, cons148, cons31, cons94, cons608)
    def replacement951(m, n, p, q, x, e, b, d, c, a):
        return -e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(a*d*q + b*c*p) + (a*d + b*c)*(m + n + S(1)), x), x)/(a*c*(m + S(1))) + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c*e*(m + S(1)))
    rule951 = ReplacementRule(pattern951, replacement951)
    rubi.add(rule951.pattern, label = 951)

    pattern952 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons71, cons148, cons31, cons612)
    def replacement952(m, n, x, e, b, d, c, a):
        return -a*e**n*Int((e*x)**(m - n)/(a + b*x**n), x)/(-a*d + b*c) + c*e**n*Int((e*x)**(m - n)/(c + d*x**n), x)/(-a*d + b*c)
    rule952 = ReplacementRule(pattern952, replacement952)
    rubi.add(rule952.pattern, label = 952)

    pattern953 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons71, cons148)
    def replacement953(m, n, x, e, b, d, c, a):
        return b*Int((e*x)**m/(a + b*x**n), x)/(-a*d + b*c) - d*Int((e*x)**m/(c + d*x**n), x)/(-a*d + b*c)
    rule953 = ReplacementRule(pattern953, replacement953)
    rubi.add(rule953.pattern, label = 953)

    pattern954 = Pattern(Integral(x_**m_/((a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71, cons613, cons614, cons615)
    def replacement954(m, n, x, b, d, c, a):
        return -a*Int(x**(m - n)/((a + b*x**n)*sqrt(c + d*x**n)), x)/b + Int(x**(m - n)/sqrt(c + d*x**n), x)/b
    rule954 = ReplacementRule(pattern954, replacement954)
    rubi.add(rule954.pattern, label = 954)

    def With955(x, b, d, c, a):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return -s*Int(S(1)/(sqrt(c + d*x**S(4))*(r - s*x**S(2))), x)/(S(2)*b) + s*Int(S(1)/(sqrt(c + d*x**S(4))*(r + s*x**S(2))), x)/(S(2)*b)
    pattern955 = Pattern(Integral(x_**S(2)/((a_ + x_**S(4)*WC('b', S(1)))*sqrt(c_ + x_**S(4)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71 )
    rule955 = ReplacementRule(pattern955, With955)
    rubi.add(rule955.pattern, label = 955)

    def With956(x, b, d, c, a):
        q = Rt(d/c, S(3))
        return -S(2)**(S(1)/3)*sqrt(S(3))*q*ArcTan(sqrt(S(3))/S(3) + S(2)**(S(2)/3)*sqrt(S(3))*(sqrt(c) - sqrt(c + d*x**S(3)))/(S(3)*sqrt(c)*q*x))/(S(18)*b*sqrt(c)) + S(2)**(S(1)/3)*sqrt(S(3))*q*ArcTan(sqrt(S(3))/S(3) + S(2)**(S(2)/3)*sqrt(S(3))*(sqrt(c) + sqrt(c + d*x**S(3)))/(S(3)*sqrt(c)*q*x))/(S(18)*b*sqrt(c)) + S(2)**(S(1)/3)*q*log(-S(2)**(S(1)/3)*q*x + S(1) - sqrt(c + d*x**S(3))/sqrt(c))/(S(12)*b*sqrt(c)) - S(2)**(S(1)/3)*q*log(-S(2)**(S(1)/3)*q*x + S(1) + sqrt(c + d*x**S(3))/sqrt(c))/(S(12)*b*sqrt(c)) + S(2)**(S(1)/3)*q*atanh(sqrt(c + d*x**S(3))/sqrt(c))/(S(18)*b*sqrt(c))
    pattern956 = Pattern(Integral(x_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71, cons616 )
    rule956 = ReplacementRule(pattern956, With956)
    rubi.add(rule956.pattern, label = 956)

    pattern957 = Pattern(Integral(x_**m_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71, cons616, cons617)
    def replacement957(m, x, b, d, c, a):
        return -a*Int(x**(m + S(-3))/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x)/b + Int(x**(m + S(-3))/sqrt(c + d*x**S(3)), x)/b
    rule957 = ReplacementRule(pattern957, replacement957)
    rubi.add(rule957.pattern, label = 957)

    pattern958 = Pattern(Integral(x_**m_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71, cons616, cons618)
    def replacement958(m, x, b, d, c, a):
        return -b*Int(x**(m + S(3))/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x)/a + Int(x**m/sqrt(c + d*x**S(3)), x)/a
    rule958 = ReplacementRule(pattern958, replacement958)
    rubi.add(rule958.pattern, label = 958)

    pattern959 = Pattern(Integral(x_**S(2)*sqrt(c_ + x_**S(4)*WC('d', S(1)))/(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons71)
    def replacement959(x, b, d, c, a):
        return d*Int(x**S(2)/sqrt(c + d*x**S(4)), x)/b + (-a*d + b*c)*Int(x**S(2)/((a + b*x**S(4))*sqrt(c + d*x**S(4))), x)/b
    rule959 = ReplacementRule(pattern959, replacement959)
    rubi.add(rule959.pattern, label = 959)

    pattern960 = Pattern(Integral(x_**WC('m', S(1))*sqrt(c_ + x_**S(3)*WC('d', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons71, cons616, cons619)
    def replacement960(m, x, b, d, c, a):
        return d*Int(x**m/sqrt(c + d*x**S(3)), x)/b + (-a*d + b*c)*Int(x**m/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x)/b
    rule960 = ReplacementRule(pattern960, replacement960)
    rubi.add(rule960.pattern, label = 960)

    pattern961 = Pattern(Integral(x_**S(2)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71, cons466, cons573, cons580)
    def replacement961(x, b, d, c, a):
        return -c*Int(sqrt(a + b*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/b + x*sqrt(a + b*x**S(2))/(b*sqrt(c + d*x**S(2)))
    rule961 = ReplacementRule(pattern961, replacement961)
    rubi.add(rule961.pattern, label = 961)

    pattern962 = Pattern(Integral(x_**n_/(sqrt(a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons71, cons620, cons621)
    def replacement962(n, x, b, d, c, a):
        return -a*Int(S(1)/(sqrt(a + b*x**n)*sqrt(c + d*x**n)), x)/b + Int(sqrt(a + b*x**n)/sqrt(c + d*x**n), x)/b
    rule962 = ReplacementRule(pattern962, replacement962)
    rubi.add(rule962.pattern, label = 962)

    def With963(m, n, q, p, x, b, d, c, a):
        k = Denominator(p)
        return a**(p + (m + S(1))/n)*k*Subst(Int(x**(k*(m + S(1))/n + S(-1))*(c - x**k*(-a*d + b*c))**q*(-b*x**k + S(1))**(-p - q + S(-1) - (m + S(1))/n), x), x, x**(n/k)*(a + b*x**n)**(-S(1)/k))/n
    pattern963 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons148, cons244, cons622, cons485 )
    rule963 = ReplacementRule(pattern963, With963)
    rubi.add(rule963.pattern, label = 963)

    pattern964 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons5, cons50, cons71, cons196, cons17)
    def replacement964(m, n, p, q, x, b, d, c, a):
        return -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q, x), x, S(1)/x)
    rule964 = ReplacementRule(pattern964, replacement964)
    rubi.add(rule964.pattern, label = 964)

    def With965(m, n, p, q, x, e, b, d, c, a):
        g = Denominator(m)
        return -g*Subst(Int(x**(-g*(m + S(1)) + S(-1))*(a + b*e**(-n)*x**(-g*n))**p*(c + d*e**(-n)*x**(-g*n))**q, x), x, (e*x)**(-S(1)/g))/e
    pattern965 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons50, cons196, cons367 )
    rule965 = ReplacementRule(pattern965, With965)
    rubi.add(rule965.pattern, label = 965)

    pattern966 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons50, cons71, cons196, cons356)
    def replacement966(m, n, p, q, x, e, b, d, c, a):
        return -(e*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q, x), x, S(1)/x)
    rule966 = ReplacementRule(pattern966, replacement966)
    rubi.add(rule966.pattern, label = 966)

    def With967(m, n, p, q, x, b, d, c, a):
        g = Denominator(n)
        return g*Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + b*x**(g*n))**p*(c + d*x**(g*n))**q, x), x, x**(S(1)/g))
    pattern967 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons21, cons5, cons50, cons71, cons489 )
    rule967 = ReplacementRule(pattern967, With967)
    rubi.add(rule967.pattern, label = 967)

    pattern968 = Pattern(Integral((e_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons50, cons71, cons489)
    def replacement968(m, n, p, q, x, e, b, d, c, a):
        return e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x)
    rule968 = ReplacementRule(pattern968, replacement968)
    rubi.add(rule968.pattern, label = 968)

    pattern969 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons50, cons71, cons541, cons23)
    def replacement969(m, n, p, q, x, b, d, c, a):
        return Subst(Int((a + b*x**(n/(m + S(1))))**p*(c + d*x**(n/(m + S(1))))**q, x), x, x**(m + S(1)))/(m + S(1))
    rule969 = ReplacementRule(pattern969, replacement969)
    rubi.add(rule969.pattern, label = 969)

    pattern970 = Pattern(Integral((e_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons71, cons541, cons23)
    def replacement970(m, n, p, q, x, e, b, d, c, a):
        return e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x)
    rule970 = ReplacementRule(pattern970, replacement970)
    rubi.add(rule970.pattern, label = 970)

    pattern971 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons71, cons402, cons137, cons576, cons608)
    def replacement971(m, n, p, q, x, e, b, d, c, a):
        return Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)) + d*x**n*(b*c*n*(p + S(1)) + (-a*d + b*c)*(m + n*(q + S(-1)) + S(1))), x), x)/(a*b*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(-a*d + b*c)/(a*b*e*n*(p + S(1)))
    rule971 = ReplacementRule(pattern971, replacement971)
    rubi.add(rule971.pattern, label = 971)

    pattern972 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons71, cons402, cons137, cons574, cons608)
    def replacement972(m, n, p, q, x, e, b, d, c, a):
        return Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m + n*(p + S(1)) + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*n*(p + S(1)))
    rule972 = ReplacementRule(pattern972, replacement972)
    rubi.add(rule972.pattern, label = 972)

    pattern973 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons50, cons71, cons13, cons137, cons608)
    def replacement973(m, n, p, q, x, e, b, d, c, a):
        return -b*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*e*n*(p + S(1))*(-a*d + b*c)) + Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c*(m + S(1)) + b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c))
    rule973 = ReplacementRule(pattern973, replacement973)
    rubi.add(rule973.pattern, label = 973)

    pattern974 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons71, cons402, cons403, cons163, cons608)
    def replacement974(m, n, p, q, x, e, b, d, c, a):
        return n*Int((e*x)**m*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(m + n*(p + q) + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + n*(p + q) + S(1)))
    rule974 = ReplacementRule(pattern974, replacement974)
    rubi.add(rule974.pattern, label = 974)

    pattern975 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons71, cons395, cons576, cons608)
    def replacement975(m, n, p, q, x, e, b, d, c, a):
        return d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*e*(m + n*(p + q) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)) + x**n*(b*c*d*n*(p + q) + d*n*(q + S(-1))*(-a*d + b*c) + d*(m + S(1))*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1)))
    rule975 = ReplacementRule(pattern975, replacement975)
    rubi.add(rule975.pattern, label = 975)

    pattern976 = Pattern(Integral(x_**m_/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons71, cons623)
    def replacement976(m, n, x, b, d, c, a):
        return -a*Int(x**(m - n)/(a + b*x**n), x)/(-a*d + b*c) + c*Int(x**(m - n)/(c + d*x**n), x)/(-a*d + b*c)
    rule976 = ReplacementRule(pattern976, replacement976)
    rubi.add(rule976.pattern, label = 976)

    pattern977 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons21, cons71)
    def replacement977(m, n, x, e, b, d, c, a):
        return b*Int((e*x)**m/(a + b*x**n), x)/(-a*d + b*c) - d*Int((e*x)**m/(c + d*x**n), x)/(-a*d + b*c)
    rule977 = ReplacementRule(pattern977, replacement977)
    rubi.add(rule977.pattern, label = 977)

    pattern978 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons71, cons624, cons625, cons626)
    def replacement978(m, n, p, q, x, e, b, d, c, a):
        return Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x)
    rule978 = ReplacementRule(pattern978, replacement978)
    rubi.add(rule978.pattern, label = 978)

    pattern979 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons585, cons586, cons587)
    def replacement979(m, n, q, p, x, b, d, c, mn, a):
        return Int(x**(m - n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x)
    rule979 = ReplacementRule(pattern979, replacement979)
    rubi.add(rule979.pattern, label = 979)

    pattern980 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons50, cons585, cons386, cons147)
    def replacement980(m, n, p, q, x, b, d, c, mn, a):
        return x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(m - n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x)
    rule980 = ReplacementRule(pattern980, replacement980)
    rubi.add(rule980.pattern, label = 980)

    pattern981 = Pattern(Integral((e_*x_)**m_*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons585)
    def replacement981(m, n, q, p, e, x, b, d, c, mn, a):
        return e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**(-n))**q, x)
    rule981 = ReplacementRule(pattern981, replacement981)
    rubi.add(rule981.pattern, label = 981)

    pattern982 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons71, cons66, cons627, cons43, cons177)
    def replacement982(m, n, p, q, x, e, b, d, c, a):
        return a**p*c**q*(e*x)**(m + S(1))*AppellF1((m + S(1))/n, -p, -q, S(1) + (m + S(1))/n, -b*x**n/a, -d*x**n/c)/(e*(m + S(1)))
    rule982 = ReplacementRule(pattern982, replacement982)
    rubi.add(rule982.pattern, label = 982)

    pattern983 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons71, cons66, cons627, cons448)
    def replacement983(m, n, p, q, x, e, b, d, c, a):
        return a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((e*x)**m*(S(1) + b*x**n/a)**p*(c + d*x**n)**q, x)
    rule983 = ReplacementRule(pattern983, replacement983)
    rubi.add(rule983.pattern, label = 983)

    pattern984 = Pattern(Integral(x_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons552, cons17, cons553)
    def replacement984(m, n, p, q, v, x, b, d, c, a):
        return Coefficient(v, x, S(1))**(-m + S(-1))*Subst(Int(SimplifyIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(x - Coefficient(v, x, S(0)))**m, x), x), x, v)
    rule984 = ReplacementRule(pattern984, replacement984)
    rubi.add(rule984.pattern, label = 984)

    pattern985 = Pattern(Integral(u_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons50, cons554)
    def replacement985(m, n, p, q, v, x, b, d, c, u, a):
        return u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x, v)/Coefficient(v, x, S(1))
    rule985 = ReplacementRule(pattern985, replacement985)
    rubi.add(rule985.pattern, label = 985)

    pattern986 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons57, cons58, cons59, cons60, cons7, cons27, cons4, cons5, cons50, cons593, cons55, cons494)
    def replacement986(n, a1, q, non2, p, b2, a2, b1, x, d, c, u):
        return Int(u*(c + d*x**n)**q*(a1*a2 + b1*b2*x**n)**p, x)
    rule986 = ReplacementRule(pattern986, replacement986)
    rubi.add(rule986.pattern, label = 986)

    pattern987 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)) + x_**WC('n2', S(1))*WC('e', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons57, cons58, cons59, cons60, cons7, cons27, cons48, cons4, cons5, cons50, cons593, cons46, cons55, cons494)
    def replacement987(n, a1, n2, q, non2, p, b2, a2, e, b1, x, d, c, u):
        return Int(u*(a1*a2 + b1*b2*x**n)**p*(c + d*x**n + e*x**(S(2)*n))**q, x)
    rule987 = ReplacementRule(pattern987, replacement987)
    rubi.add(rule987.pattern, label = 987)

    pattern988 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**p_*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons57, cons58, cons59, cons60, cons7, cons27, cons4, cons5, cons50, cons593, cons55)
    def replacement988(n, a1, q, non2, p, b2, a2, b1, x, d, c, u):
        return (a1 + b1*x**(n/S(2)))**FracPart(p)*(a2 + b2*x**(n/S(2)))**FracPart(p)*(a1*a2 + b1*b2*x**n)**(-FracPart(p))*Int(u*(c + d*x**n)**q*(a1*a2 + b1*b2*x**n)**p, x)
    rule988 = ReplacementRule(pattern988, replacement988)
    rubi.add(rule988.pattern, label = 988)

    pattern989 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)) + x_**WC('n2', S(1))*WC('e', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons57, cons58, cons59, cons60, cons7, cons27, cons48, cons4, cons5, cons50, cons593, cons46, cons55)
    def replacement989(n, a1, n2, q, non2, p, b2, a2, e, b1, x, d, c, u):
        return (a1 + b1*x**(n/S(2)))**FracPart(p)*(a2 + b2*x**(n/S(2)))**FracPart(p)*(a1*a2 + b1*b2*x**n)**(-FracPart(p))*Int(u*(a1*a2 + b1*b2*x**n)**p*(c + d*x**n + e*x**(S(2)*n))**q, x)
    rule989 = ReplacementRule(pattern989, replacement989)
    rubi.add(rule989.pattern, label = 989)

    pattern990 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons628)
    def replacement990(n, p, q, x, e, b, d, f, r, c, a):
        return Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x)
    rule990 = ReplacementRule(pattern990, replacement990)
    rubi.add(rule990.pattern, label = 990)

    pattern991 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons629)
    def replacement991(n, x, e, b, d, f, c, a):
        return (-a*f + b*e)*Int(S(1)/(a + b*x**n), x)/(-a*d + b*c) - (-c*f + d*e)*Int(S(1)/(c + d*x**n), x)/(-a*d + b*c)
    rule991 = ReplacementRule(pattern991, replacement991)
    rubi.add(rule991.pattern, label = 991)

    pattern992 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/((a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons629)
    def replacement992(n, x, e, b, d, f, c, a):
        return f*Int(S(1)/sqrt(c + d*x**n), x)/b + (-a*f + b*e)*Int(S(1)/((a + b*x**n)*sqrt(c + d*x**n)), x)/b
    rule992 = ReplacementRule(pattern992, replacement992)
    rubi.add(rule992.pattern, label = 992)

    pattern993 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/(sqrt(a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons630)
    def replacement993(n, x, e, b, d, f, c, a):
        return f*Int(sqrt(a + b*x**n)/sqrt(c + d*x**n), x)/b + (-a*f + b*e)*Int(S(1)/(sqrt(a + b*x**n)*sqrt(c + d*x**n)), x)/b
    rule993 = ReplacementRule(pattern993, replacement993)
    rubi.add(rule993.pattern, label = 993)

    pattern994 = Pattern(Integral((e_ + x_**S(2)*WC('f', S(1)))/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons466, cons573)
    def replacement994(x, e, b, d, f, c, a):
        return (-a*f + b*e)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - (-c*f + d*e)*Int(sqrt(a + b*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c)
    rule994 = ReplacementRule(pattern994, replacement994)
    rubi.add(rule994.pattern, label = 994)

    pattern995 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons402, cons137, cons403)
    def replacement995(n, q, p, x, e, b, d, f, c, a):
        return -x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*n*(p + S(1))) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(-a*f + b*e*n*(p + S(1)) + b*e) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(n*q + S(1))), x), x)/(a*b*n*(p + S(1)))
    rule995 = ReplacementRule(pattern995, replacement995)
    rubi.add(rule995.pattern, label = 995)

    pattern996 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons50, cons13, cons137)
    def replacement996(n, q, p, x, e, b, d, f, c, a):
        return -x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*n*(p + S(1))*(-a*d + b*c)) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c))
    rule996 = ReplacementRule(pattern996, replacement996)
    rubi.add(rule996.pattern, label = 996)

    pattern997 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons395, cons403, cons631)
    def replacement997(n, p, q, x, e, b, d, f, c, a):
        return f*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*(n*(p + q + S(1)) + S(1))) + Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(-a*f + b*e*n*(p + q + S(1)) + b*e) + x**n*(b*d*e*n*(p + q + S(1)) + d*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(n*(p + q + S(1)) + S(1)))
    rule997 = ReplacementRule(pattern997, replacement997)
    rubi.add(rule997.pattern, label = 997)

    pattern998 = Pattern(Integral((e_ + x_**S(4)*WC('f', S(1)))/((a_ + x_**S(4)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(4)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons153)
    def replacement998(x, e, b, d, f, c, a):
        return (-a*f + b*e)*Int((a + b*x**S(4))**(S(-3)/4), x)/(-a*d + b*c) - (-c*f + d*e)*Int((a + b*x**S(4))**(S(1)/4)/(c + d*x**S(4)), x)/(-a*d + b*c)
    rule998 = ReplacementRule(pattern998, replacement998)
    rubi.add(rule998.pattern, label = 998)

    pattern999 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons4, cons632)
    def replacement999(n, p, x, e, b, d, f, c, a):
        return f*Int((a + b*x**n)**p, x)/d + (-c*f + d*e)*Int((a + b*x**n)**p/(c + d*x**n), x)/d
    rule999 = ReplacementRule(pattern999, replacement999)
    rubi.add(rule999.pattern, label = 999)

    pattern1000 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons50, cons633)
    def replacement1000(n, p, q, x, e, b, d, f, c, a):
        return e*Int((a + b*x**n)**p*(c + d*x**n)**q, x) + f*Int(x**n*(a + b*x**n)**p*(c + d*x**n)**q, x)
    rule1000 = ReplacementRule(pattern1000, replacement1000)
    rubi.add(rule1000.pattern, label = 1000)

    pattern1001 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons153)
    def replacement1001(x, e, b, d, f, c, a):
        return b*Int(S(1)/((a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/(-a*d + b*c) - d*Int(S(1)/((c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(-a*d + b*c)
    rule1001 = ReplacementRule(pattern1001, replacement1001)
    rubi.add(rule1001.pattern, label = 1001)

    pattern1002 = Pattern(Integral(S(1)/(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons7, cons27, cons48, cons125, cons176)
    def replacement1002(x, e, d, f, c):
        return -d*Int(S(1)/((c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/c + Int(S(1)/(x**S(2)*sqrt(e + f*x**S(2))), x)/c
    rule1002 = ReplacementRule(pattern1002, replacement1002)
    rubi.add(rule1002.pattern, label = 1002)

    pattern1003 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons634, cons635, cons636)
    def replacement1003(x, e, b, d, f, c, a):
        return d*Int(sqrt(e + f*x**S(2))/sqrt(c + d*x**S(2)), x)/b + (-a*d + b*c)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/b
    rule1003 = ReplacementRule(pattern1003, replacement1003)
    rubi.add(rule1003.pattern, label = 1003)

    pattern1004 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons637)
    def replacement1004(x, e, b, d, f, c, a):
        return d*Int(sqrt(e + f*x**S(2))/sqrt(c + d*x**S(2)), x)/b + (-a*d + b*c)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/b
    rule1004 = ReplacementRule(pattern1004, replacement1004)
    rubi.add(rule1004.pattern, label = 1004)

    pattern1005 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons573, cons638, cons636)
    def replacement1005(x, e, b, d, f, c, a):
        return b*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*f + b*e) - f*Int(S(1)/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(-a*f + b*e)
    rule1005 = ReplacementRule(pattern1005, replacement1005)
    rubi.add(rule1005.pattern, label = 1005)

    pattern1006 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons581, cons177, cons178, cons639)
    def replacement1006(x, e, b, d, f, c, a):
        return EllipticPi(b*c/(a*d), asin(x*Rt(-d/c, S(2))), c*f/(d*e))/(a*sqrt(c)*sqrt(e)*Rt(-d/c, S(2)))
    rule1006 = ReplacementRule(pattern1006, replacement1006)
    rubi.add(rule1006.pattern, label = 1006)

    pattern1007 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons117)
    def replacement1007(x, e, b, d, f, c, a):
        return sqrt(S(1) + d*x**S(2)/c)*Int(S(1)/(sqrt(S(1) + d*x**S(2)/c)*(a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/sqrt(c + d*x**S(2))
    rule1007 = ReplacementRule(pattern1007, replacement1007)
    rubi.add(rule1007.pattern, label = 1007)

    pattern1008 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons573)
    def replacement1008(x, e, b, d, f, c, a):
        return c*sqrt(e + f*x**S(2))*EllipticPi(S(1) - b*c/(a*d), ArcTan(x*Rt(d/c, S(2))), -c*f/(d*e) + S(1))/(a*e*sqrt(c*(e + f*x**S(2))/(e*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2)))
    rule1008 = ReplacementRule(pattern1008, replacement1008)
    rubi.add(rule1008.pattern, label = 1008)

    pattern1009 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons581)
    def replacement1009(x, e, b, d, f, c, a):
        return d*Int(S(1)/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/b + (-a*d + b*c)*Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/b
    rule1009 = ReplacementRule(pattern1009, replacement1009)
    rubi.add(rule1009.pattern, label = 1009)

    pattern1010 = Pattern(Integral(sqrt(e_ + x_**S(2)*WC('f', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons573, cons638)
    def replacement1010(x, e, b, d, f, c, a):
        return b*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - d*Int(sqrt(e + f*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c)
    rule1010 = ReplacementRule(pattern1010, replacement1010)
    rubi.add(rule1010.pattern, label = 1010)

    pattern1011 = Pattern(Integral((e_ + x_**S(2)*WC('f', S(1)))**(S(3)/2)/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons573, cons638)
    def replacement1011(x, e, b, d, f, c, a):
        return (-a*f + b*e)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - (-c*f + d*e)*Int(sqrt(e + f*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c)
    rule1011 = ReplacementRule(pattern1011, replacement1011)
    rubi.add(rule1011.pattern, label = 1011)

    pattern1012 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons573, cons638)
    def replacement1012(x, e, b, d, f, c, a):
        return d*Int(sqrt(e + f*x**S(2))*(-a*d + S(2)*b*c + b*d*x**S(2))/sqrt(c + d*x**S(2)), x)/b**S(2) + (-a*d + b*c)**S(2)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/b**S(2)
    rule1012 = ReplacementRule(pattern1012, replacement1012)
    rubi.add(rule1012.pattern, label = 1012)

    pattern1013 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons640, cons396, cons641)
    def replacement1013(q, x, e, b, d, f, r, c, a):
        return b*(-a*f + b*e)*Int((c + d*x**S(2))**(q + S(2))*(e + f*x**S(2))**(r + S(-1))/(a + b*x**S(2)), x)/(-a*d + b*c)**S(2) - Int((c + d*x**S(2))**q*(e + f*x**S(2))**(r + S(-1))*(-a*d**S(2)*e - b*c**S(2)*f + S(2)*b*c*d*e + d**S(2)*x**S(2)*(-a*f + b*e)), x)/(-a*d + b*c)**S(2)
    rule1013 = ReplacementRule(pattern1013, replacement1013)
    rubi.add(rule1013.pattern, label = 1013)

    pattern1014 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons52, cons395, cons576)
    def replacement1014(q, x, e, b, d, f, r, c, a):
        return d*Int((c + d*x**S(2))**(q + S(-1))*(e + f*x**S(2))**r, x)/b + (-a*d + b*c)*Int((c + d*x**S(2))**(q + S(-1))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/b
    rule1014 = ReplacementRule(pattern1014, replacement1014)
    rubi.add(rule1014.pattern, label = 1014)

    pattern1015 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons52, cons395, cons396)
    def replacement1015(q, x, e, b, d, f, r, c, a):
        return b**S(2)*Int((c + d*x**S(2))**(q + S(2))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/(-a*d + b*c)**S(2) - d*Int((c + d*x**S(2))**q*(e + f*x**S(2))**r*(-a*d + S(2)*b*c + b*d*x**S(2)), x)/(-a*d + b*c)**S(2)
    rule1015 = ReplacementRule(pattern1015, replacement1015)
    rubi.add(rule1015.pattern, label = 1015)

    pattern1016 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons52, cons395, cons642)
    def replacement1016(q, x, e, b, d, f, r, c, a):
        return b*Int((c + d*x**S(2))**(q + S(1))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/(-a*d + b*c) - d*Int((c + d*x**S(2))**q*(e + f*x**S(2))**r, x)/(-a*d + b*c)
    rule1016 = ReplacementRule(pattern1016, replacement1016)
    rubi.add(rule1016.pattern, label = 1016)

    pattern1017 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1)))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons153)
    def replacement1017(x, e, b, d, f, c, a):
        return x*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*a*(a + b*x**S(2))) + d*f*Int((a - b*x**S(2))/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*b**S(2)) + (-a**S(2)*d*f + b**S(2)*c*e)*Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*b**S(2))
    rule1017 = ReplacementRule(pattern1017, replacement1017)
    rubi.add(rule1017.pattern, label = 1017)

    pattern1018 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**S(2)*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons153)
    def replacement1018(x, e, b, d, f, c, a):
        return b**S(2)*x*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*a*(a + b*x**S(2))*(-a*d + b*c)*(-a*f + b*e)) - d*f*Int((a + b*x**S(2))/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*(-a*d + b*c)*(-a*f + b*e)) + (S(3)*a**S(2)*d*f - S(2)*a*b*(c*f + d*e) + b**S(2)*c*e)*Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*(-a*d + b*c)*(-a*f + b*e))
    rule1018 = ReplacementRule(pattern1018, replacement1018)
    rubi.add(rule1018.pattern, label = 1018)

    pattern1019 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons52, cons63, cons395, cons403)
    def replacement1019(n, p, q, x, e, b, d, f, r, c, a):
        return d*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(e + f*x**n)**r, x)/b + (-a*d + b*c)*Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*(e + f*x**n)**r, x)/b
    rule1019 = ReplacementRule(pattern1019, replacement1019)
    rubi.add(rule1019.pattern, label = 1019)

    pattern1020 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons50, cons63, cons395, cons642)
    def replacement1020(n, p, q, x, e, b, d, f, r, c, a):
        return b*Int((a + b*x**n)**p*(c + d*x**n)**(q + S(1))*(e + f*x**n)**r, x)/(-a*d + b*c) - d*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(e + f*x**n)**r, x)/(-a*d + b*c)
    rule1020 = ReplacementRule(pattern1020, replacement1020)
    rubi.add(rule1020.pattern, label = 1020)

    pattern1021 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons153)
    def replacement1021(x, e, b, d, f, c, a):
        return sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))*Subst(Int(S(1)/(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)*sqrt(S(1) - x**S(2)*(-a*f + b*e)/e)), x), x, x/sqrt(a + b*x**S(2)))/(c*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2)))
    rule1021 = ReplacementRule(pattern1021, replacement1021)
    rubi.add(rule1021.pattern, label = 1021)

    pattern1022 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons153)
    def replacement1022(x, e, b, d, f, c, a):
        return a*sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))*Subst(Int(S(1)/(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)*sqrt(S(1) - x**S(2)*(-a*f + b*e)/e)*(-b*x**S(2) + S(1))), x), x, x/sqrt(a + b*x**S(2)))/(c*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2)))
    rule1022 = ReplacementRule(pattern1022, replacement1022)
    rubi.add(rule1022.pattern, label = 1022)

    pattern1023 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))**(S(3)/2)*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons153)
    def replacement1023(x, e, b, d, f, c, a):
        return sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))*Subst(Int(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)/sqrt(S(1) - x**S(2)*(-a*f + b*e)/e), x), x, x/sqrt(a + b*x**S(2)))/(a*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2)))
    rule1023 = ReplacementRule(pattern1023, replacement1023)
    rubi.add(rule1023.pattern, label = 1023)

    pattern1024 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/sqrt(e_ + x_**S(2)*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons643)
    def replacement1024(x, e, b, d, f, c, a):
        return b*c*(-c*f + d*e)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*d*f) - c*(-c*f + d*e)*Int(sqrt(a + b*x**S(2))/((c + d*x**S(2))**(S(3)/2)*sqrt(e + f*x**S(2))), x)/(S(2)*f) + d*x*sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*f*sqrt(c + d*x**S(2))) - (-a*d*f - b*c*f + b*d*e)*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*d*f)
    rule1024 = ReplacementRule(pattern1024, replacement1024)
    rubi.add(rule1024.pattern, label = 1024)

    pattern1025 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/sqrt(e_ + x_**S(2)*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons644)
    def replacement1025(x, e, b, d, f, c, a):
        return e*(-a*f + b*e)*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*(e + f*x**S(2))**(S(3)/2)), x)/(S(2)*f) + x*sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))/(S(2)*sqrt(e + f*x**S(2))) + (-a*f + b*e)*(-S(2)*c*f + d*e)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*f**S(2)) - (-a*d*f - b*c*f + b*d*e)*Int(sqrt(e + f*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(S(2)*f**S(2))
    rule1025 = ReplacementRule(pattern1025, replacement1025)
    rubi.add(rule1025.pattern, label = 1025)

    pattern1026 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/(e_ + x_**S(2)*WC('f', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons153)
    def replacement1026(x, e, b, d, f, c, a):
        return b*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/f - (-a*f + b*e)*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*(e + f*x**S(2))**(S(3)/2)), x)/f
    rule1026 = ReplacementRule(pattern1026, replacement1026)
    rubi.add(rule1026.pattern, label = 1026)

    def With1027(n, p, q, x, e, b, d, f, r, c, a):
        u = ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)
        if SumQ(u):
            return Int(u, x)
        return False
    pattern1027 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons52, cons148, CustomConstraint(With1027))
    rule1027 = ReplacementRule(pattern1027, With1027)
    rubi.add(rule1027.pattern, label = 1027)

    pattern1028 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons52, cons196)
    def replacement1028(n, p, q, x, e, b, d, f, r, c, a):
        return -Subst(Int((a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r/x**S(2), x), x, S(1)/x)
    rule1028 = ReplacementRule(pattern1028, replacement1028)
    rubi.add(rule1028.pattern, label = 1028)

    pattern1029 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons50, cons52, cons645)
    def replacement1029(n, p, q, x, e, b, d, f, r, c, a):
        return Int((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)
    rule1029 = ReplacementRule(pattern1029, replacement1029)
    rubi.add(rule1029.pattern, label = 1029)

    pattern1030 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons4, cons50, cons52, cons646, cons647, cons68, cons69)
    def replacement1030(n, w, p, q, v, e, x, b, d, f, r, c, u, a):
        return Subst(Int((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x, u)/Coefficient(u, x, S(1))
    rule1030 = ReplacementRule(pattern1030, replacement1030)
    rubi.add(rule1030.pattern, label = 1030)

    pattern1031 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons52, cons585, cons586)
    def replacement1031(n, q, p, x, e, b, d, f, r, c, mn, a):
        return Int(x**(-n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x)
    rule1031 = ReplacementRule(pattern1031, replacement1031)
    rubi.add(rule1031.pattern, label = 1031)

    pattern1032 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons50, cons585, cons38, cons648)
    def replacement1032(n, q, p, x, e, b, d, f, r, c, mn, a):
        return Int(x**(n*(p + r))*(c + d*x**(-n))**q*(a*x**(-n) + b)**p*(e*x**(-n) + f)**r, x)
    rule1032 = ReplacementRule(pattern1032, replacement1032)
    rubi.add(rule1032.pattern, label = 1032)

    pattern1033 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons50, cons52, cons585, cons386)
    def replacement1033(n, p, q, x, e, b, d, f, r, c, mn, a):
        return x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(-n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x)
    rule1033 = ReplacementRule(pattern1033, replacement1033)
    rubi.add(rule1033.pattern, label = 1033)

    pattern1034 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons652, cons653, cons654, cons655, cons4, cons5, cons50, cons52, cons649, cons650, cons651)
    def replacement1034(n, n2, p, q, e2, c, f1, e1, x, b, d, r, f2, a):
        return Int((a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x)
    rule1034 = ReplacementRule(pattern1034, replacement1034)
    rubi.add(rule1034.pattern, label = 1034)

    pattern1035 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons652, cons653, cons654, cons655, cons4, cons5, cons50, cons52, cons649, cons650)
    def replacement1035(n, n2, p, q, e2, c, f1, e1, x, b, d, r, f2, a):
        return (e1 + f1*x**(n/S(2)))**FracPart(r)*(e2 + f2*x**(n/S(2)))**FracPart(r)*(e1*e2 + f1*f2*x**n)**(-FracPart(r))*Int((a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x)
    rule1035 = ReplacementRule(pattern1035, replacement1035)
    rubi.add(rule1035.pattern, label = 1035)

    pattern1036 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons50, cons52, cons656, cons500)
    def replacement1036(m, n, q, p, g, x, e, b, d, f, r, c):
        return b**(S(1) - (m + S(1))/n)*g**m*Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n)*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n
    rule1036 = ReplacementRule(pattern1036, replacement1036)
    rubi.add(rule1036.pattern, label = 1036)

    pattern1037 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons50, cons52, cons656, cons501)
    def replacement1037(m, n, q, p, g, x, e, b, d, f, r, c):
        return b**IntPart(p)*g**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p)*(c + d*x**n)**q*(e + f*x**n)**r, x)
    rule1037 = ReplacementRule(pattern1037, replacement1037)
    rubi.add(rule1037.pattern, label = 1037)

    pattern1038 = Pattern(Integral((g_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons50, cons52, cons18)
    def replacement1038(m, n, q, p, g, x, e, b, d, f, r, c):
        return g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)
    rule1038 = ReplacementRule(pattern1038, replacement1038)
    rubi.add(rule1038.pattern, label = 1038)

    pattern1039 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons657)
    def replacement1039(m, n, p, q, g, x, e, b, d, f, r, c, a):
        return Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x)
    rule1039 = ReplacementRule(pattern1039, replacement1039)
    rubi.add(rule1039.pattern, label = 1039)

    pattern1040 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons52, cons53)
    def replacement1040(m, n, p, q, x, e, b, d, f, r, c, a):
        return Subst(Int((a + b*x)**p*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n
    rule1040 = ReplacementRule(pattern1040, replacement1040)
    rubi.add(rule1040.pattern, label = 1040)

    pattern1041 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons658, cons502)
    def replacement1041(m, n, p, q, x, e, b, d, f, r, c, a):
        return Int(x**(m + n*(p + q + r))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q*(e*x**(-n) + f)**r, x)
    rule1041 = ReplacementRule(pattern1041, replacement1041)
    rubi.add(rule1041.pattern, label = 1041)

    pattern1042 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons52, cons500)
    def replacement1042(m, n, p, q, x, e, b, d, f, r, c, a):
        return Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n
    rule1042 = ReplacementRule(pattern1042, replacement1042)
    rubi.add(rule1042.pattern, label = 1042)

    pattern1043 = Pattern(Integral((g_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons50, cons52, cons500)
    def replacement1043(m, n, p, q, g, x, e, b, d, f, r, c, a):
        return g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)
    rule1043 = ReplacementRule(pattern1043, replacement1043)
    rubi.add(rule1043.pattern, label = 1043)

    def With1044(m, n, p, q, x, e, b, d, f, r, c, a):
        k = GCD(m + S(1), n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p*(c + d*x**(n/k))**q*(e + f*x**(n/k))**r, x), x, x**k)/k
        return False
    pattern1044 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons52, cons148, cons17, CustomConstraint(With1044))
    rule1044 = ReplacementRule(pattern1044, With1044)
    rubi.add(rule1044.pattern, label = 1044)

    def With1045(m, n, p, q, g, x, e, b, d, f, r, c, a):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*g**(-n)*x**(k*n))**p*(c + d*g**(-n)*x**(k*n))**q*(e + f*g**(-n)*x**(k*n))**r, x), x, (g*x)**(S(1)/k))/g
    pattern1045 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons50, cons52, cons148, cons367 )
    rule1045 = ReplacementRule(pattern1045, With1045)
    rubi.add(rule1045.pattern, label = 1045)

    pattern1046 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons148, cons402, cons137, cons403, cons659)
    def replacement1046(m, n, q, p, g, x, e, b, d, f, c, a):
        return Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + S(1)) + (m + S(1))*(-a*f + b*e)) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(m + n*q + S(1))), x), x)/(a*b*n*(p + S(1))) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*g*n*(p + S(1)))
    rule1046 = ReplacementRule(pattern1046, replacement1046)
    rubi.add(rule1046.pattern, label = 1046)

    pattern1047 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons50, cons148, cons244, cons137, cons607)
    def replacement1047(m, n, p, q, g, x, e, b, d, f, c, a):
        return -g**n*Int((g*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(-a*f + b*e)*(m - n + S(1)) + x**n*(-b*n*(p + S(1))*(c*f - d*e) + d*(-a*f + b*e)*(m + n*q + S(1))), x), x)/(b*n*(p + S(1))*(-a*d + b*c)) + g**(n + S(-1))*(g*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(b*n*(p + S(1))*(-a*d + b*c))
    rule1047 = ReplacementRule(pattern1047, replacement1047)
    rubi.add(rule1047.pattern, label = 1047)

    pattern1048 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons50, cons148, cons13, cons137)
    def replacement1048(m, n, p, q, g, x, e, b, d, f, c, a):
        return Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(m + n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*g*n*(p + S(1))*(-a*d + b*c))
    rule1048 = ReplacementRule(pattern1048, replacement1048)
    rubi.add(rule1048.pattern, label = 1048)

    pattern1049 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons148, cons611, cons403, cons94, cons660)
    def replacement1049(m, n, q, p, g, x, e, b, d, f, c, a):
        return e*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*g*(m + S(1))) - g**(-n)*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + e*n*(a*d*q + b*c*(p + S(1))), x), x)/(a*(m + S(1)))
    rule1049 = ReplacementRule(pattern1049, replacement1049)
    rubi.add(rule1049.pattern, label = 1049)

    pattern1050 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons148, cons395, cons403, cons660)
    def replacement1050(m, n, p, q, g, x, e, b, d, f, c, a):
        return f*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*g*(m + n*(p + q + S(1)) + S(1))) + Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + x**n*(b*d*e*n*(p + q + S(1)) + d*(m + S(1))*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q + S(1)) + S(1)))
    rule1050 = ReplacementRule(pattern1050, replacement1050)
    rubi.add(rule1050.pattern, label = 1050)

    pattern1051 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons50, cons148, cons31, cons530)
    def replacement1051(m, n, p, q, g, x, e, b, d, f, c, a):
        return f*g**(n + S(-1))*(g*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*d*(m + n*(p + q + S(1)) + S(1))) - g**n*Int((g*x)**(m - n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*f*(m - n + S(1)) + x**n*(a*d*f*(m + n*q + S(1)) + b*(c*f*(m + n*p + S(1)) - d*e*(m + n*(p + q + S(1)) + S(1)))), x), x)/(b*d*(m + n*(p + q + S(1)) + S(1)))
    rule1051 = ReplacementRule(pattern1051, replacement1051)
    rubi.add(rule1051.pattern, label = 1051)

    pattern1052 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons50, cons148, cons31, cons94)
    def replacement1052(m, n, q, p, g, x, e, b, d, f, c, a):
        return e*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c*g*(m + S(1))) + g**(-n)*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*f*(m + S(1)) - b*d*e*x**n*(m + n*(p + q + S(2)) + S(1)) - e*n*(a*d*q + b*c*p) - e*(a*d + b*c)*(m + n + S(1)), x), x)/(a*c*(m + S(1)))
    rule1052 = ReplacementRule(pattern1052, replacement1052)
    rubi.add(rule1052.pattern, label = 1052)

    pattern1053 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons148)
    def replacement1053(m, n, p, g, x, e, b, d, f, c, a):
        return Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(e + f*x**n)/(c + d*x**n), x), x)
    rule1053 = ReplacementRule(pattern1053, replacement1053)
    rubi.add(rule1053.pattern, label = 1053)

    pattern1054 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons50, cons148)
    def replacement1054(m, n, p, q, g, x, e, b, d, f, c, a):
        return e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x) + e**(-n)*f*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q, x)
    rule1054 = ReplacementRule(pattern1054, replacement1054)
    rubi.add(rule1054.pattern, label = 1054)

    pattern1055 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons50, cons148, cons661)
    def replacement1055(m, n, p, q, g, x, e, b, d, f, r, c, a):
        return e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**(r + S(-1)), x) + e**(-n)*f*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**(r + S(-1)), x)
    rule1055 = ReplacementRule(pattern1055, replacement1055)
    rubi.add(rule1055.pattern, label = 1055)

    pattern1056 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons52, cons196, cons17)
    def replacement1056(m, n, p, q, x, e, b, d, f, r, c, a):
        return -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r, x), x, S(1)/x)
    rule1056 = ReplacementRule(pattern1056, replacement1056)
    rubi.add(rule1056.pattern, label = 1056)

    def With1057(m, n, q, p, g, x, e, b, d, f, r, c, a):
        k = Denominator(m)
        return -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*g**(-n)*x**(-k*n))**p*(c + d*g**(-n)*x**(-k*n))**q*(e + f*g**(-n)*x**(-k*n))**r, x), x, (g*x)**(-S(1)/k))/g
    pattern1057 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons50, cons52, cons196, cons367 )
    rule1057 = ReplacementRule(pattern1057, With1057)
    rubi.add(rule1057.pattern, label = 1057)

    pattern1058 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons50, cons52, cons196, cons356)
    def replacement1058(m, n, q, p, g, x, e, b, d, f, r, c, a):
        return -(g*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r, x), x, S(1)/x)
    rule1058 = ReplacementRule(pattern1058, replacement1058)
    rubi.add(rule1058.pattern, label = 1058)

    def With1059(m, n, p, q, x, e, b, d, f, r, c, a):
        k = Denominator(n)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n))**p*(c + d*x**(k*n))**q*(e + f*x**(k*n))**r, x), x, x**(S(1)/k))
    pattern1059 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons5, cons50, cons52, cons489 )
    rule1059 = ReplacementRule(pattern1059, With1059)
    rubi.add(rule1059.pattern, label = 1059)

    pattern1060 = Pattern(Integral((g_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons5, cons50, cons52, cons489)
    def replacement1060(m, n, p, q, g, x, e, b, d, f, r, c, a):
        return g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)
    rule1060 = ReplacementRule(pattern1060, replacement1060)
    rubi.add(rule1060.pattern, label = 1060)

    pattern1061 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons52, cons541)
    def replacement1061(m, n, p, q, x, e, b, d, f, r, c, a):
        return Subst(Int((a + b*x**(n/(m + S(1))))**p*(c + d*x**(n/(m + S(1))))**q*(e + f*x**(n/(m + S(1))))**r, x), x, x**(m + S(1)))/(m + S(1))
    rule1061 = ReplacementRule(pattern1061, replacement1061)
    rubi.add(rule1061.pattern, label = 1061)

    pattern1062 = Pattern(Integral((g_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons50, cons52, cons541)
    def replacement1062(m, n, p, q, g, x, e, b, d, f, r, c, a):
        return g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)
    rule1062 = ReplacementRule(pattern1062, replacement1062)
    rubi.add(rule1062.pattern, label = 1062)

    pattern1063 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons402, cons137, cons403, cons659)
    def replacement1063(m, n, q, p, g, x, e, b, d, f, c, a):
        return Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + S(1)) + (m + S(1))*(-a*f + b*e)) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(m + n*q + S(1))), x), x)/(a*b*n*(p + S(1))) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*g*n*(p + S(1)))
    rule1063 = ReplacementRule(pattern1063, replacement1063)
    rubi.add(rule1063.pattern, label = 1063)

    pattern1064 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons50, cons13, cons137)
    def replacement1064(m, n, p, q, g, x, e, b, d, f, c, a):
        return Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(m + n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*g*n*(p + S(1))*(-a*d + b*c))
    rule1064 = ReplacementRule(pattern1064, replacement1064)
    rubi.add(rule1064.pattern, label = 1064)

    pattern1065 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons395, cons403, cons660)
    def replacement1065(m, n, p, q, g, x, e, b, d, f, c, a):
        return f*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*g*(m + n*(p + q + S(1)) + S(1))) + Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + x**n*(b*d*e*n*(p + q + S(1)) + d*(m + S(1))*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q + S(1)) + S(1)))
    rule1065 = ReplacementRule(pattern1065, replacement1065)
    rubi.add(rule1065.pattern, label = 1065)

    pattern1066 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons380)
    def replacement1066(m, n, p, g, x, e, b, d, f, c, a):
        return Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(e + f*x**n)/(c + d*x**n), x), x)
    rule1066 = ReplacementRule(pattern1066, replacement1066)
    rubi.add(rule1066.pattern, label = 1066)

    pattern1067 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons50, cons662)
    def replacement1067(m, n, p, q, g, x, e, b, d, f, c, a):
        return e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x) + f*x**(-m)*(g*x)**m*Int(x**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q, x)
    rule1067 = ReplacementRule(pattern1067, replacement1067)
    rubi.add(rule1067.pattern, label = 1067)

    pattern1068 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons52, cons585, cons586)
    def replacement1068(m, n, q, p, x, e, b, d, f, r, c, mn, a):
        return Int(x**(m - n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x)
    rule1068 = ReplacementRule(pattern1068, replacement1068)
    rubi.add(rule1068.pattern, label = 1068)

    pattern1069 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons50, cons585, cons38, cons648)
    def replacement1069(m, n, q, p, x, e, b, d, f, r, c, mn, a):
        return Int(x**(m + n*(p + r))*(c + d*x**(-n))**q*(a*x**(-n) + b)**p*(e*x**(-n) + f)**r, x)
    rule1069 = ReplacementRule(pattern1069, replacement1069)
    rubi.add(rule1069.pattern, label = 1069)

    pattern1070 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons52, cons585, cons386)
    def replacement1070(m, n, p, q, x, e, b, d, f, r, c, mn, a):
        return x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(m - n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x)
    rule1070 = ReplacementRule(pattern1070, replacement1070)
    rubi.add(rule1070.pattern, label = 1070)

    pattern1071 = Pattern(Integral((g_*x_)**m_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons50, cons52, cons585)
    def replacement1071(m, n, q, p, g, x, e, b, d, f, r, c, mn, a):
        return g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**(-n))**q*(e + f*x**n)**r, x)
    rule1071 = ReplacementRule(pattern1071, replacement1071)
    rubi.add(rule1071.pattern, label = 1071)

    pattern1072 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons50, cons52, cons663)
    def replacement1072(m, n, p, q, g, x, e, b, d, f, r, c, a):
        return Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)
    rule1072 = ReplacementRule(pattern1072, replacement1072)
    rubi.add(rule1072.pattern, label = 1072)

    pattern1073 = Pattern(Integral(u_**WC('m', S(1))*(e_ + v_**n_*WC('f', S(1)))**WC('r', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons52, cons554)
    def replacement1073(m, n, p, q, v, e, x, b, d, f, r, c, u, a):
        return u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x, v)/Coefficient(v, x, S(1))
    rule1073 = ReplacementRule(pattern1073, replacement1073)
    rubi.add(rule1073.pattern, label = 1073)

    pattern1074 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons652, cons653, cons654, cons655, cons208, cons21, cons4, cons5, cons50, cons52, cons649, cons650, cons651)
    def replacement1074(m, n, n2, p, q, e2, c, f1, g, x, e1, b, d, r, f2, a):
        return Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x)
    rule1074 = ReplacementRule(pattern1074, replacement1074)
    rubi.add(rule1074.pattern, label = 1074)

    pattern1075 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons27, cons652, cons653, cons654, cons655, cons208, cons21, cons4, cons5, cons50, cons52, cons649, cons650)
    def replacement1075(m, n, n2, p, q, e2, c, f1, g, x, e1, b, d, r, f2, a):
        return (e1 + f1*x**(n/S(2)))**FracPart(r)*(e2 + f2*x**(n/S(2)))**FracPart(r)*(e1*e2 + f1*f2*x**n)**(-FracPart(r))*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x)
    rule1075 = ReplacementRule(pattern1075, replacement1075)
    rubi.add(rule1075.pattern, label = 1075)

    return (rubi, (rule689, rule690, rule691, rule692, rule693, rule694, rule695, rule696, rule697, rule698, rule699, rule700, rule701, rule702, rule703, rule704, rule705, rule706, rule707, rule708, rule709, rule710, rule711, rule712, rule713, rule714, rule715, rule716, rule717, rule718, rule719, rule720, rule721, rule722, rule723, rule724, rule725, rule726, rule727, rule728, rule729, rule730, rule731, rule732, rule733, rule734, rule735, rule736, rule737, rule738, rule739, rule740, rule741, rule742, rule743, rule744, rule745, rule746, rule747, rule748, rule749, rule750, rule751, rule752, rule753, rule754, rule755, rule756, rule757, rule758, rule759, rule760, rule761, rule762, rule763, rule764, rule765, rule766, rule767, rule768, rule769, rule770, rule771, rule772, rule773, rule774, rule775, rule776, rule777, rule778, rule779, rule780, rule781, rule782, rule783, rule784, rule785, rule786, rule787, rule788, rule789, rule790, rule791, rule792, rule793, rule794, rule795, rule796, rule797, rule798, rule799, rule800, rule801, rule802, rule803, rule804, rule805, rule806, rule807, rule808, rule809, rule810, rule811, rule812, rule813, rule814, rule815, rule816, rule817, rule818, rule819, rule820, rule821, rule822, rule823, rule824, rule825, rule826, rule827, rule828, rule829, rule830, rule831, rule832, rule833, rule834, rule835, rule836, rule837, rule838, rule839, rule840, rule841, rule842, rule843, rule844, rule845, rule846, rule847, rule848, rule849, rule850, rule851, rule852, rule853, rule854, rule855, rule856, rule857, rule858, rule859, rule860, rule861, rule862, rule863, rule864, rule865, rule866, rule867, rule868, rule869, rule870, rule871, rule872, rule873, rule874, rule875, rule876, rule877, rule878, rule879, rule880, rule881, rule882, rule883, rule884, rule885, rule886, rule887, rule888, rule889, rule890, rule891, rule892, rule893, rule894, rule895, rule896, rule897, rule898, rule899, rule900, rule901, rule902, rule903, rule904, rule905, rule906, rule907, rule908, rule909, rule910, rule911, rule912, rule913, rule914, rule915, rule916, rule917, rule918, rule919, rule920, rule921, rule922, rule923, rule924, rule925, rule926, rule927, rule928, rule929, rule930, rule931, rule932, rule933, rule934, rule935, rule936, rule937, rule938, rule939, rule940, rule941, rule942, rule943, rule944, rule945, rule946, rule947, rule948, rule949, rule950, rule951, rule952, rule953, rule954, rule955, rule956, rule957, rule958, rule959, rule960, rule961, rule962, rule963, rule964, rule965, rule966, rule967, rule968, rule969, rule970, rule971, rule972, rule973, rule974, rule975, rule976, rule977, rule978, rule979, rule980, rule981, rule982, rule983, rule984, rule985, rule986, rule987, rule988, rule989, rule990, rule991, rule992, rule993, rule994, rule995, rule996, rule997, rule998, rule999, rule1000, rule1001, rule1002, rule1003, rule1004, rule1005, rule1006, rule1007, rule1008, rule1009, rule1010, rule1011, rule1012, rule1013, rule1014, rule1015, rule1016, rule1017, rule1018, rule1019, rule1020, rule1021, rule1022, rule1023, rule1024, rule1025, rule1026, rule1027, rule1028, rule1029, rule1030, rule1031, rule1032, rule1033, rule1034, rule1035, rule1036, rule1037, rule1038, rule1039, rule1040, rule1041, rule1042, rule1043, rule1044, rule1045, rule1046, rule1047, rule1048, rule1049, rule1050, rule1051, rule1052, rule1053, rule1054, rule1055, rule1056, rule1057, rule1058, rule1059, rule1060, rule1061, rule1062, rule1063, rule1064, rule1065, rule1066, rule1067, rule1068, rule1069, rule1070, rule1071, rule1072, rule1073, rule1074, rule1075, ))
