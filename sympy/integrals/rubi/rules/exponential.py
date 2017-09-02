from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (Int, Set, With, Module, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcSec, ArcCsch, ArcSech, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct, SplitSum, SubstFor, SubstForAux, FresnelS, FresnelC, Erfc, Erfi, Gamma, FunctionOfTrigOfLinearQ, ElementaryFunctionQ, Complex, UnsameQ, _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify, _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum, _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux, TrigSimplifyAux, Cancel, Part, PolyLog, D, Dist)
    from sympy import Integral, S, sqrt
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)

    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]

    _UseGamma = False

def exponential(rubi):

    pattern1 = Pattern(Integral((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda m: IntegerQ(S(2)*m)), CustomConstraint(lambda : Not(SameQ(_UseGamma, True))))
    rule1 = ReplacementRule(pattern1, lambda d, g, n, F, m, b, e, c, f, x : -d*m*Int((F**(g*(e + f*x))*b)**n*(c + d*x)**(m + S(-1)), x)/(f*g*n*log(F)) + (F**(g*(e + f*x))*b)**n*(c + d*x)**m/(f*g*n*log(F)))
    rubi.add(rule1)

    pattern2 = Pattern(Integral((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: IntegerQ(S(2)*m)), CustomConstraint(lambda : Not(SameQ(_UseGamma, True))))
    rule2 = ReplacementRule(pattern2, lambda d, g, n, F, m, b, e, c, f, x : -f*g*n*Int((F**(g*(e + f*x))*b)**n*(c + d*x)**(m + S(1)), x)*log(F)/(d*(m + S(1))) + (F**(g*(e + f*x))*b)**n*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule2)

    pattern3 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda : Not(SameQ(_UseGamma, True))))
    rule3 = ReplacementRule(pattern3, lambda d, g, F, e, c, f, x : F**(g*(-c*f/d + e))*ExpIntegralEi(f*g*(c + d*x)*log(F)/d)/d)
    rubi.add(rule3)

    pattern4 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule4 = ReplacementRule(pattern4, lambda d, g, F, m, e, c, f, x : F**(g*(-c*f/d + e))*f**(-m + S(-1))*g**(-m + S(-1))*(-d)**m*Gamma(m + S(1), -f*g*(c + d*x)*log(F)/d)*log(F)**(-m + S(-1)))
    rubi.add(rule4)

    pattern5 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))/sqrt(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda : Not(SameQ(_UseGamma, True))))
    rule5 = ReplacementRule(pattern5, lambda d, g, F, e, c, f, x : S(2)*Subst(Int(F**(g*(-c*f/d + e) + f*g*x**S(2)/d), x), x, sqrt(c + d*x))/d)
    rubi.add(rule5)

    pattern6 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule6 = ReplacementRule(pattern6, lambda d, g, F, m, e, c, f, x : -F**(g*(-c*f/d + e))*(-f*g*log(F)/d)**(-IntPart(m) + S(-1))*(-f*g*(c + d*x)*log(F)/d)**(-FracPart(m))*(c + d*x)**FracPart(m)*Gamma(m + S(1), -f*g*(c + d*x)*log(F)/d)/d)
    rubi.add(rule6)

    pattern7 = Pattern(Integral((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*WC('b', S(1)))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule7 = ReplacementRule(pattern7, lambda d, g, n, F, m, b, e, c, f, x : F**(-g*n*(e + f*x))*(F**(g*(e + f*x))*b)**n*Int(F**(g*n*(e + f*x))*(c + d*x)**m, x))
    rubi.add(rule7)

    pattern8 = Pattern(Integral((a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule8 = ReplacementRule(pattern8, lambda d, g, n, F, a, p, m, b, e, c, f, x : Int(ExpandIntegrand((c + d*x)**m, (a + b*(F**(g*(e + f*x)))**n)**p, x), x))
    rubi.add(rule8)

    pattern9 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule9 = ReplacementRule(pattern9, lambda d, g, n, F, a, m, b, e, c, f, x : d*m*Int((c + d*x)**(m + S(-1))*log(a*(F**(g*(e + f*x)))**(-n)/b + S(1)), x)/(a*f*g*n*log(F)) - (c + d*x)**m*log(a*(F**(g*(e + f*x)))**(-n)/b + S(1))/(a*f*g*n*log(F)))
    rubi.add(rule9)

    pattern10 = Pattern(Integral((a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)))**p_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), )
    def With10(d, g, n, F, a, p, m, b, e, c, f, x):
        u = IntHide((a + b*(F**(g*(e + f*x)))**n)**p, x)
        return -d*m*Int(u*(c + d*x)**(m + S(-1)), x) + Dist((c + d*x)**m, u, x)
    rule10 = ReplacementRule(pattern10, lambda d, g, n, F, a, p, m, b, e, c, f, x : With10(d, g, n, F, a, p, m, b, e, c, f, x))
    rubi.add(rule10)

    pattern11 = Pattern(Integral(u_**WC('m', S(1))*((F_**(v_*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda x, u: PowerOfLinearQ(u, x)), CustomConstraint(lambda u, v, x: Not(LinearMatchQ(v, x) & PowerOfLinearMatchQ(u, x))), CustomConstraint(lambda m: IntegerQ(m)))
    rule11 = ReplacementRule(pattern11, lambda x, g, n, F, a, p, m, b, v, u : Int((a + b*(F**(g*ExpandToSum(v, x)))**n)**p*NormalizePowerOfLinear(u, x)**m, x))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(u_**WC('m', S(1))*((F_**(v_*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda x, u: PowerOfLinearQ(u, x)), CustomConstraint(lambda u, v, x: Not(LinearMatchQ(v, x) & PowerOfLinearMatchQ(u, x))), CustomConstraint(lambda m: Not(IntegerQ(m))), )
    def With12(x, g, n, F, a, p, m, b, v, u):
        uu = NormalizePowerOfLinear(u, x)
        z = Symbol('z')
        z = If(PowerQ(uu), Part(uu, 1)**(m*Part(uu, 2)), uu**m)
        return uu**m*Int(z*(a + b*(F**(g*ExpandToSum(v, x)))**n)**p, x)/z
    rule12 = ReplacementRule(pattern12, lambda x, g, n, F, a, p, m, b, v, u : With12(x, g, n, F, a, p, m, b, v, u))
    rubi.add(rule12)

    pattern13 = Pattern(Integral((a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule13 = ReplacementRule(pattern13, lambda d, g, n, F, a, p, m, b, e, c, f, x : Int((a + b*(F**(g*(e + f*x)))**n)**p*(c + d*x)**m, x))
    rubi.add(rule13)

    pattern14 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))/(a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule14 = ReplacementRule(pattern14, lambda d, g, n, F, a, m, b, e, c, f, x : -d*m*Int((c + d*x)**(m + S(-1))*log(S(1) + b*(F**(g*(e + f*x)))**n/a), x)/(b*f*g*n*log(F)) + (c + d*x)**m*log(S(1) + b*(F**(g*(e + f*x)))**n/a)/(b*f*g*n*log(F)))
    rubi.add(rule14)

    pattern15 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule15 = ReplacementRule(pattern15, lambda d, g, n, F, a, p, m, b, e, c, f, x : -d*m*Int((a + b*(F**(g*(e + f*x)))**n)**(p + S(1))*(c + d*x)**(m + S(-1)), x)/(b*f*g*n*(p + S(1))*log(F)) + (a + b*(F**(g*(e + f*x)))**n)**(p + S(1))*(c + d*x)**m/(b*f*g*n*(p + S(1))*log(F)))
    rubi.add(rule15)

    pattern16 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule16 = ReplacementRule(pattern16, lambda d, g, n, F, a, p, m, b, e, c, f, x : Int((a + b*(F**(g*(e + f*x)))**n)**p*(c + d*x)**m*(F**(g*(e + f*x)))**n, x))
    rubi.add(rule16)

    pattern17 = Pattern(Integral((G_**((x_*WC('i', S(1)) + WC('h', S(0)))*WC('j', S(1)))*WC('k', S(1)))**WC('q', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda k, x: FreeQ(k, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda i, g, n, F, q, G, j, f: ZeroQ(f*g*n*log(F) - i*j*q*log(G))), CustomConstraint(lambda i, q, g, F, n, h, G, j, e, k, f, x: NonzeroQ((G**(j*(h + i*x))*k)**q - (F**(g*(e + f*x)))**n)))
    rule17 = ReplacementRule(pattern17, lambda i, d, g, n, q, a, h, p, j, m, F, G, b, e, k, c, f, x : (G**(j*(h + i*x))*k)**q*(F**(g*(e + f*x)))**(-n)*Int((a + b*(F**(g*(e + f*x)))**n)**p*(c + d*x)**m*(F**(g*(e + f*x)))**n, x))
    rubi.add(rule17)

    pattern18 = Pattern(Integral((F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule18 = ReplacementRule(pattern18, lambda n, F, a, b, c, x : (F**(c*(a + b*x)))**n/(b*c*n*log(F)))
    rubi.add(rule18)

    pattern19 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda x, u: PolynomialQ(u, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda : SameQ(_UseGamma, True)))
    rule19 = ReplacementRule(pattern19, lambda u, F, c, v, x : Int(ExpandIntegrand(F**(c*ExpandToSum(v, x))*u, x), x))
    rubi.add(rule19)

    pattern20 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda x, u: PolynomialQ(u, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda : Not(SameQ(_UseGamma, True))))
    rule20 = ReplacementRule(pattern20, lambda u, F, c, v, x : Int(ExpandIntegrand(F**(c*ExpandToSum(v, x)), u, x), x))
    rubi.add(rule20)

    pattern21 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda w, v, x, u: LinearQ(List(u, v, w), x)), CustomConstraint(lambda u, x, F, m, c, v, w: ZeroQ(-c*(-Coefficient(u, x, S(0))*Coefficient(w, x, S(1)) + Coefficient(u, x, S(1))*Coefficient(w, x, S(0)))*Coefficient(v, x, S(1))*log(F) + (m + S(1))*Coefficient(u, x, S(1))*Coefficient(w, x, S(1)))))
    rule21 = ReplacementRule(pattern21, lambda x, w, F, m, c, v, u : F**(c*v)*u**(m + S(1))*Coefficient(w, x, S(1))/(c*Coefficient(u, x, S(1))*Coefficient(v, x, S(1))*log(F)))
    rubi.add(rule21)

    pattern22 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda x, w: PolynomialQ(w, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda x, u: PowerOfLinearQ(u, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda : SameQ(_UseGamma, True)))
    rule22 = ReplacementRule(pattern22, lambda x, w, F, m, c, v, u : Int(ExpandIntegrand(F**(c*ExpandToSum(v, x))*w*NormalizePowerOfLinear(u, x)**m, x), x))
    rubi.add(rule22)

    pattern23 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda x, w: PolynomialQ(w, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda x, u: PowerOfLinearQ(u, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda : Not(SameQ(_UseGamma, True))))
    rule23 = ReplacementRule(pattern23, lambda x, w, F, m, c, v, u : Int(ExpandIntegrand(F**(c*ExpandToSum(v, x)), w*NormalizePowerOfLinear(u, x)**m, x), x))
    rubi.add(rule23)

    pattern24 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, w: PolynomialQ(w, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda x, u: PowerOfLinearQ(u, x)), CustomConstraint(lambda m: Not(IntegerQ(m))), )
    def With24(x, w, F, m, c, v, u):
        uu = NormalizePowerOfLinear(u, x)
        z = Symbol('z')
        z = If(PowerQ(uu), Part(uu, 1)**(m*Part(uu, 2)), uu**m)
        return uu**m*Int(ExpandIntegrand(F**(c*ExpandToSum(v, x))*w*z, x), x)/z
    rule24 = ReplacementRule(pattern24, lambda x, w, F, m, c, v, u : With24(x, w, F, m, c, v, u))
    rubi.add(rule24)

    pattern25 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(e_ + (x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1))*log(x_*WC('d', S(1))))*log(x_*WC('d', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, h, e, f: ZeroQ(e - f*h*(n + S(1)))), CustomConstraint(lambda g, n, F, h, b, e, c: ZeroQ(-b*c*e*log(F) + g*h*(n + S(1)))), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule25 = ReplacementRule(pattern25, lambda d, g, n, F, a, h, b, e, c, f, x : F**(c*(a + b*x))*e*x*log(d*x)**(n + S(1))/(n + S(1)))
    rubi.add(rule25)

    pattern26 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*(e_ + (x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1))*log(x_*WC('d', S(1))))*log(x_*WC('d', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, h, m, e, f: ZeroQ(e*(m + S(1)) - f*h*(n + S(1)))), CustomConstraint(lambda g, n, F, h, b, e, c: ZeroQ(-b*c*e*log(F) + g*h*(n + S(1)))), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule26 = ReplacementRule(pattern26, lambda d, g, n, F, a, h, m, b, e, c, f, x : F**(c*(a + b*x))*e*x**(m + S(1))*log(d*x)**(n + S(1))/(n + S(1)))
    rubi.add(rule26)

    pattern27 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule27 = ReplacementRule(pattern27, lambda d, F, a, b, c, x : F**(a + b*(c + d*x))/(b*d*log(F)))
    rubi.add(rule27)

    pattern28 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b: PosQ(b)))
    rule28 = ReplacementRule(pattern28, lambda d, F, a, b, c, x : F**a*sqrt(Pi)*Erfi((c + d*x)*Rt(b*log(F), S(2)))/(S(2)*d*Rt(b*log(F), S(2))))
    rubi.add(rule28)

    pattern29 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b: NegQ(b)))
    rule29 = ReplacementRule(pattern29, lambda d, F, a, b, c, x : F**a*sqrt(Pi)*Erf((c + d*x)*Rt(-b*log(F), S(2)))/(S(2)*d*Rt(-b*log(F), S(2))))
    rubi.add(rule29)

    pattern30 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(S(2)/n)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule30 = ReplacementRule(pattern30, lambda d, n, F, a, b, c, x : F**(a + b*(c + d*x)**n)*(c + d*x)/d - b*n*Int(F**(a + b*(c + d*x)**n)*(c + d*x)**n, x)*log(F))
    rubi.add(rule30)

    pattern31 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(S(2)/n)), CustomConstraint(lambda n: Not(IntegerQ(n))), )
    def With31(d, n, F, a, b, c, x):
        k = Denominator(n)
        return k*Subst(Int(F**(a + b*x**(k*n))*x**(k + S(-1)), x), x, (c + d*x)**(1/k))/d
    rule31 = ReplacementRule(pattern31, lambda d, n, F, a, b, c, x : With31(d, n, F, a, b, c, x))
    rubi.add(rule31)

    pattern32 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(S(2)/n))))
    rule32 = ReplacementRule(pattern32, lambda d, n, F, a, b, c, x : -F**a*(-b*(c + d*x)**n*log(F))**(-S(1)/n)*(c + d*x)*Gamma(1/n, -b*(c + d*x)**n*log(F))/(d*n))
    rubi.add(rule32)

    pattern33 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))), CustomConstraint(lambda c, d, e, f: ZeroQ(-c*f + d*e)))
    rule33 = ReplacementRule(pattern33, lambda d, n, F, a, m, b, e, c, f, x : F**(a + b*(c + d*x)**n)*(c + d*x)**(-n)*(e + f*x)**n/(b*f*n*log(F)))
    rubi.add(rule33)

    pattern34 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))/(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e, f: ZeroQ(-c*f + d*e)))
    rule34 = ReplacementRule(pattern34, lambda d, n, F, a, b, e, c, f, x : F**a*ExpIntegralEi(b*(c + d*x)**n*log(F))/(f*n))
    rubi.add(rule34)

    pattern35 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(-S(2)*m + n + S(-2))))
    rule35 = ReplacementRule(pattern35, lambda d, n, F, a, m, b, c, x : Subst(Int(F**(a + b*x**S(2)), x), x, (c + d*x)**(m + S(1)))/(d*(m + S(1))))
    rubi.add(rule35)

    pattern36 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n, m: IntegerQ(S(2)*(m + S(1))/n)), CustomConstraint(lambda n, m: Less(S(0), (m + S(1))/n, S(5))), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n, m: Less(m, n, S(0)) | Less(S(0), n, m + S(1))))
    rule36 = ReplacementRule(pattern36, lambda d, n, F, a, m, b, c, x : F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n + S(1))/(b*d*n*log(F)) - (m - n + S(1))*Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n), x)/(b*n*log(F)))
    rubi.add(rule36)

    pattern37 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: IntegerQ(S(2)*(m + S(1))/n)), CustomConstraint(lambda n, m: Less(S(0), (m + S(1))/n, S(5))), CustomConstraint(lambda m: Not(RationalQ(m))), CustomConstraint(lambda n, m: SumSimplerQ(m, -n)))
    rule37 = ReplacementRule(pattern37, lambda d, n, F, a, m, b, c, x : F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n + S(1))/(b*d*n*log(F)) - (m - n + S(1))*Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n), x)/(b*n*log(F)))
    rubi.add(rule37)

    pattern38 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n, m: IntegerQ(S(2)*(m + S(1))/n)), CustomConstraint(lambda n, m: Less(S(-4), (m + S(1))/n, S(5))), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n, m: (Greater(n, S(0)) & Less(m, S(-1))) | Inequality(S(0), Less, -n, LessEqual, m + S(1))))
    rule38 = ReplacementRule(pattern38, lambda d, n, F, a, m, b, c, x : F**(a + b*(c + d*x)**n)*(c + d*x)**(m + S(1))/(d*(m + S(1))) - b*n*Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m + n), x)*log(F)/(m + S(1)))
    rubi.add(rule38)

    pattern39 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: IntegerQ(S(2)*(m + S(1))/n)), CustomConstraint(lambda n, m: Less(S(-4), (m + S(1))/n, S(5))), CustomConstraint(lambda m: Not(RationalQ(m))), CustomConstraint(lambda n, m: SumSimplerQ(m, n)))
    rule39 = ReplacementRule(pattern39, lambda d, n, F, a, m, b, c, x : F**(a + b*(c + d*x)**n)*(c + d*x)**(m + S(1))/(d*(m + S(1))) - b*n*Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m + n), x)*log(F)/(m + S(1)))
    rubi.add(rule39)

    pattern40 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n, m: IntegerQ(S(2)*(m + S(1))/n)), CustomConstraint(lambda n, m: Less(S(0), (m + S(1))/n, S(5))), CustomConstraint(lambda n: Not(IntegerQ(n))), )
    def With40(d, n, F, a, m, b, c, x):
        k = Denominator(n)
        return k*Subst(Int(F**(a + b*x**(k*n))*x**(k*(m + S(1)) + S(-1)), x), x, (c + d*x)**(1/k))/d
    rule40 = ReplacementRule(pattern40, lambda d, n, F, a, m, b, c, x : With40(d, n, F, a, m, b, c, x))
    rubi.add(rule40)

    pattern41 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e, f: ZeroQ(-c*f + d*e)), CustomConstraint(lambda n, m: IntegerQ(S(2)*(m + S(1))/n)), CustomConstraint(lambda d, f: NonzeroQ(-d + f)), CustomConstraint(lambda m: Not(IntegerQ(m))), CustomConstraint(lambda c, e: NonzeroQ(c*e)))
    rule41 = ReplacementRule(pattern41, lambda d, n, F, a, m, b, e, c, f, x : (c + d*x)**(-m)*(e + f*x)**m*Int(F**(a + b*(c + d*x)**n)*(c + d*x)**m, x))
    rubi.add(rule41)

    pattern42 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e, f: ZeroQ(-c*f + d*e)))
    rule42 = ReplacementRule(pattern42, lambda d, n, F, a, m, b, e, c, f, x : -F**a*(-b*(c + d*x)**n*log(F))**(-(m + S(1))/n)*(e + f*x)**(m + S(1))*Gamma((m + S(1))/n, -b*(c + d*x)**n*log(F))/(f*n))
    rubi.add(rule42)

    pattern43 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e, f: NonzeroQ(-c*f + d*e)), CustomConstraint(lambda m: FractionQ(m)), CustomConstraint(lambda m: Greater(m, S(1))))
    rule43 = ReplacementRule(pattern43, lambda d, F, a, m, b, e, c, f, x : F**(a + b*(c + d*x)**S(2))*f*(e + f*x)**(m + S(-1))/(S(2)*b*d**S(2)*log(F)) + (-c*f + d*e)*Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(-1)), x)/d - f**S(2)*(m + S(-1))*Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(-2)), x)/(S(2)*b*d**S(2)*log(F)))
    rubi.add(rule43)

    pattern44 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e, f: NonzeroQ(-c*f + d*e)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule44 = ReplacementRule(pattern44, lambda d, F, a, m, b, e, c, f, x : F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(1))/(f*(m + S(1))) - S(2)*b*d**S(2)*Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(2)), x)*log(F)/(f**S(2)*(m + S(1))) + S(2)*b*d*(-c*f + d*e)*Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(1)), x)*log(F)/(f**S(2)*(m + S(1))))
    rubi.add(rule44)

    pattern45 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e, f: NonzeroQ(-c*f + d*e)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(2))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule45 = ReplacementRule(pattern45, lambda d, n, F, a, m, b, e, c, f, x : F**(a + b*(c + d*x)**n)*(e + f*x)**(m + S(1))/(f*(m + S(1))) - b*d*n*Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(n + S(-1))*(e + f*x)**(m + S(1)), x)*log(F)/(f*(m + S(1))))
    rubi.add(rule45)

    pattern46 = Pattern(Integral(F_**(WC('a', S(0)) + WC('b', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e, f: NonzeroQ(-c*f + d*e)))
    rule46 = ReplacementRule(pattern46, lambda d, F, a, b, e, c, f, x : d*Int(F**(a + b/(c + d*x))/(c + d*x), x)/f - (-c*f + d*e)*Int(F**(a + b/(c + d*x))/((c + d*x)*(e + f*x)), x)/f)
    rubi.add(rule46)

    pattern47 = Pattern(Integral(F_**(WC('a', S(0)) + WC('b', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e, f: NonzeroQ(-c*f + d*e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule47 = ReplacementRule(pattern47, lambda d, F, a, m, b, e, c, f, x : F**(a + b/(c + d*x))*(e + f*x)**(m + S(1))/(f*(m + S(1))) + b*d*Int(F**(a + b/(c + d*x))*(e + f*x)**(m + S(1))/(c + d*x)**S(2), x)*log(F)/(f*(m + S(1))))
    rubi.add(rule47)

    pattern48 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))/(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e, f: NonzeroQ(-c*f + d*e)))
    rule48 = ReplacementRule(pattern48, lambda d, n, F, a, b, e, c, f, x : Int(F**(a + b*(c + d*x)**n)/(e + f*x), x))
    rubi.add(rule48)

    pattern49 = Pattern(Integral(F_**v_*u_**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda v, x: BinomialQ(v, x)), CustomConstraint(lambda v, x, u: Not(BinomialMatchQ(v, x) & LinearMatchQ(u, x))))
    rule49 = ReplacementRule(pattern49, lambda u, F, m, v, x : Int(F**ExpandToSum(v, x)*ExpandToSum(u, x)**m, x))
    rubi.add(rule49)

    pattern50 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*u_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, u: PolynomialQ(u, x)))
    rule50 = ReplacementRule(pattern50, lambda u, d, n, F, a, b, c, x : Int(ExpandLinearProduct(F**(a + b*(c + d*x)**n), u, c, d, x), x))
    rubi.add(rule50)

    pattern51 = Pattern(Integral(F_**(v_*WC('b', S(1)) + WC('a', S(0)))*WC('u', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda x, u: PolynomialQ(u, x)), CustomConstraint(lambda v, x: PowerOfLinearQ(v, x)), CustomConstraint(lambda v, x: Not(PowerOfLinearMatchQ(v, x))))
    rule51 = ReplacementRule(pattern51, lambda x, F, a, b, v, u : Int(F**(a + b*NormalizePowerOfLinear(v, x))*u, x))
    rubi.add(rule51)

    pattern52 = Pattern(Integral(F_**(WC('a', S(0)) + WC('b', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/((x_*WC('f', S(1)) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e, f: ZeroQ(-c*f + d*e)))
    rule52 = ReplacementRule(pattern52, lambda d, g, F, a, h, b, e, c, f, x : -d*Subst(Int(F**(a + b*d*x/(-c*h + d*g) - b*h/(-c*h + d*g))/x, x), x, (g + h*x)/(c + d*x))/(f*(-c*h + d*g)))
    rubi.add(rule52)

    pattern53 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, b, a, d: ZeroQ(-a*d + b*c)))
    rule53 = ReplacementRule(pattern53, lambda d, g, F, a, h, m, b, e, c, f, x : F**(b*f/d + e)*Int((g + h*x)**m, x))
    rubi.add(rule53)

    pattern54 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, b, a, d: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda c, d, h, g: ZeroQ(-c*h + d*g)))
    rule54 = ReplacementRule(pattern54, lambda d, g, F, a, h, m, b, e, c, f, x : Int(F**(-f*(-a*d + b*c)/(d*(c + d*x)) + (b*f + d*e)/d)*(g + h*x)**m, x))
    rubi.add(rule54)

    pattern55 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda c, b, a, d: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda c, d, h, g: NonzeroQ(-c*h + d*g)))
    rule55 = ReplacementRule(pattern55, lambda d, g, F, a, h, b, e, c, f, x : d*Int(F**(e + f*(a + b*x)/(c + d*x))/(c + d*x), x)/h - (-c*h + d*g)*Int(F**(e + f*(a + b*x)/(c + d*x))/((c + d*x)*(g + h*x)), x)/h)
    rubi.add(rule55)

    pattern56 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda c, b, a, d: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda c, d, h, g: NonzeroQ(-c*h + d*g)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule56 = ReplacementRule(pattern56, lambda d, g, F, a, h, m, b, e, c, f, x : F**(e + f*(a + b*x)/(c + d*x))*(g + h*x)**(m + S(1))/(h*(m + S(1))) - f*(-a*d + b*c)*Int(F**(e + f*(a + b*x)/(c + d*x))*(g + h*x)**(m + S(1))/(c + d*x)**S(2), x)*log(F)/(h*(m + S(1))))
    rubi.add(rule56)

    pattern57 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_*WC('j', S(1)) + WC('i', S(0)))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda c, d, h, g: ZeroQ(-c*h + d*g)))
    rule57 = ReplacementRule(pattern57, lambda d, i, g, F, a, h, j, b, e, c, f, x : -d*Subst(Int(F**(e - f*x*(-a*d + b*c)/(-c*j + d*i) + f*(-a*j + b*i)/(-c*j + d*i))/x, x), x, (i + j*x)/(c + d*x))/(h*(-c*j + d*i)))
    rubi.add(rule57)

    pattern58 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule58 = ReplacementRule(pattern58, lambda F, a, b, c, x : F**(a - b**S(2)/(S(4)*c))*Int(F**((b + S(2)*c*x)**S(2)/(S(4)*c)), x))
    rubi.add(rule58)

    pattern59 = Pattern(Integral(F_**v_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda v, x: Not(QuadraticMatchQ(v, x))))
    rule59 = ReplacementRule(pattern59, lambda F, v, x : Int(F**ExpandToSum(v, x), x))
    rubi.add(rule59)

    pattern60 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, b, d, e: ZeroQ(b*e - S(2)*c*d)))
    rule60 = ReplacementRule(pattern60, lambda d, F, a, b, e, c, x : F**(a + b*x + c*x**S(2))*e/(S(2)*c*log(F)))
    rubi.add(rule60)

    pattern61 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, b, d, e: ZeroQ(b*e - S(2)*c*d)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(1))))
    rule61 = ReplacementRule(pattern61, lambda d, F, a, m, b, e, c, x : F**(a + b*x + c*x**S(2))*e*(d + e*x)**(m + S(-1))/(S(2)*c*log(F)) - e**S(2)*(m + S(-1))*Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(-2)), x)/(S(2)*c*log(F)))
    rubi.add(rule61)

    pattern62 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, b, d, e: ZeroQ(b*e - S(2)*c*d)))
    rule62 = ReplacementRule(pattern62, lambda d, F, a, b, e, c, x : F**(a - b**S(2)/(S(4)*c))*ExpIntegralEi((b + S(2)*c*x)**S(2)*log(F)/(S(4)*c))/(S(2)*e))
    rubi.add(rule62)

    pattern63 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, b, d, e: ZeroQ(b*e - S(2)*c*d)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule63 = ReplacementRule(pattern63, lambda d, F, a, m, b, e, c, x : F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(1))/(e*(m + S(1))) - S(2)*c*Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(2)), x)*log(F)/(e**S(2)*(m + S(1))))
    rubi.add(rule63)

    pattern64 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, b, d, e: NonzeroQ(b*e - S(2)*c*d)))
    rule64 = ReplacementRule(pattern64, lambda d, F, a, b, e, c, x : F**(a + b*x + c*x**S(2))*e/(S(2)*c*log(F)) - (b*e - S(2)*c*d)*Int(F**(a + b*x + c*x**S(2)), x)/(S(2)*c))
    rubi.add(rule64)

    pattern65 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, b, d, e: NonzeroQ(b*e - S(2)*c*d)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(1))))
    rule65 = ReplacementRule(pattern65, lambda d, F, a, m, b, e, c, x : F**(a + b*x + c*x**S(2))*e*(d + e*x)**(m + S(-1))/(S(2)*c*log(F)) - e**S(2)*(m + S(-1))*Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(-2)), x)/(S(2)*c*log(F)) - (b*e - S(2)*c*d)*Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(-1)), x)/(S(2)*c))
    rubi.add(rule65)

    pattern66 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, b, d, e: NonzeroQ(b*e - S(2)*c*d)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule66 = ReplacementRule(pattern66, lambda d, F, a, m, b, e, c, x : F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(1))/(e*(m + S(1))) - S(2)*c*Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(2)), x)*log(F)/(e**S(2)*(m + S(1))) - (b*e - S(2)*c*d)*Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(1)), x)*log(F)/(e**S(2)*(m + S(1))))
    rubi.add(rule66)

    pattern67 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule67 = ReplacementRule(pattern67, lambda d, F, a, m, b, e, c, x : Int(F**(a + b*x + c*x**S(2))*(d + e*x)**m, x))
    rubi.add(rule67)

    pattern68 = Pattern(Integral(F_**v_*u_**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda v, x, u: Not(LinearMatchQ(u, x) & QuadraticMatchQ(v, x))))
    rule68 = ReplacementRule(pattern68, lambda u, F, m, v, x : Int(F**ExpandToSum(v, x)*ExpandToSum(u, x)**m, x))
    rubi.add(rule68)

    pattern69 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*x_**WC('m', S(1))*(F_**v_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, e, c, v, x: ZeroQ(S(2)*e*(c + d*x) - v)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: NegativeIntegerQ(n)), )
    def With69(d, n, F, a, m, b, e, c, v, x):
        u = IntHide(F**(e*(c + d*x))*(F**v*b + a)**n, x)
        return -m*Int(u*x**(m + S(-1)), x) + Dist(x**m, u, x)
    rule69 = ReplacementRule(pattern69, lambda d, n, F, a, m, b, e, c, v, x : With69(d, n, F, a, m, b, e, c, v, x))
    rubi.add(rule69)

    pattern70 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, h, a, c, m, n, e, g, G, d, F, b, x: RationalQ(m) & GreaterEqual(Abs(m), S(1))))
    def With70(d, g, n, F, a, h, G, b, e, c, f, x):
        m = FullSimplify(g*h*log(G)/(d*e*log(F)))
        return G**(-c*g*h/d + f*h)*Denominator(m)*Subst(Int(x**(Numerator(m) + S(-1))*(a + b*x**Denominator(m))**n, x), x, F**(e*(c + d*x)/Denominator(m)))/(d*e*log(F))
    rule70 = ReplacementRule(pattern70, lambda d, g, n, F, a, h, G, b, e, c, f, x : With70(d, g, n, F, a, h, G, b, e, c, f, x))
    rubi.add(rule70)

    pattern71 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, h, a, c, m, n, e, g, G, d, F, b, x: RationalQ(m) & Greater(Abs(m), S(1))))
    def With71(d, g, n, F, a, h, G, b, e, c, f, x):
        m = FullSimplify(d*e*log(F)/(g*h*log(G)))
        return Denominator(m)*Subst(Int(x**(Denominator(m) + S(-1))*(F**(c*e - d*e*f/g)*b*x**Numerator(m) + a)**n, x), x, G**(h*(f + g*x)/Denominator(m)))/(g*h*log(G))
    rule71 = ReplacementRule(pattern71, lambda d, g, n, F, a, h, G, b, e, c, f, x : With71(d, g, n, F, a, h, G, b, e, c, f, x))
    rubi.add(rule71)

    pattern72 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda d, g, F, G, h, e: Not(RationalQ(FullSimplify(g*h*log(G)/(d*e*log(F)))))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule72 = ReplacementRule(pattern72, lambda d, g, n, F, a, h, G, b, e, c, f, x : Int(G**(f*h)*G**(g*h*x)*(F**(c*e)*F**(d*e*x)*b + a)**n, x))
    rubi.add(rule72)

    pattern73 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda d, g, F, G, h, e: Not(RationalQ(FullSimplify(g*h*log(G)/(d*e*log(F)))))), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule73 = ReplacementRule(pattern73, lambda d, g, n, F, a, h, G, b, e, c, f, x : G**(h*(f + g*x))*a**n*Hypergeometric2F1(-n, g*h*log(G)/(d*e*log(F)), S(1) + g*h*log(G)/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(g*h*log(G)))
    rubi.add(rule73)

    pattern74 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, g, F, G, h, e: Not(RationalQ(FullSimplify(g*h*log(G)/(d*e*log(F)))))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule74 = ReplacementRule(pattern74, lambda d, g, n, F, a, h, G, b, e, c, f, x : G**(h*(f + g*x))*(F**(e*(c + d*x))*b + a)**(n + S(1))*Hypergeometric2F1(S(1), n + S(1) + g*h*log(G)/(d*e*log(F)), S(1) + g*h*log(G)/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(a*g*h*log(G)))
    rubi.add(rule74)

    pattern75 = Pattern(Integral(G_**(u_*WC('h', S(1)))*(F_**(v_*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda v, x, u: LinearQ(List(u, v), x)), CustomConstraint(lambda v, x, u: Not(LinearMatchQ(List(u, v), x))))
    rule75 = ReplacementRule(pattern75, lambda x, n, F, a, h, G, b, e, v, u : Int(G**(h*ExpandToSum(u, x))*(F**(e*ExpandToSum(v, x))*b + a)**n, x))
    rubi.add(rule75)

    pattern76 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda H, x: FreeQ(H, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda t, x: FreeQ(t, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, s, g, G, d, F, f, h, r, a, c, H, n, e, b, x, t: RationalQ(m)))
    def With76(d, g, n, t, a, h, F, G, b, s, H, r, e, c, f, x):
        m = FullSimplify((g*h*log(G) + s*t*log(H))/(d*e*log(F)))
        return G**(-c*g*h/d + f*h)*H**(-c*s*t/d + r*t)*Denominator(m)*Subst(Int(x**(Numerator(m) + S(-1))*(a + b*x**Denominator(m))**n, x), x, F**(e*(c + d*x)/Denominator(m)))/(d*e*log(F))
    rule76 = ReplacementRule(pattern76, lambda d, g, n, t, a, h, F, G, b, s, H, r, e, c, f, x : With76(d, g, n, t, a, h, F, G, b, s, H, r, e, c, f, x))
    rubi.add(rule76)

    pattern77 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda H, x: FreeQ(H, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda t, x: FreeQ(t, x)), CustomConstraint(lambda d, g, n, F, h, G, e: ZeroQ(d*e*n*log(F) + g*h*log(G))), CustomConstraint(lambda n: IntegerQ(n)))
    rule77 = ReplacementRule(pattern77, lambda d, g, n, t, a, h, F, G, b, s, H, r, e, c, f, x : G**(h*(-c*g/d + f))*Int(H**(t*(r + s*x))*(b + F**(-e*(c + d*x))*a)**n, x))
    rubi.add(rule77)

    pattern78 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda H, x: FreeQ(H, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda t, x: FreeQ(t, x)), CustomConstraint(lambda d, g, t, F, G, h, s, H, e: Not(RationalQ(FullSimplify((g*h*log(G) + s*t*log(H))/(d*e*log(F)))))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule78 = ReplacementRule(pattern78, lambda d, g, n, t, a, h, F, G, b, s, H, r, e, c, f, x : Int(G**(f*h)*G**(g*h*x)*H**(r*t)*H**(s*t*x)*(F**(c*e)*F**(d*e*x)*b + a)**n, x))
    rubi.add(rule78)

    pattern79 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda H, x: FreeQ(H, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda t, x: FreeQ(t, x)), CustomConstraint(lambda d, g, t, F, G, h, s, H, e: Not(RationalQ(FullSimplify((g*h*log(G) + s*t*log(H))/(d*e*log(F)))))), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule79 = ReplacementRule(pattern79, lambda d, g, t, F, a, h, n, G, b, s, H, r, e, c, f, x : G**(h*(f + g*x))*H**(t*(r + s*x))*a**n*Hypergeometric2F1(-n, (g*h*log(G) + s*t*log(H))/(d*e*log(F)), S(1) + (g*h*log(G) + s*t*log(H))/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(g*h*log(G) + s*t*log(H)))
    rubi.add(rule79)

    pattern80 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda H, x: FreeQ(H, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda t, x: FreeQ(t, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, g, t, F, G, h, s, H, e: Not(RationalQ(FullSimplify((g*h*log(G) + s*t*log(H))/(d*e*log(F)))))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule80 = ReplacementRule(pattern80, lambda d, g, t, F, a, h, n, G, b, s, H, r, e, c, f, x : G**(h*(f + g*x))*H**(t*(r + s*x))*((F**(e*(c + d*x))*b + a)/a)**(-n)*(F**(e*(c + d*x))*b + a)**n*Hypergeometric2F1(-n, (g*h*log(G) + s*t*log(H))/(d*e*log(F)), S(1) + (g*h*log(G) + s*t*log(H))/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(g*h*log(G) + s*t*log(H)))
    rubi.add(rule80)

    pattern81 = Pattern(Integral(G_**(u_*WC('h', S(1)))*H_**(w_*WC('t', S(1)))*(F_**(v_*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda H, x: FreeQ(H, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda t, x: FreeQ(t, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda w, v, x, u: LinearQ(List(u, v, w), x)), CustomConstraint(lambda w, v, x, u: Not(LinearMatchQ(List(u, v, w), x))))
    rule81 = ReplacementRule(pattern81, lambda x, t, F, a, h, n, G, w, b, H, e, v, u : Int(G**(h*ExpandToSum(u, x))*H**(t*ExpandToSum(w, x))*(F**(e*ExpandToSum(v, x))*b + a)**n, x))
    rubi.add(rule81)

    pattern82 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule82 = ReplacementRule(pattern82, lambda d, n, F, a, p, b, e, c, x : -a*n*Int(x**(n + S(-1))*(F**(e*(c + d*x))*b + a*x**n)**p, x)/(b*d*e*log(F)) + (F**(e*(c + d*x))*b + a*x**n)**(p + S(1))/(b*d*e*(p + S(1))*log(F)))
    rubi.add(rule82)

    pattern83 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*x_**WC('m', S(1))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule83 = ReplacementRule(pattern83, lambda d, n, F, a, p, m, b, e, c, x : -a*n*Int(x**(m + n + S(-1))*(F**(e*(c + d*x))*b + a*x**n)**p, x)/(b*d*e*log(F)) - m*Int(x**(m + S(-1))*(F**(e*(c + d*x))*b + a*x**n)**(p + S(1)), x)/(b*d*e*(p + S(1))*log(F)) + x**m*(F**(e*(c + d*x))*b + a*x**n)**(p + S(1))/(b*d*e*(p + S(1))*log(F)))
    rubi.add(rule83)

    pattern84 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(F_**u_*WC('b', S(1)) + F_**v_*WC('c', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda v, u: ZeroQ(-S(2)*u + v)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda c, b, a: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda m: PositiveIntegerQ(m)), )
    def With84(x, g, F, a, m, b, c, v, f, u):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int((f + g*x)**m/(S(2)*F**u*c + b - q), x)/q - S(2)*c*Int((f + g*x)**m/(S(2)*F**u*c + b + q), x)/q
    rule84 = ReplacementRule(pattern84, lambda x, g, F, a, m, b, c, v, f, u : With84(x, g, F, a, m, b, c, v, f, u))
    rubi.add(rule84)

    pattern85 = Pattern(Integral(F_**u_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(F_**u_*WC('b', S(1)) + F_**v_*WC('c', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda v, u: ZeroQ(-S(2)*u + v)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda c, b, a: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda m: PositiveIntegerQ(m)), )
    def With85(x, g, F, a, m, b, c, v, f, u):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int(F**u*(f + g*x)**m/(S(2)*F**u*c + b - q), x)/q - S(2)*c*Int(F**u*(f + g*x)**m/(S(2)*F**u*c + b + q), x)/q
    rule85 = ReplacementRule(pattern85, lambda x, g, F, a, m, b, c, v, f, u : With85(x, g, F, a, m, b, c, v, f, u))
    rubi.add(rule85)

    pattern86 = Pattern(Integral((F_**u_*WC('i', S(1)) + h_)*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(F_**u_*WC('b', S(1)) + F_**v_*WC('c', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda v, u: ZeroQ(-S(2)*u + v)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda c, b, a: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda m: PositiveIntegerQ(m)), )
    def With86(i, x, g, F, a, h, m, b, c, v, f, u):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -(-i + (-b*i + S(2)*c*h)/q)*Int((f + g*x)**m/(S(2)*F**u*c + b + q), x) + (i + (-b*i + S(2)*c*h)/q)*Int((f + g*x)**m/(S(2)*F**u*c + b - q), x)
    rule86 = ReplacementRule(pattern86, lambda i, x, g, F, a, h, m, b, c, v, f, u : With86(i, x, g, F, a, h, m, b, c, v, f, u))
    rubi.add(rule86)

    pattern87 = Pattern(Integral(x_**WC('m', S(1))/(F_**v_*WC('b', S(1)) + F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('a', S(1))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d, v, x: ZeroQ(c + d*x + v)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))), )
    def With87(d, F, a, m, b, c, v, x):
        u = IntHide(1/(F**v*b + F**(c + d*x)*a), x)
        return -m*Int(u*x**(m + S(-1)), x) + u*x**m
    rule87 = ReplacementRule(pattern87, lambda d, F, a, m, b, c, v, x : With87(d, F, a, m, b, c, v, x))
    rubi.add(rule87)

    pattern88 = Pattern(Integral(u_/(F_**v_*WC('b', S(1)) + F_**w_*WC('c', S(1)) + a_), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda x, w: LinearQ(w, x)), CustomConstraint(lambda v, w: ZeroQ(v + w)), CustomConstraint(lambda w, v, x: If(RationalQ(Coefficient(v, x, S(1))), Greater(Coefficient(v, x, S(1)), S(0)), Less(LeafCount(v), LeafCount(w)))))
    rule88 = ReplacementRule(pattern88, lambda u, x, F, a, b, c, v, w : Int(F**v*u/(F**(S(2)*v)*b + F**v*a + c), x))
    rubi.add(rule88)

    pattern89 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule89 = ReplacementRule(pattern89, lambda d, g, n, F, a, b, e, c, x : Int(ExpandIntegrand(F**(g*(d + e*x)**n), 1/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule89)

    pattern90 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule90 = ReplacementRule(pattern90, lambda d, g, n, F, a, e, c, x : Int(ExpandIntegrand(F**(g*(d + e*x)**n), 1/(a + c*x**S(2)), x), x))
    rubi.add(rule90)

    pattern91 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))*u_**WC('m', S(1))/(c_*x_**S(2) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, u: PolynomialQ(u, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule91 = ReplacementRule(pattern91, lambda u, d, g, n, F, a, m, b, e, c, x : Int(ExpandIntegrand(F**(g*(d + e*x)**n), u**m/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule91)

    pattern92 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))*u_**WC('m', S(1))/(a_ + c_*x_**S(2)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, u: PolynomialQ(u, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule92 = ReplacementRule(pattern92, lambda u, d, g, n, F, a, m, e, c, x : Int(ExpandIntegrand(F**(g*(d + e*x)**n), u**m/(a + c*x**S(2)), x), x))
    rubi.add(rule92)

    pattern93 = Pattern(Integral(F_**((x_**S(4)*WC('b', S(1)) + WC('a', S(0)))/x_**S(2)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule93 = ReplacementRule(pattern93, lambda b, a, F, x : -sqrt(Pi)*Erf((-x**S(2)*sqrt(-b*log(F)) + sqrt(-a*log(F)))/x)*Exp(-S(2)*sqrt(-a*log(F))*sqrt(-b*log(F)))/(S(4)*sqrt(-b*log(F))) + sqrt(Pi)*Erf((x**S(2)*sqrt(-b*log(F)) + sqrt(-a*log(F)))/x)*Exp(S(2)*sqrt(-a*log(F))*sqrt(-b*log(F)))/(S(4)*sqrt(-b*log(F))))
    rubi.add(rule93)

    pattern94 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('m', S(1)) + exp(x_))**n_, x_), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: Less(n, S(0))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule94 = ReplacementRule(pattern94, lambda n, m, x : m*Int(x**(m + S(-1))*(x**m + exp(x))**n, x) + Int((x**m + exp(x))**(n + S(1)), x) - (x**m + exp(x))**(n + S(1))/(n + S(1)))
    rubi.add(rule94)

    pattern95 = Pattern(Integral(log(a_ + (F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a: PositiveQ(a)))
    rule95 = ReplacementRule(pattern95, lambda d, n, F, a, b, e, c, x : Subst(Int(log(a + b*x)/x, x), x, (F**(e*(c + d*x)))**n)/(d*e*n*log(F)))
    rubi.add(rule95)

    pattern96 = Pattern(Integral(log(a_ + (F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a: Not(PositiveQ(a))))
    rule96 = ReplacementRule(pattern96, lambda d, n, F, a, b, e, c, x : -b*d*e*n*Int(x*(F**(e*(c + d*x)))**n/(a + b*(F**(e*(c + d*x)))**n), x)*log(F) + x*log(a + b*(F**(e*(c + d*x)))**n))
    rubi.add(rule96)

    pattern97 = Pattern(Integral((F_**v_*WC('a', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule97 = ReplacementRule(pattern97, lambda x, n, F, a, v, u : F**(-n*v)*(F**v*a)**n*Int(F**(n*v)*u, x))
    rubi.add(rule97)

    pattern98 = Pattern(Integral(u_, x_), CustomConstraint(lambda x, u: FunctionOfExponentialQ(u, x)), )
    def With98(x, u):
        v = FunctionOfExponential(u, x)
        return v*Subst(Int(FunctionOfExponentialFunction(u, x)/x, x), x, v)/D(v, x)
    rule98 = ReplacementRule(pattern98, lambda x, u : With98(x, u))
    rubi.add(rule98)

    pattern99 = Pattern(Integral((F_**v_*WC('a', S(1)) + F_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda v, x, w: LinearQ(List(v, w), x)))
    rule99 = ReplacementRule(pattern99, lambda x, w, F, a, n, b, v, u : Int(F**(n*v)*u*(F**ExpandToSum(-v + w, x)*b + a)**n, x))
    rubi.add(rule99)

    pattern100 = Pattern(Integral((F_**v_*WC('a', S(1)) + G_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda v, x, w: LinearQ(List(v, w), x)))
    rule100 = ReplacementRule(pattern100, lambda x, w, F, a, G, n, b, v, u : Int(F**(n*v)*u*(a + b*exp(ExpandToSum(-v*log(F) + w*log(G), x)))**n, x))
    rubi.add(rule100)

    pattern101 = Pattern(Integral((F_**v_*WC('a', S(1)) + F_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n))), CustomConstraint(lambda v, x, w: LinearQ(List(v, w), x)))
    rule101 = ReplacementRule(pattern101, lambda x, w, F, a, n, b, v, u : F**(-n*v)*(F**v*a + F**w*b)**n*(F**ExpandToSum(-v + w, x)*b + a)**(-n)*Int(F**(n*v)*u*(F**ExpandToSum(-v + w, x)*b + a)**n, x))
    rubi.add(rule101)

    pattern102 = Pattern(Integral((F_**v_*WC('a', S(1)) + G_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n))), CustomConstraint(lambda v, x, w: LinearQ(List(v, w), x)))
    rule102 = ReplacementRule(pattern102, lambda x, w, F, a, G, n, b, v, u : F**(-n*v)*(a + b*exp(ExpandToSum(-v*log(F) + w*log(G), x)))**(-n)*(F**v*a + G**w*b)**n*Int(F**(n*v)*u*(a + b*exp(ExpandToSum(-v*log(F) + w*log(G), x)))**n, x))
    rubi.add(rule102)

    pattern103 = Pattern(Integral(F_**v_*G_**w_*WC('u', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda G, x: FreeQ(G, x)), CustomConstraint(lambda v, x, w: BinomialQ(v + w, x) | (PolynomialQ(v + w, x) & LessEqual(Exponent(v + w, x), S(2)))))
    rule103 = ReplacementRule(pattern103, lambda u, x, F, G, v, w : Int(u*NormalizeIntegrand(exp(v*log(F) + w*log(G)), x), x))
    rubi.add(rule103)

    pattern104 = Pattern(Integral(F_**u_*(v_ + w_)*WC('y', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda F, z, x, y, u, w: ZeroQ(-w*y + D(z, x))))
    def With104(u, x, F, v, y, w):
        z = v*y/(D(u, x)*log(F))
        return F**u*z
    rule104 = ReplacementRule(pattern104, lambda u, x, F, v, y, w : With104(u, x, F, v, y, w))
    rubi.add(rule104)

    pattern105 = Pattern(Integral(F_**u_*v_**WC('n', S(1))*w_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, u: PolynomialQ(u, x)), CustomConstraint(lambda v, x: PolynomialQ(v, x)), CustomConstraint(lambda x, w: PolynomialQ(w, x)), CustomConstraint(lambda v, w, F, z, x, u, n: Equal(Exponent(w, x), Exponent(z, x)) & ZeroQ(w*Coefficient(z, x, Exponent(z, x)) - z*Coefficient(w, x, Exponent(w, x)))))
    def With105(u, x, n, F, v, w):
        z = v*D(u, x)*log(F) + (n + S(1))*D(v, x)
        return F**u*v**(n + S(1))*Coefficient(w, x, Exponent(w, x))/Coefficient(z, x, Exponent(z, x))
    rule105 = ReplacementRule(pattern105, lambda u, x, n, F, v, w : With105(u, x, n, F, v, w))
    rubi.add(rule105)

    return rubi
