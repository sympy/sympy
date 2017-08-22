from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (Int, Set, With, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcCsch, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, FixSimplify, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, SimpFixFactor, Simp, SimpHelp, SplitProduct, SplitSum, UnsameQ)
    from sympy import Integral, S, sqrt
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols

    A_, B_, C_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, z_ = [WC(i) for i in 'ABCabcdefghijklmnpqrtuvswxz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, jn_, mn_, non2_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'n1', 'n2', 'n3',' Pq', 'Pm', ' Px', 'Qm', 'Qr', 'jn', 'mn', 'non2']]
    p, q, r, s, mn, gcd, P, Q = symbols('p q r s mn gcd P Q')

    _UseGamma = False

def binomial_products(rubi):

    pattern1 = Pattern(Integral((x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule1 = ReplacementRule(pattern1, lambda b, p, n, x : b**IntPart(p)*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(n*p), x))
    rubi.add(rule1)

    pattern2 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n: ZeroQ(p + S(1) + 1/n)))
    rule2 = ReplacementRule(pattern2, lambda p, x, a, b, n : x*(a + b*x**n)**(p + S(1))/a)
    rubi.add(rule2)

    pattern3 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n: NegativeIntegerQ(p + S(1) + 1/n)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule3 = ReplacementRule(pattern3, lambda p, x, a, b, n : -x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))) + (n*(p + S(1)) + S(1))*Int((a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))))
    rubi.add(rule3)

    pattern4 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: NonzeroQ(S(3)*n + S(1))))
    rule4 = ReplacementRule(pattern4, lambda b, n, a, x : Int(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n), x))
    rubi.add(rule4)

    pattern5 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(0))), CustomConstraint(lambda p: IntegerQ(p)))
    rule5 = ReplacementRule(pattern5, lambda p, x, a, b, n : Int(x**(n*p)*(a*x**(-n) + b)**p, x))
    rubi.add(rule5)

    pattern6 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)))
    rule6 = ReplacementRule(pattern6, lambda p, x, a, b, n : Int(ExpandIntegrand((a + b*x**n)**p, x), x))
    rubi.add(rule6)

    pattern7 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n: IntegerQ(S(2)*p) | (Equal(n, S(2)) & IntegerQ(S(3)*p)) | (Equal(n, S(2)) & IntegerQ(S(4)*p)) | Less(Denominator(p + 1/n), Denominator(p))))
    rule7 = ReplacementRule(pattern7, lambda p, x, a, b, n : a*n*p*Int((a + b*x**n)**(p + S(-1)), x)/(n*p + S(1)) + x*(a + b*x**n)**p/(n*p + S(1)))
    rubi.add(rule7)

    pattern8 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-5)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: PosQ(b/a)), CustomConstraint(lambda a: PositiveQ(a)))
    rule8 = ReplacementRule(pattern8, lambda b, a, x : S(2)*EllipticE(ArcTan(x*Rt(b/a, S(2)))/S(2), S(2))/(a**(S(5)/4)*Rt(b/a, S(2))))
    rubi.add(rule8)

    pattern9 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-5)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: PosQ(b/a)), CustomConstraint(lambda a: ~(PositiveQ(a))))
    rule9 = ReplacementRule(pattern9, lambda b, a, x : (S(1) + b*x**S(2)/a)**(S(1)/4)*Int((S(1) + b*x**S(2)/a)**(S(-5)/4), x)/(a*(a + b*x**S(2))**(S(1)/4)))
    rubi.add(rule9)

    pattern10 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-7)/6), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule10 = ReplacementRule(pattern10, lambda b, a, x : Subst(Int((-b*x**S(2) + S(1))**(S(-1)/3), x), x, x/Sqrt(a + b*x**S(2)))/((a/(a + b*x**S(2)))**(S(2)/3)*(a + b*x**S(2))**(S(2)/3)))
    rubi.add(rule10)

    pattern11 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n: IntegerQ(S(2)*p) | (Equal(n, S(2)) & IntegerQ(S(3)*p)) | (Equal(n, S(2)) & IntegerQ(S(4)*p)) | Less(Denominator(p + 1/n), Denominator(p))))
    rule11 = ReplacementRule(pattern11, lambda p, x, a, b, n : -x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))) + (n*(p + S(1)) + S(1))*Int((a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(1/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule12 = ReplacementRule(pattern12, lambda b, a, x : Int((-x*Rt(b, S(3)) + S(2)*Rt(a, S(3)))/(x**S(2)*Rt(b, S(3))**S(2) - x*Rt(a, S(3))*Rt(b, S(3)) + Rt(a, S(3))**S(2)), x)/(S(3)*Rt(a, S(3))**S(2)) + Int(1/(x*Rt(b, S(3)) + Rt(a, S(3))), x)/(S(3)*Rt(a, S(3))**S(2)))
    rubi.add(rule12)

    pattern13 = Pattern(Integral(1/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2) + S(-3)/2)), CustomConstraint(lambda a, b: PosQ(a/b)))
    rule13 = ReplacementRule(pattern13, lambda b, n, a, x : Module(List(Set(r, Numerator(Rt(a/b, n))), Set(s, Denominator(Rt(a/b, n))), k, u), CompoundExpression(Set(u, Int((r - s*x*Cos(Pi*(S(2)*k + S(-1))/n))/(r**S(2) - S(2)*r*s*x*Cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)), Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r*Int(1/(r + s*x), x)/(a*n))))
    rubi.add(rule13)

    pattern14 = Pattern(Integral(1/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2) + S(-3)/2)), CustomConstraint(lambda a, b: NegQ(a/b)))
    rule14 = ReplacementRule(pattern14, lambda b, n, a, x : Module(List(Set(r, Numerator(Rt(-a/b, n))), Set(s, Denominator(Rt(-a/b, n))), k, u), CompoundExpression(Set(u, Int((r + s*x*Cos(Pi*(S(2)*k + S(-1))/n))/(r**S(2) + S(2)*r*s*x*Cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)), Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r*Int(1/(r - s*x), x)/(a*n))))
    rubi.add(rule14)

    pattern15 = Pattern(Integral(1/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a, b: PosQ(a/b)), CustomConstraint(lambda a, b: PositiveQ(a) | PositiveQ(b)))
    rule15 = ReplacementRule(pattern15, lambda b, a, x : ArcTan(x*Rt(b, S(2))/Rt(a, S(2)))/(Rt(a, S(2))*Rt(b, S(2))))
    rubi.add(rule15)

    pattern16 = Pattern(Integral(1/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a, b: PosQ(a/b)), CustomConstraint(lambda a, b: NegativeQ(a) | NegativeQ(b)))
    rule16 = ReplacementRule(pattern16, lambda b, a, x : -ArcTan(x*Rt(-b, S(2))/Rt(-a, S(2)))/(Rt(-a, S(2))*Rt(-b, S(2))))
    rubi.add(rule16)

    pattern17 = Pattern(Integral(1/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a, b: PosQ(a/b)))
    rule17 = ReplacementRule(pattern17, lambda b, a, x : ArcTan(x/Rt(a/b, S(2)))*Rt(a/b, S(2))/a)
    rubi.add(rule17)

    pattern18 = Pattern(Integral(1/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a, b: NegQ(a/b)), CustomConstraint(lambda a, b: NegativeQ(b) | PositiveQ(a)))
    rule18 = ReplacementRule(pattern18, lambda b, a, x : ArcTanh(x*Rt(-b, S(2))/Rt(a, S(2)))/(Rt(a, S(2))*Rt(-b, S(2))))
    rubi.add(rule18)

    pattern19 = Pattern(Integral(1/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a, b: NegQ(a/b)), CustomConstraint(lambda a, b: NegativeQ(a) | PositiveQ(b)))
    rule19 = ReplacementRule(pattern19, lambda b, a, x : -ArcTanh(x*Rt(b, S(2))/Rt(-a, S(2)))/(Rt(-a, S(2))*Rt(b, S(2))))
    rubi.add(rule19)

    pattern20 = Pattern(Integral(1/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a, b: NegQ(a/b)))
    rule20 = ReplacementRule(pattern20, lambda b, a, x : ArcTanh(x/Rt(-a/b, S(2)))*Rt(-a/b, S(2))/a)
    rubi.add(rule20)

    pattern21 = Pattern(Integral(1/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(4) + S(-1)/2)), CustomConstraint(lambda a, b: PosQ(a/b)))
    rule21 = ReplacementRule(pattern21, lambda b, n, a, x : Module(List(Set(r, Numerator(Rt(a/b, n))), Set(s, Denominator(Rt(a/b, n))), k, u, v), CompoundExpression(Set(u, Int((r - s*x*Cos(Pi*(S(2)*k + S(-1))/n))/(r**S(2) - S(2)*r*s*x*Cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x) + Int((r + s*x*Cos(Pi*(S(2)*k + S(-1))/n))/(r**S(2) + S(2)*r*s*x*Cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)), Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**S(2)*Int(1/(r**S(2) + s**S(2)*x**S(2)), x)/(a*n))))
    rubi.add(rule21)

    pattern22 = Pattern(Integral(1/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(4) + S(-1)/2)), CustomConstraint(lambda a, b: NegQ(a/b)))
    rule22 = ReplacementRule(pattern22, lambda b, n, a, x : Module(List(Set(r, Numerator(Rt(-a/b, n))), Set(s, Denominator(Rt(-a/b, n))), k, u), CompoundExpression(Set(u, Int((r - s*x*Cos(S(2)*Pi*k/n))/(r**S(2) - S(2)*r*s*x*Cos(S(2)*Pi*k/n) + s**S(2)*x**S(2)), x) + Int((r + s*x*Cos(S(2)*Pi*k/n))/(r**S(2) + S(2)*r*s*x*Cos(S(2)*Pi*k/n) + s**S(2)*x**S(2)), x)), Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**S(2)*Int(1/(r**S(2) - s**S(2)*x**S(2)), x)/(a*n))))
    rubi.add(rule22)

    pattern23 = Pattern(Integral(1/(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a, b: PositiveQ(a/b) | (AtomQ(SplitProduct(SumBaseQ, a)) & AtomQ(SplitProduct(SumBaseQ, b)) & PosQ(a/b))))
    rule23 = ReplacementRule(pattern23, lambda b, a, x : With(List(Set(r, Numerator(Rt(a/b, S(2)))), Set(s, Denominator(Rt(a/b, S(2))))), Int((r - s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*r) + Int((r + s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*r)))
    rubi.add(rule23)

    pattern24 = Pattern(Integral(1/(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a, b: ~(PositiveQ(a/b))))
    rule24 = ReplacementRule(pattern24, lambda b, a, x : With(List(Set(r, Numerator(Rt(-a/b, S(2)))), Set(s, Denominator(Rt(-a/b, S(2))))), r*Int(1/(r - s*x**S(2)), x)/(S(2)*a) + r*Int(1/(r + s*x**S(2)), x)/(S(2)*a)))
    rubi.add(rule24)

    pattern25 = Pattern(Integral(1/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(4) + S(-1))), CustomConstraint(lambda a, b: PositiveQ(a/b)))
    rule25 = ReplacementRule(pattern25, lambda b, n, a, x : With(List(Set(r, Numerator(Rt(a/b, S(4)))), Set(s, Denominator(Rt(a/b, S(4))))), r*Int((r*Sqrt(S(2)) - s*x**(n/S(4)))/(r**S(2) - r*s*x**(n/S(4))*Sqrt(S(2)) + s**S(2)*x**(n/S(2))), x)/(S(2)*a*Sqrt(S(2))) + r*Int((r*Sqrt(S(2)) + s*x**(n/S(4)))/(r**S(2) + r*s*x**(n/S(4))*Sqrt(S(2)) + s**S(2)*x**(n/S(2))), x)/(S(2)*a*Sqrt(S(2)))))
    rubi.add(rule25)

    pattern26 = Pattern(Integral(1/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(4) + S(-1))), CustomConstraint(lambda a, b: ~(PositiveQ(a/b))))
    rule26 = ReplacementRule(pattern26, lambda b, n, a, x : With(List(Set(r, Numerator(Rt(-a/b, S(2)))), Set(s, Denominator(Rt(-a/b, S(2))))), r*Int(1/(r - s*x**(n/S(2))), x)/(S(2)*a) + r*Int(1/(r + s*x**(n/S(2))), x)/(S(2)*a)))
    rubi.add(rule26)

    pattern27 = Pattern(Integral(1/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: PositiveQ(a)), CustomConstraint(lambda b: PosQ(b)))
    rule27 = ReplacementRule(pattern27, lambda b, a, x : ArcSinh(x*Rt(b, S(2))/Sqrt(a))/Rt(b, S(2)))
    rubi.add(rule27)

    pattern28 = Pattern(Integral(1/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: PositiveQ(a)), CustomConstraint(lambda b: NegQ(b)))
    rule28 = ReplacementRule(pattern28, lambda b, a, x : ArcSin(x*Rt(-b, S(2))/Sqrt(a))/Rt(-b, S(2)))
    rubi.add(rule28)

    pattern29 = Pattern(Integral(1/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: ~(PositiveQ(a))))
    rule29 = ReplacementRule(pattern29, lambda b, a, x : Subst(Int(1/(-b*x**S(2) + S(1)), x), x, x/Sqrt(a + b*x**S(2))))
    rubi.add(rule29)

    pattern30 = Pattern(Integral(1/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: PosQ(a)))
    rule30 = ReplacementRule(pattern30, lambda b, a, x : With(List(Set(r, Numer(Rt(b/a, S(3)))), Set(s, Denom(Rt(b/a, S(3))))), S(2)*S(3)**(S(3)/4)*(r*x + s)*EllipticF(ArcSin((r*x + s*(-Sqrt(S(3)) + S(1)))/(r*x + s*(Sqrt(S(3)) + S(1)))), -S(4)*Sqrt(S(3)) + S(-7))*Sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(Sqrt(S(3)) + S(1)))**S(2))*Sqrt(Sqrt(S(3)) + S(2))/(S(3)*r*Sqrt(s*(r*x + s)/(r*x + s*(Sqrt(S(3)) + S(1)))**S(2))*Sqrt(a + b*x**S(3)))))
    rubi.add(rule30)

    pattern31 = Pattern(Integral(1/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: NegQ(a)))
    rule31 = ReplacementRule(pattern31, lambda b, a, x : With(List(Set(r, Numer(Rt(b/a, S(3)))), Set(s, Denom(Rt(b/a, S(3))))), S(2)*S(3)**(S(3)/4)*(r*x + s)*EllipticF(ArcSin((r*x + s*(Sqrt(S(3)) + S(1)))/(r*x + s*(-Sqrt(S(3)) + S(1)))), S(4)*Sqrt(S(3)) + S(-7))*Sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(-Sqrt(S(3)) + S(1)))**S(2))*Sqrt(-Sqrt(S(3)) + S(2))/(S(3)*r*Sqrt(-s*(r*x + s)/(r*x + s*(-Sqrt(S(3)) + S(1)))**S(2))*Sqrt(a + b*x**S(3)))))
    rubi.add(rule31)

    pattern32 = Pattern(Integral(1/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule32 = ReplacementRule(pattern32, lambda b, a, x : With(List(Set(q, Rt(b/a, S(4)))), (q**S(2)*x**S(2) + S(1))*EllipticF(S(2)*ArcTan(q*x), S(1)/2)*Sqrt((a + b*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))/(S(2)*q*Sqrt(a + b*x**S(4)))))
    rubi.add(rule32)

    pattern33 = Pattern(Integral(1/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: NegQ(b/a)), CustomConstraint(lambda a: PositiveQ(a)))
    rule33 = ReplacementRule(pattern33, lambda b, a, x : EllipticF(ArcSin(x*Rt(-b, S(4))/Rt(a, S(4))), S(-1))/(Rt(a, S(4))*Rt(-b, S(4))))
    rubi.add(rule33)

    pattern34 = Pattern(Integral(1/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: NegativeQ(a)), CustomConstraint(lambda b: PositiveQ(b)))
    rule34 = ReplacementRule(pattern34, lambda b, a, x : With(List(Set(q, Rt(-a*b, S(2)))), Condition(EllipticF(ArcSin(x/Sqrt((a + q*x**S(2))/(S(2)*q))), S(1)/2)*Sqrt((a + q*x**S(2))/q)*Sqrt(-a + q*x**S(2))/(Sqrt(S(2))*Sqrt(-a)*Sqrt(a + b*x**S(4))), IntegerQ(q))))
    rubi.add(rule34)

    pattern35 = Pattern(Integral(1/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: NegativeQ(a)), CustomConstraint(lambda b: PositiveQ(b)))
    rule35 = ReplacementRule(pattern35, lambda b, a, x : With(List(Set(q, Rt(-a*b, S(2)))), EllipticF(ArcSin(x/Sqrt((a + q*x**S(2))/(S(2)*q))), S(1)/2)*Sqrt((a + q*x**S(2))/q)*Sqrt((a - q*x**S(2))/(a + q*x**S(2)))/(Sqrt(S(2))*Sqrt(a/(a + q*x**S(2)))*Sqrt(a + b*x**S(4)))))
    rubi.add(rule35)

    pattern36 = Pattern(Integral(1/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: NegQ(b/a)), CustomConstraint(lambda a: ~(PositiveQ(a))))
    rule36 = ReplacementRule(pattern36, lambda b, a, x : Int(1/Sqrt(S(1) + b*x**S(4)/a), x)*Sqrt(S(1) + b*x**S(4)/a)/Sqrt(a + b*x**S(4)))
    rubi.add(rule36)

    pattern37 = Pattern(Integral(1/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule37 = ReplacementRule(pattern37, lambda b, a, x : With(List(Set(r, Numer(Rt(b/a, S(3)))), Set(s, Denom(Rt(b/a, S(3))))), S(3)**(S(3)/4)*x*(r*x**S(2) + s)*EllipticF(ArcCos((r*x**S(2)*(-Sqrt(S(3)) + S(1)) + s)/(r*x**S(2)*(Sqrt(S(3)) + S(1)) + s)), Sqrt(S(3))/S(4) + S(1)/2)*Sqrt((r**S(2)*x**S(4) - r*s*x**S(2) + s**S(2))/(r*x**S(2)*(Sqrt(S(3)) + S(1)) + s)**S(2))/(S(6)*s*Sqrt(r*x**S(2)*(r*x**S(2) + s)/(r*x**S(2)*(Sqrt(S(3)) + S(1)) + s)**S(2))*Sqrt(a + b*x**S(6)))))
    rubi.add(rule37)

    pattern38 = Pattern(Integral(1/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule38 = ReplacementRule(pattern38, lambda b, a, x : Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/Sqrt(a + b*x**S(8)), x)/S(2) + Int((x**S(2)*Rt(b/a, S(4)) + S(1))/Sqrt(a + b*x**S(8)), x)/S(2))
    rubi.add(rule38)

    pattern39 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule39 = ReplacementRule(pattern39, lambda b, a, x : -a*Int((a + b*x**S(2))**(S(-5)/4), x) + S(2)*x/(a + b*x**S(2))**(S(1)/4))
    rubi.add(rule39)

    pattern40 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: NegQ(b/a)), CustomConstraint(lambda a: PositiveQ(a)))
    rule40 = ReplacementRule(pattern40, lambda b, a, x : S(2)*EllipticE(ArcSin(x*Rt(-b/a, S(2)))/S(2), S(2))/(a**(S(1)/4)*Rt(-b/a, S(2))))
    rubi.add(rule40)

    pattern41 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: NegQ(b/a)), CustomConstraint(lambda a: ~(PositiveQ(a))))
    rule41 = ReplacementRule(pattern41, lambda b, a, x : (S(1) + b*x**S(2)/a)**(S(1)/4)*Int((S(1) + b*x**S(2)/a)**(S(-1)/4), x)/(a + b*x**S(2))**(S(1)/4))
    rubi.add(rule41)

    pattern42 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: PositiveQ(a)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule42 = ReplacementRule(pattern42, lambda b, a, x : S(2)*EllipticF(ArcTan(x*Rt(b/a, S(2)))/S(2), S(2))/(a**(S(3)/4)*Rt(b/a, S(2))))
    rubi.add(rule42)

    pattern43 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: PositiveQ(a)), CustomConstraint(lambda b, a: NegQ(b/a)))
    rule43 = ReplacementRule(pattern43, lambda b, a, x : S(2)*EllipticF(ArcSin(x*Rt(-b/a, S(2)))/S(2), S(2))/(a**(S(3)/4)*Rt(-b/a, S(2))))
    rubi.add(rule43)

    pattern44 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: ~(PositiveQ(a))))
    rule44 = ReplacementRule(pattern44, lambda b, a, x : (S(1) + b*x**S(2)/a)**(S(3)/4)*Int((S(1) + b*x**S(2)/a)**(S(-3)/4), x)/(a + b*x**S(2))**(S(3)/4))
    rubi.add(rule44)

    pattern45 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/3), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule45 = ReplacementRule(pattern45, lambda b, a, x : S(3)*Sqrt(b*x**S(2))*Subst(Int(x/Sqrt(-a + x**S(3)), x), x, (a + b*x**S(2))**(S(1)/3))/(S(2)*b*x))
    rubi.add(rule45)

    pattern46 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-2)/3), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule46 = ReplacementRule(pattern46, lambda b, a, x : S(3)*Sqrt(b*x**S(2))*Subst(Int(1/Sqrt(-a + x**S(3)), x), x, (a + b*x**S(2))**(S(1)/3))/(S(2)*b*x))
    rubi.add(rule46)

    pattern47 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**(S(-3)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule47 = ReplacementRule(pattern47, lambda b, a, x : x**S(3)*(a/(b*x**S(4)) + S(1))**(S(3)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(3)/4)), x)/(a + b*x**S(4))**(S(3)/4))
    rubi.add(rule47)

    pattern48 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/6), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule48 = ReplacementRule(pattern48, lambda b, a, x : -a*Int((a + b*x**S(2))**(S(-7)/6), x)/S(2) + S(3)*x/(S(2)*(a + b*x**S(2))**(S(1)/6)))
    rubi.add(rule48)

    pattern49 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(S(-1), p, S(0))), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda p, n: IntegerQ(p + 1/n)))
    rule49 = ReplacementRule(pattern49, lambda p, x, a, b, n : a**(p + 1/n)*Subst(Int((-b*x**n + S(1))**(-p + S(-1) - S(1)/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule49)

    pattern50 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(S(-1), p, S(0))), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda p, n: Less(Denominator(p + 1/n), Denominator(p))))
    rule50 = ReplacementRule(pattern50, lambda p, x, a, b, n : (a/(a + b*x**n))**(p + 1/n)*(a + b*x**n)**(p + 1/n)*Subst(Int((-b*x**n + S(1))**(-p + S(-1) - S(1)/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule50)

    pattern51 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule51 = ReplacementRule(pattern51, lambda p, x, a, b, n : -Subst(Int((a + b*x**(-n))**p/x**S(2), x), x, 1/x))
    rubi.add(rule51)

    pattern52 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: FractionQ(n)))
    rule52 = ReplacementRule(pattern52, lambda p, x, a, b, n : With(List(Set(k, Denominator(n))), k*Subst(Int(x**(k + S(-1))*(a + b*x**(k*n))**p, x), x, x**(1/k))))
    rubi.add(rule52)

    pattern53 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule53 = ReplacementRule(pattern53, lambda p, x, a, b, n : Int(ExpandIntegrand((a + b*x**n)**p, x), x))
    rubi.add(rule53)

    pattern54 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: ~(PositiveIntegerQ(p))), CustomConstraint(lambda n: ~(IntegerQ(1/n))), CustomConstraint(lambda p, n: ~(NegativeIntegerQ(Simplify(p + 1/n)))), CustomConstraint(lambda a, p: IntegerQ(p) | PositiveQ(a)))
    rule54 = ReplacementRule(pattern54, lambda p, x, a, b, n : a**p*x*Hypergeometric2F1(-p, 1/n, S(1) + 1/n, -b*x**n/a))
    rubi.add(rule54)

    pattern55 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: ~(PositiveIntegerQ(p))), CustomConstraint(lambda n: ~(IntegerQ(1/n))), CustomConstraint(lambda p, n: ~(NegativeIntegerQ(Simplify(p + 1/n)))), CustomConstraint(lambda a, p: ~(IntegerQ(p) | PositiveQ(a))))
    rule55 = ReplacementRule(pattern55, lambda p, x, a, b, n : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((S(1) + b*x**n/a)**p, x))
    rubi.add(rule55)

    pattern56 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda x, u: NonzeroQ(u - x)))
    rule56 = ReplacementRule(pattern56, lambda p, x, a, u, b, n : Subst(Int((a + b*x**n)**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule56)

    pattern57 = Pattern(Integral((x_**n_*WC('b1', S(1)) + WC('a1', S(0)))**WC('p', S(1))*(x_**n_*WC('b2', S(1)) + WC('a2', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p, a1, a2: IntegerQ(p) | (PositiveQ(a1) & PositiveQ(a2))))
    rule57 = ReplacementRule(pattern57, lambda b2, x, n, b1, a1, p, a2 : Int((a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule57)

    pattern58 = Pattern(Integral((a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n: IntegerQ(S(2)*p) | Less(Denominator(p + 1/n), Denominator(p))))
    rule58 = ReplacementRule(pattern58, lambda b2, n, x, b1, a1, p, a2 : S(2)*a1*a2*n*p*Int((a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(S(2)*n*p + S(1)) + x*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(S(2)*n*p + S(1)))
    rubi.add(rule58)

    pattern59 = Pattern(Integral((a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n: IntegerQ(S(2)*p) | Less(Denominator(p + 1/n), Denominator(p))))
    rule59 = ReplacementRule(pattern59, lambda p, b2, x, b1, a1, n, a2 : -x*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*n*(p + S(1))) + (S(2)*n*(p + S(1)) + S(1))*Int((a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))))
    rubi.add(rule59)

    pattern60 = Pattern(Integral((a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: NegativeIntegerQ(S(2)*n)))
    rule60 = ReplacementRule(pattern60, lambda p, b2, x, b1, a1, n, a2 : -Subst(Int((a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p/x**S(2), x), x, 1/x))
    rubi.add(rule60)

    pattern61 = Pattern(Integral((a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: FractionQ(S(2)*n)))
    rule61 = ReplacementRule(pattern61, lambda p, b2, x, b1, a1, n, a2 : With(List(Set(k, Denominator(S(2)*n))), k*Subst(Int(x**(k + S(-1))*(a1 + b1*x**(k*n))**p*(a2 + b2*x**(k*n))**p, x), x, x**(1/k))))
    rubi.add(rule61)

    pattern62 = Pattern(Integral((x_**n_*WC('b1', S(1)) + WC('a1', S(0)))**p_*(x_**n_*WC('b2', S(1)) + WC('a2', S(0)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule62 = ReplacementRule(pattern62, lambda p, b2, x, b1, a1, n, a2 : (a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p))*Int((a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule62)

    pattern63 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p, a1, a2: IntegerQ(p) | (PositiveQ(a1) & PositiveQ(a2))))
    rule63 = ReplacementRule(pattern63, lambda p, b2, x, n, b1, a1, c, m, a2 : Int((c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule63)

    pattern64 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, m: IntegerQ(m) | PositiveQ(c)), CustomConstraint(lambda n, m: IntegerQ(Simplify((m + S(1))/n))))
    rule64 = ReplacementRule(pattern64, lambda p, x, n, b, c, m : b**(-Simplify((m + S(1))/n) + S(1))*c**m*Subst(Int((b*x)**(p + Simplify((m + S(1))/n) + S(-1)), x), x, x**n)/n)
    rubi.add(rule64)

    pattern65 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, m: IntegerQ(m) | PositiveQ(c)), CustomConstraint(lambda n, m: ~(IntegerQ(Simplify((m + S(1))/n)))))
    rule65 = ReplacementRule(pattern65, lambda p, x, b, c, n, m : b**IntPart(p)*c**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p), x))
    rubi.add(rule65)

    pattern66 = Pattern(Integral((c_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule66 = ReplacementRule(pattern66, lambda p, x, b, c, n, m : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(b*x**n)**p, x))
    rubi.add(rule66)

    pattern67 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: NegQ(n)))
    rule67 = ReplacementRule(pattern67, lambda p, x, a, b, n, m : Int(x**(m + n*p)*(a*x**(-n) + b)**p, x))
    rubi.add(rule67)

    pattern68 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, m: ZeroQ(p + S(1) + (m + S(1))/n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule68 = ReplacementRule(pattern68, lambda p, x, n, a, b, c, m : (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))))
    rubi.add(rule68)

    pattern69 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p, n, m: ZeroQ(p + S(1) + (m + S(1))/(S(2)*n))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule69 = ReplacementRule(pattern69, lambda p, b2, x, n, b1, a1, c, m, a2 : (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))))
    rubi.add(rule69)

    pattern70 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: IntegerQ(Simplify((m + S(1))/n))))
    rule70 = ReplacementRule(pattern70, lambda p, x, a, b, n, m : Subst(Int(x**(Simplify((m + S(1))/n) + S(-1))*(a + b*x)**p, x), x, x**n)/n)
    rubi.add(rule70)

    pattern71 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n, m: IntegerQ(Simplify((m + S(1))/(S(2)*n)))))
    rule71 = ReplacementRule(pattern71, lambda p, b2, x, b1, a1, n, m, a2 : Subst(Int(x**(Simplify((m + S(1))/n) + S(-1))*(a1 + b1*x)**p*(a2 + b2*x)**p, x), x, x**n)/n)
    rubi.add(rule71)

    pattern72 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: IntegerQ(Simplify((m + S(1))/n))))
    rule72 = ReplacementRule(pattern72, lambda p, x, n, a, b, c, m : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule72)

    pattern73 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n, m: IntegerQ(Simplify((m + S(1))/(S(2)*n)))))
    rule73 = ReplacementRule(pattern73, lambda p, b2, x, n, b1, a1, c, m, a2 : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule73)

    pattern74 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule74 = ReplacementRule(pattern74, lambda x, n, a, b, c, p, m : Int(ExpandIntegrand((c*x)**m*(a + b*x**n)**p, x), x))
    rubi.add(rule74)

    pattern75 = Pattern(Integral(x_**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, m: NegativeIntegerQ(Simplify(p + S(1) + (m + S(1))/n))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule75 = ReplacementRule(pattern75, lambda p, x, a, b, n, m : -b*(m + n*(p + S(1)) + S(1))*Int(x**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + x**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*(m + S(1))))
    rubi.add(rule75)

    pattern76 = Pattern(Integral(x_**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p, n, m: NegativeIntegerQ(Simplify(p + S(1) + (m + S(1))/(S(2)*n)))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule76 = ReplacementRule(pattern76, lambda p, b2, x, b1, a1, n, m, a2 : -b1*b2*(m + S(2)*n*(p + S(1)) + S(1))*Int(x**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + x**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*(m + S(1))))
    rubi.add(rule76)

    pattern77 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, m: NegativeIntegerQ(Simplify(p + S(1) + (m + S(1))/n))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule77 = ReplacementRule(pattern77, lambda p, x, n, a, b, c, m : (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))))
    rubi.add(rule77)

    pattern78 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p, n, m: NegativeIntegerQ(Simplify(p + S(1) + (m + S(1))/(S(2)*n)))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule78 = ReplacementRule(pattern78, lambda p, b2, x, n, b1, a1, c, m, a2 : (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))))
    rubi.add(rule78)

    pattern79 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule79 = ReplacementRule(pattern79, lambda p, x, a, b, n, m : With(List(Set(k, GCD(m + S(1), n))), Condition(Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p, x), x, x**k)/k, Unequal(k, S(1)))))
    rubi.add(rule79)

    pattern80 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule80 = ReplacementRule(pattern80, lambda p, b2, x, b1, a1, n, m, a2 : With(List(Set(k, GCD(m + S(1), S(2)*n))), Condition(Subst(Int(x**(S(-1) + (m + S(1))/k)*(a1 + b1*x**(n/k))**p*(a2 + b2*x**(n/k))**p, x), x, x**k)/k, Unequal(k, S(1)))))
    rubi.add(rule80)

    pattern81 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p, n, m: ~(NegativeIntegerQ((m + n*p + n + S(1))/n))), CustomConstraint(lambda p, n, x, a, b, c, m: IntBinomialQ(a, b, c, n, m, p, x)))
    rule81 = ReplacementRule(pattern81, lambda p, x, n, a, b, c, m : -b*c**(-n)*n*p*Int((c*x)**(m + n)*(a + b*x**n)**(p + S(-1)), x)/(m + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + S(1))))
    rubi.add(rule81)

    pattern82 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n, m: NonzeroQ(m + S(2)*n*p + S(1))), CustomConstraint(lambda p, b2, n, x, b1, a1, c, m, a2: IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x)))
    rule82 = ReplacementRule(pattern82, lambda p, b2, x, n, b1, a1, c, m, a2 : S(2)*a1*a2*n*p*Int((c*x)**m*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(2)*n*p + S(1)) + (c*x)**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(c*(m + S(2)*n*p + S(1))))
    rubi.add(rule82)

    pattern83 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n*p + S(1))), CustomConstraint(lambda p, n, x, a, b, c, m: IntBinomialQ(a, b, c, n, m, p, x)))
    rule83 = ReplacementRule(pattern83, lambda p, x, n, a, b, c, m : a*n*p*Int((c*x)**m*(a + b*x**n)**(p + S(-1)), x)/(m + n*p + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + n*p + S(1))))
    rubi.add(rule83)

    pattern84 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule84 = ReplacementRule(pattern84, lambda b, a, x : x*(a/(b*x**S(4)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(5)/4)), x)/(b*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule84)

    pattern85 = Pattern(Integral(x_**m_/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: PosQ(b/a)), CustomConstraint(lambda m: PositiveIntegerQ(m/S(4) + S(-1)/2)))
    rule85 = ReplacementRule(pattern85, lambda b, a, x, m : -a*(m + S(-3))*Int(x**(m + S(-4))/(a + b*x**S(4))**(S(5)/4), x)/(b*(m + S(-4))) + x**(m + S(-3))/(b*(a + b*x**S(4))**(S(1)/4)*(m + S(-4))))
    rubi.add(rule85)

    pattern86 = Pattern(Integral(x_**m_/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: PosQ(b/a)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(4) + S(-1)/2)))
    rule86 = ReplacementRule(pattern86, lambda b, a, x, m : -b*m*Int(x**(m + S(4))/(a + b*x**S(4))**(S(5)/4), x)/(a*(m + S(1))) + x**(m + S(1))/(a*(a + b*x**S(4))**(S(1)/4)*(m + S(1))))
    rubi.add(rule86)

    pattern87 = Pattern(Integral(sqrt(x_*WC('c', S(1)))/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule87 = ReplacementRule(pattern87, lambda b, c, x, a : (a/(b*x**S(2)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(2)*(a/(b*x**S(2)) + S(1))**(S(5)/4)), x)*Sqrt(c*x)/(b*(a + b*x**S(2))**(S(1)/4)))
    rubi.add(rule87)

    pattern88 = Pattern(Integral((x_*WC('c', S(1)))**m_/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, a: PosQ(b/a)), CustomConstraint(lambda m: IntegerQ(S(2)*m)), CustomConstraint(lambda m: Greater(m, S(3)/2)))
    rule88 = ReplacementRule(pattern88, lambda x, a, b, c, m : -S(2)*a*c**S(2)*(m + S(-1))*Int((c*x)**(m + S(-2))/(a + b*x**S(2))**(S(5)/4), x)/(b*(S(2)*m + S(-3))) + S(2)*c*(c*x)**(m + S(-1))/(b*(a + b*x**S(2))**(S(1)/4)*(S(2)*m + S(-3))))
    rubi.add(rule88)

    pattern89 = Pattern(Integral((x_*WC('c', S(1)))**m_/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, a: PosQ(b/a)), CustomConstraint(lambda m: IntegerQ(S(2)*m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule89 = ReplacementRule(pattern89, lambda x, a, b, c, m : -b*(S(2)*m + S(1))*Int((c*x)**(m + S(2))/(a + b*x**S(2))**(S(5)/4), x)/(S(2)*a*c**S(2)*(m + S(1))) + (c*x)**(m + S(1))/(a*c*(a + b*x**S(2))**(S(1)/4)*(m + S(1))))
    rubi.add(rule89)

    pattern90 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: NegQ(b/a)))
    rule90 = ReplacementRule(pattern90, lambda b, a, x : -Int(S(1)/(x**S(2)*(a + b*x**S(4))**(S(1)/4)), x)/b - S(1)/(b*x*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule90)

    pattern91 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, m: Greater(m + S(1), n)), CustomConstraint(lambda p, n, m: ~(NegativeIntegerQ((m + n*(p + S(1)) + S(1))/n))), CustomConstraint(lambda p, n, x, a, b, c, m: IntBinomialQ(a, b, c, n, m, p, x)))
    rule91 = ReplacementRule(pattern91, lambda p, x, n, a, b, c, m : -c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**(p + S(1)), x)/(b*n*(p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))))
    rubi.add(rule91)

    pattern92 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, m: Greater(m + S(1), S(2)*n)), CustomConstraint(lambda p, n, m: ~(NegativeIntegerQ((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*n)))), CustomConstraint(lambda p, b2, n, x, b1, a1, c, m, a2: IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x)))
    rule92 = ReplacementRule(pattern92, lambda p, b2, x, n, b1, a1, c, m, a2 : -c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*b1*b2*n*(p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*b1*b2*n*(p + S(1))))
    rubi.add(rule92)

    pattern93 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n, x, a, b, c, m: IntBinomialQ(a, b, c, n, m, p, x)))
    rule93 = ReplacementRule(pattern93, lambda p, x, n, a, b, c, m : (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))))
    rubi.add(rule93)

    pattern94 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, b2, n, x, b1, a1, c, m, a2: IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x)))
    rule94 = ReplacementRule(pattern94, lambda p, b2, x, n, b1, a1, c, m, a2 : (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))))
    rubi.add(rule94)

    pattern95 = Pattern(Integral(x_/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule95 = ReplacementRule(pattern95, lambda b, a, x : Int((x*Rt(b, S(3)) + Rt(a, S(3)))/(x**S(2)*Rt(b, S(3))**S(2) - x*Rt(a, S(3))*Rt(b, S(3)) + Rt(a, S(3))**S(2)), x)/(S(3)*Rt(a, S(3))*Rt(b, S(3))) - Int(1/(x*Rt(b, S(3)) + Rt(a, S(3))), x)/(S(3)*Rt(a, S(3))*Rt(b, S(3))))
    rubi.add(rule95)

    pattern96 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2) + S(-1)/2)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n, m: Less(m, n + S(-1))), CustomConstraint(lambda a, b: PosQ(a/b)))
    rule96 = ReplacementRule(pattern96, lambda x, a, b, n, m : Module(List(Set(r, Numerator(Rt(a/b, n))), Set(s, Denominator(Rt(a/b, n))), k, u), CompoundExpression(Set(u, Int((r*Cos(Pi*m*(S(2)*k + S(-1))/n) - s*x*Cos(Pi*(S(2)*k + S(-1))*(m + S(1))/n))/(r**S(2) - S(2)*r*s*x*Cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)), Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) - s**(-m)*(-r)**(m + S(1))*Int(1/(r + s*x), x)/(a*n))))
    rubi.add(rule96)

    pattern97 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n/S(2) + S(-1)/2)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n, m: Less(m, n + S(-1))), CustomConstraint(lambda a, b: NegQ(a/b)))
    rule97 = ReplacementRule(pattern97, lambda x, a, b, n, m : Module(List(Set(r, Numerator(Rt(-a/b, n))), Set(s, Denominator(Rt(-a/b, n))), k, u), CompoundExpression(Set(u, Int((r*Cos(Pi*m*(S(2)*k + S(-1))/n) + s*x*Cos(Pi*(S(2)*k + S(-1))*(m + S(1))/n))/(r**S(2) + S(2)*r*s*x*Cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)), -Dist(S(2)*s**(-m)*(-r)**(m + S(1))/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r**(m + S(1))*s**(-m)*Int(1/(r - s*x), x)/(a*n))))
    rubi.add(rule97)

    pattern98 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n/S(4) + S(-1)/2)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n, m: Less(m, n + S(-1))), CustomConstraint(lambda a, b: PosQ(a/b)))
    rule98 = ReplacementRule(pattern98, lambda x, a, b, n, m : Module(List(Set(r, Numerator(Rt(a/b, n))), Set(s, Denominator(Rt(a/b, n))), k, u), CompoundExpression(Set(u, Int((r*Cos(Pi*m*(S(2)*k + S(-1))/n) - s*x*Cos(Pi*(S(2)*k + S(-1))*(m + S(1))/n))/(r**S(2) - S(2)*r*s*x*Cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x) + Int((r*Cos(Pi*m*(S(2)*k + S(-1))/n) + s*x*Cos(Pi*(S(2)*k + S(-1))*(m + S(1))/n))/(r**S(2) + S(2)*r*s*x*Cos(Pi*(S(2)*k + S(-1))/n) + s**S(2)*x**S(2)), x)), S(2)*(S(-1))**(m/S(2))*r**(m + S(2))*s**(-m)*Int(1/(r**S(2) + s**S(2)*x**S(2)), x)/(a*n) + Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x))))
    rubi.add(rule98)

    pattern99 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n/S(4) + S(-1)/2)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n, m: Less(m, n + S(-1))), CustomConstraint(lambda a, b: NegQ(a/b)))
    rule99 = ReplacementRule(pattern99, lambda x, a, b, n, m : Module(List(Set(r, Numerator(Rt(-a/b, n))), Set(s, Denominator(Rt(-a/b, n))), k, u), CompoundExpression(Set(u, Int((r*Cos(S(2)*Pi*k*m/n) - s*x*Cos(S(2)*Pi*k*(m + S(1))/n))/(r**S(2) - S(2)*r*s*x*Cos(S(2)*Pi*k/n) + s**S(2)*x**S(2)), x) + Int((r*Cos(S(2)*Pi*k*m/n) + s*x*Cos(S(2)*Pi*k*(m + S(1))/n))/(r**S(2) + S(2)*r*s*x*Cos(S(2)*Pi*k/n) + s**S(2)*x**S(2)), x)), Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**(m + S(2))*s**(-m)*Int(1/(r**S(2) - s**S(2)*x**S(2)), x)/(a*n))))
    rubi.add(rule99)

    pattern100 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a, b: PositiveQ(a/b) | (AtomQ(SplitProduct(SumBaseQ, a)) & AtomQ(SplitProduct(SumBaseQ, b)) & PosQ(a/b))))
    rule100 = ReplacementRule(pattern100, lambda b, a, x : With(List(Set(r, Numerator(Rt(a/b, S(2)))), Set(s, Denominator(Rt(a/b, S(2))))), -Int((r - s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*s) + Int((r + s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*s)))
    rubi.add(rule100)

    pattern101 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a, b: ~(PositiveQ(a/b))))
    rule101 = ReplacementRule(pattern101, lambda b, a, x : With(List(Set(r, Numerator(Rt(-a/b, S(2)))), Set(s, Denominator(Rt(-a/b, S(2))))), -s*Int(1/(r - s*x**S(2)), x)/(S(2)*b) + s*Int(1/(r + s*x**S(2)), x)/(S(2)*b)))
    rubi.add(rule101)

    pattern102 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n/S(4))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n, m: Less(m, n + S(-1))), CustomConstraint(lambda a, b: PositiveQ(a/b)))
    rule102 = ReplacementRule(pattern102, lambda x, a, b, n, m : With(List(Set(r, Numerator(Rt(a/b, S(4)))), Set(s, Denominator(Rt(a/b, S(4))))), s**S(3)*Int(x**(m - n/S(4))/(r**S(2) - r*s*x**(n/S(4))*Sqrt(S(2)) + s**S(2)*x**(n/S(2))), x)/(S(2)*b*r*Sqrt(S(2))) - s**S(3)*Int(x**(m - n/S(4))/(r**S(2) + r*s*x**(n/S(4))*Sqrt(S(2)) + s**S(2)*x**(n/S(2))), x)/(S(2)*b*r*Sqrt(S(2)))))
    rubi.add(rule102)

    pattern103 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n/S(4))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n, m: Less(m, n/S(2))), CustomConstraint(lambda a, b: ~(PositiveQ(a/b))))
    rule103 = ReplacementRule(pattern103, lambda x, a, b, n, m : With(List(Set(r, Numerator(Rt(-a/b, S(2)))), Set(s, Denominator(Rt(-a/b, S(2))))), r*Int(x**m/(r - s*x**(n/S(2))), x)/(S(2)*a) + r*Int(x**m/(r + s*x**(n/S(2))), x)/(S(2)*a)))
    rubi.add(rule103)

    pattern104 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n/S(4))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n, m: Inequality(n/S(2), LessEqual, m, Less, n)), CustomConstraint(lambda a, b: ~(PositiveQ(a/b))))
    rule104 = ReplacementRule(pattern104, lambda x, a, b, n, m : With(List(Set(r, Numerator(Rt(-a/b, S(2)))), Set(s, Denominator(Rt(-a/b, S(2))))), -s*Int(x**(m - n/S(2))/(r - s*x**(n/S(2))), x)/(S(2)*b) + s*Int(x**(m - n/S(2))/(r + s*x**(n/S(2))), x)/(S(2)*b)))
    rubi.add(rule104)

    pattern105 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)), CustomConstraint(lambda n, m: Greater(m, S(2)*n + S(-1))))
    rule105 = ReplacementRule(pattern105, lambda x, a, b, n, m : Int(PolynomialDivide(x**m, a + b*x**n, x), x))
    rubi.add(rule105)

    pattern106 = Pattern(Integral(x_/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: PosQ(a)))
    rule106 = ReplacementRule(pattern106, lambda b, a, x : With(List(Set(r, Numer(Rt(b/a, S(3)))), Set(s, Denom(Rt(b/a, S(3))))), s*Int(1/Sqrt(a + b*x**S(3)), x)*Sqrt(S(2))/(r*Sqrt(Sqrt(S(3)) + S(2))) + Int((r*x + s*(-Sqrt(S(3)) + S(1)))/Sqrt(a + b*x**S(3)), x)/r))
    rubi.add(rule106)

    pattern107 = Pattern(Integral(x_/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: NegQ(a)))
    rule107 = ReplacementRule(pattern107, lambda b, a, x : With(List(Set(r, Numer(Rt(b/a, S(3)))), Set(s, Denom(Rt(b/a, S(3))))), -s*Int(1/Sqrt(a + b*x**S(3)), x)*Sqrt(S(2))/(r*Sqrt(-Sqrt(S(3)) + S(2))) + Int((r*x + s*(Sqrt(S(3)) + S(1)))/Sqrt(a + b*x**S(3)), x)/r))
    rubi.add(rule107)

    pattern108 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule108 = ReplacementRule(pattern108, lambda b, a, x : With(List(Set(q, Rt(b/a, S(2)))), -Int((-q*x**S(2) + S(1))/Sqrt(a + b*x**S(4)), x)/q + Int(1/Sqrt(a + b*x**S(4)), x)/q))
    rubi.add(rule108)

    pattern109 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda a: NegativeQ(a)), CustomConstraint(lambda b: PositiveQ(b)))
    rule109 = ReplacementRule(pattern109, lambda b, a, x : With(List(Set(q, Rt(-b/a, S(2)))), -Int((-q*x**S(2) + S(1))/Sqrt(a + b*x**S(4)), x)/q + Int(1/Sqrt(a + b*x**S(4)), x)/q))
    rubi.add(rule109)

    pattern110 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: NegQ(b/a)))
    rule110 = ReplacementRule(pattern110, lambda b, a, x : With(List(Set(q, Rt(-b/a, S(2)))), Int((q*x**S(2) + S(1))/Sqrt(a + b*x**S(4)), x)/q - Int(1/Sqrt(a + b*x**S(4)), x)/q))
    rubi.add(rule110)

    pattern111 = Pattern(Integral(x_**S(4)/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule111 = ReplacementRule(pattern111, lambda b, a, x : With(List(Set(r, Numer(Rt(b/a, S(3)))), Set(s, Denom(Rt(b/a, S(3))))), s**S(2)*(Sqrt(S(3)) + S(-1))*Int(1/Sqrt(a + b*x**S(6)), x)/(S(2)*r**S(2)) - Int((-S(2)*r**S(2)*x**S(4) + s**S(2)*(Sqrt(S(3)) + S(-1)))/Sqrt(a + b*x**S(6)), x)/(S(2)*r**S(2))))
    rubi.add(rule111)

    pattern112 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule112 = ReplacementRule(pattern112, lambda b, a, x : -Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/Sqrt(a + b*x**S(8)), x)/(S(2)*Rt(b/a, S(4))) + Int((x**S(2)*Rt(b/a, S(4)) + S(1))/Sqrt(a + b*x**S(8)), x)/(S(2)*Rt(b/a, S(4))))
    rubi.add(rule112)

    pattern113 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule113 = ReplacementRule(pattern113, lambda b, a, x : -a*Int(x**S(2)/(a + b*x**S(4))**(S(5)/4), x)/S(2) + x**S(3)/(S(2)*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule113)

    pattern114 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: NegQ(b/a)))
    rule114 = ReplacementRule(pattern114, lambda b, a, x : a*Int(S(1)/(x**S(2)*(a + b*x**S(4))**(S(1)/4)), x)/(S(2)*b) + (a + b*x**S(4))**(S(3)/4)/(S(2)*b*x))
    rubi.add(rule114)

    pattern115 = Pattern(Integral(S(1)/(x_**S(2)*(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule115 = ReplacementRule(pattern115, lambda b, a, x : -b*Int(x**S(2)/(a + b*x**S(4))**(S(5)/4), x) - S(1)/(x*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule115)

    pattern116 = Pattern(Integral(S(1)/(x_**S(2)*(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda b, a: NegQ(b/a)))
    rule116 = ReplacementRule(pattern116, lambda b, a, x : x*(a/(b*x**S(4)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(1)/4)), x)/(a + b*x**S(4))**(S(1)/4))
    rubi.add(rule116)

    pattern117 = Pattern(Integral(sqrt(c_*x_)/(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule117 = ReplacementRule(pattern117, lambda b, c, x, a : -a*Int(Sqrt(c*x)/(a + b*x**S(2))**(S(5)/4), x)/S(2) + x*Sqrt(c*x)/(a + b*x**S(2))**(S(1)/4))
    rubi.add(rule117)

    pattern118 = Pattern(Integral(sqrt(c_*x_)/(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, a: NegQ(b/a)))
    rule118 = ReplacementRule(pattern118, lambda b, c, x, a : a*c**S(2)*Int(S(1)/((c*x)**(S(3)/2)*(a + b*x**S(2))**(S(1)/4)), x)/(S(2)*b) + c*(a + b*x**S(2))**(S(3)/4)/(b*Sqrt(c*x)))
    rubi.add(rule118)

    pattern119 = Pattern(Integral(S(1)/((x_*WC('c', S(1)))**(S(3)/2)*(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule119 = ReplacementRule(pattern119, lambda b, c, x, a : -b*Int(Sqrt(c*x)/(a + b*x**S(2))**(S(5)/4), x)/c**S(2) - S(2)/(c*(a + b*x**S(2))**(S(1)/4)*Sqrt(c*x)))
    rubi.add(rule119)

    pattern120 = Pattern(Integral(S(1)/((x_*WC('c', S(1)))**(S(3)/2)*(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, a: NegQ(b/a)))
    rule120 = ReplacementRule(pattern120, lambda b, c, x, a : (a/(b*x**S(2)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(2)*(a/(b*x**S(2)) + S(1))**(S(1)/4)), x)*Sqrt(c*x)/(c**S(2)*(a + b*x**S(2))**(S(1)/4)))
    rubi.add(rule120)

    pattern121 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n, m: Greater(m, n + S(-1))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n*p + S(1))), CustomConstraint(lambda p, n, x, a, b, c, m: IntBinomialQ(a, b, c, n, m, p, x)))
    rule121 = ReplacementRule(pattern121, lambda p, x, n, a, b, c, m : -a*c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**p, x)/(b*(m + n*p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*p + S(1))))
    rubi.add(rule121)

    pattern122 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, m: SumSimplerQ(m, -n)), CustomConstraint(lambda p, n, m: NonzeroQ(m + n*p + S(1))), CustomConstraint(lambda p, n, m: NegativeIntegerQ(Simplify(p + (m + S(1))/n))))
    rule122 = ReplacementRule(pattern122, lambda p, x, n, a, b, c, m : -a*c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**p, x)/(b*(m + n*p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*p + S(1))))
    rubi.add(rule122)

    pattern123 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n, m: Greater(m, S(2)*n + S(-1))), CustomConstraint(lambda p, n, m: NonzeroQ(m + S(2)*n*p + S(1))), CustomConstraint(lambda p, b2, n, x, b1, a1, c, m, a2: IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x)))
    rule123 = ReplacementRule(pattern123, lambda p, b2, x, n, b1, a1, c, m, a2 : -a1*a2*c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(b1*b2*(m + S(2)*n*p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(b1*b2*(m + S(2)*n*p + S(1))))
    rubi.add(rule123)

    pattern124 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda n, m: SumSimplerQ(m, -S(2)*n)), CustomConstraint(lambda p, n, m: NonzeroQ(m + S(2)*n*p + S(1))), CustomConstraint(lambda p, n, m: NegativeIntegerQ(Simplify(p + (m + S(1))/(S(2)*n)))))
    rule124 = ReplacementRule(pattern124, lambda p, b2, x, n, b1, a1, c, m, a2 : -a1*a2*c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(b1*b2*(m + S(2)*n*p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(b1*b2*(m + S(2)*n*p + S(1))))
    rubi.add(rule124)

    pattern125 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p, n, x, a, b, c, m: IntBinomialQ(a, b, c, n, m, p, x)))
    rule125 = ReplacementRule(pattern125, lambda p, x, n, a, b, c, m : -b*c**(-n)*(m + n*(p + S(1)) + S(1))*Int((c*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))))
    rubi.add(rule125)

    pattern126 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, m: SumSimplerQ(m, n)), CustomConstraint(lambda p, n, m: NegativeIntegerQ(Simplify(p + (m + S(1))/n))))
    rule126 = ReplacementRule(pattern126, lambda p, x, n, a, b, c, m : -b*c**(-n)*(m + n*(p + S(1)) + S(1))*Int((c*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))))
    rubi.add(rule126)

    pattern127 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p, b2, n, x, b1, a1, c, m, a2: IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x)))
    rule127 = ReplacementRule(pattern127, lambda p, b2, x, n, b1, a1, c, m, a2 : -b1*b2*c**(-S(2)*n)*(m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))))
    rubi.add(rule127)

    pattern128 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda n, m: SumSimplerQ(m, S(2)*n)), CustomConstraint(lambda p, n, m: NegativeIntegerQ(Simplify(p + (m + S(1))/(S(2)*n)))))
    rule128 = ReplacementRule(pattern128, lambda p, b2, x, n, b1, a1, c, m, a2 : -b1*b2*c**(-S(2)*n)*(m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))))
    rubi.add(rule128)

    pattern129 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)), CustomConstraint(lambda p, n, x, a, b, c, m: IntBinomialQ(a, b, c, n, m, p, x)))
    rule129 = ReplacementRule(pattern129, lambda p, x, n, a, b, c, m : With(List(Set(k, Denominator(m))), k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*c**(-n)*x**(k*n))**p, x), x, (c*x)**(1/k))/c))
    rubi.add(rule129)

    pattern130 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda m: FractionQ(m)), CustomConstraint(lambda p, b2, n, x, b1, a1, c, m, a2: IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x)))
    rule130 = ReplacementRule(pattern130, lambda p, b2, x, n, b1, a1, c, m, a2 : With(List(Set(k, Denominator(m))), k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a1 + b1*c**(-n)*x**(k*n))**p*(a2 + b2*c**(-n)*x**(k*n))**p, x), x, (c*x)**(1/k))/c))
    rubi.add(rule130)

    pattern131 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(S(-1), p, S(0))), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda n, p, m: IntegersQ(m, p + (m + S(1))/n)))
    rule131 = ReplacementRule(pattern131, lambda p, x, a, b, n, m : a**(p + (m + S(1))/n)*Subst(Int(x**m*(-b*x**n + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule131)

    pattern132 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(S(-1), p, S(0))), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda n, p, m: IntegersQ(m, p + (m + S(1))/(S(2)*n))))
    rule132 = ReplacementRule(pattern132, lambda p, b2, x, b1, a1, n, m, a2 : (a1*a2)**(p + (m + S(1))/(S(2)*n))*Subst(Int(x**m*(-b1*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x*(a1 + b1*x**n)**(-S(1)/(S(2)*n))*(a2 + b2*x**n)**(-S(1)/(S(2)*n))))
    rubi.add(rule132)

    pattern133 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(S(-1), p, S(0))), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda n, p, m: Less(Denominator(p + (m + S(1))/n), Denominator(p))))
    rule133 = ReplacementRule(pattern133, lambda p, x, a, b, n, m : (a/(a + b*x**n))**(p + (m + S(1))/n)*(a + b*x**n)**(p + (m + S(1))/n)*Subst(Int(x**m*(-b*x**n + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule133)

    pattern134 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: PositiveIntegerQ(S(2)*n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(S(-1), p, S(0))), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda n, p, m: Less(Denominator(p + (m + S(1))/(S(2)*n)), Denominator(p))))
    rule134 = ReplacementRule(pattern134, lambda p, b2, x, b1, a1, n, m, a2 : (a1/(a1 + b1*x**n))**(p + (m + S(1))/(S(2)*n))*(a2/(a2 + b2*x**n))**(p + (m + S(1))/(S(2)*n))*(a1 + b1*x**n)**(p + (m + S(1))/(S(2)*n))*(a2 + b2*x**n)**(p + (m + S(1))/(S(2)*n))*Subst(Int(x**m*(-b1*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x*(a1 + b1*x**n)**(-S(1)/(S(2)*n))*(a2 + b2*x**n)**(-S(1)/(S(2)*n))))
    rubi.add(rule134)

    pattern135 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule135 = ReplacementRule(pattern135, lambda p, x, a, b, n, m : -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p, x), x, 1/x))
    rubi.add(rule135)

    pattern136 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: NegativeIntegerQ(S(2)*n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule136 = ReplacementRule(pattern136, lambda p, b2, x, b1, a1, n, m, a2 : -Subst(Int(x**(-m + S(-2))*(a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p, x), x, 1/x))
    rubi.add(rule136)

    pattern137 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)))
    rule137 = ReplacementRule(pattern137, lambda p, x, n, a, b, c, m : With(List(Set(k, Denominator(m))), -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*c**(-n)*x**(-k*n))**p, x), x, (c*x)**(-S(1)/k))/c))
    rubi.add(rule137)

    pattern138 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: NegativeIntegerQ(S(2)*n)), CustomConstraint(lambda m: FractionQ(m)))
    rule138 = ReplacementRule(pattern138, lambda p, b2, x, n, b1, a1, c, m, a2 : With(List(Set(k, Denominator(m))), -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a1 + b1*c**(-n)*x**(-k*n))**p*(a2 + b2*c**(-n)*x**(-k*n))**p, x), x, (c*x)**(-S(1)/k))/c))
    rubi.add(rule138)

    pattern139 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: ~(RationalQ(m))))
    rule139 = ReplacementRule(pattern139, lambda p, x, n, a, b, c, m : -(c*x)**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p, x), x, 1/x))
    rubi.add(rule139)

    pattern140 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: NegativeIntegerQ(S(2)*n)), CustomConstraint(lambda m: ~(RationalQ(m))))
    rule140 = ReplacementRule(pattern140, lambda p, b2, x, n, b1, a1, c, m, a2 : -(c*x)**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p, x), x, 1/x))
    rubi.add(rule140)

    pattern141 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: FractionQ(n)))
    rule141 = ReplacementRule(pattern141, lambda p, x, a, b, n, m : With(List(Set(k, Denominator(n))), k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n))**p, x), x, x**(1/k))))
    rubi.add(rule141)

    pattern142 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: FractionQ(S(2)*n)))
    rule142 = ReplacementRule(pattern142, lambda p, b2, x, b1, a1, n, m, a2 : With(List(Set(k, Denominator(S(2)*n))), k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a1 + b1*x**(k*n))**p*(a2 + b2*x**(k*n))**p, x), x, x**(1/k))))
    rubi.add(rule142)

    pattern143 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: FractionQ(n)))
    rule143 = ReplacementRule(pattern143, lambda p, x, n, a, b, c, m : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule143)

    pattern144 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n: FractionQ(S(2)*n)))
    rule144 = ReplacementRule(pattern144, lambda p, b2, x, n, b1, a1, c, m, a2 : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule144)

    pattern145 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: IntegerQ(Simplify(n/(m + S(1))))), CustomConstraint(lambda n: ~(IntegerQ(n))))
    rule145 = ReplacementRule(pattern145, lambda p, x, a, b, n, m : Subst(Int((a + b*x**Simplify(n/(m + S(1))))**p, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule145)

    pattern146 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n, m: IntegerQ(Simplify(S(2)*n/(m + S(1))))), CustomConstraint(lambda n: ~(IntegerQ(S(2)*n))))
    rule146 = ReplacementRule(pattern146, lambda p, b2, x, b1, a1, n, m, a2 : Subst(Int((a1 + b1*x**Simplify(n/(m + S(1))))**p*(a2 + b2*x**Simplify(n/(m + S(1))))**p, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule146)

    pattern147 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: IntegerQ(Simplify(n/(m + S(1))))), CustomConstraint(lambda n: ~(IntegerQ(n))))
    rule147 = ReplacementRule(pattern147, lambda p, x, n, a, b, c, m : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule147)

    pattern148 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n, m: IntegerQ(Simplify(S(2)*n/(m + S(1))))), CustomConstraint(lambda n: ~(IntegerQ(S(2)*n))))
    rule148 = ReplacementRule(pattern148, lambda p, b2, x, n, b1, a1, c, m, a2 : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule148)

    pattern149 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, n, m: ZeroQ(p + (m + S(1))/n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule149 = ReplacementRule(pattern149, lambda p, x, a, b, n, m : -b*n*p*Int(x**(m + n)*(a + b*x**n)**(p + S(-1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*x**n)**p/(m + S(1)))
    rubi.add(rule149)

    pattern150 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p, n, m: ZeroQ(p + (m + S(1))/(S(2)*n))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule150 = ReplacementRule(pattern150, lambda p, b2, x, b1, a1, n, m, a2 : -S(2)*b1*b2*n*p*Int(x**(m + n)*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(1)) + x**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(m + S(1)))
    rubi.add(rule150)

    pattern151 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, n, m: ZeroQ(p + (m + S(1))/n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule151 = ReplacementRule(pattern151, lambda p, x, n, a, b, c, m : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule151)

    pattern152 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p, n, m: ZeroQ(p + (m + S(1))/(S(2)*n))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule152 = ReplacementRule(pattern152, lambda p, b2, x, n, b1, a1, c, m, a2 : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule152)

    pattern153 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, p, m: IntegerQ(p + Simplify((m + S(1))/n))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n*p + S(1))))
    rule153 = ReplacementRule(pattern153, lambda p, x, n, a, b, c, m : a*n*p*Int((c*x)**m*(a + b*x**n)**(p + S(-1)), x)/(m + n*p + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + n*p + S(1))))
    rubi.add(rule153)

    pattern154 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n, p, m: IntegerQ(p + Simplify((m + S(1))/(S(2)*n)))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n, m: NonzeroQ(m + S(2)*n*p + S(1))))
    rule154 = ReplacementRule(pattern154, lambda p, b2, x, n, b1, a1, c, m, a2 : S(2)*a1*a2*n*p*Int((c*x)**m*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(2)*n*p + S(1)) + (c*x)**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(c*(m + S(2)*n*p + S(1))))
    rubi.add(rule154)

    pattern155 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, p, m: IntegerQ(p + Simplify((m + S(1))/n))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(S(-1), p, S(0))))
    rule155 = ReplacementRule(pattern155, lambda p, x, a, b, n, m : With(List(Set(k, Denominator(p))), a**(p + Simplify((m + S(1))/n))*k*Subst(Int(x**(k*Simplify((m + S(1))/n) + S(-1))*(-b*x**k + S(1))**(-p - Simplify((m + S(1))/n) + S(-1)), x), x, x**(n/k)*(a + b*x**n)**(-S(1)/k))/n))
    rubi.add(rule155)

    pattern156 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n, p, m: IntegerQ(p + Simplify((m + S(1))/(S(2)*n)))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(S(-1), p, S(0))))
    rule156 = ReplacementRule(pattern156, lambda p, b2, x, b1, a1, n, m, a2 : With(List(Set(k, Denominator(p))), k*(a1*a2)**(p + Simplify((m + S(1))/(S(2)*n)))*Subst(Int(x**(k*Simplify((m + S(1))/(S(2)*n)) + S(-1))*(-b1*x**k + S(1))**(-p - Simplify((m + S(1))/(S(2)*n)) + S(-1))*(-b2*x**k + S(1))**(-p - Simplify((m + S(1))/(S(2)*n)) + S(-1)), x), x, x**(S(2)*n/k)*(a1 + b1*x**n)**(-S(1)/k)*(a2 + b2*x**n)**(-S(1)/k))/(S(2)*n)))
    rubi.add(rule156)

    pattern157 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, p, m: IntegerQ(p + Simplify((m + S(1))/n))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(S(-1), p, S(0))))
    rule157 = ReplacementRule(pattern157, lambda p, x, n, a, b, c, m : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule157)

    pattern158 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n, p, m: IntegerQ(p + Simplify((m + S(1))/(S(2)*n)))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(S(-1), p, S(0))))
    rule158 = ReplacementRule(pattern158, lambda p, b2, x, n, b1, a1, c, m, a2 : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule158)

    pattern159 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, p, m: IntegerQ(p + Simplify((m + S(1))/n))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule159 = ReplacementRule(pattern159, lambda p, x, n, a, b, c, m : (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))))
    rubi.add(rule159)

    pattern160 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n, p, m: IntegerQ(p + Simplify((m + S(1))/n))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule160 = ReplacementRule(pattern160, lambda p, b2, x, n, b1, a1, c, m, a2 : (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))))
    rubi.add(rule160)

    pattern161 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: FractionQ(Simplify((m + S(1))/n))), CustomConstraint(lambda n, m: SumSimplerQ(m, -n)))
    rule161 = ReplacementRule(pattern161, lambda x, a, b, n, m : With(List(Set(mn, Simplify(m - n))), -a*Int(x**mn/(a + b*x**n), x)/b + x**(mn + S(1))/(b*(mn + S(1)))))
    rubi.add(rule161)

    pattern162 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: FractionQ(Simplify((m + S(1))/n))), CustomConstraint(lambda n, m: SumSimplerQ(m, n)))
    rule162 = ReplacementRule(pattern162, lambda x, a, b, n, m : -b*Int(x**Simplify(m + n)/(a + b*x**n), x)/a + x**(m + S(1))/(a*(m + S(1))))
    rubi.add(rule162)

    pattern163 = Pattern(Integral((c_*x_)**m_/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: FractionQ(Simplify((m + S(1))/n))), CustomConstraint(lambda n, m: SumSimplerQ(m, n) | SumSimplerQ(m, -n)))
    rule163 = ReplacementRule(pattern163, lambda x, n, a, b, c, m : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m/(a + b*x**n), x))
    rubi.add(rule163)

    pattern164 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: ~(PositiveIntegerQ(p))), CustomConstraint(lambda a, p: NegativeIntegerQ(p) | PositiveQ(a)))
    rule164 = ReplacementRule(pattern164, lambda p, x, n, a, b, c, m : a**p*(c*x)**(m + S(1))*Hypergeometric2F1(-p, (m + S(1))/n, S(1) + (m + S(1))/n, -b*x**n/a)/(c*(m + S(1))))
    rubi.add(rule164)

    pattern165 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: ~(PositiveIntegerQ(p))), CustomConstraint(lambda a, p: ~(NegativeIntegerQ(p) | PositiveQ(a))))
    rule165 = ReplacementRule(pattern165, lambda p, x, n, a, b, c, m : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((c*x)**m*(S(1) + b*x**n/a)**p, x))
    rubi.add(rule165)

    pattern166 = Pattern(Integral(x_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda v, x: NonzeroQ(v - x)))
    rule166 = ReplacementRule(pattern166, lambda n, x, a, v, b, p, m : Coefficient(v, x, S(1))**(-m + S(-1))*Subst(Int(SimplifyIntegrand((a + b*x**n)**p*(x - Coefficient(v, x, S(0)))**m, x), x), x, v))
    rubi.add(rule166)

    pattern167 = Pattern(Integral(u_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, v, u: LinearPairQ(u, v, x)))
    rule167 = ReplacementRule(pattern167, lambda n, x, a, v, u, b, p, m : u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p, x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule167)

    pattern168 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule168 = ReplacementRule(pattern168, lambda p, b2, x, n, b1, a1, c, m, a2 : (a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p))*Int((c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule168)

    pattern169 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)))
    rule169 = ReplacementRule(pattern169, lambda q, x, n, a, d, b, c, p : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule169)

    pattern170 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, q: IntegersQ(p, q)), CustomConstraint(lambda n: NegQ(n)))
    rule170 = ReplacementRule(pattern170, lambda q, x, n, a, d, b, c, p : Int(x**(n*(p + q))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q, x))
    rubi.add(rule170)

    pattern171 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule171 = ReplacementRule(pattern171, lambda q, x, n, a, d, b, c, p : -Subst(Int((a + b*x**(-n))**p*(c + d*x**(-n))**q/x**S(2), x), x, 1/x))
    rubi.add(rule171)

    pattern172 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: FractionQ(n)))
    rule172 = ReplacementRule(pattern172, lambda q, x, n, a, d, b, c, p : With(List(Set(g, Denominator(n))), g*Subst(Int(x**(g + S(-1))*(a + b*x**(g*n))**p*(c + d*x**(g*n))**q, x), x, x**(1/g))))
    rubi.add(rule172)

    pattern173 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_/(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n: ZeroQ(n*p + S(1))), CustomConstraint(lambda n: IntegerQ(n)))
    rule173 = ReplacementRule(pattern173, lambda p, x, a, d, b, c, n : Subst(Int(1/(c - x**n*(-a*d + b*c)), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule173)

    pattern174 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n, q: ZeroQ(n*(p + q + S(1)) + S(1))), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule174 = ReplacementRule(pattern174, lambda p, x, q, a, d, b, c, n : -c*q*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1)), x)/(a*(p + S(1))) - x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*n*(p + S(1))))
    rubi.add(rule174)

    pattern175 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n, q: ZeroQ(n*(p + q + S(1)) + S(1))), CustomConstraint(lambda p: NegativeIntegerQ(p)))
    rule175 = ReplacementRule(pattern175, lambda p, x, q, a, d, b, c, n : a**p*c**(-p + S(-1))*x*(c + d*x**n)**(-S(1)/n)*Hypergeometric2F1(1/n, -p, S(1) + 1/n, x**n*(a*d - b*c)/(a*(c + d*x**n))))
    rubi.add(rule175)

    pattern176 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n, q: ZeroQ(n*(p + q + S(1)) + S(1))))
    rule176 = ReplacementRule(pattern176, lambda p, x, q, a, d, b, c, n : x*(c*(a + b*x**n)/(a*(c + d*x**n)))**(-p)*(a + b*x**n)**p*(c + d*x**n)**(-p - S(1)/n)*Hypergeometric2F1(1/n, -p, S(1) + 1/n, x**n*(a*d - b*c)/(a*(c + d*x**n)))/c)
    rubi.add(rule176)

    pattern177 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n, q: ZeroQ(n*(p + q + S(2)) + S(1))), CustomConstraint(lambda q, a, d, b, c, p: ZeroQ(a*d*(p + S(1)) + b*c*(q + S(1)))))
    rule177 = ReplacementRule(pattern177, lambda p, x, q, a, d, b, c, n : x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c))
    rubi.add(rule177)

    pattern178 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n, q: ZeroQ(n*(p + q + S(2)) + S(1))), CustomConstraint(lambda p, q: (RationalQ(p) & Less(p, S(-1))) | ~(RationalQ(q) & Less(q, S(-1)))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule178 = ReplacementRule(pattern178, lambda p, x, q, a, d, b, c, n : -b*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*n*(p + S(1))*(-a*d + b*c)) + (b*c + n*(p + S(1))*(-a*d + b*c))*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q, x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule178)

    pattern179 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n, a, d, b, c: ZeroQ(a*d - b*c*(n*(p + S(1)) + S(1)))))
    rule179 = ReplacementRule(pattern179, lambda x, n, a, d, b, c, p : c*x*(a + b*x**n)**(p + S(1))/a)
    rubi.add(rule179)

    pattern180 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n: NegativeIntegerQ(p + 1/n) | (RationalQ(p) & Less(p, S(-1)))))
    rule180 = ReplacementRule(pattern180, lambda p, x, a, d, b, c, n : x*(a + b*x**n)**(p + S(1))*(a*d - b*c)/(a*b*n*(p + S(1))) - (a*d - b*c*(n*(p + S(1)) + S(1)))*Int((a + b*x**n)**(p + S(1)), x)/(a*b*n*(p + S(1))))
    rubi.add(rule180)

    pattern181 = Pattern(Integral((c_ + x_**n_*WC('d', S(1)))/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(0))))
    rule181 = ReplacementRule(pattern181, lambda x, a, d, b, c, n : c*x/a - (-a*d + b*c)*Int(1/(a*x**(-n) + b), x)/a)
    rubi.add(rule181)

    pattern182 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n: NonzeroQ(n*(p + S(1)) + S(1))))
    rule182 = ReplacementRule(pattern182, lambda p, x, a, d, b, c, n : d*x*(a + b*x**n)**(p + S(1))/(b*(n*(p + S(1)) + S(1))) - (a*d - b*c*(n*(p + S(1)) + S(1)))*Int((a + b*x**n)**p, x)/(b*(n*(p + S(1)) + S(1))))
    rubi.add(rule182)

    pattern183 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)), CustomConstraint(lambda q: NegativeIntegerQ(q)), CustomConstraint(lambda p, q: GreaterEqual(p, -q)))
    rule183 = ReplacementRule(pattern183, lambda p, x, q, a, d, b, c, n : Int(PolynomialDivide((a + b*x**n)**p, (c + d*x**n)**(-q), x), x))
    rubi.add(rule183)

    pattern184 = Pattern(Integral(S(1)/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)))
    rule184 = ReplacementRule(pattern184, lambda x, a, d, b, c, n : b*Int(1/(a + b*x**n), x)/(-a*d + b*c) - d*Int(1/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule184)

    pattern185 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/3)*(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda b, c, d, a: ZeroQ(S(3)*a*d + b*c)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule185 = ReplacementRule(pattern185, lambda x, a, d, b, c : Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(-x*Rt(b/a, S(2)) + Sqrt(S(3)))), x)*Sqrt(S(3))/(S(2)*c) + Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(x*Rt(b/a, S(2)) + Sqrt(S(3)))), x)*Sqrt(S(3))/(S(2)*c))
    rubi.add(rule185)

    pattern186 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/3)*(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda b, c, d, a: ZeroQ(S(3)*a*d + b*c)), CustomConstraint(lambda b, a: NegQ(b/a)))
    rule186 = ReplacementRule(pattern186, lambda x, a, d, b, c : Int((-x*Rt(-b/a, S(2)) + S(3))/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/S(6) + Int((x*Rt(-b/a, S(2)) + S(3))/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/S(6))
    rubi.add(rule186)

    pattern187 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(2)/3)/(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda b, c, d, a: ZeroQ(S(3)*a*d + b*c)))
    rule187 = ReplacementRule(pattern187, lambda x, a, d, b, c : b*Int((a + b*x**S(2))**(S(-1)/3), x)/d - (-a*d + b*c)*Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/d)
    rubi.add(rule187)

    pattern188 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)*(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)))
    rule188 = ReplacementRule(pattern188, lambda x, a, d, b, c : Sqrt(-b*x**S(2)/a)*Subst(Int(S(1)/((a + b*x)**(S(1)/4)*(c + d*x)*Sqrt(-b*x/a)), x), x, x**S(2))/(S(2)*x))
    rubi.add(rule188)

    pattern189 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)))
    rule189 = ReplacementRule(pattern189, lambda x, a, d, b, c : Sqrt(-b*x**S(2)/a)*Subst(Int(S(1)/((a + b*x)**(S(3)/4)*(c + d*x)*Sqrt(-b*x/a)), x), x, x**S(2))/(S(2)*x))
    rubi.add(rule189)

    pattern190 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**WC('p', S(1))/(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p: Equal(p, S(1)/2) | Equal(Denominator(p), S(4))))
    rule190 = ReplacementRule(pattern190, lambda x, a, d, b, c, p : b*Int((a + b*x**S(2))**(p + S(-1)), x)/d - (-a*d + b*c)*Int((a + b*x**S(2))**(p + S(-1))/(c + d*x**S(2)), x)/d)
    rubi.add(rule190)

    pattern191 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**p_/(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Equal(Denominator(p), S(4))), CustomConstraint(lambda p: Equal(p, S(-7)/4) | Equal(p, S(-5)/4)))
    rule191 = ReplacementRule(pattern191, lambda x, a, d, b, c, p : b*Int((a + b*x**S(2))**p, x)/(-a*d + b*c) - d*Int((a + b*x**S(2))**(p + S(1))/(c + d*x**S(2)), x)/(-a*d + b*c))
    rubi.add(rule191)

    pattern192 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: ZeroQ(a*d + b*c)), CustomConstraint(lambda a, b: PosQ(a*b)))
    rule192 = ReplacementRule(pattern192, lambda x, a, d, b, c : a*Subst(Int(1/(-S(4)*a*b*x**S(4) + S(1)), x), x, x/Sqrt(a + b*x**S(4)))/c)
    rubi.add(rule192)

    pattern193 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: ZeroQ(a*d + b*c)), CustomConstraint(lambda a, b: NegQ(a*b)))
    rule193 = ReplacementRule(pattern193, lambda x, a, d, b, c : With(List(Set(q, Rt(-a*b, S(4)))), a*ArcTan(q*x*(a + q**S(2)*x**S(2))/(a*Sqrt(a + b*x**S(4))))/(S(2)*c*q) + a*ArcTanh(q*x*(a - q**S(2)*x**S(2))/(a*Sqrt(a + b*x**S(4))))/(S(2)*c*q)))
    rubi.add(rule193)

    pattern194 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)))
    rule194 = ReplacementRule(pattern194, lambda x, a, d, b, c : b*Int(1/Sqrt(a + b*x**S(4)), x)/d - (-a*d + b*c)*Int(S(1)/((c + d*x**S(4))*Sqrt(a + b*x**S(4))), x)/d)
    rubi.add(rule194)

    pattern195 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)/(c_ + x_**S(4)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)))
    rule195 = ReplacementRule(pattern195, lambda x, a, d, b, c : Sqrt(a/(a + b*x**S(4)))*Sqrt(a + b*x**S(4))*Subst(Int(S(1)/((c - x**S(4)*(-a*d + b*c))*Sqrt(-b*x**S(4) + S(1))), x), x, x/(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule195)

    pattern196 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**p_/(c_ + x_**S(4)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Equal(p, S(3)/4) | Equal(p, S(5)/4)))
    rule196 = ReplacementRule(pattern196, lambda x, a, d, b, c, p : b*Int((a + b*x**S(4))**(p + S(-1)), x)/d - (-a*d + b*c)*Int((a + b*x**S(4))**(p + S(-1))/(c + d*x**S(4)), x)/d)
    rubi.add(rule196)

    pattern197 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('b', S(1)))*(c_ + x_**S(4)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)))
    rule197 = ReplacementRule(pattern197, lambda x, a, d, b, c : Int(S(1)/((-x**S(2)*Rt(-d/c, S(2)) + S(1))*Sqrt(a + b*x**S(4))), x)/(S(2)*c) + Int(S(1)/((x**S(2)*Rt(-d/c, S(2)) + S(1))*Sqrt(a + b*x**S(4))), x)/(S(2)*c))
    rubi.add(rule197)

    pattern198 = Pattern(Integral(S(1)/((a_ + x_**S(4)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(4)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)))
    rule198 = ReplacementRule(pattern198, lambda x, a, d, b, c : b*Int((a + b*x**S(4))**(S(-3)/4), x)/(-a*d + b*c) - d*Int((a + b*x**S(4))**(S(1)/4)/(c + d*x**S(4)), x)/(-a*d + b*c))
    rubi.add(rule198)

    pattern199 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, a: PosQ(b/a)), CustomConstraint(lambda c, d: PosQ(d/c)))
    rule199 = ReplacementRule(pattern199, lambda x, a, d, b, c : EllipticE(ArcTan(x*Rt(d/c, S(2))), S(1) - b*c/(a*d))*Sqrt(a + b*x**S(2))/(c*Rt(d/c, S(2))*Sqrt(c*(a + b*x**S(2))/(a*(c + d*x**S(2))))*Sqrt(c + d*x**S(2))))
    rubi.add(rule199)

    pattern200 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda q: Less(S(0), q, S(1))), CustomConstraint(lambda p, n, q, x, a, d, b, c: IntBinomialQ(a, b, c, d, n, p, q, x)))
    rule200 = ReplacementRule(pattern200, lambda p, x, q, a, d, b, c, n : -x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*n*(p + S(1))) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(n*(p + S(1)) + S(1)) + d*x**n*(n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1))))
    rubi.add(rule200)

    pattern201 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda q: Greater(q, S(1))), CustomConstraint(lambda p, n, q, x, a, d, b, c: IntBinomialQ(a, b, c, d, n, p, q, x)))
    rule201 = ReplacementRule(pattern201, lambda p, x, q, a, d, b, c, n : x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(a*d - b*c)/(a*b*n*(p + S(1))) - Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(a*d - b*c*(n*(p + S(1)) + S(1))) + d*x**n*(a*d*(n*(q + S(-1)) + S(1)) - b*c*(n*(p + q) + S(1))), x), x)/(a*b*n*(p + S(1))))
    rubi.add(rule201)

    pattern202 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, q: ~(IntegerQ(q) & Less(q, S(-1)) & ~(IntegerQ(p)))), CustomConstraint(lambda p, n, q, x, a, d, b, c: IntBinomialQ(a, b, c, d, n, p, q, x)))
    rule202 = ReplacementRule(pattern202, lambda p, x, q, a, d, b, c, n : -b*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*n*(p + S(1))*(-a*d + b*c)) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c + b*d*x**n*(n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule202)

    pattern203 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, q: IntegersQ(p, q)), CustomConstraint(lambda p, q: Greater(p + q, S(0))))
    rule203 = ReplacementRule(pattern203, lambda p, x, q, a, d, b, c, n : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule203)

    pattern204 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: Greater(q, S(1))), CustomConstraint(lambda p, n, q: NonzeroQ(n*(p + q) + S(1))), CustomConstraint(lambda p: ~(IntegerQ(p) & Greater(p, S(1)))), CustomConstraint(lambda p, n, q, x, a, d, b, c: IntBinomialQ(a, b, c, d, n, p, q, x)))
    rule204 = ReplacementRule(pattern204, lambda p, x, q, a, d, b, c, n : d*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*(n*(p + q) + S(1))) + Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(-a*d + b*c*(n*(p + q) + S(1))) + d*x**n*(-a*d*(n*(q + S(-1)) + S(1)) + b*c*(n*(p + S(2)*q + S(-1)) + S(1))), x), x)/(b*(n*(p + q) + S(1))))
    rubi.add(rule204)

    pattern205 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n, q, x, a, d, b, c: IntBinomialQ(a, b, c, d, n, p, q, x)))
    rule205 = ReplacementRule(pattern205, lambda p, x, q, a, d, b, c, n : n*Int((a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(n*(p + q) + S(1)) + x*(a + b*x**n)**p*(c + d*x**n)**q/(n*(p + q) + S(1)))
    rubi.add(rule205)

    pattern206 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: PosQ(d/c)), CustomConstraint(lambda b, a: PosQ(b/a)), CustomConstraint(lambda b, a, c, d: ~(SimplerSqrtQ(b/a, d/c))))
    rule206 = ReplacementRule(pattern206, lambda x, a, d, b, c : EllipticF(ArcTan(x*Rt(d/c, S(2))), S(1) - b*c/(a*d))*Sqrt(a + b*x**S(2))/(a*Rt(d/c, S(2))*Sqrt(c*(a + b*x**S(2))/(a*(c + d*x**S(2))))*Sqrt(c + d*x**S(2))))
    rubi.add(rule206)

    pattern207 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NegQ(d/c)), CustomConstraint(lambda c: PositiveQ(c)), CustomConstraint(lambda a: PositiveQ(a)), CustomConstraint(lambda b, a, c, d: ~(NegQ(b/a) & SimplerSqrtQ(-b/a, -d/c))))
    rule207 = ReplacementRule(pattern207, lambda x, a, d, b, c : EllipticF(ArcSin(x*Rt(-d/c, S(2))), b*c/(a*d))/(Rt(-d/c, S(2))*Sqrt(a)*Sqrt(c)))
    rubi.add(rule207)

    pattern208 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NegQ(d/c)), CustomConstraint(lambda c: PositiveQ(c)), CustomConstraint(lambda a, b, c, d: PositiveQ(a - b*c/d)))
    rule208 = ReplacementRule(pattern208, lambda x, a, d, b, c : -EllipticF(ArcCos(x*Rt(-d/c, S(2))), b*c/(-a*d + b*c))/(Rt(-d/c, S(2))*Sqrt(c)*Sqrt(a - b*c/d)))
    rubi.add(rule208)

    pattern209 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ~(PositiveQ(c))))
    rule209 = ReplacementRule(pattern209, lambda x, a, d, b, c : Int(S(1)/(Sqrt(S(1) + d*x**S(2)/c)*Sqrt(a + b*x**S(2))), x)*Sqrt(S(1) + d*x**S(2)/c)/Sqrt(c + d*x**S(2)))
    rubi.add(rule209)

    pattern210 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: PosQ(d/c)), CustomConstraint(lambda b, a: PosQ(b/a)))
    rule210 = ReplacementRule(pattern210, lambda x, a, d, b, c : a*Int(S(1)/(Sqrt(a + b*x**S(2))*Sqrt(c + d*x**S(2))), x) + b*Int(x**S(2)/(Sqrt(a + b*x**S(2))*Sqrt(c + d*x**S(2))), x))
    rubi.add(rule210)

    pattern211 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: PosQ(d/c)), CustomConstraint(lambda b, a: NegQ(b/a)))
    rule211 = ReplacementRule(pattern211, lambda x, a, d, b, c : b*Int(Sqrt(c + d*x**S(2))/Sqrt(a + b*x**S(2)), x)/d - (-a*d + b*c)*Int(S(1)/(Sqrt(a + b*x**S(2))*Sqrt(c + d*x**S(2))), x)/d)
    rubi.add(rule211)

    pattern212 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NegQ(d/c)), CustomConstraint(lambda c: PositiveQ(c)), CustomConstraint(lambda a: PositiveQ(a)))
    rule212 = ReplacementRule(pattern212, lambda x, a, d, b, c : EllipticE(ArcSin(x*Rt(-d/c, S(2))), b*c/(a*d))*Sqrt(a)/(Rt(-d/c, S(2))*Sqrt(c)))
    rubi.add(rule212)

    pattern213 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NegQ(d/c)), CustomConstraint(lambda c: PositiveQ(c)), CustomConstraint(lambda a, b, c, d: PositiveQ(a - b*c/d)))
    rule213 = ReplacementRule(pattern213, lambda x, a, d, b, c : -EllipticE(ArcCos(x*Rt(-d/c, S(2))), b*c/(-a*d + b*c))*Sqrt(a - b*c/d)/(Rt(-d/c, S(2))*Sqrt(c)))
    rubi.add(rule213)

    pattern214 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NegQ(d/c)), CustomConstraint(lambda c: PositiveQ(c)), CustomConstraint(lambda a: ~(PositiveQ(a))))
    rule214 = ReplacementRule(pattern214, lambda x, a, d, b, c : Int(Sqrt(S(1) + b*x**S(2)/a)/Sqrt(c + d*x**S(2)), x)*Sqrt(a + b*x**S(2))/Sqrt(S(1) + b*x**S(2)/a))
    rubi.add(rule214)

    pattern215 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NegQ(d/c)), CustomConstraint(lambda c: ~(PositiveQ(c))))
    rule215 = ReplacementRule(pattern215, lambda x, a, d, b, c : Int(Sqrt(a + b*x**S(2))/Sqrt(S(1) + d*x**S(2)/c), x)*Sqrt(S(1) + d*x**S(2)/c)/Sqrt(c + d*x**S(2)))
    rubi.add(rule215)

    pattern216 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule216 = ReplacementRule(pattern216, lambda p, x, q, a, d, b, c, n : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule216)

    pattern217 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: NonzeroQ(n + S(1))), CustomConstraint(lambda a: PositiveQ(a)), CustomConstraint(lambda c: PositiveQ(c)))
    rule217 = ReplacementRule(pattern217, lambda p, x, q, a, d, b, c, n : a**p*c**q*x*AppellF1(1/n, -p, -q, S(1) + 1/n, -b*x**n/a, -d*x**n/c))
    rubi.add(rule217)

    pattern218 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: NonzeroQ(n + S(1))), CustomConstraint(lambda a: ~(PositiveQ(a))))
    rule218 = ReplacementRule(pattern218, lambda p, x, q, a, d, b, c, n : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((S(1) + b*x**n/a)**p*(c + d*x**n)**q, x))
    rubi.add(rule218)

    pattern219 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)), CustomConstraint(lambda q: IntegerQ(q)), CustomConstraint(lambda p, n: PosQ(n) | ~(IntegerQ(p))))
    rule219 = ReplacementRule(pattern219, lambda p, q, x, mn, a, d, b, c, n : Int(x**(-n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule219)

    pattern220 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)), CustomConstraint(lambda q: ~(IntegerQ(q))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule220 = ReplacementRule(pattern220, lambda p, x, q, mn, a, d, b, c, n : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(-n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule220)

    pattern221 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda x, u: NonzeroQ(u - x)))
    rule221 = ReplacementRule(pattern221, lambda q, n, x, a, u, d, b, c, p : Subst(Int((a + b*x**n)**p*(c + d*x**n)**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule221)

    pattern222 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1)), x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda x, v, u: PseudoBinomialPairQ(u, v, x)))
    rule222 = ReplacementRule(pattern222, lambda x, q, v, u, p : Int(NormalizePseudoBinomial(u, x)**p*NormalizePseudoBinomial(v, x)**q, x))
    rubi.add(rule222)

    pattern223 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*x_**WC('m', S(1)), x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p, m: IntegersQ(p, m/p)), CustomConstraint(lambda x, v, u, p, m: PseudoBinomialPairQ(u*x**(m/p), v, x)))
    rule223 = ReplacementRule(pattern223, lambda q, x, v, u, p, m : Int(NormalizePseudoBinomial(v, x)**q*NormalizePseudoBinomial(u*x**(m/p), x)**p, x))
    rubi.add(rule223)

    pattern224 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda m, e: IntegerQ(m) | PositiveQ(e)), CustomConstraint(lambda n, m: IntegerQ(Simplify((m + S(1))/n))))
    rule224 = ReplacementRule(pattern224, lambda p, q, x, d, b, c, e, n, m : b**(-Simplify((m + S(1))/n) + S(1))*e**m*Subst(Int((b*x)**(p + Simplify((m + S(1))/n) + S(-1))*(c + d*x)**q, x), x, x**n)/n)
    rubi.add(rule224)

    pattern225 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda m, e: IntegerQ(m) | PositiveQ(e)), CustomConstraint(lambda n, m: ~(IntegerQ(Simplify((m + S(1))/n)))))
    rule225 = ReplacementRule(pattern225, lambda p, q, x, d, b, c, e, n, m : b**IntPart(p)*e**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p)*(c + d*x**n)**q, x))
    rubi.add(rule225)

    pattern226 = Pattern(Integral((e_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule226 = ReplacementRule(pattern226, lambda p, q, x, d, b, c, e, n, m : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule226)

    pattern227 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))))
    rule227 = ReplacementRule(pattern227, lambda q, x, n, a, d, b, c, p, m : Subst(Int((a + b*x)**p*(c + d*x)**q, x), x, x**n)/n)
    rubi.add(rule227)

    pattern228 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, q: IntegersQ(p, q)), CustomConstraint(lambda n: NegQ(n)))
    rule228 = ReplacementRule(pattern228, lambda q, x, n, a, d, b, c, p, m : Int(x**(m + n*(p + q))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q, x))
    rubi.add(rule228)

    pattern229 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, m: IntegerQ(Simplify((m + S(1))/n))))
    rule229 = ReplacementRule(pattern229, lambda q, x, n, a, d, b, c, p, m : Subst(Int(x**(Simplify((m + S(1))/n) + S(-1))*(a + b*x)**p*(c + d*x)**q, x), x, x**n)/n)
    rubi.add(rule229)

    pattern230 = Pattern(Integral((e_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, m: IntegerQ(Simplify((m + S(1))/n))))
    rule230 = ReplacementRule(pattern230, lambda q, x, n, a, d, b, c, e, p, m : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule230)

    pattern231 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)))
    rule231 = ReplacementRule(pattern231, lambda q, x, n, a, d, b, c, e, p, m : Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule231)

    pattern232 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n, a, d, b, c, m: ZeroQ(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule232 = ReplacementRule(pattern232, lambda x, n, a, d, b, c, e, p, m : c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))))
    rubi.add(rule232)

    pattern233 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda non2, n: ZeroQ(-n/S(2) + non2)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p, b2, n, b1, d, a1, c, m, a2: ZeroQ(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule233 = ReplacementRule(pattern233, lambda b2, x, n, non2, b1, d, c, e, a1, p, m, a2 : c*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(a1*a2*e*(m + S(1))))
    rubi.add(rule233)

    pattern234 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n, m: ZeroQ(m + n*(p + S(1)) + S(1))), CustomConstraint(lambda n, e: IntegerQ(n) | PositiveQ(e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n, m: (Greater(n, S(0)) & Less(m, S(-1))) | (Less(n, S(0)) & Greater(m + n, S(-1)))))
    rule234 = ReplacementRule(pattern234, lambda x, n, a, d, b, c, e, p, m : d*e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p, x) + c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))))
    rubi.add(rule234)

    pattern235 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n, m: ZeroQ(m + n*(p + S(1)) + S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule235 = ReplacementRule(pattern235, lambda x, n, a, d, b, c, e, p, m : d*Int((e*x)**m*(a + b*x**n)**(p + S(1)), x)/b + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*e*(m + S(1))))
    rubi.add(rule235)

    pattern236 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, e: IntegerQ(n) | PositiveQ(e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n, m: (Greater(n, S(0)) & Less(m, S(-1))) | (Less(n, S(0)) & Greater(m + n, S(-1)))), CustomConstraint(lambda p: ~(IntegerQ(p) & Less(p, S(-1)))))
    rule236 = ReplacementRule(pattern236, lambda x, n, a, d, b, c, e, p, m : c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))) + e**(-n)*(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))))
    rubi.add(rule236)

    pattern237 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda non2, n: ZeroQ(-n/S(2) + non2)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n, e: IntegerQ(n) | PositiveQ(e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n, m: (Greater(n, S(0)) & Less(m, S(-1))) | (Less(n, S(0)) & Greater(m + n, S(-1)))), CustomConstraint(lambda p: ~(IntegerQ(p) & Less(p, S(-1)))))
    rule237 = ReplacementRule(pattern237, lambda b2, x, n, non2, b1, d, c, e, a1, p, m, a2 : c*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(a1*a2*e*(m + S(1))) + e**(-n)*(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**(m + n)*(a1 + b1*x**(n/S(2)))**p*(a2 + b2*x**(n/S(2)))**p, x)/(a1*a2*(m + S(1))))
    rubi.add(rule237)

    pattern238 = Pattern(Integral(x_**m_*(a_ + x_**S(2)*WC('b', S(1)))**p_*(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m/S(2))), CustomConstraint(lambda p, m: IntegerQ(p) | Equal(m + S(2)*p + S(1), S(0))))
    rule238 = ReplacementRule(pattern238, lambda x, a, d, b, c, p, m : b**(-m/S(2) + S(-1))*x*(-a)**(m/S(2) + S(-1))*(a + b*x**S(2))**(p + S(1))*(-a*d + b*c)/(S(2)*(p + S(1))) + b**(-m/S(2) + S(-1))*Int((a + b*x**S(2))**(p + S(1))*ExpandToSum(S(2)*b*x**S(2)*(p + S(1))*Together((b**(m/S(2))*x**(m + S(-2))*(c + d*x**S(2)) - (-a)**(m/S(2) + S(-1))*(-a*d + b*c))/(a + b*x**S(2))) - (-a)**(m/S(2) + S(-1))*(-a*d + b*c), x), x)/(S(2)*(p + S(1))))
    rubi.add(rule238)

    pattern239 = Pattern(Integral(x_**m_*(a_ + x_**S(2)*WC('b', S(1)))**p_*(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2))), CustomConstraint(lambda p, m: IntegerQ(p) | Equal(m + S(2)*p + S(1), S(0))))
    rule239 = ReplacementRule(pattern239, lambda x, a, d, b, c, p, m : b**(-m/S(2) + S(-1))*x*(-a)**(m/S(2) + S(-1))*(a + b*x**S(2))**(p + S(1))*(-a*d + b*c)/(S(2)*(p + S(1))) + b**(-m/S(2) + S(-1))*Int(x**m*(a + b*x**S(2))**(p + S(1))*ExpandToSum(S(2)*b*(p + S(1))*Together((b**(m/S(2))*(c + d*x**S(2)) - x**(-m + S(2))*(-a)**(m/S(2) + S(-1))*(-a*d + b*c))/(a + b*x**S(2))) - x**(-m)*(-a)**(m/S(2) + S(-1))*(-a*d + b*c), x), x)/(S(2)*(p + S(1))))
    rubi.add(rule239)

    pattern240 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, p, m: IntegerQ(p) | ~(RationalQ(m)) | (PositiveIntegerQ(n) & NegativeIntegerQ(p + S(1)/2) & LessEqual(S(-1), m, -n*(p + S(1))))))
    rule240 = ReplacementRule(pattern240, lambda x, n, a, d, b, c, e, p, m : -(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*b*n*(p + S(1))) + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(a*d - b*c)/(a*b*e*n*(p + S(1))))
    rubi.add(rule240)

    pattern241 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda non2, n: ZeroQ(-n/S(2) + non2)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, p, m: IntegerQ(p) | ~(RationalQ(m)) | (PositiveIntegerQ(n) & NegativeIntegerQ(p + S(1)/2) & LessEqual(S(-1), m, -n*(p + S(1))))))
    rule241 = ReplacementRule(pattern241, lambda b2, x, n, non2, b1, d, c, e, a1, p, m, a2 : -(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1)), x)/(a1*a2*b1*b2*n*(p + S(1))) + (e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))*(a1*a2*d - b1*b2*c)/(a1*a2*b1*b2*e*n*(p + S(1))))
    rubi.add(rule241)

    pattern242 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n, m: NonzeroQ(m + n*(p + S(1)) + S(1))))
    rule242 = ReplacementRule(pattern242, lambda x, n, a, d, b, c, e, p, m : d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(b*e*(m + n*(p + S(1)) + S(1))) - (a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a + b*x**n)**p, x)/(b*(m + n*(p + S(1)) + S(1))))
    rubi.add(rule242)

    pattern243 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda non2, n: ZeroQ(-n/S(2) + non2)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p, n, m: NonzeroQ(m + n*(p + S(1)) + S(1))))
    rule243 = ReplacementRule(pattern243, lambda b2, x, n, non2, b1, d, c, e, a1, p, m, a2 : d*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(b1*b2*e*(m + n*(p + S(1)) + S(1))) - (a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a1 + b1*x**(n/S(2)))**p*(a2 + b2*x**(n/S(2)))**p, x)/(b1*b2*(m + n*(p + S(1)) + S(1))))
    rubi.add(rule243)

    pattern244 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_/(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m) | ~(RationalQ(m)) | PositiveIntegerQ(S(2)*m + S(2))))
    rule244 = ReplacementRule(pattern244, lambda p, x, a, d, b, c, e, n, m : Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p/(c + d*x**n), x), x))
    rubi.add(rule244)

    pattern245 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda n: Greater(n, S(0))))
    rule245 = ReplacementRule(pattern245, lambda p, x, a, d, b, c, e, n, m : c**S(2)*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))) - e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*Simp(-a*d**S(2)*x**n*(m + S(1)) + b*c**S(2)*n*(p + S(1)) + c*(m + S(1))*(-S(2)*a*d + b*c), x), x)/(a*(m + S(1))))
    rubi.add(rule245)

    pattern246 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule246 = ReplacementRule(pattern246, lambda p, x, a, d, b, c, e, n, m : Int((e*x)**m*(a + b*x**n)**(p + S(1))*Simp(a*b*d**S(2)*n*x**n*(p + S(1)) + b**S(2)*c**S(2)*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)**S(2), x), x)/(a*b**S(2)*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)**S(2)/(a*b**S(2)*e*n*(p + S(1))))
    rubi.add(rule246)

    pattern247 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, n, m: NonzeroQ(m + n*(p + S(2)) + S(1))))
    rule247 = ReplacementRule(pattern247, lambda p, x, a, d, b, c, e, n, m : d**S(2)*e**(-n + S(-1))*(e*x)**(m + n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*(p + S(2)) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*Simp(b*c**S(2)*(m + n*(p + S(2)) + S(1)) + d*x**n*(S(2)*b*c*n*(p + S(1)) + (-a*d + S(2)*b*c)*(m + n + S(1))), x), x)/(b*(m + n*(p + S(2)) + S(1))))
    rubi.add(rule247)

    pattern248 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule248 = ReplacementRule(pattern248, lambda p, x, q, a, d, b, c, n, m : With(List(Set(k, GCD(m + S(1), n))), Condition(Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p*(c + d*x**(n/k))**q, x), x, x**k)/k, Unequal(k, S(1)))))
    rubi.add(rule248)

    pattern249 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)), CustomConstraint(lambda p: IntegerQ(p)))
    rule249 = ReplacementRule(pattern249, lambda p, x, q, a, d, b, c, e, n, m : With(List(Set(k, Denominator(m))), k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*e**(-n)*x**(k*n))**p*(c + d*e**(-n)*x**(k*n))**q, x), x, (e*x)**(1/k))/e))
    rubi.add(rule249)

    pattern250 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, q, m: RationalQ(m, p, q)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda n, m: Greater(m - n + S(1), S(0))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule250 = ReplacementRule(pattern250, lambda p, x, q, a, d, b, c, e, n, m : -e**n*Int((e*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m - n + S(1)) + d*x**n*(m + n*(q + S(-1)) + S(1)), x), x)/(b*n*(p + S(1))) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*n*(p + S(1))))
    rubi.add(rule250)

    pattern251 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda q: Greater(q, S(1))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule251 = ReplacementRule(pattern251, lambda p, x, q, a, d, b, c, e, n, m : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)) + d*x**n*(b*c*n*(p + S(1)) + (-a*d + b*c)*(m + n*(q + S(-1)) + S(1))), x), x)/(a*b*n*(p + S(1))) + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(a*d - b*c)/(a*b*e*n*(p + S(1))))
    rubi.add(rule251)

    pattern252 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda q: Less(S(0), q, S(1))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule252 = ReplacementRule(pattern252, lambda p, x, q, a, d, b, c, e, n, m : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m + n*(p + S(1)) + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*n*(p + S(1))))
    rubi.add(rule252)

    pattern253 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, m: Greater(m - n + S(1), n)), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule253 = ReplacementRule(pattern253, lambda p, x, q, a, d, b, c, e, n, m : -a*e**(S(2)*n + S(-1))*(e*x)**(m - S(2)*n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*n*(p + S(1))*(-a*d + b*c)) + e**(S(2)*n)*Int((e*x)**(m - S(2)*n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(a*c*(m - S(2)*n + S(1)) + x**n*(a*d*(m + n*q - n + S(1)) + b*c*n*(p + S(1))), x), x)/(b*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule253)

    pattern254 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, m: Inequality(n, GreaterEqual, m - n + S(1), Greater, S(0))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule254 = ReplacementRule(pattern254, lambda p, x, q, a, d, b, c, e, n, m : -e**n*Int((e*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m - n + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(n*(p + S(1))*(-a*d + b*c)) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule254)

    pattern255 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule255 = ReplacementRule(pattern255, lambda p, x, q, a, d, b, c, e, n, m : -b*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*e*n*(p + S(1))*(-a*d + b*c)) + Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c*(m + S(1)) + b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule255)

    pattern256 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, q, m: RationalQ(m, p, q)), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule256 = ReplacementRule(pattern256, lambda p, x, q, a, d, b, c, e, n, m : -e**(-n)*n*Int((e*x)**(m + n)*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*d*q + b*c*p + b*d*x**n*(p + q), x), x)/(m + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + S(1))))
    rubi.add(rule256)

    pattern257 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda q, m: RationalQ(m, q)), CustomConstraint(lambda q: Greater(q, S(1))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule257 = ReplacementRule(pattern257, lambda p, x, q, a, d, b, c, e, n, m : c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(a*e*(m + S(1))) - e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*n*(a*d*(q + S(-1)) + b*c*(p + S(1))) + c*(m + S(1))*(-a*d + b*c) + d*x**n*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)), x), x)/(a*(m + S(1))))
    rubi.add(rule257)

    pattern258 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda q, m: RationalQ(m, q)), CustomConstraint(lambda q: Less(S(0), q, S(1))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule258 = ReplacementRule(pattern258, lambda p, x, q, a, d, b, c, e, n, m : -e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(b*c*(m + S(1)) + d*x**n*(b*n*(p + q + S(1)) + b*(m + S(1))) + n*(a*d*q + b*c*(p + S(1))), x), x)/(a*(m + S(1))) + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*(m + S(1))))
    rubi.add(rule258)

    pattern259 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule259 = ReplacementRule(pattern259, lambda p, x, q, a, d, b, c, e, n, m : n*Int((e*x)**m*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(m + n*(p + q) + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + n*(p + q) + S(1))))
    rubi.add(rule259)

    pattern260 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: Greater(q, S(1))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule260 = ReplacementRule(pattern260, lambda p, x, q, a, d, b, c, e, n, m : d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*e*(m + n*(p + q) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)) + x**n*(b*c*d*n*(p + q) + d*n*(q + S(-1))*(-a*d + b*c) + d*(m + S(1))*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1))))
    rubi.add(rule260)

    pattern261 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda q, m: RationalQ(m, q)), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda n, m: Greater(m - n + S(1), S(0))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule261 = ReplacementRule(pattern261, lambda p, x, q, a, d, b, c, e, n, m : -e**n*Int((e*x)**(m - n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(a*c*(m - n + S(1)) + x**n*(a*d*(m - n + S(1)) - n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1))) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*(m + n*(p + q) + S(1))))
    rubi.add(rule261)

    pattern262 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n, m: Greater(m - n + S(1), n)), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule262 = ReplacementRule(pattern262, lambda p, x, q, a, d, b, c, e, n, m : -e**(S(2)*n)*Int((e*x)**(m - S(2)*n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*(m - S(2)*n + S(1)) + x**n*(a*d*(m + n*(q + S(-1)) + S(1)) + b*c*(m + n*(p + S(-1)) + S(1))), x), x)/(b*d*(m + n*(p + q) + S(1))) + e**(S(2)*n + S(-1))*(e*x)**(m - S(2)*n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*d*(m + n*(p + q) + S(1))))
    rubi.add(rule262)

    pattern263 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule263 = ReplacementRule(pattern263, lambda p, x, q, a, d, b, c, e, n, m : -e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(a*d*q + b*c*p) + (a*d + b*c)*(m + n + S(1)), x), x)/(a*c*(m + S(1))) + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c*e*(m + S(1))))
    rubi.add(rule263)

    pattern264 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n, m: LessEqual(n, m, S(2)*n + S(-1))))
    rule264 = ReplacementRule(pattern264, lambda x, a, d, b, c, e, n, m : -a*e**n*Int((e*x)**(m - n)/(a + b*x**n), x)/(-a*d + b*c) + c*e**n*Int((e*x)**(m - n)/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule264)

    pattern265 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule265 = ReplacementRule(pattern265, lambda x, a, d, b, c, e, n, m : b*Int((e*x)**m/(a + b*x**n), x)/(-a*d + b*c) - d*Int((e*x)**m/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule265)

    pattern266 = Pattern(Integral(x_**m_/((a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, m: IntegersQ(m/S(2), n/S(2))), CustomConstraint(lambda n, m: Less(S(0), m - n + S(1), n)), CustomConstraint(lambda n: LessEqual(n, S(4))))
    rule266 = ReplacementRule(pattern266, lambda x, a, d, b, c, n, m : -a*Int(x**(m - n)/((a + b*x**n)*Sqrt(c + d*x**n)), x)/b + Int(x**(m - n)/Sqrt(c + d*x**n), x)/b)
    rubi.add(rule266)

    pattern267 = Pattern(Integral(x_**S(2)/((a_ + x_**S(4)*WC('b', S(1)))*sqrt(c_ + x_**S(4)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)))
    rule267 = ReplacementRule(pattern267, lambda x, a, d, b, c : With(List(Set(r, Numerator(Rt(-a/b, S(2)))), Set(s, Denominator(Rt(-a/b, S(2))))), -s*Int(S(1)/((r - s*x**S(2))*Sqrt(c + d*x**S(4))), x)/(S(2)*b) + s*Int(S(1)/((r + s*x**S(2))*Sqrt(c + d*x**S(4))), x)/(S(2)*b)))
    rubi.add(rule267)

    pattern268 = Pattern(Integral(x_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda b, c, d, a: ZeroQ(-a*d + S(4)*b*c)))
    rule268 = ReplacementRule(pattern268, lambda x, a, d, b, c : With(List(Set(q, Rt(d/c, S(3)))), -S(2)**(S(1)/3)*q*ArcTan(1/Sqrt(S(3)) + S(2)**(S(2)/3)*(Sqrt(c) - Sqrt(c + d*x**S(3)))/(q*x*Sqrt(S(3))*Sqrt(c)))/(S(6)*b*Sqrt(S(3))*Sqrt(c)) + S(2)**(S(1)/3)*q*ArcTan(1/Sqrt(S(3)) + S(2)**(S(2)/3)*(Sqrt(c) + Sqrt(c + d*x**S(3)))/(q*x*Sqrt(S(3))*Sqrt(c)))/(S(6)*b*Sqrt(S(3))*Sqrt(c)) + S(2)**(S(1)/3)*q*ArcTanh(Sqrt(c + d*x**S(3))/Sqrt(c))/(S(18)*b*Sqrt(c)) + S(2)**(S(1)/3)*q*Log(-S(2)**(S(1)/3)*q*x + S(1) - Sqrt(c + d*x**S(3))/Sqrt(c))/(S(12)*b*Sqrt(c)) - S(2)**(S(1)/3)*q*Log(-S(2)**(S(1)/3)*q*x + S(1) + Sqrt(c + d*x**S(3))/Sqrt(c))/(S(12)*b*Sqrt(c))))
    rubi.add(rule268)

    pattern269 = Pattern(Integral(x_**m_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda b, c, d, a: ZeroQ(-a*d + S(4)*b*c)), CustomConstraint(lambda m: PositiveIntegerQ(m/S(3) + S(-1)/3)))
    rule269 = ReplacementRule(pattern269, lambda x, a, d, b, c, m : -a*Int(x**(m + S(-3))/((a + b*x**S(3))*Sqrt(c + d*x**S(3))), x)/b + Int(x**(m + S(-3))/Sqrt(c + d*x**S(3)), x)/b)
    rubi.add(rule269)

    pattern270 = Pattern(Integral(x_**m_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda b, c, d, a: ZeroQ(-a*d + S(4)*b*c)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(3) + S(-1)/3)))
    rule270 = ReplacementRule(pattern270, lambda x, a, d, b, c, m : -b*Int(x**(m + S(3))/((a + b*x**S(3))*Sqrt(c + d*x**S(3))), x)/a + Int(x**m/Sqrt(c + d*x**S(3)), x)/a)
    rubi.add(rule270)

    pattern271 = Pattern(Integral(x_**S(2)*sqrt(c_ + x_**S(4)*WC('d', S(1)))/(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)))
    rule271 = ReplacementRule(pattern271, lambda x, a, d, b, c : d*Int(x**S(2)/Sqrt(c + d*x**S(4)), x)/b + (-a*d + b*c)*Int(x**S(2)/((a + b*x**S(4))*Sqrt(c + d*x**S(4))), x)/b)
    rubi.add(rule271)

    pattern272 = Pattern(Integral(x_**WC('m', S(1))*sqrt(c_ + x_**S(3)*WC('d', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda b, c, d, a: ZeroQ(-a*d + S(4)*b*c)), CustomConstraint(lambda m: IntegerQ(m/S(3) + S(-1)/3)))
    rule272 = ReplacementRule(pattern272, lambda x, a, d, b, c, m : d*Int(x**m/Sqrt(c + d*x**S(3)), x)/b + (-a*d + b*c)*Int(x**m/((a + b*x**S(3))*Sqrt(c + d*x**S(3))), x)/b)
    rubi.add(rule272)

    pattern273 = Pattern(Integral(x_**S(2)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda b, a: PosQ(b/a)), CustomConstraint(lambda c, d: PosQ(d/c)), CustomConstraint(lambda b, a, c, d: ~(SimplerSqrtQ(b/a, d/c))))
    rule273 = ReplacementRule(pattern273, lambda x, a, d, b, c : -c*Int(Sqrt(a + b*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/b + x*Sqrt(a + b*x**S(2))/(b*Sqrt(c + d*x**S(2))))
    rubi.add(rule273)

    pattern274 = Pattern(Integral(x_**n_/(sqrt(a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: EqQ(n, S(2)) | EqQ(n, S(4))), CustomConstraint(lambda a, d, b, c, n: ~(EqQ(n, S(2)) & SimplerSqrtQ(-b/a, -d/c))))
    rule274 = ReplacementRule(pattern274, lambda x, a, d, b, c, n : -a*Int(S(1)/(Sqrt(a + b*x**n)*Sqrt(c + d*x**n)), x)/b + Int(Sqrt(a + b*x**n)/Sqrt(c + d*x**n), x)/b)
    rubi.add(rule274)

    pattern275 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda n, p, q, m: IntegersQ(p + (m + S(1))/n, q)), CustomConstraint(lambda p: Less(S(-1), p, S(0))))
    rule275 = ReplacementRule(pattern275, lambda p, x, q, a, d, b, c, n, m : With(List(Set(k, Denominator(p))), a**(p + (m + S(1))/n)*k*Subst(Int(x**(k*(m + S(1))/n + S(-1))*(c - x**k*(-a*d + b*c))**q*(-b*x**k + S(1))**(-p - q + S(-1) - (m + S(1))/n), x), x, x**(n/k)*(a + b*x**n)**(-S(1)/k))/n))
    rubi.add(rule275)

    pattern276 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule276 = ReplacementRule(pattern276, lambda p, x, q, a, d, b, c, n, m : -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q, x), x, 1/x))
    rubi.add(rule276)

    pattern277 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)))
    rule277 = ReplacementRule(pattern277, lambda p, x, q, a, d, b, c, e, n, m : With(List(Set(g, Denominator(m))), -g*Subst(Int(x**(-g*(m + S(1)) + S(-1))*(a + b*e**(-n)*x**(-g*n))**p*(c + d*e**(-n)*x**(-g*n))**q, x), x, (e*x)**(-S(1)/g))/e))
    rubi.add(rule277)

    pattern278 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: ~(RationalQ(m))))
    rule278 = ReplacementRule(pattern278, lambda p, x, q, a, d, b, c, e, n, m : -(e*x)**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q, x), x, 1/x))
    rubi.add(rule278)

    pattern279 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: FractionQ(n)))
    rule279 = ReplacementRule(pattern279, lambda p, x, q, a, d, b, c, n, m : With(List(Set(g, Denominator(n))), g*Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + b*x**(g*n))**p*(c + d*x**(g*n))**q, x), x, x**(1/g))))
    rubi.add(rule279)

    pattern280 = Pattern(Integral((e_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n: FractionQ(n)))
    rule280 = ReplacementRule(pattern280, lambda p, x, q, a, d, b, c, e, n, m : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule280)

    pattern281 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, m: IntegerQ(Simplify(n/(m + S(1))))), CustomConstraint(lambda n: ~(IntegerQ(n))))
    rule281 = ReplacementRule(pattern281, lambda p, x, q, a, d, b, c, n, m : Subst(Int((a + b*x**Simplify(n/(m + S(1))))**p*(c + d*x**Simplify(n/(m + S(1))))**q, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule281)

    pattern282 = Pattern(Integral((e_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, m: IntegerQ(Simplify(n/(m + S(1))))), CustomConstraint(lambda n: ~(IntegerQ(n))))
    rule282 = ReplacementRule(pattern282, lambda p, x, q, a, d, b, c, e, n, m : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule282)

    pattern283 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda q: Greater(q, S(1))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule283 = ReplacementRule(pattern283, lambda p, x, q, a, d, b, c, e, n, m : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)) + d*x**n*(b*c*n*(p + S(1)) + (-a*d + b*c)*(m + n*(q + S(-1)) + S(1))), x), x)/(a*b*n*(p + S(1))) + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(a*d - b*c)/(a*b*e*n*(p + S(1))))
    rubi.add(rule283)

    pattern284 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda q: Less(S(0), q, S(1))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule284 = ReplacementRule(pattern284, lambda p, x, q, a, d, b, c, e, n, m : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m + n*(p + S(1)) + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*n*(p + S(1))))
    rubi.add(rule284)

    pattern285 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule285 = ReplacementRule(pattern285, lambda p, x, q, a, d, b, c, e, n, m : -b*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*e*n*(p + S(1))*(-a*d + b*c)) + Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c*(m + S(1)) + b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule285)

    pattern286 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule286 = ReplacementRule(pattern286, lambda p, x, q, a, d, b, c, e, n, m : n*Int((e*x)**m*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(m + n*(p + q) + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + n*(p + q) + S(1))))
    rubi.add(rule286)

    pattern287 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: Greater(q, S(1))), CustomConstraint(lambda p, n, q, x, a, d, b, e, c, m: IntBinomialQ(a, b, c, d, e, m, n, p, q, x)))
    rule287 = ReplacementRule(pattern287, lambda p, x, q, a, d, b, c, e, n, m : d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*e*(m + n*(p + q) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)) + x**n*(b*c*d*n*(p + q) + d*n*(q + S(-1))*(-a*d + b*c) + d*(m + S(1))*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1))))
    rubi.add(rule287)

    pattern288 = Pattern(Integral(x_**m_/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, m: ZeroQ(m - n) | ZeroQ(m - S(2)*n + S(1))))
    rule288 = ReplacementRule(pattern288, lambda x, a, d, b, c, n, m : -a*Int(x**(m - n)/(a + b*x**n), x)/(-a*d + b*c) + c*Int(x**(m - n)/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule288)

    pattern289 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)))
    rule289 = ReplacementRule(pattern289, lambda x, a, d, b, c, e, n, m : b*Int((e*x)**m/(a + b*x**n), x)/(-a*d + b*c) - d*Int((e*x)**m/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule289)

    pattern290 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, q, m: IntegersQ(m, p, q)), CustomConstraint(lambda p: GreaterEqual(p, S(-2))), CustomConstraint(lambda m, q: GreaterEqual(q, S(-2)) | (Equal(q, S(-3)) & IntegerQ(m/S(2) + S(-1)/2))))
    rule290 = ReplacementRule(pattern290, lambda p, x, q, a, d, b, c, e, n, m : Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule290)

    pattern291 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)), CustomConstraint(lambda q: IntegerQ(q)), CustomConstraint(lambda p, n: PosQ(n) | ~(IntegerQ(p))))
    rule291 = ReplacementRule(pattern291, lambda p, q, x, mn, a, d, b, c, n, m : Int(x**(m - n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule291)

    pattern292 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)), CustomConstraint(lambda q: ~(IntegerQ(q))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule292 = ReplacementRule(pattern292, lambda p, x, q, mn, a, d, b, c, n, m : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(m - n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule292)

    pattern293 = Pattern(Integral((e_*x_)**m_*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)))
    rule293 = ReplacementRule(pattern293, lambda p, q, x, mn, a, d, b, c, e, n, m : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**(-n))**q, x))
    rubi.add(rule293)

    pattern294 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: NonzeroQ(m - n + S(1))), CustomConstraint(lambda a: PositiveQ(a)), CustomConstraint(lambda c: PositiveQ(c)))
    rule294 = ReplacementRule(pattern294, lambda p, x, q, a, d, b, c, e, n, m : a**p*c**q*(e*x)**(m + S(1))*AppellF1((m + S(1))/n, -p, -q, S(1) + (m + S(1))/n, -b*x**n/a, -d*x**n/c)/(e*(m + S(1))))
    rubi.add(rule294)

    pattern295 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, c, d, a: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: NonzeroQ(m - n + S(1))), CustomConstraint(lambda a: ~(PositiveQ(a))))
    rule295 = ReplacementRule(pattern295, lambda p, x, q, a, d, b, c, e, n, m : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((e*x)**m*(S(1) + b*x**n/a)**p*(c + d*x**n)**q, x))
    rubi.add(rule295)

    pattern296 = Pattern(Integral(x_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda v, x: NonzeroQ(v - x)))
    rule296 = ReplacementRule(pattern296, lambda q, n, x, a, v, d, b, c, p, m : Coefficient(v, x, S(1))**(-m + S(-1))*Subst(Int(SimplifyIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(x - Coefficient(v, x, S(0)))**m, x), x), x, v))
    rubi.add(rule296)

    pattern297 = Pattern(Integral(u_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda x, v, u: LinearPairQ(u, v, x)))
    rule297 = ReplacementRule(pattern297, lambda q, n, x, a, v, u, d, b, c, p, m : u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule297)

    pattern298 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda non2, n: ZeroQ(-n/S(2) + non2)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p, a1, a2: IntegerQ(p) | (PositiveQ(a1) & PositiveQ(a2))))
    rule298 = ReplacementRule(pattern298, lambda p, b2, q, x, non2, u, d, b1, c, a1, n, a2 : Int(u*(c + d*x**n)**q*(a1*a2 + b1*b2*x**n)**p, x))
    rubi.add(rule298)

    pattern299 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)) + x_**WC('n2', S(1))*WC('e', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda non2, n: ZeroQ(-n/S(2) + non2)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda p, a1, a2: IntegerQ(p) | (PositiveQ(a1) & PositiveQ(a2))))
    rule299 = ReplacementRule(pattern299, lambda p, b2, n2, q, x, non2, u, d, b1, c, e, a1, n, a2 : Int(u*(a1*a2 + b1*b2*x**n)**p*(c + d*x**n + e*x**(S(2)*n))**q, x))
    rubi.add(rule299)

    pattern300 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**p_*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda non2, n: ZeroQ(-n/S(2) + non2)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)))
    rule300 = ReplacementRule(pattern300, lambda p, b2, q, x, non2, u, d, b1, c, a1, n, a2 : (a1 + b1*x**(n/S(2)))**FracPart(p)*(a2 + b2*x**(n/S(2)))**FracPart(p)*(a1*a2 + b1*b2*x**n)**(-FracPart(p))*Int(u*(c + d*x**n)**q*(a1*a2 + b1*b2*x**n)**p, x))
    rubi.add(rule300)

    pattern301 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)) + x_**WC('n2', S(1))*WC('e', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda non2, n: ZeroQ(-n/S(2) + non2)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda b2, a1, b1, a2: ZeroQ(a1*b2 + a2*b1)))
    rule301 = ReplacementRule(pattern301, lambda p, b2, n2, q, x, non2, u, d, b1, c, e, a1, n, a2 : (a1 + b1*x**(n/S(2)))**FracPart(p)*(a2 + b2*x**(n/S(2)))**FracPart(p)*(a1*a2 + b1*b2*x**n)**(-FracPart(p))*Int(u*(a1*a2 + b1*b2*x**n)**p*(c + d*x**n + e*x**(S(2)*n))**q, x))
    rubi.add(rule301)

    pattern302 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, q, r: PositiveIntegerQ(p, q, r)))
    rule302 = ReplacementRule(pattern302, lambda f, q, x, n, a, d, b, c, e, r, p : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x))
    rubi.add(rule302)

    pattern303 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule303 = ReplacementRule(pattern303, lambda f, x, a, d, b, c, e, n : (-a*f + b*e)*Int(1/(a + b*x**n), x)/(-a*d + b*c) - (-c*f + d*e)*Int(1/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule303)

    pattern304 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/((a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule304 = ReplacementRule(pattern304, lambda f, x, a, d, b, c, e, n : f*Int(1/Sqrt(c + d*x**n), x)/b + (-a*f + b*e)*Int(S(1)/((a + b*x**n)*Sqrt(c + d*x**n)), x)/b)
    rubi.add(rule304)

    pattern305 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/(sqrt(a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, d, b, c, n: ~(ZeroQ(n + S(-2)) & ((PosQ(b/a) & PosQ(d/c)) | (NegQ(b/a) & (PosQ(d/c) | (PositiveQ(a) & (~(PositiveQ(c)) | SimplerSqrtQ(-b/a, -d/c)))))))))
    rule305 = ReplacementRule(pattern305, lambda f, x, a, d, b, c, e, n : f*Int(Sqrt(a + b*x**n)/Sqrt(c + d*x**n), x)/b + (-a*f + b*e)*Int(S(1)/(Sqrt(a + b*x**n)*Sqrt(c + d*x**n)), x)/b)
    rubi.add(rule305)

    pattern306 = Pattern(Integral((e_ + x_**S(2)*WC('f', S(1)))/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: PosQ(b/a)), CustomConstraint(lambda c, d: PosQ(d/c)))
    rule306 = ReplacementRule(pattern306, lambda f, x, a, d, b, e, c : (-a*f + b*e)*Int(S(1)/(Sqrt(a + b*x**S(2))*Sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - (-c*f + d*e)*Int(Sqrt(a + b*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c))
    rubi.add(rule306)

    pattern307 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda q: Greater(q, S(0))))
    rule307 = ReplacementRule(pattern307, lambda p, f, q, x, a, d, b, c, e, n : x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(a*f - b*e)/(a*b*n*(p + S(1))) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(-a*f + b*e*n*(p + S(1)) + b*e) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(n*q + S(1))), x), x)/(a*b*n*(p + S(1))))
    rubi.add(rule307)

    pattern308 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule308 = ReplacementRule(pattern308, lambda p, f, q, x, a, d, b, c, e, n : x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(a*f - b*e)/(a*n*(p + S(1))*(-a*d + b*c)) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule308)

    pattern309 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda p, n, q: NonzeroQ(n*(p + q + S(1)) + S(1))))
    rule309 = ReplacementRule(pattern309, lambda f, q, x, n, a, d, b, c, e, p : f*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*(n*(p + q + S(1)) + S(1))) + Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(-a*f + b*e*n*(p + q + S(1)) + b*e) + x**n*(b*d*e*n*(p + q + S(1)) + d*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(n*(p + q + S(1)) + S(1))))
    rubi.add(rule309)

    pattern310 = Pattern(Integral((e_ + x_**S(4)*WC('f', S(1)))/((a_ + x_**S(4)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(4)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule310 = ReplacementRule(pattern310, lambda f, x, a, d, b, e, c : (-a*f + b*e)*Int((a + b*x**S(4))**(S(-3)/4), x)/(-a*d + b*c) - (-c*f + d*e)*Int((a + b*x**S(4))**(S(1)/4)/(c + d*x**S(4)), x)/(-a*d + b*c))
    rubi.add(rule310)

    pattern311 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule311 = ReplacementRule(pattern311, lambda p, f, x, a, d, b, c, e, n : f*Int((a + b*x**n)**p, x)/d + (-c*f + d*e)*Int((a + b*x**n)**p/(c + d*x**n), x)/d)
    rubi.add(rule311)

    pattern312 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)))
    rule312 = ReplacementRule(pattern312, lambda f, q, x, n, a, d, b, c, e, p : e*Int((a + b*x**n)**p*(c + d*x**n)**q, x) + f*Int(x**n*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule312)

    pattern313 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule313 = ReplacementRule(pattern313, lambda f, x, a, d, b, e, c : b*Int(S(1)/((a + b*x**S(2))*Sqrt(e + f*x**S(2))), x)/(-a*d + b*c) - d*Int(S(1)/((c + d*x**S(2))*Sqrt(e + f*x**S(2))), x)/(-a*d + b*c))
    rubi.add(rule313)

    pattern314 = Pattern(Integral(S(1)/(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, f, e, d: NonzeroQ(-c*f + d*e)))
    rule314 = ReplacementRule(pattern314, lambda f, x, d, e, c : -d*Int(S(1)/((c + d*x**S(2))*Sqrt(e + f*x**S(2))), x)/c + Int(S(1)/(x**S(2)*Sqrt(e + f*x**S(2))), x)/c)
    rubi.add(rule314)

    pattern315 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d: PositiveQ(d/c)), CustomConstraint(lambda f, e: PositiveQ(f/e)), CustomConstraint(lambda c, f, e, d: ~(SimplerSqrtQ(d/c, f/e))))
    rule315 = ReplacementRule(pattern315, lambda f, x, a, d, b, e, c : d*Int(Sqrt(e + f*x**S(2))/Sqrt(c + d*x**S(2)), x)/b + (-a*d + b*c)*Int(Sqrt(e + f*x**S(2))/((a + b*x**S(2))*Sqrt(c + d*x**S(2))), x)/b)
    rubi.add(rule315)

    pattern316 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda f, c, e, d: ~(SimplerSqrtQ(-f/e, -d/c))))
    rule316 = ReplacementRule(pattern316, lambda f, x, a, d, b, e, c : d*Int(Sqrt(e + f*x**S(2))/Sqrt(c + d*x**S(2)), x)/b + (-a*d + b*c)*Int(Sqrt(e + f*x**S(2))/((a + b*x**S(2))*Sqrt(c + d*x**S(2))), x)/b)
    rubi.add(rule316)

    pattern317 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d: PosQ(d/c)), CustomConstraint(lambda f, e: PosQ(f/e)), CustomConstraint(lambda c, f, e, d: ~(SimplerSqrtQ(d/c, f/e))))
    rule317 = ReplacementRule(pattern317, lambda f, x, a, d, b, e, c : b*Int(Sqrt(e + f*x**S(2))/((a + b*x**S(2))*Sqrt(c + d*x**S(2))), x)/(-a*f + b*e) - f*Int(S(1)/(Sqrt(c + d*x**S(2))*Sqrt(e + f*x**S(2))), x)/(-a*f + b*e))
    rubi.add(rule317)

    pattern318 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d: NegQ(d/c)), CustomConstraint(lambda c: PositiveQ(c)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda f, c, e, d: ~(NegQ(f/e) & SimplerSqrtQ(-f/e, -d/c))))
    rule318 = ReplacementRule(pattern318, lambda f, x, a, d, b, e, c : EllipticPi(b*c/(a*d), ArcSin(x*Rt(-d/c, S(2))), c*f/(d*e))/(a*Rt(-d/c, S(2))*Sqrt(c)*Sqrt(e)))
    rubi.add(rule318)

    pattern319 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c: ~(PositiveQ(c))))
    rule319 = ReplacementRule(pattern319, lambda f, x, a, d, b, e, c : Int(S(1)/((a + b*x**S(2))*Sqrt(S(1) + d*x**S(2)/c)*Sqrt(e + f*x**S(2))), x)*Sqrt(S(1) + d*x**S(2)/c)/Sqrt(c + d*x**S(2)))
    rubi.add(rule319)

    pattern320 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d: PosQ(d/c)))
    rule320 = ReplacementRule(pattern320, lambda f, x, a, d, b, e, c : c*EllipticPi(S(1) - b*c/(a*d), ArcTan(x*Rt(d/c, S(2))), -c*f/(d*e) + S(1))*Sqrt(e + f*x**S(2))/(a*e*Rt(d/c, S(2))*Sqrt(c*(e + f*x**S(2))/(e*(c + d*x**S(2))))*Sqrt(c + d*x**S(2))))
    rubi.add(rule320)

    pattern321 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d: NegQ(d/c)))
    rule321 = ReplacementRule(pattern321, lambda f, x, a, d, b, e, c : d*Int(S(1)/(Sqrt(c + d*x**S(2))*Sqrt(e + f*x**S(2))), x)/b + (-a*d + b*c)*Int(S(1)/((a + b*x**S(2))*Sqrt(c + d*x**S(2))*Sqrt(e + f*x**S(2))), x)/b)
    rubi.add(rule321)

    pattern322 = Pattern(Integral(sqrt(e_ + x_**S(2)*WC('f', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d: PosQ(d/c)), CustomConstraint(lambda f, e: PosQ(f/e)))
    rule322 = ReplacementRule(pattern322, lambda f, x, a, d, b, e, c : b*Int(Sqrt(e + f*x**S(2))/((a + b*x**S(2))*Sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - d*Int(Sqrt(e + f*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c))
    rubi.add(rule322)

    pattern323 = Pattern(Integral((e_ + x_**S(2)*WC('f', S(1)))**(S(3)/2)/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d: PosQ(d/c)), CustomConstraint(lambda f, e: PosQ(f/e)))
    rule323 = ReplacementRule(pattern323, lambda f, x, a, d, b, e, c : (-a*f + b*e)*Int(Sqrt(e + f*x**S(2))/((a + b*x**S(2))*Sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - (-c*f + d*e)*Int(Sqrt(e + f*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c))
    rubi.add(rule323)

    pattern324 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d: PosQ(d/c)), CustomConstraint(lambda f, e: PosQ(f/e)))
    rule324 = ReplacementRule(pattern324, lambda f, x, a, d, b, e, c : d*Int((-a*d + S(2)*b*c + b*d*x**S(2))*Sqrt(e + f*x**S(2))/Sqrt(c + d*x**S(2)), x)/b**S(2) + (-a*d + b*c)**S(2)*Int(Sqrt(e + f*x**S(2))/((a + b*x**S(2))*Sqrt(c + d*x**S(2))), x)/b**S(2))
    rubi.add(rule324)

    pattern325 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda q, r: RationalQ(q, r)), CustomConstraint(lambda q: Less(q, S(-1))), CustomConstraint(lambda r: Greater(r, S(1))))
    rule325 = ReplacementRule(pattern325, lambda f, x, q, a, d, b, e, r, c : b*(-a*f + b*e)*Int((c + d*x**S(2))**(q + S(2))*(e + f*x**S(2))**(r + S(-1))/(a + b*x**S(2)), x)/(-a*d + b*c)**S(2) - Int((c + d*x**S(2))**q*(e + f*x**S(2))**(r + S(-1))*(-a*d**S(2)*e - b*c**S(2)*f + S(2)*b*c*d*e + d**S(2)*x**S(2)*(-a*f + b*e)), x)/(-a*d + b*c)**S(2))
    rubi.add(rule325)

    pattern326 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: Greater(q, S(1))))
    rule326 = ReplacementRule(pattern326, lambda f, x, q, a, d, b, e, r, c : d*Int((c + d*x**S(2))**(q + S(-1))*(e + f*x**S(2))**r, x)/b + (-a*d + b*c)*Int((c + d*x**S(2))**(q + S(-1))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/b)
    rubi.add(rule326)

    pattern327 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: Less(q, S(-1))))
    rule327 = ReplacementRule(pattern327, lambda f, x, q, a, d, b, e, r, c : b**S(2)*Int((c + d*x**S(2))**(q + S(2))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/(-a*d + b*c)**S(2) - d*Int((c + d*x**S(2))**q*(e + f*x**S(2))**r*(-a*d + S(2)*b*c + b*d*x**S(2)), x)/(-a*d + b*c)**S(2))
    rubi.add(rule327)

    pattern328 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: LessEqual(q, S(-1))))
    rule328 = ReplacementRule(pattern328, lambda f, x, q, a, d, b, e, r, c : b*Int((c + d*x**S(2))**(q + S(1))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/(-a*d + b*c) - d*Int((c + d*x**S(2))**q*(e + f*x**S(2))**r, x)/(-a*d + b*c))
    rubi.add(rule328)

    pattern329 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule329 = ReplacementRule(pattern329, lambda f, x, a, d, b, e, c : x*Sqrt(c + d*x**S(2))*Sqrt(e + f*x**S(2))/(S(2)*a*(a + b*x**S(2))) + d*f*Int((a - b*x**S(2))/(Sqrt(c + d*x**S(2))*Sqrt(e + f*x**S(2))), x)/(S(2)*a*b**S(2)) + (-a**S(2)*d*f + b**S(2)*c*e)*Int(S(1)/((a + b*x**S(2))*Sqrt(c + d*x**S(2))*Sqrt(e + f*x**S(2))), x)/(S(2)*a*b**S(2)))
    rubi.add(rule329)

    pattern330 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**S(2)*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule330 = ReplacementRule(pattern330, lambda f, x, a, d, b, e, c : b**S(2)*x*Sqrt(c + d*x**S(2))*Sqrt(e + f*x**S(2))/(S(2)*a*(a + b*x**S(2))*(-a*d + b*c)*(-a*f + b*e)) - d*f*Int((a + b*x**S(2))/(Sqrt(c + d*x**S(2))*Sqrt(e + f*x**S(2))), x)/(S(2)*a*(-a*d + b*c)*(-a*f + b*e)) + (S(3)*a**S(2)*d*f - S(2)*a*b*(c*f + d*e) + b**S(2)*c*e)*Int(S(1)/((a + b*x**S(2))*Sqrt(c + d*x**S(2))*Sqrt(e + f*x**S(2))), x)/(S(2)*a*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule330)

    pattern331 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: Greater(q, S(0))))
    rule331 = ReplacementRule(pattern331, lambda p, f, x, q, a, d, b, c, e, r, n : d*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(e + f*x**n)**r, x)/b + (-a*d + b*c)*Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*(e + f*x**n)**r, x)/b)
    rubi.add(rule331)

    pattern332 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: LessEqual(q, S(-1))))
    rule332 = ReplacementRule(pattern332, lambda p, f, x, q, a, d, b, c, e, r, n : b*Int((a + b*x**n)**p*(c + d*x**n)**(q + S(1))*(e + f*x**n)**r, x)/(-a*d + b*c) - d*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(e + f*x**n)**r, x)/(-a*d + b*c))
    rubi.add(rule332)

    pattern333 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule333 = ReplacementRule(pattern333, lambda f, x, a, d, b, e, c : Sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*Sqrt(c + d*x**S(2))*Subst(Int(S(1)/(Sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)*Sqrt(S(1) - x**S(2)*(-a*f + b*e)/e)), x), x, x/Sqrt(a + b*x**S(2)))/(c*Sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*Sqrt(e + f*x**S(2))))
    rubi.add(rule333)

    pattern334 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule334 = ReplacementRule(pattern334, lambda f, x, a, d, b, e, c : a*Sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*Sqrt(c + d*x**S(2))*Subst(Int(S(1)/((-b*x**S(2) + S(1))*Sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)*Sqrt(S(1) - x**S(2)*(-a*f + b*e)/e)), x), x, x/Sqrt(a + b*x**S(2)))/(c*Sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*Sqrt(e + f*x**S(2))))
    rubi.add(rule334)

    pattern335 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))**(S(3)/2)*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule335 = ReplacementRule(pattern335, lambda f, x, a, d, b, e, c : Sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*Sqrt(c + d*x**S(2))*Subst(Int(Sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)/Sqrt(S(1) - x**S(2)*(-a*f + b*e)/e), x), x, x/Sqrt(a + b*x**S(2)))/(a*Sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*Sqrt(e + f*x**S(2))))
    rubi.add(rule335)

    pattern336 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/sqrt(e_ + x_**S(2)*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, f, e, d: PosQ((-c*f + d*e)/c)))
    rule336 = ReplacementRule(pattern336, lambda f, x, a, d, b, e, c : b*c*(-c*f + d*e)*Int(S(1)/(Sqrt(a + b*x**S(2))*Sqrt(c + d*x**S(2))*Sqrt(e + f*x**S(2))), x)/(S(2)*d*f) - c*(-c*f + d*e)*Int(Sqrt(a + b*x**S(2))/((c + d*x**S(2))**(S(3)/2)*Sqrt(e + f*x**S(2))), x)/(S(2)*f) + d*x*Sqrt(a + b*x**S(2))*Sqrt(e + f*x**S(2))/(S(2)*f*Sqrt(c + d*x**S(2))) - (-a*d*f - b*c*f + b*d*e)*Int(Sqrt(c + d*x**S(2))/(Sqrt(a + b*x**S(2))*Sqrt(e + f*x**S(2))), x)/(S(2)*d*f))
    rubi.add(rule336)

    pattern337 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/sqrt(e_ + x_**S(2)*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, f, e, d: NegQ((-c*f + d*e)/c)))
    rule337 = ReplacementRule(pattern337, lambda f, x, a, d, b, e, c : e*(-a*f + b*e)*Int(Sqrt(c + d*x**S(2))/((e + f*x**S(2))**(S(3)/2)*Sqrt(a + b*x**S(2))), x)/(S(2)*f) + x*Sqrt(a + b*x**S(2))*Sqrt(c + d*x**S(2))/(S(2)*Sqrt(e + f*x**S(2))) + (-a*f + b*e)*(-S(2)*c*f + d*e)*Int(S(1)/(Sqrt(a + b*x**S(2))*Sqrt(c + d*x**S(2))*Sqrt(e + f*x**S(2))), x)/(S(2)*f**S(2)) - (-a*d*f - b*c*f + b*d*e)*Int(Sqrt(e + f*x**S(2))/(Sqrt(a + b*x**S(2))*Sqrt(c + d*x**S(2))), x)/(S(2)*f**S(2)))
    rubi.add(rule337)

    pattern338 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/(e_ + x_**S(2)*WC('f', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule338 = ReplacementRule(pattern338, lambda f, x, a, d, b, e, c : b*Int(Sqrt(c + d*x**S(2))/(Sqrt(a + b*x**S(2))*Sqrt(e + f*x**S(2))), x)/f - (-a*f + b*e)*Int(Sqrt(c + d*x**S(2))/((e + f*x**S(2))**(S(3)/2)*Sqrt(a + b*x**S(2))), x)/f)
    rubi.add(rule338)

    pattern339 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule339 = ReplacementRule(pattern339, lambda p, f, x, q, a, d, b, c, e, r, n : With(List(Set(u, ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))), Condition(Int(u, x), SumQ(u))))
    rubi.add(rule339)

    pattern340 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule340 = ReplacementRule(pattern340, lambda p, f, x, q, a, d, b, c, e, r, n : -Subst(Int((a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r/x**S(2), x), x, 1/x))
    rubi.add(rule340)

    pattern341 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)))
    rule341 = ReplacementRule(pattern341, lambda f, q, x, n, a, d, b, c, e, r, p : Int((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule341)

    pattern342 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda v, u: ZeroQ(u - v)), CustomConstraint(lambda w, u: ZeroQ(u - w)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda x, u: NonzeroQ(u - x)))
    rule342 = ReplacementRule(pattern342, lambda f, q, n, x, a, w, v, u, d, b, c, e, r, p : Subst(Int((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule342)

    pattern343 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)), CustomConstraint(lambda q: IntegerQ(q)))
    rule343 = ReplacementRule(pattern343, lambda p, f, q, x, mn, a, d, b, c, e, r, n : Int(x**(-n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule343)

    pattern344 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda r: IntegerQ(r)))
    rule344 = ReplacementRule(pattern344, lambda p, f, q, x, mn, a, d, b, c, e, r, n : Int(x**(n*(p + r))*(c + d*x**(-n))**q*(a*x**(-n) + b)**p*(e*x**(-n) + f)**r, x))
    rubi.add(rule344)

    pattern345 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)), CustomConstraint(lambda q: ~(IntegerQ(q))))
    rule345 = ReplacementRule(pattern345, lambda p, f, x, q, mn, a, d, b, c, e, r, n : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(-n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule345)

    pattern346 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda f1, x: FreeQ(f1, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f2, x: FreeQ(f2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n2, n: ZeroQ(-n/S(2) + n2)), CustomConstraint(lambda e2, f1, f2, e1: ZeroQ(e1*f2 + e2*f1)), CustomConstraint(lambda e2, r, e1: IntegerQ(r) | (PositiveQ(e1) & PositiveQ(e2))))
    rule346 = ReplacementRule(pattern346, lambda n2, q, x, n, e2, a, d, f2, b, c, r, e1, f1, p : Int((a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule346)

    pattern347 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda f1, x: FreeQ(f1, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f2, x: FreeQ(f2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n2, n: ZeroQ(-n/S(2) + n2)), CustomConstraint(lambda e2, f1, f2, e1: ZeroQ(e1*f2 + e2*f1)))
    rule347 = ReplacementRule(pattern347, lambda n2, q, x, n, e2, a, d, f2, b, c, r, e1, f1, p : (e1 + f1*x**(n/S(2)))**FracPart(r)*(e2 + f2*x**(n/S(2)))**FracPart(r)*(e1*e2 + f1*f2*x**n)**(-FracPart(r))*Int((a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule347)

    pattern348 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda g, m: IntegerQ(m) | PositiveQ(g)), CustomConstraint(lambda n, m: IntegerQ(Simplify((m + S(1))/n))))
    rule348 = ReplacementRule(pattern348, lambda p, f, q, x, d, b, c, e, r, g, n, m : b**(-Simplify((m + S(1))/n) + S(1))*g**m*Subst(Int((b*x)**(p + Simplify((m + S(1))/n) + S(-1))*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n)
    rubi.add(rule348)

    pattern349 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda g, m: IntegerQ(m) | PositiveQ(g)), CustomConstraint(lambda n, m: ~(IntegerQ(Simplify((m + S(1))/n)))))
    rule349 = ReplacementRule(pattern349, lambda p, f, q, x, d, b, c, e, r, g, n, m : b**IntPart(p)*g**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p)*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule349)

    pattern350 = Pattern(Integral((g_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule350 = ReplacementRule(pattern350, lambda p, f, q, x, d, b, c, e, r, g, n, m : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule350)

    pattern351 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, q, r: PositiveIntegerQ(p + S(2), q, r)))
    rule351 = ReplacementRule(pattern351, lambda f, q, x, n, a, d, b, c, e, r, g, p, m : Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x))
    rubi.add(rule351)

    pattern352 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))))
    rule352 = ReplacementRule(pattern352, lambda f, q, x, n, a, d, b, c, e, r, p, m : Subst(Int((a + b*x)**p*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n)
    rubi.add(rule352)

    pattern353 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, q, r: IntegersQ(p, q, r)), CustomConstraint(lambda n: NegQ(n)))
    rule353 = ReplacementRule(pattern353, lambda f, q, x, n, a, d, b, c, e, r, p, m : Int(x**(m + n*(p + q + r))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q*(e*x**(-n) + f)**r, x))
    rubi.add(rule353)

    pattern354 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n, m: IntegerQ(Simplify((m + S(1))/n))))
    rule354 = ReplacementRule(pattern354, lambda f, q, x, n, a, d, b, c, e, r, p, m : Subst(Int(x**(Simplify((m + S(1))/n) + S(-1))*(a + b*x)**p*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n)
    rubi.add(rule354)

    pattern355 = Pattern(Integral((g_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n, m: IntegerQ(Simplify((m + S(1))/n))))
    rule355 = ReplacementRule(pattern355, lambda f, q, x, n, a, d, b, c, e, r, g, p, m : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule355)

    pattern356 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule356 = ReplacementRule(pattern356, lambda f, q, x, n, a, d, b, c, e, r, p, m : With(List(Set(k, GCD(m + S(1), n))), Condition(Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p*(c + d*x**(n/k))**q*(e + f*x**(n/k))**r, x), x, x**k)/k, Unequal(k, S(1)))))
    rubi.add(rule356)

    pattern357 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)))
    rule357 = ReplacementRule(pattern357, lambda p, f, x, q, a, d, b, c, e, r, g, n, m : With(List(Set(k, Denominator(m))), k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*g**(-n)*x**(k*n))**p*(c + d*g**(-n)*x**(k*n))**q*(e + f*g**(-n)*x**(k*n))**r, x), x, (g*x)**(1/k))/g))
    rubi.add(rule357)

    pattern358 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda f, q, a, d, b, e, c: ~(Equal(q, S(1)) & SimplerQ(-a*d + b*c, -a*f + b*e))))
    rule358 = ReplacementRule(pattern358, lambda p, f, q, x, a, d, b, c, e, g, n, m : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + S(1)) + (m + S(1))*(-a*f + b*e)) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(m + n*q + S(1))), x), x)/(a*b*n*(p + S(1))) + (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(a*f - b*e)/(a*b*g*n*(p + S(1))))
    rubi.add(rule358)

    pattern359 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, m: Greater(m - n + S(1), S(0))))
    rule359 = ReplacementRule(pattern359, lambda p, f, x, q, a, d, b, c, e, g, n, m : -g**n*Int((g*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(-a*f + b*e)*(m - n + S(1)) + x**n*(-b*n*(p + S(1))*(c*f - d*e) + d*(-a*f + b*e)*(m + n*q + S(1))), x), x)/(b*n*(p + S(1))*(-a*d + b*c)) + g**(n + S(-1))*(g*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(b*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule359)

    pattern360 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule360 = ReplacementRule(pattern360, lambda p, f, x, q, a, d, b, c, e, g, n, m : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(m + n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)) + (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(a*f - b*e)/(a*g*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule360)

    pattern361 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda q, m: RationalQ(m, q)), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda f, x, q, d, c, e, n: ~(Equal(q, S(1)) & SimplerQ(e + f*x**n, c + d*x**n))))
    rule361 = ReplacementRule(pattern361, lambda f, q, x, n, a, d, b, c, e, g, p, m : e*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*g*(m + S(1))) - g**(-n)*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + e*n*(a*d*q + b*c*(p + S(1))), x), x)/(a*(m + S(1))))
    rubi.add(rule361)

    pattern362 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda f, x, q, d, c, e, n: ~(Equal(q, S(1)) & SimplerQ(e + f*x**n, c + d*x**n))))
    rule362 = ReplacementRule(pattern362, lambda f, q, x, n, a, d, b, c, e, g, p, m : f*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*g*(m + n*(p + q + S(1)) + S(1))) + Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + x**n*(b*d*e*n*(p + q + S(1)) + d*(m + S(1))*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q + S(1)) + S(1))))
    rubi.add(rule362)

    pattern363 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n, m: Greater(m, n + S(-1))))
    rule363 = ReplacementRule(pattern363, lambda f, q, x, n, a, d, b, c, e, g, p, m : f*g**(n + S(-1))*(g*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*d*(m + n*(p + q + S(1)) + S(1))) - g**n*Int((g*x)**(m - n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*f*(m - n + S(1)) + x**n*(a*d*f*(m + n*q + S(1)) + b*(c*f*(m + n*p + S(1)) - d*e*(m + n*(p + q + S(1)) + S(1)))), x), x)/(b*d*(m + n*(p + q + S(1)) + S(1))))
    rubi.add(rule363)

    pattern364 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule364 = ReplacementRule(pattern364, lambda f, q, x, n, a, d, b, c, e, g, p, m : e*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c*g*(m + S(1))) + g**(-n)*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*f*(m + S(1)) - b*d*e*x**n*(m + n*(p + q + S(2)) + S(1)) - e*n*(a*d*q + b*c*p) - e*(a*d + b*c)*(m + n + S(1)), x), x)/(a*c*(m + S(1))))
    rubi.add(rule364)

    pattern365 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule365 = ReplacementRule(pattern365, lambda p, f, x, a, d, b, c, e, g, n, m : Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(e + f*x**n)/(c + d*x**n), x), x))
    rubi.add(rule365)

    pattern366 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule366 = ReplacementRule(pattern366, lambda f, q, x, n, a, d, b, c, e, g, p, m : e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x) + e**(-n)*f*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule366)

    pattern367 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda r: PositiveIntegerQ(r)))
    rule367 = ReplacementRule(pattern367, lambda f, q, x, n, a, d, b, c, e, r, g, p, m : e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**(r + S(-1)), x) + e**(-n)*f*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**(r + S(-1)), x))
    rubi.add(rule367)

    pattern368 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule368 = ReplacementRule(pattern368, lambda f, q, x, n, a, d, b, c, e, r, p, m : -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r, x), x, 1/x))
    rubi.add(rule368)

    pattern369 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)))
    rule369 = ReplacementRule(pattern369, lambda f, q, x, n, a, d, b, c, e, r, g, p, m : With(List(Set(k, Denominator(m))), -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*g**(-n)*x**(-k*n))**p*(c + d*g**(-n)*x**(-k*n))**q*(e + f*g**(-n)*x**(-k*n))**r, x), x, (g*x)**(-S(1)/k))/g))
    rubi.add(rule369)

    pattern370 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: ~(RationalQ(m))))
    rule370 = ReplacementRule(pattern370, lambda f, q, x, n, a, d, b, c, e, r, g, p, m : -(g*x)**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r, x), x, 1/x))
    rubi.add(rule370)

    pattern371 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n: FractionQ(n)))
    rule371 = ReplacementRule(pattern371, lambda f, q, x, n, a, d, b, c, e, r, p, m : With(List(Set(k, Denominator(n))), k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n))**p*(c + d*x**(k*n))**q*(e + f*x**(k*n))**r, x), x, x**(1/k))))
    rubi.add(rule371)

    pattern372 = Pattern(Integral((g_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n: FractionQ(n)))
    rule372 = ReplacementRule(pattern372, lambda f, q, x, n, a, d, b, c, e, r, g, p, m : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule372)

    pattern373 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n, m: IntegerQ(Simplify(n/(m + S(1))))))
    rule373 = ReplacementRule(pattern373, lambda f, q, x, n, a, d, b, c, e, r, p, m : Subst(Int((a + b*x**Simplify(n/(m + S(1))))**p*(c + d*x**Simplify(n/(m + S(1))))**q*(e + f*x**Simplify(n/(m + S(1))))**r, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule373)

    pattern374 = Pattern(Integral((g_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n, m: IntegerQ(Simplify(n/(m + S(1))))))
    rule374 = ReplacementRule(pattern374, lambda f, q, x, n, a, d, b, c, e, r, g, p, m : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule374)

    pattern375 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, q: RationalQ(p, q)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda f, q, a, d, b, e, c: ~(Equal(q, S(1)) & SimplerQ(-a*d + b*c, -a*f + b*e))))
    rule375 = ReplacementRule(pattern375, lambda p, f, q, x, a, d, b, c, e, g, n, m : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + S(1)) + (m + S(1))*(-a*f + b*e)) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(m + n*q + S(1))), x), x)/(a*b*n*(p + S(1))) + (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(a*f - b*e)/(a*b*g*n*(p + S(1))))
    rubi.add(rule375)

    pattern376 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule376 = ReplacementRule(pattern376, lambda p, f, x, q, a, d, b, c, e, g, n, m : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(m + n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)) + (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(a*f - b*e)/(a*g*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule376)

    pattern377 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q: RationalQ(q)), CustomConstraint(lambda q: Greater(q, S(0))), CustomConstraint(lambda f, x, q, d, c, e, n: ~(Equal(q, S(1)) & SimplerQ(e + f*x**n, c + d*x**n))))
    rule377 = ReplacementRule(pattern377, lambda f, q, x, n, a, d, b, c, e, g, p, m : f*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*g*(m + n*(p + q + S(1)) + S(1))) + Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + x**n*(b*d*e*n*(p + q + S(1)) + d*(m + S(1))*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q + S(1)) + S(1))))
    rubi.add(rule377)

    pattern378 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule378 = ReplacementRule(pattern378, lambda p, f, x, a, d, b, c, e, g, n, m : Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(e + f*x**n)/(c + d*x**n), x), x))
    rubi.add(rule378)

    pattern379 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)))
    rule379 = ReplacementRule(pattern379, lambda p, f, x, q, a, d, b, c, e, g, n, m : e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x) + f*x**(-m)*(g*x)**m*Int(x**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule379)

    pattern380 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)), CustomConstraint(lambda q: IntegerQ(q)))
    rule380 = ReplacementRule(pattern380, lambda p, f, q, x, mn, a, d, b, c, e, r, n, m : Int(x**(m - n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule380)

    pattern381 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda r: IntegerQ(r)))
    rule381 = ReplacementRule(pattern381, lambda p, f, q, x, mn, a, d, b, c, e, r, n, m : Int(x**(m + n*(p + r))*(c + d*x**(-n))**q*(a*x**(-n) + b)**p*(e*x**(-n) + f)**r, x))
    rubi.add(rule381)

    pattern382 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)), CustomConstraint(lambda q: ~(IntegerQ(q))))
    rule382 = ReplacementRule(pattern382, lambda p, f, x, q, mn, a, d, b, c, e, r, n, m : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(m - n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule382)

    pattern383 = Pattern(Integral((g_*x_)**m_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda mn, n: EqQ(mn, -n)))
    rule383 = ReplacementRule(pattern383, lambda p, f, q, x, mn, a, d, b, c, e, r, g, n, m : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**(-n))**q*(e + f*x**n)**r, x))
    rubi.add(rule383)

    pattern384 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)))
    rule384 = ReplacementRule(pattern384, lambda f, q, x, n, a, d, b, c, e, r, g, p, m : Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule384)

    pattern385 = Pattern(Integral(u_**WC('m', S(1))*(e_ + v_**n_*WC('f', S(1)))**WC('r', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda x, v, u: LinearPairQ(u, v, x)))
    rule385 = ReplacementRule(pattern385, lambda f, q, n, x, a, v, u, d, b, c, e, r, p, m : u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule385)

    pattern386 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda f1, x: FreeQ(f1, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f2, x: FreeQ(f2, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n2, n: ZeroQ(-n/S(2) + n2)), CustomConstraint(lambda e2, f1, f2, e1: ZeroQ(e1*f2 + e2*f1)), CustomConstraint(lambda e2, r, e1: IntegerQ(r) | (PositiveQ(e1) & PositiveQ(e2))))
    rule386 = ReplacementRule(pattern386, lambda n2, q, x, n, e2, a, d, f2, b, c, r, e1, f1, g, p, m : Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule386)

    pattern387 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda f1, x: FreeQ(f1, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f2, x: FreeQ(f2, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda n2, n: ZeroQ(-n/S(2) + n2)), CustomConstraint(lambda e2, f1, f2, e1: ZeroQ(e1*f2 + e2*f1)))
    rule387 = ReplacementRule(pattern387, lambda n2, q, x, n, e2, a, d, f2, b, c, r, e1, f1, g, p, m : (e1 + f1*x**(n/S(2)))**FracPart(r)*(e2 + f2*x**(n/S(2)))**FracPart(r)*(e1*e2 + f1*f2*x**n)**(-FracPart(r))*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule387)

    return rubi
