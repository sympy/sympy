
from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (SubstFor, Dist, Int, Set, With, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcCsch, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct, SplitSum, Complex, UnsameQ, _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify, _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum, _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux, TrigSimplifyAux)
    from sympy import Integral, S, sqrt
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)


    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3',' Pq', 'Pm', ' Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]
    d, p, q, r, s, v, mn, gcd, P, Q, lst = symbols('d p q r s v mn gcd P Q lst')

    _UseGamma = False

def miscellaneous_trig(rubi):

    pattern1 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule1 = ReplacementRule(pattern1, lambda m, d, n, u, c, a, x, b : (c*Tan(a + b*x))**m*(d*Cos(a + b*x))**m*(d*Sin(a + b*x))**(-m)*Int((d*Cos(a + b*x))**(-m)*(d*Sin(a + b*x))**(m + n)*ActivateTrig(u), x))
    rubi.add(rule1)

    pattern2 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule2 = ReplacementRule(pattern2, lambda m, d, n, u, c, a, x, b : (c*Tan(a + b*x))**m*(d*Cos(a + b*x))**m*(d*Sin(a + b*x))**(-m)*Int((d*Cos(a + b*x))**(-m + n)*(d*Sin(a + b*x))**m*ActivateTrig(u), x))
    rubi.add(rule2)

    pattern3 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule3 = ReplacementRule(pattern3, lambda m, d, n, u, c, a, x, b : (c*Cot(a + b*x))**m*(d*Cos(a + b*x))**(-m)*(d*Sin(a + b*x))**m*Int((d*Cos(a + b*x))**m*(d*Sin(a + b*x))**(-m + n)*ActivateTrig(u), x))
    rubi.add(rule3)

    pattern4 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule4 = ReplacementRule(pattern4, lambda m, d, n, u, c, a, x, b : (c*Cot(a + b*x))**m*(d*Cos(a + b*x))**(-m)*(d*Sin(a + b*x))**m*Int((d*Cos(a + b*x))**(m + n)*(d*Sin(a + b*x))**(-m)*ActivateTrig(u), x))
    rubi.add(rule4)

    pattern5 = Pattern(Integral(u_*(WC('c', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule5 = ReplacementRule(pattern5, lambda m, d, n, u, c, a, x, b : (c*Csc(a + b*x))**m*(d*Sin(a + b*x))**m*Int((d*Sin(a + b*x))**(-m + n)*ActivateTrig(u), x))
    rubi.add(rule5)

    pattern6 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule6 = ReplacementRule(pattern6, lambda m, b, c, a, x, u : (c*Cos(a + b*x))**m*(c*Sin(a + b*x))**(-m)*(c*Tan(a + b*x))**m*Int((c*Cos(a + b*x))**(-m)*(c*Sin(a + b*x))**m*ActivateTrig(u), x))
    rubi.add(rule6)

    pattern7 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule7 = ReplacementRule(pattern7, lambda a, m, u, c, x, b : (c*Cos(a + b*x))**(-m)*(c*Cot(a + b*x))**m*(c*Sin(a + b*x))**m*Int((c*Cos(a + b*x))**m*(c*Sin(a + b*x))**(-m)*ActivateTrig(u), x))
    rubi.add(rule7)

    pattern8 = Pattern(Integral(u_*(WC('c', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule8 = ReplacementRule(pattern8, lambda m, b, c, a, x, u : (c*Cos(a + b*x))**m*(c*Sec(a + b*x))**m*Int((c*Cos(a + b*x))**(-m)*ActivateTrig(u), x))
    rubi.add(rule8)

    pattern9 = Pattern(Integral(u_*(WC('c', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule9 = ReplacementRule(pattern9, lambda a, m, u, c, x, b : (c*Csc(a + b*x))**m*(c*Sin(a + b*x))**m*Int((c*Sin(a + b*x))**(-m)*ActivateTrig(u), x))
    rubi.add(rule9)

    pattern10 = Pattern(Integral(u_*(WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule10 = ReplacementRule(pattern10, lambda n, u, A, c, a, x, b, B : c*Int((c*Sin(a + b*x))**(n + S(-1))*(A*Sin(a + b*x) + B)*ActivateTrig(u), x))
    rubi.add(rule10)

    pattern11 = Pattern(Integral(u_*(WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule11 = ReplacementRule(pattern11, lambda n, u, A, c, a, x, b, B : c*Int((c*Cos(a + b*x))**(n + S(-1))*(A*Cos(a + b*x) + B)*ActivateTrig(u), x))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(u_*(A_ + WC('B', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule12 = ReplacementRule(pattern12, lambda a, u, A, x, b, B : Int((A*Sin(a + b*x) + B)*ActivateTrig(u)/Sin(a + b*x), x))
    rubi.add(rule12)

    pattern13 = Pattern(Integral(u_*(A_ + WC('B', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule13 = ReplacementRule(pattern13, lambda b, A, a, x, u, B : Int((A*Cos(a + b*x) + B)*ActivateTrig(u)/Cos(a + b*x), x))
    rubi.add(rule13)

    pattern14 = Pattern(Integral((WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule14 = ReplacementRule(pattern14, lambda n, u, C, A, c, a, x, b, B : c**S(2)*Int((c*Sin(a + b*x))**(n + S(-2))*(A*Sin(a + b*x)**S(2) + B*Sin(a + b*x) + C)*ActivateTrig(u), x))
    rubi.add(rule14)

    pattern15 = Pattern(Integral((WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule15 = ReplacementRule(pattern15, lambda n, b, C, A, c, a, x, u, B : c**S(2)*Int((c*Cos(a + b*x))**(n + S(-2))*(A*Cos(a + b*x)**S(2) + B*Cos(a + b*x) + C)*ActivateTrig(u), x))
    rubi.add(rule15)

    pattern16 = Pattern(Integral((WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule16 = ReplacementRule(pattern16, lambda n, b, C, A, c, a, x, u : c**S(2)*Int((c*Sin(a + b*x))**(n + S(-2))*(A*Sin(a + b*x)**S(2) + C)*ActivateTrig(u), x))
    rubi.add(rule16)

    pattern17 = Pattern(Integral((WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule17 = ReplacementRule(pattern17, lambda n, b, C, A, c, a, x, u : c**S(2)*Int((c*Cos(a + b*x))**(n + S(-2))*(A*Cos(a + b*x)**S(2) + C)*ActivateTrig(u), x))
    rubi.add(rule17)

    pattern18 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule18 = ReplacementRule(pattern18, lambda u, C, A, a, x, b, B : Int((A*Sin(a + b*x)**S(2) + B*Sin(a + b*x) + C)*ActivateTrig(u)/Sin(a + b*x)**S(2), x))
    rubi.add(rule18)

    pattern19 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule19 = ReplacementRule(pattern19, lambda u, C, A, a, x, b, B : Int((A*Cos(a + b*x)**S(2) + B*Cos(a + b*x) + C)*ActivateTrig(u)/Cos(a + b*x)**S(2), x))
    rubi.add(rule19)

    pattern20 = Pattern(Integral(u_*(A_ + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule20 = ReplacementRule(pattern20, lambda a, u, C, A, x, b : Int((A*Sin(a + b*x)**S(2) + C)*ActivateTrig(u)/Sin(a + b*x)**S(2), x))
    rubi.add(rule20)

    pattern21 = Pattern(Integral(u_*(A_ + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownSineIntegrandQ(u, x)))
    rule21 = ReplacementRule(pattern21, lambda b, C, A, a, x, u : Int((A*Cos(a + b*x)**S(2) + C)*ActivateTrig(u)/Cos(a + b*x)**S(2), x))
    rubi.add(rule21)

    pattern22 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)))
    rule22 = ReplacementRule(pattern22, lambda u, C, A, a, x, b, B : Int((A*Sin(a + b*x) + B*Sin(a + b*x)**S(2) + C)*ActivateTrig(u)/Sin(a + b*x), x))
    rubi.add(rule22)

    pattern23 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)))
    rule23 = ReplacementRule(pattern23, lambda u, C, A, a, x, b, B : Int((A*Cos(a + b*x) + B*Cos(a + b*x)**S(2) + C)*ActivateTrig(u)/Cos(a + b*x), x))
    rubi.add(rule23)

    pattern24 = Pattern(Integral(u_*(WC('A', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n1, n: ZeroQ(-n + n1 + S(-1))), CustomConstraint(lambda n, n2: ZeroQ(-n + n2 + S(-2))))
    rule24 = ReplacementRule(pattern24, lambda n1, n2, n, u, C, A, a, x, b, B : Int((A + B*Sin(a + b*x) + C*Sin(a + b*x)**S(2))*ActivateTrig(u)*Sin(a + b*x)**n, x))
    rubi.add(rule24)

    pattern25 = Pattern(Integral(u_*(WC('A', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n1, n: ZeroQ(-n + n1 + S(-1))), CustomConstraint(lambda n, n2: ZeroQ(-n + n2 + S(-2))))
    rule25 = ReplacementRule(pattern25, lambda n1, n2, n, u, C, A, a, x, b, B : Int((A + B*Cos(a + b*x) + C*Cos(a + b*x)**S(2))*ActivateTrig(u)*Cos(a + b*x)**n, x))
    rubi.add(rule25)

    pattern26 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownTangentIntegrandQ(u, x)))
    rule26 = ReplacementRule(pattern26, lambda m, d, n, u, c, a, x, b : (c*Cot(a + b*x))**m*(d*Tan(a + b*x))**m*Int((d*Tan(a + b*x))**(-m + n)*ActivateTrig(u), x))
    rubi.add(rule26)

    pattern27 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownCotangentIntegrandQ(u, x)))
    rule27 = ReplacementRule(pattern27, lambda m, d, n, u, c, a, x, b : (c*Tan(a + b*x))**m*(d*Cos(a + b*x))**m*(d*Sin(a + b*x))**(-m)*Int((d*Cos(a + b*x))**(-m + n)*(d*Sin(a + b*x))**m*ActivateTrig(u), x))
    rubi.add(rule27)

    pattern28 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda u, x: KnownTangentIntegrandQ(u, x)))
    rule28 = ReplacementRule(pattern28, lambda a, m, u, c, x, b : (c*Cot(a + b*x))**m*(c*Tan(a + b*x))**m*Int((c*Tan(a + b*x))**(-m)*ActivateTrig(u), x))
    rubi.add(rule28)

    pattern29 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda u, x: KnownCotangentIntegrandQ(u, x)))
    rule29 = ReplacementRule(pattern29, lambda m, b, c, a, x, u : (c*Cot(a + b*x))**m*(c*Tan(a + b*x))**m*Int((c*Cot(a + b*x))**(-m)*ActivateTrig(u), x))
    rubi.add(rule29)

    pattern30 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownTangentIntegrandQ(u, x)))
    rule30 = ReplacementRule(pattern30, lambda n, u, A, c, a, x, b, B : c*Int((c*Tan(a + b*x))**(n + S(-1))*(A*Tan(a + b*x) + B)*ActivateTrig(u), x))
    rubi.add(rule30)

    pattern31 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownCotangentIntegrandQ(u, x)))
    rule31 = ReplacementRule(pattern31, lambda n, u, A, c, a, x, b, B : c*Int((c*Cot(a + b*x))**(n + S(-1))*(A*Cot(a + b*x) + B)*ActivateTrig(u), x))
    rubi.add(rule31)

    pattern32 = Pattern(Integral(u_*(A_ + WC('B', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda u, x: KnownTangentIntegrandQ(u, x)))
    rule32 = ReplacementRule(pattern32, lambda a, u, A, x, b, B : Int((A*Tan(a + b*x) + B)*ActivateTrig(u)/Tan(a + b*x), x))
    rubi.add(rule32)

    pattern33 = Pattern(Integral(u_*(A_ + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda u, x: KnownCotangentIntegrandQ(u, x)))
    rule33 = ReplacementRule(pattern33, lambda b, A, a, x, u, B : Int((A*Cot(a + b*x) + B)*ActivateTrig(u)/Cot(a + b*x), x))
    rubi.add(rule33)

    pattern34 = Pattern(Integral((WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownTangentIntegrandQ(u, x)))
    rule34 = ReplacementRule(pattern34, lambda n, u, C, A, c, a, x, b, B : c**S(2)*Int((c*Tan(a + b*x))**(n + S(-2))*(A*Tan(a + b*x)**S(2) + B*Tan(a + b*x) + C)*ActivateTrig(u), x))
    rubi.add(rule34)

    pattern35 = Pattern(Integral((WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownCotangentIntegrandQ(u, x)))
    rule35 = ReplacementRule(pattern35, lambda n, b, C, A, c, a, x, u, B : c**S(2)*Int((c*Cot(a + b*x))**(n + S(-2))*(A*Cot(a + b*x)**S(2) + B*Cot(a + b*x) + C)*ActivateTrig(u), x))
    rubi.add(rule35)

    pattern36 = Pattern(Integral((WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownTangentIntegrandQ(u, x)))
    rule36 = ReplacementRule(pattern36, lambda n, b, C, A, c, a, x, u : c**S(2)*Int((c*Tan(a + b*x))**(n + S(-2))*(A*Tan(a + b*x)**S(2) + C)*ActivateTrig(u), x))
    rubi.add(rule36)

    pattern37 = Pattern(Integral((WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownCotangentIntegrandQ(u, x)))
    rule37 = ReplacementRule(pattern37, lambda n, b, C, A, c, a, x, u : c**S(2)*Int((c*Cot(a + b*x))**(n + S(-2))*(A*Cot(a + b*x)**S(2) + C)*ActivateTrig(u), x))
    rubi.add(rule37)

    pattern38 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownTangentIntegrandQ(u, x)))
    rule38 = ReplacementRule(pattern38, lambda u, C, A, a, x, b, B : Int((A*Tan(a + b*x)**S(2) + B*Tan(a + b*x) + C)*ActivateTrig(u)/Tan(a + b*x)**S(2), x))
    rubi.add(rule38)

    pattern39 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownCotangentIntegrandQ(u, x)))
    rule39 = ReplacementRule(pattern39, lambda u, C, A, a, x, b, B : Int((A*Cot(a + b*x)**S(2) + B*Cot(a + b*x) + C)*ActivateTrig(u)/Cot(a + b*x)**S(2), x))
    rubi.add(rule39)

    pattern40 = Pattern(Integral(u_*(A_ + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownTangentIntegrandQ(u, x)))
    rule40 = ReplacementRule(pattern40, lambda a, u, C, A, x, b : Int((A*Tan(a + b*x)**S(2) + C)*ActivateTrig(u)/Tan(a + b*x)**S(2), x))
    rubi.add(rule40)

    pattern41 = Pattern(Integral(u_*(A_ + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownCotangentIntegrandQ(u, x)))
    rule41 = ReplacementRule(pattern41, lambda b, C, A, a, x, u : Int((A*Cot(a + b*x)**S(2) + C)*ActivateTrig(u)/Cot(a + b*x)**S(2), x))
    rubi.add(rule41)

    pattern42 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)))
    rule42 = ReplacementRule(pattern42, lambda u, C, A, a, x, b, B : Int((A*Tan(a + b*x) + B*Tan(a + b*x)**S(2) + C)*ActivateTrig(u)/Tan(a + b*x), x))
    rubi.add(rule42)

    pattern43 = Pattern(Integral(u_*(WC('A', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n1, n: ZeroQ(-n + n1 + S(-1))), CustomConstraint(lambda n, n2: ZeroQ(-n + n2 + S(-2))))
    rule43 = ReplacementRule(pattern43, lambda n1, n2, n, u, C, A, a, x, b, B : Int((A + B*Tan(a + b*x) + C*Tan(a + b*x)**S(2))*ActivateTrig(u)*Tan(a + b*x)**n, x))
    rubi.add(rule43)

    pattern44 = Pattern(Integral(u_*(WC('A', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n1, n: ZeroQ(-n + n1 + S(-1))), CustomConstraint(lambda n, n2: ZeroQ(-n + n2 + S(-2))))
    rule44 = ReplacementRule(pattern44, lambda n1, n2, n, u, C, A, a, x, b, B : Int((A + B*Cot(a + b*x) + C*Cot(a + b*x)**S(2))*ActivateTrig(u)*Cot(a + b*x)**n, x))
    rubi.add(rule44)

    pattern45 = Pattern(Integral(u_*(WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule45 = ReplacementRule(pattern45, lambda m, d, n, u, c, a, x, b : (c*Sin(a + b*x))**m*(d*Csc(a + b*x))**m*Int((d*Csc(a + b*x))**(-m + n)*ActivateTrig(u), x))
    rubi.add(rule45)

    pattern46 = Pattern(Integral(u_*(WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule46 = ReplacementRule(pattern46, lambda m, d, n, u, c, a, x, b : (c*Cos(a + b*x))**m*(d*Sec(a + b*x))**m*Int((d*Sec(a + b*x))**(-m + n)*ActivateTrig(u), x))
    rubi.add(rule46)

    pattern47 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule47 = ReplacementRule(pattern47, lambda m, d, n, u, c, a, x, b : (c*Tan(a + b*x))**m*(d*Csc(a + b*x))**m*(d*Sec(a + b*x))**(-m)*Int((d*Csc(a + b*x))**(-m)*(d*Sec(a + b*x))**(m + n)*ActivateTrig(u), x))
    rubi.add(rule47)

    pattern48 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule48 = ReplacementRule(pattern48, lambda m, d, n, u, c, a, x, b : (c*Tan(a + b*x))**m*(d*Csc(a + b*x))**m*(d*Sec(a + b*x))**(-m)*Int((d*Csc(a + b*x))**(-m + n)*(d*Sec(a + b*x))**m*ActivateTrig(u), x))
    rubi.add(rule48)

    pattern49 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule49 = ReplacementRule(pattern49, lambda m, d, n, u, c, a, x, b : (c*Cot(a + b*x))**m*(d*Csc(a + b*x))**(-m)*(d*Sec(a + b*x))**m*Int((d*Csc(a + b*x))**m*(d*Sec(a + b*x))**(-m + n)*ActivateTrig(u), x))
    rubi.add(rule49)

    pattern50 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule50 = ReplacementRule(pattern50, lambda m, d, n, u, c, a, x, b : (c*Cot(a + b*x))**m*(d*Csc(a + b*x))**(-m)*(d*Sec(a + b*x))**m*Int((d*Csc(a + b*x))**(m + n)*(d*Sec(a + b*x))**(-m)*ActivateTrig(u), x))
    rubi.add(rule50)

    pattern51 = Pattern(Integral(u_*(WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule51 = ReplacementRule(pattern51, lambda m, b, c, a, x, u : (c*Csc(a + b*x))**m*(c*Sin(a + b*x))**m*Int((c*Csc(a + b*x))**(-m)*ActivateTrig(u), x))
    rubi.add(rule51)

    pattern52 = Pattern(Integral(u_*(WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule52 = ReplacementRule(pattern52, lambda a, m, u, c, x, b : (c*Cos(a + b*x))**m*(c*Sec(a + b*x))**m*Int((c*Sec(a + b*x))**(-m)*ActivateTrig(u), x))
    rubi.add(rule52)

    pattern53 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule53 = ReplacementRule(pattern53, lambda m, b, c, a, x, u : (c*Csc(a + b*x))**m*(c*Sec(a + b*x))**(-m)*(c*Tan(a + b*x))**m*Int((c*Csc(a + b*x))**(-m)*(c*Sec(a + b*x))**m*ActivateTrig(u), x))
    rubi.add(rule53)

    pattern54 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule54 = ReplacementRule(pattern54, lambda a, m, u, c, x, b : (c*Cot(a + b*x))**m*(c*Csc(a + b*x))**(-m)*(c*Sec(a + b*x))**m*Int((c*Csc(a + b*x))**m*(c*Sec(a + b*x))**(-m)*ActivateTrig(u), x))
    rubi.add(rule54)

    pattern55 = Pattern(Integral(u_*(WC('c', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule55 = ReplacementRule(pattern55, lambda n, u, A, c, a, x, b, B : c*Int((c*Sec(a + b*x))**(n + S(-1))*(A*Sec(a + b*x) + B)*ActivateTrig(u), x))
    rubi.add(rule55)

    pattern56 = Pattern(Integral(u_*(WC('c', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule56 = ReplacementRule(pattern56, lambda n, u, A, c, a, x, b, B : c*Int((c*Csc(a + b*x))**(n + S(-1))*(A*Csc(a + b*x) + B)*ActivateTrig(u), x))
    rubi.add(rule56)

    pattern57 = Pattern(Integral(u_*(A_ + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule57 = ReplacementRule(pattern57, lambda a, u, A, x, b, B : Int((A*Sec(a + b*x) + B)*ActivateTrig(u)/Sec(a + b*x), x))
    rubi.add(rule57)

    pattern58 = Pattern(Integral(u_*(A_ + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule58 = ReplacementRule(pattern58, lambda b, A, a, x, u, B : Int((A*Csc(a + b*x) + B)*ActivateTrig(u)/Csc(a + b*x), x))
    rubi.add(rule58)

    pattern59 = Pattern(Integral((WC('c', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule59 = ReplacementRule(pattern59, lambda n, u, C, A, c, a, x, b, B : c**S(2)*Int((c*Sec(a + b*x))**(n + S(-2))*(A*Sec(a + b*x)**S(2) + B*Sec(a + b*x) + C)*ActivateTrig(u), x))
    rubi.add(rule59)

    pattern60 = Pattern(Integral((WC('c', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule60 = ReplacementRule(pattern60, lambda n, b, C, A, c, a, x, u, B : c**S(2)*Int((c*Csc(a + b*x))**(n + S(-2))*(A*Csc(a + b*x)**S(2) + B*Csc(a + b*x) + C)*ActivateTrig(u), x))
    rubi.add(rule60)

    pattern61 = Pattern(Integral((WC('c', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule61 = ReplacementRule(pattern61, lambda n, b, C, A, c, a, x, u : c**S(2)*Int((c*Sec(a + b*x))**(n + S(-2))*(A*Sec(a + b*x)**S(2) + C)*ActivateTrig(u), x))
    rubi.add(rule61)

    pattern62 = Pattern(Integral((WC('c', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule62 = ReplacementRule(pattern62, lambda n, b, C, A, c, a, x, u : c**S(2)*Int((c*Csc(a + b*x))**(n + S(-2))*(A*Csc(a + b*x)**S(2) + C)*ActivateTrig(u), x))
    rubi.add(rule62)

    pattern63 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule63 = ReplacementRule(pattern63, lambda u, C, A, a, x, b, B : Int((A*Sec(a + b*x)**S(2) + B*Sec(a + b*x) + C)*ActivateTrig(u)/Sec(a + b*x)**S(2), x))
    rubi.add(rule63)

    pattern64 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule64 = ReplacementRule(pattern64, lambda u, C, A, a, x, b, B : Int((A*Csc(a + b*x)**S(2) + B*Csc(a + b*x) + C)*ActivateTrig(u)/Csc(a + b*x)**S(2), x))
    rubi.add(rule64)

    pattern65 = Pattern(Integral(u_*(A_ + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule65 = ReplacementRule(pattern65, lambda a, u, C, A, x, b : Int((A*Sec(a + b*x)**S(2) + C)*ActivateTrig(u)/Sec(a + b*x)**S(2), x))
    rubi.add(rule65)

    pattern66 = Pattern(Integral(u_*(A_ + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda u, x: KnownSecantIntegrandQ(u, x)))
    rule66 = ReplacementRule(pattern66, lambda b, C, A, a, x, u : Int((A*Csc(a + b*x)**S(2) + C)*ActivateTrig(u)/Csc(a + b*x)**S(2), x))
    rubi.add(rule66)

    pattern67 = Pattern(Integral(u_*(WC('A', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n1, n: ZeroQ(-n + n1 + S(-1))), CustomConstraint(lambda n, n2: ZeroQ(-n + n2 + S(-2))))
    rule67 = ReplacementRule(pattern67, lambda n1, n2, n, u, C, A, a, x, b, B : Int((A + B*Sec(a + b*x) + C*Sec(a + b*x)**S(2))*ActivateTrig(u)*Sec(a + b*x)**n, x))
    rubi.add(rule67)

    pattern68 = Pattern(Integral(u_*(WC('A', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n1, n: ZeroQ(-n + n1 + S(-1))), CustomConstraint(lambda n, n2: ZeroQ(-n + n2 + S(-2))))
    rule68 = ReplacementRule(pattern68, lambda n1, n2, n, u, C, A, a, x, b, B : Int((A + B*Csc(a + b*x) + C*Csc(a + b*x)**S(2))*ActivateTrig(u)*Csc(a + b*x)**n, x))
    rubi.add(rule68)

    pattern69 = Pattern(Integral(sin(x_*WC('b', S(1)) + WC('a', S(0)))*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, d: NonzeroQ(b**S(2) - d**S(2))))
    rule69 = ReplacementRule(pattern69, lambda a, d, c, x, b : -Sin(a + c + x*(b + d))/(S(2)*b + S(2)*d) + Sin(a - c + x*(b - d))/(S(2)*b - S(2)*d))
    rubi.add(rule69)

    pattern70 = Pattern(Integral(cos(x_*WC('b', S(1)) + WC('a', S(0)))*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, d: NonzeroQ(b**S(2) - d**S(2))))
    rule70 = ReplacementRule(pattern70, lambda a, d, c, x, b : Sin(a + c + x*(b + d))/(S(2)*b + S(2)*d) + Sin(a - c + x*(b - d))/(S(2)*b - S(2)*d))
    rubi.add(rule70)

    pattern71 = Pattern(Integral(sin(x_*WC('b', S(1)) + WC('a', S(0)))*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, d: NonzeroQ(b**S(2) - d**S(2))))
    rule71 = ReplacementRule(pattern71, lambda a, d, c, x, b : -Cos(a + c + x*(b + d))/(S(2)*b + S(2)*d) - Cos(a - c + x*(b - d))/(S(2)*b - S(2)*d))
    rubi.add(rule71)

    pattern72 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: PositiveIntegerQ(p/S(2))))
    rule72 = ReplacementRule(pattern72, lambda a, d, p, c, g, x, b : Int((g*Sin(c + d*x))**p, x)/S(2) + Int((g*Sin(c + d*x))**p*Cos(c + d*x), x)/S(2))
    rubi.add(rule72)

    pattern73 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: PositiveIntegerQ(p/S(2))))
    rule73 = ReplacementRule(pattern73, lambda a, d, p, c, g, x, b : Int((g*Sin(c + d*x))**p, x)/S(2) - Int((g*Sin(c + d*x))**p*Cos(c + d*x), x)/S(2))
    rubi.add(rule73)

    pattern74 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: IntegerQ(p)))
    rule74 = ReplacementRule(pattern74, lambda m, d, p, c, a, e, x, b : S(2)**p*e**(-p)*Int((e*Cos(a + b*x))**(m + p)*Sin(a + b*x)**p, x))
    rubi.add(rule74)

    pattern75 = Pattern(Integral((WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: IntegerQ(p)))
    rule75 = ReplacementRule(pattern75, lambda d, n, f, p, c, a, x, b : S(2)**p*f**(-p)*Int((f*Sin(a + b*x))**(n + p)*Cos(a + b*x)**p, x))
    rubi.add(rule75)

    pattern76 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: ZeroQ(m + p + S(-1))))
    rule76 = ReplacementRule(pattern76, lambda a, m, d, p, c, g, e, x, b : e**S(2)*(e*Cos(a + b*x))**(m + S(-2))*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))))
    rubi.add(rule76)

    pattern77 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: ZeroQ(m + p + S(-1))))
    rule77 = ReplacementRule(pattern77, lambda a, m, d, p, c, g, e, x, b : -e**S(2)*(e*Sin(a + b*x))**(m + S(-2))*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))))
    rubi.add(rule77)

    pattern78 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))))
    rule78 = ReplacementRule(pattern78, lambda m, d, p, g, a, e, c, x, b : -(e*Cos(a + b*x))**m*(g*Sin(c + d*x))**(p + S(1))/(b*g*m))
    rubi.add(rule78)

    pattern79 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))))
    rule79 = ReplacementRule(pattern79, lambda m, d, p, g, a, e, c, x, b : (e*Sin(a + b*x))**m*(g*Sin(c + d*x))**(p + S(1))/(b*g*m))
    rubi.add(rule79)

    pattern80 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda m: Greater(m, S(2))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, m: Equal(p, S(-3)/2) | Greater(m, S(3))), CustomConstraint(lambda p, m: IntegersQ(S(2)*m, S(2)*p)))
    rule80 = ReplacementRule(pattern80, lambda a, m, d, p, c, g, e, x, b : e**S(4)*(m + p + S(-1))*Int((e*Cos(a + b*x))**(m + S(-4))*(g*Sin(c + d*x))**(p + S(2)), x)/(S(4)*g**S(2)*(p + S(1))) + e**S(2)*(e*Cos(a + b*x))**(m + S(-2))*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))))
    rubi.add(rule80)

    pattern81 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda m: Greater(m, S(2))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, m: Equal(p, S(-3)/2) | Greater(m, S(3))), CustomConstraint(lambda p, m: IntegersQ(S(2)*m, S(2)*p)))
    rule81 = ReplacementRule(pattern81, lambda a, m, d, p, c, g, e, x, b : e**S(4)*(m + p + S(-1))*Int((e*Sin(a + b*x))**(m + S(-4))*(g*Sin(c + d*x))**(p + S(2)), x)/(S(4)*g**S(2)*(p + S(1))) - e**S(2)*(e*Sin(a + b*x))**(m + S(-2))*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))))
    rubi.add(rule81)

    pattern82 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p, m: Equal(m, S(2)) | Less(p, S(-2))), CustomConstraint(lambda p, m: IntegersQ(S(2)*m, S(2)*p)))
    rule82 = ReplacementRule(pattern82, lambda a, m, d, p, c, g, e, x, b : e**S(2)*(m + S(2)*p + S(2))*Int((e*Cos(a + b*x))**(m + S(-2))*(g*Sin(c + d*x))**(p + S(2)), x)/(S(4)*g**S(2)*(p + S(1))) + (e*Cos(a + b*x))**m*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))))
    rubi.add(rule82)

    pattern83 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p, m: Equal(m, S(2)) | Less(p, S(-2))), CustomConstraint(lambda p, m: IntegersQ(S(2)*m, S(2)*p)))
    rule83 = ReplacementRule(pattern83, lambda a, m, d, p, c, g, e, x, b : e**S(2)*(m + S(2)*p + S(2))*Int((e*Sin(a + b*x))**(m + S(-2))*(g*Sin(c + d*x))**(p + S(2)), x)/(S(4)*g**S(2)*(p + S(1))) - (e*Sin(a + b*x))**m*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))))
    rubi.add(rule83)

    pattern84 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p)), CustomConstraint(lambda p, m: IntegersQ(S(2)*m, S(2)*p)))
    rule84 = ReplacementRule(pattern84, lambda a, m, d, p, c, g, e, x, b : e**S(2)*(m + p + S(-1))*Int((e*Cos(a + b*x))**(m + S(-2))*(g*Sin(c + d*x))**p, x)/(m + S(2)*p) + e**S(2)*(e*Cos(a + b*x))**(m + S(-2))*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + S(2)*p)))
    rubi.add(rule84)

    pattern85 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p)), CustomConstraint(lambda p, m: IntegersQ(S(2)*m, S(2)*p)))
    rule85 = ReplacementRule(pattern85, lambda a, m, d, p, c, g, e, x, b : e**S(2)*(m + p + S(-1))*Int((e*Sin(a + b*x))**(m + S(-2))*(g*Sin(c + d*x))**p, x)/(m + S(2)*p) - e**S(2)*(e*Sin(a + b*x))**(m + S(-2))*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + S(2)*p)))
    rubi.add(rule85)

    pattern86 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p, m: NonzeroQ(m + p + S(1))), CustomConstraint(lambda p, m: IntegersQ(S(2)*m, S(2)*p)))
    rule86 = ReplacementRule(pattern86, lambda a, m, d, p, c, g, e, x, b : (m + S(2)*p + S(2))*Int((e*Cos(a + b*x))**(m + S(2))*(g*Sin(c + d*x))**p, x)/(e**S(2)*(m + p + S(1))) - (e*Cos(a + b*x))**m*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + p + S(1))))
    rubi.add(rule86)

    pattern87 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p, m: NonzeroQ(m + p + S(1))), CustomConstraint(lambda p, m: IntegersQ(S(2)*m, S(2)*p)))
    rule87 = ReplacementRule(pattern87, lambda a, m, d, p, c, g, e, x, b : (m + S(2)*p + S(2))*Int((e*Sin(a + b*x))**(m + S(2))*(g*Sin(c + d*x))**p, x)/(e**S(2)*(m + p + S(1))) + (e*Sin(a + b*x))**m*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + p + S(1))))
    rubi.add(rule87)

    pattern88 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule88 = ReplacementRule(pattern88, lambda a, d, p, c, g, x, b : S(2)*g*p*Int((g*Sin(c + d*x))**(p + S(-1))*Sin(a + b*x), x)/(S(2)*p + S(1)) + S(2)*(g*Sin(c + d*x))**p*Sin(a + b*x)/(d*(S(2)*p + S(1))))
    rubi.add(rule88)

    pattern89 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule89 = ReplacementRule(pattern89, lambda a, d, p, c, g, x, b : S(2)*g*p*Int((g*Sin(c + d*x))**(p + S(-1))*Cos(a + b*x), x)/(S(2)*p + S(1)) - S(2)*(g*Sin(c + d*x))**p*Cos(a + b*x)/(d*(S(2)*p + S(1))))
    rubi.add(rule89)

    pattern90 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule90 = ReplacementRule(pattern90, lambda a, d, p, c, g, x, b : (S(2)*p + S(3))*Int((g*Sin(c + d*x))**(p + S(1))*Sin(a + b*x), x)/(S(2)*g*(p + S(1))) + (g*Sin(c + d*x))**(p + S(1))*Cos(a + b*x)/(S(2)*b*g*(p + S(1))))
    rubi.add(rule90)

    pattern91 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule91 = ReplacementRule(pattern91, lambda a, d, p, c, g, x, b : (S(2)*p + S(3))*Int((g*Sin(c + d*x))**(p + S(1))*Cos(a + b*x), x)/(S(2)*g*(p + S(1))) - (g*Sin(c + d*x))**(p + S(1))*Sin(a + b*x)/(S(2)*b*g*(p + S(1))))
    rubi.add(rule91)

    pattern92 = Pattern(Integral(cos(x_*WC('b', S(1)) + WC('a', S(0)))/sqrt(sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)))
    rule92 = ReplacementRule(pattern92, lambda a, d, c, x, b : -ArcSin(Cos(a + b*x) - Sin(a + b*x))/d + Log(Cos(a + b*x) + Sin(a + b*x) + Sqrt(Sin(c + d*x)))/d)
    rubi.add(rule92)

    pattern93 = Pattern(Integral(sin(x_*WC('b', S(1)) + WC('a', S(0)))/sqrt(sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)))
    rule93 = ReplacementRule(pattern93, lambda a, d, c, x, b : -ArcSin(Cos(a + b*x) - Sin(a + b*x))/d - Log(Cos(a + b*x) + Sin(a + b*x) + Sqrt(Sin(c + d*x)))/d)
    rubi.add(rule93)

    pattern94 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_/cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule94 = ReplacementRule(pattern94, lambda a, d, p, c, g, x, b : S(2)*g*Int((g*Sin(c + d*x))**(p + S(-1))*Sin(a + b*x), x))
    rubi.add(rule94)

    pattern95 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_/sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule95 = ReplacementRule(pattern95, lambda a, d, p, c, g, x, b : S(2)*g*Int((g*Sin(c + d*x))**(p + S(-1))*Cos(a + b*x), x))
    rubi.add(rule95)

    pattern96 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule96 = ReplacementRule(pattern96, lambda m, d, p, g, a, e, c, x, b : (e*Cos(a + b*x))**(-p)*(g*Sin(c + d*x))**p*Int((e*Cos(a + b*x))**(m + p)*Sin(a + b*x)**p, x)*Sin(a + b*x)**(-p))
    rubi.add(rule96)

    pattern97 = Pattern(Integral((WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule97 = ReplacementRule(pattern97, lambda d, n, f, p, g, a, c, x, b : (f*Sin(a + b*x))**(-p)*(g*Sin(c + d*x))**p*Cos(a + b*x)**(-p)*Int((f*Sin(a + b*x))**(n + p)*Cos(a + b*x)**p, x))
    rubi.add(rule97)

    pattern98 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: PositiveIntegerQ(p/S(2))))
    rule98 = ReplacementRule(pattern98, lambda a, d, p, c, g, x, b : Int((g*Sin(c + d*x))**p, x)/S(4) - Int((g*Sin(c + d*x))**p*Cos(c + d*x)**S(2), x)/S(4))
    rubi.add(rule98)

    pattern99 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: IntegerQ(p)))
    rule99 = ReplacementRule(pattern99, lambda m, d, f, n, p, c, a, e, x, b : S(2)**p*e**(-p)*f**(-p)*Int((e*Cos(a + b*x))**(m + p)*(f*Sin(a + b*x))**(n + p), x))
    rubi.add(rule99)

    pattern100 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: ZeroQ(m + p + S(1))))
    rule100 = ReplacementRule(pattern100, lambda m, d, f, n, p, g, a, e, c, x, b : e*(e*Cos(a + b*x))**(m + S(-1))*(f*Sin(a + b*x))**(n + S(1))*(g*Sin(c + d*x))**p/(b*f*(n + p + S(1))))
    rubi.add(rule100)

    pattern101 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: ZeroQ(m + p + S(1))))
    rule101 = ReplacementRule(pattern101, lambda m, d, f, n, p, g, a, e, c, x, b : -e*(e*Sin(a + b*x))**(m + S(-1))*(f*Cos(a + b*x))**(n + S(1))*(g*Sin(c + d*x))**p/(b*f*(n + p + S(1))))
    rubi.add(rule101)

    pattern102 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, n, m: ZeroQ(m + n + S(2)*p + S(2))), CustomConstraint(lambda p, m: NonzeroQ(m + p + S(1))))
    rule102 = ReplacementRule(pattern102, lambda m, d, f, n, p, g, a, e, c, x, b : -(e*Cos(a + b*x))**(m + S(1))*(f*Sin(a + b*x))**(n + S(1))*(g*Sin(c + d*x))**p/(b*e*f*(m + p + S(1))))
    rubi.add(rule102)

    pattern103 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda m: Greater(m, S(3))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n: NonzeroQ(n + p + S(1))), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule103 = ReplacementRule(pattern103, lambda m, d, f, n, p, g, a, e, c, x, b : e**S(4)*(m + p + S(-1))*Int((e*Cos(a + b*x))**(m + S(-4))*(f*Sin(a + b*x))**n*(g*Sin(c + d*x))**(p + S(2)), x)/(S(4)*g**S(2)*(n + p + S(1))) + e**S(2)*(e*Cos(a + b*x))**(m + S(-2))*(f*Sin(a + b*x))**n*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))))
    rubi.add(rule103)

    pattern104 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda m: Greater(m, S(3))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n: NonzeroQ(n + p + S(1))), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule104 = ReplacementRule(pattern104, lambda m, d, f, n, p, g, a, e, c, x, b : e**S(4)*(m + p + S(-1))*Int((e*Sin(a + b*x))**(m + S(-4))*(f*Cos(a + b*x))**n*(g*Sin(c + d*x))**(p + S(2)), x)/(S(4)*g**S(2)*(n + p + S(1))) - e**S(2)*(e*Sin(a + b*x))**(m + S(-2))*(f*Cos(a + b*x))**n*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))))
    rubi.add(rule104)

    pattern105 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n + S(2)*p + S(2))), CustomConstraint(lambda p, n: NonzeroQ(n + p + S(1))), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)), CustomConstraint(lambda p, m: Equal(m, S(2)) | Equal(m, S(3)) | Less(p, S(-2))))
    rule105 = ReplacementRule(pattern105, lambda m, d, n, f, p, g, a, e, c, x, b : e**S(2)*(m + n + S(2)*p + S(2))*Int((e*Cos(a + b*x))**(m + S(-2))*(f*Sin(a + b*x))**n*(g*Sin(c + d*x))**(p + S(2)), x)/(S(4)*g**S(2)*(n + p + S(1))) + (e*Cos(a + b*x))**m*(f*Sin(a + b*x))**n*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))))
    rubi.add(rule105)

    pattern106 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n + S(2)*p + S(2))), CustomConstraint(lambda p, n: NonzeroQ(n + p + S(1))), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)), CustomConstraint(lambda p, m: Equal(m, S(2)) | Equal(m, S(3)) | Less(p, S(-2))))
    rule106 = ReplacementRule(pattern106, lambda m, d, n, f, p, g, a, e, c, x, b : e**S(2)*(m + n + S(2)*p + S(2))*Int((e*Sin(a + b*x))**(m + S(-2))*(f*Cos(a + b*x))**n*(g*Sin(c + d*x))**(p + S(2)), x)/(S(4)*g**S(2)*(n + p + S(1))) - (e*Sin(a + b*x))**m*(f*Cos(a + b*x))**n*(g*Sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))))
    rubi.add(rule106)

    pattern107 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, n: NonzeroQ(n + p + S(1))), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule107 = ReplacementRule(pattern107, lambda m, d, f, n, p, g, a, e, c, x, b : e**S(2)*(m + p + S(-1))*Int((e*Cos(a + b*x))**(m + S(-2))*(f*Sin(a + b*x))**(n + S(2))*(g*Sin(c + d*x))**p, x)/(f**S(2)*(n + p + S(1))) + e*(e*Cos(a + b*x))**(m + S(-1))*(f*Sin(a + b*x))**(n + S(1))*(g*Sin(c + d*x))**p/(b*f*(n + p + S(1))))
    rubi.add(rule107)

    pattern108 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, n: NonzeroQ(n + p + S(1))), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule108 = ReplacementRule(pattern108, lambda m, d, f, n, p, g, a, e, c, x, b : e**S(2)*(m + p + S(-1))*Int((e*Sin(a + b*x))**(m + S(-2))*(f*Cos(a + b*x))**(n + S(2))*(g*Sin(c + d*x))**p, x)/(f**S(2)*(n + p + S(1))) - e*(e*Sin(a + b*x))**(m + S(-1))*(f*Cos(a + b*x))**(n + S(1))*(g*Sin(c + d*x))**p/(b*f*(n + p + S(1))))
    rubi.add(rule108)

    pattern109 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n + S(2)*p)), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule109 = ReplacementRule(pattern109, lambda m, d, n, f, p, g, a, e, c, x, b : e**S(2)*(m + p + S(-1))*Int((e*Cos(a + b*x))**(m + S(-2))*(f*Sin(a + b*x))**n*(g*Sin(c + d*x))**p, x)/(m + n + S(2)*p) + e*(e*Cos(a + b*x))**(m + S(-1))*(f*Sin(a + b*x))**(n + S(1))*(g*Sin(c + d*x))**p/(b*f*(m + n + S(2)*p)))
    rubi.add(rule109)

    pattern110 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n + S(2)*p)), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule110 = ReplacementRule(pattern110, lambda m, d, n, f, p, g, a, e, c, x, b : e**S(2)*(m + p + S(-1))*Int((e*Sin(a + b*x))**(m + S(-2))*(f*Cos(a + b*x))**n*(g*Sin(c + d*x))**p, x)/(m + n + S(2)*p) - e*(e*Sin(a + b*x))**(m + S(-1))*(f*Cos(a + b*x))**(n + S(1))*(g*Sin(c + d*x))**p/(b*f*(m + n + S(2)*p)))
    rubi.add(rule110)

    pattern111 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, n, m: RationalQ(m, n, p)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n + S(2)*p)), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule111 = ReplacementRule(pattern111, lambda m, d, n, f, p, g, a, e, c, x, b : S(2)*f*g*(n + p + S(-1))*Int((e*Cos(a + b*x))**(m + S(1))*(f*Sin(a + b*x))**(n + S(-1))*(g*Sin(c + d*x))**(p + S(-1)), x)/(e*(m + n + S(2)*p)) - f*(e*Cos(a + b*x))**(m + S(1))*(f*Sin(a + b*x))**(n + S(-1))*(g*Sin(c + d*x))**p/(b*e*(m + n + S(2)*p)))
    rubi.add(rule111)

    pattern112 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, n, m: RationalQ(m, n, p)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n + S(2)*p)), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule112 = ReplacementRule(pattern112, lambda m, d, n, f, p, g, a, e, c, x, b : S(2)*f*g*(n + p + S(-1))*Int((e*Sin(a + b*x))**(m + S(1))*(f*Cos(a + b*x))**(n + S(-1))*(g*Sin(c + d*x))**(p + S(-1)), x)/(e*(m + n + S(2)*p)) + f*(e*Sin(a + b*x))**(m + S(1))*(f*Cos(a + b*x))**(n + S(-1))*(g*Sin(c + d*x))**p/(b*e*(m + n + S(2)*p)))
    rubi.add(rule112)

    pattern113 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, n, m: RationalQ(m, n, p)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n + S(2)*p + S(2))), CustomConstraint(lambda p, m: NonzeroQ(m + p + S(1))), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule113 = ReplacementRule(pattern113, lambda m, d, n, f, p, g, a, e, c, x, b : f*(m + n + S(2)*p + S(2))*Int((e*Cos(a + b*x))**(m + S(1))*(f*Sin(a + b*x))**(n + S(-1))*(g*Sin(c + d*x))**(p + S(1)), x)/(S(2)*e*g*(m + p + S(1))) - (e*Cos(a + b*x))**(m + S(1))*(f*Sin(a + b*x))**(n + S(1))*(g*Sin(c + d*x))**p/(b*e*f*(m + p + S(1))))
    rubi.add(rule113)

    pattern114 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda p, n, m: RationalQ(m, n, p)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n + S(2)*p + S(2))), CustomConstraint(lambda p, m: NonzeroQ(m + p + S(1))), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule114 = ReplacementRule(pattern114, lambda m, d, n, f, p, g, a, e, c, x, b : f*(m + n + S(2)*p + S(2))*Int((e*Sin(a + b*x))**(m + S(1))*(f*Cos(a + b*x))**(n + S(-1))*(g*Sin(c + d*x))**(p + S(1)), x)/(S(2)*e*g*(m + p + S(1))) + (e*Sin(a + b*x))**(m + S(1))*(f*Cos(a + b*x))**(n + S(1))*(g*Sin(c + d*x))**p/(b*e*f*(m + p + S(1))))
    rubi.add(rule114)

    pattern115 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n + S(2)*p + S(2))), CustomConstraint(lambda p, m: NonzeroQ(m + p + S(1))), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule115 = ReplacementRule(pattern115, lambda m, d, n, f, p, g, a, e, c, x, b : (m + n + S(2)*p + S(2))*Int((e*Cos(a + b*x))**(m + S(2))*(f*Sin(a + b*x))**n*(g*Sin(c + d*x))**p, x)/(e**S(2)*(m + p + S(1))) - (e*Cos(a + b*x))**(m + S(1))*(f*Sin(a + b*x))**(n + S(1))*(g*Sin(c + d*x))**p/(b*e*f*(m + p + S(1))))
    rubi.add(rule115)

    pattern116 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p, n, m: NonzeroQ(m + n + S(2)*p + S(2))), CustomConstraint(lambda p, m: NonzeroQ(m + p + S(1))), CustomConstraint(lambda p, n, m: IntegersQ(S(2)*m, S(2)*n, S(2)*p)))
    rule116 = ReplacementRule(pattern116, lambda m, d, n, f, p, g, a, e, c, x, b : (m + n + S(2)*p + S(2))*Int((e*Sin(a + b*x))**(m + S(2))*(f*Cos(a + b*x))**n*(g*Sin(c + d*x))**p, x)/(e**S(2)*(m + p + S(1))) + (e*Sin(a + b*x))**(m + S(1))*(f*Cos(a + b*x))**(n + S(1))*(g*Sin(c + d*x))**p/(b*e*f*(m + p + S(1))))
    rubi.add(rule116)

    pattern117 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: ZeroQ(S(-2) + d/b)), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule117 = ReplacementRule(pattern117, lambda m, d, f, n, p, g, a, e, c, x, b : (e*Cos(a + b*x))**(-p)*(f*Sin(a + b*x))**(-p)*(g*Sin(c + d*x))**p*Int((e*Cos(a + b*x))**(m + p)*(f*Sin(a + b*x))**(n + p), x))
    rubi.add(rule117)

    pattern118 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda m, b, d: ZeroQ(-Abs(m + S(2)) + d/b)))
    rule118 = ReplacementRule(pattern118, lambda m, d, c, a, e, x, b : (e*Cos(a + b*x))**(m + S(1))*(-m + S(-2))*Cos((a + b*x)*(m + S(1)))/(d*e*(m + S(1))))
    rubi.add(rule118)

    pattern119 = Pattern(Integral((F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_)**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda F: InertTrigQ(F)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule119 = ReplacementRule(pattern119, lambda d, n, x, p, c, a, F, b : Int(Expand((a + b*F(c + d*x)**n)**p, x), x))
    rubi.add(rule119)

    pattern120 = Pattern(Integral(1/(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda F: InertTrigQ(F)), CustomConstraint(lambda n: EvenQ(n)), CustomConstraint(lambda n: Greater(n, S(2))))
    rule120 = ReplacementRule(pattern120, lambda d, n, x, c, a, F, b : Dist(S(2)/(a*n), Sum(Int(1/(S(1) - (S(-1))**(-S(4)*k/n)*F(c + d*x)**S(2)/Rt(-a/b, n/S(2))), x), List(k, S(1), n/S(2))), x))
    rubi.add(rule120)

    pattern121 = Pattern(Integral(1/(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda F: InertTrigQ(F)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda n: Greater(n, S(2))))
    rule121 = ReplacementRule(pattern121, lambda d, n, x, c, a, F, b : Int(ExpandTrig(1/(a + b*F(c + d*x)**n), x), x))
    rubi.add(rule121)

    pattern122 = Pattern(Integral(G_**(x_*WC('d', S(1)) + WC('c', S(0)))/(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda G, F: InertTrigQ(F, G)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(2))))
    rule122 = ReplacementRule(pattern122, lambda G, m, d, n, x, c, a, F, b : Int(ExpandTrig(G(c + d*x)**m, 1/(a + b*F(c + d*x)**n), x), x))
    rubi.add(rule122)

    pattern123 = Pattern(Integral((F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('a', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda F: InertTrigQ(F)), CustomConstraint(lambda n: ~(IntegerQ(n))), CustomConstraint(lambda p: IntegerQ(p)))
    rule123 = ReplacementRule(pattern123, lambda a, d, n, p, c, x, F : With(List(Set(v, ActivateTrig(F(c + d*x)))), a**IntPart(n)*(a*v**p)**FracPart(n)*(v/NonfreeFactors(v, x))**(p*IntPart(n))*Int(NonfreeFactors(v, x)**(n*p), x)*NonfreeFactors(v, x)**(-p*FracPart(n))))
    rubi.add(rule123)

    pattern124 = Pattern(Integral(((F_*(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**p_*WC('a', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda F: InertTrigQ(F)), CustomConstraint(lambda n: ~(IntegerQ(n))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule124 = ReplacementRule(pattern124, lambda d, n, x, p, c, a, F, b : With(List(Set(v, ActivateTrig(F(c + d*x)))), a**IntPart(n)*(a*(b*v)**p)**FracPart(n)*(b*v)**(-p*FracPart(n))*Int((b*v)**(n*p), x)))
    rubi.add(rule124)

    pattern125 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda F: SameQ(F, cos) | SameQ(F, Cos)))
    rule125 = ReplacementRule(pattern125, lambda b, x, c, a, F, u : With(List(Set(d, FreeFactors(Sin(c*(a + b*x)), x))), Condition(d*Subst(Int(SubstFor(S(1), Sin(c*(a + b*x))/d, u, x), x), x, Sin(c*(a + b*x))/d)/(b*c), FunctionOfQ(Sin(c*(a + b*x))/d, u, x, True))))
    rubi.add(rule125)

    pattern126 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda F: SameQ(F, sin) | SameQ(F, Sin)))
    rule126 = ReplacementRule(pattern126, lambda b, x, c, a, F, u : With(List(Set(d, FreeFactors(Cos(c*(a + b*x)), x))), Condition(-d*Subst(Int(SubstFor(S(1), Cos(c*(a + b*x))/d, u, x), x), x, Cos(c*(a + b*x))/d)/(b*c), FunctionOfQ(Cos(c*(a + b*x))/d, u, x, True))))
    rubi.add(rule126)

    pattern127 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda F: SameQ(F, cot) | SameQ(F, Cot)))
    rule127 = ReplacementRule(pattern127, lambda b, x, c, a, F, u : With(List(Set(d, FreeFactors(Sin(c*(a + b*x)), x))), Condition(Subst(Int(SubstFor(1/x, Sin(c*(a + b*x))/d, u, x), x), x, Sin(c*(a + b*x))/d)/(b*c), FunctionOfQ(Sin(c*(a + b*x))/d, u, x, True))))
    rubi.add(rule127)

    pattern128 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda F: SameQ(F, tan) | SameQ(F, Tan)))
    rule128 = ReplacementRule(pattern128, lambda b, x, c, a, F, u : With(List(Set(d, FreeFactors(Cos(c*(a + b*x)), x))), Condition(-Subst(Int(SubstFor(1/x, Cos(c*(a + b*x))/d, u, x), x), x, Cos(c*(a + b*x))/d)/(b*c), FunctionOfQ(Cos(c*(a + b*x))/d, u, x, True))))
    rubi.add(rule128)

    pattern129 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda F: SameQ(F, sec) | SameQ(F, Sec)))
    rule129 = ReplacementRule(pattern129, lambda b, x, c, a, F, u : With(List(Set(d, FreeFactors(Tan(c*(a + b*x)), x))), Condition(d*Subst(Int(SubstFor(S(1), Tan(c*(a + b*x))/d, u, x), x), x, Tan(c*(a + b*x))/d)/(b*c), FunctionOfQ(Tan(c*(a + b*x))/d, u, x, True))))
    rubi.add(rule129)

    pattern130 = Pattern(Integral(u_/cos((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda u: NonsumQ(u)))
    rule130 = ReplacementRule(pattern130, lambda a, b, c, x, u : With(List(Set(d, FreeFactors(Tan(c*(a + b*x)), x))), Condition(d*Subst(Int(SubstFor(S(1), Tan(c*(a + b*x))/d, u, x), x), x, Tan(c*(a + b*x))/d)/(b*c), FunctionOfQ(Tan(c*(a + b*x))/d, u, x, True))))
    rubi.add(rule130)

    pattern131 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda F: SameQ(F, csc) | SameQ(F, Csc)))
    rule131 = ReplacementRule(pattern131, lambda b, x, c, a, F, u : With(List(Set(d, FreeFactors(Cot(c*(a + b*x)), x))), Condition(-d*Subst(Int(SubstFor(S(1), Cot(c*(a + b*x))/d, u, x), x), x, Cot(c*(a + b*x))/d)/(b*c), FunctionOfQ(Cot(c*(a + b*x))/d, u, x, True))))
    rubi.add(rule131)

    pattern132 = Pattern(Integral(u_/sin((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda u: NonsumQ(u)))
    rule132 = ReplacementRule(pattern132, lambda a, b, c, x, u : With(List(Set(d, FreeFactors(Cot(c*(a + b*x)), x))), Condition(-d*Subst(Int(SubstFor(S(1), Cot(c*(a + b*x))/d, u, x), x), x, Cot(c*(a + b*x))/d)/(b*c), FunctionOfQ(Cot(c*(a + b*x))/d, u, x, True))))
    rubi.add(rule132)

    pattern133 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda F: SameQ(F, cot) | SameQ(F, Cot)))
    rule133 = ReplacementRule(pattern133, lambda n, b, x, c, a, F, u : With(List(Set(d, FreeFactors(Tan(c*(a + b*x)), x))), Condition(d**(-n + S(1))*Subst(Int(SubstFor(x**(-n)/(d**S(2)*x**S(2) + S(1)), Tan(c*(a + b*x))/d, u, x), x), x, Tan(c*(a + b*x))/d)/(b*c), TryPureTanSubst(ActivateTrig(u)*Cot(c*(a + b*x))**n, x) & FunctionOfQ(Tan(c*(a + b*x))/d, u, x, True))))
    rubi.add(rule133)

    pattern134 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda F: SameQ(F, tan) | SameQ(F, Tan)))
    rule134 = ReplacementRule(pattern134, lambda n, b, x, c, a, F, u : With(List(Set(d, FreeFactors(Cot(c*(a + b*x)), x))), Condition(-d**(-n + S(1))*Subst(Int(SubstFor(x**(-n)/(d**S(2)*x**S(2) + S(1)), Cot(c*(a + b*x))/d, u, x), x), x, Cot(c*(a + b*x))/d)/(b*c), TryPureTanSubst(ActivateTrig(u)*Tan(c*(a + b*x))**n, x) & FunctionOfQ(Cot(c*(a + b*x))/d, u, x, True))))
    rubi.add(rule134)

    pattern135 = Pattern(Integral(u_, x_))
    rule135 = ReplacementRule(pattern135, lambda u, x : With(List(Set(v, FunctionOfTrig(u, x))), Condition(With(List(Set(d, FreeFactors(Cot(v), x))), Dist(-d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(1/(d**S(2)*x**S(2) + S(1)), Cot(v)/d, u, x), x), x, Cot(v)/d), x)), ~(FalseQ(v)) & TryPureTanSubst(ActivateTrig(u), x) & FunctionOfQ(NonfreeFactors(Cot(v), x), u, x, True))))
    rubi.add(rule135)

    pattern136 = Pattern(Integral(u_, x_))
    rule136 = ReplacementRule(pattern136, lambda u, x : With(List(Set(v, FunctionOfTrig(u, x))), Condition(With(List(Set(d, FreeFactors(Tan(v), x))), Dist(d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(1/(d**S(2)*x**S(2) + S(1)), Tan(v)/d, u, x), x), x, Tan(v)/d), x)), ~(FalseQ(v)) & TryPureTanSubst(ActivateTrig(u), x) & FunctionOfQ(NonfreeFactors(Tan(v), x), u, x, True))))
    rubi.add(rule136)

    pattern137 = Pattern(Integral(F_**(x_*WC('b', S(1)) + WC('a', S(0)))*G_**(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda F: SameQ(F, cos) | SameQ(F, sin)), CustomConstraint(lambda G: SameQ(G, cos) | SameQ(G, sin)), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)))
    rule137 = ReplacementRule(pattern137, lambda q, G, d, x, p, c, a, F, b : Int(ExpandTrigReduce(ActivateTrig(F(a + b*x)**p*G(c + d*x)**q), x), x))
    rubi.add(rule137)

    pattern138 = Pattern(Integral(F_**(x_*WC('b', S(1)) + WC('a', S(0)))*G_**(x_*WC('d', S(1)) + WC('c', S(0)))*H_**(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda F: SameQ(F, cos) | SameQ(F, sin)), CustomConstraint(lambda G: SameQ(G, cos) | SameQ(G, sin)), CustomConstraint(lambda H: SameQ(H, cos) | SameQ(H, sin)), CustomConstraint(lambda p, q, r: PositiveIntegerQ(p, q, r)))
    rule138 = ReplacementRule(pattern138, lambda H, q, r, G, d, f, x, p, c, a, e, F, b : Int(ExpandTrigReduce(ActivateTrig(F(a + b*x)**p*G(c + d*x)**q*H(e + f*x)**r), x), x))
    rubi.add(rule138)

    pattern139 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda F: SameQ(F, cos) | SameQ(F, Cos)))
    rule139 = ReplacementRule(pattern139, lambda b, x, c, a, F, u : With(List(Set(d, FreeFactors(Sin(c*(a + b*x)), x))), Condition(d*Subst(Int(SubstFor(S(1), Sin(c*(a + b*x))/d, u, x), x), x, Sin(c*(a + b*x))/d)/(b*c), FunctionOfQ(Sin(c*(a + b*x))/d, u, x))))
    rubi.add(rule139)

    pattern140 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda F: SameQ(F, sin) | SameQ(F, Sin)))
    rule140 = ReplacementRule(pattern140, lambda b, x, c, a, F, u : With(List(Set(d, FreeFactors(Cos(c*(a + b*x)), x))), Condition(-d*Subst(Int(SubstFor(S(1), Cos(c*(a + b*x))/d, u, x), x), x, Cos(c*(a + b*x))/d)/(b*c), FunctionOfQ(Cos(c*(a + b*x))/d, u, x))))
    rubi.add(rule140)

    pattern141 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda F: SameQ(F, cot) | SameQ(F, Cot)))
    rule141 = ReplacementRule(pattern141, lambda b, x, c, a, F, u : With(List(Set(d, FreeFactors(Sin(c*(a + b*x)), x))), Condition(Subst(Int(SubstFor(1/x, Sin(c*(a + b*x))/d, u, x), x), x, Sin(c*(a + b*x))/d)/(b*c), FunctionOfQ(Sin(c*(a + b*x))/d, u, x))))
    rubi.add(rule141)

    pattern142 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda F: SameQ(F, tan) | SameQ(F, Tan)))
    rule142 = ReplacementRule(pattern142, lambda b, x, c, a, F, u : With(List(Set(d, FreeFactors(Cos(c*(a + b*x)), x))), Condition(-Subst(Int(SubstFor(1/x, Cos(c*(a + b*x))/d, u, x), x), x, Cos(c*(a + b*x))/d)/(b*c), FunctionOfQ(Cos(c*(a + b*x))/d, u, x))))
    rubi.add(rule142)

    pattern143 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda F: SameQ(F, cos) | SameQ(F, Cos)))
    rule143 = ReplacementRule(pattern143, lambda n, b, x, c, a, F, u : With(List(Set(d, FreeFactors(Sin(c*(a + b*x)), x))), Condition(d*Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), Sin(c*(a + b*x))/d, u, x), x), x, Sin(c*(a + b*x))/d)/(b*c), FunctionOfQ(Sin(c*(a + b*x))/d, u, x))))
    rubi.add(rule143)

    pattern144 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda F: SameQ(F, sec) | SameQ(F, Sec)))
    rule144 = ReplacementRule(pattern144, lambda n, b, x, c, a, F, u : With(List(Set(d, FreeFactors(Sin(c*(a + b*x)), x))), Condition(d*Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(-n/S(2) + S(-1)/2), Sin(c*(a + b*x))/d, u, x), x), x, Sin(c*(a + b*x))/d)/(b*c), FunctionOfQ(Sin(c*(a + b*x))/d, u, x))))
    rubi.add(rule144)

    pattern145 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda F: SameQ(F, sin) | SameQ(F, Sin)))
    rule145 = ReplacementRule(pattern145, lambda n, b, x, c, a, F, u : With(List(Set(d, FreeFactors(Cos(c*(a + b*x)), x))), Condition(-d*Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), Cos(c*(a + b*x))/d, u, x), x), x, Cos(c*(a + b*x))/d)/(b*c), FunctionOfQ(Cos(c*(a + b*x))/d, u, x))))
    rubi.add(rule145)

    pattern146 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda F: SameQ(F, csc) | SameQ(F, Csc)))
    rule146 = ReplacementRule(pattern146, lambda n, b, x, c, a, F, u : With(List(Set(d, FreeFactors(Cos(c*(a + b*x)), x))), Condition(-d*Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(-n/S(2) + S(-1)/2), Cos(c*(a + b*x))/d, u, x), x), x, Cos(c*(a + b*x))/d)/(b*c), FunctionOfQ(Cos(c*(a + b*x))/d, u, x))))
    rubi.add(rule146)

    pattern147 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda F: SameQ(F, cot) | SameQ(F, Cot)))
    rule147 = ReplacementRule(pattern147, lambda n, b, x, c, a, F, u : With(List(Set(d, FreeFactors(Sin(c*(a + b*x)), x))), Condition(d**(-n + S(1))*Subst(Int(SubstFor(x**(-n)*(-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), Sin(c*(a + b*x))/d, u, x), x), x, Sin(c*(a + b*x))/d)/(b*c), FunctionOfQ(Sin(c*(a + b*x))/d, u, x))))
    rubi.add(rule147)

    pattern148 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda F: SameQ(F, tan) | SameQ(F, Tan)))
    rule148 = ReplacementRule(pattern148, lambda n, b, x, c, a, F, u : With(List(Set(d, FreeFactors(Cos(c*(a + b*x)), x))), Condition(-d**(-n + S(1))*Subst(Int(SubstFor(x**(-n)*(-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), Cos(c*(a + b*x))/d, u, x), x), x, Cos(c*(a + b*x))/d)/(b*c), FunctionOfQ(Cos(c*(a + b*x))/d, u, x))))
    rubi.add(rule148)

    pattern149 = Pattern(Integral(u_*(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*WC('d', S(1)) + v_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda v, x: NFreeQ(v, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda F: SameQ(F, cos) | SameQ(F, Cos)))
    rule149 = ReplacementRule(pattern149, lambda v, d, n, x, u, c, a, F, b : With(List(Set(e, FreeFactors(Sin(c*(a + b*x)), x))), Condition(d*Int(ActivateTrig(u)*Cos(c*(a + b*x))**n, x) + Int(ActivateTrig(u*v), x), FunctionOfQ(Sin(c*(a + b*x))/e, u, x))))
    rubi.add(rule149)

    pattern150 = Pattern(Integral(u_*(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*WC('d', S(1)) + v_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda v, x: NFreeQ(v, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda F: SameQ(F, sin) | SameQ(F, Sin)))
    rule150 = ReplacementRule(pattern150, lambda v, d, n, x, u, c, a, F, b : With(List(Set(e, FreeFactors(Cos(c*(a + b*x)), x))), Condition(d*Int(ActivateTrig(u)*Sin(c*(a + b*x))**n, x) + Int(ActivateTrig(u*v), x), FunctionOfQ(Cos(c*(a + b*x))/e, u, x))))
    rubi.add(rule150)

    pattern151 = Pattern(Integral(u_, x_))
    rule151 = ReplacementRule(pattern151, lambda u, x : With(List(Set(v, FunctionOfTrig(u, x))), Condition(With(List(Set(d, FreeFactors(Sin(v), x))), Dist(d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(S(1), Sin(v)/d, u/Cos(v), x), x), x, Sin(v)/d), x)), ~(FalseQ(v)) & FunctionOfQ(NonfreeFactors(Sin(v), x), u/Cos(v), x))))
    rubi.add(rule151)

    pattern152 = Pattern(Integral(u_, x_))
    rule152 = ReplacementRule(pattern152, lambda u, x : With(List(Set(v, FunctionOfTrig(u, x))), Condition(With(List(Set(d, FreeFactors(Cos(v), x))), Dist(-d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(S(1), Cos(v)/d, u/Sin(v), x), x), x, Cos(v)/d), x)), ~(FalseQ(v)) & FunctionOfQ(NonfreeFactors(Cos(v), x), u/Sin(v), x))))
    rubi.add(rule152)

    pattern153 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c: ZeroQ(b - c)))
    rule153 = ReplacementRule(pattern153, lambda d, u, p, c, a, e, x, b : (a + c)**p*Int(ActivateTrig(u), x))
    rubi.add(rule153)

    pattern154 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))*sec(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c: ZeroQ(b + c)))
    rule154 = ReplacementRule(pattern154, lambda d, b, p, c, a, e, x, u : (a + c)**p*Int(ActivateTrig(u), x))
    rubi.add(rule154)

    pattern155 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cot(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))*csc(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c: ZeroQ(b + c)))
    rule155 = ReplacementRule(pattern155, lambda d, p, u, c, a, e, x, b : (a + c)**p*Int(ActivateTrig(u), x))
    rubi.add(rule155)

    pattern156 = Pattern(Integral(u_/y_, x_), CustomConstraint(lambda u: ~(InertTrigFreeQ(u))))
    rule156 = ReplacementRule(pattern156, lambda u, y, x : With(List(Set(q, DerivativeDivides(ActivateTrig(y), ActivateTrig(u), x))), Condition(q*Log(RemoveContent(ActivateTrig(y), x)), ~(FalseQ(q)))))
    rubi.add(rule156)

    pattern157 = Pattern(Integral(u_/(w_*y_), x_), CustomConstraint(lambda u: ~(InertTrigFreeQ(u))))
    rule157 = ReplacementRule(pattern157, lambda x, u, y, w : With(List(Set(q, DerivativeDivides(ActivateTrig(w*y), ActivateTrig(u), x))), Condition(q*Log(RemoveContent(ActivateTrig(w*y), x)), ~(FalseQ(q)))))
    rubi.add(rule157)

    pattern158 = Pattern(Integral(u_*y_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u: ~(InertTrigFreeQ(u))))
    rule158 = ReplacementRule(pattern158, lambda u, m, y, x : With(List(Set(q, DerivativeDivides(ActivateTrig(y), ActivateTrig(u), x))), Condition(q*ActivateTrig(y**(m + S(1)))/(m + S(1)), ~(FalseQ(q)))))
    rubi.add(rule158)

    pattern159 = Pattern(Integral(u_*y_**WC('m', S(1))*z_**WC('n', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u: ~(InertTrigFreeQ(u))))
    rule159 = ReplacementRule(pattern159, lambda m, z, n, y, x, u : With(List(Set(q, DerivativeDivides(ActivateTrig(y*z), ActivateTrig(u*z**(-m + n)), x))), Condition(q*ActivateTrig(y**(m + S(1))*z**(m + S(1)))/(m + S(1)), ~(FalseQ(q)))))
    rubi.add(rule159)

    pattern160 = Pattern(Integral((F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('a', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda F: InertTrigQ(F)), CustomConstraint(lambda n: ~(IntegerQ(n))), CustomConstraint(lambda p: IntegerQ(p)))
    rule160 = ReplacementRule(pattern160, lambda d, n, x, p, c, a, F, u : With(List(Set(v, ActivateTrig(F(c + d*x)))), a**IntPart(n)*(a*v**p)**FracPart(n)*(v/NonfreeFactors(v, x))**(p*IntPart(n))*Int(ActivateTrig(u)*NonfreeFactors(v, x)**(n*p), x)*NonfreeFactors(v, x)**(-p*FracPart(n))))
    rubi.add(rule160)

    pattern161 = Pattern(Integral(((F_*(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**p_*WC('a', S(1)))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda F: InertTrigQ(F)), CustomConstraint(lambda n: ~(IntegerQ(n))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule161 = ReplacementRule(pattern161, lambda d, n, b, x, p, c, a, F, u : With(List(Set(v, ActivateTrig(F(c + d*x)))), a**IntPart(n)*(a*(b*v)**p)**FracPart(n)*(b*v)**(-p*FracPart(n))*Int((b*v)**(n*p)*ActivateTrig(u), x)))
    rubi.add(rule161)

    pattern162 = Pattern(Integral(u_, x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)))
    rule162 = ReplacementRule(pattern162, lambda u, x : With(List(Set(v, FunctionOfTrig(u, x))), Condition(With(List(Set(d, FreeFactors(Tan(v), x))), Dist(d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(1/(d**S(2)*x**S(2) + S(1)), Tan(v)/d, u, x), x), x, Tan(v)/d), x)), ~(FalseQ(v)) & FunctionOfQ(NonfreeFactors(Tan(v), x), u, x))))
    rubi.add(rule162)

    pattern163 = Pattern(Integral((WC('a', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)) + WC('b', S(1))*sec(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, n: IntegersQ(n, p)))
    rule163 = ReplacementRule(pattern163, lambda d, n, b, p, c, a, x, u : Int((a*Sin(c + d*x)**n + b)**p*ActivateTrig(u)*Sec(c + d*x)**(n*p), x))
    rubi.add(rule163)

    pattern164 = Pattern(Integral((WC('a', S(1))*cot(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)) + WC('b', S(1))*csc(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, n: IntegersQ(n, p)))
    rule164 = ReplacementRule(pattern164, lambda d, n, u, p, c, a, x, b : Int((a*Cos(c + d*x)**n + b)**p*ActivateTrig(u)*Csc(c + d*x)**(n*p), x))
    rubi.add(rule164)

    pattern165 = Pattern(Integral(u_*(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*a_ + F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda F: InertTrigQ(F)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda q, p: PosQ(-p + q)))
    rule165 = ReplacementRule(pattern165, lambda q, d, n, x, p, u, c, a, F, b : Int(ActivateTrig(u*(a + b*F(c + d*x)**(-p + q))**n*F(c + d*x)**(n*p)), x))
    rubi.add(rule165)

    pattern166 = Pattern(Integral(u_*(F_**(x_*WC('e', S(1)) + WC('d', S(0)))*a_ + F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1)) + F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda F: InertTrigQ(F)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda q, p: PosQ(-p + q)), CustomConstraint(lambda r, p: PosQ(-p + r)))
    rule166 = ReplacementRule(pattern166, lambda q, r, d, n, x, p, u, c, a, e, F, b : Int(ActivateTrig(u*(a + b*F(d + e*x)**(-p + q) + c*F(d + e*x)**(-p + r))**n*F(d + e*x)**(n*p)), x))
    rubi.add(rule166)

    pattern167 = Pattern(Integral(u_*(F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1)) + F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1)) + a_)**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda F: InertTrigQ(F)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda p: NegQ(p)))
    rule167 = ReplacementRule(pattern167, lambda q, d, n, x, p, u, c, a, e, F, b : Int(ActivateTrig(u*(a*F(d + e*x)**(-p) + b + c*F(d + e*x)**(-p + q))**n*F(d + e*x)**(n*p)), x))
    rubi.add(rule167)

    pattern168 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) + b**S(2))))
    rule168 = ReplacementRule(pattern168, lambda d, n, u, c, a, x, b : Int((a*exp(-a*(c + d*x)/b))**n*ActivateTrig(u), x))
    rubi.add(rule168)

    pattern169 = Pattern(Integral(u_, x_), CustomConstraint(lambda u: TrigSimplifyQ(u)))
    rule169 = ReplacementRule(pattern169, lambda u, x : Int(TrigSimplify(u), x))
    rubi.add(rule169)

    pattern170 = Pattern(Integral((a_*v_)**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda v: ~(InertTrigFreeQ(v))))
    rule170 = ReplacementRule(pattern170, lambda a, p, u, x, v : With(List(Set(uu, ActivateTrig(u)), Set(vv, ActivateTrig(v))), a**IntPart(p)*vv**(-FracPart(p))*(a*vv)**FracPart(p)*Int(uu*vv**p, x)))
    rubi.add(rule170)

    pattern171 = Pattern(Integral((v_**m_)**p_*WC('u', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda v: ~(InertTrigFreeQ(v))))
    rule171 = ReplacementRule(pattern171, lambda m, p, u, x, v : With(List(Set(uu, ActivateTrig(u)), Set(vv, ActivateTrig(v))), vv**(-m*FracPart(p))*(vv**m)**FracPart(p)*Int(uu*vv**(m*p), x)))
    rubi.add(rule171)

    pattern172 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda v, w: ~(InertTrigFreeQ(v)) | ~(InertTrigFreeQ(w))))
    rule172 = ReplacementRule(pattern172, lambda m, n, u, p, w, x, v : With(List(Set(uu, ActivateTrig(u)), Set(vv, ActivateTrig(v)), Set(ww, ActivateTrig(w))), vv**(-m*FracPart(p))*ww**(-n*FracPart(p))*(vv**m*ww**n)**FracPart(p)*Int(uu*vv**(m*p)*ww**(n*p), x)))
    rubi.add(rule172)

    pattern173 = Pattern(Integral(u_, x_), CustomConstraint(lambda u: ~(InertTrigFreeQ(u))))
    rule173 = ReplacementRule(pattern173, lambda u, x : With(List(Set(v, ExpandTrig(u, x))), Condition(Int(v, x), SumQ(v))))
    rubi.add(rule173)

    pattern174 = Pattern(Integral(u_, x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda u, x: ~(FalseQ(FunctionOfTrig(u, x)))))
    rule174 = ReplacementRule(pattern174, lambda u, x : With(List(Set(w, Block(List(Set(ShowSteps, False)), Int(SubstFor(1/(x**S(2)*FreeFactors(Tan(FunctionOfTrig(u, x)/S(2)), x)**S(2) + S(1)), Tan(FunctionOfTrig(u, x)/S(2))/FreeFactors(Tan(FunctionOfTrig(u, x)/S(2)), x), u, x), x)))), Condition(Module(List(Set(v, FunctionOfTrig(u, x)), d), CompoundExpression(Set(d, FreeFactors(Tan(v/S(2)), x)), Dist(S(2)*d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(1/(d**S(2)*x**S(2) + S(1)), Tan(v/S(2))/d, u, x), x), x, Tan(v/S(2))/d), x))))))
    rubi.add(rule174)

    pattern175 = Pattern(Integral(u_, x_), CustomConstraint(lambda u: ~(InertTrigFreeQ(u))))
    rule175 = ReplacementRule(pattern175, lambda u, x : With(List(Set(v, ActivateTrig(u))), Int(v, x)))
    rubi.add(rule175)

    pattern176 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cos(x_*WC('b', S(1)) + WC('a', S(0)))*Sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule176 = ReplacementRule(pattern176, lambda m, d, n, c, a, x, b : -d*m*Int((c + d*x)**(m + S(-1))*Sin(a + b*x)**(n + S(1)), x)/(b*(n + S(1))) + (c + d*x)**m*Sin(a + b*x)**(n + S(1))/(b*(n + S(1))))
    rubi.add(rule176)

    pattern177 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule177 = ReplacementRule(pattern177, lambda m, d, n, c, a, x, b : d*m*Int((c + d*x)**(m + S(-1))*Cos(a + b*x)**(n + S(1)), x)/(b*(n + S(1))) - (c + d*x)**m*Cos(a + b*x)**(n + S(1))/(b*(n + S(1))))
    rubi.add(rule177)

    pattern178 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*Sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)))
    rule178 = ReplacementRule(pattern178, lambda m, d, n, p, c, a, x, b : Int(ExpandTrigReduce((c + d*x)**m, Cos(a + b*x)**p*Sin(a + b*x)**n, x), x))
    rubi.add(rule178)

    pattern179 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)))
    rule179 = ReplacementRule(pattern179, lambda m, d, n, p, c, a, x, b : -Int((c + d*x)**m*Sin(a + b*x)**n*Tan(a + b*x)**(p + S(-2)), x) + Int((c + d*x)**m*Sin(a + b*x)**(n + S(-2))*Tan(a + b*x)**p, x))
    rubi.add(rule179)

    pattern180 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Cot(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)))
    rule180 = ReplacementRule(pattern180, lambda m, d, n, p, c, a, x, b : -Int((c + d*x)**m*Cos(a + b*x)**n*Cot(a + b*x)**(p + S(-2)), x) + Int((c + d*x)**m*Cos(a + b*x)**(n + S(-2))*Cot(a + b*x)**p, x))
    rubi.add(rule180)

    pattern181 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: SameQ(p, S(1))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule181 = ReplacementRule(pattern181, lambda m, d, n, p, c, a, x, b : -d*m*Int((c + d*x)**(m + S(-1))*Sec(a + b*x)**n, x)/(b*n) + (c + d*x)**m*Sec(a + b*x)**n/(b*n))
    rubi.add(rule181)

    pattern182 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cot(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*Csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: SameQ(p, S(1))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule182 = ReplacementRule(pattern182, lambda m, d, n, p, c, a, x, b : d*m*Int((c + d*x)**(m + S(-1))*Csc(a + b*x)**n, x)/(b*n) - (c + d*x)**m*Csc(a + b*x)**n/(b*n))
    rubi.add(rule182)

    pattern183 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Sec(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*Tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule183 = ReplacementRule(pattern183, lambda m, d, n, c, a, x, b : -d*m*Int((c + d*x)**(m + S(-1))*Tan(a + b*x)**(n + S(1)), x)/(b*(n + S(1))) + (c + d*x)**m*Tan(a + b*x)**(n + S(1))/(b*(n + S(1))))
    rubi.add(rule183)

    pattern184 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cot(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Csc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule184 = ReplacementRule(pattern184, lambda m, d, n, c, a, x, b : d*m*Int((c + d*x)**(m + S(-1))*Cot(a + b*x)**(n + S(1)), x)/(b*(n + S(1))) - (c + d*x)**m*Cot(a + b*x)**(n + S(1))/(b*(n + S(1))))
    rubi.add(rule184)

    pattern185 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Sec(x_*WC('b', S(1)) + WC('a', S(0)))*Tan(x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: PositiveIntegerQ(p/S(2))))
    rule185 = ReplacementRule(pattern185, lambda m, d, p, c, a, x, b : -Int((c + d*x)**m*Sec(a + b*x)*Tan(a + b*x)**(p + S(-2)), x) + Int((c + d*x)**m*Sec(a + b*x)**S(3)*Tan(a + b*x)**(p + S(-2)), x))
    rubi.add(rule185)

    pattern186 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Tan(x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p/S(2))))
    rule186 = ReplacementRule(pattern186, lambda m, d, n, p, c, a, x, b : -Int((c + d*x)**m*Sec(a + b*x)**n*Tan(a + b*x)**(p + S(-2)), x) + Int((c + d*x)**m*Sec(a + b*x)**(n + S(2))*Tan(a + b*x)**(p + S(-2)), x))
    rubi.add(rule186)

    pattern187 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cot(x_*WC('b', S(1)) + WC('a', S(0)))**p_*Csc(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: PositiveIntegerQ(p/S(2))))
    rule187 = ReplacementRule(pattern187, lambda m, d, p, c, a, x, b : -Int((c + d*x)**m*Cot(a + b*x)**(p + S(-2))*Csc(a + b*x), x) + Int((c + d*x)**m*Cot(a + b*x)**(p + S(-2))*Csc(a + b*x)**S(3), x))
    rubi.add(rule187)

    pattern188 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cot(x_*WC('b', S(1)) + WC('a', S(0)))**p_*Csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p/S(2))))
    rule188 = ReplacementRule(pattern188, lambda m, d, n, p, c, a, x, b : -Int((c + d*x)**m*Cot(a + b*x)**(p + S(-2))*Csc(a + b*x)**n, x) + Int((c + d*x)**m*Cot(a + b*x)**(p + S(-2))*Csc(a + b*x)**(n + S(2)), x))
    rubi.add(rule188)

    pattern189 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, n: EvenQ(n) | OddQ(p)))
    rule189 = ReplacementRule(pattern189, lambda m, d, n, p, c, a, x, b : Module(List(Set(u, IntHide(Sec(a + b*x)**n*Tan(a + b*x)**p, x))), -d*m*Int(u*(c + d*x)**(m + S(-1)), x) + Dist((c + d*x)**m, u, x)))
    rubi.add(rule189)

    pattern190 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cot(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*Csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, n: EvenQ(n) | OddQ(p)))
    rule190 = ReplacementRule(pattern190, lambda m, d, n, p, c, a, x, b : Module(List(Set(u, IntHide(Cot(a + b*x)**p*Csc(a + b*x)**n, x))), -d*m*Int(u*(c + d*x)**(m + S(-1)), x) + Dist((c + d*x)**m, u, x)))
    rubi.add(rule190)

    pattern191 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n: IntegerQ(n)))
    rule191 = ReplacementRule(pattern191, lambda m, d, n, c, a, x, b : S(2)**n*Int((c + d*x)**m*Csc(S(2)*a + S(2)*b*x)**n, x))
    rubi.add(rule191)

    pattern192 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, n: IntegersQ(n, p)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p, n: Unequal(n, p)))
    rule192 = ReplacementRule(pattern192, lambda m, d, n, p, c, a, x, b : Module(List(Set(u, IntHide(Csc(a + b*x)**n*Sec(a + b*x)**p, x))), -d*m*Int(u*(c + d*x)**(m + S(-1)), x) + Dist((c + d*x)**m, u, x)))
    rubi.add(rule192)

    pattern193 = Pattern(Integral(F_**v_*G_**w_*u_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda F: TrigQ(F)), CustomConstraint(lambda G: TrigQ(G)), CustomConstraint(lambda v, w: ZeroQ(v - w)), CustomConstraint(lambda v, x, u, w: LinearQ(List(u, v, w), x)), CustomConstraint(lambda v, x, u, w: ~(LinearMatchQ(List(u, v, w), x))))
    rule193 = ReplacementRule(pattern193, lambda v, G, m, n, p, w, x, F, u : Int(ExpandToSum(u, x)**m*F(ExpandToSum(v, x))**n*G(ExpandToSum(v, x))**p, x))
    rubi.add(rule193)

    pattern194 = Pattern(Integral((a_ + Sin(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule194 = ReplacementRule(pattern194, lambda a, m, d, f, n, c, x, e, b : -f*m*Int((a + b*Sin(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) + (a + b*Sin(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule194)

    pattern195 = Pattern(Integral((a_ + Cos(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule195 = ReplacementRule(pattern195, lambda a, m, d, n, f, c, x, e, b : f*m*Int((a + b*Cos(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) - (a + b*Cos(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule195)

    pattern196 = Pattern(Integral((a_ + Tan(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Sec(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule196 = ReplacementRule(pattern196, lambda a, m, d, f, n, c, x, e, b : -f*m*Int((a + b*Tan(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) + (a + b*Tan(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule196)

    pattern197 = Pattern(Integral((a_ + Cot(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Csc(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule197 = ReplacementRule(pattern197, lambda a, m, d, n, f, c, x, e, b : f*m*Int((a + b*Cot(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) - (a + b*Cot(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule197)

    pattern198 = Pattern(Integral((a_ + Sec(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Sec(x_*WC('d', S(1)) + WC('c', S(0)))*Tan(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule198 = ReplacementRule(pattern198, lambda a, m, d, f, n, c, x, e, b : -f*m*Int((a + b*Sec(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) + (a + b*Sec(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule198)

    pattern199 = Pattern(Integral((a_ + Csc(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cot(x_*WC('d', S(1)) + WC('c', S(0)))*Csc(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule199 = ReplacementRule(pattern199, lambda a, m, d, n, f, c, x, e, b : f*m*Int((a + b*Csc(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) - (a + b*Csc(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule199)

    pattern200 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*Sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)), CustomConstraint(lambda m: IntegerQ(m)))
    rule200 = ReplacementRule(pattern200, lambda q, m, d, f, p, c, a, e, x, b : Int(ExpandTrigReduce((e + f*x)**m, Sin(a + b*x)**p*Sin(c + d*x)**q, x), x))
    rubi.add(rule200)

    pattern201 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*Cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)), CustomConstraint(lambda m: IntegerQ(m)))
    rule201 = ReplacementRule(pattern201, lambda q, m, d, f, p, c, a, e, x, b : Int(ExpandTrigReduce((e + f*x)**m, Cos(a + b*x)**p*Cos(c + d*x)**q, x), x))
    rubi.add(rule201)

    pattern202 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*Sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)))
    rule202 = ReplacementRule(pattern202, lambda q, m, d, f, p, c, a, e, x, b : Int(ExpandTrigReduce((e + f*x)**m, Cos(c + d*x)**q*Sin(a + b*x)**p, x), x))
    rubi.add(rule202)

    pattern203 = Pattern(Integral(F_**(x_*WC('b', S(1)) + WC('a', S(0)))*G_**(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda F: MemberQ(List(Sin, Cos), F)), CustomConstraint(lambda G: MemberQ(List(Sec, Csc), G)), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)), CustomConstraint(lambda b, d, c, a: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: PositiveIntegerQ(b/d + S(-1))))
    rule203 = ReplacementRule(pattern203, lambda q, G, m, d, f, F, p, c, a, e, x, b : Int(ExpandTrigExpand((e + f*x)**m*G(c + d*x)**q, F, c + d*x, p, b/d, x), x))
    rubi.add(rule203)

    pattern204 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Sin(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, F, b, c: NonzeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2))))
    rule204 = ReplacementRule(pattern204, lambda d, x, c, a, e, F, b : F**(c*(a + b*x))*b*c*Log(F)*Sin(d + e*x)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)) - F**(c*(a + b*x))*e*Cos(d + e*x)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)))
    rubi.add(rule204)

    pattern205 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cos(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, F, b, c: NonzeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2))))
    rule205 = ReplacementRule(pattern205, lambda d, F, c, a, e, x, b : F**(c*(a + b*x))*b*c*Cos(d + e*x)*Log(F)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)) + F**(c*(a + b*x))*e*Sin(d + e*x)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)))
    rubi.add(rule205)

    pattern206 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, c, e, F, b: NonzeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule206 = ReplacementRule(pattern206, lambda d, n, x, c, a, e, F, b : F**(c*(a + b*x))*b*c*Log(F)*Sin(d + e*x)**n/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2)) - F**(c*(a + b*x))*e*n*Cos(d + e*x)*Sin(d + e*x)**(n + S(-1))/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2)) + e**S(2)*n*(n + S(-1))*Int(F**(c*(a + b*x))*Sin(d + e*x)**(n + S(-2)), x)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2)))
    rubi.add(rule206)

    pattern207 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, c, e, F, b: NonzeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*m**S(2))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(1))))
    rule207 = ReplacementRule(pattern207, lambda m, d, F, c, a, e, x, b : F**(c*(a + b*x))*b*c*Cos(d + e*x)**m*Log(F)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*m**S(2)) + F**(c*(a + b*x))*e*m*Cos(d + e*x)**(m + S(-1))*Sin(d + e*x)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*m**S(2)) + e**S(2)*m*(m + S(-1))*Int(F**(c*(a + b*x))*Cos(d + e*x)**(m + S(-2)), x)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*m**S(2)))
    rubi.add(rule207)

    pattern208 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, c, e, F, b: ZeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(2))**S(2))), CustomConstraint(lambda n: NonzeroQ(n + S(1))), CustomConstraint(lambda n: NonzeroQ(n + S(2))))
    rule208 = ReplacementRule(pattern208, lambda d, n, x, c, a, e, F, b : -F**(c*(a + b*x))*b*c*Log(F)*Sin(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))) + F**(c*(a + b*x))*Cos(d + e*x)*Sin(d + e*x)**(n + S(1))/(e*(n + S(1))))
    rubi.add(rule208)

    pattern209 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, c, e, F, b: ZeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(2))**S(2))), CustomConstraint(lambda n: NonzeroQ(n + S(1))), CustomConstraint(lambda n: NonzeroQ(n + S(2))))
    rule209 = ReplacementRule(pattern209, lambda d, n, F, c, a, e, x, b : -F**(c*(a + b*x))*b*c*Cos(d + e*x)**(n + S(2))*Log(F)/(e**S(2)*(n + S(1))*(n + S(2))) - F**(c*(a + b*x))*Cos(d + e*x)**(n + S(1))*Sin(d + e*x)/(e*(n + S(1))))
    rubi.add(rule209)

    pattern210 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, c, e, F, b: NonzeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(2))**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule210 = ReplacementRule(pattern210, lambda d, n, x, c, a, e, F, b : -F**(c*(a + b*x))*b*c*Log(F)*Sin(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))) + F**(c*(a + b*x))*Cos(d + e*x)*Sin(d + e*x)**(n + S(1))/(e*(n + S(1))) + (b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(2))**S(2))*Int(F**(c*(a + b*x))*Sin(d + e*x)**(n + S(2)), x)/(e**S(2)*(n + S(1))*(n + S(2))))
    rubi.add(rule210)

    pattern211 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, c, e, F, b: NonzeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(2))**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule211 = ReplacementRule(pattern211, lambda d, n, F, c, a, e, x, b : -F**(c*(a + b*x))*b*c*Cos(d + e*x)**(n + S(2))*Log(F)/(e**S(2)*(n + S(1))*(n + S(2))) - F**(c*(a + b*x))*Cos(d + e*x)**(n + S(1))*Sin(d + e*x)/(e*(n + S(1))) + (b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(2))**S(2))*Int(F**(c*(a + b*x))*Cos(d + e*x)**(n + S(2)), x)/(e**S(2)*(n + S(1))*(n + S(2))))
    rubi.add(rule211)

    pattern212 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(IntegerQ(n))))
    rule212 = ReplacementRule(pattern212, lambda d, n, x, c, a, e, F, b : (exp(S(2)*ImaginaryI*(d + e*x)) + S(-1))**(-n)*Int(F**(c*(a + b*x))*(exp(S(2)*ImaginaryI*(d + e*x)) + S(-1))**n*exp(-ImaginaryI*n*(d + e*x)), x)*Sin(d + e*x)**n*exp(ImaginaryI*n*(d + e*x)))
    rubi.add(rule212)

    pattern213 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(IntegerQ(n))))
    rule213 = ReplacementRule(pattern213, lambda d, n, F, c, a, e, x, b : (exp(S(2)*ImaginaryI*(d + e*x)) + S(1))**(-n)*Cos(d + e*x)**n*Int(F**(c*(a + b*x))*(exp(S(2)*ImaginaryI*(d + e*x)) + S(1))**n*exp(-ImaginaryI*n*(d + e*x)), x)*exp(ImaginaryI*n*(d + e*x)))
    rubi.add(rule213)

    pattern214 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule214 = ReplacementRule(pattern214, lambda d, n, x, c, a, e, F, b : ImaginaryI**n*Int(ExpandIntegrand(F**(c*(a + b*x))*(-exp(S(2)*ImaginaryI*(d + e*x)) + S(1))**n*(exp(S(2)*ImaginaryI*(d + e*x)) + S(1))**(-n), x), x))
    rubi.add(rule214)

    pattern215 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cot(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule215 = ReplacementRule(pattern215, lambda d, n, F, c, a, e, x, b : (-ImaginaryI)**n*Int(ExpandIntegrand(F**(c*(a + b*x))*(-exp(S(2)*ImaginaryI*(d + e*x)) + S(1))**(-n)*(exp(S(2)*ImaginaryI*(d + e*x)) + S(1))**n, x), x))
    rubi.add(rule215)

    pattern216 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Sec(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, c, e, F, b: NonzeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule216 = ReplacementRule(pattern216, lambda d, n, x, c, a, e, F, b : F**(c*(a + b*x))*b*c*Log(F)*Sec(d + e*x)**n/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2)) - F**(c*(a + b*x))*e*n*Sec(d + e*x)**(n + S(1))*Sin(d + e*x)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2)) + e**S(2)*n*(n + S(1))*Int(F**(c*(a + b*x))*Sec(d + e*x)**(n + S(2)), x)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2)))
    rubi.add(rule216)

    pattern217 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Csc(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, c, e, F, b: NonzeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule217 = ReplacementRule(pattern217, lambda d, n, F, c, a, e, x, b : F**(c*(a + b*x))*b*c*Csc(d + e*x)**n*Log(F)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2)) + F**(c*(a + b*x))*e*n*Cos(d + e*x)*Csc(d + e*x)**(n + S(1))/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2)) + e**S(2)*n*(n + S(1))*Int(F**(c*(a + b*x))*Csc(d + e*x)**(n + S(2)), x)/(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*n**S(2)))
    rubi.add(rule217)

    pattern218 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Sec(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, c, e, F, b: ZeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))), CustomConstraint(lambda n: NonzeroQ(n + S(-1))), CustomConstraint(lambda n: NonzeroQ(n + S(-2))))
    rule218 = ReplacementRule(pattern218, lambda d, n, x, c, a, e, F, b : -F**(c*(a + b*x))*b*c*Log(F)*Sec(d + e*x)**(n + S(-2))/(e**S(2)*(n + S(-2))*(n + S(-1))) + F**(c*(a + b*x))*Sec(d + e*x)**(n + S(-1))*Sin(d + e*x)/(e*(n + S(-1))))
    rubi.add(rule218)

    pattern219 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Csc(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, c, e, F, b: ZeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))), CustomConstraint(lambda n: NonzeroQ(n + S(-1))), CustomConstraint(lambda n: NonzeroQ(n + S(-2))))
    rule219 = ReplacementRule(pattern219, lambda d, n, F, c, a, e, x, b : -F**(c*(a + b*x))*b*c*Csc(d + e*x)**(n + S(-2))*Log(F)/(e**S(2)*(n + S(-2))*(n + S(-1))) + F**(c*(a + b*x))*Cos(d + e*x)*Csc(d + e*x)**(n + S(-1))/(e*(n + S(-1))))
    rubi.add(rule219)

    pattern220 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Sec(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, c, e, F, b: NonzeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda n: Unequal(n, S(2))))
    rule220 = ReplacementRule(pattern220, lambda d, n, x, c, a, e, F, b : -F**(c*(a + b*x))*b*c*Log(F)*Sec(d + e*x)**(n + S(-2))/(e**S(2)*(n + S(-2))*(n + S(-1))) + F**(c*(a + b*x))*Sec(d + e*x)**(n + S(-1))*Sin(d + e*x)/(e*(n + S(-1))) + (b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))*Int(F**(c*(a + b*x))*Sec(d + e*x)**(n + S(-2)), x)/(e**S(2)*(n + S(-2))*(n + S(-1))))
    rubi.add(rule220)

    pattern221 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Csc(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, c, e, F, b: NonzeroQ(b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda n: Unequal(n, S(2))))
    rule221 = ReplacementRule(pattern221, lambda d, n, F, c, a, e, x, b : -F**(c*(a + b*x))*b*c*Csc(d + e*x)**(n + S(-2))*Log(F)/(e**S(2)*(n + S(-2))*(n + S(-1))) - F**(c*(a + b*x))*Cos(d + e*x)*Csc(d + e*x)**(n + S(-1))/(e*(n + S(-1))) + (b**S(2)*c**S(2)*Log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))*Int(F**(c*(a + b*x))*Csc(d + e*x)**(n + S(-2)), x)/(e**S(2)*(n + S(-2))*(n + S(-1))))
    rubi.add(rule221)

    pattern222 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Sec(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule222 = ReplacementRule(pattern222, lambda d, n, x, c, a, e, F, b : S(2)**n*F**(c*(a + b*x))*Hypergeometric2F1(n, -ImaginaryI*b*c*Log(F)/(S(2)*e) + n/S(2), -ImaginaryI*b*c*Log(F)/(S(2)*e) + n/S(2) + S(1), -exp(S(2)*ImaginaryI*(d + e*x)))*exp(ImaginaryI*n*(d + e*x))/(ImaginaryI*e*n + b*c*Log(F)))
    rubi.add(rule222)

    pattern223 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Csc(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule223 = ReplacementRule(pattern223, lambda d, n, F, c, a, e, x, b : F**(c*(a + b*x))*(-S(2)*ImaginaryI)**n*Hypergeometric2F1(n, -ImaginaryI*b*c*Log(F)/(S(2)*e) + n/S(2), -ImaginaryI*b*c*Log(F)/(S(2)*e) + n/S(2) + S(1), exp(S(2)*ImaginaryI*(d + e*x)))*exp(ImaginaryI*n*(d + e*x))/(ImaginaryI*e*n + b*c*Log(F)))
    rubi.add(rule223)

    pattern224 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Sec(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: ~(IntegerQ(n))))
    rule224 = ReplacementRule(pattern224, lambda d, n, x, c, a, e, F, b : (exp(S(2)*ImaginaryI*(d + e*x)) + S(1))**n*Int(SimplifyIntegrand(F**(c*(a + b*x))*(exp(S(2)*ImaginaryI*(d + e*x)) + S(1))**(-n)*exp(ImaginaryI*n*(d + e*x)), x), x)*Sec(d + e*x)**n*exp(-ImaginaryI*n*(d + e*x)))
    rubi.add(rule224)

    pattern225 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Csc(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: ~(IntegerQ(n))))
    rule225 = ReplacementRule(pattern225, lambda d, n, F, c, a, e, x, b : (S(1) - exp(-S(2)*ImaginaryI*(d + e*x)))**n*Csc(d + e*x)**n*Int(SimplifyIntegrand(F**(c*(a + b*x))*(S(1) - exp(-S(2)*ImaginaryI*(d + e*x)))**(-n)*exp(-ImaginaryI*n*(d + e*x)), x), x)*exp(ImaginaryI*n*(d + e*x)))
    rubi.add(rule225)

    pattern226 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + Sin(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda f, g: ZeroQ(f**S(2) - g**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule226 = ReplacementRule(pattern226, lambda a, d, n, f, x, g, c, e, F, b : S(2)**n*f**n*Int(F**(c*(a + b*x))*Cos(-Pi*f/(S(4)*g) + d/S(2) + e*x/S(2))**(S(2)*n), x))
    rubi.add(rule226)

    pattern227 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + Cos(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda f, g: ZeroQ(f - g)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule227 = ReplacementRule(pattern227, lambda a, d, n, f, F, c, g, e, x, b : S(2)**n*f**n*Int(F**(c*(a + b*x))*Cos(d/S(2) + e*x/S(2))**(S(2)*n), x))
    rubi.add(rule227)

    pattern228 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + Cos(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda f, g: ZeroQ(f + g)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule228 = ReplacementRule(pattern228, lambda a, d, n, f, F, c, g, e, x, b : S(2)**n*f**n*Int(F**(c*(a + b*x))*Sin(d/S(2) + e*x/S(2))**(S(2)*n), x))
    rubi.add(rule228)

    pattern229 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + Sin(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1)))**WC('n', S(1))*Cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda f, g: ZeroQ(f**S(2) - g**S(2))), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Equal(m + n, S(0))))
    rule229 = ReplacementRule(pattern229, lambda m, d, n, f, F, c, a, e, g, x, b : g**n*Int(F**(c*(a + b*x))*Tan(Pi*f/(S(4)*g) - d/S(2) - e*x/S(2))**m, x))
    rubi.add(rule229)

    pattern230 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + Cos(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1)))**WC('n', S(1))*Sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda f, g: ZeroQ(f - g)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Equal(m + n, S(0))))
    rule230 = ReplacementRule(pattern230, lambda a, m, d, n, f, F, c, g, e, x, b : f**n*Int(F**(c*(a + b*x))*Tan(d/S(2) + e*x/S(2))**m, x))
    rubi.add(rule230)

    pattern231 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + Cos(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1)))**WC('n', S(1))*Sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda f, g: ZeroQ(f + g)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Equal(m + n, S(0))))
    rule231 = ReplacementRule(pattern231, lambda a, m, d, n, f, F, c, g, e, x, b : f**n*Int(F**(c*(a + b*x))*Cot(d/S(2) + e*x/S(2))**m, x))
    rubi.add(rule231)

    pattern232 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(h_ + Cos(x_*WC('e', S(1)) + WC('d', S(0)))*WC('i', S(1)))/(f_ + Sin(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda f, g: ZeroQ(f**S(2) - g**S(2))), CustomConstraint(lambda h, i: ZeroQ(h**S(2) - i**S(2))), CustomConstraint(lambda f, h, g, i: ZeroQ(-f*i + g*h)))
    rule232 = ReplacementRule(pattern232, lambda d, f, i, x, h, c, a, e, g, F, b : S(2)*i*Int(F**(c*(a + b*x))*Cos(d + e*x)/(f + g*Sin(d + e*x)), x) + Int(F**(c*(a + b*x))*(h - i*Cos(d + e*x))/(f + g*Sin(d + e*x)), x))
    rubi.add(rule232)

    pattern233 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(h_ + Sin(x_*WC('e', S(1)) + WC('d', S(0)))*WC('i', S(1)))/(f_ + Cos(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda f, g: ZeroQ(f**S(2) - g**S(2))), CustomConstraint(lambda h, i: ZeroQ(h**S(2) - i**S(2))), CustomConstraint(lambda f, h, g, i: ZeroQ(f*i + g*h)))
    rule233 = ReplacementRule(pattern233, lambda a, d, f, F, i, h, c, g, e, x, b : S(2)*i*Int(F**(c*(a + b*x))*Sin(d + e*x)/(f + g*Cos(d + e*x)), x) + Int(F**(c*(a + b*x))*(h - i*Sin(d + e*x))/(f + g*Cos(d + e*x)), x))
    rubi.add(rule233)

    pattern234 = Pattern(Integral(F_**(u_*WC('c', S(1)))*G_**v_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda G: TrigQ(G)), CustomConstraint(lambda v, u, x: LinearQ(List(u, v), x)), CustomConstraint(lambda v, u, x: ~(LinearMatchQ(List(u, v), x))))
    rule234 = ReplacementRule(pattern234, lambda G, v, n, c, x, F, u : Int(F**(c*ExpandToSum(u, x))*G(ExpandToSum(v, x))**n, x))
    rubi.add(rule234)

    pattern235 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*Sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule235 = ReplacementRule(pattern235, lambda m, d, n, x, c, a, e, F, b : Module(List(Set(u, IntHide(F**(c*(a + b*x))*Sin(d + e*x)**n, x))), -m*Int(u*x**(m + S(-1)), x) + Dist(x**m, u, x)))
    rubi.add(rule235)

    pattern236 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*Cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule236 = ReplacementRule(pattern236, lambda m, d, n, F, c, a, e, x, b : Module(List(Set(u, IntHide(F**(c*(a + b*x))*Cos(d + e*x)**n, x))), -m*Int(u*x**(m + S(-1)), x) + Dist(x**m, u, x)))
    rubi.add(rule236)

    pattern237 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cos(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*Sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule237 = ReplacementRule(pattern237, lambda a, m, d, n, f, F, c, g, e, x, b : Int(ExpandTrigReduce(F**(c*(a + b*x)), Cos(f + g*x)**n*Sin(d + e*x)**m, x), x))
    rubi.add(rule237)

    pattern238 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('p', S(1))*Cos(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*Sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, n, m: PositiveIntegerQ(m, n, p)))
    rule238 = ReplacementRule(pattern238, lambda a, m, d, n, f, F, p, c, g, e, x, b : Int(ExpandTrigReduce(F**(c*(a + b*x))*x**p, Cos(f + g*x)**n*Sin(d + e*x)**m, x), x))
    rubi.add(rule238)

    pattern239 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*G_**(x_*WC('e', S(1)) + WC('d', S(0)))*H_**(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)), CustomConstraint(lambda G: TrigQ(G)), CustomConstraint(lambda H: TrigQ(H)))
    rule239 = ReplacementRule(pattern239, lambda H, G, m, d, n, x, c, a, e, F, b : Int(ExpandTrigToExp(F**(c*(a + b*x)), G(d + e*x)**m*H(d + e*x)**n, x), x))
    rubi.add(rule239)

    pattern240 = Pattern(Integral(F_**u_*Sin(v_)**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda u, x: LinearQ(u, x) | PolyQ(u, x, S(2))), CustomConstraint(lambda v, x: LinearQ(v, x) | PolyQ(v, x, S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule240 = ReplacementRule(pattern240, lambda n, u, x, F, v : Int(ExpandTrigToExp(F**u, Sin(v)**n, x), x))
    rubi.add(rule240)

    pattern241 = Pattern(Integral(F_**u_*Cos(v_)**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda u, x: LinearQ(u, x) | PolyQ(u, x, S(2))), CustomConstraint(lambda v, x: LinearQ(v, x) | PolyQ(v, x, S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule241 = ReplacementRule(pattern241, lambda v, n, x, F, u : Int(ExpandTrigToExp(F**u, Cos(v)**n, x), x))
    rubi.add(rule241)

    pattern242 = Pattern(Integral(F_**u_*Cos(v_)**WC('n', S(1))*Sin(v_)**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda u, x: LinearQ(u, x) | PolyQ(u, x, S(2))), CustomConstraint(lambda v, x: LinearQ(v, x) | PolyQ(v, x, S(2))), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule242 = ReplacementRule(pattern242, lambda m, n, u, x, F, v : Int(ExpandTrigToExp(F**u, Cos(v)**n*Sin(v)**m, x), x))
    rubi.add(rule242)

    pattern243 = Pattern(Integral(Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b: ZeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule243 = ReplacementRule(pattern243, lambda a, n, p, c, x, b : x*(p + S(2))*Sin(a + b*Log(c*x**n))**(p + S(2))/(p + S(1)) + x*Cot(a + b*Log(c*x**n))*Sin(a + b*Log(c*x**n))**(p + S(2))/(b*n*(p + S(1))))
    rubi.add(rule243)

    pattern244 = Pattern(Integral(Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b: ZeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule244 = ReplacementRule(pattern244, lambda a, n, p, c, x, b : x*(p + S(2))*Cos(a + b*Log(c*x**n))**(p + S(2))/(p + S(1)) - x*Cos(a + b*Log(c*x**n))**(p + S(2))*Tan(a + b*Log(c*x**n))/(b*n*(p + S(1))))
    rubi.add(rule244)

    pattern245 = Pattern(Integral(Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda p, n, b: ZeroQ(b**S(2)*n**S(2)*p**S(2) + S(1))))
    rule245 = ReplacementRule(pattern245, lambda n, p, c, a, x, b : Int(ExpandIntegrand((-(c*x**n)**(S(1)/(n*p))*exp(-a*b*n*p)/(S(2)*b*n*p) + (c*x**n)**(-S(1)/(n*p))*exp(a*b*n*p)/(S(2)*b*n*p))**p, x), x))
    rubi.add(rule245)

    pattern246 = Pattern(Integral(Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda p, n, b: ZeroQ(b**S(2)*n**S(2)*p**S(2) + S(1))))
    rule246 = ReplacementRule(pattern246, lambda n, p, c, a, x, b : Int(ExpandIntegrand((-(c*x**n)**(S(1)/(n*p))*exp(-a*b*n*p)/S(2) + (c*x**n)**(-S(1)/(n*p))*exp(a*b*n*p)/S(2))**p, x), x))
    rubi.add(rule246)

    pattern247 = Pattern(Integral(Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, b: NonzeroQ(b**S(2)*n**S(2) + S(1))))
    rule247 = ReplacementRule(pattern247, lambda a, n, c, x, b : -b*n*x*Cos(a + b*Log(c*x**n))/(b**S(2)*n**S(2) + S(1)) + x*Sin(a + b*Log(c*x**n))/(b**S(2)*n**S(2) + S(1)))
    rubi.add(rule247)

    pattern248 = Pattern(Integral(Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, b: NonzeroQ(b**S(2)*n**S(2) + S(1))))
    rule248 = ReplacementRule(pattern248, lambda a, n, c, x, b : b*n*x*Sin(a + b*Log(c*x**n))/(b**S(2)*n**S(2) + S(1)) + x*Cos(a + b*Log(c*x**n))/(b**S(2)*n**S(2) + S(1)))
    rubi.add(rule248)

    pattern249 = Pattern(Integral(Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(1))))
    rule249 = ReplacementRule(pattern249, lambda a, n, p, c, x, b : b**S(2)*n**S(2)*p*(p + S(-1))*Int(Sin(a + b*Log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*p**S(2) + S(1)) - b*n*p*x*Cos(a + b*Log(c*x**n))*Sin(a + b*Log(c*x**n))**(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + S(1)) + x*Sin(a + b*Log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)))
    rubi.add(rule249)

    pattern250 = Pattern(Integral(Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(1))))
    rule250 = ReplacementRule(pattern250, lambda a, n, p, c, x, b : b**S(2)*n**S(2)*p*(p + S(-1))*Int(Cos(a + b*Log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*p**S(2) + S(1)) + b*n*p*x*Cos(a + b*Log(c*x**n))**(p + S(-1))*Sin(a + b*Log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + S(1)) + x*Cos(a + b*Log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)))
    rubi.add(rule250)

    pattern251 = Pattern(Integral(Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + S(1))))
    rule251 = ReplacementRule(pattern251, lambda a, n, p, c, x, b : x*Cot(a + b*Log(c*x**n))*Sin(a + b*Log(c*x**n))**(p + S(2))/(b*n*(p + S(1))) - x*Sin(a + b*Log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))) + (b**S(2)*n**S(2)*(p + S(2))**S(2) + S(1))*Int(Sin(a + b*Log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))))
    rubi.add(rule251)

    pattern252 = Pattern(Integral(Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + S(1))))
    rule252 = ReplacementRule(pattern252, lambda a, n, p, c, x, b : -x*Cos(a + b*Log(c*x**n))**(p + S(2))*Tan(a + b*Log(c*x**n))/(b*n*(p + S(1))) - x*Cos(a + b*Log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))) + (b**S(2)*n**S(2)*(p + S(2))**S(2) + S(1))*Int(Cos(a + b*Log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))))
    rubi.add(rule252)

    pattern253 = Pattern(Integral(Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(1))))
    rule253 = ReplacementRule(pattern253, lambda a, n, p, c, x, b : x*(-S(2)*(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(2))**(-p)*(-ImaginaryI*(c*x**n)**(ImaginaryI*b)*exp(ImaginaryI*a) + ImaginaryI*(c*x**n)**(-ImaginaryI*b)*exp(-ImaginaryI*a))**p*Hypergeometric2F1(-p, (-ImaginaryI*b*n*p + S(1))/(S(2)*ImaginaryI*b*n), S(1) + (-ImaginaryI*b*n*p + S(1))/(S(2)*ImaginaryI*b*n), (c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a))/(-ImaginaryI*b*n*p + S(1)))
    rubi.add(rule253)

    pattern254 = Pattern(Integral(Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(1))))
    rule254 = ReplacementRule(pattern254, lambda a, n, p, c, x, b : x*((c*x**n)**(ImaginaryI*b)*exp(ImaginaryI*a) + (c*x**n)**(-ImaginaryI*b)*exp(-ImaginaryI*a))**p*(S(2)*(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(2))**(-p)*Hypergeometric2F1(-p, (-ImaginaryI*b*n*p + S(1))/(S(2)*ImaginaryI*b*n), S(1) + (-ImaginaryI*b*n*p + S(1))/(S(2)*ImaginaryI*b*n), -(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a))/(-ImaginaryI*b*n*p + S(1)))
    rubi.add(rule254)

    pattern255 = Pattern(Integral(x_**WC('m', S(1))*Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b, m: ZeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + (m + S(1))**S(2))), CustomConstraint(lambda p: NonzeroQ(p + S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule255 = ReplacementRule(pattern255, lambda m, n, p, c, a, x, b : x**(m + S(1))*(p + S(2))*Sin(a + b*Log(c*x**n))**(p + S(2))/((m + S(1))*(p + S(1))) + x**(m + S(1))*Cot(a + b*Log(c*x**n))*Sin(a + b*Log(c*x**n))**(p + S(2))/(b*n*(p + S(1))))
    rubi.add(rule255)

    pattern256 = Pattern(Integral(x_**WC('m', S(1))*Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b, m: ZeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + (m + S(1))**S(2))), CustomConstraint(lambda p: NonzeroQ(p + S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule256 = ReplacementRule(pattern256, lambda m, n, p, c, a, x, b : x**(m + S(1))*(p + S(2))*Cos(a + b*Log(c*x**n))**(p + S(2))/((m + S(1))*(p + S(1))) - x**(m + S(1))*Cos(a + b*Log(c*x**n))**(p + S(2))*Tan(a + b*Log(c*x**n))/(b*n*(p + S(1))))
    rubi.add(rule256)

    pattern257 = Pattern(Integral(x_**WC('m', S(1))*Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda p, n, b, m: ZeroQ(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2))))
    rule257 = ReplacementRule(pattern257, lambda m, n, p, c, a, x, b : S(2)**(-p)*Int(ExpandIntegrand(x**m*((c*x**n)**((-m + S(-1))/(n*p))*(m + S(1))*exp(a*b*n*p/(m + S(1)))/(b*n*p) - (c*x**n)**((m + S(1))/(n*p))*(m + S(1))*exp(-a*b*n*p/(m + S(1)))/(b*n*p))**p, x), x))
    rubi.add(rule257)

    pattern258 = Pattern(Integral(x_**WC('m', S(1))*Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda p, n, b, m: ZeroQ(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2))))
    rule258 = ReplacementRule(pattern258, lambda m, n, p, c, a, x, b : S(2)**(-p)*Int(ExpandIntegrand(x**m*((c*x**n)**((-m + S(-1))/(n*p))*exp(a*b*n*p/(m + S(1))) - (c*x**n)**((m + S(1))/(n*p))*exp(-a*b*n*p/(m + S(1))))**p, x), x))
    rubi.add(rule258)

    pattern259 = Pattern(Integral(x_**WC('m', S(1))*Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, n, b: NonzeroQ(b**S(2)*n**S(2) + (m + S(1))**S(2))))
    rule259 = ReplacementRule(pattern259, lambda m, n, c, a, x, b : -b*n*x**(m + S(1))*Cos(a + b*Log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)) + x**(m + S(1))*(m + S(1))*Sin(a + b*Log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)))
    rubi.add(rule259)

    pattern260 = Pattern(Integral(x_**WC('m', S(1))*Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, n, b: NonzeroQ(b**S(2)*n**S(2) + (m + S(1))**S(2))))
    rule260 = ReplacementRule(pattern260, lambda m, n, c, a, x, b : b*n*x**(m + S(1))*Sin(a + b*Log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)) + x**(m + S(1))*(m + S(1))*Cos(a + b*Log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)))
    rubi.add(rule260)

    pattern261 = Pattern(Integral(x_**WC('m', S(1))*Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2))))
    rule261 = ReplacementRule(pattern261, lambda m, n, p, c, a, x, b : b**S(2)*n**S(2)*p*(p + S(-1))*Int(x**m*Sin(a + b*Log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)) - b*n*p*x**(m + S(1))*Cos(a + b*Log(c*x**n))*Sin(a + b*Log(c*x**n))**(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)) + x**(m + S(1))*(m + S(1))*Sin(a + b*Log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)))
    rubi.add(rule261)

    pattern262 = Pattern(Integral(x_**WC('m', S(1))*Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2))))
    rule262 = ReplacementRule(pattern262, lambda m, n, p, c, a, x, b : b**S(2)*n**S(2)*p*(p + S(-1))*Int(x**m*Cos(a + b*Log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)) + b*n*p*x**(m + S(1))*Cos(a + b*Log(c*x**n))**(p + S(-1))*Sin(a + b*Log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)) + x**(m + S(1))*(m + S(1))*Cos(a + b*Log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)))
    rubi.add(rule262)

    pattern263 = Pattern(Integral(x_**WC('m', S(1))*Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + (m + S(1))**S(2))))
    rule263 = ReplacementRule(pattern263, lambda m, n, p, c, a, x, b : x**(m + S(1))*Cot(a + b*Log(c*x**n))*Sin(a + b*Log(c*x**n))**(p + S(2))/(b*n*(p + S(1))) - x**(m + S(1))*(m + S(1))*Sin(a + b*Log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))) + (b**S(2)*n**S(2)*(p + S(2))**S(2) + (m + S(1))**S(2))*Int(x**m*Sin(a + b*Log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))))
    rubi.add(rule263)

    pattern264 = Pattern(Integral(x_**WC('m', S(1))*Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + (m + S(1))**S(2))))
    rule264 = ReplacementRule(pattern264, lambda m, n, p, c, a, x, b : -x**(m + S(1))*Cos(a + b*Log(c*x**n))**(p + S(2))*Tan(a + b*Log(c*x**n))/(b*n*(p + S(1))) - x**(m + S(1))*(m + S(1))*Cos(a + b*Log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))) + (b**S(2)*n**S(2)*(p + S(2))**S(2) + (m + S(1))**S(2))*Int(x**m*Cos(a + b*Log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))))
    rubi.add(rule264)

    pattern265 = Pattern(Integral(x_**WC('m', S(1))*Sin(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2))))
    rule265 = ReplacementRule(pattern265, lambda m, n, p, c, a, x, b : x**(m + S(1))*(-S(2)*(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(2))**(-p)*(-ImaginaryI*(c*x**n)**(ImaginaryI*b)*exp(ImaginaryI*a) + ImaginaryI*(c*x**n)**(-ImaginaryI*b)*exp(-ImaginaryI*a))**p*Hypergeometric2F1(-p, (-ImaginaryI*b*n*p + m + S(1))/(S(2)*ImaginaryI*b*n), S(1) + (-ImaginaryI*b*n*p + m + S(1))/(S(2)*ImaginaryI*b*n), (c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a))/(-ImaginaryI*b*n*p + m + S(1)))
    rubi.add(rule265)

    pattern266 = Pattern(Integral(x_**WC('m', S(1))*Cos(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2))))
    rule266 = ReplacementRule(pattern266, lambda m, n, p, c, a, x, b : x**(m + S(1))*((c*x**n)**(ImaginaryI*b)*exp(ImaginaryI*a) + (c*x**n)**(-ImaginaryI*b)*exp(-ImaginaryI*a))**p*(S(2)*(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(2))**(-p)*Hypergeometric2F1(-p, (-ImaginaryI*b*n*p + m + S(1))/(S(2)*ImaginaryI*b*n), S(1) + (-ImaginaryI*b*n*p + m + S(1))/(S(2)*ImaginaryI*b*n), -(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a))/(-ImaginaryI*b*n*p + m + S(1)))
    rubi.add(rule266)

    pattern267 = Pattern(Integral(Sec(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, b: ZeroQ(b**S(2)*n**S(2) + S(1))))
    rule267 = ReplacementRule(pattern267, lambda a, n, c, x, b : S(2)*Int((c*x**n)**(1/n)/((c*x**n)**(S(2)/n) + exp(S(2)*a*b*n)), x)*exp(a*b*n))
    rubi.add(rule267)

    pattern268 = Pattern(Integral(Csc(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, b: ZeroQ(b**S(2)*n**S(2) + S(1))))
    rule268 = ReplacementRule(pattern268, lambda a, n, c, x, b : S(2)*b*n*Int((c*x**n)**(1/n)/(-(c*x**n)**(S(2)/n) + exp(S(2)*a*b*n)), x)*exp(a*b*n))
    rubi.add(rule268)

    pattern269 = Pattern(Integral(Sec(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b: ZeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(-1))))
    rule269 = ReplacementRule(pattern269, lambda a, n, p, c, x, b : x*(p + S(-2))*Sec(a + b*Log(c*x**n))**(p + S(-2))/(p + S(-1)) + x*Sec(a + b*Log(c*x**n))**(p + S(-2))*Tan(a + b*Log(c*x**n))/(b*n*(p + S(-1))))
    rubi.add(rule269)

    pattern270 = Pattern(Integral(Csc(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b: ZeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(-1))))
    rule270 = ReplacementRule(pattern270, lambda a, n, p, c, x, b : x*(p + S(-2))*Csc(a + b*Log(c*x**n))**(p + S(-2))/(p + S(-1)) - x*Cot(a + b*Log(c*x**n))*Csc(a + b*Log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))))
    rubi.add(rule270)

    pattern271 = Pattern(Integral(Sec(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p: Unequal(p, S(2))), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(1))))
    rule271 = ReplacementRule(pattern271, lambda a, n, p, c, x, b : x*Sec(a + b*Log(c*x**n))**(p + S(-2))*Tan(a + b*Log(c*x**n))/(b*n*(p + S(-1))) - x*Sec(a + b*Log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))) + (b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(1))*Int(Sec(a + b*Log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))))
    rubi.add(rule271)

    pattern272 = Pattern(Integral(Csc(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p: Unequal(p, S(2))), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(1))))
    rule272 = ReplacementRule(pattern272, lambda a, n, p, c, x, b : -x*Cot(a + b*Log(c*x**n))*Csc(a + b*Log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))) - x*Csc(a + b*Log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))) + (b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(1))*Int(Csc(a + b*Log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))))
    rubi.add(rule272)

    pattern273 = Pattern(Integral(Sec(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(1))))
    rule273 = ReplacementRule(pattern273, lambda a, n, p, c, x, b : b**S(2)*n**S(2)*p*(p + S(1))*Int(Sec(a + b*Log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*p**S(2) + S(1)) - b*n*p*x*Sec(a + b*Log(c*x**n))**(p + S(1))*Sin(a + b*Log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + S(1)) + x*Sec(a + b*Log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)))
    rubi.add(rule273)

    pattern274 = Pattern(Integral(Csc(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(1))))
    rule274 = ReplacementRule(pattern274, lambda a, n, p, c, x, b : b**S(2)*n**S(2)*p*(p + S(1))*Int(Csc(a + b*Log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*p**S(2) + S(1)) + b*n*p*x*Cos(a + b*Log(c*x**n))*Csc(a + b*Log(c*x**n))**(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + S(1)) + x*Csc(a + b*Log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)))
    rubi.add(rule274)

    pattern275 = Pattern(Integral(Sec(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(1))))
    rule275 = ReplacementRule(pattern275, lambda n, p, c, a, x, b : x*((c*x**n)**(ImaginaryI*b)*exp(ImaginaryI*a)/((c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(1)))**p*(S(2)*(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(2))**p*Hypergeometric2F1(p, (ImaginaryI*b*n*p + S(1))/(S(2)*ImaginaryI*b*n), S(1) + (ImaginaryI*b*n*p + S(1))/(S(2)*ImaginaryI*b*n), -(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a))/(ImaginaryI*b*n*p + S(1)))
    rubi.add(rule275)

    pattern276 = Pattern(Integral(Csc(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(1))))
    rule276 = ReplacementRule(pattern276, lambda n, p, c, a, x, b : x*(-ImaginaryI*(c*x**n)**(ImaginaryI*b)*exp(ImaginaryI*a)/(-(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(1)))**p*(-S(2)*(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(2))**p*Hypergeometric2F1(p, (ImaginaryI*b*n*p + S(1))/(S(2)*ImaginaryI*b*n), S(1) + (ImaginaryI*b*n*p + S(1))/(S(2)*ImaginaryI*b*n), (c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a))/(ImaginaryI*b*n*p + S(1)))
    rubi.add(rule276)

    pattern277 = Pattern(Integral(x_**WC('m', S(1))*Sec(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, n, b: ZeroQ(b**S(2)*n**S(2) + (m + S(1))**S(2))))
    rule277 = ReplacementRule(pattern277, lambda m, n, c, a, x, b : S(2)*Int(x**m*(c*x**n)**((m + S(1))/n)/((c*x**n)**(S(2)*(m + S(1))/n) + exp(S(2)*a*b*n/(m + S(1)))), x)*exp(a*b*n/(m + S(1))))
    rubi.add(rule277)

    pattern278 = Pattern(Integral(x_**WC('m', S(1))*Csc(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, n, b: ZeroQ(b**S(2)*n**S(2) + (m + S(1))**S(2))))
    rule278 = ReplacementRule(pattern278, lambda m, n, c, a, x, b : S(2)*b*n*Int(x**m*(c*x**n)**((m + S(1))/n)/(-(c*x**n)**(S(2)*(m + S(1))/n) + exp(S(2)*a*b*n/(m + S(1)))), x)*exp(a*b*n/(m + S(1)))/(m + S(1)))
    rubi.add(rule278)

    pattern279 = Pattern(Integral(x_**WC('m', S(1))*Sec(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b, m: ZeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(-1))))
    rule279 = ReplacementRule(pattern279, lambda m, n, p, c, a, x, b : x**(m + S(1))*(p + S(-2))*Sec(a + b*Log(c*x**n))**(p + S(-2))/((m + S(1))*(p + S(-1))) + x**(m + S(1))*Sec(a + b*Log(c*x**n))**(p + S(-2))*Tan(a + b*Log(c*x**n))/(b*n*(p + S(-1))))
    rubi.add(rule279)

    pattern280 = Pattern(Integral(x_**WC('m', S(1))*Csc(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b, m: ZeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(-1))))
    rule280 = ReplacementRule(pattern280, lambda m, n, p, c, a, x, b : x**(m + S(1))*(p + S(-2))*Csc(a + b*Log(c*x**n))**(p + S(-2))/((m + S(1))*(p + S(-1))) - x**(m + S(1))*Cot(a + b*Log(c*x**n))*Csc(a + b*Log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))))
    rubi.add(rule280)

    pattern281 = Pattern(Integral(x_**WC('m', S(1))*Sec(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p: Unequal(p, S(2))), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))))
    rule281 = ReplacementRule(pattern281, lambda m, n, p, c, a, x, b : x**(m + S(1))*Sec(a + b*Log(c*x**n))**(p + S(-2))*Tan(a + b*Log(c*x**n))/(b*n*(p + S(-1))) - x**(m + S(1))*(m + S(1))*Sec(a + b*Log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))) + (b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))*Int(x**m*Sec(a + b*Log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))))
    rubi.add(rule281)

    pattern282 = Pattern(Integral(x_**WC('m', S(1))*Csc(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p: Unequal(p, S(2))), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))))
    rule282 = ReplacementRule(pattern282, lambda m, n, p, c, a, x, b : -x**(m + S(1))*Cot(a + b*Log(c*x**n))*Csc(a + b*Log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))) - x**(m + S(1))*(m + S(1))*Csc(a + b*Log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))) + (b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))*Int(x**m*Csc(a + b*Log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))))
    rubi.add(rule282)

    pattern283 = Pattern(Integral(x_**WC('m', S(1))*Sec(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2))))
    rule283 = ReplacementRule(pattern283, lambda m, n, p, c, a, x, b : b**S(2)*n**S(2)*p*(p + S(1))*Int(x**m*Sec(a + b*Log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)) - b*n*p*x**(m + S(1))*Sec(a + b*Log(c*x**n))**(p + S(1))*Sin(a + b*Log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)) + x**(m + S(1))*(m + S(1))*Sec(a + b*Log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)))
    rubi.add(rule283)

    pattern284 = Pattern(Integral(x_**WC('m', S(1))*Csc(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2))))
    rule284 = ReplacementRule(pattern284, lambda m, n, p, c, a, x, b : b**S(2)*n**S(2)*p*(p + S(1))*Int(x**m*Csc(a + b*Log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)) + b*n*p*x**(m + S(1))*Cos(a + b*Log(c*x**n))*Csc(a + b*Log(c*x**n))**(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)) + x**(m + S(1))*(m + S(1))*Csc(a + b*Log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)))
    rubi.add(rule284)

    pattern285 = Pattern(Integral(x_**WC('m', S(1))*Sec(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2))))
    rule285 = ReplacementRule(pattern285, lambda m, n, p, c, a, x, b : x**(m + S(1))*((c*x**n)**(ImaginaryI*b)*exp(ImaginaryI*a)/((c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(1)))**p*(S(2)*(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(2))**p*Hypergeometric2F1(p, (ImaginaryI*b*n*p + m + S(1))/(S(2)*ImaginaryI*b*n), S(1) + (ImaginaryI*b*n*p + m + S(1))/(S(2)*ImaginaryI*b*n), -(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a))/(ImaginaryI*b*n*p + m + S(1)))
    rubi.add(rule285)

    pattern286 = Pattern(Integral(x_**WC('m', S(1))*Csc(Log(x_**WC('n', S(1))*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, n, b, m: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2))))
    rule286 = ReplacementRule(pattern286, lambda m, n, p, c, a, x, b : x**(m + S(1))*(-ImaginaryI*(c*x**n)**(ImaginaryI*b)*exp(ImaginaryI*a)/(-(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(1)))**p*(-S(2)*(c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a) + S(2))**p*Hypergeometric2F1(p, (ImaginaryI*b*n*p + m + S(1))/(S(2)*ImaginaryI*b*n), S(1) + (ImaginaryI*b*n*p + m + S(1))/(S(2)*ImaginaryI*b*n), (c*x**n)**(S(2)*ImaginaryI*b)*exp(S(2)*ImaginaryI*a))/(ImaginaryI*b*n*p + m + S(1)))
    rubi.add(rule286)

    pattern287 = Pattern(Integral(Log(x_*WC('b', S(1)))**WC('p', S(1))*Sin(x_*Log(x_*WC('b', S(1)))**WC('p', S(1))*WC('a', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule287 = ReplacementRule(pattern287, lambda p, x, b, a : -p*Int(Log(b*x)**(p + S(-1))*Sin(a*x*Log(b*x)**p), x) - Cos(a*x*Log(b*x)**p)/a)
    rubi.add(rule287)

    pattern288 = Pattern(Integral(Cos(x_*Log(x_*WC('b', S(1)))**WC('p', S(1))*WC('a', S(1)))*Log(x_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule288 = ReplacementRule(pattern288, lambda p, x, b, a : -p*Int(Cos(a*x*Log(b*x)**p)*Log(b*x)**(p + S(-1)), x) + Sin(a*x*Log(b*x)**p)/a)
    rubi.add(rule288)

    pattern289 = Pattern(Integral(Log(x_*WC('b', S(1)))**WC('p', S(1))*Sin(x_**n_*Log(x_*WC('b', S(1)))**WC('p', S(1))*WC('a', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule289 = ReplacementRule(pattern289, lambda a, n, p, x, b : -p*Int(Log(b*x)**(p + S(-1))*Sin(a*x**n*Log(b*x)**p), x)/n - x**(-n + S(1))*Cos(a*x**n*Log(b*x)**p)/(a*n) - (n + S(-1))*Int(x**(-n)*Cos(a*x**n*Log(b*x)**p), x)/(a*n))
    rubi.add(rule289)

    pattern290 = Pattern(Integral(Cos(x_**n_*Log(x_*WC('b', S(1)))**WC('p', S(1))*WC('a', S(1)))*Log(x_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule290 = ReplacementRule(pattern290, lambda a, n, p, x, b : -p*Int(Cos(a*x**n*Log(b*x)**p)*Log(b*x)**(p + S(-1)), x)/n + x**(-n + S(1))*Sin(a*x**n*Log(b*x)**p)/(a*n) + (n + S(-1))*Int(x**(-n)*Sin(a*x**n*Log(b*x)**p), x)/(a*n))
    rubi.add(rule290)

    pattern291 = Pattern(Integral(x_**WC('m', S(1))*Log(x_*WC('b', S(1)))**WC('p', S(1))*Sin(x_**WC('n', S(1))*Log(x_*WC('b', S(1)))**WC('p', S(1))*WC('a', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule291 = ReplacementRule(pattern291, lambda m, n, p, a, x, b : -p*Int(x**m*Log(b*x)**(p + S(-1))*Sin(a*x**n*Log(b*x)**p), x)/n - Cos(a*x**n*Log(b*x)**p)/(a*n))
    rubi.add(rule291)

    pattern292 = Pattern(Integral(x_**WC('m', S(1))*Cos(x_**WC('n', S(1))*Log(x_*WC('b', S(1)))**WC('p', S(1))*WC('a', S(1)))*Log(x_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule292 = ReplacementRule(pattern292, lambda m, n, p, a, x, b : -p*Int(x**m*Cos(a*x**n*Log(b*x)**p)*Log(b*x)**(p + S(-1)), x)/n + Sin(a*x**n*Log(b*x)**p)/(a*n))
    rubi.add(rule292)

    pattern293 = Pattern(Integral(x_**WC('m', S(1))*Log(x_*WC('b', S(1)))**WC('p', S(1))*Sin(x_**WC('n', S(1))*Log(x_*WC('b', S(1)))**WC('p', S(1))*WC('a', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n, m: NonzeroQ(m - n + S(1))))
    rule293 = ReplacementRule(pattern293, lambda m, n, p, a, x, b : -p*Int(x**m*Log(b*x)**(p + S(-1))*Sin(a*x**n*Log(b*x)**p), x)/n - x**(m - n + S(1))*Cos(a*x**n*Log(b*x)**p)/(a*n) + (m - n + S(1))*Int(x**(m - n)*Cos(a*x**n*Log(b*x)**p), x)/(a*n))
    rubi.add(rule293)

    pattern294 = Pattern(Integral(x_**m_*Cos(x_**WC('n', S(1))*Log(x_*WC('b', S(1)))**WC('p', S(1))*WC('a', S(1)))*Log(x_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n, m: NonzeroQ(m - n + S(1))))
    rule294 = ReplacementRule(pattern294, lambda a, m, n, p, x, b : -p*Int(x**m*Cos(a*x**n*Log(b*x)**p)*Log(b*x)**(p + S(-1)), x)/n + x**(m - n + S(1))*Sin(a*x**n*Log(b*x)**p)/(a*n) - (m - n + S(1))*Int(x**(m - n)*Sin(a*x**n*Log(b*x)**p), x)/(a*n))
    rubi.add(rule294)

    pattern295 = Pattern(Integral(Sin(WC('a', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule295 = ReplacementRule(pattern295, lambda a, d, n, c, x : -Subst(Int(Sin(a*x)**n/x**S(2), x), x, 1/(c + d*x))/d)
    rubi.add(rule295)

    pattern296 = Pattern(Integral(Cos(WC('a', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule296 = ReplacementRule(pattern296, lambda a, d, n, c, x : -Subst(Int(Cos(a*x)**n/x**S(2), x), x, 1/(c + d*x))/d)
    rubi.add(rule296)

    pattern297 = Pattern(Integral(Sin((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda b, d, c, a: NonzeroQ(-a*d + b*c)))
    rule297 = ReplacementRule(pattern297, lambda d, n, c, a, e, x, b : -Subst(Int(Sin(b*e/d - e*x*(-a*d + b*c)/d)**n/x**S(2), x), x, 1/(c + d*x))/d)
    rubi.add(rule297)

    pattern298 = Pattern(Integral(Cos((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda b, d, c, a: NonzeroQ(-a*d + b*c)))
    rule298 = ReplacementRule(pattern298, lambda d, n, c, a, e, x, b : -Subst(Int(Cos(b*e/d - e*x*(-a*d + b*c)/d)**n/x**S(2), x), x, 1/(c + d*x))/d)
    rubi.add(rule298)

    pattern299 = Pattern(Integral(Sin(u_)**WC('n', S(1)), x_), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda u, x: QuotientOfLinearsQ(u, x)))
    rule299 = ReplacementRule(pattern299, lambda n, u, x : Module(List(Set(lst, QuotientOfLinearsParts(u, x))), Int(Sin((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**n, x)))
    rubi.add(rule299)

    pattern300 = Pattern(Integral(Cos(u_)**WC('n', S(1)), x_), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda u, x: QuotientOfLinearsQ(u, x)))
    rule300 = ReplacementRule(pattern300, lambda n, u, x : Module(List(Set(lst, QuotientOfLinearsParts(u, x))), Int(Cos((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**n, x)))
    rubi.add(rule300)

    pattern301 = Pattern(Integral(Sin(v_)**WC('p', S(1))*Sin(w_)**WC('q', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda v, w: ZeroQ(v - w)))
    rule301 = ReplacementRule(pattern301, lambda q, p, u, w, x, v : Int(u*Sin(v)**(p + q), x))
    rubi.add(rule301)

    pattern302 = Pattern(Integral(Cos(v_)**WC('p', S(1))*Cos(w_)**WC('q', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda v, w: ZeroQ(v - w)))
    rule302 = ReplacementRule(pattern302, lambda q, p, u, w, x, v : Int(u*Cos(v)**(p + q), x))
    rubi.add(rule302)

    pattern303 = Pattern(Integral(Sin(v_)**WC('p', S(1))*Sin(w_)**WC('q', S(1)), x_), CustomConstraint(lambda v, w, x: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(Cancel(v/w), x))), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)))
    rule303 = ReplacementRule(pattern303, lambda q, p, w, x, v : Int(ExpandTrigReduce(Sin(v)**p*Sin(w)**q, x), x))
    rubi.add(rule303)

    pattern304 = Pattern(Integral(Cos(v_)**WC('p', S(1))*Cos(w_)**WC('q', S(1)), x_), CustomConstraint(lambda v, w, x: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(Cancel(v/w), x))), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)))
    rule304 = ReplacementRule(pattern304, lambda q, p, w, x, v : Int(ExpandTrigReduce(Cos(v)**p*Cos(w)**q, x), x))
    rubi.add(rule304)

    pattern305 = Pattern(Integral(x_**WC('m', S(1))*Sin(v_)**WC('p', S(1))*Sin(w_)**WC('q', S(1)), x_), CustomConstraint(lambda p, q, m: PositiveIntegerQ(m, p, q)), CustomConstraint(lambda v, w, x: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(Cancel(v/w), x))))
    rule305 = ReplacementRule(pattern305, lambda q, m, p, w, x, v : Int(ExpandTrigReduce(x**m, Sin(v)**p*Sin(w)**q, x), x))
    rubi.add(rule305)

    pattern306 = Pattern(Integral(x_**WC('m', S(1))*Cos(v_)**WC('p', S(1))*Cos(w_)**WC('q', S(1)), x_), CustomConstraint(lambda p, q, m: PositiveIntegerQ(m, p, q)), CustomConstraint(lambda v, w, x: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(Cancel(v/w), x))))
    rule306 = ReplacementRule(pattern306, lambda q, m, p, w, x, v : Int(ExpandTrigReduce(x**m, Cos(v)**p*Cos(w)**q, x), x))
    rubi.add(rule306)

    pattern307 = Pattern(Integral(Cos(w_)**WC('p', S(1))*Sin(v_)**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda v, w: ZeroQ(v - w)), CustomConstraint(lambda p: IntegerQ(p)))
    rule307 = ReplacementRule(pattern307, lambda p, u, w, x, v : S(2)**(-p)*Int(u*Sin(S(2)*v)**p, x))
    rubi.add(rule307)

    pattern308 = Pattern(Integral(Cos(w_)**WC('q', S(1))*Sin(v_)**WC('p', S(1)), x_), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)), CustomConstraint(lambda v, w, x: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(Cancel(v/w), x))))
    rule308 = ReplacementRule(pattern308, lambda q, p, w, x, v : Int(ExpandTrigReduce(Cos(w)**q*Sin(v)**p, x), x))
    rubi.add(rule308)

    pattern309 = Pattern(Integral(x_**WC('m', S(1))*Cos(w_)**WC('q', S(1))*Sin(v_)**WC('p', S(1)), x_), CustomConstraint(lambda p, q, m: PositiveIntegerQ(m, p, q)), CustomConstraint(lambda v, w, x: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(Cancel(v/w), x))))
    rule309 = ReplacementRule(pattern309, lambda q, m, p, w, x, v : Int(ExpandTrigReduce(x**m, Cos(w)**q*Sin(v)**p, x), x))
    rubi.add(rule309)

    pattern310 = Pattern(Integral(Sin(v_)*Tan(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule310 = ReplacementRule(pattern310, lambda n, v, x, w : Cos(v - w)*Int(Sec(w)*Tan(w)**(n + S(-1)), x) - Int(Cos(v)*Tan(w)**(n + S(-1)), x))
    rubi.add(rule310)

    pattern311 = Pattern(Integral(Cos(v_)*Cot(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule311 = ReplacementRule(pattern311, lambda n, v, x, w : Cos(v - w)*Int(Cot(w)**(n + S(-1))*Csc(w), x) - Int(Cot(w)**(n + S(-1))*Sin(v), x))
    rubi.add(rule311)

    pattern312 = Pattern(Integral(Cot(w_)**WC('n', S(1))*Sin(v_), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule312 = ReplacementRule(pattern312, lambda v, n, w, x : Int(Cos(v)*Cot(w)**(n + S(-1)), x) + Int(Cot(w)**(n + S(-1))*Csc(w), x)*Sin(v - w))
    rubi.add(rule312)

    pattern313 = Pattern(Integral(Cos(v_)*Tan(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule313 = ReplacementRule(pattern313, lambda n, v, x, w : -Int(Sec(w)*Tan(w)**(n + S(-1)), x)*Sin(v - w) + Int(Sin(v)*Tan(w)**(n + S(-1)), x))
    rubi.add(rule313)

    pattern314 = Pattern(Integral(Sec(w_)**WC('n', S(1))*Sin(v_), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule314 = ReplacementRule(pattern314, lambda v, n, w, x : Cos(v - w)*Int(Sec(w)**(n + S(-1))*Tan(w), x) + Int(Sec(w)**(n + S(-1)), x)*Sin(v - w))
    rubi.add(rule314)

    pattern315 = Pattern(Integral(Cos(v_)*Csc(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule315 = ReplacementRule(pattern315, lambda n, v, x, w : Cos(v - w)*Int(Cot(w)*Csc(w)**(n + S(-1)), x) - Int(Csc(w)**(n + S(-1)), x)*Sin(v - w))
    rubi.add(rule315)

    pattern316 = Pattern(Integral(Csc(w_)**WC('n', S(1))*Sin(v_), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule316 = ReplacementRule(pattern316, lambda v, n, w, x : Cos(v - w)*Int(Csc(w)**(n + S(-1)), x) + Int(Cot(w)*Csc(w)**(n + S(-1)), x)*Sin(v - w))
    rubi.add(rule316)

    pattern317 = Pattern(Integral(Cos(v_)*Sec(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule317 = ReplacementRule(pattern317, lambda n, v, x, w : Cos(v - w)*Int(Sec(w)**(n + S(-1)), x) - Int(Sec(w)**(n + S(-1))*Tan(w), x)*Sin(v - w))
    rubi.add(rule317)

    pattern318 = Pattern(Integral((a_ + Cos(x_*WC('d', S(1)) + WC('c', S(0)))*Sin(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule318 = ReplacementRule(pattern318, lambda a, m, d, f, n, c, x, e, b : Int((a + b*Sin(S(2)*c + S(2)*d*x)/S(2))**n*(e + f*x)**m, x))
    rubi.add(rule318)

    pattern319 = Pattern(Integral(x_**WC('m', S(1))*(a_ + Sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, a: NonzeroQ(a + b)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: Less(n, S(0))), CustomConstraint(lambda n, m: Equal(n, S(-1)) | (Equal(m, S(1)) & Equal(n, S(-2)))))
    rule319 = ReplacementRule(pattern319, lambda a, m, d, n, c, x, b : S(2)**(-n)*Int(x**m*(S(2)*a - b*Cos(S(2)*c + S(2)*d*x) + b)**n, x))
    rubi.add(rule319)

    pattern320 = Pattern(Integral(x_**WC('m', S(1))*(a_ + Cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, a: NonzeroQ(a + b)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: Less(n, S(0))), CustomConstraint(lambda n, m: Equal(n, S(-1)) | (Equal(m, S(1)) & Equal(n, S(-2)))))
    rule320 = ReplacementRule(pattern320, lambda a, m, d, n, c, x, b : S(2)**(-n)*Int(x**m*(S(2)*a + b*Cos(S(2)*c + S(2)*d*x) + b)**n, x))
    rubi.add(rule320)

    pattern321 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Sin((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p: RationalQ(p)))
    rule321 = ReplacementRule(pattern321, lambda m, d, f, n, p, c, a, e, x, b : d**(-m + S(-1))*Subst(Int((-c*f + d*e + f*x)**m*Sin(a + b*x**n)**p, x), x, c + d*x))
    rubi.add(rule321)

    pattern322 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cos((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p: RationalQ(p)))
    rule322 = ReplacementRule(pattern322, lambda m, d, f, n, p, c, a, e, x, b : d**(-m + S(-1))*Subst(Int((-c*f + d*e + f*x)**m*Cos(a + b*x**n)**p, x), x, c + d*x))
    rubi.add(rule322)

    pattern323 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(Cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*WC('b', S(1)) + Sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*WC('c', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda b, a: NonzeroQ(a + b)), CustomConstraint(lambda c, a: NonzeroQ(a + c)))
    rule323 = ReplacementRule(pattern323, lambda m, d, f, g, a, e, c, x, b : S(2)*Int((f + g*x)**m/(S(2)*a + b + c + (b - c)*Cos(S(2)*d + S(2)*e*x)), x))
    rubi.add(rule323)

    pattern324 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*Sec(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(b_ + Tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*WC('c', S(1))), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule324 = ReplacementRule(pattern324, lambda m, d, f, g, c, e, x, b : S(2)*Int((f + g*x)**m/(b + c + (b - c)*Cos(S(2)*d + S(2)*e*x)), x))
    rubi.add(rule324)

    pattern325 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*Sec(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(Sec(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*WC('a', S(1)) + Tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*WC('c', S(1)) + WC('b', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda b, a: NonzeroQ(a + b)), CustomConstraint(lambda c, a: NonzeroQ(a + c)))
    rule325 = ReplacementRule(pattern325, lambda m, d, f, g, a, e, c, x, b : S(2)*Int((f + g*x)**m/(S(2)*a + b + c + (b - c)*Cos(S(2)*d + S(2)*e*x)), x))
    rubi.add(rule325)

    pattern326 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*Csc(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(c_ + Cot(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule326 = ReplacementRule(pattern326, lambda m, d, f, g, x, e, c, b : S(2)*Int((f + g*x)**m/(b + c + (b - c)*Cos(S(2)*d + S(2)*e*x)), x))
    rubi.add(rule326)

    pattern327 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*Csc(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(Cot(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*WC('b', S(1)) + Csc(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*WC('a', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda b, a: NonzeroQ(a + b)), CustomConstraint(lambda c, a: NonzeroQ(a + c)))
    rule327 = ReplacementRule(pattern327, lambda m, d, f, c, a, e, g, x, b : S(2)*Int((f + g*x)**m/(S(2)*a + b + c + (b - c)*Cos(S(2)*d + S(2)*e*x)), x))
    rubi.add(rule327)

    pattern328 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cos(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + Sin(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda b, a: PosQ(a**S(2) - b**S(2))))
    rule328 = ReplacementRule(pattern328, lambda a, m, d, f, c, x, e, b : -ImaginaryI*(e + f*x)**(m + S(1))/(b*f*(m + S(1))) + Int((e + f*x)**m*exp(ImaginaryI*(c + d*x))/(-ImaginaryI*b*exp(ImaginaryI*(c + d*x)) + a - Rt(a**S(2) - b**S(2), S(2))), x) + Int((e + f*x)**m*exp(ImaginaryI*(c + d*x))/(-ImaginaryI*b*exp(ImaginaryI*(c + d*x)) + a + Rt(a**S(2) - b**S(2), S(2))), x))
    rubi.add(rule328)

    pattern329 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Sin(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + Cos(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda b, a: PosQ(a**S(2) - b**S(2))))
    rule329 = ReplacementRule(pattern329, lambda a, m, d, f, c, x, e, b : -ImaginaryI*Int((e + f*x)**m*exp(ImaginaryI*(c + d*x))/(a + b*exp(ImaginaryI*(c + d*x)) - Rt(a**S(2) - b**S(2), S(2))), x) - ImaginaryI*Int((e + f*x)**m*exp(ImaginaryI*(c + d*x))/(a + b*exp(ImaginaryI*(c + d*x)) + Rt(a**S(2) - b**S(2), S(2))), x) + ImaginaryI*(e + f*x)**(m + S(1))/(b*f*(m + S(1))))
    rubi.add(rule329)

    pattern330 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cos(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + Sin(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda b, a: NegQ(a**S(2) - b**S(2))))
    rule330 = ReplacementRule(pattern330, lambda a, m, d, f, c, x, e, b : ImaginaryI*Int((e + f*x)**m*exp(ImaginaryI*(c + d*x))/(ImaginaryI*a + b*exp(ImaginaryI*(c + d*x)) - Rt(-a**S(2) + b**S(2), S(2))), x) + ImaginaryI*Int((e + f*x)**m*exp(ImaginaryI*(c + d*x))/(ImaginaryI*a + b*exp(ImaginaryI*(c + d*x)) + Rt(-a**S(2) + b**S(2), S(2))), x) - ImaginaryI*(e + f*x)**(m + S(1))/(b*f*(m + S(1))))
    rubi.add(rule330)

    pattern331 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Sin(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + Cos(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda b, a: NegQ(a**S(2) - b**S(2))))
    rule331 = ReplacementRule(pattern331, lambda a, m, d, f, c, x, e, b : ImaginaryI*(e + f*x)**(m + S(1))/(b*f*(m + S(1))) + Int((e + f*x)**m*exp(ImaginaryI*(c + d*x))/(ImaginaryI*a + ImaginaryI*b*exp(ImaginaryI*(c + d*x)) - Rt(-a**S(2) + b**S(2), S(2))), x) + Int((e + f*x)**m*exp(ImaginaryI*(c + d*x))/(ImaginaryI*a + ImaginaryI*b*exp(ImaginaryI*(c + d*x)) + Rt(-a**S(2) + b**S(2), S(2))), x))
    rubi.add(rule331)

    pattern332 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + Sin(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))))
    rule332 = ReplacementRule(pattern332, lambda a, m, d, f, n, c, x, e, b : -Int((e + f*x)**m*Cos(c + d*x)**(n + S(-2))*Sin(c + d*x), x)/b + Int((e + f*x)**m*Cos(c + d*x)**(n + S(-2)), x)/a)
    rubi.add(rule332)

    pattern333 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + Cos(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))))
    rule333 = ReplacementRule(pattern333, lambda a, m, d, f, n, c, x, e, b : -Int((e + f*x)**m*Cos(c + d*x)*Sin(c + d*x)**(n + S(-2)), x)/b + Int((e + f*x)**m*Sin(c + d*x)**(n + S(-2)), x)/a)
    rubi.add(rule333)

    pattern334 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + Sin(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))))
    rule334 = ReplacementRule(pattern334, lambda a, m, d, f, n, c, x, e, b : a*Int((e + f*x)**m*Cos(c + d*x)**(n + S(-2)), x)/b**S(2) - Int((e + f*x)**m*Cos(c + d*x)**(n + S(-2))*Sin(c + d*x), x)/b - (a**S(2) - b**S(2))*Int((e + f*x)**m*Cos(c + d*x)**(n + S(-2))/(a + b*Sin(c + d*x)), x)/b**S(2))
    rubi.add(rule334)

    pattern335 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + Cos(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))))
    rule335 = ReplacementRule(pattern335, lambda a, m, d, f, n, c, x, e, b : a*Int((e + f*x)**m*Sin(c + d*x)**(n + S(-2)), x)/b**S(2) - Int((e + f*x)**m*Cos(c + d*x)*Sin(c + d*x)**(n + S(-2)), x)/b - (a**S(2) - b**S(2))*Int((e + f*x)**m*Sin(c + d*x)**(n + S(-2))/(a + b*Cos(c + d*x)), x)/b**S(2))
    rubi.add(rule335)

    pattern336 = Pattern(Integral((A_ + Sin(x_*WC('d', S(1)) + WC('c', S(0)))*WC('B', S(1)))*(x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + Sin(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda B, A, b, a: ZeroQ(A*a - B*b)))
    rule336 = ReplacementRule(pattern336, lambda a, d, f, A, c, x, e, b, B : B*f*Int(Cos(c + d*x)/(a + b*Sin(c + d*x)), x)/(a*d) - B*(e + f*x)*Cos(c + d*x)/(a*d*(a + b*Sin(c + d*x))))
    rubi.add(rule336)

    pattern337 = Pattern(Integral((A_ + Cos(x_*WC('d', S(1)) + WC('c', S(0)))*WC('B', S(1)))*(x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + Cos(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda B, A, b, a: ZeroQ(A*a - B*b)))
    rule337 = ReplacementRule(pattern337, lambda a, d, f, A, c, x, e, b, B : -B*f*Int(Sin(c + d*x)/(a + b*Cos(c + d*x)), x)/(a*d) + B*(e + f*x)*Sin(c + d*x)/(a*d*(a + b*Cos(c + d*x))))
    rubi.add(rule337)

    pattern338 = Pattern(Integral((a_ + Tan(v_)*WC('b', S(1)))**WC('n', S(1))*Sec(v_)**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Equal(m + n, S(0))), CustomConstraint(lambda m: OddQ(m)))
    rule338 = ReplacementRule(pattern338, lambda m, n, b, a, x, v : Int((a*Cos(v) + b*Sin(v))**n, x))
    rubi.add(rule338)

    pattern339 = Pattern(Integral((a_ + Cot(v_)*WC('b', S(1)))**WC('n', S(1))*Csc(v_)**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Equal(m + n, S(0))), CustomConstraint(lambda m: OddQ(m)))
    rule339 = ReplacementRule(pattern339, lambda m, n, b, a, x, v : Int((a*Sin(v) + b*Cos(v))**n, x))
    rubi.add(rule339)

    pattern340 = Pattern(Integral(Sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*Sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule340 = ReplacementRule(pattern340, lambda m, d, n, b, c, a, x, u : Int(ExpandTrigReduce(u, Sin(a + b*x)**m*Sin(c + d*x)**n, x), x))
    rubi.add(rule340)

    pattern341 = Pattern(Integral(Cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*Cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule341 = ReplacementRule(pattern341, lambda m, d, n, u, c, a, x, b : Int(ExpandTrigReduce(u, Cos(a + b*x)**m*Cos(c + d*x)**n, x), x))
    rubi.add(rule341)

    pattern342 = Pattern(Integral(Sec(c_ + x_*WC('d', S(1)))*Sec(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, d: ZeroQ(b**S(2) - d**S(2))), CustomConstraint(lambda b, d, c, a: NonzeroQ(-a*d + b*c)))
    rule342 = ReplacementRule(pattern342, lambda a, d, c, x, b : Csc((-a*d + b*c)/b)*Int(Tan(c + d*x), x) - Csc((-a*d + b*c)/d)*Int(Tan(a + b*x), x))
    rubi.add(rule342)

    pattern343 = Pattern(Integral(Csc(c_ + x_*WC('d', S(1)))*Csc(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, d: ZeroQ(b**S(2) - d**S(2))), CustomConstraint(lambda b, d, c, a: NonzeroQ(-a*d + b*c)))
    rule343 = ReplacementRule(pattern343, lambda a, d, c, x, b : Csc((-a*d + b*c)/b)*Int(Cot(a + b*x), x) - Csc((-a*d + b*c)/d)*Int(Cot(c + d*x), x))
    rubi.add(rule343)

    pattern344 = Pattern(Integral(Tan(c_ + x_*WC('d', S(1)))*Tan(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, d: ZeroQ(b**S(2) - d**S(2))), CustomConstraint(lambda b, d, c, a: NonzeroQ(-a*d + b*c)))
    rule344 = ReplacementRule(pattern344, lambda a, d, c, x, b : -b*x/d + b*Cos((-a*d + b*c)/d)*Int(Sec(a + b*x)*Sec(c + d*x), x)/d)
    rubi.add(rule344)

    pattern345 = Pattern(Integral(Cot(c_ + x_*WC('d', S(1)))*Cot(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, d: ZeroQ(b**S(2) - d**S(2))), CustomConstraint(lambda b, d, c, a: NonzeroQ(-a*d + b*c)))
    rule345 = ReplacementRule(pattern345, lambda a, d, c, x, b : -b*x/d + Cos((-a*d + b*c)/d)*Int(Csc(a + b*x)*Csc(c + d*x), x))
    rubi.add(rule345)

    pattern346 = Pattern(Integral((Cos(v_)*WC('a', S(1)) + Sin(v_)*WC('b', S(1)))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) + b**S(2))))
    rule346 = ReplacementRule(pattern346, lambda n, b, u, a, x, v : Int(u*(a*exp(-a*v/b))**n, x))
    rubi.add(rule346)

    return rubi
