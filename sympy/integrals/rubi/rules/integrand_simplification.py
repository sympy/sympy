
from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (Int, Set, With, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcCsch, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct, SplitSum, Complex, UnsameQ, _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify, _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum, _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux, TrigSimplifyAux)
    from sympy import Integral, S, sqrt
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)


    A_, B_, C_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, z_ = [WC(i) for i in 'ABCabcdefghijklmnpqrtuvswxz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, jn_, mn_, non2_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'n1', 'n2', 'n3',' Pq', 'Pm', ' Px', 'Qm', 'Qr', 'jn', 'mn', 'non2']]
    p, q, r, s, mn, gcd, P, Q = symbols('p q r s mn gcd P Q')

    _UseGamma = False

def integrand_simplification(rubi):

    pattern1 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a: ZeroQ(a)))
    rule1 = ReplacementRule(pattern1, lambda u, a, p, x, b, n : Int(u*(b*x**n)**p, x))
    rubi.add(rule1)

    pattern2 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b: ZeroQ(b)))
    rule2 = ReplacementRule(pattern2, lambda u, a, p, x, b, n : Int(a**p*u, x))
    rubi.add(rule2)

    pattern3 = Pattern(Integral((a_ + x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, j: ZeroQ(j - S(2)*n)), CustomConstraint(lambda a: ZeroQ(a)))
    rule3 = ReplacementRule(pattern3, lambda u, a, j, p, x, b, c, n : Int(u*(b*x**n + c*x**(S(2)*n))**p, x))
    rubi.add(rule3)

    pattern4 = Pattern(Integral((x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, j: ZeroQ(j - S(2)*n)), CustomConstraint(lambda b: ZeroQ(b)))
    rule4 = ReplacementRule(pattern4, lambda u, a, j, p, x, b, c, n : Int(u*(a + c*x**(S(2)*n))**p, x))
    rubi.add(rule4)

    pattern5 = Pattern(Integral((x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, j: ZeroQ(j - S(2)*n)), CustomConstraint(lambda c: ZeroQ(c)))
    rule5 = ReplacementRule(pattern5, lambda u, a, j, p, x, b, c, n : Int(u*(a + b*x**n)**p, x))
    rubi.add(rule5)

    pattern6 = Pattern(Integral(Pm_**p_*WC('u', S(1)), x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, Pm: PolyQ(Pm, x)), CustomConstraint(lambda p: ~(RationalQ(p))), CustomConstraint(lambda p: RationalQ(Simplify(p))))
    rule6 = ReplacementRule(pattern6, lambda p, u, x, Pm : Int(Pm**Simplify(p)*u, x))
    rubi.add(rule6)

    pattern7 = Pattern(Integral(a_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)))
    rule7 = ReplacementRule(pattern7, lambda a, x : a*x)
    rubi.add(rule7)

    pattern8 = Pattern(Integral(a_*(b_ + x_*WC('c', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule8 = ReplacementRule(pattern8, lambda b, c, a, x : a*(b + c*x)**S(2)/(S(2)*c))
    rubi.add(rule8)

    pattern9 = Pattern(Integral(-u_, x_))
    rule9 = ReplacementRule(pattern9, lambda u, x : Int(u, x)*I)
    rubi.add(rule9)

    pattern10 = Pattern(Integral(u_*Complex(S(0), a_), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda a: EqQ(a**S(2), S(1))))
    rule10 = ReplacementRule(pattern10, lambda a, u, x : Int(u, x)*Complex(I, a))
    rubi.add(rule10)

    pattern11 = Pattern(Integral(a_*u_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)))
    rule11 = ReplacementRule(pattern11, lambda a, u, x : a*Int(u, x))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(u_, x_), CustomConstraint(lambda u: SumQ(u)))
    rule12 = ReplacementRule(pattern12, lambda u, x : IntSum(u, x))
    rubi.add(rule12)

    pattern13 = Pattern(Integral(v_**WC('m', S(1))*(b_*v_)**n_*WC('u', S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule13 = ReplacementRule(pattern13, lambda v, u, m, x, b, n : b**(-m)*Int(u*(b*v)**(m + n), x))
    rubi.add(rule13)

    pattern14 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda n: PositiveIntegerQ(n + S(1)/2)), CustomConstraint(lambda n, m: IntegerQ(m + n)))
    rule14 = ReplacementRule(pattern14, lambda v, u, m, a, x, b, n : a**(m + S(1)/2)*b**(n + S(-1)/2)*Int(u*v**(m + n), x)*Sqrt(b*v)/Sqrt(a*v))
    rubi.add(rule14)

    pattern15 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda n: NegativeIntegerQ(n + S(-1)/2)), CustomConstraint(lambda n, m: IntegerQ(m + n)))
    rule15 = ReplacementRule(pattern15, lambda v, u, m, a, x, b, n : a**(m + S(-1)/2)*b**(n + S(1)/2)*Int(u*v**(m + n), x)*Sqrt(a*v)/Sqrt(b*v))
    rubi.add(rule15)

    pattern16 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda n: ~(IntegerQ(n))), CustomConstraint(lambda n, m: IntegerQ(m + n)))
    rule16 = ReplacementRule(pattern16, lambda v, u, m, a, x, b, n : a**(m + n)*(a*v)**(-n)*(b*v)**n*Int(u*v**(m + n), x))
    rubi.add(rule16)

    pattern17 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda n: ~(IntegerQ(n))), CustomConstraint(lambda n, m: ~(IntegerQ(m + n))))
    rule17 = ReplacementRule(pattern17, lambda v, u, m, a, x, b, n : a**(-IntPart(n))*b**IntPart(n)*(a*v)**(-FracPart(n))*(b*v)**FracPart(n)*Int(u*(a*v)**(m + n), x))
    rubi.add(rule17)

    pattern18 = Pattern(Integral((a_ + v_*WC('b', S(1)))**WC('m', S(1))*(c_ + v_*WC('d', S(1)))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, a, d: ZeroQ(-a*d + b*c)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d, a, x, b, c, n: ~(IntegerQ(n)) | SimplerQ(c + d*x, a + b*x)))
    rule18 = ReplacementRule(pattern18, lambda v, u, m, d, a, x, b, c, n : (b/d)**m*Int(u*(c + d*v)**(m + n), x))
    rubi.add(rule18)

    pattern19 = Pattern(Integral((a_ + v_*WC('b', S(1)))**m_*(c_ + v_*WC('d', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, a, d: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, d: PositiveQ(b/d)), CustomConstraint(lambda n, m: ~(IntegerQ(m) | IntegerQ(n))))
    rule19 = ReplacementRule(pattern19, lambda v, u, m, d, a, x, b, c, n : (b/d)**m*Int(u*(c + d*v)**(m + n), x))
    rubi.add(rule19)

    pattern20 = Pattern(Integral((a_ + v_*WC('b', S(1)))**m_*(c_ + v_*WC('d', S(1)))**n_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, c, a, d: ZeroQ(-a*d + b*c)), CustomConstraint(lambda b, n, m, d: ~(IntegerQ(m) | IntegerQ(n) | PositiveQ(b/d))))
    rule20 = ReplacementRule(pattern20, lambda v, u, m, d, a, x, b, c, n : (a + b*v)**m*(c + d*v)**(-m)*Int(u*(c + d*v)**(m + n), x))
    rubi.add(rule20)

    pattern21 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_**S(2)*WC('c', S(1)) + v_*WC('b', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: LessEqual(m, S(-1))))
    rule21 = ReplacementRule(pattern21, lambda v, u, m, a, x, b, c : Int(u*(a*v)**(m + S(1))*(b + c*v), x)/a)
    rubi.add(rule21)

    pattern22 = Pattern(Integral((a_ + v_*WC('b', S(1)))**m_*(v_**S(2)*WC('C', S(1)) + v_*WC('B', S(1)) + WC('A', S(0)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda A, a, C, b, B: ZeroQ(A*b**S(2) - B*a*b + C*a**S(2))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: LessEqual(m, S(-1))))
    rule22 = ReplacementRule(pattern22, lambda v, A, u, m, a, x, C, b, B : Int(u*(a + b*v)**(m + S(1))*Simp(B*b - C*a + C*b*v, x), x)/b**S(2))
    rubi.add(rule22)

    pattern23 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('m', S(1))*(c_ + x_**WC('q', S(1))*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, q: ZeroQ(n + q)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda b, c, a, d: ZeroQ(a*c - b*d)), CustomConstraint(lambda n, m: ~(IntegerQ(m) & NegQ(n))))
    rule23 = ReplacementRule(pattern23, lambda u, m, d, a, p, x, b, c, n, q : (d/a)**p*Int(u*x**(-n*p)*(a + b*x**n)**(m + p), x))
    rubi.add(rule23)

    pattern24 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('m', S(1))*(c_ + x_**j_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, j: ZeroQ(j - S(2)*n)), CustomConstraint(lambda p, m: ZeroQ(m + p)), CustomConstraint(lambda b, c, a, d: ZeroQ(a**S(2)*d + b**S(2)*c)), CustomConstraint(lambda a: PositiveQ(a)), CustomConstraint(lambda d: NegativeQ(d)))
    rule24 = ReplacementRule(pattern24, lambda u, m, d, a, j, p, x, b, c, n : (-b**S(2)/d)**m*Int(u*(a - b*x**n)**(-m), x))
    rubi.add(rule24)

    pattern25 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, c, a: ZeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule25 = ReplacementRule(pattern25, lambda u, a, p, x, b, c : Int(u*Cancel(c**(-p)*(b/S(2) + c*x)**(S(2)*p)), x))
    rubi.add(rule25)

    pattern26 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda b, c, a: ZeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule26 = ReplacementRule(pattern26, lambda u, a, p, x, b, c, n, n2 : c**(-p)*Int(u*(b/S(2) + c*x**n)**(S(2)*p), x))
    rubi.add(rule26)

    pattern27 = Pattern(Integral((d_ + x_*WC('e', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, c, e, d: ZeroQ(-b*e + S(2)*c*d)))
    rule27 = ReplacementRule(pattern27, lambda d, a, p, e, x, b, c : d*Subst(Int(x**p, x), x, a + b*x + c*x**S(2))/b)
    rubi.add(rule27)

    pattern28 = Pattern(Integral((x_**WC('p', S(1))*WC('a', S(1)) + x_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p, q: PosQ(-p + q)))
    rule28 = ReplacementRule(pattern28, lambda u, m, a, p, x, b, q : Int(u*x**(m*p)*(a + b*x**(-p + q))**m, x))
    rubi.add(rule28)

    pattern29 = Pattern(Integral((x_**WC('p', S(1))*WC('a', S(1)) + x_**WC('q', S(1))*WC('b', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p, q: PosQ(-p + q)), CustomConstraint(lambda p, r: PosQ(-p + r)))
    rule29 = ReplacementRule(pattern29, lambda u, m, a, p, r, x, b, c, q : Int(u*x**(m*p)*(a + b*x**(-p + q) + c*x**(-p + r))**m, x))
    rubi.add(rule29)

    pattern30 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))))
    rule30 = ReplacementRule(pattern30, lambda m, a, x, b, n : Log(RemoveContent(a + b*x**n, x))/(b*n))
    rubi.add(rule30)

    pattern31 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule31 = ReplacementRule(pattern31, lambda m, a, p, x, b, n : (a + b*x**n)**(p + S(1))/(b*n*(p + S(1))))
    rubi.add(rule31)

    pattern32 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**p_, x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b1, b2, a2, a1: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda n, m: ZeroQ(m - S(2)*n + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule32 = ReplacementRule(pattern32, lambda m, a2, b1, p, b2, x, n, a1 : (a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*b1*b2*n*(p + S(1))))
    rubi.add(rule32)

    pattern33 = Pattern(Integral(Qm_*(Pm_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, Pm: PolyQ(Pm, x)), CustomConstraint(lambda Qm, x: PolyQ(Qm, x)))
    rule33 = ReplacementRule(pattern33, lambda a, Qm, p, x, Pm, b, n : Module(List(Set(m, Expon(Pm, x))), Condition(Coeff(Qm, x, m + S(-1))*Subst(Int((a + b*x**n)**p, x), x, Pm)/(m*Coeff(Pm, x, m)), Equal(Expon(Qm, x), m + S(-1)) & ZeroQ(-Qm*m*Coeff(Pm, x, m) + Coeff(Qm, x, m + S(-1))*D(Pm, x)))))
    rubi.add(rule33)

    pattern34 = Pattern(Integral(Qm_*(Pm_**WC('n', S(1))*WC('b', S(1)) + Pm_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda x, Pm: PolyQ(Pm, x)), CustomConstraint(lambda Qm, x: PolyQ(Qm, x)))
    rule34 = ReplacementRule(pattern34, lambda a, Qm, p, x, Pm, b, c, n, n2 : Module(List(Set(m, Expon(Pm, x))), Condition(Coeff(Qm, x, m + S(-1))*Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, Pm)/(m*Coeff(Pm, x, m)), Equal(Expon(Qm, x), m + S(-1)) & ZeroQ(-Qm*m*Coeff(Pm, x, m) + Coeff(Qm, x, m + S(-1))*D(Pm, x)))))
    rubi.add(rule34)

    pattern35 = Pattern(Integral(Pq_**m_*Qr_**p_*WC('u', S(1)), x_), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda x, Qr: PolyQ(Qr, x)))
    rule35 = ReplacementRule(pattern35, lambda Pq, u, m, Qr, p, x : Module(List(Set(gcd, PolyGCD(Pq, Qr, x))), Condition(Int(gcd**(m + p)*u*PolynomialQuotient(Pq, gcd, x)**m*PolynomialQuotient(Qr, gcd, x)**p, x), NonzeroQ(gcd + S(-1)))))
    rubi.add(rule35)

    pattern36 = Pattern(Integral(Pq_*Qr_**p_*WC('u', S(1)), x_), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda x, Qr: PolyQ(Qr, x)))
    rule36 = ReplacementRule(pattern36, lambda Pq, u, Qr, p, x : Module(List(Set(gcd, PolyGCD(Pq, Qr, x))), Condition(Int(gcd**(p + S(1))*u*PolynomialQuotient(Pq, gcd, x)*PolynomialQuotient(Qr, gcd, x)**p, x), NonzeroQ(gcd + S(-1)))))
    rubi.add(rule36)

    return rubi
