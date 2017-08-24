
from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (Gamma, Part, Module, Int, Set, With, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcCsch, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct, SplitSum, Complex, UnsameQ, _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify, _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum, _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux, TrigSimplifyAux)
    from sympy import Integral, S, sqrt
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)


    A_, B_, C_, F_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, jn_, mn_, non2_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'n1', 'n2', 'n3',' Pq', 'Pm', ' Px', 'Qm', 'Qr', 'jn', 'mn', 'non2']]
    p, q, r, s, mn, gcd, P, Q, lst = symbols('p q r s mn gcd P Q lst')

    _UseGamma = False

def miscellaneous_integration(rubi):

    pattern1 = Pattern(Integral(((d_*(x_*WC('b', S(1)) + WC('a', S(0))))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda q: ~(IntegerQ(q))))
    rule1 = ReplacementRule(pattern1, lambda b, d, q, a, p, c, u, x : (c*(d*(a + b*x))**p)**q*(a + b*x)**(-p*q)*Int(u*(a + b*x)**(p*q), x))
    rubi.add(rule1)

    pattern2 = Pattern(Integral((((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('d', S(1)))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda q: ~(IntegerQ(q))))
    rule2 = ReplacementRule(pattern2, lambda b, d, q, n, a, p, c, u, x : (c*(d*(a + b*x)**n)**p)**q*(a + b*x)**(-n*p*q)*Int(u*(a + b*x)**(n*p*q), x))
    rubi.add(rule2)

    pattern3 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda g, e: ZeroQ(e + g)), CustomConstraint(lambda f, d: ZeroQ(d + f + S(-2))), CustomConstraint(lambda d, A, f, e, C: ZeroQ(A*e**S(2) + C*d*f)), CustomConstraint(lambda B, e, d, C: ZeroQ(-B*e + S(2)*C*(d + S(-1)))))
    rule3 = ReplacementRule(pattern3, lambda b, d, n, g, a, A, f, c, x, e, B, F, C : g*Subst(Int((a + b*F(c*x))**n/x, x), x, Sqrt(d + e*x)/Sqrt(f + g*x))/C)
    rubi.add(rule3)

    pattern4 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + S(1))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda g, e: ZeroQ(e + g)), CustomConstraint(lambda e, A, C: ZeroQ(A*e**S(2) + C)))
    rule4 = ReplacementRule(pattern4, lambda b, n, g, a, A, c, e, x, F, C : g*Subst(Int((a + b*F(c*x))**n/x, x), x, Sqrt(e*x + S(1))/Sqrt(g*x + S(1)))/C)
    rubi.add(rule4)

    pattern5 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda g, e: ZeroQ(e + g)), CustomConstraint(lambda f, d: ZeroQ(d + f + S(-2))), CustomConstraint(lambda d, A, f, e, C: ZeroQ(A*e**S(2) + C*d*f)), CustomConstraint(lambda B, e, d, C: ZeroQ(-B*e + S(2)*C*(d + S(-1)))))
    rule5 = ReplacementRule(pattern5, lambda b, d, n, g, a, A, f, c, x, e, B, F, C : g*Subst(Int((F**(c*x)*b + a)**n/x, x), x, Sqrt(d + e*x)/Sqrt(f + g*x))/C)
    rubi.add(rule5)

    pattern6 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda g, e: ZeroQ(e + g)), CustomConstraint(lambda e, A, C: ZeroQ(A*e**S(2) + C)))
    rule6 = ReplacementRule(pattern6, lambda b, n, g, a, A, c, e, x, F, C : g*Subst(Int((F**(c*x)*b + a)**n/x, x), x, Sqrt(e*x + S(1))/Sqrt(g*x + S(1)))/C)
    rubi.add(rule6)

    pattern7 = Pattern(Integral(u_/y_, x_))
    rule7 = ReplacementRule(pattern7, lambda x, y, u : Module(List(Set(q, DerivativeDivides(y, u, x))), Condition(q*Log(RemoveContent(y, x)), ~(FalseQ(q)))))
    rubi.add(rule7)

    pattern8 = Pattern(Integral(u_/(w_*y_), x_))
    rule8 = ReplacementRule(pattern8, lambda x, w, y, u : Module(List(Set(q, DerivativeDivides(w*y, u, x))), Condition(q*Log(RemoveContent(w*y, x)), ~(FalseQ(q)))))
    rubi.add(rule8)

    pattern9 = Pattern(Integral(u_*y_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule9 = ReplacementRule(pattern9, lambda m, x, y, u : Module(List(Set(q, DerivativeDivides(y, u, x))), Condition(q*y**(m + S(1))/(m + S(1)), ~(FalseQ(q)))))
    rubi.add(rule9)

    pattern10 = Pattern(Integral(u_*y_**WC('m', S(1))*z_**WC('n', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule10 = ReplacementRule(pattern10, lambda n, z, m, y, u, x : Module(List(Set(q, DerivativeDivides(y*z, u*z**(-m + n), x))), Condition(q*y**(m + S(1))*z**(m + S(1))/(m + S(1)), ~(FalseQ(q)))))
    rubi.add(rule10)

    pattern11 = Pattern(Integral(u_, x_))
    rule11 = ReplacementRule(pattern11, lambda x, u : Module(List(Set(v, SimplifyIntegrand(u, x))), Condition(Int(v, x), SimplerIntegrandQ(v, u, x))))
    rubi.add(rule11)

    pattern12 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda b, f, d, e: ZeroQ(b*e**S(2) - d*f**S(2))))
    rule12 = ReplacementRule(pattern12, lambda b, d, n, a, f, c, m, e, u, x : (a*e**S(2) - c*f**S(2))**m*Int(ExpandIntegrand(u*(e*Sqrt(a + b*x**n) - f*Sqrt(c + d*x**n))**(-m), x), x))
    rubi.add(rule12)

    pattern13 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda a, f, c, e: ZeroQ(a*e**S(2) - c*f**S(2))))
    rule13 = ReplacementRule(pattern13, lambda b, d, n, a, f, c, m, e, u, x : (b*e**S(2) - d*f**S(2))**m*Int(ExpandIntegrand(u*x**(m*n)*(e*Sqrt(a + b*x**n) - f*Sqrt(c + d*x**n))**(-m), x), x))
    rubi.add(rule13)

    pattern14 = Pattern(Integral(u_**WC('m', S(1))*w_*(u_**n_*WC('a', S(1)) + v_)**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: ~(PositiveQ(n))), CustomConstraint(lambda v, x: NFreeQ(v, x)))
    rule14 = ReplacementRule(pattern14, lambda n, a, w, p, x, m, u, v : Int(u**(m + n*p)*w*(a + u**(-n)*v)**p, x))
    rubi.add(rule14)

    pattern15 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda v, y: ZeroQ(-v + y)))
    rule15 = ReplacementRule(pattern15, lambda b, d, n, a, c, x, m, y, u, v : Module(List(Set(q, DerivativeDivides(y, u, x))), Condition(q*Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, y), ~(FalseQ(q)))))
    rubi.add(rule15)

    pattern16 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda v, y: ZeroQ(-v + y)), CustomConstraint(lambda w, y: ZeroQ(-w + y)))
    rule16 = ReplacementRule(pattern16, lambda b, d, n, a, w, f, c, p, x, m, e, y, u, v : Module(List(Set(q, DerivativeDivides(y, u, x))), Condition(q*Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x, y), ~(FalseQ(q)))))
    rubi.add(rule16)

    pattern17 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(z_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda v, y: ZeroQ(-v + y)), CustomConstraint(lambda w, y: ZeroQ(-w + y)), CustomConstraint(lambda y, z: ZeroQ(y - z)))
    rule17 = ReplacementRule(pattern17, lambda b, d, q, n, g, a, w, f, c, p, z, m, e, x, y, u, v, h : Module(List(Set(r, DerivativeDivides(y, u, x))), Condition(r*Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x, y), ~(FalseQ(r)))))
    rubi.add(rule17)

    pattern18 = Pattern(Integral((a_ + y_**n_*WC('b', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule18 = ReplacementRule(pattern18, lambda b, n, a, y, u, x : Module(List(Set(q, DerivativeDivides(y, u, x))), Condition(a*Int(u, x) + b*q*Subst(Int(x**n, x), x, y), ~(FalseQ(q)))))
    rubi.add(rule18)

    pattern19 = Pattern(Integral((y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule19 = ReplacementRule(pattern19, lambda b, n, a, p, y, u, x : Module(List(Set(q, DerivativeDivides(y, u, x))), Condition(q*Subst(Int((a + b*x**n)**p, x), x, y), ~(FalseQ(q)))))
    rubi.add(rule19)

    pattern20 = Pattern(Integral(v_**WC('m', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule20 = ReplacementRule(pattern20, lambda b, n, a, p, x, m, y, u, v : Module(List(q, r), Condition(q*r*Subst(Int(x**m*(a + b*x**n)**p, x), x, y), ~(FalseQ(Set(q, DerivativeDivides(y, u, x)))) & ~(FalseQ(Set(r, Divides(y**m, v**m, x)))))))
    rubi.add(rule20)

    pattern21 = Pattern(Integral((v_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda v, y: ZeroQ(-v + y)))
    rule21 = ReplacementRule(pattern21, lambda b, n, a, n2, p, c, x, y, u, v : Module(List(Set(q, DerivativeDivides(y, u, x))), Condition(q*Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, y), ~(FalseQ(q)))))
    rubi.add(rule21)

    pattern22 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda v, y: ZeroQ(-v + y)), CustomConstraint(lambda w, y: ZeroQ(-w + y)))
    rule22 = ReplacementRule(pattern22, lambda b, n, a, n2, w, p, c, A, x, y, v, u, B : Module(List(Set(q, DerivativeDivides(y, u, x))), Condition(q*Subst(Int((A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), ~(FalseQ(q)))))
    rubi.add(rule22)

    pattern23 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda w, y: ZeroQ(-w + y)))
    rule23 = ReplacementRule(pattern23, lambda n, a, n2, w, p, c, A, x, y, u, B : Module(List(Set(q, DerivativeDivides(y, u, x))), Condition(q*Subst(Int((A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), ~(FalseQ(q)))))
    rubi.add(rule23)

    pattern24 = Pattern(Integral(v_**WC('m', S(1))*(w_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda w, y: ZeroQ(-w + y)))
    rule24 = ReplacementRule(pattern24, lambda b, n, a, n2, w, p, c, x, m, y, u, v : Module(List(q, r), Condition(q*r*Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), ~(FalseQ(Set(q, DerivativeDivides(y, u, x)))) & ~(FalseQ(Set(r, Divides(y**m, v**m, x)))))))
    rubi.add(rule24)

    pattern25 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda v, y: ZeroQ(-v + y)), CustomConstraint(lambda w, y: ZeroQ(-w + y)))
    rule25 = ReplacementRule(pattern25, lambda b, n, a, n2, w, p, c, A, z, m, x, y, v, u, B : Module(List(q, r), Condition(q*r*Subst(Int(x**m*(A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), ~(FalseQ(Set(q, DerivativeDivides(y, u, x)))) & ~(FalseQ(Set(r, Divides(y**m, z**m, x)))))))
    rubi.add(rule25)

    pattern26 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda w, y: ZeroQ(-w + y)))
    rule26 = ReplacementRule(pattern26, lambda n, a, n2, w, p, c, A, z, m, x, y, u, B : Module(List(q, r), Condition(q*r*Subst(Int(x**m*(A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), ~(FalseQ(Set(q, DerivativeDivides(y, u, x)))) & ~(FalseQ(Set(r, Divides(y**m, z**m, x)))))))
    rubi.add(rule26)

    pattern27 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda v, y: ZeroQ(-v + y)))
    rule27 = ReplacementRule(pattern27, lambda b, d, n, a, p, c, x, m, y, u, v : Module(List(Set(q, DerivativeDivides(y, u, x))), Condition(q*Subst(Int((a + b*x**n)**m*(c + d*x**n)**p, x), x, y), ~(FalseQ(q)))))
    rubi.add(rule27)

    pattern28 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('q', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda v, y: ZeroQ(-v + y)), CustomConstraint(lambda w, y: ZeroQ(-w + y)))
    rule28 = ReplacementRule(pattern28, lambda b, d, q, n, a, w, p, c, f, x, m, e, y, u, v : Module(List(Set(r, DerivativeDivides(y, u, x))), Condition(r*Subst(Int((a + b*x**n)**m*(c + d*x**n)**p*(e + f*x**n)**q, x), x, y), ~(FalseQ(r)))))
    rubi.add(rule28)

    pattern29 = Pattern(Integral(F_**v_*u_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)))
    rule29 = ReplacementRule(pattern29, lambda v, x, F, u : Module(List(Set(q, DerivativeDivides(v, u, x))), Condition(F**v*q/Log(F), ~(FalseQ(q)))))
    rubi.add(rule29)

    pattern30 = Pattern(Integral(F_**v_*u_*w_**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda v, w: ZeroQ(-v + w)))
    rule30 = ReplacementRule(pattern30, lambda w, x, m, u, v, F : Module(List(Set(q, DerivativeDivides(v, u, x))), Condition(q*Subst(Int(F**x*x**m, x), x, v), ~(FalseQ(q)))))
    rubi.add(rule30)

    pattern31 = Pattern(Integral(u_*(a_ + v_**WC('p', S(1))*w_**WC('p', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: IntegerQ(p)))
    rule31 = ReplacementRule(pattern31, lambda b, a, w, p, x, m, u, v : Module(List(Set(c, Simplify(u/(v*D(w, x) + w*D(v, x))))), Condition(c*Subst(Int((a + b*x**p)**m, x), x, v*w))))
    rubi.add(rule31)

    pattern32 = Pattern(Integral(u_*v_**WC('r', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda r, p, q: ZeroQ(p - q*(r + S(1)))), CustomConstraint(lambda r: NonzeroQ(r + S(1))), CustomConstraint(lambda r, p: IntegerQ(p/(r + S(1)))))
    rule32 = ReplacementRule(pattern32, lambda b, q, a, w, p, x, m, r, u, v : Module(List(Set(c, Simplify(u/(p*w*D(v, x) + q*v*D(w, x))))), Condition(c*p*Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w)/(r + S(1)))))
    rubi.add(rule32)

    pattern33 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda s, r, p, q: ZeroQ(p*(s + S(1)) - q*(r + S(1)))), CustomConstraint(lambda r: NonzeroQ(r + S(1))), CustomConstraint(lambda r, p: IntegerQ(p/(r + S(1)))))
    rule33 = ReplacementRule(pattern33, lambda b, q, a, w, p, x, s, r, m, u, v : Module(List(Set(c, Simplify(u/(p*w*D(v, x) + q*v*D(w, x))))), Condition(c*p*Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w**(s + S(1)))/(r + S(1)))))
    rubi.add(rule33)

    pattern34 = Pattern(Integral(u_*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda m, p, q: ZeroQ(p + q*(m*p + S(1)))), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m)))
    rule34 = ReplacementRule(pattern34, lambda b, q, a, w, p, x, m, u, v : Module(List(Set(c, Simplify(u/(p*w*D(v, x) - q*v*D(w, x))))), Condition(c*p*Subst(Int((a*x**p + b)**m, x), x, v*w**(m*q + S(1))))))
    rubi.add(rule34)

    pattern35 = Pattern(Integral(u_*v_**WC('r', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda m, r, p, q: ZeroQ(p + q*(m*p + r + S(1)))), CustomConstraint(lambda q: IntegerQ(q)), CustomConstraint(lambda m: IntegerQ(m)))
    rule35 = ReplacementRule(pattern35, lambda b, q, a, w, p, x, m, r, u, v : Module(List(Set(c, Simplify(u/(p*w*D(v, x) - q*v*D(w, x))))), Condition(-c*q*Subst(Int((a + b*x**q)**m, x), x, v**(m*p + r + S(1))*w))))
    rubi.add(rule35)

    pattern36 = Pattern(Integral(u_*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda s, m, p, q: ZeroQ(p*(s + S(1)) + q*(m*p + S(1)))), CustomConstraint(lambda s: NonzeroQ(s + S(1))), CustomConstraint(lambda s, q: IntegerQ(q/(s + S(1)))), CustomConstraint(lambda m: IntegerQ(m)))
    rule36 = ReplacementRule(pattern36, lambda b, q, a, w, p, x, s, m, u, v : Module(List(Set(c, Simplify(u/(p*w*D(v, x) - q*v*D(w, x))))), Condition(-c*q*Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + S(1))*w**(s + S(1)))/(s + S(1)))))
    rubi.add(rule36)

    pattern37 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda q, p, s, m, r: ZeroQ(p*(s + S(1)) + q*(m*p + r + S(1)))), CustomConstraint(lambda s: NonzeroQ(s + S(1))), CustomConstraint(lambda s, q: IntegerQ(q/(s + S(1)))), CustomConstraint(lambda m: IntegerQ(m)))
    rule37 = ReplacementRule(pattern37, lambda b, q, a, w, p, x, s, r, m, u, v : Module(List(Set(c, Simplify(u/(p*w*D(v, x) - q*v*D(w, x))))), Condition(-c*q*Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + r + S(1))*w**(s + S(1)))/(s + S(1)))))
    rubi.add(rule37)

    pattern38 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda x, m, u: FunctionOfQ(x**(m + S(1)), u, x)))
    rule38 = ReplacementRule(pattern38, lambda m, x, u : Subst(Int(SubstFor(x**(m + S(1)), u, x), x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule38)

    pattern39 = Pattern(Integral(u_, x_))
    rule39 = ReplacementRule(pattern39, lambda x, u : Module(List(Set(lst, SubstForFractionalPowerOfLinear(u, x))), Condition(Part(lst, S(2))*Part(lst, S(4))*Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(1/Part(lst, S(2)))), ~(FalseQ(lst)) & SubstForFractionalPowerQ(u, Part(lst, S(3)), x))))
    rubi.add(rule39)

    pattern40 = Pattern(Integral(u_, x_))
    rule40 = ReplacementRule(pattern40, lambda x, u : Module(List(Set(lst, SubstForFractionalPowerOfQuotientOfLinears(u, x))), Condition(Part(lst, S(2))*Part(lst, S(4))*Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(1/Part(lst, S(2)))), ~(FalseQ(lst)))))
    rubi.add(rule40)

    pattern41 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*z_**WC('q', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda v, x: NFreeQ(v, x)), CustomConstraint(lambda x, w: NFreeQ(w, x)), CustomConstraint(lambda x, z: NFreeQ(z, x)))
    rule41 = ReplacementRule(pattern41, lambda q, n, a, w, p, z, x, m, u, v : a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*z**(-q*FracPart(p))*(a*v**m*w**n*z**q)**FracPart(p)*Int(u*v**(m*p)*w**(n*p)*z**(p*q), x))
    rubi.add(rule41)

    pattern42 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda v, x: NFreeQ(v, x)), CustomConstraint(lambda x, w: NFreeQ(w, x)))
    rule42 = ReplacementRule(pattern42, lambda n, a, w, p, x, m, u, v : a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*(a*v**m*w**n)**FracPart(p)*Int(u*v**(m*p)*w**(n*p), x))
    rubi.add(rule42)

    pattern43 = Pattern(Integral((v_**WC('m', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda v, x: NFreeQ(v, x)), CustomConstraint(lambda m, a: ~(EqQ(a, S(1)) & EqQ(m, S(1)))), CustomConstraint(lambda v, m, x: ~(EqQ(m, S(1)) & EqQ(v, x))))
    rule43 = ReplacementRule(pattern43, lambda a, p, x, m, u, v : a**IntPart(p)*v**(-m*FracPart(p))*(a*v**m)**FracPart(p)*Int(u*v**(m*p), x))
    rubi.add(rule43)

    pattern44 = Pattern(Integral((x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda x, u: ~(RationalFunctionQ(u, x))))
    rule44 = ReplacementRule(pattern44, lambda b, n, a, p, u, x : FullSimplify(x**(-n/S(2))*Sqrt(a + b*x**n)/Sqrt(a*x**(-n) + b))*Int(u*x**(n*p)*(a*x**(-n) + b)**p, x))
    rubi.add(rule44)

    pattern45 = Pattern(Integral((v_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda v, x: BinomialQ(v, x)), CustomConstraint(lambda v, x: ~(LinearQ(v, x))))
    rule45 = ReplacementRule(pattern45, lambda b, n, a, p, x, u, v : v**(-n*FracPart(p))*(a + b*v**n)**FracPart(p)*(a*v**(-n) + b)**(-FracPart(p))*Int(u*v**(n*p)*(a*v**(-n) + b)**p, x))
    rubi.add(rule45)

    pattern46 = Pattern(Integral((v_**n_*x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda v, x: BinomialQ(v, x)))
    rule46 = ReplacementRule(pattern46, lambda b, n, a, p, x, m, u, v : v**(-n*FracPart(p))*(a + b*v**n*x**m)**FracPart(p)*(a*v**(-n) + b*x**m)**(-FracPart(p))*Int(u*v**(n*p)*(a*v**(-n) + b*x**m)**p, x))
    rubi.add(rule46)

    pattern47 = Pattern(Integral((x_**WC('r', S(1))*WC('a', S(1)) + x_**WC('s', S(1))*WC('b', S(1)))**m_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda s, r: PosQ(-r + s)))
    rule47 = ReplacementRule(pattern47, lambda b, a, s, r, m, u, x : With(List(Set(v, x**(-r*FracPart(m))*(a + b*x**(-r + s))**(-FracPart(m))*(a*x**r + b*x**s)**FracPart(m))), Condition(v*Int(u*x**(m*r)*(a + b*x**(-r + s))**m, x), ~(EqQ(Simplify(v), S(1))))))
    rubi.add(rule47)

    pattern48 = Pattern(Integral(u_/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule48 = ReplacementRule(pattern48, lambda b, n, a, u, x : With(List(Set(v, RationalFunctionExpand(u/(a + b*x**n), x))), Condition(Int(v, x), SumQ(v))))
    rubi.add(rule48)

    pattern49 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda b, a, c: ZeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda x, u: ~(AlgebraicFunctionQ(u, x))))
    rule49 = ReplacementRule(pattern49, lambda b, n, a, n2, p, c, u, x : S(4)**(-p)*c**(-p)*Int(u*(b + S(2)*c*x**n)**(S(2)*p), x))
    rubi.add(rule49)

    pattern50 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda b, a, c: ZeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda x, u: ~(AlgebraicFunctionQ(u, x))))
    rule50 = ReplacementRule(pattern50, lambda b, n, a, n2, p, c, u, x : (b + S(2)*c*x**n)**(-S(2)*p)*(a + b*x**n + c*x**(S(2)*n))**p*Int(u*(b + S(2)*c*x**n)**(S(2)*p), x))
    rubi.add(rule50)

    pattern51 = Pattern(Integral(u_/(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule51 = ReplacementRule(pattern51, lambda b, n, a, n2, c, u, x : With(List(Set(v, RationalFunctionExpand(u/(a + b*x**n + c*x**(S(2)*n)), x))), Condition(Int(v, x), SumQ(v))))
    rubi.add(rule51)

    pattern52 = Pattern(Integral(WC('u', S(1))/(x_**WC('m', S(1))*WC('a', S(1)) + sqrt(x_**n_*WC('c', S(1)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule52 = ReplacementRule(pattern52, lambda b, n, a, c, m, u, x : Int(u*(a*x**m - b*Sqrt(c*x**n))/(a**S(2)*x**(S(2)*m) - b**S(2)*c*x**n), x))
    rubi.add(rule52)

    pattern53 = Pattern(Integral(u_, x_))
    rule53 = ReplacementRule(pattern53, lambda x, u : Module(List(Set(lst, FunctionOfLinear(u, x))), Condition(Dist(1/Part(lst, S(3)), Subst(Int(Part(lst, S(1)), x), x, x*Part(lst, S(3)) + Part(lst, S(2))), x), ~(FalseQ(lst)))))
    rubi.add(rule53)

    pattern54 = Pattern(Integral(u_/x_, x_), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda x, u: ~(RationalFunctionQ(u, x))))
    rule54 = ReplacementRule(pattern54, lambda x, u : Module(List(Set(lst, PowerVariableExpn(u, S(0), x))), Condition(Subst(Int(NormalizeIntegrand(Simplify(Part(lst, S(1))/x), x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2)))/Part(lst, S(2)), ~(FalseQ(lst)) & NonzeroQ(Part(lst, S(2))))))
    rubi.add(rule54)

    pattern55 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Unequal(m, S(-1))), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda m, x, u: Greater(m, S(0)) | ~(AlgebraicFunctionQ(u, x))))
    rule55 = ReplacementRule(pattern55, lambda m, x, u : Module(List(Set(lst, PowerVariableExpn(u, m + S(1), x))), Condition(Subst(Int(NormalizeIntegrand(Simplify(Part(lst, S(1))/x), x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2)))/Part(lst, S(2)), ~(FalseQ(lst)) & NonzeroQ(-m + Part(lst, S(2)) + S(-1)))))
    rubi.add(rule55)

    pattern56 = Pattern(Integral(u_*x_**m_, x_), CustomConstraint(lambda m: FractionQ(m)))
    rule56 = ReplacementRule(pattern56, lambda x, m, u : Module(List(Set(k, Denominator(m))), k*Subst(Int(x**(k*(m + S(1)) + S(-1))*ReplaceAll(u, Rule(x, x**k)), x), x, x**(1/k))))
    rubi.add(rule56)

    pattern57 = Pattern(Integral(u_, x_), CustomConstraint(lambda x, u: EulerIntegrandQ(u, x)))
    rule57 = ReplacementRule(pattern57, lambda x, u : Module(List(Set(lst, FunctionOfSquareRootOfQuadratic(u, x))), Condition(S(2)*Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(2))), ~(FalseQ(lst)))))
    rubi.add(rule57)

    pattern58 = Pattern(Integral(1/(a_ + v_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule58 = ReplacementRule(pattern58, lambda v, b, a, x : Int(Together(1/(-v/Rt(-a/b, S(2)) + S(1))), x)/(S(2)*a) + Int(Together(1/(v/Rt(-a/b, S(2)) + S(1))), x)/(S(2)*a))
    rubi.add(rule58)

    pattern59 = Pattern(Integral(1/(a_ + v_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: EvenQ(n)), CustomConstraint(lambda n: Greater(n, S(2))))
    rule59 = ReplacementRule(pattern59, lambda b, n, a, x, v : Dist(S(2)/(a*n), Sum(Int(Together(1/(S(1) - (S(-1))**(-S(4)*k/n)*v**S(2)/Rt(-a/b, n/S(2)))), x), List(k, S(1), n/S(2))), x))
    rubi.add(rule59)

    pattern60 = Pattern(Integral(1/(a_ + v_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule60 = ReplacementRule(pattern60, lambda b, n, a, x, v : Dist(S(1)/(a*n), Sum(Int(Together(1/(S(1) - (S(-1))**(-S(2)*k/n)*v/Rt(-a/b, n))), x), List(k, S(1), n)), x))
    rubi.add(rule60)

    pattern61 = Pattern(Integral(v_/(a_ + u_**WC('n', S(1))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda v, x, u: PolynomialInQ(v, u, x)))
    rule61 = ReplacementRule(pattern61, lambda b, n, a, x, u, v : Int(ReplaceAll(ExpandIntegrand(PolynomialInSubst(v, u, x)/(a + b*x**n), x), Rule(x, u)), x))
    rubi.add(rule61)

    pattern62 = Pattern(Integral(u_, x_))
    rule62 = ReplacementRule(pattern62, lambda x, u : Module(List(Set(v, NormalizeIntegrand(u, x))), Condition(Int(v, x), UnsameQ(v, u))))
    rubi.add(rule62)

    pattern63 = Pattern(Integral(u_, x_))
    rule63 = ReplacementRule(pattern63, lambda x, u : Module(List(Set(v, ExpandIntegrand(u, x))), Condition(Int(v, x), SumQ(v))))
    rubi.add(rule63)

    pattern64 = Pattern(Integral((x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda a, d: ZeroQ(a + d)), CustomConstraint(lambda b, c: ZeroQ(b + c)), CustomConstraint(lambda m, n: ZeroQ(m + n)), CustomConstraint(lambda p, q: ZeroQ(p + q)))
    rule64 = ReplacementRule(pattern64, lambda b, d, q, n, a, p, c, m, u, x : x**(-m*p)*(a + b*x**m)**p*(c + d*x**n)**q*Int(u*x**(m*p), x))
    rubi.add(rule64)

    pattern65 = Pattern(Integral(u_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda b, a, c: ZeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule65 = ReplacementRule(pattern65, lambda b, n, n2, a, p, c, u, x : (S(4)*c)**(-p + S(1)/2)*Int(u*(b + S(2)*c*x**n)**(S(2)*p), x)*Sqrt(a + b*x**n + c*x**(S(2)*n))/(b + S(2)*c*x**n))
    rubi.add(rule65)

    pattern66 = Pattern(Integral(u_, x_))
    rule66 = ReplacementRule(pattern66, lambda x, u : Module(List(Set(lst, SubstForFractionalPowerOfLinear(u, x))), Condition(Part(lst, S(2))*Part(lst, S(4))*Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(1/Part(lst, S(2)))), ~(FalseQ(lst)))))
    rubi.add(rule66)

    pattern67 = Pattern(Integral(u_, x_))
    rule67 = ReplacementRule(pattern67, lambda x, u : Int(u, x))
    rubi.add(rule67)

    return rubi
