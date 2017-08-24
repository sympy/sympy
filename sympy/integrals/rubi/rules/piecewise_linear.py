
from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (D, Module, Int, Set, With, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcCsch, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct, SplitSum, Complex, UnsameQ, _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify, _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum, _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux, TrigSimplifyAux)
    from sympy import Integral, S, sqrt
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)


    A_, B_, C_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, z_ = [WC(i) for i in 'ABCabcdefghijklmnpqrtuvswxz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, jn_, mn_, non2_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'n1', 'n2', 'n3',' Pq', 'Pm', ' Px', 'Qm', 'Qr', 'jn', 'mn', 'non2']]
    c, p, q, r, s, mn, gcd, P, Q = symbols('c p q r s mn gcd P Q')

    _UseGamma = False

def piecewise_linear(rubi):

    pattern1 = Pattern(Integral(u_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda u, x: PiecewiseLinearQ(u, x)))
    rule1 = ReplacementRule(pattern1, lambda m, u, x : Module(List(Set(c, Simplify(D(u, x)))), Subst(Int(x**m, x), x, u)/c))
    rubi.add(rule1)

    pattern2 = Pattern(Integral(v_/u_, x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)))
    rule2 = ReplacementRule(pattern2, lambda u, x, v : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(b*x/a - (-a*v + b*u)*Int(1/u, x)/a, NonzeroQ(-a*v + b*u))))
    rubi.add(rule2)

    pattern3 = Pattern(Integral(v_**n_/u_, x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda n: Unequal(n, S(1))))
    rule3 = ReplacementRule(pattern3, lambda u, x, v, n : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(-(-a*v + b*u)*Int(v**(n + S(-1))/u, x)/a + v**n/(a*n), NonzeroQ(-a*v + b*u))))
    rubi.add(rule3)

    pattern4 = Pattern(Integral(S(1)/(u_*v_), x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)))
    rule4 = ReplacementRule(pattern4, lambda u, x, v : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(-a*Int(1/u, x)/(-a*v + b*u) + b*Int(1/v, x)/(-a*v + b*u), NonzeroQ(-a*v + b*u))))
    rubi.add(rule4)

    pattern5 = Pattern(Integral(S(1)/(u_*sqrt(v_)), x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)))
    rule5 = ReplacementRule(pattern5, lambda u, x, v : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(S(2)*ArcTan(Sqrt(v)/Rt((-a*v + b*u)/a, S(2)))/(a*Rt((-a*v + b*u)/a, S(2))), NonzeroQ(-a*v + b*u) & PosQ((-a*v + b*u)/a))))
    rubi.add(rule5)

    pattern6 = Pattern(Integral(S(1)/(u_*sqrt(v_)), x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)))
    rule6 = ReplacementRule(pattern6, lambda u, x, v : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(-S(2)*ArcTanh(Sqrt(v)/Rt((a*v - b*u)/a, S(2)))/(a*Rt((a*v - b*u)/a, S(2))), NonzeroQ(-a*v + b*u) & NegQ((-a*v + b*u)/a))))
    rubi.add(rule6)

    pattern7 = Pattern(Integral(v_**n_/u_, x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule7 = ReplacementRule(pattern7, lambda u, x, v, n : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(-a*Int(v**(n + S(1))/u, x)/(-a*v + b*u) + v**(n + S(1))/((n + S(1))*(-a*v + b*u)), NonzeroQ(-a*v + b*u))))
    rubi.add(rule7)

    pattern8 = Pattern(Integral(v_**n_/u_, x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda n: ~(IntegerQ(n))))
    rule8 = ReplacementRule(pattern8, lambda u, x, v, n : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(v**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), -a*v/(-a*v + b*u))/((n + S(1))*(-a*v + b*u)), NonzeroQ(-a*v + b*u))))
    rubi.add(rule8)

    pattern9 = Pattern(Integral(S(1)/(sqrt(u_)*sqrt(v_)), x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)))
    rule9 = ReplacementRule(pattern9, lambda u, x, v : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(S(2)*ArcTanh(Rt(a*b, S(2))*Sqrt(u)/(a*Sqrt(v)))/Rt(a*b, S(2)), PosQ(a*b) & NonzeroQ(-a*v + b*u))))
    rubi.add(rule9)

    pattern10 = Pattern(Integral(S(1)/(sqrt(u_)*sqrt(v_)), x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)))
    rule10 = ReplacementRule(pattern10, lambda u, x, v : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(S(2)*ArcTan(Rt(-a*b, S(2))*Sqrt(u)/(a*Sqrt(v)))/Rt(-a*b, S(2)), NegQ(a*b) & NonzeroQ(-a*v + b*u))))
    rubi.add(rule10)

    pattern11 = Pattern(Integral(u_**m_*v_**n_, x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda m, n: ZeroQ(m + n + S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule11 = ReplacementRule(pattern11, lambda x, v, m, u, n : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(-u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u)), NonzeroQ(-a*v + b*u))))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(u_**m_*v_**WC('n', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda m, n: (NegativeIntegerQ(m) & ~(IntegerQ(n))) | (PositiveIntegerQ(n) & ~(IntegerQ(m))) | (LessEqual(n, m) & PositiveIntegerQ(n, m)) | (Greater(n, S(0)) & Less(m, S(-1)) & RationalQ(m, n) & ~(IntegerQ(m + n) & Less(m + n + S(2), S(0)) & (FractionQ(m) | GreaterEqual(m + S(2)*n + S(1), S(0)))))))
    rule12 = ReplacementRule(pattern12, lambda x, v, m, u, n : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(-b*n*Int(u**(m + S(1))*v**(n + S(-1)), x)/(a*(m + S(1))) + u**(m + S(1))*v**n/(a*(m + S(1))), NonzeroQ(-a*v + b*u))))
    rubi.add(rule12)

    pattern13 = Pattern(Integral(u_**m_*v_**WC('n', S(1)), x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda m, n: NonzeroQ(m + n + S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m, n: NonzeroQ(m + n + S(1))), CustomConstraint(lambda m, n: ~(PositiveIntegerQ(m) & (~(IntegerQ(n)) | Less(S(0), m, n)))), CustomConstraint(lambda m, n: ~(IntegerQ(m + n) & Less(m + n + S(2), S(0)))))
    rule13 = ReplacementRule(pattern13, lambda x, v, m, u, n : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(-n*(-a*v + b*u)*Int(u**m*v**(n + S(-1)), x)/(a*(m + n + S(1))) + u**(m + S(1))*v**n/(a*(m + n + S(1))), NonzeroQ(-a*v + b*u))))
    rubi.add(rule13)

    pattern14 = Pattern(Integral(u_**m_*v_**n_, x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda m, n: NonzeroQ(m + n + S(1))), CustomConstraint(lambda n: ~(RationalQ(n))), CustomConstraint(lambda n: SumSimplerQ(n, S(-1))))
    rule14 = ReplacementRule(pattern14, lambda x, v, m, u, n : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(-n*(-a*v + b*u)*Int(u**m*v**Simplify(n + S(-1)), x)/(a*(m + n + S(1))) + u**(m + S(1))*v**n/(a*(m + n + S(1))), NonzeroQ(-a*v + b*u))))
    rubi.add(rule14)

    pattern15 = Pattern(Integral(u_**m_*v_**n_, x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda m, n: NonzeroQ(m + n + S(2))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule15 = ReplacementRule(pattern15, lambda x, v, m, u, n : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(b*(m + n + S(2))*Int(u**(m + S(1))*v**n, x)/((m + S(1))*(-a*v + b*u)) - u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u)), NonzeroQ(-a*v + b*u))))
    rubi.add(rule15)

    pattern16 = Pattern(Integral(u_**m_*v_**n_, x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda m: ~(RationalQ(m))), CustomConstraint(lambda m: SumSimplerQ(m, S(1))))
    rule16 = ReplacementRule(pattern16, lambda x, v, m, u, n : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(b*(m + n + S(2))*Int(u**Simplify(m + S(1))*v**n, x)/((m + S(1))*(-a*v + b*u)) - u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u)), NonzeroQ(-a*v + b*u))))
    rubi.add(rule16)

    pattern17 = Pattern(Integral(u_**m_*v_**n_, x_), CustomConstraint(lambda u, x, v: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda m: ~(IntegerQ(m))), CustomConstraint(lambda n: ~(IntegerQ(n))))
    rule17 = ReplacementRule(pattern17, lambda x, v, m, u, n : Module(List(Set(a, Simplify(D(u, x))), Set(b, Simplify(D(v, x)))), Condition(u**m*v**(n + S(1))*(b*u/(-a*v + b*u))**(-m)*Hypergeometric2F1(-m, n + S(1), n + S(2), -a*v/(-a*v + b*u))/(b*(n + S(1))), NonzeroQ(-a*v + b*u))))
    rubi.add(rule17)

    pattern18 = Pattern(Integral(u_**WC('n', S(1))*Log(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda u, x: PiecewiseLinearQ(u, x)), CustomConstraint(lambda u, x: ~(LinearQ(u, x))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule18 = ReplacementRule(pattern18, lambda x, b, a, u, n : Module(List(Set(c, Simplify(D(u, x)))), -Int(u**n, x) - c*n*Int(u**(n + S(-1))*(a + b*x)*Log(a + b*x), x)/b + u**n*(a + b*x)*Log(a + b*x)/b))
    rubi.add(rule18)

    pattern19 = Pattern(Integral(u_**WC('n', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*Log(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda u, x: PiecewiseLinearQ(u, x)), CustomConstraint(lambda u, x: ~(LinearQ(u, x))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule19 = ReplacementRule(pattern19, lambda x, b, a, m, n, u : Module(List(Set(c, Simplify(D(u, x)))), -Int(u**n*(a + b*x)**m, x)/(m + S(1)) - c*n*Int(u**(n + S(-1))*(a + b*x)**(m + S(1))*Log(a + b*x), x)/(b*(m + S(1))) + u**n*(a + b*x)**(m + S(1))*Log(a + b*x)/(b*(m + S(1)))))
    rubi.add(rule19)

    return rubi
