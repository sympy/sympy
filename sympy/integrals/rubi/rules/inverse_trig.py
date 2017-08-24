
from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (ArcSec, ArcSech, PolyLog, Int, Set, With, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcCsch, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct, SplitSum, Complex, UnsameQ, _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify, _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum, _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux, TrigSimplifyAux)
    from sympy import Integral, S, sqrt
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)


    A_, B_, C_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, z_ = [WC(i) for i in 'ABCabcdefghijklmnpqrtuvswxz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, jn_, mn_, non2_, RFx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'n1', 'n2', 'n3',' Pq', 'Pm', ' Px', 'Qm', 'Qr', 'jn', 'mn', 'non2', 'RFx']]
    p, q, r, s, mn, gcd, P, Q = symbols('p q r s mn gcd P Q')

    _UseGamma = False

def inverse_trig(rubi):

    pattern1 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule1 = ReplacementRule(pattern1, lambda n, c, a, x, b : -b*c*n*Int(x*(a + b*ArcSin(c*x))**(n + S(-1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x) + x*(a + b*ArcSin(c*x))**n)
    rubi.add(rule1)

    pattern2 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule2 = ReplacementRule(pattern2, lambda n, c, a, x, b : b*c*n*Int(x*(a + b*ArcCos(c*x))**(n + S(-1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x) + x*(a + b*ArcCos(c*x))**n)
    rubi.add(rule2)

    pattern3 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule3 = ReplacementRule(pattern3, lambda n, c, a, x, b : c*Int(x*(a + b*ArcSin(c*x))**(n + S(1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*(n + S(1))) + (a + b*ArcSin(c*x))**(n + S(1))*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule3)

    pattern4 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule4 = ReplacementRule(pattern4, lambda n, c, a, x, b : -c*Int(x*(a + b*ArcCos(c*x))**(n + S(1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*(n + S(1))) - (a + b*ArcCos(c*x))**(n + S(1))*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule4)

    pattern5 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule5 = ReplacementRule(pattern5, lambda n, c, a, x, b : Subst(Int(x**n*Cos(a/b - x/b), x), x, a + b*ArcSin(c*x))/(b*c))
    rubi.add(rule5)

    pattern6 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule6 = ReplacementRule(pattern6, lambda n, c, a, x, b : Subst(Int(x**n*Sin(a/b - x/b), x), x, a + b*ArcCos(c*x))/(b*c))
    rubi.add(rule6)

    pattern7 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule7 = ReplacementRule(pattern7, lambda n, c, a, x, b : Subst(Int((a + b*x)**n/Tan(x), x), x, ArcSin(c*x)))
    rubi.add(rule7)

    pattern8 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule8 = ReplacementRule(pattern8, lambda n, c, a, x, b : -Subst(Int((a + b*x)**n/Cot(x), x), x, ArcCos(c*x)))
    rubi.add(rule8)

    pattern9 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule9 = ReplacementRule(pattern9, lambda n, c, d, a, m, x, b : -b*c*n*Int((d*x)**(m + S(1))*(a + b*ArcSin(c*x))**(n + S(-1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(d*(m + S(1))) + (d*x)**(m + S(1))*(a + b*ArcSin(c*x))**n/(d*(m + S(1))))
    rubi.add(rule9)

    pattern10 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule10 = ReplacementRule(pattern10, lambda n, c, d, a, m, x, b : b*c*n*Int((d*x)**(m + S(1))*(a + b*ArcCos(c*x))**(n + S(-1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(d*(m + S(1))) + (d*x)**(m + S(1))*(a + b*ArcCos(c*x))**n/(d*(m + S(1))))
    rubi.add(rule10)

    pattern11 = Pattern(Integral(x_**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule11 = ReplacementRule(pattern11, lambda n, c, a, m, x, b : -b*c*n*Int(x**(m + S(1))*(a + b*ArcSin(c*x))**(n + S(-1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcSin(c*x))**n/(m + S(1)))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(x_**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule12 = ReplacementRule(pattern12, lambda n, c, a, m, x, b : b*c*n*Int(x**(m + S(1))*(a + b*ArcCos(c*x))**(n + S(-1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcCos(c*x))**n/(m + S(1)))
    rubi.add(rule12)

    pattern13 = Pattern(Integral(x_**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Inequality(S(-2), LessEqual, n, Less, S(-1))))
    rule13 = ReplacementRule(pattern13, lambda n, c, a, m, x, b : -c**(-m + S(-1))*Subst(Int(ExpandTrigReduce((a + b*x)**(n + S(1)), (m - (m + S(1))*Sin(x)**S(2))*Sin(x)**(m + S(-1)), x), x), x, ArcSin(c*x))/(b*(n + S(1))) + x**m*(a + b*ArcSin(c*x))**(n + S(1))*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule13)

    pattern14 = Pattern(Integral(x_**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Inequality(S(-2), LessEqual, n, Less, S(-1))))
    rule14 = ReplacementRule(pattern14, lambda n, c, a, m, x, b : -c**(-m + S(-1))*Subst(Int(ExpandTrigReduce((a + b*x)**(n + S(1)), (m - (m + S(1))*Cos(x)**S(2))*Cos(x)**(m + S(-1)), x), x), x, ArcCos(c*x))/(b*(n + S(1))) - x**m*(a + b*ArcCos(c*x))**(n + S(1))*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule14)

    pattern15 = Pattern(Integral(x_**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-2))))
    rule15 = ReplacementRule(pattern15, lambda n, c, a, m, x, b : c*(m + S(1))*Int(x**(m + S(1))*(a + b*ArcSin(c*x))**(n + S(1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*(n + S(1))) - m*Int(x**(m + S(-1))*(a + b*ArcSin(c*x))**(n + S(1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*c*(n + S(1))) + x**m*(a + b*ArcSin(c*x))**(n + S(1))*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule15)

    pattern16 = Pattern(Integral(x_**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-2))))
    rule16 = ReplacementRule(pattern16, lambda n, c, a, m, x, b : -c*(m + S(1))*Int(x**(m + S(1))*(a + b*ArcCos(c*x))**(n + S(1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*(n + S(1))) + m*Int(x**(m + S(-1))*(a + b*ArcCos(c*x))**(n + S(1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*c*(n + S(1))) - x**m*(a + b*ArcCos(c*x))**(n + S(1))*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule16)

    pattern17 = Pattern(Integral(x_**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule17 = ReplacementRule(pattern17, lambda n, c, a, m, x, b : c**(-m + S(-1))*Subst(Int((a + b*x)**n*Cos(x)*Sin(x)**m, x), x, ArcSin(c*x)))
    rubi.add(rule17)

    pattern18 = Pattern(Integral(x_**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule18 = ReplacementRule(pattern18, lambda n, c, a, m, x, b : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*Cos(x)**m*Sin(x), x), x, ArcCos(c*x)))
    rubi.add(rule18)

    pattern19 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule19 = ReplacementRule(pattern19, lambda n, c, d, a, m, x, b : Int((d*x)**m*(a + b*ArcSin(c*x))**n, x))
    rubi.add(rule19)

    pattern20 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule20 = ReplacementRule(pattern20, lambda n, c, d, a, m, x, b : Int((d*x)**m*(a + b*ArcCos(c*x))**n, x))
    rubi.add(rule20)

    pattern21 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule21 = ReplacementRule(pattern21, lambda c, d, a, x, e, b : Log(a + b*ArcSin(c*x))/(b*c*Sqrt(d)))
    rubi.add(rule21)

    pattern22 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule22 = ReplacementRule(pattern22, lambda c, d, a, x, e, b : -Log(a + b*ArcCos(c*x))/(b*c*Sqrt(d)))
    rubi.add(rule22)

    pattern23 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule23 = ReplacementRule(pattern23, lambda n, c, d, a, x, e, b : (a + b*ArcSin(c*x))**(n + S(1))/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule23)

    pattern24 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule24 = ReplacementRule(pattern24, lambda n, c, d, a, x, e, b : -(a + b*ArcCos(c*x))**(n + S(1))/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule24)

    pattern25 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: ~(PositiveQ(d))))
    rule25 = ReplacementRule(pattern25, lambda n, c, d, a, x, e, b : Int((a + b*ArcSin(c*x))**n/Sqrt(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(-c**S(2)*x**S(2) + S(1))/Sqrt(d + e*x**S(2)))
    rubi.add(rule25)

    pattern26 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: ~(PositiveQ(d))))
    rule26 = ReplacementRule(pattern26, lambda n, c, d, a, x, e, b : Int((a + b*ArcCos(c*x))**n/Sqrt(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(-c**S(2)*x**S(2) + S(1))/Sqrt(d + e*x**S(2)))
    rubi.add(rule26)

    pattern27 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule27 = ReplacementRule(pattern27, lambda p, c, d, a, x, e, b : With(List(Set(u, IntHide((d + e*x**S(2))**p, x))), -b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcSin(c*x), u, x)))
    rubi.add(rule27)

    pattern28 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule28 = ReplacementRule(pattern28, lambda p, c, d, a, x, e, b : With(List(Set(u, IntHide((d + e*x**S(2))**p, x))), b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcCos(c*x), u, x)))
    rubi.add(rule28)

    pattern29 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule29 = ReplacementRule(pattern29, lambda n, c, d, a, x, e, b : -b*c*n*Int(x*(a + b*ArcSin(c*x))**(n + S(-1)), x)*Sqrt(d + e*x**S(2))/(S(2)*Sqrt(-c**S(2)*x**S(2) + S(1))) + x*(a + b*ArcSin(c*x))**n*Sqrt(d + e*x**S(2))/S(2) + Int((a + b*ArcSin(c*x))**n/Sqrt(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(d + e*x**S(2))/(S(2)*Sqrt(-c**S(2)*x**S(2) + S(1))))
    rubi.add(rule29)

    pattern30 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule30 = ReplacementRule(pattern30, lambda n, c, d, a, x, e, b : b*c*n*Int(x*(a + b*ArcCos(c*x))**(n + S(-1)), x)*Sqrt(d + e*x**S(2))/(S(2)*Sqrt(-c**S(2)*x**S(2) + S(1))) + x*(a + b*ArcCos(c*x))**n*Sqrt(d + e*x**S(2))/S(2) + Int((a + b*ArcCos(c*x))**n/Sqrt(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(d + e*x**S(2))/(S(2)*Sqrt(-c**S(2)*x**S(2) + S(1))))
    rubi.add(rule30)

    pattern31 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))))
    rule31 = ReplacementRule(pattern31, lambda p, n, c, d, a, x, e, b : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*ArcSin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(S(2)*p + S(1)) + S(2)*d*p*Int((a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule31)

    pattern32 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))))
    rule32 = ReplacementRule(pattern32, lambda p, n, c, d, a, x, e, b : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*ArcCos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(S(2)*p + S(1)) + S(2)*d*p*Int((a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule32)

    pattern33 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda d: PositiveQ(d)))
    rule33 = ReplacementRule(pattern33, lambda n, c, d, a, x, e, b : -b*c*n*Int(x*(a + b*ArcSin(c*x))**(n + S(-1))/(d + e*x**S(2)), x)/Sqrt(d) + x*(a + b*ArcSin(c*x))**n/(d*Sqrt(d + e*x**S(2))))
    rubi.add(rule33)

    pattern34 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda d: PositiveQ(d)))
    rule34 = ReplacementRule(pattern34, lambda n, c, d, a, x, e, b : b*c*n*Int(x*(a + b*ArcCos(c*x))**(n + S(-1))/(d + e*x**S(2)), x)/Sqrt(d) + x*(a + b*ArcCos(c*x))**n/(d*Sqrt(d + e*x**S(2))))
    rubi.add(rule34)

    pattern35 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule35 = ReplacementRule(pattern35, lambda n, c, d, a, x, e, b : -b*c*n*Int(x*(a + b*ArcSin(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(-c**S(2)*x**S(2) + S(1))/(d*Sqrt(d + e*x**S(2))) + x*(a + b*ArcSin(c*x))**n/(d*Sqrt(d + e*x**S(2))))
    rubi.add(rule35)

    pattern36 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule36 = ReplacementRule(pattern36, lambda n, c, d, a, x, e, b : b*c*n*Int(x*(a + b*ArcCos(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(-c**S(2)*x**S(2) + S(1))/(d*Sqrt(d + e*x**S(2))) + x*(a + b*ArcCos(c*x))**n/(d*Sqrt(d + e*x**S(2))))
    rubi.add(rule36)

    pattern37 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule37 = ReplacementRule(pattern37, lambda p, n, c, d, a, x, e, b : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*ArcSin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*(p + S(1))) - x*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule37)

    pattern38 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule38 = ReplacementRule(pattern38, lambda p, n, c, d, a, x, e, b : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*ArcCos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*(p + S(1))) - x*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule38)

    pattern39 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule39 = ReplacementRule(pattern39, lambda n, c, d, a, x, e, b : Subst(Int((a + b*x)**n*Sec(x), x), x, ArcSin(c*x))/(c*d))
    rubi.add(rule39)

    pattern40 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule40 = ReplacementRule(pattern40, lambda n, c, d, a, x, e, b : -Subst(Int((a + b*x)**n*Csc(x), x), x, ArcCos(c*x))/(c*d))
    rubi.add(rule40)

    pattern41 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule41 = ReplacementRule(pattern41, lambda p, n, c, d, a, x, e, b : c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(S(2)*p + S(1))*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*ArcSin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*(n + S(1))) + (a + b*ArcSin(c*x))**(n + S(1))*(d + e*x**S(2))**p*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule41)

    pattern42 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule42 = ReplacementRule(pattern42, lambda p, n, c, d, a, x, e, b : -c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(S(2)*p + S(1))*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*ArcCos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*(n + S(1))) - (a + b*ArcCos(c*x))**(n + S(1))*(d + e*x**S(2))**p*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule42)

    pattern43 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)), CustomConstraint(lambda d, p: IntegerQ(p) | PositiveQ(d)))
    rule43 = ReplacementRule(pattern43, lambda p, n, c, d, a, x, e, b : d**p*Subst(Int((a + b*x)**n*Cos(x)**(S(2)*p + S(1)), x), x, ArcSin(c*x))/c)
    rubi.add(rule43)

    pattern44 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)), CustomConstraint(lambda d, p: IntegerQ(p) | PositiveQ(d)))
    rule44 = ReplacementRule(pattern44, lambda p, n, c, d, a, x, e, b : -d**p*Subst(Int((a + b*x)**n*Sin(x)**(S(2)*p + S(1)), x), x, ArcCos(c*x))/c)
    rubi.add(rule44)

    pattern45 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)), CustomConstraint(lambda d, p: ~(IntegerQ(p) | PositiveQ(d))))
    rule45 = ReplacementRule(pattern45, lambda p, n, c, d, a, x, e, b : d**(p + S(-1)/2)*Int((a + b*ArcSin(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x)*Sqrt(d + e*x**S(2))/Sqrt(-c**S(2)*x**S(2) + S(1)))
    rubi.add(rule45)

    pattern46 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)), CustomConstraint(lambda d, p: ~(IntegerQ(p) | PositiveQ(d))))
    rule46 = ReplacementRule(pattern46, lambda p, n, c, d, a, x, e, b : d**(p + S(-1)/2)*Int((a + b*ArcCos(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x)*Sqrt(d + e*x**S(2))/Sqrt(-c**S(2)*x**S(2) + S(1)))
    rubi.add(rule46)

    pattern47 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)))
    rule47 = ReplacementRule(pattern47, lambda p, c, d, a, x, e, b : With(List(Set(u, IntHide((d + e*x**S(2))**p, x))), -b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcSin(c*x), u, x)))
    rubi.add(rule47)

    pattern48 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)))
    rule48 = ReplacementRule(pattern48, lambda p, c, d, a, x, e, b : With(List(Set(u, IntHide((d + e*x**S(2))**p, x))), b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcCos(c*x), u, x)))
    rubi.add(rule48)

    pattern49 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p, n: PositiveIntegerQ(n) | Greater(p, S(0))))
    rule49 = ReplacementRule(pattern49, lambda p, n, c, d, a, x, e, b : Int(ExpandIntegrand((a + b*ArcSin(c*x))**n, (d + e*x**S(2))**p, x), x))
    rubi.add(rule49)

    pattern50 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p, n: PositiveIntegerQ(n) | Greater(p, S(0))))
    rule50 = ReplacementRule(pattern50, lambda p, n, c, d, a, x, e, b : Int(ExpandIntegrand((a + b*ArcCos(c*x))**n, (d + e*x**S(2))**p, x), x))
    rubi.add(rule50)

    pattern51 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule51 = ReplacementRule(pattern51, lambda p, n, c, d, a, x, e, b : Int((a + b*ArcSin(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule51)

    pattern52 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule52 = ReplacementRule(pattern52, lambda p, n, c, d, a, x, e, b : Int((a + b*ArcCos(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule52)

    pattern53 = Pattern(Integral((d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, f, g, e: ZeroQ(d*g + e*f)), CustomConstraint(lambda f, g, c: ZeroQ(c**S(2)*f**S(2) - g**S(2))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule53 = ReplacementRule(pattern53, lambda p, n, c, d, a, f, x, e, b, g : (d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p))*Int((a + b*ArcSin(c*x))**n*(d*f + e*g*x**S(2))**p, x))
    rubi.add(rule53)

    pattern54 = Pattern(Integral((d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, f, g, e: ZeroQ(d*g + e*f)), CustomConstraint(lambda f, g, c: ZeroQ(c**S(2)*f**S(2) - g**S(2))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule54 = ReplacementRule(pattern54, lambda p, n, c, d, a, f, x, e, b, g : (d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p))*Int((a + b*ArcCos(c*x))**n*(d*f + e*g*x**S(2))**p, x))
    rubi.add(rule54)

    pattern55 = Pattern(Integral(x_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule55 = ReplacementRule(pattern55, lambda n, c, d, a, x, e, b : -Subst(Int((a + b*x)**n*Tan(x), x), x, ArcSin(c*x))/e)
    rubi.add(rule55)

    pattern56 = Pattern(Integral(x_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule56 = ReplacementRule(pattern56, lambda n, c, d, a, x, e, b : Subst(Int((a + b*x)**n*Cot(x), x), x, ArcCos(c*x))/e)
    rubi.add(rule56)

    pattern57 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule57 = ReplacementRule(pattern57, lambda p, n, c, d, a, x, e, b : b*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*ArcSin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) + (a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule57)

    pattern58 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule58 = ReplacementRule(pattern58, lambda p, n, c, d, a, x, e, b : -b*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*ArcCos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) + (a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule58)

    pattern59 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule59 = ReplacementRule(pattern59, lambda n, c, d, a, x, e, b : Subst(Int((a + b*x)**n/(Cos(x)*Sin(x)), x), x, ArcSin(c*x))/d)
    rubi.add(rule59)

    pattern60 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule60 = ReplacementRule(pattern60, lambda n, c, d, a, x, e, b : -Subst(Int((a + b*x)**n/(Cos(x)*Sin(x)), x), x, ArcCos(c*x))/d)
    rubi.add(rule60)

    pattern61 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule61 = ReplacementRule(pattern61, lambda p, n, c, d, a, f, x, m, e, b : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*ArcSin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + (f*x)**(m + S(1))*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule61)

    pattern62 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule62 = ReplacementRule(pattern62, lambda p, n, c, d, a, f, x, m, e, b : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*ArcCos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + (f*x)**(m + S(1))*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule62)

    pattern63 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule63 = ReplacementRule(pattern63, lambda p, c, d, a, x, e, b : -b*c*d**p*Int((-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(S(2)*p) + d*Int((a + b*ArcSin(c*x))*(d + e*x**S(2))**(p + S(-1))/x, x) + (a + b*ArcSin(c*x))*(d + e*x**S(2))**p/(S(2)*p))
    rubi.add(rule63)

    pattern64 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule64 = ReplacementRule(pattern64, lambda p, c, d, a, x, e, b : b*c*d**p*Int((-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(S(2)*p) + d*Int((a + b*ArcCos(c*x))*(d + e*x**S(2))**(p + S(-1))/x, x) + (a + b*ArcCos(c*x))*(d + e*x**S(2))**p/(S(2)*p))
    rubi.add(rule64)

    pattern65 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2) + S(1)/2)))
    rule65 = ReplacementRule(pattern65, lambda p, c, d, a, f, x, m, e, b : -b*c*d**p*Int((f*x)**(m + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*ArcSin(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*ArcSin(c*x))*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule65)

    pattern66 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2) + S(1)/2)))
    rule66 = ReplacementRule(pattern66, lambda p, c, d, a, f, x, m, e, b : b*c*d**p*Int((f*x)**(m + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*ArcCos(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*ArcCos(c*x))*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule66)

    pattern67 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule67 = ReplacementRule(pattern67, lambda p, c, d, a, f, x, m, e, b : With(List(Set(u, IntHide((f*x)**m*(d + e*x**S(2))**p, x))), -b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcSin(c*x), u, x)))
    rubi.add(rule67)

    pattern68 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule68 = ReplacementRule(pattern68, lambda p, c, d, a, f, x, m, e, b : With(List(Set(u, IntHide((f*x)**m*(d + e*x**S(2))**p, x))), b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcCos(c*x), u, x)))
    rubi.add(rule68)

    pattern69 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda d: PositiveQ(d)))
    rule69 = ReplacementRule(pattern69, lambda p, c, d, a, x, m, e, b : With(List(Set(u, IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x))), -b*c*d**p*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(d**p*(a + b*ArcSin(c*x)), u, x)))
    rubi.add(rule69)

    pattern70 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda d: PositiveQ(d)))
    rule70 = ReplacementRule(pattern70, lambda p, c, d, a, x, m, e, b : With(List(Set(u, IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x))), b*c*d**p*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(d**p*(a + b*ArcCos(c*x)), u, x)))
    rubi.add(rule70)

    pattern71 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)))
    rule71 = ReplacementRule(pattern71, lambda p, c, d, a, x, m, e, b : With(List(Set(u, IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x))), -b*c*d**(p + S(-1)/2)*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x)*Sqrt(d + e*x**S(2))/Sqrt(-c**S(2)*x**S(2) + S(1)) + (a + b*ArcSin(c*x))*Int(x**m*(d + e*x**S(2))**p, x)))
    rubi.add(rule71)

    pattern72 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)))
    rule72 = ReplacementRule(pattern72, lambda p, c, d, a, x, m, e, b : With(List(Set(u, IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x))), b*c*d**(p + S(-1)/2)*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x)*Sqrt(d + e*x**S(2))/Sqrt(-c**S(2)*x**S(2) + S(1)) + (a + b*ArcCos(c*x))*Int(x**m*(d + e*x**S(2))**p, x)))
    rubi.add(rule72)

    pattern73 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule73 = ReplacementRule(pattern73, lambda n, c, d, a, f, x, m, e, b : -b*c*n*Int((f*x)**(m + S(1))*(a + b*ArcSin(c*x))**(n + S(-1)), x)*Sqrt(d + e*x**S(2))/(f*(m + S(1))*Sqrt(-c**S(2)*x**S(2) + S(1))) + c**S(2)*Int((f*x)**(m + S(2))*(a + b*ArcSin(c*x))**n/Sqrt(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(d + e*x**S(2))/(f**S(2)*(m + S(1))*Sqrt(-c**S(2)*x**S(2) + S(1))) + (f*x)**(m + S(1))*(a + b*ArcSin(c*x))**n*Sqrt(d + e*x**S(2))/(f*(m + S(1))))
    rubi.add(rule73)

    pattern74 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule74 = ReplacementRule(pattern74, lambda n, c, d, a, f, x, m, e, b : b*c*n*Int((f*x)**(m + S(1))*(a + b*ArcCos(c*x))**(n + S(-1)), x)*Sqrt(d + e*x**S(2))/(f*(m + S(1))*Sqrt(-c**S(2)*x**S(2) + S(1))) + c**S(2)*Int((f*x)**(m + S(2))*(a + b*ArcCos(c*x))**n/Sqrt(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(d + e*x**S(2))/(f**S(2)*(m + S(1))*Sqrt(-c**S(2)*x**S(2) + S(1))) + (f*x)**(m + S(1))*(a + b*ArcCos(c*x))**n*Sqrt(d + e*x**S(2))/(f*(m + S(1))))
    rubi.add(rule74)

    pattern75 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule75 = ReplacementRule(pattern75, lambda p, n, c, d, a, f, x, m, e, b : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*ArcSin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule75)

    pattern76 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule76 = ReplacementRule(pattern76, lambda p, n, c, d, a, f, x, m, e, b : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*ArcCos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule76)

    pattern77 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: ~(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule77 = ReplacementRule(pattern77, lambda n, c, d, a, f, x, m, e, b : -b*c*n*Int((f*x)**(m + S(1))*(a + b*ArcSin(c*x))**(n + S(-1)), x)*Sqrt(d + e*x**S(2))/(f*(m + S(2))*Sqrt(-c**S(2)*x**S(2) + S(1))) + Int((f*x)**m*(a + b*ArcSin(c*x))**n/Sqrt(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(d + e*x**S(2))/((m + S(2))*Sqrt(-c**S(2)*x**S(2) + S(1))) + (f*x)**(m + S(1))*(a + b*ArcSin(c*x))**n*Sqrt(d + e*x**S(2))/(f*(m + S(2))))
    rubi.add(rule77)

    pattern78 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: ~(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule78 = ReplacementRule(pattern78, lambda n, c, d, a, f, x, m, e, b : b*c*n*Int((f*x)**(m + S(1))*(a + b*ArcCos(c*x))**(n + S(-1)), x)*Sqrt(d + e*x**S(2))/(f*(m + S(2))*Sqrt(-c**S(2)*x**S(2) + S(1))) + Int((f*x)**m*(a + b*ArcCos(c*x))**n/Sqrt(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(d + e*x**S(2))/((m + S(2))*Sqrt(-c**S(2)*x**S(2) + S(1))) + (f*x)**(m + S(1))*(a + b*ArcCos(c*x))**n*Sqrt(d + e*x**S(2))/(f*(m + S(2))))
    rubi.add(rule78)

    pattern79 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: ~(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule79 = ReplacementRule(pattern79, lambda p, n, c, d, a, f, x, m, e, b : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*ArcSin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(2)*p + S(1))) + S(2)*d*p*Int((f*x)**m*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(m + S(2)*p + S(1)) + (f*x)**(m + S(1))*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(2)*p + S(1))))
    rubi.add(rule79)

    pattern80 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: ~(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule80 = ReplacementRule(pattern80, lambda p, n, c, d, a, f, x, m, e, b : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*ArcCos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(2)*p + S(1))) + S(2)*d*p*Int((f*x)**m*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(m + S(2)*p + S(1)) + (f*x)**(m + S(1))*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(2)*p + S(1))))
    rubi.add(rule80)

    pattern81 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule81 = ReplacementRule(pattern81, lambda p, n, c, d, a, f, x, m, e, b : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*ArcSin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + c**S(2)*(m + S(2)*p + S(3))*Int((f*x)**(m + S(2))*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**p, x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule81)

    pattern82 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule82 = ReplacementRule(pattern82, lambda p, n, c, d, a, f, x, m, e, b : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*ArcCos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + c**S(2)*(m + S(2)*p + S(3))*Int((f*x)**(m + S(2))*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**p, x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule82)

    pattern83 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule83 = ReplacementRule(pattern83, lambda p, n, c, d, a, f, x, m, e, b : b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*ArcSin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) - f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*e*(p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule83)

    pattern84 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule84 = ReplacementRule(pattern84, lambda p, n, c, d, a, f, x, m, e, b : -b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*ArcCos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) - f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*e*(p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule84)

    pattern85 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: ~(RationalQ(m) & Greater(m, S(1)))), CustomConstraint(lambda p, n, m: IntegerQ(m) | IntegerQ(p) | Equal(n, S(1))))
    rule85 = ReplacementRule(pattern85, lambda p, n, c, d, a, f, x, m, e, b : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*ArcSin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*f*(p + S(1))) + (m + S(2)*p + S(3))*Int((f*x)**m*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))) - (f*x)**(m + S(1))*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*f*(p + S(1))))
    rubi.add(rule85)

    pattern86 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: ~(RationalQ(m) & Greater(m, S(1)))), CustomConstraint(lambda p, n, m: IntegerQ(m) | IntegerQ(p) | Equal(n, S(1))))
    rule86 = ReplacementRule(pattern86, lambda p, n, c, d, a, f, x, m, e, b : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*ArcCos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*f*(p + S(1))) + (m + S(2)*p + S(3))*Int((f*x)**m*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))) - (f*x)**(m + S(1))*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*f*(p + S(1))))
    rubi.add(rule86)

    pattern87 = Pattern(Integral((x_*WC('f', S(1)))**m_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule87 = ReplacementRule(pattern87, lambda n, c, d, a, f, x, m, e, b : b*f*n*Int((f*x)**(m + S(-1))*(a + b*ArcSin(c*x))**(n + S(-1)), x)*Sqrt(-c**S(2)*x**S(2) + S(1))/(c*m*Sqrt(d + e*x**S(2))) + f*(f*x)**(m + S(-1))*(a + b*ArcSin(c*x))**n*Sqrt(d + e*x**S(2))/(e*m) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*ArcSin(c*x))**n/Sqrt(d + e*x**S(2)), x)/(c**S(2)*m))
    rubi.add(rule87)

    pattern88 = Pattern(Integral((x_*WC('f', S(1)))**m_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule88 = ReplacementRule(pattern88, lambda n, c, d, a, f, x, m, e, b : -b*f*n*Int((f*x)**(m + S(-1))*(a + b*ArcCos(c*x))**(n + S(-1)), x)*Sqrt(-c**S(2)*x**S(2) + S(1))/(c*m*Sqrt(d + e*x**S(2))) + f*(f*x)**(m + S(-1))*(a + b*ArcCos(c*x))**n*Sqrt(d + e*x**S(2))/(e*m) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*ArcCos(c*x))**n/Sqrt(d + e*x**S(2)), x)/(c**S(2)*m))
    rubi.add(rule88)

    pattern89 = Pattern(Integral(x_**m_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule89 = ReplacementRule(pattern89, lambda n, c, d, a, x, m, e, b : c**(-m + S(-1))*Subst(Int((a + b*x)**n*Sin(x)**m, x), x, ArcSin(c*x))/Sqrt(d))
    rubi.add(rule89)

    pattern90 = Pattern(Integral(x_**m_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule90 = ReplacementRule(pattern90, lambda n, c, d, a, x, m, e, b : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*Cos(x)**m, x), x, ArcCos(c*x))/Sqrt(d))
    rubi.add(rule90)

    pattern91 = Pattern(Integral((x_*WC('f', S(1)))**m_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule91 = ReplacementRule(pattern91, lambda c, d, a, f, x, m, e, b : -b*c*(f*x)**(m + S(2))*HypergeometricPFQ(List(S(1), m/S(2) + S(1), m/S(2) + S(1)), List(m/S(2) + S(3)/2, m/S(2) + S(2)), c**S(2)*x**S(2))/(f**S(2)*(m + S(1))*(m + S(2))*Sqrt(d)) + (f*x)**(m + S(1))*(a + b*ArcSin(c*x))*Hypergeometric2F1(S(1)/2, m/S(2) + S(1)/2, m/S(2) + S(3)/2, c**S(2)*x**S(2))/(f*(m + S(1))*Sqrt(d)))
    rubi.add(rule91)

    pattern92 = Pattern(Integral((x_*WC('f', S(1)))**m_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule92 = ReplacementRule(pattern92, lambda c, d, a, f, x, m, e, b : b*c*(f*x)**(m + S(2))*HypergeometricPFQ(List(S(1), m/S(2) + S(1), m/S(2) + S(1)), List(m/S(2) + S(3)/2, m/S(2) + S(2)), c**S(2)*x**S(2))/(f**S(2)*(m + S(1))*(m + S(2))*Sqrt(d)) + (f*x)**(m + S(1))*(a + b*ArcCos(c*x))*Hypergeometric2F1(S(1)/2, m/S(2) + S(1)/2, m/S(2) + S(3)/2, c**S(2)*x**S(2))/(f*(m + S(1))*Sqrt(d)))
    rubi.add(rule92)

    pattern93 = Pattern(Integral((x_*WC('f', S(1)))**m_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda d: ~(PositiveQ(d))), CustomConstraint(lambda n, m: IntegerQ(m) | Equal(n, S(1))))
    rule93 = ReplacementRule(pattern93, lambda n, c, d, a, f, x, m, e, b : Int((f*x)**m*(a + b*ArcSin(c*x))**n/Sqrt(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(-c**S(2)*x**S(2) + S(1))/Sqrt(d + e*x**S(2)))
    rubi.add(rule93)

    pattern94 = Pattern(Integral((x_*WC('f', S(1)))**m_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda d: ~(PositiveQ(d))), CustomConstraint(lambda n, m: IntegerQ(m) | Equal(n, S(1))))
    rule94 = ReplacementRule(pattern94, lambda n, c, d, a, f, x, m, e, b : Int((f*x)**m*(a + b*ArcCos(c*x))**n/Sqrt(-c**S(2)*x**S(2) + S(1)), x)*Sqrt(-c**S(2)*x**S(2) + S(1))/Sqrt(d + e*x**S(2)))
    rubi.add(rule94)

    pattern95 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule95 = ReplacementRule(pattern95, lambda p, n, c, d, a, f, x, m, e, b : b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*ArcSin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(c*(m + S(2)*p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(e*(m + S(2)*p + S(1))) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**p, x)/(c**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule95)

    pattern96 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule96 = ReplacementRule(pattern96, lambda p, n, c, d, a, f, x, m, e, b : -b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*ArcCos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(c*(m + S(2)*p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(e*(m + S(2)*p + S(1))) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**p, x)/(c**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule96)

    pattern97 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(1))))
    rule97 = ReplacementRule(pattern97, lambda p, n, c, d, a, f, m, x, e, b : -d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*ArcSin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) + (f*x)**m*(a + b*ArcSin(c*x))**(n + S(1))*(d + e*x**S(2))**p*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule97)

    pattern98 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(1))))
    rule98 = ReplacementRule(pattern98, lambda p, n, c, d, a, f, m, x, e, b : d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*ArcCos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) - (f*x)**m*(a + b*ArcCos(c*x))**(n + S(1))*(d + e*x**S(2))**p*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule98)

    pattern99 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda d: PositiveQ(d)))
    rule99 = ReplacementRule(pattern99, lambda n, c, d, a, f, m, x, e, b : -f*m*Int((f*x)**(m + S(-1))*(a + b*ArcSin(c*x))**(n + S(1)), x)/(b*c*(n + S(1))*Sqrt(d)) + (f*x)**m*(a + b*ArcSin(c*x))**(n + S(1))/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule99)

    pattern100 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda d: PositiveQ(d)))
    rule100 = ReplacementRule(pattern100, lambda n, c, d, a, f, m, x, e, b : f*m*Int((f*x)**(m + S(-1))*(a + b*ArcCos(c*x))**(n + S(1)), x)/(b*c*(n + S(1))*Sqrt(d)) - (f*x)**m*(a + b*ArcCos(c*x))**(n + S(1))/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule100)

    pattern101 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(-3))), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)))
    rule101 = ReplacementRule(pattern101, lambda p, n, c, d, a, f, m, x, e, b : c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*(m + S(2)*p + S(1))*Int((f*x)**(m + S(1))*(a + b*ArcSin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*f*(n + S(1))) - d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*ArcSin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) + (f*x)**m*(a + b*ArcSin(c*x))**(n + S(1))*(d + e*x**S(2))**p*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule101)

    pattern102 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(-3))), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)))
    rule102 = ReplacementRule(pattern102, lambda p, n, c, d, a, f, m, x, e, b : -c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*(m + S(2)*p + S(1))*Int((f*x)**(m + S(1))*(a + b*ArcCos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*f*(n + S(1))) + d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*ArcCos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) - (f*x)**m*(a + b*ArcCos(c*x))**(n + S(1))*(d + e*x**S(2))**p*Sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule102)

    pattern103 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, p: IntegerQ(p) | PositiveQ(d)))
    rule103 = ReplacementRule(pattern103, lambda p, n, c, d, a, m, x, e, b : c**(-m + S(-1))*d**p*Subst(Int((a + b*x)**n*Cos(x)**(S(2)*p + S(1))*Sin(x)**m, x), x, ArcSin(c*x)))
    rubi.add(rule103)

    pattern104 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, p: IntegerQ(p) | PositiveQ(d)))
    rule104 = ReplacementRule(pattern104, lambda p, n, c, d, a, m, x, e, b : -c**(-m + S(-1))*d**p*Subst(Int((a + b*x)**n*Cos(x)**m*Sin(x)**(S(2)*p + S(1)), x), x, ArcCos(c*x)))
    rubi.add(rule104)

    pattern105 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, p: ~(IntegerQ(p) | PositiveQ(d))))
    rule105 = ReplacementRule(pattern105, lambda p, n, c, d, a, m, x, e, b : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x**m*(a + b*ArcSin(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x))
    rubi.add(rule105)

    pattern106 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, p: ~(IntegerQ(p) | PositiveQ(d))))
    rule106 = ReplacementRule(pattern106, lambda p, n, c, d, a, m, x, e, b : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x**m*(a + b*ArcCos(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x))
    rubi.add(rule106)

    pattern107 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda m: ~(PositiveIntegerQ(m/S(2) + S(1)/2))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Less(S(-3), m, S(0))))
    rule107 = ReplacementRule(pattern107, lambda p, n, c, d, a, f, x, m, e, b : Int(ExpandIntegrand((a + b*ArcSin(c*x))**n/Sqrt(d + e*x**S(2)), (f*x)**m*(d + e*x**S(2))**(p + S(1)/2), x), x))
    rubi.add(rule107)

    pattern108 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda m: ~(PositiveIntegerQ(m/S(2) + S(1)/2))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Less(S(-3), m, S(0))))
    rule108 = ReplacementRule(pattern108, lambda p, n, c, d, a, f, x, m, e, b : Int(ExpandIntegrand((a + b*ArcCos(c*x))**n/Sqrt(d + e*x**S(2)), (f*x)**m*(d + e*x**S(2))**(p + S(1)/2), x), x))
    rubi.add(rule108)

    pattern109 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule109 = ReplacementRule(pattern109, lambda p, c, d, a, x, e, b : -b*c*Int((d + e*x**S(2))**(p + S(1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(S(2)*e*(p + S(1))) + (a + b*ArcSin(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule109)

    pattern110 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule110 = ReplacementRule(pattern110, lambda p, c, d, a, x, e, b : b*c*Int((d + e*x**S(2))**(p + S(1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(S(2)*e*(p + S(1))) + (a + b*ArcCos(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule110)

    pattern111 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (LessEqual(m + p, S(0)) & PositiveIntegerQ(m/S(2) + S(-1)/2))))
    rule111 = ReplacementRule(pattern111, lambda p, c, d, a, f, m, x, e, b : With(List(Set(u, IntHide((f*x)**m*(d + e*x**S(2))**p, x))), -b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcSin(c*x), u, x)))
    rubi.add(rule111)

    pattern112 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (LessEqual(m + p, S(0)) & PositiveIntegerQ(m/S(2) + S(-1)/2))))
    rule112 = ReplacementRule(pattern112, lambda p, c, d, a, f, m, x, e, b : With(List(Set(u, IntHide((f*x)**m*(d + e*x**S(2))**p, x))), b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcCos(c*x), u, x)))
    rubi.add(rule112)

    pattern113 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m)))
    rule113 = ReplacementRule(pattern113, lambda p, n, c, d, a, f, m, x, e, b : Int(ExpandIntegrand((a + b*ArcSin(c*x))**n, (f*x)**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule113)

    pattern114 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, c, e: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m)))
    rule114 = ReplacementRule(pattern114, lambda p, n, c, d, a, f, m, x, e, b : Int(ExpandIntegrand((a + b*ArcCos(c*x))**n, (f*x)**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule114)

    pattern115 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule115 = ReplacementRule(pattern115, lambda p, n, c, d, a, f, m, x, e, b : Int((f*x)**m*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule115)

    pattern116 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule116 = ReplacementRule(pattern116, lambda p, n, c, d, a, f, m, x, e, b : Int((f*x)**m*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule116)

    pattern117 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, f, g, e: ZeroQ(d*g + e*f)), CustomConstraint(lambda f, g, c: ZeroQ(c**S(2)*f**S(2) - g**S(2))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule117 = ReplacementRule(pattern117, lambda p, n, c, d, a, f, m, x, e, h, b, g : (d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p))*Int((h*x)**m*(a + b*ArcSin(c*x))**n*(d*f + e*g*x**S(2))**p, x))
    rubi.add(rule117)

    pattern118 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, f, g, e: ZeroQ(d*g + e*f)), CustomConstraint(lambda f, g, c: ZeroQ(c**S(2)*f**S(2) - g**S(2))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule118 = ReplacementRule(pattern118, lambda p, n, c, d, a, f, m, x, e, h, b, g : (d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p))*Int((h*x)**m*(a + b*ArcCos(c*x))**n*(d*f + e*g*x**S(2))**p, x))
    rubi.add(rule118)

    pattern119 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule119 = ReplacementRule(pattern119, lambda n, c, d, a, x, e, b : Subst(Int((a + b*x)**n*Cos(x)/(c*d + e*Sin(x)), x), x, ArcSin(c*x)))
    rubi.add(rule119)

    pattern120 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule120 = ReplacementRule(pattern120, lambda n, c, d, a, x, e, b : -Subst(Int((a + b*x)**n*Sin(x)/(c*d + e*Cos(x)), x), x, ArcCos(c*x)))
    rubi.add(rule120)

    pattern121 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule121 = ReplacementRule(pattern121, lambda n, c, d, a, m, x, e, b : -b*c*n*Int((a + b*ArcSin(c*x))**(n + S(-1))*(d + e*x)**(m + S(1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(e*(m + S(1))) + (a + b*ArcSin(c*x))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))))
    rubi.add(rule121)

    pattern122 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule122 = ReplacementRule(pattern122, lambda n, c, d, a, m, x, e, b : b*c*n*Int((a + b*ArcCos(c*x))**(n + S(-1))*(d + e*x)**(m + S(1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x)/(e*(m + S(1))) + (a + b*ArcCos(c*x))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))))
    rubi.add(rule122)

    pattern123 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule123 = ReplacementRule(pattern123, lambda n, c, d, a, m, x, e, b : Int(ExpandIntegrand((a + b*ArcSin(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule123)

    pattern124 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule124 = ReplacementRule(pattern124, lambda n, c, d, a, m, x, e, b : Int(ExpandIntegrand((a + b*ArcCos(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule124)

    pattern125 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule125 = ReplacementRule(pattern125, lambda n, c, d, a, m, x, e, b : c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*d + e*Sin(x))**m*Cos(x), x), x, ArcSin(c*x)))
    rubi.add(rule125)

    pattern126 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule126 = ReplacementRule(pattern126, lambda n, c, d, a, m, x, e, b : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*d + e*Cos(x))**m*Sin(x), x), x, ArcCos(c*x)))
    rubi.add(rule126)

    pattern127 = Pattern(Integral(Px_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)))
    rule127 = ReplacementRule(pattern127, lambda c, a, Px, x, b : With(List(Set(u, IntHide(Px, x))), -b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcSin(c*x), u, x)))
    rubi.add(rule127)

    pattern128 = Pattern(Integral(Px_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)))
    rule128 = ReplacementRule(pattern128, lambda c, a, Px, x, b : With(List(Set(u, IntHide(Px, x))), b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcCos(c*x), u, x)))
    rubi.add(rule128)

    pattern129 = Pattern(Integral(Px_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)))
    rule129 = ReplacementRule(pattern129, lambda n, c, a, Px, x, b : Int(ExpandIntegrand(Px*(a + b*ArcSin(c*x))**n, x), x))
    rubi.add(rule129)

    pattern130 = Pattern(Integral(Px_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)))
    rule130 = ReplacementRule(pattern130, lambda n, c, a, Px, x, b : Int(ExpandIntegrand(Px*(a + b*ArcCos(c*x))**n, x), x))
    rubi.add(rule130)

    pattern131 = Pattern(Integral(Px_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)))
    rule131 = ReplacementRule(pattern131, lambda c, d, a, Px, m, x, e, b : With(List(Set(u, IntHide(Px*(d + e*x)**m, x))), -b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcSin(c*x), u, x)))
    rubi.add(rule131)

    pattern132 = Pattern(Integral(Px_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)))
    rule132 = ReplacementRule(pattern132, lambda c, d, a, Px, m, x, e, b : With(List(Set(u, IntHide(Px*(d + e*x)**m, x))), b*c*Int(SimplifyIntegrand(u/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcCos(c*x), u, x)))
    rubi.add(rule132)

    pattern133 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda p, m: Less(m + p + S(1), S(0))))
    rule133 = ReplacementRule(pattern133, lambda p, n, c, d, a, f, x, m, e, b, g : With(List(Set(u, IntHide((d + e*x)**m*(f + g*x)**p, x))), -b*c*n*Int(SimplifyIntegrand(u*(a + b*ArcSin(c*x))**(n + S(-1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist((a + b*ArcSin(c*x))**n, u, x)))
    rubi.add(rule133)

    pattern134 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda p, m: Less(m + p + S(1), S(0))))
    rule134 = ReplacementRule(pattern134, lambda p, n, c, d, a, f, x, m, e, b, g : With(List(Set(u, IntHide((d + e*x)**m*(f + g*x)**p, x))), b*c*n*Int(SimplifyIntegrand(u*(a + b*ArcCos(c*x))**(n + S(-1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist((a + b*ArcCos(c*x))**n, u, x)))
    rubi.add(rule134)

    pattern135 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_*(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))/(d_ + x_*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)), CustomConstraint(lambda d, g, h, e: ZeroQ(-S(2)*d*h + e*g)))
    rule135 = ReplacementRule(pattern135, lambda p, n, c, d, a, f, x, e, h, b, g : With(List(Set(u, IntHide((f + g*x + h*x**S(2))**p/(d + e*x)**S(2), x))), -b*c*n*Int(SimplifyIntegrand(u*(a + b*ArcSin(c*x))**(n + S(-1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist((a + b*ArcSin(c*x))**n, u, x)))
    rubi.add(rule135)

    pattern136 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_*(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))/(d_ + x_*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)), CustomConstraint(lambda d, g, h, e: ZeroQ(-S(2)*d*h + e*g)))
    rule136 = ReplacementRule(pattern136, lambda p, n, c, d, a, f, x, e, h, b, g : With(List(Set(u, IntHide((f + g*x + h*x**S(2))**p/(d + e*x)**S(2), x))), b*c*n*Int(SimplifyIntegrand(u*(a + b*ArcCos(c*x))**(n + S(-1))/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist((a + b*ArcCos(c*x))**n, u, x)))
    rubi.add(rule136)

    pattern137 = Pattern(Integral(Px_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule137 = ReplacementRule(pattern137, lambda n, c, d, a, Px, m, x, e, b : Int(ExpandIntegrand(Px*(a + b*ArcSin(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule137)

    pattern138 = Pattern(Integral(Px_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule138 = ReplacementRule(pattern138, lambda n, c, d, a, Px, m, x, e, b : Int(ExpandIntegrand(Px*(a + b*ArcCos(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule138)

    pattern139 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p, m: Greater(m, S(3)) | Less(m, -S(2)*p + S(-1))))
    rule139 = ReplacementRule(pattern139, lambda p, c, d, a, f, m, x, e, b, g : With(List(Set(u, IntHide((d + e*x**S(2))**p*(f + g*x)**m, x))), -b*c*Int(Dist(1/Sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x) + Dist(a + b*ArcSin(c*x), u, x)))
    rubi.add(rule139)

    pattern140 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p, m: Greater(m, S(3)) | Less(m, -S(2)*p + S(-1))))
    rule140 = ReplacementRule(pattern140, lambda p, c, d, a, f, m, x, e, b, g : With(List(Set(u, IntHide((d + e*x**S(2))**p*(f + g*x)**m, x))), b*c*Int(Dist(1/Sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x) + Dist(a + b*ArcCos(c*x), u, x)))
    rubi.add(rule140)

    pattern141 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p, n, m: Equal(m, S(1)) | Greater(p, S(0)) | (Equal(m, S(2)) & Less(p, S(-2))) | (Equal(n, S(1)) & Greater(p, S(-1)))))
    rule141 = ReplacementRule(pattern141, lambda p, n, c, d, a, f, m, x, e, b, g : Int(ExpandIntegrand((a + b*ArcSin(c*x))**n*(d + e*x**S(2))**p, (f + g*x)**m, x), x))
    rubi.add(rule141)

    pattern142 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p, n, m: Equal(m, S(1)) | Greater(p, S(0)) | (Equal(m, S(2)) & Less(p, S(-2))) | (Equal(n, S(1)) & Greater(p, S(-1)))))
    rule142 = ReplacementRule(pattern142, lambda p, n, c, d, a, f, m, x, e, b, g : Int(ExpandIntegrand((a + b*ArcCos(c*x))**n*(d + e*x**S(2))**p, (f + g*x)**m, x), x))
    rubi.add(rule142)

    pattern143 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(f_ + x_*WC('g', S(1)))**m_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule143 = ReplacementRule(pattern143, lambda n, c, d, a, f, x, m, e, b, g : (a + b*ArcSin(c*x))**(n + S(1))*(d + e*x**S(2))*(f + g*x)**m/(b*c*(n + S(1))*Sqrt(d)) - Int((a + b*ArcSin(c*x))**(n + S(1))*(f + g*x)**(m + S(-1))*(d*g*m + S(2)*e*f*x + e*g*x**S(2)*(m + S(2))), x)/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule143)

    pattern144 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(f_ + x_*WC('g', S(1)))**m_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule144 = ReplacementRule(pattern144, lambda n, c, d, a, f, x, m, e, b, g : -(a + b*ArcCos(c*x))**(n + S(1))*(d + e*x**S(2))*(f + g*x)**m/(b*c*(n + S(1))*Sqrt(d)) + Int((a + b*ArcCos(c*x))**(n + S(1))*(f + g*x)**(m + S(-1))*(d*g*m + S(2)*e*f*x + e*g*x**S(2)*(m + S(2))), x)/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule144)

    pattern145 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule145 = ReplacementRule(pattern145, lambda p, n, c, d, a, f, m, x, e, b, g : Int(ExpandIntegrand((a + b*ArcSin(c*x))**n*Sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(-1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule145)

    pattern146 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule146 = ReplacementRule(pattern146, lambda p, n, c, d, a, f, m, x, e, b, g : Int(ExpandIntegrand((a + b*ArcCos(c*x))**n*Sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(-1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule146)

    pattern147 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule147 = ReplacementRule(pattern147, lambda p, n, c, d, a, f, m, x, e, b, g : (a + b*ArcSin(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m/(b*c*(n + S(1))*Sqrt(d)) - Int(ExpandIntegrand((a + b*ArcSin(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), (d + e*x**S(2))**(p + S(-1)/2)*(d*g*m + e*f*x*(S(2)*p + S(1)) + e*g*x**S(2)*(m + S(2)*p + S(1))), x), x)/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule147)

    pattern148 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule148 = ReplacementRule(pattern148, lambda p, n, c, d, a, f, m, x, e, b, g : -(a + b*ArcCos(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m/(b*c*(n + S(1))*Sqrt(d)) + Int(ExpandIntegrand((a + b*ArcCos(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), (d + e*x**S(2))**(p + S(-1)/2)*(d*g*m + e*f*x*(S(2)*p + S(1)) + e*g*x**S(2)*(m + S(2)*p + S(1))), x), x)/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule148)

    pattern149 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule149 = ReplacementRule(pattern149, lambda n, c, d, a, f, m, x, e, b, g : -g*m*Int((a + b*ArcSin(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), x)/(b*c*(n + S(1))*Sqrt(d)) + (a + b*ArcSin(c*x))**(n + S(1))*(f + g*x)**m/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule149)

    pattern150 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule150 = ReplacementRule(pattern150, lambda n, c, d, a, f, m, x, e, b, g : g*m*Int((a + b*ArcCos(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), x)/(b*c*(n + S(1))*Sqrt(d)) - (a + b*ArcCos(c*x))**(n + S(1))*(f + g*x)**m/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule150)

    pattern151 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n, m: PositiveIntegerQ(n) | Greater(m, S(0))))
    rule151 = ReplacementRule(pattern151, lambda n, c, d, a, f, m, x, e, b, g : c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*f + g*Sin(x))**m, x), x, ArcSin(c*x))/Sqrt(d))
    rubi.add(rule151)

    pattern152 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n, m: PositiveIntegerQ(n) | Greater(m, S(0))))
    rule152 = ReplacementRule(pattern152, lambda n, c, d, a, f, m, x, e, b, g : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*f + g*Cos(x))**m, x), x, ArcCos(c*x))/Sqrt(d))
    rubi.add(rule152)

    pattern153 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule153 = ReplacementRule(pattern153, lambda p, n, c, d, a, f, m, x, e, b, g : Int(ExpandIntegrand((a + b*ArcSin(c*x))**n/Sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule153)

    pattern154 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule154 = ReplacementRule(pattern154, lambda p, n, c, d, a, f, m, x, e, b, g : Int(ExpandIntegrand((a + b*ArcCos(c*x))**n/Sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule154)

    pattern155 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: ~(PositiveQ(d))))
    rule155 = ReplacementRule(pattern155, lambda p, n, c, d, a, f, m, x, e, b, g : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*ArcSin(c*x))**n*(f + g*x)**m*(-c**S(2)*x**S(2) + S(1))**p, x))
    rubi.add(rule155)

    pattern156 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: ~(PositiveQ(d))))
    rule156 = ReplacementRule(pattern156, lambda p, n, c, d, a, f, m, x, e, b, g : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*ArcCos(c*x))**n*(f + g*x)**m*(-c**S(2)*x**S(2) + S(1))**p, x))
    rubi.add(rule156)

    pattern157 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule157 = ReplacementRule(pattern157, lambda n, c, d, a, f, m, x, e, h, b, g : -g*m*Int((a + b*ArcSin(c*x))**(n + S(1))/(f + g*x), x)/(b*c*(n + S(1))*Sqrt(d)) + (a + b*ArcSin(c*x))**(n + S(1))*Log(h*(f + g*x)**m)/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule157)

    pattern158 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule158 = ReplacementRule(pattern158, lambda n, c, d, a, f, m, x, e, h, b, g : g*m*Int((a + b*ArcCos(c*x))**(n + S(1))/(f + g*x), x)/(b*c*(n + S(1))*Sqrt(d)) - (a + b*ArcCos(c*x))**(n + S(1))*Log(h*(f + g*x)**m)/(b*c*(n + S(1))*Sqrt(d)))
    rubi.add(rule158)

    pattern159 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: ~(PositiveQ(d))))
    rule159 = ReplacementRule(pattern159, lambda p, n, c, d, a, f, m, x, e, h, b, g : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*ArcSin(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p*Log(h*(f + g*x)**m), x))
    rubi.add(rule159)

    pattern160 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: ~(PositiveQ(d))))
    rule160 = ReplacementRule(pattern160, lambda p, n, c, d, a, f, m, x, e, h, b, g : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*ArcCos(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p*Log(h*(f + g*x)**m), x))
    rubi.add(rule160)

    pattern161 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**m_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m + S(1)/2)))
    rule161 = ReplacementRule(pattern161, lambda c, d, a, f, x, m, e, b, g : With(List(Set(u, IntHide((d + e*x)**m*(f + g*x)**m, x))), -b*c*Int(Dist(1/Sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x) + Dist(a + b*ArcSin(c*x), u, x)))
    rubi.add(rule161)

    pattern162 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**m_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m + S(1)/2)))
    rule162 = ReplacementRule(pattern162, lambda c, d, a, f, x, m, e, b, g : With(List(Set(u, IntHide((d + e*x)**m*(f + g*x)**m, x))), b*c*Int(Dist(1/Sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x) + Dist(a + b*ArcCos(c*x), u, x)))
    rubi.add(rule162)

    pattern163 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule163 = ReplacementRule(pattern163, lambda n, c, d, a, f, m, x, e, b, g : Int(ExpandIntegrand((a + b*ArcSin(c*x))**n*(d + e*x)**m*(f + g*x)**m, x), x))
    rubi.add(rule163)

    pattern164 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule164 = ReplacementRule(pattern164, lambda n, c, d, a, f, m, x, e, b, g : Int(ExpandIntegrand((a + b*ArcCos(c*x))**n*(d + e*x)**m*(f + g*x)**m, x), x))
    rubi.add(rule164)

    pattern165 = Pattern(Integral(u_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule165 = ReplacementRule(pattern165, lambda c, u, a, x, b : With(List(Set(v, IntHide(u, x))), Condition(-b*c*Int(SimplifyIntegrand(v/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcSin(c*x), v, x), InverseFunctionFreeQ(v, x))))
    rubi.add(rule165)

    pattern166 = Pattern(Integral(u_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule166 = ReplacementRule(pattern166, lambda c, u, a, x, b : With(List(Set(v, IntHide(u, x))), Condition(b*c*Int(SimplifyIntegrand(v/Sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcCos(c*x), v, x), InverseFunctionFreeQ(v, x))))
    rubi.add(rule166)

    pattern167 = Pattern(Integral(Px_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule167 = ReplacementRule(pattern167, lambda p, n, c, d, a, Px, x, e, b : With(List(Set(u, ExpandIntegrand(Px*(a + b*ArcSin(c*x))**n*(d + e*x**S(2))**p, x))), Condition(Int(u, x), SumQ(u))))
    rubi.add(rule167)

    pattern168 = Pattern(Integral(Px_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule168 = ReplacementRule(pattern168, lambda p, n, c, d, a, Px, x, e, b : With(List(Set(u, ExpandIntegrand(Px*(a + b*ArcCos(c*x))**n*(d + e*x**S(2))**p, x))), Condition(Int(u, x), SumQ(u))))
    rubi.add(rule168)

    pattern169 = Pattern(Integral((f_ + (d_ + x_**S(2)*WC('e', S(1)))**p_*WC('g', S(1)))**WC('m', S(1))*(ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*WC('Px', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda n, m: IntegersQ(m, n)))
    rule169 = ReplacementRule(pattern169, lambda p, n, c, d, Px, a, f, m, x, e, b, g : With(List(Set(u, ExpandIntegrand(Px*(a + b*ArcSin(c*x))**n*(f + g*(d + e*x**S(2))**p)**m, x))), Condition(Int(u, x), SumQ(u))))
    rubi.add(rule169)

    pattern170 = Pattern(Integral((f_ + (d_ + x_**S(2)*WC('e', S(1)))**p_*WC('g', S(1)))**WC('m', S(1))*(ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*WC('Px', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda n, m: IntegersQ(m, n)))
    rule170 = ReplacementRule(pattern170, lambda p, n, c, d, Px, a, f, m, x, e, b, g : With(List(Set(u, ExpandIntegrand(Px*(a + b*ArcCos(c*x))**n*(f + g*(d + e*x**S(2))**p)**m, x))), Condition(Int(u, x), SumQ(u))))
    rubi.add(rule170)

    pattern171 = Pattern(Integral(RFx_*ArcSin(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule171 = ReplacementRule(pattern171, lambda RFx, n, c, x : With(List(Set(u, ExpandIntegrand(ArcSin(c*x)**n, RFx, x))), Condition(Int(u, x), SumQ(u))))
    rubi.add(rule171)

    pattern172 = Pattern(Integral(RFx_*ArcCos(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule172 = ReplacementRule(pattern172, lambda RFx, n, c, x : With(List(Set(u, ExpandIntegrand(ArcCos(c*x)**n, RFx, x))), Condition(Int(u, x), SumQ(u))))
    rubi.add(rule172)

    pattern173 = Pattern(Integral(RFx_*(a_ + ArcSin(x_*WC('c', S(1)))*WC('b', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule173 = ReplacementRule(pattern173, lambda RFx, n, c, a, x, b : Int(ExpandIntegrand(RFx*(a + b*ArcSin(c*x))**n, x), x))
    rubi.add(rule173)

    pattern174 = Pattern(Integral(RFx_*(a_ + ArcCos(x_*WC('c', S(1)))*WC('b', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule174 = ReplacementRule(pattern174, lambda RFx, n, c, a, x, b : Int(ExpandIntegrand(RFx*(a + b*ArcCos(c*x))**n, x), x))
    rubi.add(rule174)

    pattern175 = Pattern(Integral(RFx_*(d_ + x_**S(2)*WC('e', S(1)))**p_*ArcSin(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule175 = ReplacementRule(pattern175, lambda p, n, c, d, x, e, RFx : With(List(Set(u, ExpandIntegrand((d + e*x**S(2))**p*ArcSin(c*x)**n, RFx, x))), Condition(Int(u, x), SumQ(u))))
    rubi.add(rule175)

    pattern176 = Pattern(Integral(RFx_*(d_ + x_**S(2)*WC('e', S(1)))**p_*ArcCos(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule176 = ReplacementRule(pattern176, lambda p, n, c, d, x, e, RFx : With(List(Set(u, ExpandIntegrand((d + e*x**S(2))**p*ArcCos(c*x)**n, RFx, x))), Condition(Int(u, x), SumQ(u))))
    rubi.add(rule176)

    pattern177 = Pattern(Integral(RFx_*(a_ + ArcSin(x_*WC('c', S(1)))*WC('b', S(1)))**WC('n', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule177 = ReplacementRule(pattern177, lambda RFx, n, p, c, d, a, x, e, b : Int(ExpandIntegrand((d + e*x**S(2))**p, RFx*(a + b*ArcSin(c*x))**n, x), x))
    rubi.add(rule177)

    pattern178 = Pattern(Integral(RFx_*(a_ + ArcCos(x_*WC('c', S(1)))*WC('b', S(1)))**WC('n', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule178 = ReplacementRule(pattern178, lambda RFx, n, p, c, d, a, x, e, b : Int(ExpandIntegrand((d + e*x**S(2))**p, RFx*(a + b*ArcCos(c*x))**n, x), x))
    rubi.add(rule178)

    pattern179 = Pattern(Integral((ArcSin(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule179 = ReplacementRule(pattern179, lambda n, c, u, a, x, b : Int(u*(a + b*ArcSin(c*x))**n, x))
    rubi.add(rule179)

    pattern180 = Pattern(Integral((ArcCos(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule180 = ReplacementRule(pattern180, lambda n, c, u, a, x, b : Int(u*(a + b*ArcCos(c*x))**n, x))
    rubi.add(rule180)

    pattern181 = Pattern(Integral((ArcSin(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule181 = ReplacementRule(pattern181, lambda n, c, d, a, x, b : Subst(Int((a + b*ArcSin(x))**n, x), x, c + d*x)/d)
    rubi.add(rule181)

    pattern182 = Pattern(Integral((ArcCos(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule182 = ReplacementRule(pattern182, lambda n, c, d, a, x, b : Subst(Int((a + b*ArcCos(x))**n, x), x, c + d*x)/d)
    rubi.add(rule182)

    pattern183 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcSin(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule183 = ReplacementRule(pattern183, lambda n, c, d, a, f, m, x, e, b : Subst(Int((a + b*ArcSin(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule183)

    pattern184 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcCos(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule184 = ReplacementRule(pattern184, lambda n, c, d, a, f, m, x, e, b : Subst(Int((a + b*ArcCos(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule184)

    pattern185 = Pattern(Integral((ArcSin(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, B, c, A: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda d, C, B, c: ZeroQ(-B*d + S(2)*C*c)))
    rule185 = ReplacementRule(pattern185, lambda p, n, B, c, d, a, C, x, A, b : Subst(Int((a + b*ArcSin(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule185)

    pattern186 = Pattern(Integral((ArcCos(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, B, c, A: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda d, C, B, c: ZeroQ(-B*d + S(2)*C*c)))
    rule186 = ReplacementRule(pattern186, lambda p, n, B, c, d, a, C, x, A, b : Subst(Int((a + b*ArcCos(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule186)

    pattern187 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcSin(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, B, c, A: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda d, C, B, c: ZeroQ(-B*d + S(2)*C*c)))
    rule187 = ReplacementRule(pattern187, lambda p, n, B, c, d, a, f, C, m, x, e, A, b : Subst(Int((a + b*ArcSin(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule187)

    pattern188 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcCos(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, B, c, A: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda d, C, B, c: ZeroQ(-B*d + S(2)*C*c)))
    rule188 = ReplacementRule(pattern188, lambda p, n, B, c, d, a, f, C, m, x, e, A, b : Subst(Int((a + b*ArcCos(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule188)

    pattern189 = Pattern(Integral(sqrt(ArcSin(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))))
    rule189 = ReplacementRule(pattern189, lambda c, d, a, x, b : x*(-c*Sin(a/(S(2)*b)) + Cos(a/(S(2)*b)))*FresnelS(Sqrt(c/(Pi*b))*Sqrt(a + b*ArcSin(c + d*x**S(2))))*Sqrt(Pi)/((-c*Sin(ArcSin(c + d*x**S(2))/S(2)) + Cos(ArcSin(c + d*x**S(2))/S(2)))*Sqrt(c/b)) - x*(c*Sin(a/(S(2)*b)) + Cos(a/(S(2)*b)))*FresnelC(Sqrt(c/(Pi*b))*Sqrt(a + b*ArcSin(c + d*x**S(2))))*Sqrt(Pi)/((-c*Sin(ArcSin(c + d*x**S(2))/S(2)) + Cos(ArcSin(c + d*x**S(2))/S(2)))*Sqrt(c/b)) + x*Sqrt(a + b*ArcSin(c + d*x**S(2))))
    rubi.add(rule189)

    pattern190 = Pattern(Integral(sqrt(ArcCos(x_**S(2)*WC('d', S(1)) + S(1))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule190 = ReplacementRule(pattern190, lambda d, a, b, x : S(2)*Cos(a/(S(2)*b))*FresnelS(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(1))))*Sin(ArcCos(d*x**S(2) + S(1))/S(2))*Sqrt(Pi)/(d*x*Sqrt(1/b)) - S(2)*FresnelC(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(1))))*Sin(a/(S(2)*b))*Sin(ArcCos(d*x**S(2) + S(1))/S(2))*Sqrt(Pi)/(d*x*Sqrt(1/b)) - S(2)*Sin(ArcCos(d*x**S(2) + S(1))/S(2))**S(2)*Sqrt(a + b*ArcCos(d*x**S(2) + S(1)))/(d*x))
    rubi.add(rule190)

    pattern191 = Pattern(Integral(sqrt(ArcCos(x_**S(2)*WC('d', S(1)) + S(-1))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule191 = ReplacementRule(pattern191, lambda d, a, b, x : -S(2)*Cos(a/(S(2)*b))*Cos(ArcCos(d*x**S(2) + S(-1))/S(2))*FresnelC(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(-1))))*Sqrt(Pi)/(d*x*Sqrt(1/b)) + S(2)*Cos(ArcCos(d*x**S(2) + S(-1))/S(2))**S(2)*Sqrt(a + b*ArcCos(d*x**S(2) + S(-1)))/(d*x) - S(2)*Cos(ArcCos(d*x**S(2) + S(-1))/S(2))*FresnelS(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(-1))))*Sin(a/(S(2)*b))*Sqrt(Pi)/(d*x*Sqrt(1/b)))
    rubi.add(rule191)

    pattern192 = Pattern(Integral((ArcSin(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule192 = ReplacementRule(pattern192, lambda n, c, d, a, x, b : -S(4)*b**S(2)*n*(n + S(-1))*Int((a + b*ArcSin(c + d*x**S(2)))**(n + S(-2)), x) + S(2)*b*n*(a + b*ArcSin(c + d*x**S(2)))**(n + S(-1))*Sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(d*x) + x*(a + b*ArcSin(c + d*x**S(2)))**n)
    rubi.add(rule192)

    pattern193 = Pattern(Integral((ArcCos(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule193 = ReplacementRule(pattern193, lambda n, c, d, a, x, b : -S(4)*b**S(2)*n*(n + S(-1))*Int((a + b*ArcCos(c + d*x**S(2)))**(n + S(-2)), x) - S(2)*b*n*(a + b*ArcCos(c + d*x**S(2)))**(n + S(-1))*Sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(d*x) + x*(a + b*ArcCos(c + d*x**S(2)))**n)
    rubi.add(rule193)

    pattern194 = Pattern(Integral(1/(ArcSin(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))))
    rule194 = ReplacementRule(pattern194, lambda c, d, a, x, b : -x*(c*Cos(a/(S(2)*b)) - Sin(a/(S(2)*b)))*CosIntegral(c*(a + b*ArcSin(c + d*x**S(2)))/(S(2)*b))/(S(2)*b*(-c*Sin(ArcSin(c + d*x**S(2))/S(2)) + Cos(ArcSin(c + d*x**S(2))/S(2)))) - x*(c*Cos(a/(S(2)*b)) + Sin(a/(S(2)*b)))*SinIntegral(c*(a + b*ArcSin(c + d*x**S(2)))/(S(2)*b))/(S(2)*b*(-c*Sin(ArcSin(c + d*x**S(2))/S(2)) + Cos(ArcSin(c + d*x**S(2))/S(2)))))
    rubi.add(rule194)

    pattern195 = Pattern(Integral(1/(ArcCos(x_**S(2)*WC('d', S(1)) + S(1))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule195 = ReplacementRule(pattern195, lambda d, a, b, x : x*Cos(a/(S(2)*b))*CosIntegral((a + b*ArcCos(d*x**S(2) + S(1)))/(S(2)*b))/(b*Sqrt(S(2))*Sqrt(-d*x**S(2))) + x*Sin(a/(S(2)*b))*SinIntegral((a + b*ArcCos(d*x**S(2) + S(1)))/(S(2)*b))/(b*Sqrt(S(2))*Sqrt(-d*x**S(2))))
    rubi.add(rule195)

    pattern196 = Pattern(Integral(1/(ArcCos(x_**S(2)*WC('d', S(1)) + S(-1))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule196 = ReplacementRule(pattern196, lambda d, a, b, x : -x*Cos(a/(S(2)*b))*SinIntegral((a + b*ArcCos(d*x**S(2) + S(-1)))/(S(2)*b))/(b*Sqrt(S(2))*Sqrt(d*x**S(2))) + x*CosIntegral((a + b*ArcCos(d*x**S(2) + S(-1)))/(S(2)*b))*Sin(a/(S(2)*b))/(b*Sqrt(S(2))*Sqrt(d*x**S(2))))
    rubi.add(rule196)

    pattern197 = Pattern(Integral(1/sqrt(ArcSin(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))))
    rule197 = ReplacementRule(pattern197, lambda c, d, a, x, b : -x*(-c*Sin(a/(S(2)*b)) + Cos(a/(S(2)*b)))*FresnelC(Sqrt(a + b*ArcSin(c + d*x**S(2)))/(Sqrt(Pi)*Sqrt(b*c)))*Sqrt(Pi)/((-c*Sin(ArcSin(c + d*x**S(2))/S(2)) + Cos(ArcSin(c + d*x**S(2))/S(2)))*Sqrt(b*c)) - x*(c*Sin(a/(S(2)*b)) + Cos(a/(S(2)*b)))*FresnelS(Sqrt(a + b*ArcSin(c + d*x**S(2)))/(Sqrt(Pi)*Sqrt(b*c)))*Sqrt(Pi)/((-c*Sin(ArcSin(c + d*x**S(2))/S(2)) + Cos(ArcSin(c + d*x**S(2))/S(2)))*Sqrt(b*c)))
    rubi.add(rule197)

    pattern198 = Pattern(Integral(1/sqrt(ArcCos(x_**S(2)*WC('d', S(1)) + S(1))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule198 = ReplacementRule(pattern198, lambda d, a, b, x : -S(2)*Cos(a/(S(2)*b))*FresnelC(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(1))))*Sin(ArcCos(d*x**S(2) + S(1))/S(2))*Sqrt(Pi/b)/(d*x) - S(2)*FresnelS(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(1))))*Sin(a/(S(2)*b))*Sin(ArcCos(d*x**S(2) + S(1))/S(2))*Sqrt(Pi/b)/(d*x))
    rubi.add(rule198)

    pattern199 = Pattern(Integral(1/sqrt(ArcCos(x_**S(2)*WC('d', S(1)) + S(-1))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule199 = ReplacementRule(pattern199, lambda d, a, b, x : -S(2)*Cos(a/(S(2)*b))*Cos(ArcCos(d*x**S(2) + S(-1))/S(2))*FresnelS(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(-1))))*Sqrt(Pi/b)/(d*x) + S(2)*Cos(ArcCos(d*x**S(2) + S(-1))/S(2))*FresnelC(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(-1))))*Sin(a/(S(2)*b))*Sqrt(Pi/b)/(d*x))
    rubi.add(rule199)

    pattern200 = Pattern(Integral((ArcSin(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**(S(-3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))))
    rule200 = ReplacementRule(pattern200, lambda c, d, a, x, b : x*(c/b)**(S(3)/2)*(-c*Sin(a/(S(2)*b)) + Cos(a/(S(2)*b)))*FresnelS(Sqrt(c/(Pi*b))*Sqrt(a + b*ArcSin(c + d*x**S(2))))*Sqrt(Pi)/(-c*Sin(ArcSin(c + d*x**S(2))/S(2)) + Cos(ArcSin(c + d*x**S(2))/S(2))) - x*(c/b)**(S(3)/2)*(c*Sin(a/(S(2)*b)) + Cos(a/(S(2)*b)))*FresnelC(Sqrt(c/(Pi*b))*Sqrt(a + b*ArcSin(c + d*x**S(2))))*Sqrt(Pi)/(-c*Sin(ArcSin(c + d*x**S(2))/S(2)) + Cos(ArcSin(c + d*x**S(2))/S(2))) - Sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(b*d*x*Sqrt(a + b*ArcSin(c + d*x**S(2)))))
    rubi.add(rule200)

    pattern201 = Pattern(Integral((ArcCos(x_**S(2)*WC('d', S(1)) + S(1))*WC('b', S(1)) + WC('a', S(0)))**(S(-3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule201 = ReplacementRule(pattern201, lambda d, a, b, x : S(2)*(1/b)**(S(3)/2)*Cos(a/(S(2)*b))*FresnelS(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(1))))*Sin(ArcCos(d*x**S(2) + S(1))/S(2))*Sqrt(Pi)/(d*x) - S(2)*(1/b)**(S(3)/2)*FresnelC(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(1))))*Sin(a/(S(2)*b))*Sin(ArcCos(d*x**S(2) + S(1))/S(2))*Sqrt(Pi)/(d*x) + Sqrt(-d**S(2)*x**S(4) - S(2)*d*x**S(2))/(b*d*x*Sqrt(a + b*ArcCos(d*x**S(2) + S(1)))))
    rubi.add(rule201)

    pattern202 = Pattern(Integral((ArcCos(x_**S(2)*WC('d', S(1)) + S(-1))*WC('b', S(1)) + WC('a', S(0)))**(S(-3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule202 = ReplacementRule(pattern202, lambda d, a, b, x : -S(2)*(1/b)**(S(3)/2)*Cos(a/(S(2)*b))*Cos(ArcCos(d*x**S(2) + S(-1))/S(2))*FresnelC(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(-1))))*Sqrt(Pi)/(d*x) - S(2)*(1/b)**(S(3)/2)*Cos(ArcCos(d*x**S(2) + S(-1))/S(2))*FresnelS(Sqrt(S(1)/(Pi*b))*Sqrt(a + b*ArcCos(d*x**S(2) + S(-1))))*Sin(a/(S(2)*b))*Sqrt(Pi)/(d*x) + Sqrt(-d**S(2)*x**S(4) + S(2)*d*x**S(2))/(b*d*x*Sqrt(a + b*ArcCos(d*x**S(2) + S(-1)))))
    rubi.add(rule202)

    pattern203 = Pattern(Integral((ArcSin(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**(S(-2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))))
    rule203 = ReplacementRule(pattern203, lambda c, d, a, x, b : -Sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(S(2)*b*d*x*(a + b*ArcSin(c + d*x**S(2)))) + x*(-c*Sin(a/(S(2)*b)) + Cos(a/(S(2)*b)))*SinIntegral(c*(a + b*ArcSin(c + d*x**S(2)))/(S(2)*b))/(S(4)*b**S(2)*(-c*Sin(ArcSin(c + d*x**S(2))/S(2)) + Cos(ArcSin(c + d*x**S(2))/S(2)))) - x*(c*Sin(a/(S(2)*b)) + Cos(a/(S(2)*b)))*CosIntegral(c*(a + b*ArcSin(c + d*x**S(2)))/(S(2)*b))/(S(4)*b**S(2)*(-c*Sin(ArcSin(c + d*x**S(2))/S(2)) + Cos(ArcSin(c + d*x**S(2))/S(2)))))
    rubi.add(rule203)

    pattern204 = Pattern(Integral((ArcCos(x_**S(2)*WC('d', S(1)) + S(1))*WC('b', S(1)) + WC('a', S(0)))**(S(-2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule204 = ReplacementRule(pattern204, lambda d, a, b, x : Sqrt(-d**S(2)*x**S(4) - S(2)*d*x**S(2))/(S(2)*b*d*x*(a + b*ArcCos(d*x**S(2) + S(1)))) - x*Cos(a/(S(2)*b))*SinIntegral((a + b*ArcCos(d*x**S(2) + S(1)))/(S(2)*b))/(S(2)*b**S(2)*Sqrt(S(2))*Sqrt(-d*x**S(2))) + x*CosIntegral((a + b*ArcCos(d*x**S(2) + S(1)))/(S(2)*b))*Sin(a/(S(2)*b))/(S(2)*b**S(2)*Sqrt(S(2))*Sqrt(-d*x**S(2))))
    rubi.add(rule204)

    pattern205 = Pattern(Integral((ArcCos(x_**S(2)*WC('d', S(1)) + S(-1))*WC('b', S(1)) + WC('a', S(0)))**(S(-2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule205 = ReplacementRule(pattern205, lambda d, a, b, x : Sqrt(-d**S(2)*x**S(4) + S(2)*d*x**S(2))/(S(2)*b*d*x*(a + b*ArcCos(d*x**S(2) + S(-1)))) - x*Cos(a/(S(2)*b))*CosIntegral((a + b*ArcCos(d*x**S(2) + S(-1)))/(S(2)*b))/(S(2)*b**S(2)*Sqrt(S(2))*Sqrt(d*x**S(2))) - x*Sin(a/(S(2)*b))*SinIntegral((a + b*ArcCos(d*x**S(2) + S(-1)))/(S(2)*b))/(S(2)*b**S(2)*Sqrt(S(2))*Sqrt(d*x**S(2))))
    rubi.add(rule205)

    pattern206 = Pattern(Integral((ArcSin(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule206 = ReplacementRule(pattern206, lambda n, c, d, a, x, b : (a + b*ArcSin(c + d*x**S(2)))**(n + S(1))*Sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(S(2)*b*d*x*(n + S(1))) + x*(a + b*ArcSin(c + d*x**S(2)))**(n + S(2))/(S(4)*b**S(2)*(n + S(1))*(n + S(2))) - Int((a + b*ArcSin(c + d*x**S(2)))**(n + S(2)), x)/(S(4)*b**S(2)*(n + S(1))*(n + S(2))))
    rubi.add(rule206)

    pattern207 = Pattern(Integral((ArcCos(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule207 = ReplacementRule(pattern207, lambda n, c, d, a, x, b : -(a + b*ArcCos(c + d*x**S(2)))**(n + S(1))*Sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(S(2)*b*d*x*(n + S(1))) + x*(a + b*ArcCos(c + d*x**S(2)))**(n + S(2))/(S(4)*b**S(2)*(n + S(1))*(n + S(2))) - Int((a + b*ArcCos(c + d*x**S(2)))**(n + S(2)), x)/(S(4)*b**S(2)*(n + S(1))*(n + S(2))))
    rubi.add(rule207)

    pattern208 = Pattern(Integral(ArcSin(x_**p_*WC('a', S(1)))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule208 = ReplacementRule(pattern208, lambda a, p, n, x : Subst(Int(x**n*Cot(x), x), x, ArcSin(a*x**p))/p)
    rubi.add(rule208)

    pattern209 = Pattern(Integral(ArcCos(x_**p_*WC('a', S(1)))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule209 = ReplacementRule(pattern209, lambda a, p, n, x : -Subst(Int(x**n*Tan(x), x), x, ArcCos(a*x**p))/p)
    rubi.add(rule209)

    pattern210 = Pattern(Integral(ArcSin(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule210 = ReplacementRule(pattern210, lambda n, c, u, a, m, x, b : Int(u*ArcCsc(a/c + b*x**n/c)**m, x))
    rubi.add(rule210)

    pattern211 = Pattern(Integral(ArcCos(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule211 = ReplacementRule(pattern211, lambda n, c, u, a, m, x, b : Int(u*ArcSec(a/c + b*x**n/c)**m, x))
    rubi.add(rule211)

    pattern212 = Pattern(Integral(ArcSin(sqrt(x_**S(2)*WC('b', S(1)) + S(1)))**WC('n', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule212 = ReplacementRule(pattern212, lambda b, n, x : Sqrt(-b*x**S(2))*Subst(Int(ArcSin(x)**n/Sqrt(-x**S(2) + S(1)), x), x, Sqrt(b*x**S(2) + S(1)))/(b*x))
    rubi.add(rule212)

    pattern213 = Pattern(Integral(ArcCos(sqrt(x_**S(2)*WC('b', S(1)) + S(1)))**WC('n', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule213 = ReplacementRule(pattern213, lambda b, n, x : Sqrt(-b*x**S(2))*Subst(Int(ArcCos(x)**n/Sqrt(-x**S(2) + S(1)), x), x, Sqrt(b*x**S(2) + S(1)))/(b*x))
    rubi.add(rule213)

    pattern214 = Pattern(Integral(f_**(ArcSin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*WC('c', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule214 = ReplacementRule(pattern214, lambda n, c, u, a, f, x, b : Subst(Int(f**(c*x**n)*Cos(x)*ReplaceAll(u, Rule(x, -a/b + Sin(x)/b)), x), x, ArcSin(a + b*x))/b)
    rubi.add(rule214)

    pattern215 = Pattern(Integral(f_**(ArcCos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*WC('c', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule215 = ReplacementRule(pattern215, lambda n, c, u, a, f, x, b : -Subst(Int(f**(c*x**n)*ReplaceAll(u, Rule(x, -a/b + Cos(x)/b))*Sin(x), x), x, ArcCos(a + b*x))/b)
    rubi.add(rule215)

    pattern216 = Pattern(Integral(ArcSin(x_**S(2)*WC('a', S(1)) + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c: EqQ(b**S(2)*c, S(1))))
    rule216 = ReplacementRule(pattern216, lambda c, d, a, x, b : x*ArcSin(a*x**S(2) + b*Sqrt(c + d*x**S(2))) - x*Int(x*(S(2)*a*Sqrt(c + d*x**S(2)) + b*d)/(Sqrt(c + d*x**S(2))*Sqrt(a**S(2)*x**S(2) + S(2)*a*b*Sqrt(c + d*x**S(2)) + b**S(2)*d)), x)*Sqrt(a**S(2)*x**S(2) + S(2)*a*b*Sqrt(c + d*x**S(2)) + b**S(2)*d)/Sqrt(-x**S(2)*(a**S(2)*x**S(2) + S(2)*a*b*Sqrt(c + d*x**S(2)) + b**S(2)*d)))
    rubi.add(rule216)

    pattern217 = Pattern(Integral(ArcCos(x_**S(2)*WC('a', S(1)) + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c: EqQ(b**S(2)*c, S(1))))
    rule217 = ReplacementRule(pattern217, lambda c, d, a, x, b : x*ArcCos(a*x**S(2) + b*Sqrt(c + d*x**S(2))) + x*Int(x*(S(2)*a*Sqrt(c + d*x**S(2)) + b*d)/(Sqrt(c + d*x**S(2))*Sqrt(a**S(2)*x**S(2) + S(2)*a*b*Sqrt(c + d*x**S(2)) + b**S(2)*d)), x)*Sqrt(a**S(2)*x**S(2) + S(2)*a*b*Sqrt(c + d*x**S(2)) + b**S(2)*d)/Sqrt(-x**S(2)*(a**S(2)*x**S(2) + S(2)*a*b*Sqrt(c + d*x**S(2)) + b**S(2)*d)))
    rubi.add(rule217)

    pattern218 = Pattern(Integral(ArcSin(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda u, x: ~(FunctionOfExponentialQ(u, x))))
    rule218 = ReplacementRule(pattern218, lambda u, x : x*ArcSin(u) - Int(SimplifyIntegrand(x*D(u, x)/Sqrt(-u**S(2) + S(1)), x), x))
    rubi.add(rule218)

    pattern219 = Pattern(Integral(ArcCos(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda u, x: ~(FunctionOfExponentialQ(u, x))))
    rule219 = ReplacementRule(pattern219, lambda u, x : x*ArcCos(u) + Int(SimplifyIntegrand(x*D(u, x)/Sqrt(-u**S(2) + S(1)), x), x))
    rubi.add(rule219)

    pattern220 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(ArcSin(u_)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, d, u, m, x: ~(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, x: ~(FunctionOfExponentialQ(u, x))))
    rule220 = ReplacementRule(pattern220, lambda c, d, a, u, m, x, b : -b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/Sqrt(-u**S(2) + S(1)), x), x)/(d*(m + S(1))) + (a + b*ArcSin(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule220)

    pattern221 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(ArcCos(u_)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, d, u, m, x: ~(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, x: ~(FunctionOfExponentialQ(u, x))))
    rule221 = ReplacementRule(pattern221, lambda c, d, a, u, m, x, b : b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/Sqrt(-u**S(2) + S(1)), x), x)/(d*(m + S(1))) + (a + b*ArcCos(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule221)

    pattern222 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule222 = ReplacementRule(pattern222, lambda n, c, a, x, b : -b*c*n*Int(x*(a + b*ArcTan(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x) + x*(a + b*ArcTan(c*x))**n)
    rubi.add(rule222)

    pattern223 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule223 = ReplacementRule(pattern223, lambda n, c, a, x, b : b*c*n*Int(x*(a + b*ArcCot(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x) + x*(a + b*ArcCot(c*x))**n)
    rubi.add(rule223)

    pattern224 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(PositiveIntegerQ(n))))
    rule224 = ReplacementRule(pattern224, lambda n, c, a, x, b : Int((a + b*ArcTan(c*x))**n, x))
    rubi.add(rule224)

    pattern225 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(PositiveIntegerQ(n))))
    rule225 = ReplacementRule(pattern225, lambda n, c, a, x, b : Int((a + b*ArcCot(c*x))**n, x))
    rubi.add(rule225)

    pattern226 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule226 = ReplacementRule(pattern226, lambda n, c, d, a, x, e, b : b*c*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*Log(S(2)*d/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x)/e - (a + b*ArcTan(c*x))**n*Log(S(2)*d/(d + e*x))/e)
    rubi.add(rule226)

    pattern227 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule227 = ReplacementRule(pattern227, lambda n, c, d, a, x, e, b : -b*c*n*Int((a + b*ArcCot(c*x))**(n + S(-1))*Log(S(2)*d/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x)/e - (a + b*ArcCot(c*x))**n*Log(S(2)*d/(d + e*x))/e)
    rubi.add(rule227)

    pattern228 = Pattern(Integral(ArcTan(x_*WC('c', S(1)))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: PositiveQ(ImaginaryI*c*d/e + S(1))), CustomConstraint(lambda d, c, e: NegativeQ(ImaginaryI*c*d/e + S(-1))))
    rule228 = ReplacementRule(pattern228, lambda d, c, x, e : ImaginaryI*PolyLog(S(2), Simp(ImaginaryI*c*(d + e*x)/(ImaginaryI*c*d - e), x))/(S(2)*e) - ImaginaryI*PolyLog(S(2), Simp(ImaginaryI*c*(d + e*x)/(ImaginaryI*c*d + e), x))/(S(2)*e) - ArcTan(c*d/e)*Log(d + e*x)/e)
    rubi.add(rule228)

    pattern229 = Pattern(Integral(ArcTan(x_*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule229 = ReplacementRule(pattern229, lambda d, c, x, e : ImaginaryI*Int(Log(-ImaginaryI*c*x + S(1))/(d + e*x), x)/S(2) - ImaginaryI*Int(Log(ImaginaryI*c*x + S(1))/(d + e*x), x)/S(2))
    rubi.add(rule229)

    pattern230 = Pattern(Integral(ArcCot(x_*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule230 = ReplacementRule(pattern230, lambda d, c, x, e : ImaginaryI*Int(Log(-ImaginaryI/(c*x) + S(1))/(d + e*x), x)/S(2) - ImaginaryI*Int(Log(ImaginaryI/(c*x) + S(1))/(d + e*x), x)/S(2))
    rubi.add(rule230)

    pattern231 = Pattern(Integral((a_ + ArcTan(x_*WC('c', S(1)))*WC('b', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule231 = ReplacementRule(pattern231, lambda c, d, a, x, e, b : a*Log(RemoveContent(d + e*x, x))/e + b*Int(ArcTan(c*x)/(d + e*x), x))
    rubi.add(rule231)

    pattern232 = Pattern(Integral((a_ + ArcCot(x_*WC('c', S(1)))*WC('b', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule232 = ReplacementRule(pattern232, lambda c, d, a, x, e, b : a*Log(RemoveContent(d + e*x, x))/e + b*Int(ArcCot(c*x)/(d + e*x), x))
    rubi.add(rule232)

    pattern233 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule233 = ReplacementRule(pattern233, lambda p, c, d, a, x, e, b : -b*c*Int((d + e*x)**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x)/(e*(p + S(1))) + (a + b*ArcTan(c*x))*(d + e*x)**(p + S(1))/(e*(p + S(1))))
    rubi.add(rule233)

    pattern234 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule234 = ReplacementRule(pattern234, lambda p, c, d, a, x, e, b : b*c*Int((d + e*x)**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x)/(e*(p + S(1))) + (a + b*ArcCot(c*x))*(d + e*x)**(p + S(1))/(e*(p + S(1))))
    rubi.add(rule234)

    pattern235 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule235 = ReplacementRule(pattern235, lambda n, c, a, x, b : -S(2)*b*c*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*ArcTanh(-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))/(c**S(2)*x**S(2) + S(1)), x) + S(2)*(a + b*ArcTan(c*x))**n*ArcTanh(-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1)))
    rubi.add(rule235)

    pattern236 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule236 = ReplacementRule(pattern236, lambda n, c, a, x, b : S(2)*b*c*n*Int((a + b*ArcCot(c*x))**(n + S(-1))*ArcCoth(-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))/(c**S(2)*x**S(2) + S(1)), x) + S(2)*(a + b*ArcCot(c*x))**n*ArcCoth(-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1)))
    rubi.add(rule236)

    pattern237 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule237 = ReplacementRule(pattern237, lambda n, c, a, m, x, b : -b*c*n*Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcTan(c*x))**n/(m + S(1)))
    rubi.add(rule237)

    pattern238 = Pattern(Integral(x_**WC('m', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule238 = ReplacementRule(pattern238, lambda n, c, a, m, x, b : b*c*n*Int(x**(m + S(1))*(a + b*ArcCot(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcCot(c*x))**n/(m + S(1)))
    rubi.add(rule238)

    pattern239 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)))
    rule239 = ReplacementRule(pattern239, lambda p, n, c, d, a, x, e, b : Int(ExpandIntegrand((a + b*ArcTan(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule239)

    pattern240 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)))
    rule240 = ReplacementRule(pattern240, lambda p, n, c, d, a, x, e, b : Int(ExpandIntegrand((a + b*ArcCot(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule240)

    pattern241 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule241 = ReplacementRule(pattern241, lambda p, n, c, d, a, x, e, b : Int((a + b*ArcTan(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule241)

    pattern242 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule242 = ReplacementRule(pattern242, lambda p, n, c, d, a, x, e, b : Int((a + b*ArcCot(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule242)

    pattern243 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule243 = ReplacementRule(pattern243, lambda n, c, d, a, m, x, e, b : -d*Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**n/(d + e*x), x)/e + Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**n, x)/e)
    rubi.add(rule243)

    pattern244 = Pattern(Integral(x_**WC('m', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule244 = ReplacementRule(pattern244, lambda n, c, d, a, m, x, e, b : -d*Int(x**(m + S(-1))*(a + b*ArcCot(c*x))**n/(d + e*x), x)/e + Int(x**(m + S(-1))*(a + b*ArcCot(c*x))**n, x)/e)
    rubi.add(rule244)

    pattern245 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*(d_ + x_*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule245 = ReplacementRule(pattern245, lambda n, c, d, a, x, e, b : -b*c*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*Log(S(2)*e*x/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x)/d + (a + b*ArcTan(c*x))**n*Log(S(2)*e*x/(d + e*x))/d)
    rubi.add(rule245)

    pattern246 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*(d_ + x_*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule246 = ReplacementRule(pattern246, lambda n, c, d, a, x, e, b : b*c*n*Int((a + b*ArcCot(c*x))**(n + S(-1))*Log(S(2)*e*x/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x)/d + (a + b*ArcCot(c*x))**n*Log(S(2)*e*x/(d + e*x))/d)
    rubi.add(rule246)

    pattern247 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule247 = ReplacementRule(pattern247, lambda n, c, d, a, x, m, e, b : -e*Int(x**(m + S(1))*(a + b*ArcTan(c*x))**n/(d + e*x), x)/d + Int(x**m*(a + b*ArcTan(c*x))**n, x)/d)
    rubi.add(rule247)

    pattern248 = Pattern(Integral(x_**m_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule248 = ReplacementRule(pattern248, lambda n, c, d, a, x, m, e, b : -e*Int(x**(m + S(1))*(a + b*ArcCot(c*x))**n/(d + e*x), x)/d + Int(x**m*(a + b*ArcCot(c*x))**n, x)/d)
    rubi.add(rule248)

    pattern249 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda a, p, m: IntegerQ(m) | NonzeroQ(a) | Greater(p, S(0))))
    rule249 = ReplacementRule(pattern249, lambda p, n, c, d, a, m, x, e, b : Int(ExpandIntegrand(x**m*(a + b*ArcTan(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule249)

    pattern250 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda a, p, m: IntegerQ(m) | NonzeroQ(a) | Greater(p, S(0))))
    rule250 = ReplacementRule(pattern250, lambda p, n, c, d, a, m, x, e, b : Int(ExpandIntegrand(x**m*(a + b*ArcCot(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule250)

    pattern251 = Pattern(Integral(x_**WC('m', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule251 = ReplacementRule(pattern251, lambda p, n, c, d, a, m, x, e, b : Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule251)

    pattern252 = Pattern(Integral(x_**WC('m', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule252 = ReplacementRule(pattern252, lambda p, n, c, d, a, m, x, e, b : Int(x**m*(a + b*ArcCot(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule252)

    pattern253 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule253 = ReplacementRule(pattern253, lambda p, c, d, a, x, e, b : -b*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*ArcTan(c*x))*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule253)

    pattern254 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule254 = ReplacementRule(pattern254, lambda p, c, d, a, x, e, b : b*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*ArcCot(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*ArcCot(c*x))*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule254)

    pattern255 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule255 = ReplacementRule(pattern255, lambda p, n, c, d, a, x, e, b : b**S(2)*d*n*(n + S(-1))*Int((a + b*ArcTan(c*x))**(n + S(-2))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p*(S(2)*p + S(1))) - b*n*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule255)

    pattern256 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule256 = ReplacementRule(pattern256, lambda p, n, c, d, a, x, e, b : b**S(2)*d*n*(n + S(-1))*Int((a + b*ArcCot(c*x))**(n + S(-2))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p*(S(2)*p + S(1))) + b*n*(a + b*ArcCot(c*x))**(n + S(-1))*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*ArcCot(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule256)

    pattern257 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)))
    rule257 = ReplacementRule(pattern257, lambda c, d, a, x, e, b : Log(RemoveContent(a + b*ArcTan(c*x), x))/(b*c*d))
    rubi.add(rule257)

    pattern258 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)))
    rule258 = ReplacementRule(pattern258, lambda c, d, a, x, e, b : -Log(RemoveContent(a + b*ArcCot(c*x), x))/(b*c*d))
    rubi.add(rule258)

    pattern259 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule259 = ReplacementRule(pattern259, lambda n, c, d, a, x, e, b : (a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule259)

    pattern260 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule260 = ReplacementRule(pattern260, lambda n, c, d, a, x, e, b : -(a + b*ArcCot(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule260)

    pattern261 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule261 = ReplacementRule(pattern261, lambda c, d, a, x, e, b : ImaginaryI*b*PolyLog(S(2), -ImaginaryI*Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/(c*Sqrt(d)) - ImaginaryI*b*PolyLog(S(2), ImaginaryI*Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/(c*Sqrt(d)) - S(2)*ImaginaryI*(a + b*ArcTan(c*x))*ArcTan(Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/(c*Sqrt(d)))
    rubi.add(rule261)

    pattern262 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule262 = ReplacementRule(pattern262, lambda c, d, a, x, e, b : -ImaginaryI*b*PolyLog(S(2), -ImaginaryI*Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/(c*Sqrt(d)) + ImaginaryI*b*PolyLog(S(2), ImaginaryI*Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/(c*Sqrt(d)) - S(2)*ImaginaryI*(a + b*ArcCot(c*x))*ArcTan(Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/(c*Sqrt(d)))
    rubi.add(rule262)

    pattern263 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule263 = ReplacementRule(pattern263, lambda n, c, d, a, x, e, b : Subst(Int((a + b*x)**n*Sec(x), x), x, ArcTan(c*x))/(c*Sqrt(d)))
    rubi.add(rule263)

    pattern264 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule264 = ReplacementRule(pattern264, lambda n, c, d, a, x, e, b : -x*Sqrt(S(1) + S(1)/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*Csc(x), x), x, ArcCot(c*x))/Sqrt(d + e*x**S(2)))
    rubi.add(rule264)

    pattern265 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: ~(PositiveQ(d))))
    rule265 = ReplacementRule(pattern265, lambda n, c, d, a, x, e, b : Int((a + b*ArcTan(c*x))**n/Sqrt(c**S(2)*x**S(2) + S(1)), x)*Sqrt(c**S(2)*x**S(2) + S(1))/Sqrt(d + e*x**S(2)))
    rubi.add(rule265)

    pattern266 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: ~(PositiveQ(d))))
    rule266 = ReplacementRule(pattern266, lambda n, c, d, a, x, e, b : Int((a + b*ArcCot(c*x))**n/Sqrt(c**S(2)*x**S(2) + S(1)), x)*Sqrt(c**S(2)*x**S(2) + S(1))/Sqrt(d + e*x**S(2)))
    rubi.add(rule266)

    pattern267 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule267 = ReplacementRule(pattern267, lambda n, c, d, a, x, e, b : -b*c*n*Int(x*(a + b*ArcTan(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/S(2) + x*(a + b*ArcTan(c*x))**n/(S(2)*d*(d + e*x**S(2))) + (a + b*ArcTan(c*x))**(n + S(1))/(S(2)*b*c*d**S(2)*(n + S(1))))
    rubi.add(rule267)

    pattern268 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule268 = ReplacementRule(pattern268, lambda n, c, d, a, x, e, b : b*c*n*Int(x*(a + b*ArcCot(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/S(2) + x*(a + b*ArcCot(c*x))**n/(S(2)*d*(d + e*x**S(2))) - (a + b*ArcCot(c*x))**(n + S(1))/(S(2)*b*c*d**S(2)*(n + S(1))))
    rubi.add(rule268)

    pattern269 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)))
    rule269 = ReplacementRule(pattern269, lambda c, d, a, x, e, b : b/(c*d*Sqrt(d + e*x**S(2))) + x*(a + b*ArcTan(c*x))/(d*Sqrt(d + e*x**S(2))))
    rubi.add(rule269)

    pattern270 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)))
    rule270 = ReplacementRule(pattern270, lambda c, d, a, x, e, b : -b/(c*d*Sqrt(d + e*x**S(2))) + x*(a + b*ArcCot(c*x))/(d*Sqrt(d + e*x**S(2))))
    rubi.add(rule270)

    pattern271 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule271 = ReplacementRule(pattern271, lambda p, c, d, a, x, e, b : b*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule271)

    pattern272 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule272 = ReplacementRule(pattern272, lambda p, c, d, a, x, e, b : -b*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*ArcCot(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*ArcCot(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule272)

    pattern273 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule273 = ReplacementRule(pattern273, lambda n, c, d, a, x, e, b : -b**S(2)*n*(n + S(-1))*Int((a + b*ArcTan(c*x))**(n + S(-2))/(d + e*x**S(2))**(S(3)/2), x) + b*n*(a + b*ArcTan(c*x))**(n + S(-1))/(c*d*Sqrt(d + e*x**S(2))) + x*(a + b*ArcTan(c*x))**n/(d*Sqrt(d + e*x**S(2))))
    rubi.add(rule273)

    pattern274 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule274 = ReplacementRule(pattern274, lambda n, c, d, a, x, e, b : -b**S(2)*n*(n + S(-1))*Int((a + b*ArcCot(c*x))**(n + S(-2))/(d + e*x**S(2))**(S(3)/2), x) - b*n*(a + b*ArcCot(c*x))**(n + S(-1))/(c*d*Sqrt(d + e*x**S(2))) + x*(a + b*ArcCot(c*x))**n/(d*Sqrt(d + e*x**S(2))))
    rubi.add(rule274)

    pattern275 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule275 = ReplacementRule(pattern275, lambda p, n, c, d, a, x, e, b : -b**S(2)*n*(n + S(-1))*Int((a + b*ArcTan(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/(S(4)*(p + S(1))**S(2)) + b*n*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule275)

    pattern276 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule276 = ReplacementRule(pattern276, lambda p, n, c, d, a, x, e, b : -b**S(2)*n*(n + S(-1))*Int((a + b*ArcCot(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/(S(4)*(p + S(1))**S(2)) - b*n*(a + b*ArcCot(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*ArcCot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule276)

    pattern277 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))))
    rule277 = ReplacementRule(pattern277, lambda p, n, c, d, a, x, e, b : -S(2)*c*(p + S(1))*Int(x*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) + (a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule277)

    pattern278 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))))
    rule278 = ReplacementRule(pattern278, lambda p, n, c, d, a, x, e, b : S(2)*c*(p + S(1))*Int(x*(a + b*ArcCot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) - (a + b*ArcCot(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule278)

    pattern279 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda d, p: IntegerQ(p) | PositiveQ(d)))
    rule279 = ReplacementRule(pattern279, lambda p, n, c, d, a, x, e, b : d**p*Subst(Int((a + b*x)**n*Cos(x)**(-S(2)*p + S(-2)), x), x, ArcTan(c*x))/c)
    rubi.add(rule279)

    pattern280 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda d, p: ~(IntegerQ(p) | PositiveQ(d))))
    rule280 = ReplacementRule(pattern280, lambda p, n, c, d, a, x, e, b : d**(p + S(1)/2)*Int((a + b*ArcTan(c*x))**n*(c**S(2)*x**S(2) + S(1))**p, x)*Sqrt(c**S(2)*x**S(2) + S(1))/Sqrt(d + e*x**S(2)))
    rubi.add(rule280)

    pattern281 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule281 = ReplacementRule(pattern281, lambda p, n, c, d, a, x, e, b : -d**p*Subst(Int((a + b*x)**n*Sin(x)**(-S(2)*p + S(-2)), x), x, ArcCot(c*x))/c)
    rubi.add(rule281)

    pattern282 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule282 = ReplacementRule(pattern282, lambda p, n, c, d, a, x, e, b : -d**(p + S(1)/2)*x*Sqrt((c**S(2)*x**S(2) + S(1))/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*Sin(x)**(-S(2)*p + S(-2)), x), x, ArcCot(c*x))/Sqrt(d + e*x**S(2)))
    rubi.add(rule282)

    pattern283 = Pattern(Integral(ArcTan(x_*WC('c', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule283 = ReplacementRule(pattern283, lambda d, c, x, e : ImaginaryI*Int(Log(-ImaginaryI*c*x + S(1))/(d + e*x**S(2)), x)/S(2) - ImaginaryI*Int(Log(ImaginaryI*c*x + S(1))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule283)

    pattern284 = Pattern(Integral(ArcCot(x_*WC('c', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule284 = ReplacementRule(pattern284, lambda d, c, x, e : ImaginaryI*Int(Log(-ImaginaryI/(c*x) + S(1))/(d + e*x**S(2)), x)/S(2) - ImaginaryI*Int(Log(ImaginaryI/(c*x) + S(1))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule284)

    pattern285 = Pattern(Integral((a_ + ArcTan(x_*WC('c', S(1)))*WC('b', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule285 = ReplacementRule(pattern285, lambda c, d, a, x, e, b : a*Int(1/(d + e*x**S(2)), x) + b*Int(ArcTan(c*x)/(d + e*x**S(2)), x))
    rubi.add(rule285)

    pattern286 = Pattern(Integral((a_ + ArcCot(x_*WC('c', S(1)))*WC('b', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule286 = ReplacementRule(pattern286, lambda c, d, a, x, e, b : a*Int(1/(d + e*x**S(2)), x) + b*Int(ArcCot(c*x)/(d + e*x**S(2)), x))
    rubi.add(rule286)

    pattern287 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p) | NegativeIntegerQ(p + S(1)/2)))
    rule287 = ReplacementRule(pattern287, lambda p, c, d, a, x, e, b : With(List(Set(u, IntHide((d + e*x**S(2))**p, x))), -b*c*Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcTan(c*x), u, x)))
    rubi.add(rule287)

    pattern288 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p) | NegativeIntegerQ(p + S(1)/2)))
    rule288 = ReplacementRule(pattern288, lambda p, c, d, a, x, e, b : With(List(Set(u, IntHide((d + e*x**S(2))**p, x))), b*c*Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcCot(c*x), u, x)))
    rubi.add(rule288)

    pattern289 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule289 = ReplacementRule(pattern289, lambda p, n, c, d, a, x, e, b : Int(ExpandIntegrand((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule289)

    pattern290 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule290 = ReplacementRule(pattern290, lambda p, n, c, d, a, x, e, b : Int(ExpandIntegrand((a + b*ArcCot(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule290)

    pattern291 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule291 = ReplacementRule(pattern291, lambda p, n, c, d, a, x, e, b : Int((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule291)

    pattern292 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule292 = ReplacementRule(pattern292, lambda p, n, c, d, a, x, e, b : Int((a + b*ArcCot(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule292)

    pattern293 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule293 = ReplacementRule(pattern293, lambda n, c, d, a, x, m, e, b : -d*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n, x)/e)
    rubi.add(rule293)

    pattern294 = Pattern(Integral(x_**m_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule294 = ReplacementRule(pattern294, lambda n, c, d, a, x, m, e, b : -d*Int(x**(m + S(-2))*(a + b*ArcCot(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*ArcCot(c*x))**n, x)/e)
    rubi.add(rule294)

    pattern295 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule295 = ReplacementRule(pattern295, lambda n, c, d, a, x, m, e, b : -e*Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*ArcTan(c*x))**n, x)/d)
    rubi.add(rule295)

    pattern296 = Pattern(Integral(x_**m_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule296 = ReplacementRule(pattern296, lambda n, c, d, a, x, m, e, b : -e*Int(x**(m + S(2))*(a + b*ArcCot(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*ArcCot(c*x))**n, x)/d)
    rubi.add(rule296)

    pattern297 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule297 = ReplacementRule(pattern297, lambda n, c, d, a, x, e, b : -ImaginaryI*(a + b*ArcTan(c*x))**(n + S(1))/(b*e*(n + S(1))) - Int((a + b*ArcTan(c*x))**n/(ImaginaryI - c*x), x)/(c*d))
    rubi.add(rule297)

    pattern298 = Pattern(Integral(x_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule298 = ReplacementRule(pattern298, lambda n, c, d, a, x, e, b : ImaginaryI*(a + b*ArcCot(c*x))**(n + S(1))/(b*e*(n + S(1))) - Int((a + b*ArcCot(c*x))**n/(ImaginaryI - c*x), x)/(c*d))
    rubi.add(rule298)

    pattern299 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: ~(PositiveIntegerQ(n))), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule299 = ReplacementRule(pattern299, lambda n, c, d, a, x, e, b : x*(a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(n + S(1))) - Int((a + b*ArcTan(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))))
    rubi.add(rule299)

    pattern300 = Pattern(Integral(x_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: ~(PositiveIntegerQ(n))), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule300 = ReplacementRule(pattern300, lambda n, c, d, a, x, e, b : -x*(a + b*ArcCot(c*x))**(n + S(1))/(b*c*d*(n + S(1))) + Int((a + b*ArcCot(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))))
    rubi.add(rule300)

    pattern301 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule301 = ReplacementRule(pattern301, lambda n, c, d, a, x, m, e, b : -d*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n, x)/e)
    rubi.add(rule301)

    pattern302 = Pattern(Integral(x_**m_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule302 = ReplacementRule(pattern302, lambda n, c, d, a, x, m, e, b : -d*Int(x**(m + S(-2))*(a + b*ArcCot(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*ArcCot(c*x))**n, x)/e)
    rubi.add(rule302)

    pattern303 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule303 = ReplacementRule(pattern303, lambda n, c, d, a, x, e, b : ImaginaryI*Int((a + b*ArcTan(c*x))**n/(x*(ImaginaryI + c*x)), x)/d - ImaginaryI*(a + b*ArcTan(c*x))**(n + S(1))/(b*d*(n + S(1))))
    rubi.add(rule303)

    pattern304 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule304 = ReplacementRule(pattern304, lambda n, c, d, a, x, e, b : ImaginaryI*Int((a + b*ArcCot(c*x))**n/(x*(ImaginaryI + c*x)), x)/d + ImaginaryI*(a + b*ArcCot(c*x))**(n + S(1))/(b*d*(n + S(1))))
    rubi.add(rule304)

    pattern305 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule305 = ReplacementRule(pattern305, lambda n, c, d, a, x, m, e, b : -e*Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*ArcTan(c*x))**n, x)/d)
    rubi.add(rule305)

    pattern306 = Pattern(Integral(x_**m_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule306 = ReplacementRule(pattern306, lambda n, c, d, a, x, m, e, b : -e*Int(x**(m + S(2))*(a + b*ArcCot(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*ArcCot(c*x))**n, x)/d)
    rubi.add(rule306)

    pattern307 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule307 = ReplacementRule(pattern307, lambda n, c, d, a, x, m, e, b : -m*Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))) + x**m*(a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule307)

    pattern308 = Pattern(Integral(x_**m_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule308 = ReplacementRule(pattern308, lambda n, c, d, a, x, m, e, b : m*Int(x**(m + S(-1))*(a + b*ArcCot(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))) - x**m*(a + b*ArcCot(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule308)

    pattern309 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda a, m: ~(NonzeroQ(a) & Equal(m, S(1)))))
    rule309 = ReplacementRule(pattern309, lambda c, d, a, m, x, e, b : Int(ExpandIntegrand(a + b*ArcTan(c*x), x**m/(d + e*x**S(2)), x), x))
    rubi.add(rule309)

    pattern310 = Pattern(Integral(x_**WC('m', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda a, m: ~(NonzeroQ(a) & Equal(m, S(1)))))
    rule310 = ReplacementRule(pattern310, lambda c, d, a, m, x, e, b : Int(ExpandIntegrand(a + b*ArcCot(c*x), x**m/(d + e*x**S(2)), x), x))
    rubi.add(rule310)

    pattern311 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule311 = ReplacementRule(pattern311, lambda p, n, c, d, a, x, e, b : -b*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(S(2)*c*(p + S(1))) + (a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule311)

    pattern312 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule312 = ReplacementRule(pattern312, lambda p, n, c, d, a, x, e, b : b*n*Int((a + b*ArcCot(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(S(2)*c*(p + S(1))) + (a + b*ArcCot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule312)

    pattern313 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule313 = ReplacementRule(pattern313, lambda n, c, d, a, x, e, b : x*(a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(d + e*x**S(2))*(n + S(1))) - S(4)*Int(x*(a + b*ArcTan(c*x))**(n + S(2))/(d + e*x**S(2))**S(2), x)/(b**S(2)*(n + S(1))*(n + S(2))) - (a + b*ArcTan(c*x))**(n + S(2))*(-c**S(2)*x**S(2) + S(1))/(b**S(2)*e*(d + e*x**S(2))*(n + S(1))*(n + S(2))))
    rubi.add(rule313)

    pattern314 = Pattern(Integral(x_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule314 = ReplacementRule(pattern314, lambda n, c, d, a, x, e, b : -x*(a + b*ArcCot(c*x))**(n + S(1))/(b*c*d*(d + e*x**S(2))*(n + S(1))) - S(4)*Int(x*(a + b*ArcCot(c*x))**(n + S(2))/(d + e*x**S(2))**S(2), x)/(b**S(2)*(n + S(1))*(n + S(2))) - (a + b*ArcCot(c*x))**(n + S(2))*(-c**S(2)*x**S(2) + S(1))/(b**S(2)*e*(d + e*x**S(2))*(n + S(1))*(n + S(2))))
    rubi.add(rule314)

    pattern315 = Pattern(Integral(x_**S(2)*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-5)/2)))
    rule315 = ReplacementRule(pattern315, lambda p, c, d, a, x, e, b : -b*(d + e*x**S(2))**(p + S(1))/(S(4)*c**S(3)*d*(p + S(1))**S(2)) + x*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*c**S(2)*d*(p + S(1))) - Int((a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*c**S(2)*d*(p + S(1))))
    rubi.add(rule315)

    pattern316 = Pattern(Integral(x_**S(2)*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-5)/2)))
    rule316 = ReplacementRule(pattern316, lambda p, c, d, a, x, e, b : b*(d + e*x**S(2))**(p + S(1))/(S(4)*c**S(3)*d*(p + S(1))**S(2)) + x*(a + b*ArcCot(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*c**S(2)*d*(p + S(1))) - Int((a + b*ArcCot(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*c**S(2)*d*(p + S(1))))
    rubi.add(rule316)

    pattern317 = Pattern(Integral(x_**S(2)*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule317 = ReplacementRule(pattern317, lambda n, c, d, a, x, e, b : b*n*Int(x*(a + b*ArcTan(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/(S(2)*c) - x*(a + b*ArcTan(c*x))**n/(S(2)*c**S(2)*d*(d + e*x**S(2))) + (a + b*ArcTan(c*x))**(n + S(1))/(S(2)*b*c**S(3)*d**S(2)*(n + S(1))))
    rubi.add(rule317)

    pattern318 = Pattern(Integral(x_**S(2)*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule318 = ReplacementRule(pattern318, lambda n, c, d, a, x, e, b : -b*n*Int(x*(a + b*ArcCot(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/(S(2)*c) - x*(a + b*ArcCot(c*x))**n/(S(2)*c**S(2)*d*(d + e*x**S(2))) - (a + b*ArcCot(c*x))**(n + S(1))/(S(2)*b*c**S(3)*d**S(2)*(n + S(1))))
    rubi.add(rule318)

    pattern319 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule319 = ReplacementRule(pattern319, lambda p, c, d, a, x, m, e, b : b*x**m*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) - x**(m + S(-1))*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) + (m + S(-1))*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule319)

    pattern320 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule320 = ReplacementRule(pattern320, lambda p, c, d, a, x, m, e, b : -b*x**m*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) - x**(m + S(-1))*(a + b*ArcCot(c*x))*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) + (m + S(-1))*Int(x**(m + S(-2))*(a + b*ArcCot(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule320)

    pattern321 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule321 = ReplacementRule(pattern321, lambda p, n, c, d, a, x, m, e, b : -b**S(2)*n*(n + S(-1))*Int(x**m*(a + b*ArcTan(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/m**S(2) + b*n*x**m*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) - x**(m + S(-1))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) + (m + S(-1))*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule321)

    pattern322 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule322 = ReplacementRule(pattern322, lambda p, n, c, d, a, x, m, e, b : -b**S(2)*n*(n + S(-1))*Int(x**m*(a + b*ArcCot(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/m**S(2) - b*n*x**m*(a + b*ArcCot(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) - x**(m + S(-1))*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) + (m + S(-1))*Int(x**(m + S(-2))*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule322)

    pattern323 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule323 = ReplacementRule(pattern323, lambda p, n, c, d, a, m, x, e, b : -m*Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) + x**m*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule323)

    pattern324 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule324 = ReplacementRule(pattern324, lambda p, n, c, d, a, m, x, e, b : m*Int(x**(m + S(-1))*(a + b*ArcCot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) - x**m*(a + b*ArcCot(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule324)

    pattern325 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule325 = ReplacementRule(pattern325, lambda p, n, c, d, a, m, x, e, b : -b*c*n*Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*(m + S(1))))
    rubi.add(rule325)

    pattern326 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule326 = ReplacementRule(pattern326, lambda p, n, c, d, a, m, x, e, b : b*c*n*Int(x**(m + S(1))*(a + b*ArcCot(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*(m + S(1))))
    rubi.add(rule326)

    pattern327 = Pattern(Integral(x_**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: NonzeroQ(m + S(2))))
    rule327 = ReplacementRule(pattern327, lambda c, d, a, x, m, e, b : -b*c*d*Int(x**(m + S(1))/Sqrt(d + e*x**S(2)), x)/(m + S(2)) + d*Int(x**m*(a + b*ArcTan(c*x))/Sqrt(d + e*x**S(2)), x)/(m + S(2)) + x**(m + S(1))*(a + b*ArcTan(c*x))*Sqrt(d + e*x**S(2))/(m + S(2)))
    rubi.add(rule327)

    pattern328 = Pattern(Integral(x_**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: NonzeroQ(m + S(2))))
    rule328 = ReplacementRule(pattern328, lambda c, d, a, x, m, e, b : b*c*d*Int(x**(m + S(1))/Sqrt(d + e*x**S(2)), x)/(m + S(2)) + d*Int(x**m*(a + b*ArcCot(c*x))/Sqrt(d + e*x**S(2)), x)/(m + S(2)) + x**(m + S(1))*(a + b*ArcCot(c*x))*Sqrt(d + e*x**S(2))/(m + S(2)))
    rubi.add(rule328)

    pattern329 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p: Greater(p, S(1))))
    rule329 = ReplacementRule(pattern329, lambda p, n, c, d, a, x, m, e, b : Int(ExpandIntegrand(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule329)

    pattern330 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p: Greater(p, S(1))))
    rule330 = ReplacementRule(pattern330, lambda p, n, c, d, a, x, m, e, b : Int(ExpandIntegrand(x**m*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule330)

    pattern331 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, n, m: RationalQ(m) | (IntegerQ(p) & EqQ(n, S(1)))))
    rule331 = ReplacementRule(pattern331, lambda p, n, c, d, a, x, m, e, b : c**S(2)*d*Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x) + d*Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x))
    rubi.add(rule331)

    pattern332 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, n, m: RationalQ(m) | (IntegerQ(p) & EqQ(n, S(1)))))
    rule332 = ReplacementRule(pattern332, lambda p, n, c, d, a, x, m, e, b : c**S(2)*d*Int(x**(m + S(2))*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x) + d*Int(x**m*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x))
    rubi.add(rule332)

    pattern333 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule333 = ReplacementRule(pattern333, lambda n, c, d, a, x, m, e, b : -b*n*Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(-1))/Sqrt(d + e*x**S(2)), x)/(c*m) - (m + S(-1))*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n/Sqrt(d + e*x**S(2)), x)/(c**S(2)*m) + x**(m + S(-1))*(a + b*ArcTan(c*x))**n*Sqrt(d + e*x**S(2))/(c**S(2)*d*m))
    rubi.add(rule333)

    pattern334 = Pattern(Integral(x_**m_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule334 = ReplacementRule(pattern334, lambda n, c, d, a, x, m, e, b : b*n*Int(x**(m + S(-1))*(a + b*ArcCot(c*x))**(n + S(-1))/Sqrt(d + e*x**S(2)), x)/(c*m) - (m + S(-1))*Int(x**(m + S(-2))*(a + b*ArcCot(c*x))**n/Sqrt(d + e*x**S(2)), x)/(c**S(2)*m) + x**(m + S(-1))*(a + b*ArcCot(c*x))**n*Sqrt(d + e*x**S(2))/(c**S(2)*d*m))
    rubi.add(rule334)

    pattern335 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule335 = ReplacementRule(pattern335, lambda c, d, a, x, e, b : ImaginaryI*b*PolyLog(S(2), -Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/Sqrt(d) - ImaginaryI*b*PolyLog(S(2), Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/Sqrt(d) - S(2)*(a + b*ArcTan(c*x))*ArcTanh(Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/Sqrt(d))
    rubi.add(rule335)

    pattern336 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule336 = ReplacementRule(pattern336, lambda c, d, a, x, e, b : -ImaginaryI*b*PolyLog(S(2), -Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/Sqrt(d) + ImaginaryI*b*PolyLog(S(2), Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/Sqrt(d) - S(2)*(a + b*ArcCot(c*x))*ArcTanh(Sqrt(ImaginaryI*c*x + S(1))/Sqrt(-ImaginaryI*c*x + S(1)))/Sqrt(d))
    rubi.add(rule336)

    pattern337 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule337 = ReplacementRule(pattern337, lambda n, c, d, a, x, e, b : Subst(Int((a + b*x)**n*Csc(x), x), x, ArcTan(c*x))/Sqrt(d))
    rubi.add(rule337)

    pattern338 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule338 = ReplacementRule(pattern338, lambda n, c, d, a, x, e, b : -c*x*Sqrt(S(1) + S(1)/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*Sec(x), x), x, ArcCot(c*x))/Sqrt(d + e*x**S(2)))
    rubi.add(rule338)

    pattern339 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: ~(PositiveQ(d))))
    rule339 = ReplacementRule(pattern339, lambda n, c, d, a, x, e, b : Int((a + b*ArcTan(c*x))**n/(x*Sqrt(c**S(2)*x**S(2) + S(1))), x)*Sqrt(c**S(2)*x**S(2) + S(1))/Sqrt(d + e*x**S(2)))
    rubi.add(rule339)

    pattern340 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: ~(PositiveQ(d))))
    rule340 = ReplacementRule(pattern340, lambda n, c, d, a, x, e, b : Int((a + b*ArcCot(c*x))**n/(x*Sqrt(c**S(2)*x**S(2) + S(1))), x)*Sqrt(c**S(2)*x**S(2) + S(1))/Sqrt(d + e*x**S(2)))
    rubi.add(rule340)

    pattern341 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule341 = ReplacementRule(pattern341, lambda n, c, d, a, x, e, b : b*c*n*Int((a + b*ArcTan(c*x))**(n + S(-1))/(x*Sqrt(d + e*x**S(2))), x) - (a + b*ArcTan(c*x))**n*Sqrt(d + e*x**S(2))/(d*x))
    rubi.add(rule341)

    pattern342 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule342 = ReplacementRule(pattern342, lambda n, c, d, a, x, e, b : -b*c*n*Int((a + b*ArcCot(c*x))**(n + S(-1))/(x*Sqrt(d + e*x**S(2))), x) - (a + b*ArcCot(c*x))**n*Sqrt(d + e*x**S(2))/(d*x))
    rubi.add(rule342)

    pattern343 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: Unequal(m, S(-2))))
    rule343 = ReplacementRule(pattern343, lambda n, c, d, a, x, m, e, b : -b*c*n*Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(-1))/Sqrt(d + e*x**S(2)), x)/(m + S(1)) - c**S(2)*(m + S(2))*Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n/Sqrt(d + e*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcTan(c*x))**n*Sqrt(d + e*x**S(2))/(d*(m + S(1))))
    rubi.add(rule343)

    pattern344 = Pattern(Integral(x_**m_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: Unequal(m, S(-2))))
    rule344 = ReplacementRule(pattern344, lambda n, c, d, a, x, m, e, b : b*c*n*Int(x**(m + S(1))*(a + b*ArcCot(c*x))**(n + S(-1))/Sqrt(d + e*x**S(2)), x)/(m + S(1)) - c**S(2)*(m + S(2))*Int(x**(m + S(2))*(a + b*ArcCot(c*x))**n/Sqrt(d + e*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcCot(c*x))**n*Sqrt(d + e*x**S(2))/(d*(m + S(1))))
    rubi.add(rule344)

    pattern345 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule345 = ReplacementRule(pattern345, lambda p, n, c, d, a, x, m, e, b : -d*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x)/e + Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/e)
    rubi.add(rule345)

    pattern346 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule346 = ReplacementRule(pattern346, lambda p, n, c, d, a, x, m, e, b : -d*Int(x**(m + S(-2))*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**p, x)/e + Int(x**(m + S(-2))*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/e)
    rubi.add(rule346)

    pattern347 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Less(m, S(0))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule347 = ReplacementRule(pattern347, lambda p, n, c, d, a, x, m, e, b : -e*Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x)/d + Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/d)
    rubi.add(rule347)

    pattern348 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Less(m, S(0))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule348 = ReplacementRule(pattern348, lambda p, n, c, d, a, x, m, e, b : -e*Int(x**(m + S(2))*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**p, x)/d + Int(x**m*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/d)
    rubi.add(rule348)

    pattern349 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n, m: RationalQ(m, n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(2))))
    rule349 = ReplacementRule(pattern349, lambda p, n, c, d, a, m, x, e, b : -c*(m + S(2)*p + S(2))*Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) - m*Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) + x**m*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule349)

    pattern350 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n, m: RationalQ(m, n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(2))))
    rule350 = ReplacementRule(pattern350, lambda p, n, c, d, a, m, x, e, b : c*(m + S(2)*p + S(2))*Int(x**(m + S(1))*(a + b*ArcCot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) + m*Int(x**(m + S(-1))*(a + b*ArcCot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) - x**m*(a + b*ArcCot(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule350)

    pattern351 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda d, p: IntegerQ(p) | PositiveQ(d)))
    rule351 = ReplacementRule(pattern351, lambda p, n, c, d, a, m, x, e, b : c**(-m + S(-1))*d**p*Subst(Int((a + b*x)**n*Cos(x)**(-m - S(2)*p + S(-2))*Sin(x)**m, x), x, ArcTan(c*x)))
    rubi.add(rule351)

    pattern352 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda d, p: ~(IntegerQ(p) | PositiveQ(d))))
    rule352 = ReplacementRule(pattern352, lambda p, n, c, d, a, m, x, e, b : d**(p + S(1)/2)*Int(x**m*(a + b*ArcTan(c*x))**n*(c**S(2)*x**S(2) + S(1))**p, x)*Sqrt(c**S(2)*x**S(2) + S(1))/Sqrt(d + e*x**S(2)))
    rubi.add(rule352)

    pattern353 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda p: IntegerQ(p)))
    rule353 = ReplacementRule(pattern353, lambda p, n, c, d, a, m, x, e, b : -c**(-m + S(-1))*d**p*Subst(Int((a + b*x)**n*Cos(x)**m*Sin(x)**(-m - S(2)*p + S(-2)), x), x, ArcCot(c*x)))
    rubi.add(rule353)

    pattern354 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule354 = ReplacementRule(pattern354, lambda p, n, c, d, a, m, x, e, b : -c**(-m)*d**(p + S(1)/2)*x*Sqrt((c**S(2)*x**S(2) + S(1))/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*Cos(x)**m*Sin(x)**(-m - S(2)*p + S(-2)), x), x, ArcCot(c*x))/Sqrt(d + e*x**S(2)))
    rubi.add(rule354)

    pattern355 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule355 = ReplacementRule(pattern355, lambda p, c, d, a, x, e, b : -b*c*Int((d + e*x**S(2))**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x)/(S(2)*e*(p + S(1))) + (a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule355)

    pattern356 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule356 = ReplacementRule(pattern356, lambda p, c, d, a, x, e, b : b*c*Int((d + e*x**S(2))**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x)/(S(2)*e*(p + S(1))) + (a + b*ArcCot(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule356)

    pattern357 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & ~(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & ~(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & ~(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))))
    rule357 = ReplacementRule(pattern357, lambda p, c, d, a, m, x, e, b : With(List(Set(u, IntHide(x**m*(d + e*x**S(2))**p, x))), -b*c*Int(SimplifyIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcTan(c*x), u, x)))
    rubi.add(rule357)

    pattern358 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & ~(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & ~(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & ~(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))))
    rule358 = ReplacementRule(pattern358, lambda p, c, d, a, m, x, e, b : With(List(Set(u, IntHide(x**m*(d + e*x**S(2))**p, x))), b*c*Int(SimplifyIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcCot(c*x), u, x)))
    rubi.add(rule358)

    pattern359 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (IntegerQ(m) & Less(p, S(-1)) & Unequal(m, S(1)))))
    rule359 = ReplacementRule(pattern359, lambda p, n, c, d, a, m, x, e, b : Int(ExpandIntegrand((a + b*ArcTan(c*x))**n, x**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule359)

    pattern360 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (IntegerQ(m) & Less(p, S(-1)) & Unequal(m, S(1)))))
    rule360 = ReplacementRule(pattern360, lambda p, n, c, d, a, m, x, e, b : Int(ExpandIntegrand((a + b*ArcCot(c*x))**n, x**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule360)

    pattern361 = Pattern(Integral(x_**WC('m', S(1))*(a_ + ArcTan(x_*WC('c', S(1)))*WC('b', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule361 = ReplacementRule(pattern361, lambda p, c, d, a, m, x, e, b : a*Int(x**m*(d + e*x**S(2))**p, x) + b*Int(x**m*(d + e*x**S(2))**p*ArcTan(c*x), x))
    rubi.add(rule361)

    pattern362 = Pattern(Integral(x_**WC('m', S(1))*(a_ + ArcCot(x_*WC('c', S(1)))*WC('b', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule362 = ReplacementRule(pattern362, lambda p, c, d, a, m, x, e, b : a*Int(x**m*(d + e*x**S(2))**p, x) + b*Int(x**m*(d + e*x**S(2))**p*ArcCot(c*x), x))
    rubi.add(rule362)

    pattern363 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule363 = ReplacementRule(pattern363, lambda p, n, c, d, a, m, x, e, b : Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule363)

    pattern364 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule364 = ReplacementRule(pattern364, lambda p, n, c, d, a, m, x, e, b : Int(x**m*(a + b*ArcCot(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule364)

    pattern365 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*ArcTanh(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule365 = ReplacementRule(pattern365, lambda n, c, u, a, d, x, e, b : -Int((a + b*ArcTan(c*x))**n*Log(-u + S(1))/(d + e*x**S(2)), x)/S(2) + Int((a + b*ArcTan(c*x))**n*Log(u + S(1))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule365)

    pattern366 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*ArcCoth(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule366 = ReplacementRule(pattern366, lambda n, c, u, a, d, x, e, b : -Int((a + b*ArcCot(c*x))**n*Log(SimplifyIntegrand(S(1) - S(1)/u, x))/(d + e*x**S(2)), x)/S(2) + Int((a + b*ArcCot(c*x))**n*Log(SimplifyIntegrand(S(1) + 1/u, x))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule366)

    pattern367 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*ArcTanh(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule367 = ReplacementRule(pattern367, lambda n, c, u, a, d, x, e, b : -Int((a + b*ArcTan(c*x))**n*Log(-u + S(1))/(d + e*x**S(2)), x)/S(2) + Int((a + b*ArcTan(c*x))**n*Log(u + S(1))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule367)

    pattern368 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*ArcCoth(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule368 = ReplacementRule(pattern368, lambda n, c, u, a, d, x, e, b : -Int((a + b*ArcCot(c*x))**n*Log(SimplifyIntegrand(S(1) - S(1)/u, x))/(d + e*x**S(2)), x)/S(2) + Int((a + b*ArcCot(c*x))**n*Log(SimplifyIntegrand(S(1) + 1/u, x))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule368)

    pattern369 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ((-u + S(1))**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule369 = ReplacementRule(pattern369, lambda n, c, u, a, d, x, e, b : -ImaginaryI*b*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(S(2), Together(-u + S(1)))/(d + e*x**S(2)), x)/S(2) + ImaginaryI*(a + b*ArcTan(c*x))**n*PolyLog(S(2), Together(-u + S(1)))/(S(2)*c*d))
    rubi.add(rule369)

    pattern370 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ((-u + S(1))**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule370 = ReplacementRule(pattern370, lambda n, c, u, a, d, x, e, b : ImaginaryI*b*n*Int((a + b*ArcCot(c*x))**(n + S(-1))*PolyLog(S(2), Together(-u + S(1)))/(d + e*x**S(2)), x)/S(2) + ImaginaryI*(a + b*ArcCot(c*x))**n*PolyLog(S(2), Together(-u + S(1)))/(S(2)*c*d))
    rubi.add(rule370)

    pattern371 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ((-u + S(1))**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule371 = ReplacementRule(pattern371, lambda n, c, u, a, d, x, e, b : ImaginaryI*b*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(S(2), Together(-u + S(1)))/(d + e*x**S(2)), x)/S(2) - ImaginaryI*(a + b*ArcTan(c*x))**n*PolyLog(S(2), Together(-u + S(1)))/(S(2)*c*d))
    rubi.add(rule371)

    pattern372 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*Log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ((-u + S(1))**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule372 = ReplacementRule(pattern372, lambda n, c, u, a, d, x, e, b : -ImaginaryI*b*n*Int((a + b*ArcCot(c*x))**(n + S(-1))*PolyLog(S(2), Together(-u + S(1)))/(d + e*x**S(2)), x)/S(2) - ImaginaryI*(a + b*ArcCot(c*x))**n*PolyLog(S(2), Together(-u + S(1)))/(S(2)*c*d))
    rubi.add(rule372)

    pattern373 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule373 = ReplacementRule(pattern373, lambda p, n, c, d, a, u, x, e, b : ImaginaryI*b*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) - ImaginaryI*(a + b*ArcTan(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule373)

    pattern374 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule374 = ReplacementRule(pattern374, lambda p, n, c, d, a, u, x, e, b : -ImaginaryI*b*n*Int((a + b*ArcCot(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) - ImaginaryI*(a + b*ArcCot(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule374)

    pattern375 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule375 = ReplacementRule(pattern375, lambda p, n, c, d, a, u, x, e, b : -ImaginaryI*b*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) + ImaginaryI*(a + b*ArcTan(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule375)

    pattern376 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule376 = ReplacementRule(pattern376, lambda p, n, c, d, a, u, x, e, b : ImaginaryI*b*n*Int((a + b*ArcCot(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) + ImaginaryI*(a + b*ArcCot(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule376)

    pattern377 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)))
    rule377 = ReplacementRule(pattern377, lambda c, d, a, x, e, b : (-Log(a + b*ArcCot(c*x)) + Log(a + b*ArcTan(c*x)))/(b*c*d*(S(2)*a + b*ArcCot(c*x) + b*ArcTan(c*x))))
    rubi.add(rule377)

    pattern378 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Inequality(S(0), Less, n, LessEqual, m)))
    rule378 = ReplacementRule(pattern378, lambda n, c, d, a, m, x, e, b : n*Int((a + b*ArcCot(c*x))**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(-1))/(d + e*x**S(2)), x)/(m + S(1)) - (a + b*ArcCot(c*x))**(m + S(1))*(a + b*ArcTan(c*x))**n/(b*c*d*(m + S(1))))
    rubi.add(rule378)

    pattern379 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Less(S(0), n, m)))
    rule379 = ReplacementRule(pattern379, lambda n, c, d, a, m, x, e, b : n*Int((a + b*ArcCot(c*x))**(n + S(-1))*(a + b*ArcTan(c*x))**(m + S(1))/(d + e*x**S(2)), x)/(m + S(1)) + (a + b*ArcCot(c*x))**n*(a + b*ArcTan(c*x))**(m + S(1))/(b*c*d*(m + S(1))))
    rubi.add(rule379)

    pattern380 = Pattern(Integral(ArcTan(x_*WC('a', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda d, a, n, c: ~(Equal(n, S(2)) & ZeroQ(-a**S(2)*c + d))))
    rule380 = ReplacementRule(pattern380, lambda n, c, d, a, x : ImaginaryI*Int(Log(-ImaginaryI*a*x + S(1))/(c + d*x**n), x)/S(2) - ImaginaryI*Int(Log(ImaginaryI*a*x + S(1))/(c + d*x**n), x)/S(2))
    rubi.add(rule380)

    pattern381 = Pattern(Integral(ArcCot(x_*WC('a', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda d, a, n, c: ~(Equal(n, S(2)) & ZeroQ(-a**S(2)*c + d))))
    rule381 = ReplacementRule(pattern381, lambda n, c, d, a, x : ImaginaryI*Int(Log(-ImaginaryI/(a*x) + S(1))/(c + d*x**n), x)/S(2) - ImaginaryI*Int(Log(ImaginaryI/(a*x) + S(1))/(c + d*x**n), x)/S(2))
    rubi.add(rule381)

    pattern382 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(Log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)))
    rule382 = ReplacementRule(pattern382, lambda c, d, a, f, x, e, b, g : -b*c*Int(x*(d + e*Log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x) - S(2)*e*g*Int(x**S(2)*(a + b*ArcTan(c*x))/(f + g*x**S(2)), x) + x*(a + b*ArcTan(c*x))*(d + e*Log(f + g*x**S(2))))
    rubi.add(rule382)

    pattern383 = Pattern(Integral((ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(Log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)))
    rule383 = ReplacementRule(pattern383, lambda c, d, a, f, x, e, b, g : b*c*Int(x*(d + e*Log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x) - S(2)*e*g*Int(x**S(2)*(a + b*ArcCot(c*x))/(f + g*x**S(2)), x) + x*(a + b*ArcCot(c*x))*(d + e*Log(f + g*x**S(2))))
    rubi.add(rule383)

    pattern384 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(Log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2))))
    rule384 = ReplacementRule(pattern384, lambda c, d, a, f, m, x, e, b, g : -b*c*Int(x**(m + S(1))*(d + e*Log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) - S(2)*e*g*Int(x**(m + S(2))*(a + b*ArcTan(c*x))/(f + g*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcTan(c*x))*(d + e*Log(f + g*x**S(2)))/(m + S(1)))
    rubi.add(rule384)

    pattern385 = Pattern(Integral(x_**WC('m', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(Log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2))))
    rule385 = ReplacementRule(pattern385, lambda c, d, a, f, m, x, e, b, g : b*c*Int(x**(m + S(1))*(d + e*Log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) - S(2)*e*g*Int(x**(m + S(2))*(a + b*ArcCot(c*x))/(f + g*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcCot(c*x))*(d + e*Log(f + g*x**S(2)))/(m + S(1)))
    rubi.add(rule385)

    pattern386 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(Log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m/S(2) + S(1)/2)))
    rule386 = ReplacementRule(pattern386, lambda c, d, a, f, m, x, e, b, g : With(List(Set(u, IntHide(x**m*(d + e*Log(f + g*x**S(2))), x))), -b*c*Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcTan(c*x), u, x)))
    rubi.add(rule386)

    pattern387 = Pattern(Integral(x_**WC('m', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(Log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m/S(2) + S(1)/2)))
    rule387 = ReplacementRule(pattern387, lambda c, d, a, f, m, x, e, b, g : With(List(Set(u, IntHide(x**m*(d + e*Log(f + g*x**S(2))), x))), b*c*Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcCot(c*x), u, x)))
    rubi.add(rule387)

    pattern388 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(Log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Unequal(m, S(-1))))
    rule388 = ReplacementRule(pattern388, lambda c, d, a, f, m, x, e, b, g : With(List(Set(u, IntHide(x**m*(a + b*ArcTan(c*x)), x))), -S(2)*e*g*Int(ExpandIntegrand(u*x/(f + g*x**S(2)), x), x) + Dist(d + e*Log(f + g*x**S(2)), u, x)))
    rubi.add(rule388)

    pattern389 = Pattern(Integral(x_**WC('m', S(1))*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(Log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Unequal(m, S(-1))))
    rule389 = ReplacementRule(pattern389, lambda c, d, a, f, m, x, e, b, g : With(List(Set(u, IntHide(x**m*(a + b*ArcCot(c*x)), x))), -S(2)*e*g*Int(ExpandIntegrand(u*x/(f + g*x**S(2)), x), x) + Dist(d + e*Log(f + g*x**S(2)), u, x)))
    rubi.add(rule389)

    pattern390 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**S(2)*(Log(f_ + x_**S(2)*WC('g', S(1)))*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda f, g, c: ZeroQ(-c**S(2)*f + g)))
    rule390 = ReplacementRule(pattern390, lambda c, d, a, f, x, e, b, g : b*c*e*Int(x**S(2)*(a + b*ArcTan(c*x))/(c**S(2)*x**S(2) + S(1)), x) - b*Int((a + b*ArcTan(c*x))*(d + e*Log(f + g*x**S(2))), x)/c - e*x**S(2)*(a + b*ArcTan(c*x))**S(2)/S(2) + (a + b*ArcTan(c*x))**S(2)*(d + e*Log(f + g*x**S(2)))*(f + g*x**S(2))/(S(2)*g))
    rubi.add(rule390)

    pattern391 = Pattern(Integral(x_*(ArcCot(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**S(2)*(Log(f_ + x_**S(2)*WC('g', S(1)))*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda f, g, c: ZeroQ(-c**S(2)*f + g)))
    rule391 = ReplacementRule(pattern391, lambda c, d, a, f, x, e, b, g : -b*c*e*Int(x**S(2)*(a + b*ArcCot(c*x))/(c**S(2)*x**S(2) + S(1)), x) + b*Int((a + b*ArcCot(c*x))*(d + e*Log(f + g*x**S(2))), x)/c - e*x**S(2)*(a + b*ArcCot(c*x))**S(2)/S(2) + (a + b*ArcCot(c*x))**S(2)*(d + e*Log(f + g*x**S(2)))*(f + g*x**S(2))/(S(2)*g))
    rubi.add(rule391)

    pattern392 = Pattern(Integral(exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: OddQ(ImaginaryI*n)))
    rule392 = ReplacementRule(pattern392, lambda a, n, x : Int((-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + S(1)/2)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + S(1)/2)/Sqrt(a**S(2)*x**S(2) + S(1)), x))
    rubi.add(rule392)

    pattern393 = Pattern(Integral(x_**WC('m', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: OddQ(ImaginaryI*n)))
    rule393 = ReplacementRule(pattern393, lambda a, n, m, x : Int(x**m*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + S(1)/2)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + S(1)/2)/Sqrt(a**S(2)*x**S(2) + S(1)), x))
    rubi.add(rule393)

    pattern394 = Pattern(Integral(exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(OddQ(ImaginaryI*n))))
    rule394 = ReplacementRule(pattern394, lambda a, n, x : Int((-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule394)

    pattern395 = Pattern(Integral(x_**WC('m', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(OddQ(ImaginaryI*n))))
    rule395 = ReplacementRule(pattern395, lambda a, n, m, x : Int(x**m*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule395)

    pattern396 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*c**S(2) + d**S(2))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)))
    rule396 = ReplacementRule(pattern396, lambda p, n, c, u, a, d, x : c**p*Int(u*(S(1) + d*x/c)**p*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule396)

    pattern397 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*c**S(2) + d**S(2))), CustomConstraint(lambda p, c: ~(IntegerQ(p) | PositiveQ(c))))
    rule397 = ReplacementRule(pattern397, lambda p, n, c, u, a, d, x : Int(u*(c + d*x)**p*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule397)

    pattern398 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule398 = ReplacementRule(pattern398, lambda p, n, c, u, a, d, x : d**p*Int(u*x**(-p)*(c*x/d + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule398)

    pattern399 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(n_*ArcTanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))), CustomConstraint(lambda c: PositiveQ(c)))
    rule399 = ReplacementRule(pattern399, lambda p, n, c, u, a, d, x : (S(-1))**(n/S(2))*c**p*Int(u*(S(1) - S(1)/(ImaginaryI*a*x))**(ImaginaryI*n/S(2))*(S(1) + S(1)/(ImaginaryI*a*x))**(-ImaginaryI*n/S(2))*(S(1) + d/(c*x))**p, x))
    rubi.add(rule399)

    pattern400 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))), CustomConstraint(lambda c: ~(PositiveQ(c))))
    rule400 = ReplacementRule(pattern400, lambda p, n, c, u, a, d, x : Int(u*(c + d/x)**p*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule400)

    pattern401 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule401 = ReplacementRule(pattern401, lambda p, n, c, u, a, d, x : x**p*(c + d/x)**p*(c*x/d + S(1))**(-p)*Int(u*x**(-p)*(c*x/d + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule401)

    pattern402 = Pattern(Integral(exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n))))
    rule402 = ReplacementRule(pattern402, lambda n, c, d, a, x : (a*x + n)*exp(n*ArcTan(a*x))/(a*c*(n**S(2) + S(1))*Sqrt(c + d*x**S(2))))
    rubi.add(rule402)

    pattern403 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n))), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) + S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule403 = ReplacementRule(pattern403, lambda p, n, c, d, a, x : S(2)*(p + S(1))*(S(2)*p + S(3))*Int((c + d*x**S(2))**(p + S(1))*exp(n*ArcTan(a*x)), x)/(c*(n**S(2) + S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(-S(2)*a*x*(p + S(1)) + n)*exp(n*ArcTan(a*x))/(a*c*(n**S(2) + S(4)*(p + S(1))**S(2))))
    rubi.add(rule403)

    pattern404 = Pattern(Integral(exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)))
    rule404 = ReplacementRule(pattern404, lambda n, c, d, a, x : exp(n*ArcTan(a*x))/(a*c*n))
    rubi.add(rule404)

    pattern405 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2) + S(1)/2)), CustomConstraint(lambda p, n: ~(IntegerQ(-ImaginaryI*n/S(2) + p))))
    rule405 = ReplacementRule(pattern405, lambda p, n, c, d, a, x : c**p*Int((a**S(2)*x**S(2) + S(1))**(-ImaginaryI*n/S(2) + p)*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n), x))
    rubi.add(rule405)

    pattern406 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)))
    rule406 = ReplacementRule(pattern406, lambda p, n, c, d, a, x : c**p*Int((-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + p), x))
    rubi.add(rule406)

    pattern407 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: ~(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: PositiveIntegerQ(ImaginaryI*n/S(2))))
    rule407 = ReplacementRule(pattern407, lambda p, n, c, d, a, x : c**(ImaginaryI*n/S(2))*Int((c + d*x**S(2))**(-ImaginaryI*n/S(2) + p)*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n), x))
    rubi.add(rule407)

    pattern408 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: ~(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: NegativeIntegerQ(ImaginaryI*n/S(2))))
    rule408 = ReplacementRule(pattern408, lambda p, n, c, d, a, x : c**(-ImaginaryI*n/S(2))*Int((c + d*x**S(2))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n), x))
    rubi.add(rule408)

    pattern409 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: ~(IntegerQ(p) | PositiveQ(c))))
    rule409 = ReplacementRule(pattern409, lambda p, n, c, d, a, x : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(a**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule409)

    pattern410 = Pattern(Integral(x_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n))))
    rule410 = ReplacementRule(pattern410, lambda n, c, d, a, x : (a*n*x + S(-1))*exp(n*ArcTan(a*x))/(d*(n**S(2) + S(1))*Sqrt(c + d*x**S(2))))
    rubi.add(rule410)

    pattern411 = Pattern(Integral(x_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule411 = ReplacementRule(pattern411, lambda p, n, c, d, a, x : -a*c*n*Int((c + d*x**S(2))**p*exp(n*ArcTan(a*x)), x)/(S(2)*d*(p + S(1))) + (c + d*x**S(2))**(p + S(1))*exp(n*ArcTan(a*x))/(S(2)*d*(p + S(1))))
    rubi.add(rule411)

    pattern412 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, n: ZeroQ(n**S(2) - S(2)*p + S(-2))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n))))
    rule412 = ReplacementRule(pattern412, lambda p, n, c, d, a, x : (c + d*x**S(2))**(p + S(1))*(a*n*x + S(-1))*exp(n*ArcTan(a*x))/(a*d*n*(n**S(2) + S(1))))
    rubi.add(rule412)

    pattern413 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n))), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) + S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule413 = ReplacementRule(pattern413, lambda p, n, c, d, a, x : (n**S(2) - S(2)*p + S(-2))*Int((c + d*x**S(2))**(p + S(1))*exp(n*ArcTan(a*x)), x)/(d*(n**S(2) + S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) - n)*exp(n*ArcTan(a*x))/(a*d*(n**S(2) + S(4)*(p + S(1))**S(2))))
    rubi.add(rule413)

    pattern414 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2) + S(1)/2)), CustomConstraint(lambda p, n: ~(IntegerQ(-ImaginaryI*n/S(2) + p))))
    rule414 = ReplacementRule(pattern414, lambda p, n, c, d, a, m, x : c**p*Int(x**m*(a**S(2)*x**S(2) + S(1))**(-ImaginaryI*n/S(2) + p)*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n), x))
    rubi.add(rule414)

    pattern415 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)))
    rule415 = ReplacementRule(pattern415, lambda p, n, c, d, a, m, x : c**p*Int(x**m*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + p), x))
    rubi.add(rule415)

    pattern416 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: ~(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: PositiveIntegerQ(ImaginaryI*n/S(2))))
    rule416 = ReplacementRule(pattern416, lambda p, n, c, d, a, m, x : c**(ImaginaryI*n/S(2))*Int(x**m*(c + d*x**S(2))**(-ImaginaryI*n/S(2) + p)*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n), x))
    rubi.add(rule416)

    pattern417 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: ~(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: NegativeIntegerQ(ImaginaryI*n/S(2))))
    rule417 = ReplacementRule(pattern417, lambda p, n, c, d, a, m, x : c**(-ImaginaryI*n/S(2))*Int(x**m*(c + d*x**S(2))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n), x))
    rubi.add(rule417)

    pattern418 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: ~(IntegerQ(p) | PositiveQ(c))))
    rule418 = ReplacementRule(pattern418, lambda p, n, c, d, a, m, x : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(a**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x**m*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule418)

    pattern419 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)))
    rule419 = ReplacementRule(pattern419, lambda p, n, c, d, a, u, x : c**p*Int(u*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + p), x))
    rubi.add(rule419)

    pattern420 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))))
    rule420 = ReplacementRule(pattern420, lambda p, n, c, d, a, u, x : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(-ImaginaryI*a*x + S(1))**(-FracPart(p))*(ImaginaryI*a*x + S(1))**(-FracPart(p))*Int(u*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + p), x))
    rubi.add(rule420)

    pattern421 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: ~(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))))
    rule421 = ReplacementRule(pattern421, lambda p, n, c, d, a, u, x : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(a**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(u*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule421)

    pattern422 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda p: IntegerQ(p)))
    rule422 = ReplacementRule(pattern422, lambda p, n, c, u, a, d, x : d**p*Int(u*x**(-S(2)*p)*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule422)

    pattern423 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))), CustomConstraint(lambda c: PositiveQ(c)))
    rule423 = ReplacementRule(pattern423, lambda p, n, c, u, a, d, x : c**p*Int(u*(-ImaginaryI/(a*x) + S(1))**p*(ImaginaryI/(a*x) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule423)

    pattern424 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))), CustomConstraint(lambda c: ~(PositiveQ(c))))
    rule424 = ReplacementRule(pattern424, lambda p, n, c, u, a, d, x : x**(S(2)*p)*(c + d/x**S(2))**p*(-ImaginaryI*a*x + S(1))**(-p)*(ImaginaryI*a*x + S(1))**(-p)*Int(u*x**(-S(2)*p)*(-ImaginaryI*a*x + S(1))**p*(ImaginaryI*a*x + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule424)

    pattern425 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda p: ~(IntegerQ(p))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))))
    rule425 = ReplacementRule(pattern425, lambda p, n, c, u, a, d, x : x**(S(2)*p)*(c + d/x**S(2))**p*(a**S(2)*x**S(2) + S(1))**(-p)*Int(u*x**(-S(2)*p)*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule425)

    pattern426 = Pattern(Integral(exp(ArcTan((a_ + x_*WC('b', S(1)))*WC('c', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule426 = ReplacementRule(pattern426, lambda n, c, a, x, b : Int((-ImaginaryI*a*c - ImaginaryI*b*c*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule426)

    pattern427 = Pattern(Integral(x_**m_*exp(n_*ArcTan((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda n: RationalQ(ImaginaryI*n)), CustomConstraint(lambda n: Less(S(-1), ImaginaryI*n, S(1))))
    rule427 = ReplacementRule(pattern427, lambda n, c, a, x, m, b : S(4)*ImaginaryI**(-m)*b**(-m + S(-1))*c**(-m + S(-1))*Subst(Int(x**(S(2)/(ImaginaryI*n))*(x**(S(2)/(ImaginaryI*n)) + S(1))**(-m + S(-2))*(-ImaginaryI*a*c - x**(S(2)/(ImaginaryI*n))*(ImaginaryI*a*c + S(1)) + S(1))**m, x), x, (-ImaginaryI*c*(a + b*x) + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*c*(a + b*x) + S(1))**(-ImaginaryI*n/S(2)))/n)
    rubi.add(rule427)

    pattern428 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*exp(ArcTan((a_ + x_*WC('b', S(1)))*WC('c', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule428 = ReplacementRule(pattern428, lambda n, c, d, a, m, x, e, b : Int((d + e*x)**m*(-ImaginaryI*a*c - ImaginaryI*b*c*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule428)

    pattern429 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(a_ + x_*WC('b', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, b, a, e: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda a, b, c, e: ZeroQ(b**S(2)*c - e*(a**S(2) + S(1)))), CustomConstraint(lambda a, p, c: IntegerQ(p) | PositiveQ(c/(a**S(2) + S(1)))))
    rule429 = ReplacementRule(pattern429, lambda p, n, c, u, d, a, x, e, b : (c/(a**S(2) + S(1)))**p*Int(u*(-ImaginaryI*a - ImaginaryI*b*x + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a + ImaginaryI*b*x + S(1))**(-ImaginaryI*n/S(2) + p), x))
    rubi.add(rule429)

    pattern430 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(a_ + x_*WC('b', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, b, a, e: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda a, b, c, e: ZeroQ(b**S(2)*c - e*(a**S(2) + S(1)))), CustomConstraint(lambda a, p, c: ~(IntegerQ(p) | PositiveQ(c/(a**S(2) + S(1))))))
    rule430 = ReplacementRule(pattern430, lambda p, n, c, u, d, a, x, e, b : (c + d*x + e*x**S(2))**p*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**(-p)*Int(u*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule430)

    pattern431 = Pattern(Integral(WC('u', S(1))*exp(ArcTan(WC('c', S(1))/(x_*WC('b', S(1)) + WC('a', S(0))))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule431 = ReplacementRule(pattern431, lambda n, c, u, a, x, b : Int(u*exp(n*ArcCot(a/c + b*x/c)), x))
    rubi.add(rule431)

    pattern432 = Pattern(Integral(WC('u', S(1))*exp(n_*ArcCot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))))
    rule432 = ReplacementRule(pattern432, lambda u, a, n, x : (S(-1))**(ImaginaryI*n/S(2))*Int(u*exp(-n*ArcTan(a*x)), x))
    rubi.add(rule432)

    pattern433 = Pattern(Integral(exp(n_*ArcCot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: OddQ(ImaginaryI*n)))
    rule433 = ReplacementRule(pattern433, lambda a, n, x : -Subst(Int((-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + S(1)/2)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + S(1)/2)/(x**S(2)*Sqrt(S(1) + x**S(2)/a**S(2))), x), x, 1/x))
    rubi.add(rule433)

    pattern434 = Pattern(Integral(x_**WC('m', S(1))*exp(n_*ArcCot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: OddQ(ImaginaryI*n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule434 = ReplacementRule(pattern434, lambda a, n, m, x : -Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + S(1)/2)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + S(1)/2)/Sqrt(S(1) + x**S(2)/a**S(2)), x), x, 1/x))
    rubi.add(rule434)

    pattern435 = Pattern(Integral(exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n))))
    rule435 = ReplacementRule(pattern435, lambda a, n, x : -Subst(Int((-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2))/x**S(2), x), x, 1/x))
    rubi.add(rule435)

    pattern436 = Pattern(Integral(x_**WC('m', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n))), CustomConstraint(lambda m: IntegerQ(m)))
    rule436 = ReplacementRule(pattern436, lambda a, n, m, x : -Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(n/S(2))*(ImaginaryI*x/a + S(1))**(-n/S(2)), x), x, 1/x))
    rubi.add(rule436)

    pattern437 = Pattern(Integral(x_**m_*exp(n_*ArcCot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: OddQ(ImaginaryI*n)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule437 = ReplacementRule(pattern437, lambda a, m, n, x : -x**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + S(1)/2)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + S(1)/2)/Sqrt(S(1) + x**S(2)/a**S(2)), x), x, 1/x))
    rubi.add(rule437)

    pattern438 = Pattern(Integral(x_**m_*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule438 = ReplacementRule(pattern438, lambda a, m, n, x : -Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(n/S(2))*(ImaginaryI*x/a + S(1))**(-n/S(2)), x), x, 1/x))
    rubi.add(rule438)

    pattern439 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*c**S(2) + d**S(2))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p: IntegerQ(p)))
    rule439 = ReplacementRule(pattern439, lambda p, n, c, u, a, d, x : d**p*Int(u*x**p*(c/(d*x) + S(1))**p*exp(n*ArcCot(a*x)), x))
    rubi.add(rule439)

    pattern440 = Pattern(Integral((c_ + x_*WC('d', S(1)))**p_*WC('u', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*c**S(2) + d**S(2))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule440 = ReplacementRule(pattern440, lambda p, n, c, u, a, d, x : x**(-p)*(c + d*x)**p*(c/(d*x) + S(1))**(-p)*Int(u*x**p*(c/(d*x) + S(1))**p*exp(n*ArcCot(a*x)), x))
    rubi.add(rule440)

    pattern441 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)))
    rule441 = ReplacementRule(pattern441, lambda p, n, c, d, a, x : -c**p*Subst(Int((S(1) + d*x/c)**p*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2))/x**S(2), x), x, 1/x))
    rubi.add(rule441)

    pattern442 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda m: IntegerQ(m)))
    rule442 = ReplacementRule(pattern442, lambda p, n, c, d, a, m, x : -c**p*Subst(Int(x**(-m + S(-2))*(S(1) + d*x/c)**p*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2)), x), x, 1/x))
    rubi.add(rule442)

    pattern443 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: ~(IntegerQ(p) | PositiveQ(c))))
    rule443 = ReplacementRule(pattern443, lambda p, n, c, d, a, x : (S(1) + d/(c*x))**(-p)*(c + d/x)**p*Int((S(1) + d/(c*x))**p*exp(n*ArcCot(a*x)), x))
    rubi.add(rule443)

    pattern444 = Pattern(Integral(x_**m_*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule444 = ReplacementRule(pattern444, lambda p, n, c, d, a, x, m : -c**p*x**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(S(1) + d*x/c)**p*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2)), x), x, 1/x))
    rubi.add(rule444)

    pattern445 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: ~(IntegerQ(p) | PositiveQ(c))))
    rule445 = ReplacementRule(pattern445, lambda p, n, c, u, a, d, x : (S(1) + d/(c*x))**(-p)*(c + d/x)**p*Int(u*(S(1) + d/(c*x))**p*exp(n*ArcCot(a*x)), x))
    rubi.add(rule445)

    pattern446 = Pattern(Integral(exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)))
    rule446 = ReplacementRule(pattern446, lambda n, c, d, a, x : -exp(n*ArcCot(a*x))/(a*c*n))
    rubi.add(rule446)

    pattern447 = Pattern(Integral(exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: ~(OddQ(ImaginaryI*n))))
    rule447 = ReplacementRule(pattern447, lambda n, c, d, a, x : (a*x - n)*exp(n*ArcCot(a*x))/(a*c*(n**S(2) + S(1))*Sqrt(c + d*x**S(2))))
    rubi.add(rule447)

    pattern448 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) + S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p, n: ~(IntegerQ(p) & EvenQ(ImaginaryI*n))), CustomConstraint(lambda p, n: ~(~(IntegerQ(p)) & OddQ(ImaginaryI*n))))
    rule448 = ReplacementRule(pattern448, lambda p, n, c, d, a, x : S(2)*(p + S(1))*(S(2)*p + S(3))*Int((c + d*x**S(2))**(p + S(1))*exp(n*ArcCot(a*x)), x)/(c*(n**S(2) + S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(-S(2)*a*x*(p + S(1)) - n)*exp(n*ArcCot(a*x))/(a*c*(n**S(2) + S(4)*(p + S(1))**S(2))))
    rubi.add(rule448)

    pattern449 = Pattern(Integral(x_*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: ~(OddQ(ImaginaryI*n))))
    rule449 = ReplacementRule(pattern449, lambda n, c, d, a, x : (-a*n*x + S(-1))*exp(n*ArcCot(a*x))/(a**S(2)*c*(n**S(2) + S(1))*Sqrt(c + d*x**S(2))))
    rubi.add(rule449)

    pattern450 = Pattern(Integral(x_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: LessEqual(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) + S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p, n: ~(IntegerQ(p) & EvenQ(ImaginaryI*n))), CustomConstraint(lambda p, n: ~(~(IntegerQ(p)) & OddQ(ImaginaryI*n))))
    rule450 = ReplacementRule(pattern450, lambda p, n, c, d, a, x : n*(S(2)*p + S(3))*Int((c + d*x**S(2))**(p + S(1))*exp(n*ArcCot(a*x)), x)/(a*c*(n**S(2) + S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(-a*n*x + S(2)*p + S(2))*exp(n*ArcCot(a*x))/(a**S(2)*c*(n**S(2) + S(4)*(p + S(1))**S(2))))
    rubi.add(rule450)

    pattern451 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, n: ZeroQ(n**S(2) - S(2)*p + S(-2))), CustomConstraint(lambda n: NonzeroQ(n**S(2) + S(1))))
    rule451 = ReplacementRule(pattern451, lambda p, n, c, d, a, x : (c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*ArcCot(a*x))/(a**S(3)*c*n**S(2)*(n**S(2) + S(1))))
    rubi.add(rule451)

    pattern452 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: LessEqual(p, S(-1))), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) - S(2)*p + S(-2))), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) + S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p, n: ~(IntegerQ(p) & EvenQ(ImaginaryI*n))), CustomConstraint(lambda p, n: ~(~(IntegerQ(p)) & OddQ(ImaginaryI*n))))
    rule452 = ReplacementRule(pattern452, lambda p, n, c, d, a, x : (n**S(2) - S(2)*p + S(-2))*Int((c + d*x**S(2))**(p + S(1))*exp(n*ArcCot(a*x)), x)/(a**S(2)*c*(n**S(2) + S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*ArcCot(a*x))/(a**S(3)*c*(n**S(2) + S(4)*(p + S(1))**S(2))))
    rubi.add(rule452)

    pattern453 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p, m: LessEqual(S(3), m, -S(2)*p + S(-2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule453 = ReplacementRule(pattern453, lambda p, n, c, d, a, m, x : -a**(-m + S(-1))*c**p*Subst(Int(Cos(x)**(-S(2)*p + S(-2))*Cot(x)**(m + S(2)*p + S(2))*exp(n*x), x), x, ArcCot(a*x)))
    rubi.add(rule453)

    pattern454 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p: IntegerQ(p)))
    rule454 = ReplacementRule(pattern454, lambda p, n, c, u, a, d, x : d**p*Int(u*x**(S(2)*p)*(S(1) + S(1)/(a**S(2)*x**S(2)))**p*exp(n*ArcCot(a*x)), x))
    rubi.add(rule454)

    pattern455 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*WC('u', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p: ~(IntegerQ(p))))
    rule455 = ReplacementRule(pattern455, lambda p, n, c, u, a, d, x : x**(-S(2)*p)*(S(1) + S(1)/(a**S(2)*x**S(2)))**(-p)*(c + d*x**S(2))**p*Int(u*x**(S(2)*p)*(S(1) + S(1)/(a**S(2)*x**S(2)))**p*exp(n*ArcCot(a*x)), x))
    rubi.add(rule455)

    pattern456 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*WC('u', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda p, n: IntegersQ(S(2)*p, ImaginaryI*n/S(2) + p)))
    rule456 = ReplacementRule(pattern456, lambda p, n, c, u, a, d, x : c**p*(ImaginaryI*a)**(-S(2)*p)*Int(u*x**(-S(2)*p)*(ImaginaryI*a*x + S(-1))**(-ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + p), x))
    rubi.add(rule456)

    pattern457 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda p, n: ~(IntegersQ(S(2)*p, ImaginaryI*n/S(2) + p))))
    rule457 = ReplacementRule(pattern457, lambda p, n, c, d, a, x : -c**p*Subst(Int((-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + p)/x**S(2), x), x, 1/x))
    rubi.add(rule457)

    pattern458 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda p, n: ~(IntegersQ(S(2)*p, ImaginaryI*n/S(2) + p))), CustomConstraint(lambda m: IntegerQ(m)))
    rule458 = ReplacementRule(pattern458, lambda p, n, c, d, a, m, x : -c**p*Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + p), x), x, 1/x))
    rubi.add(rule458)

    pattern459 = Pattern(Integral(x_**m_*(c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda p, n: ~(IntegersQ(S(2)*p, ImaginaryI*n/S(2) + p))), CustomConstraint(lambda m: ~(IntegerQ(m))))
    rule459 = ReplacementRule(pattern459, lambda p, n, c, d, a, x, m : -c**p*x**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + p), x), x, 1/x))
    rubi.add(rule459)

    pattern460 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(ArcCot(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, a, c: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: ~(IntegerQ(p) | PositiveQ(c))))
    rule460 = ReplacementRule(pattern460, lambda p, n, c, u, a, d, x : (S(1) + S(1)/(a**S(2)*x**S(2)))**(-p)*(c + d/x**S(2))**p*Int(u*(S(1) + S(1)/(a**S(2)*x**S(2)))**p*exp(n*ArcCot(a*x)), x))
    rubi.add(rule460)

    pattern461 = Pattern(Integral(WC('u', S(1))*exp(n_*ArcCot((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))))
    rule461 = ReplacementRule(pattern461, lambda n, c, u, a, x, b : (S(-1))**(ImaginaryI*n/S(2))*Int(u*exp(-n*ArcTan(c*(a + b*x))), x))
    rubi.add(rule461)

    pattern462 = Pattern(Integral(exp(ArcCot((a_ + x_*WC('b', S(1)))*WC('c', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))))
    rule462 = ReplacementRule(pattern462, lambda n, c, a, x, b : (ImaginaryI*c*(a + b*x))**(ImaginaryI*n/S(2))*(S(1) + S(1)/(ImaginaryI*c*(a + b*x)))**(ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(-ImaginaryI*n/S(2))*Int((ImaginaryI*a*c + ImaginaryI*b*c*x + S(-1))**(-ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(ImaginaryI*n/S(2)), x))
    rubi.add(rule462)

    pattern463 = Pattern(Integral(x_**m_*exp(n_*ArcCoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda n: RationalQ(ImaginaryI*n)), CustomConstraint(lambda n: Less(S(-1), ImaginaryI*n, S(1))))
    rule463 = ReplacementRule(pattern463, lambda n, c, a, x, m, b : S(4)*ImaginaryI**(-m)*b**(-m + S(-1))*c**(-m + S(-1))*Subst(Int(x**(S(2)/(ImaginaryI*n))*(x**(S(2)/(ImaginaryI*n)) + S(-1))**(-m + S(-2))*(ImaginaryI*a*c + x**(S(2)/(ImaginaryI*n))*(-ImaginaryI*a*c + S(1)) + S(1))**m, x), x, (S(1) - S(1)/(ImaginaryI*c*(a + b*x)))**(-ImaginaryI*n/S(2))*(S(1) + S(1)/(ImaginaryI*c*(a + b*x)))**(ImaginaryI*n/S(2)))/n)
    rubi.add(rule463)

    pattern464 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*exp(ArcCoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))))
    rule464 = ReplacementRule(pattern464, lambda n, c, d, a, m, x, e, b : (ImaginaryI*c*(a + b*x))**(ImaginaryI*n/S(2))*(S(1) + S(1)/(ImaginaryI*c*(a + b*x)))**(ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(-ImaginaryI*n/S(2))*Int((d + e*x)**m*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(-1))**(-ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(ImaginaryI*n/S(2)), x))
    rubi.add(rule464)

    pattern465 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcCot(a_ + x_*WC('b', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda d, b, a, e: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda a, b, c, e: ZeroQ(b**S(2)*c - e*(a**S(2) + S(1)))), CustomConstraint(lambda a, p, c: IntegerQ(p) | PositiveQ(c/(a**S(2) + S(1)))))
    rule465 = ReplacementRule(pattern465, lambda p, n, c, u, d, a, x, e, b : (c/(a**S(2) + S(1)))**p*((ImaginaryI*a + ImaginaryI*b*x + S(1))/(ImaginaryI*a + ImaginaryI*b*x))**(ImaginaryI*n/S(2))*((ImaginaryI*a + ImaginaryI*b*x)/(ImaginaryI*a + ImaginaryI*b*x + S(1)))**(ImaginaryI*n/S(2))*(-ImaginaryI*a - ImaginaryI*b*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a + ImaginaryI*b*x + S(-1))**(-ImaginaryI*n/S(2))*Int(u*(-ImaginaryI*a - ImaginaryI*b*x + S(1))**(-ImaginaryI*n/S(2) + p)*(ImaginaryI*a + ImaginaryI*b*x + S(1))**(ImaginaryI*n/S(2) + p), x))
    rubi.add(rule465)

    pattern466 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcCot(a_ + x_*WC('b', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: ~(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda d, b, a, e: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda a, b, c, e: ZeroQ(b**S(2)*c - e*(a**S(2) + S(1)))), CustomConstraint(lambda a, p, c: ~(IntegerQ(p) | PositiveQ(c/(a**S(2) + S(1))))))
    rule466 = ReplacementRule(pattern466, lambda p, n, c, u, d, a, x, e, b : (c + d*x + e*x**S(2))**p*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**(-p)*Int(u*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**p*exp(n*ArcCot(a*x)), x))
    rubi.add(rule466)

    pattern467 = Pattern(Integral(WC('u', S(1))*exp(ArcCot(WC('c', S(1))/(x_*WC('b', S(1)) + WC('a', S(0))))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule467 = ReplacementRule(pattern467, lambda n, c, u, a, x, b : Int(u*exp(n*ArcTan(a/c + b*x/c)), x))
    rubi.add(rule467)

    pattern468 = Pattern(Integral((ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule468 = ReplacementRule(pattern468, lambda n, c, d, a, x, b : Subst(Int((a + b*ArcTan(x))**n, x), x, c + d*x)/d)
    rubi.add(rule468)

    pattern469 = Pattern(Integral((ArcCot(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule469 = ReplacementRule(pattern469, lambda n, c, d, a, x, b : Subst(Int((a + b*ArcCot(x))**n, x), x, c + d*x)/d)
    rubi.add(rule469)

    pattern470 = Pattern(Integral((ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(PositiveIntegerQ(n))))
    rule470 = ReplacementRule(pattern470, lambda n, c, d, a, x, b : Int((a + b*ArcTan(c + d*x))**n, x))
    rubi.add(rule470)

    pattern471 = Pattern(Integral((ArcCot(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(PositiveIntegerQ(n))))
    rule471 = ReplacementRule(pattern471, lambda n, c, d, a, x, b : Int((a + b*ArcCot(c + d*x))**n, x))
    rubi.add(rule471)

    pattern472 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule472 = ReplacementRule(pattern472, lambda n, c, d, a, f, m, x, e, b : Subst(Int((a + b*ArcTan(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule472)

    pattern473 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcCot(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule473 = ReplacementRule(pattern473, lambda n, c, d, a, f, m, x, e, b : Subst(Int((a + b*ArcCot(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule473)

    pattern474 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**m_*(ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(PositiveIntegerQ(n))))
    rule474 = ReplacementRule(pattern474, lambda n, c, d, a, f, x, m, e, b : Int((a + b*ArcTan(c + d*x))**n*(e + f*x)**m, x))
    rubi.add(rule474)

    pattern475 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**m_*(ArcCot(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(PositiveIntegerQ(n))))
    rule475 = ReplacementRule(pattern475, lambda n, c, d, a, f, x, m, e, b : Int((a + b*ArcCot(c + d*x))**n*(e + f*x)**m, x))
    rubi.add(rule475)

    pattern476 = Pattern(Integral((ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, B, c, A: ZeroQ(-S(2)*A*c*d + B*(c**S(2) + S(1)))), CustomConstraint(lambda d, C, B, c: ZeroQ(-B*d + S(2)*C*c)))
    rule476 = ReplacementRule(pattern476, lambda p, n, B, c, d, a, C, x, A, b : Subst(Int((a + b*ArcTan(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule476)

    pattern477 = Pattern(Integral((ArcCot(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, B, c, A: ZeroQ(-S(2)*A*c*d + B*(c**S(2) + S(1)))), CustomConstraint(lambda d, C, B, c: ZeroQ(-B*d + S(2)*C*c)))
    rule477 = ReplacementRule(pattern477, lambda p, n, B, c, d, a, C, x, A, b : Subst(Int((a + b*ArcCot(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule477)

    pattern478 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, B, c, A: ZeroQ(-S(2)*A*c*d + B*(c**S(2) + S(1)))), CustomConstraint(lambda d, C, B, c: ZeroQ(-B*d + S(2)*C*c)))
    rule478 = ReplacementRule(pattern478, lambda p, n, B, c, d, a, f, C, m, x, e, A, b : Subst(Int((a + b*ArcTan(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule478)

    pattern479 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcCot(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda d, B, c, A: ZeroQ(-S(2)*A*c*d + B*(c**S(2) + S(1)))), CustomConstraint(lambda d, C, B, c: ZeroQ(-B*d + S(2)*C*c)))
    rule479 = ReplacementRule(pattern479, lambda p, n, B, c, d, a, f, C, m, x, e, A, b : Subst(Int((a + b*ArcCot(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule479)

    pattern480 = Pattern(Integral(ArcTan(a_ + x_*WC('b', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: RationalQ(n)))
    rule480 = ReplacementRule(pattern480, lambda n, c, d, a, x, b : ImaginaryI*Int(Log(-ImaginaryI*a - ImaginaryI*b*x + S(1))/(c + d*x**n), x)/S(2) - ImaginaryI*Int(Log(ImaginaryI*a + ImaginaryI*b*x + S(1))/(c + d*x**n), x)/S(2))
    rubi.add(rule480)

    pattern481 = Pattern(Integral(ArcCot(a_ + x_*WC('b', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: RationalQ(n)))
    rule481 = ReplacementRule(pattern481, lambda n, c, d, a, x, b : ImaginaryI*Int(Log((-ImaginaryI + a + b*x)/(a + b*x))/(c + d*x**n), x)/S(2) - ImaginaryI*Int(Log((ImaginaryI + a + b*x)/(a + b*x))/(c + d*x**n), x)/S(2))
    rubi.add(rule481)

    pattern482 = Pattern(Integral(ArcTan(a_ + x_*WC('b', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(RationalQ(n))))
    rule482 = ReplacementRule(pattern482, lambda n, c, d, a, x, b : Int(ArcTan(a + b*x)/(c + d*x**n), x))
    rubi.add(rule482)

    pattern483 = Pattern(Integral(ArcCot(a_ + x_*WC('b', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: ~(RationalQ(n))))
    rule483 = ReplacementRule(pattern483, lambda n, c, d, a, x, b : Int(ArcCot(a + b*x)/(c + d*x**n), x))
    rubi.add(rule483)

    pattern484 = Pattern(Integral(ArcTan(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule484 = ReplacementRule(pattern484, lambda a, b, n, x : -b*n*Int(x**n/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x) + x*ArcTan(a + b*x**n))
    rubi.add(rule484)

    pattern485 = Pattern(Integral(ArcCot(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule485 = ReplacementRule(pattern485, lambda a, b, n, x : b*n*Int(x**n/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x) + x*ArcCot(a + b*x**n))
    rubi.add(rule485)

    pattern486 = Pattern(Integral(ArcTan(x_**n_*WC('b', S(1)) + WC('a', S(0)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule486 = ReplacementRule(pattern486, lambda a, b, x, n : ImaginaryI*Int(Log(-ImaginaryI*a - ImaginaryI*b*x**n + S(1))/x, x)/S(2) - ImaginaryI*Int(Log(ImaginaryI*a + ImaginaryI*b*x**n + S(1))/x, x)/S(2))
    rubi.add(rule486)

    pattern487 = Pattern(Integral(ArcCot(x_**n_*WC('b', S(1)) + WC('a', S(0)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule487 = ReplacementRule(pattern487, lambda a, b, x, n : ImaginaryI*Int(Log(-ImaginaryI/(a + b*x**n) + S(1))/x, x)/S(2) - ImaginaryI*Int(Log(ImaginaryI/(a + b*x**n) + S(1))/x, x)/S(2))
    rubi.add(rule487)

    pattern488 = Pattern(Integral(x_**WC('m', S(1))*ArcTan(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda m: Unequal(m + S(1), S(0))), CustomConstraint(lambda n, m: Unequal(m + S(1), n)))
    rule488 = ReplacementRule(pattern488, lambda n, a, m, x, b : -b*n*Int(x**(m + n)/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x)/(m + S(1)) + x**(m + S(1))*ArcTan(a + b*x**n)/(m + S(1)))
    rubi.add(rule488)

    pattern489 = Pattern(Integral(x_**WC('m', S(1))*ArcCot(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda m: Unequal(m + S(1), S(0))), CustomConstraint(lambda n, m: Unequal(m + S(1), n)))
    rule489 = ReplacementRule(pattern489, lambda n, a, m, x, b : b*n*Int(x**(m + n)/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x)/(m + S(1)) + x**(m + S(1))*ArcCot(a + b*x**n)/(m + S(1)))
    rubi.add(rule489)

    pattern490 = Pattern(Integral(ArcTan(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule490 = ReplacementRule(pattern490, lambda c, d, a, f, x, b : ImaginaryI*Int(Log(-ImaginaryI*a - ImaginaryI*b*f**(c + d*x) + S(1)), x)/S(2) - ImaginaryI*Int(Log(ImaginaryI*a + ImaginaryI*b*f**(c + d*x) + S(1)), x)/S(2))
    rubi.add(rule490)

    pattern491 = Pattern(Integral(ArcCot(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule491 = ReplacementRule(pattern491, lambda c, d, a, f, x, b : ImaginaryI*Int(Log(-ImaginaryI/(a + b*f**(c + d*x)) + S(1)), x)/S(2) - ImaginaryI*Int(Log(ImaginaryI/(a + b*f**(c + d*x)) + S(1)), x)/S(2))
    rubi.add(rule491)

    pattern492 = Pattern(Integral(x_**WC('m', S(1))*ArcTan(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule492 = ReplacementRule(pattern492, lambda c, d, a, f, m, x, b : ImaginaryI*Int(x**m*Log(-ImaginaryI*a - ImaginaryI*b*f**(c + d*x) + S(1)), x)/S(2) - ImaginaryI*Int(x**m*Log(ImaginaryI*a + ImaginaryI*b*f**(c + d*x) + S(1)), x)/S(2))
    rubi.add(rule492)

    pattern493 = Pattern(Integral(x_**WC('m', S(1))*ArcCot(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule493 = ReplacementRule(pattern493, lambda c, d, a, f, m, x, b : ImaginaryI*Int(x**m*Log(-ImaginaryI/(a + b*f**(c + d*x)) + S(1)), x)/S(2) - ImaginaryI*Int(x**m*Log(ImaginaryI/(a + b*f**(c + d*x)) + S(1)), x)/S(2))
    rubi.add(rule493)

    pattern494 = Pattern(Integral(ArcTan(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule494 = ReplacementRule(pattern494, lambda n, c, u, a, m, x, b : Int(u*ArcCot(a/c + b*x**n/c)**m, x))
    rubi.add(rule494)

    pattern495 = Pattern(Integral(ArcCot(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule495 = ReplacementRule(pattern495, lambda n, c, u, a, m, x, b : Int(u*ArcTan(a/c + b*x**n/c)**m, x))
    rubi.add(rule495)

    pattern496 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0)))*ArcTan(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))))
    rule496 = ReplacementRule(pattern496, lambda a, b, c, x : Log(ArcTan(c*x/Sqrt(a + b*x**S(2))))/c)
    rubi.add(rule496)

    pattern497 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0)))*ArcCot(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))))
    rule497 = ReplacementRule(pattern497, lambda a, b, c, x : -Log(ArcCot(c*x/Sqrt(a + b*x**S(2))))/c)
    rubi.add(rule497)

    pattern498 = Pattern(Integral(ArcTan(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule498 = ReplacementRule(pattern498, lambda c, a, m, x, b : ArcTan(c*x/Sqrt(a + b*x**S(2)))**(m + S(1))/(c*(m + S(1))))
    rubi.add(rule498)

    pattern499 = Pattern(Integral(ArcCot(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule499 = ReplacementRule(pattern499, lambda c, a, m, x, b : -ArcCot(c*x/Sqrt(a + b*x**S(2)))**(m + S(1))/(c*(m + S(1))))
    rubi.add(rule499)

    pattern500 = Pattern(Integral(ArcTan(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))), CustomConstraint(lambda d, b, a, e: ZeroQ(-a*e + b*d)))
    rule500 = ReplacementRule(pattern500, lambda c, d, a, m, x, e, b : Int(ArcTan(c*x/Sqrt(a + b*x**S(2)))**m/Sqrt(a + b*x**S(2)), x)*Sqrt(a + b*x**S(2))/Sqrt(d + e*x**S(2)))
    rubi.add(rule500)

    pattern501 = Pattern(Integral(ArcCot(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))), CustomConstraint(lambda d, b, a, e: ZeroQ(-a*e + b*d)))
    rule501 = ReplacementRule(pattern501, lambda c, d, a, m, x, e, b : Int(ArcCot(c*x/Sqrt(a + b*x**S(2)))**m/Sqrt(a + b*x**S(2)), x)*Sqrt(a + b*x**S(2))/Sqrt(d + e*x**S(2)))
    rubi.add(rule501)

    pattern502 = Pattern(Integral(ArcTan(v_ + sqrt(w_)*WC('s', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda s: ZeroQ(s**S(2) + S(-1))), CustomConstraint(lambda w, v: ZeroQ(-v**S(2) + w + S(-1))))
    rule502 = ReplacementRule(pattern502, lambda u, s, x, w, v : Pi*s*Int(u, x)/S(4) + Int(u*ArcTan(v), x)/S(2))
    rubi.add(rule502)

    pattern503 = Pattern(Integral(ArcCot(v_ + sqrt(w_)*WC('s', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda s: ZeroQ(s**S(2) + S(-1))), CustomConstraint(lambda w, v: ZeroQ(-v**S(2) + w + S(-1))))
    rule503 = ReplacementRule(pattern503, lambda u, s, x, w, v : Pi*s*Int(u, x)/S(4) - Int(u*ArcTan(v), x)/S(2))
    rubi.add(rule503)

    pattern504 = Pattern(Integral(ArcTan(Tan(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: ZeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule504 = ReplacementRule(pattern504, lambda c, d, a, x, b : -ImaginaryI*b*Int(x/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*ArcTan(c + d*Tan(a + b*x)))
    rubi.add(rule504)

    pattern505 = Pattern(Integral(ArcCot(Tan(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: ZeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule505 = ReplacementRule(pattern505, lambda c, d, a, x, b : ImaginaryI*b*Int(x/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*ArcCot(c + d*Tan(a + b*x)))
    rubi.add(rule505)

    pattern506 = Pattern(Integral(ArcTan(Cot(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: ZeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule506 = ReplacementRule(pattern506, lambda c, d, a, x, b : -ImaginaryI*b*Int(x/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*ArcTan(c + d*Cot(a + b*x)))
    rubi.add(rule506)

    pattern507 = Pattern(Integral(ArcCot(Cot(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: ZeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule507 = ReplacementRule(pattern507, lambda c, d, a, x, b : ImaginaryI*b*Int(x/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*ArcCot(c + d*Cot(a + b*x)))
    rubi.add(rule507)

    pattern508 = Pattern(Integral(ArcTan(Tan(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: NonzeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule508 = ReplacementRule(pattern508, lambda c, d, a, x, b : b*(-ImaginaryI*c - d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c + d + (-ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) - b*(ImaginaryI*c + d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c - d + (ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*ArcTan(c + d*Tan(a + b*x)))
    rubi.add(rule508)

    pattern509 = Pattern(Integral(ArcCot(Tan(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: NonzeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule509 = ReplacementRule(pattern509, lambda c, d, a, x, b : -b*(-ImaginaryI*c - d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c + d + (-ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + b*(ImaginaryI*c + d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c - d + (ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*ArcCot(c + d*Tan(a + b*x)))
    rubi.add(rule509)

    pattern510 = Pattern(Integral(ArcTan(Cot(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: NonzeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule510 = ReplacementRule(pattern510, lambda c, d, a, x, b : -b*(-ImaginaryI*c + d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c - d - (-ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + b*(ImaginaryI*c - d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c + d - (ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*ArcTan(c + d*Cot(a + b*x)))
    rubi.add(rule510)

    pattern511 = Pattern(Integral(ArcCot(Cot(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: NonzeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule511 = ReplacementRule(pattern511, lambda c, d, a, x, b : b*(-ImaginaryI*c + d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c - d - (-ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) - b*(ImaginaryI*c - d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c + d - (ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*ArcCot(c + d*Cot(a + b*x)))
    rubi.add(rule511)

    pattern512 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(Tan(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: ZeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule512 = ReplacementRule(pattern512, lambda c, d, a, f, m, x, e, b : -ImaginaryI*b*Int((e + f*x)**(m + S(1))/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*Tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule512)

    pattern513 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcCot(Tan(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: ZeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule513 = ReplacementRule(pattern513, lambda c, d, a, f, m, x, e, b : ImaginaryI*b*Int((e + f*x)**(m + S(1))/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcCot(c + d*Tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule513)

    pattern514 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(Cot(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: ZeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule514 = ReplacementRule(pattern514, lambda c, d, a, f, m, x, e, b : -ImaginaryI*b*Int((e + f*x)**(m + S(1))/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*Cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule514)

    pattern515 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcCot(Cot(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: ZeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule515 = ReplacementRule(pattern515, lambda c, d, a, f, m, x, e, b : ImaginaryI*b*Int((e + f*x)**(m + S(1))/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcCot(c + d*Cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule515)

    pattern516 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(Tan(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: NonzeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule516 = ReplacementRule(pattern516, lambda c, d, a, f, m, x, e, b : b*(-ImaginaryI*c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c + d + (-ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) - b*(ImaginaryI*c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c - d + (ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*Tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule516)

    pattern517 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcCot(Tan(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: NonzeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule517 = ReplacementRule(pattern517, lambda c, d, a, f, m, x, e, b : -b*(-ImaginaryI*c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c + d + (-ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + b*(ImaginaryI*c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c - d + (ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcCot(c + d*Tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule517)

    pattern518 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(Cot(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: NonzeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule518 = ReplacementRule(pattern518, lambda c, d, a, f, m, x, e, b : -b*(-ImaginaryI*c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c - d - (-ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + b*(ImaginaryI*c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c + d - (ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*Cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule518)

    pattern519 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcCot(Cot(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: NonzeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule519 = ReplacementRule(pattern519, lambda c, d, a, f, m, x, e, b : b*(-ImaginaryI*c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c - d - (-ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) - b*(ImaginaryI*c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c + d - (ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcCot(c + d*Cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule519)

    pattern520 = Pattern(Integral(ArcTan(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule520 = ReplacementRule(pattern520, lambda a, b, x : -b*Int(x*Sech(S(2)*a + S(2)*b*x), x) + x*ArcTan(Tanh(a + b*x)))
    rubi.add(rule520)

    pattern521 = Pattern(Integral(ArcCot(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule521 = ReplacementRule(pattern521, lambda a, b, x : b*Int(x*Sech(S(2)*a + S(2)*b*x), x) + x*ArcCot(Tanh(a + b*x)))
    rubi.add(rule521)

    pattern522 = Pattern(Integral(ArcTan(Coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule522 = ReplacementRule(pattern522, lambda a, b, x : b*Int(x*Sech(S(2)*a + S(2)*b*x), x) + x*ArcTan(Coth(a + b*x)))
    rubi.add(rule522)

    pattern523 = Pattern(Integral(ArcCot(Coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule523 = ReplacementRule(pattern523, lambda a, b, x : -b*Int(x*Sech(S(2)*a + S(2)*b*x), x) + x*ArcCot(Coth(a + b*x)))
    rubi.add(rule523)

    pattern524 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule524 = ReplacementRule(pattern524, lambda a, f, m, x, e, b : -b*Int((e + f*x)**(m + S(1))*Sech(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(Tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule524)

    pattern525 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcCot(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule525 = ReplacementRule(pattern525, lambda a, f, m, x, e, b : b*Int((e + f*x)**(m + S(1))*Sech(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcCot(Tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule525)

    pattern526 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(Coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule526 = ReplacementRule(pattern526, lambda a, f, m, x, e, b : b*Int((e + f*x)**(m + S(1))*Sech(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(Coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule526)

    pattern527 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcCot(Coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule527 = ReplacementRule(pattern527, lambda a, f, m, x, e, b : -b*Int((e + f*x)**(m + S(1))*Sech(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcCot(Coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule527)

    pattern528 = Pattern(Integral(ArcTan(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: ZeroQ((c - d)**S(2) + S(1))))
    rule528 = ReplacementRule(pattern528, lambda c, d, a, x, b : -b*Int(x/(c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*ArcTan(c + d*Tanh(a + b*x)))
    rubi.add(rule528)

    pattern529 = Pattern(Integral(ArcCot(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: ZeroQ((c - d)**S(2) + S(1))))
    rule529 = ReplacementRule(pattern529, lambda c, d, a, x, b : b*Int(x/(c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*ArcCot(c + d*Tanh(a + b*x)))
    rubi.add(rule529)

    pattern530 = Pattern(Integral(ArcTan(Coth(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: ZeroQ((c - d)**S(2) + S(1))))
    rule530 = ReplacementRule(pattern530, lambda c, d, a, x, b : -b*Int(x/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*ArcTan(c + d*Coth(a + b*x)))
    rubi.add(rule530)

    pattern531 = Pattern(Integral(ArcCot(Coth(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: ZeroQ((c - d)**S(2) + S(1))))
    rule531 = ReplacementRule(pattern531, lambda c, d, a, x, b : b*Int(x/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*ArcCot(c + d*Coth(a + b*x)))
    rubi.add(rule531)

    pattern532 = Pattern(Integral(ArcTan(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: NonzeroQ((c - d)**S(2) + S(1))))
    rule532 = ReplacementRule(pattern532, lambda c, d, a, x, b : ImaginaryI*b*(ImaginaryI - c - d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d + (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x) - ImaginaryI*b*(ImaginaryI + c + d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d + (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x) + x*ArcTan(c + d*Tanh(a + b*x)))
    rubi.add(rule532)

    pattern533 = Pattern(Integral(ArcCot(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: NonzeroQ((c - d)**S(2) + S(1))))
    rule533 = ReplacementRule(pattern533, lambda c, d, a, x, b : -ImaginaryI*b*(ImaginaryI - c - d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d + (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x) + ImaginaryI*b*(ImaginaryI + c + d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d + (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x) + x*ArcCot(c + d*Tanh(a + b*x)))
    rubi.add(rule533)

    pattern534 = Pattern(Integral(ArcTan(Coth(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: NonzeroQ((c - d)**S(2) + S(1))))
    rule534 = ReplacementRule(pattern534, lambda c, d, a, x, b : -ImaginaryI*b*(ImaginaryI - c - d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d - (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x) + ImaginaryI*b*(ImaginaryI + c + d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d - (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x) + x*ArcTan(c + d*Coth(a + b*x)))
    rubi.add(rule534)

    pattern535 = Pattern(Integral(ArcCot(Coth(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, c: NonzeroQ((c - d)**S(2) + S(1))))
    rule535 = ReplacementRule(pattern535, lambda c, d, a, x, b : ImaginaryI*b*(ImaginaryI - c - d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d - (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x) - ImaginaryI*b*(ImaginaryI + c + d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d - (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x) + x*ArcCot(c + d*Coth(a + b*x)))
    rubi.add(rule535)

    pattern536 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: ZeroQ((c - d)**S(2) + S(1))))
    rule536 = ReplacementRule(pattern536, lambda c, d, a, f, m, x, e, b : -b*Int((e + f*x)**(m + S(1))/(c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*Tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule536)

    pattern537 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcCot(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: ZeroQ((c - d)**S(2) + S(1))))
    rule537 = ReplacementRule(pattern537, lambda c, d, a, f, m, x, e, b : b*Int((e + f*x)**(m + S(1))/(c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcCot(c + d*Tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule537)

    pattern538 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(Coth(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: ZeroQ((c - d)**S(2) + S(1))))
    rule538 = ReplacementRule(pattern538, lambda c, d, a, f, m, x, e, b : -b*Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*Coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule538)

    pattern539 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcCot(Coth(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: ZeroQ((c - d)**S(2) + S(1))))
    rule539 = ReplacementRule(pattern539, lambda c, d, a, f, m, x, e, b : b*Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcCot(c + d*Coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule539)

    pattern540 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: NonzeroQ((c - d)**S(2) + S(1))))
    rule540 = ReplacementRule(pattern540, lambda c, d, a, f, m, x, e, b : ImaginaryI*b*(ImaginaryI - c - d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d + (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) - ImaginaryI*b*(ImaginaryI + c + d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d + (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*Tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule540)

    pattern541 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcCot(Tanh(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: NonzeroQ((c - d)**S(2) + S(1))))
    rule541 = ReplacementRule(pattern541, lambda c, d, a, f, m, x, e, b : -ImaginaryI*b*(ImaginaryI - c - d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d + (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + ImaginaryI*b*(ImaginaryI + c + d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d + (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcCot(c + d*Tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule541)

    pattern542 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(Coth(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: NonzeroQ((c - d)**S(2) + S(1))))
    rule542 = ReplacementRule(pattern542, lambda c, d, a, f, m, x, e, b : -ImaginaryI*b*(ImaginaryI - c - d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d - (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + ImaginaryI*b*(ImaginaryI + c + d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d - (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*Coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule542)

    pattern543 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcCot(Coth(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, c: NonzeroQ((c - d)**S(2) + S(1))))
    rule543 = ReplacementRule(pattern543, lambda c, d, a, f, m, x, e, b : ImaginaryI*b*(ImaginaryI - c - d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d - (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) - ImaginaryI*b*(ImaginaryI + c + d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d - (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcCot(c + d*Coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule543)

    pattern544 = Pattern(Integral(ArcTan(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)))
    rule544 = ReplacementRule(pattern544, lambda u, x : x*ArcTan(u) - Int(SimplifyIntegrand(x*D(u, x)/(u**S(2) + S(1)), x), x))
    rubi.add(rule544)

    pattern545 = Pattern(Integral(ArcCot(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)))
    rule545 = ReplacementRule(pattern545, lambda u, x : x*ArcCot(u) + Int(SimplifyIntegrand(x*D(u, x)/(u**S(2) + S(1)), x), x))
    rubi.add(rule545)

    pattern546 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(ArcTan(u_)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, d, u, m, x: ~(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, m, x: FalseQ(PowerVariableExpn(u, m + S(1), x))))
    rule546 = ReplacementRule(pattern546, lambda c, d, a, u, m, x, b : -b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u**S(2) + S(1)), x), x)/(d*(m + S(1))) + (a + b*ArcTan(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule546)

    pattern547 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(ArcCot(u_)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, d, u, m, x: ~(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, m, x: FalseQ(PowerVariableExpn(u, m + S(1), x))))
    rule547 = ReplacementRule(pattern547, lambda c, d, a, u, m, x, b : b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u**S(2) + S(1)), x), x)/(d*(m + S(1))) + (a + b*ArcCot(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule547)

    pattern548 = Pattern(Integral(ArcTan(v_)*Log(w_)/(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda x, v: LinearQ(v, x)), CustomConstraint(lambda w, x: LinearQ(w, x)), CustomConstraint(lambda a, b, x, v: ZeroQ(Simplify(D(v/(a + b*x), x)))), CustomConstraint(lambda a, b, w, x: ZeroQ(Simplify(D(w/(a + b*x), x)))))
    rule548 = ReplacementRule(pattern548, lambda a, x, w, v, b : ImaginaryI*Int(Log(w)*Log(-ImaginaryI*v + S(1))/(a + b*x), x)/S(2) - ImaginaryI*Int(Log(w)*Log(ImaginaryI*v + S(1))/(a + b*x), x)/S(2))
    rubi.add(rule548)

    pattern549 = Pattern(Integral(ArcTan(v_)*Log(w_), x_), CustomConstraint(lambda x, v: InverseFunctionFreeQ(v, x)), CustomConstraint(lambda w, x: InverseFunctionFreeQ(w, x)))
    rule549 = ReplacementRule(pattern549, lambda w, x, v : x*ArcTan(v)*Log(w) - Int(SimplifyIntegrand(x*ArcTan(v)*D(w, x)/w, x), x) - Int(SimplifyIntegrand(x*D(v, x)*Log(w)/(v**S(2) + S(1)), x), x))
    rubi.add(rule549)

    pattern550 = Pattern(Integral(ArcCot(v_)*Log(w_), x_), CustomConstraint(lambda x, v: InverseFunctionFreeQ(v, x)), CustomConstraint(lambda w, x: InverseFunctionFreeQ(w, x)))
    rule550 = ReplacementRule(pattern550, lambda w, x, v : x*ArcCot(v)*Log(w) - Int(SimplifyIntegrand(x*ArcCot(v)*D(w, x)/w, x), x) + Int(SimplifyIntegrand(x*D(v, x)*Log(w)/(v**S(2) + S(1)), x), x))
    rubi.add(rule550)

    pattern551 = Pattern(Integral(u_*ArcTan(v_)*Log(w_), x_), CustomConstraint(lambda x, v: InverseFunctionFreeQ(v, x)), CustomConstraint(lambda w, x: InverseFunctionFreeQ(w, x)))
    rule551 = ReplacementRule(pattern551, lambda u, w, x, v : With(List(Set(z, IntHide(u, x))), Condition(Dist(ArcTan(v)*Log(w), z, x) - Int(SimplifyIntegrand(z*ArcTan(v)*D(w, x)/w, x), x) - Int(SimplifyIntegrand(z*D(v, x)*Log(w)/(v**S(2) + S(1)), x), x), InverseFunctionFreeQ(z, x))))
    rubi.add(rule551)

    pattern552 = Pattern(Integral(u_*ArcCot(v_)*Log(w_), x_), CustomConstraint(lambda x, v: InverseFunctionFreeQ(v, x)), CustomConstraint(lambda w, x: InverseFunctionFreeQ(w, x)))
    rule552 = ReplacementRule(pattern552, lambda u, w, x, v : With(List(Set(z, IntHide(u, x))), Condition(Dist(ArcCot(v)*Log(w), z, x) - Int(SimplifyIntegrand(z*ArcCot(v)*D(w, x)/w, x), x) + Int(SimplifyIntegrand(z*D(v, x)*Log(w)/(v**S(2) + S(1)), x), x), InverseFunctionFreeQ(z, x))))
    rubi.add(rule552)

    pattern553 = Pattern(Integral(ArcSec(x_*WC('c', S(1))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule553 = ReplacementRule(pattern553, lambda c, x : x*ArcSec(c*x) - Int(S(1)/(x*Sqrt(S(1) - S(1)/(c**S(2)*x**S(2)))), x)/c)
    rubi.add(rule553)

    pattern554 = Pattern(Integral(ArcCsc(x_*WC('c', S(1))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule554 = ReplacementRule(pattern554, lambda c, x : x*ArcCsc(c*x) + Int(S(1)/(x*Sqrt(S(1) - S(1)/(c**S(2)*x**S(2)))), x)/c)
    rubi.add(rule554)

    pattern555 = Pattern(Integral((ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule555 = ReplacementRule(pattern555, lambda n, c, a, x, b : Subst(Int((a + b*x)**n*Sec(x)*Tan(x), x), x, ArcSec(c*x))/c)
    rubi.add(rule555)

    pattern556 = Pattern(Integral((ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule556 = ReplacementRule(pattern556, lambda n, c, a, x, b : -Subst(Int((a + b*x)**n*Cot(x)*Csc(x), x), x, ArcCsc(c*x))/c)
    rubi.add(rule556)

    pattern557 = Pattern(Integral((ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule557 = ReplacementRule(pattern557, lambda a, b, c, x : -Subst(Int((a + b*ArcCos(x/c))/x, x), x, 1/x))
    rubi.add(rule557)

    pattern558 = Pattern(Integral((ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule558 = ReplacementRule(pattern558, lambda a, b, c, x : -Subst(Int((a + b*ArcSin(x/c))/x, x), x, 1/x))
    rubi.add(rule558)

    pattern559 = Pattern(Integral(x_**WC('m', S(1))*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule559 = ReplacementRule(pattern559, lambda c, a, m, x, b : -b*Int(x**(m + S(-1))/Sqrt(S(1) - S(1)/(c**S(2)*x**S(2))), x)/(c*(m + S(1))) + x**(m + S(1))*(a + b*ArcSec(c*x))/(m + S(1)))
    rubi.add(rule559)

    pattern560 = Pattern(Integral(x_**WC('m', S(1))*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule560 = ReplacementRule(pattern560, lambda c, a, m, x, b : b*Int(x**(m + S(-1))/Sqrt(S(1) - S(1)/(c**S(2)*x**S(2))), x)/(c*(m + S(1))) + x**(m + S(1))*(a + b*ArcCsc(c*x))/(m + S(1)))
    rubi.add(rule560)

    pattern561 = Pattern(Integral(x_**WC('m', S(1))*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule561 = ReplacementRule(pattern561, lambda n, c, a, m, x, b : c**(-m + S(-1))*Subst(Int((a + b*x)**n*Sec(x)**(m + S(1))*Tan(x), x), x, ArcSec(c*x)))
    rubi.add(rule561)

    pattern562 = Pattern(Integral(x_**WC('m', S(1))*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule562 = ReplacementRule(pattern562, lambda n, c, a, m, x, b : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*Cot(x)*Csc(x)**(m + S(1)), x), x, ArcCsc(c*x)))
    rubi.add(rule562)

    pattern563 = Pattern(Integral(x_**WC('m', S(1))*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule563 = ReplacementRule(pattern563, lambda n, c, a, m, x, b : Int(x**m*(a + b*ArcSec(c*x))**n, x))
    rubi.add(rule563)

    pattern564 = Pattern(Integral(x_**WC('m', S(1))*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule564 = ReplacementRule(pattern564, lambda n, c, a, m, x, b : Int(x**m*(a + b*ArcCsc(c*x))**n, x))
    rubi.add(rule564)

    pattern565 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)))
    rule565 = ReplacementRule(pattern565, lambda p, c, d, a, x, e, b : With(List(Set(u, IntHide((d + e*x**S(2))**p, x))), -b*c*x*Int(SimplifyIntegrand(u/(x*Sqrt(c**S(2)*x**S(2) + S(-1))), x), x)/Sqrt(c**S(2)*x**S(2)) + Dist(a + b*ArcSec(c*x), u, x)))
    rubi.add(rule565)

    pattern566 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)))
    rule566 = ReplacementRule(pattern566, lambda p, c, d, a, x, e, b : With(List(Set(u, IntHide((d + e*x**S(2))**p, x))), b*c*x*Int(SimplifyIntegrand(u/(x*Sqrt(c**S(2)*x**S(2) + S(-1))), x), x)/Sqrt(c**S(2)*x**S(2)) + Dist(a + b*ArcCsc(c*x), u, x)))
    rubi.add(rule566)

    pattern567 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)))
    rule567 = ReplacementRule(pattern567, lambda p, n, c, d, a, x, e, b : -Subst(Int(x**(-S(2)*p + S(-2))*(a + b*ArcCos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule567)

    pattern568 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)))
    rule568 = ReplacementRule(pattern568, lambda p, n, c, d, a, x, e, b : -Subst(Int(x**(-S(2)*p + S(-2))*(a + b*ArcSin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule568)

    pattern569 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule569 = ReplacementRule(pattern569, lambda p, n, c, d, a, x, e, b : -Sqrt(x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*ArcCos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule569)

    pattern570 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule570 = ReplacementRule(pattern570, lambda p, n, c, d, a, x, e, b : -Sqrt(x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*ArcSin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule570)

    pattern571 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d, e: ~(Negative(d) & PositiveQ(e))))
    rule571 = ReplacementRule(pattern571, lambda p, n, c, d, a, x, e, b : -Sqrt(d + e*x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*ArcCos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*Sqrt(d/x**S(2) + e)))
    rubi.add(rule571)

    pattern572 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d, e: ~(Negative(d) & PositiveQ(e))))
    rule572 = ReplacementRule(pattern572, lambda p, n, c, d, a, x, e, b : -Sqrt(d + e*x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*ArcSin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*Sqrt(d/x**S(2) + e)))
    rubi.add(rule572)

    pattern573 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule573 = ReplacementRule(pattern573, lambda p, n, c, d, a, x, e, b : Int((a + b*ArcSec(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule573)

    pattern574 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule574 = ReplacementRule(pattern574, lambda p, n, c, d, a, x, e, b : Int((a + b*ArcCsc(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule574)

    pattern575 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule575 = ReplacementRule(pattern575, lambda p, c, d, a, x, e, b : -b*c*x*Int((d + e*x**S(2))**(p + S(1))/(x*Sqrt(c**S(2)*x**S(2) + S(-1))), x)/(S(2)*e*(p + S(1))*Sqrt(c**S(2)*x**S(2))) + (a + b*ArcSec(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule575)

    pattern576 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule576 = ReplacementRule(pattern576, lambda p, c, d, a, x, e, b : b*c*x*Int((d + e*x**S(2))**(p + S(1))/(x*Sqrt(c**S(2)*x**S(2) + S(-1))), x)/(S(2)*e*(p + S(1))*Sqrt(c**S(2)*x**S(2))) + (a + b*ArcCsc(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule576)

    pattern577 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & ~(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & ~(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & ~(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))))
    rule577 = ReplacementRule(pattern577, lambda p, c, d, a, m, x, e, b : With(List(Set(u, IntHide(x**m*(d + e*x**S(2))**p, x))), -b*c*x*Int(SimplifyIntegrand(u/(x*Sqrt(c**S(2)*x**S(2) + S(-1))), x), x)/Sqrt(c**S(2)*x**S(2)) + Dist(a + b*ArcSec(c*x), u, x)))
    rubi.add(rule577)

    pattern578 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & ~(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & ~(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & ~(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))))
    rule578 = ReplacementRule(pattern578, lambda p, c, d, a, m, x, e, b : With(List(Set(u, IntHide(x**m*(d + e*x**S(2))**p, x))), b*c*x*Int(SimplifyIntegrand(u/(x*Sqrt(c**S(2)*x**S(2) + S(-1))), x), x)/Sqrt(c**S(2)*x**S(2)) + Dist(a + b*ArcCsc(c*x), u, x)))
    rubi.add(rule578)

    pattern579 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, m: IntegersQ(m, p)))
    rule579 = ReplacementRule(pattern579, lambda p, n, c, d, a, m, x, e, b : -Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*ArcCos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule579)

    pattern580 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, m: IntegersQ(m, p)))
    rule580 = ReplacementRule(pattern580, lambda p, n, c, d, a, m, x, e, b : -Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*ArcSin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule580)

    pattern581 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule581 = ReplacementRule(pattern581, lambda p, n, c, d, a, m, x, e, b : -Sqrt(x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*ArcCos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule581)

    pattern582 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule582 = ReplacementRule(pattern582, lambda p, n, c, d, a, m, x, e, b : -Sqrt(x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*ArcSin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule582)

    pattern583 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d, e: ~(Negative(d) & PositiveQ(e))))
    rule583 = ReplacementRule(pattern583, lambda p, n, c, d, a, m, x, e, b : -Sqrt(d + e*x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*ArcCos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*Sqrt(d/x**S(2) + e)))
    rubi.add(rule583)

    pattern584 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda d, c, e: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d, e: ~(Negative(d) & PositiveQ(e))))
    rule584 = ReplacementRule(pattern584, lambda p, n, c, d, a, m, x, e, b : -Sqrt(d + e*x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*ArcSin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*Sqrt(d/x**S(2) + e)))
    rubi.add(rule584)

    pattern585 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcSec(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule585 = ReplacementRule(pattern585, lambda p, n, c, d, a, m, x, e, b : Int(x**m*(a + b*ArcSec(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule585)

    pattern586 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcCsc(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule586 = ReplacementRule(pattern586, lambda p, n, c, d, a, m, x, e, b : Int(x**m*(a + b*ArcCsc(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule586)

    pattern587 = Pattern(Integral(ArcSec(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule587 = ReplacementRule(pattern587, lambda a, b, x : -Int(S(1)/((a + b*x)*Sqrt(S(1) - S(1)/(a + b*x)**S(2))), x) + (a + b*x)*ArcSec(a + b*x)/b)
    rubi.add(rule587)

    pattern588 = Pattern(Integral(ArcCsc(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule588 = ReplacementRule(pattern588, lambda a, b, x : Int(S(1)/((a + b*x)*Sqrt(S(1) - S(1)/(a + b*x)**S(2))), x) + (a + b*x)*ArcCsc(a + b*x)/b)
    rubi.add(rule588)

    pattern589 = Pattern(Integral(ArcSec(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule589 = ReplacementRule(pattern589, lambda a, b, n, x : Subst(Int(x**n*Sec(x)*Tan(x), x), x, ArcSec(a + b*x))/b)
    rubi.add(rule589)

    pattern590 = Pattern(Integral(ArcCsc(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule590 = ReplacementRule(pattern590, lambda a, b, n, x : -Subst(Int(x**n*Cot(x)*Csc(x), x), x, ArcCsc(a + b*x))/b)
    rubi.add(rule590)

    pattern591 = Pattern(Integral(ArcSec(a_ + x_*WC('b', S(1)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule591 = ReplacementRule(pattern591, lambda a, b, x : -ImaginaryI*PolyLog(S(2), (-Sqrt(-a**S(2) + S(1)) + S(1))*exp(ImaginaryI*ArcSec(a + b*x))/a) - ImaginaryI*PolyLog(S(2), (Sqrt(-a**S(2) + S(1)) + S(1))*exp(ImaginaryI*ArcSec(a + b*x))/a) + ImaginaryI*PolyLog(S(2), -exp(S(2)*ImaginaryI*ArcSec(a + b*x)))/S(2) + ArcSec(a + b*x)*Log(S(1) - (-Sqrt(-a**S(2) + S(1)) + S(1))*exp(ImaginaryI*ArcSec(a + b*x))/a) + ArcSec(a + b*x)*Log(S(1) - (Sqrt(-a**S(2) + S(1)) + S(1))*exp(ImaginaryI*ArcSec(a + b*x))/a) - ArcSec(a + b*x)*Log(exp(S(2)*ImaginaryI*ArcSec(a + b*x)) + S(1)))
    rubi.add(rule591)

    pattern592 = Pattern(Integral(ArcCsc(a_ + x_*WC('b', S(1)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule592 = ReplacementRule(pattern592, lambda a, b, x : ImaginaryI*ArcCsc(a + b*x)**S(2) + ImaginaryI*PolyLog(S(2), ImaginaryI*(-Sqrt(-a**S(2) + S(1)) + S(1))*exp(-ImaginaryI*ArcCsc(a + b*x))/a) + ImaginaryI*PolyLog(S(2), ImaginaryI*(Sqrt(-a**S(2) + S(1)) + S(1))*exp(-ImaginaryI*ArcCsc(a + b*x))/a) + ImaginaryI*PolyLog(S(2), exp(S(2)*ImaginaryI*ArcCsc(a + b*x)))/S(2) + ArcCsc(a + b*x)*Log(-ImaginaryI*(-Sqrt(-a**S(2) + S(1)) + S(1))*exp(-ImaginaryI*ArcCsc(a + b*x))/a + S(1)) + ArcCsc(a + b*x)*Log(-ImaginaryI*(Sqrt(-a**S(2) + S(1)) + S(1))*exp(-ImaginaryI*ArcCsc(a + b*x))/a + S(1)) - ArcCsc(a + b*x)*Log(-exp(S(2)*ImaginaryI*ArcCsc(a + b*x)) + S(1)))
    rubi.add(rule592)

    pattern593 = Pattern(Integral(x_**WC('m', S(1))*ArcSec(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule593 = ReplacementRule(pattern593, lambda a, b, m, x : b**(-m + S(-1))*(b**(m + S(1))*x**(m + S(1)) - (-a)**(m + S(1)))*ArcSec(a + b*x)/(m + S(1)) - b**(-m + S(-1))*Subst(Int(x**(-m + S(-1))*((-a*x)**(m + S(1)) - (-a*x + S(1))**(m + S(1)))/Sqrt(-x**S(2) + S(1)), x), x, 1/(a + b*x))/(m + S(1)))
    rubi.add(rule593)

    pattern594 = Pattern(Integral(x_**WC('m', S(1))*ArcCsc(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule594 = ReplacementRule(pattern594, lambda a, b, m, x : b**(-m + S(-1))*(b**(m + S(1))*x**(m + S(1)) - (-a)**(m + S(1)))*ArcCsc(a + b*x)/(m + S(1)) + b**(-m + S(-1))*Subst(Int(x**(-m + S(-1))*((-a*x)**(m + S(1)) - (-a*x + S(1))**(m + S(1)))/Sqrt(-x**S(2) + S(1)), x), x, 1/(a + b*x))/(m + S(1)))
    rubi.add(rule594)

    pattern595 = Pattern(Integral(x_**WC('m', S(1))*ArcSec(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule595 = ReplacementRule(pattern595, lambda n, a, m, x, b : b**(-m + S(-1))*Subst(Int(x**n*(-a + Sec(x))**m*Sec(x)*Tan(x), x), x, ArcSec(a + b*x)))
    rubi.add(rule595)

    pattern596 = Pattern(Integral(x_**WC('m', S(1))*ArcCsc(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule596 = ReplacementRule(pattern596, lambda n, a, m, x, b : -b**(-m + S(-1))*Subst(Int(x**n*(-a + Csc(x))**m*Cot(x)*Csc(x), x), x, ArcCsc(a + b*x)))
    rubi.add(rule596)

    pattern597 = Pattern(Integral(ArcSec(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule597 = ReplacementRule(pattern597, lambda n, c, u, a, m, x, b : Int(u*ArcCos(a/c + b*x**n/c)**m, x))
    rubi.add(rule597)

    pattern598 = Pattern(Integral(ArcCsc(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule598 = ReplacementRule(pattern598, lambda n, c, u, a, m, x, b : Int(u*ArcSin(a/c + b*x**n/c)**m, x))
    rubi.add(rule598)

    pattern599 = Pattern(Integral(f_**(ArcSec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*WC('c', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule599 = ReplacementRule(pattern599, lambda n, c, u, a, f, x, b : Subst(Int(f**(c*x**n)*ReplaceAll(u, Rule(x, -a/b + Sec(x)/b))*Sec(x)*Tan(x), x), x, ArcSec(a + b*x))/b)
    rubi.add(rule599)

    pattern600 = Pattern(Integral(f_**(ArcCsc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*WC('c', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule600 = ReplacementRule(pattern600, lambda n, c, u, a, f, x, b : -Subst(Int(f**(c*x**n)*Cot(x)*Csc(x)*ReplaceAll(u, Rule(x, -a/b + Csc(x)/b)), x), x, ArcCsc(a + b*x))/b)
    rubi.add(rule600)

    pattern601 = Pattern(Integral(ArcSec(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda u, x: ~(FunctionOfExponentialQ(u, x))))
    rule601 = ReplacementRule(pattern601, lambda u, x : -u*Int(SimplifyIntegrand(x*D(u, x)/(u*Sqrt(u**S(2) + S(-1))), x), x)/Sqrt(u**S(2)) + x*ArcSec(u))
    rubi.add(rule601)

    pattern602 = Pattern(Integral(ArcCsc(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda u, x: ~(FunctionOfExponentialQ(u, x))))
    rule602 = ReplacementRule(pattern602, lambda u, x : u*Int(SimplifyIntegrand(x*D(u, x)/(u*Sqrt(u**S(2) + S(-1))), x), x)/Sqrt(u**S(2)) + x*ArcCsc(u))
    rubi.add(rule602)

    pattern603 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(ArcSec(u_)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, d, u, m, x: ~(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, x: ~(FunctionOfExponentialQ(u, x))))
    rule603 = ReplacementRule(pattern603, lambda c, d, a, u, m, x, b : -b*u*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u*Sqrt(u**S(2) + S(-1))), x), x)/(d*(m + S(1))*Sqrt(u**S(2))) + (a + b*ArcSec(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule603)

    pattern604 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(ArcCsc(u_)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, d, u, m, x: ~(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, x: ~(FunctionOfExponentialQ(u, x))))
    rule604 = ReplacementRule(pattern604, lambda c, d, a, u, m, x, b : b*u*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u*Sqrt(u**S(2) + S(-1))), x), x)/(d*(m + S(1))*Sqrt(u**S(2))) + (a + b*ArcCsc(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule604)

    return rubi
