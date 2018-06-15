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

def inverse_hyperbolic(rubi):

    pattern1 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule1 = ReplacementRule(pattern1, lambda x, b, c, a, n : -b*c*n*Int(x*(a + b*asinh(c*x))**(n + S(-1))/sqrt(c**S(2)*x**S(2) + S(1)), x) + x*(a + b*asinh(c*x))**n)
    rubi.add(rule1)

    pattern2 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule2 = ReplacementRule(pattern2, lambda x, b, c, a, n : -b*c*n*Int(x*(a + b*acosh(c*x))**(n + S(-1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x) + x*(a + b*acosh(c*x))**n)
    rubi.add(rule2)

    pattern3 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule3 = ReplacementRule(pattern3, lambda x, b, c, a, n : -c*Int(x*(a + b*asinh(c*x))**(n + S(1))/sqrt(c**S(2)*x**S(2) + S(1)), x)/(b*(n + S(1))) + (a + b*asinh(c*x))**(n + S(1))*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule3)

    pattern4 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule4 = ReplacementRule(pattern4, lambda x, b, c, a, n : -c*Int(x*(a + b*acosh(c*x))**(n + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(b*(n + S(1))) + (a + b*acosh(c*x))**(n + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))))
    rubi.add(rule4)

    pattern5 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule5 = ReplacementRule(pattern5, lambda x, b, c, a, n : Subst(Int(x**n*Cosh(a/b - x/b), x), x, a + b*asinh(c*x))/(b*c))
    rubi.add(rule5)

    pattern6 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule6 = ReplacementRule(pattern6, lambda x, b, c, a, n : -Subst(Int(x**n*sinh(a/b - x/b), x), x, a + b*acosh(c*x))/(b*c))
    rubi.add(rule6)

    pattern7 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule7 = ReplacementRule(pattern7, lambda x, b, c, a, n : Subst(Int((a + b*x)**n/tanh(x), x), x, asinh(c*x)))
    rubi.add(rule7)

    pattern8 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule8 = ReplacementRule(pattern8, lambda x, b, c, a, n : Subst(Int((a + b*x)**n/coth(x), x), x, acosh(c*x)))
    rubi.add(rule8)

    pattern9 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule9 = ReplacementRule(pattern9, lambda x, b, m, c, a, d, n : -b*c*n*Int((d*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))/sqrt(c**S(2)*x**S(2) + S(1)), x)/(d*(m + S(1))) + (d*x)**(m + S(1))*(a + b*asinh(c*x))**n/(d*(m + S(1))))
    rubi.add(rule9)

    pattern10 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule10 = ReplacementRule(pattern10, lambda x, b, m, c, a, d, n : -b*c*n*Int((d*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(d*(m + S(1))) + (d*x)**(m + S(1))*(a + b*acosh(c*x))**n/(d*(m + S(1))))
    rubi.add(rule10)

    pattern11 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule11 = ReplacementRule(pattern11, lambda x, b, m, c, a, n : -b*c*n*Int(x**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))/sqrt(c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*asinh(c*x))**n/(m + S(1)))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule12 = ReplacementRule(pattern12, lambda x, b, m, c, a, n : -b*c*n*Int(x**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(m + S(1)) + x**(m + S(1))*(a + b*acosh(c*x))**n/(m + S(1)))
    rubi.add(rule12)

    pattern13 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Inequality(S(-2), LessEqual, n, Less, S(-1))))
    rule13 = ReplacementRule(pattern13, lambda x, b, m, c, a, n : -c**(-m + S(-1))*Subst(Int(ExpandTrigReduce((a + b*x)**(n + S(1)), (m + (m + S(1))*sinh(x)**S(2))*sinh(x)**(m + S(-1)), x), x), x, asinh(c*x))/(b*(n + S(1))) + x**m*(a + b*asinh(c*x))**(n + S(1))*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule13)

    pattern14 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Inequality(S(-2), LessEqual, n, Less, S(-1))))
    rule14 = ReplacementRule(pattern14, lambda x, b, m, c, a, n : c**(-m + S(-1))*Subst(Int(ExpandTrigReduce((a + b*x)**(n + S(1))*(m - (m + S(1))*Cosh(x)**S(2))*Cosh(x)**(m + S(-1)), x), x), x, acosh(c*x))/(b*(n + S(1))) + x**m*(a + b*acosh(c*x))**(n + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))))
    rubi.add(rule14)

    pattern15 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-2))))
    rule15 = ReplacementRule(pattern15, lambda x, b, m, c, a, n : -c*(m + S(1))*Int(x**(m + S(1))*(a + b*asinh(c*x))**(n + S(1))/sqrt(c**S(2)*x**S(2) + S(1)), x)/(b*(n + S(1))) - m*Int(x**(m + S(-1))*(a + b*asinh(c*x))**(n + S(1))/sqrt(c**S(2)*x**S(2) + S(1)), x)/(b*c*(n + S(1))) + x**m*(a + b*asinh(c*x))**(n + S(1))*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule15)

    pattern16 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-2))))
    rule16 = ReplacementRule(pattern16, lambda x, b, m, c, a, n : -c*(m + S(1))*Int(x**(m + S(1))*(a + b*acosh(c*x))**(n + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(b*(n + S(1))) + m*Int(x**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(b*c*(n + S(1))) + x**m*(a + b*acosh(c*x))**(n + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))))
    rubi.add(rule16)

    pattern17 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule17 = ReplacementRule(pattern17, lambda x, b, m, c, a, n : c**(-m + S(-1))*Subst(Int((a + b*x)**n*Cosh(x)*sinh(x)**m, x), x, asinh(c*x)))
    rubi.add(rule17)

    pattern18 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule18 = ReplacementRule(pattern18, lambda x, b, m, c, a, n : c**(-m + S(-1))*Subst(Int((a + b*x)**n*Cosh(x)**m*sinh(x), x), x, acosh(c*x)))
    rubi.add(rule18)

    pattern19 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule19 = ReplacementRule(pattern19, lambda x, b, m, c, a, d, n : Int((d*x)**m*(a + b*asinh(c*x))**n, x))
    rubi.add(rule19)

    pattern20 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule20 = ReplacementRule(pattern20, lambda x, b, m, c, a, d, n : Int((d*x)**m*(a + b*acosh(c*x))**n, x))
    rubi.add(rule20)

    pattern21 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule21 = ReplacementRule(pattern21, lambda x, b, c, a, d, e : log(a + b*asinh(c*x))/(b*c*sqrt(d)))
    rubi.add(rule21)

    pattern22 = Pattern(Integral(S(1)/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)))
    rule22 = ReplacementRule(pattern22, lambda e2, e1, d1, x, b, c, a, d2 : log(a + b*acosh(c*x))/(b*c*sqrt(-d1*d2)))
    rubi.add(rule22)

    pattern23 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule23 = ReplacementRule(pattern23, lambda x, b, c, a, d, e, n : (a + b*asinh(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule23)

    pattern24 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule24 = ReplacementRule(pattern24, lambda e2, e1, d1, x, b, c, a, d2, n : (a + b*acosh(c*x))**(n + S(1))/(b*c*sqrt(-d1*d2)*(n + S(1))))
    rubi.add(rule24)

    pattern25 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule25 = ReplacementRule(pattern25, lambda x, b, c, a, d, e, n : sqrt(c**S(2)*x**S(2) + S(1))*Int((a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule25)

    pattern26 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda d2, d1: Not(NegativeQ(d2) & PositiveQ(d1))))
    rule26 = ReplacementRule(pattern26, lambda e2, e1, d1, x, b, c, a, d2, n : sqrt(c*x + S(-1))*sqrt(c*x + S(1))*Int((a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)))
    rubi.add(rule26)

    pattern27 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), )
    def With27(p, x, b, c, a, d, e):
        u = IntHide((d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asinh(c*x), u, x)
    rule27 = ReplacementRule(pattern27, lambda p, x, b, c, a, d, e : With27(p, x, b, c, a, d, e))
    rubi.add(rule27)

    pattern28 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), )
    def With28(p, x, b, c, a, d, e):
        u = IntHide((d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Dist(a + b*acosh(c*x), u, x)
    rule28 = ReplacementRule(pattern28, lambda p, x, b, c, a, d, e : With28(p, x, b, c, a, d, e))
    rubi.add(rule28)

    pattern29 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p: IntegerQ(p)))
    rule29 = ReplacementRule(pattern29, lambda p, x, b, c, a, d, e, n : -b*c*n*(-d)**p*Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(S(2)*p + S(1)) + S(2)*d*p*Int((a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule29)

    pattern30 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule30 = ReplacementRule(pattern30, lambda x, b, c, a, d, e, n : -b*c*n*sqrt(d + e*x**S(2))*Int(x*(a + b*asinh(c*x))**(n + S(-1)), x)/(S(2)*sqrt(c**S(2)*x**S(2) + S(1))) + x*(a + b*asinh(c*x))**n*sqrt(d + e*x**S(2))/S(2) + sqrt(d + e*x**S(2))*Int((a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x)/(S(2)*sqrt(c**S(2)*x**S(2) + S(1))))
    rubi.add(rule30)

    pattern31 = Pattern(Integral(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule31 = ReplacementRule(pattern31, lambda e2, e1, d1, x, b, c, a, d2, n : -b*c*n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int(x*(a + b*acosh(c*x))**(n + S(-1)), x)/(S(2)*sqrt(c*x + S(-1))*sqrt(c*x + S(1))) + x*(a + b*acosh(c*x))**n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/S(2) - sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int((a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(S(2)*sqrt(c*x + S(-1))*sqrt(c*x + S(1))))
    rubi.add(rule31)

    pattern32 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))))
    rule32 = ReplacementRule(pattern32, lambda p, x, b, c, a, d, e, n : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(S(2)*p + S(1)) + S(2)*d*p*Int((a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule32)

    pattern33 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule33 = ReplacementRule(pattern33, lambda p, e2, e1, d1, x, b, c, a, d2, n : -b*c*n*(-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x)/((S(2)*p + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))) + S(2)*d1*d2*p*Int((a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(-1))*(d2 + e2*x)**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p/(S(2)*p + S(1)))
    rubi.add(rule33)

    pattern34 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))))
    rule34 = ReplacementRule(pattern34, lambda p, e2, e1, d1, x, b, c, a, d2, n : -b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(S(2)*p + S(1)) + S(2)*d1*d2*p*Int((a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(-1))*(d2 + e2*x)**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p/(S(2)*p + S(1)))
    rubi.add(rule34)

    pattern35 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule35 = ReplacementRule(pattern35, lambda x, b, c, a, d, e, n : -b*c*n*sqrt(c**S(2)*x**S(2) + S(1))*Int(x*(a + b*asinh(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x)/(d*sqrt(d + e*x**S(2))) + x*(a + b*asinh(c*x))**n/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule35)

    pattern36 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/((d1_ + x_*WC('e1', S(1)))**(S(3)/2)*(d2_ + x_*WC('e2', S(1)))**(S(3)/2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule36 = ReplacementRule(pattern36, lambda e2, e1, d1, x, b, c, a, d2, n : b*c*n*sqrt(c*x + S(-1))*sqrt(c*x + S(1))*Int(x*(a + b*acosh(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x)/(d1*d2*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)) + x*(a + b*acosh(c*x))**n/(d1*d2*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)))
    rubi.add(rule36)

    pattern37 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: IntegerQ(p)))
    rule37 = ReplacementRule(pattern37, lambda p, x, b, c, a, d, e, n : -b*c*n*(-d)**p*Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(S(2)*p + S(2)) - x*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule37)

    pattern38 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule38 = ReplacementRule(pattern38, lambda p, x, b, c, a, d, e, n : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*(p + S(1))) - x*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule38)

    pattern39 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)))
    rule39 = ReplacementRule(pattern39, lambda p, e2, e1, d1, x, b, c, a, d2, n : -b*c*n*(-d1*d2)**(p + S(1)/2)*sqrt(c*x + S(-1))*sqrt(c*x + S(1))*Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x)/(S(2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*(p + S(1))) - x*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*d1*d2*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x)/(S(2)*d1*d2*(p + S(1))))
    rubi.add(rule39)

    pattern40 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule40 = ReplacementRule(pattern40, lambda p, e2, e1, d1, x, b, c, a, d2, n : -b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(S(2)*(p + S(1))) - x*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*d1*d2*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x)/(S(2)*d1*d2*(p + S(1))))
    rubi.add(rule40)

    pattern41 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule41 = ReplacementRule(pattern41, lambda x, b, c, a, d, e, n : Subst(Int((a + b*x)**n*sech(x), x), x, asinh(c*x))/(c*d))
    rubi.add(rule41)

    pattern42 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule42 = ReplacementRule(pattern42, lambda x, b, c, a, d, e, n : -Subst(Int((a + b*x)**n*csch(x), x), x, acosh(c*x))/(c*d))
    rubi.add(rule42)

    pattern43 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p: IntegerQ(p)))
    rule43 = ReplacementRule(pattern43, lambda p, x, b, c, a, d, e, n : -c*(-d)**p*(S(2)*p + S(1))*Int(x*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(b*(n + S(1))) + (-d)**p*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2)/(b*c*(n + S(1))))
    rubi.add(rule43)

    pattern44 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule44 = ReplacementRule(pattern44, lambda p, x, b, c, a, d, e, n : -c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(S(2)*p + S(1))*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*asinh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*(n + S(1))) + (a + b*asinh(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule44)

    pattern45 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule45 = ReplacementRule(pattern45, lambda p, e2, e1, d1, x, b, c, a, d2, n : -c*(-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*(S(2)*p + S(1))*Int(x*(a + b*acosh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x)/(b*(n + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))) + (a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**p*(d2 + e2*x)**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))))
    rubi.add(rule45)

    pattern46 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule46 = ReplacementRule(pattern46, lambda p, e2, e1, d1, x, b, c, a, d2, n : -c*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(S(2)*p + S(1))*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int(x*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(b*(n + S(1))) + (a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**p*(d2 + e2*x)**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))))
    rubi.add(rule46)

    pattern47 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)), CustomConstraint(lambda d, p: IntegerQ(p) | PositiveQ(d)))
    rule47 = ReplacementRule(pattern47, lambda p, x, b, c, a, d, e, n : d**p*Subst(Int((a + b*x)**n*Cosh(x)**(S(2)*p + S(1)), x), x, asinh(c*x))/c)
    rubi.add(rule47)

    pattern48 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule48 = ReplacementRule(pattern48, lambda p, x, b, c, a, d, e, n : (-d)**p*Subst(Int((a + b*x)**n*sinh(x)**(S(2)*p + S(1)), x), x, acosh(c*x))/c)
    rubi.add(rule48)

    pattern49 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda d2, d1: NegativeQ(d2) & PositiveQ(d1)))
    rule49 = ReplacementRule(pattern49, lambda p, e2, e1, d1, x, b, c, a, d2, n : (-d1*d2)**p*Subst(Int((a + b*x)**n*sinh(x)**(S(2)*p + S(1)), x), x, acosh(c*x))/c)
    rubi.add(rule49)

    pattern50 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)), CustomConstraint(lambda d, p: Not(IntegerQ(p) | PositiveQ(d))))
    rule50 = ReplacementRule(pattern50, lambda p, x, b, c, a, d, e, n : d**(p + S(-1)/2)*sqrt(d + e*x**S(2))*Int((a + b*asinh(c*x))**n*(c**S(2)*x**S(2) + S(1))**p, x)/sqrt(c**S(2)*x**S(2) + S(1)))
    rubi.add(rule50)

    pattern51 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)), CustomConstraint(lambda d2, d1: Not(NegativeQ(d2) & PositiveQ(d1))))
    rule51 = ReplacementRule(pattern51, lambda p, e2, e1, d1, x, b, c, a, d2, n : (-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int((a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p, x)/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))))
    rubi.add(rule51)

    pattern52 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: NonzeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With52(p, x, b, c, a, d, e):
        u = IntHide((d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asinh(c*x), u, x)
    rule52 = ReplacementRule(pattern52, lambda p, x, b, c, a, d, e : With52(p, x, b, c, a, d, e))
    rubi.add(rule52)

    pattern53 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With53(p, x, b, c, a, d, e):
        u = IntHide((d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Dist(a + b*acosh(c*x), u, x)
    rule53 = ReplacementRule(pattern53, lambda p, x, b, c, a, d, e : With53(p, x, b, c, a, d, e))
    rubi.add(rule53)

    pattern54 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: NonzeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p, n: PositiveIntegerQ(n) | Greater(p, S(0))))
    rule54 = ReplacementRule(pattern54, lambda p, x, b, c, a, d, e, n : Int(ExpandIntegrand((a + b*asinh(c*x))**n, (d + e*x**S(2))**p, x), x))
    rubi.add(rule54)

    pattern55 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p, n: PositiveIntegerQ(n) | Greater(p, S(0))))
    rule55 = ReplacementRule(pattern55, lambda p, x, b, c, a, d, e, n : Int(ExpandIntegrand((a + b*acosh(c*x))**n, (d + e*x**S(2))**p, x), x))
    rubi.add(rule55)

    pattern56 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule56 = ReplacementRule(pattern56, lambda p, x, b, c, a, d, e, n : Int((a + b*asinh(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule56)

    pattern57 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: IntegerQ(p)))
    rule57 = ReplacementRule(pattern57, lambda p, x, b, c, a, d, e, n : Int((a + b*acosh(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule57)

    pattern58 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule58 = ReplacementRule(pattern58, lambda p, e2, e1, d1, x, b, c, a, d2, n : Int((a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x))
    rubi.add(rule58)

    pattern59 = Pattern(Integral((d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda g, f, d, e: ZeroQ(d*g + e*f)), CustomConstraint(lambda g, c, f: ZeroQ(c**S(2)*f**S(2) + g**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule59 = ReplacementRule(pattern59, lambda p, x, b, g, c, f, a, d, e, n : (d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p))*Int((a + b*asinh(c*x))**n*(d*f + e*g*x**S(2))**p, x))
    rubi.add(rule59)

    pattern60 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule60 = ReplacementRule(pattern60, lambda p, x, b, c, a, d, e, n : (-d)**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p, x))
    rubi.add(rule60)

    pattern61 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule61 = ReplacementRule(pattern61, lambda x, b, c, a, d, e, n : Subst(Int((a + b*x)**n*tanh(x), x), x, asinh(c*x))/e)
    rubi.add(rule61)

    pattern62 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule62 = ReplacementRule(pattern62, lambda x, b, c, a, d, e, n : Subst(Int((a + b*x)**n*coth(x), x), x, acosh(c*x))/e)
    rubi.add(rule62)

    pattern63 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))), CustomConstraint(lambda p: IntegerQ(p)))
    rule63 = ReplacementRule(pattern63, lambda p, x, b, c, a, d, e, n : -b*n*(-d)**p*Int((a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) + (a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule63)

    pattern64 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule64 = ReplacementRule(pattern64, lambda p, x, b, c, a, d, e, n : -b*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) + (a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule64)

    pattern65 = Pattern(Integral(x_*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)))
    rule65 = ReplacementRule(pattern65, lambda p, e2, e1, d1, x, b, c, a, d2, n : -b*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) + (a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*e1*e2*(p + S(1))))
    rubi.add(rule65)

    pattern66 = Pattern(Integral(x_*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule66 = ReplacementRule(pattern66, lambda p, e2, e1, d1, x, b, c, a, d2, n : -b*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) + (a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*e1*e2*(p + S(1))))
    rubi.add(rule66)

    pattern67 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule67 = ReplacementRule(pattern67, lambda x, b, c, a, d, e, n : Subst(Int((a + b*x)**n/(Cosh(x)*sinh(x)), x), x, asinh(c*x))/d)
    rubi.add(rule67)

    pattern68 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule68 = ReplacementRule(pattern68, lambda x, b, c, a, d, e, n : -Subst(Int((a + b*x)**n/(Cosh(x)*sinh(x)), x), x, acosh(c*x))/d)
    rubi.add(rule68)

    pattern69 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda p: IntegerQ(p)))
    rule69 = ReplacementRule(pattern69, lambda p, x, b, m, c, f, a, d, e, n : b*c*n*(-d)**p*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule69)

    pattern70 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule70 = ReplacementRule(pattern70, lambda p, x, b, m, c, f, a, d, e, n : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + (f*x)**(m + S(1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule70)

    pattern71 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)))
    rule71 = ReplacementRule(pattern71, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x)/(f*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(d1*d2*f*(m + S(1))))
    rubi.add(rule71)

    pattern72 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule72 = ReplacementRule(pattern72, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(d1*d2*f*(m + S(1))))
    rubi.add(rule72)

    pattern73 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule73 = ReplacementRule(pattern73, lambda p, x, b, c, a, d, e : -b*c*d**p*Int((c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(S(2)*p) + d*Int((a + b*asinh(c*x))*(d + e*x**S(2))**(p + S(-1))/x, x) + (a + b*asinh(c*x))*(d + e*x**S(2))**p/(S(2)*p))
    rubi.add(rule73)

    pattern74 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule74 = ReplacementRule(pattern74, lambda p, x, b, c, a, d, e : -b*c*(-d)**p*Int((c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(S(2)*p) + d*Int((a + b*acosh(c*x))*(d + e*x**S(2))**(p + S(-1))/x, x) + (a + b*acosh(c*x))*(d + e*x**S(2))**p/(S(2)*p))
    rubi.add(rule74)

    pattern75 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2) + S(1)/2)))
    rule75 = ReplacementRule(pattern75, lambda p, x, b, m, c, f, a, d, e : -b*c*d**p*Int((f*x)**(m + S(1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*asinh(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*asinh(c*x))*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule75)

    pattern76 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2) + S(1)/2)))
    rule76 = ReplacementRule(pattern76, lambda p, x, b, m, c, f, a, d, e : -b*c*(-d)**p*Int((f*x)**(m + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*acosh(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule76)

    pattern77 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), )
    def With77(p, x, b, m, c, f, a, d, e):
        u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asinh(c*x), u, x)
    rule77 = ReplacementRule(pattern77, lambda p, x, b, m, c, f, a, d, e : With77(p, x, b, m, c, f, a, d, e))
    rubi.add(rule77)

    pattern78 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), )
    def With78(p, x, b, m, c, f, a, d, e):
        u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Dist(a + b*acosh(c*x), u, x)
    rule78 = ReplacementRule(pattern78, lambda p, x, b, m, c, f, a, d, e : With78(p, x, b, m, c, f, a, d, e))
    rubi.add(rule78)

    pattern79 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda d: PositiveQ(d)), )
    def With79(p, x, b, m, c, a, d, e):
        u = IntHide(x**m*(c**S(2)*x**S(2) + S(1))**p, x)
        return -b*c*d**p*Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Dist(d**p*(a + b*asinh(c*x)), u, x)
    rule79 = ReplacementRule(pattern79, lambda p, x, b, m, c, a, d, e : With79(p, x, b, m, c, a, d, e))
    rubi.add(rule79)

    pattern80 = Pattern(Integral(x_**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), )
    def With80(p, e2, e1, d1, x, b, m, c, a, d2):
        u = IntHide(x**m*(c*x + S(-1))**p*(c*x + S(1))**p, x)
        return -b*c*(-d1*d2)**p*Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Dist((-d1*d2)**p*(a + b*acosh(c*x)), u, x)
    rule80 = ReplacementRule(pattern80, lambda p, e2, e1, d1, x, b, m, c, a, d2 : With80(p, e2, e1, d1, x, b, m, c, a, d2))
    rubi.add(rule80)

    pattern81 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)), )
    def With81(p, x, b, m, c, a, d, e):
        u = IntHide(x**m*(c**S(2)*x**S(2) + S(1))**p, x)
        return -b*c*d**(p + S(-1)/2)*sqrt(d + e*x**S(2))*Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x)/sqrt(c**S(2)*x**S(2) + S(1)) + (a + b*asinh(c*x))*Int(x**m*(d + e*x**S(2))**p, x)
    rule81 = ReplacementRule(pattern81, lambda p, x, b, m, c, a, d, e : With81(p, x, b, m, c, a, d, e))
    rubi.add(rule81)

    pattern82 = Pattern(Integral(x_**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)), )
    def With82(p, e2, e1, d1, x, b, m, c, a, d2):
        u = IntHide(x**m*(c*x + S(-1))**p*(c*x + S(1))**p, x)
        return -b*c*(-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x)/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))) + (a + b*acosh(c*x))*Int(x**m*(d1 + e1*x)**p*(d2 + e2*x)**p, x)
    rule82 = ReplacementRule(pattern82, lambda p, e2, e1, d1, x, b, m, c, a, d2 : With82(p, e2, e1, d1, x, b, m, c, a, d2))
    rubi.add(rule82)

    pattern83 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p: IntegerQ(p)))
    rule83 = ReplacementRule(pattern83, lambda p, x, b, m, c, f, a, d, e, n : -b*c*n*(-d)**p*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule83)

    pattern84 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule84 = ReplacementRule(pattern84, lambda x, b, m, c, f, a, d, e, n : -b*c*n*sqrt(d + e*x**S(2))*Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1)), x)/(f*(m + S(1))*sqrt(c**S(2)*x**S(2) + S(1))) - c**S(2)*sqrt(d + e*x**S(2))*Int((f*x)**(m + S(2))*(a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x)/(f**S(2)*(m + S(1))*sqrt(c**S(2)*x**S(2) + S(1))) + (f*x)**(m + S(1))*(a + b*asinh(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(1))))
    rubi.add(rule84)

    pattern85 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule85 = ReplacementRule(pattern85, lambda e2, e1, d1, x, b, m, c, f, a, d2, n : -b*c*n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1)), x)/(f*(m + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))) - c**S(2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(f**S(2)*(m + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(f*(m + S(1))))
    rubi.add(rule85)

    pattern86 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule86 = ReplacementRule(pattern86, lambda p, x, b, m, c, f, a, d, e, n : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule86)

    pattern87 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule87 = ReplacementRule(pattern87, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : -b*c*n*(-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x)/(f*(m + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))) - S(2)*e1*e2*p*Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(-1))*(d2 + e2*x)**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p/(f*(m + S(1))))
    rubi.add(rule87)

    pattern88 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Not(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule88 = ReplacementRule(pattern88, lambda p, x, b, m, c, f, a, d, e, n : -b*c*n*(-d)**p*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(f*(m + S(2)*p + S(1))) + S(2)*d*p*Int((f*x)**m*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(m + S(2)*p + S(1)) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(2)*p + S(1))))
    rubi.add(rule88)

    pattern89 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Not(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule89 = ReplacementRule(pattern89, lambda x, b, m, c, f, a, d, e, n : -b*c*n*sqrt(d + e*x**S(2))*Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1)), x)/(f*(m + S(2))*sqrt(c**S(2)*x**S(2) + S(1))) + sqrt(d + e*x**S(2))*Int((f*x)**m*(a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x)/((m + S(2))*sqrt(c**S(2)*x**S(2) + S(1))) + (f*x)**(m + S(1))*(a + b*asinh(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(2))))
    rubi.add(rule89)

    pattern90 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Not(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule90 = ReplacementRule(pattern90, lambda e2, e1, d1, x, b, m, c, f, a, d2, n : -b*c*n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1)), x)/(f*(m + S(2))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))) - sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int((f*x)**m*(a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/((m + S(2))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(f*(m + S(2))))
    rubi.add(rule90)

    pattern91 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Not(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule91 = ReplacementRule(pattern91, lambda p, x, b, m, c, f, a, d, e, n : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(2)*p + S(1))) + S(2)*d*p*Int((f*x)**m*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(m + S(2)*p + S(1)) + (f*x)**(m + S(1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(2)*p + S(1))))
    rubi.add(rule91)

    pattern92 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Not(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule92 = ReplacementRule(pattern92, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : -b*c*n*(-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x)/(f*sqrt(c*x + S(-1))*sqrt(c*x + S(1))*(m + S(2)*p + S(1))) + S(2)*d1*d2*p*Int((f*x)**m*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(-1))*(d2 + e2*x)**(p + S(-1)), x)/(m + S(2)*p + S(1)) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p/(f*(m + S(2)*p + S(1))))
    rubi.add(rule92)

    pattern93 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p)))
    rule93 = ReplacementRule(pattern93, lambda p, x, b, m, c, f, a, d, e, n : b*c*n*(-d)**p*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + c**S(2)*(m + S(2)*p + S(3))*Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p, x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule93)

    pattern94 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule94 = ReplacementRule(pattern94, lambda p, x, b, m, c, f, a, d, e, n : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) - c**S(2)*(m + S(2)*p + S(3))*Int((f*x)**(m + S(2))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p, x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule94)

    pattern95 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)))
    rule95 = ReplacementRule(pattern95, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x)/(f*(m + S(1))) + c**S(2)*(m + S(2)*p + S(3))*Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(d1*d2*f*(m + S(1))))
    rubi.add(rule95)

    pattern96 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule96 = ReplacementRule(pattern96, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + c**S(2)*(m + S(2)*p + S(3))*Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(d1*d2*f*(m + S(1))))
    rubi.add(rule96)

    pattern97 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p: IntegerQ(p)))
    rule97 = ReplacementRule(pattern97, lambda p, x, b, m, c, f, a, d, e, n : -b*f*n*(-d)**p*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) - f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*e*(p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule97)

    pattern98 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule98 = ReplacementRule(pattern98, lambda p, x, b, m, c, f, a, d, e, n : -b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) - f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*e*(p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule98)

    pattern99 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)))
    rule99 = ReplacementRule(pattern99, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : -b*f*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) - f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x)/(S(2)*e1*e2*(p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*e1*e2*(p + S(1))))
    rubi.add(rule99)

    pattern100 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule100 = ReplacementRule(pattern100, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : -b*f*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) - f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x)/(S(2)*e1*e2*(p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*e1*e2*(p + S(1))))
    rubi.add(rule100)

    pattern101 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Not(RationalQ(m) & Greater(m, S(1)))), CustomConstraint(lambda p: IntegerQ(p)))
    rule101 = ReplacementRule(pattern101, lambda p, x, b, m, c, f, a, d, e, n : -b*c*n*(-d)**p*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(S(2)*f*(p + S(1))) + (m + S(2)*p + S(3))*Int((f*x)**m*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))) - (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*f*(p + S(1))))
    rubi.add(rule101)

    pattern102 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Not(RationalQ(m) & Greater(m, S(1)))), CustomConstraint(lambda n, p, m: IntegerQ(m) | IntegerQ(p) | Equal(n, S(1))))
    rule102 = ReplacementRule(pattern102, lambda p, x, b, m, c, f, a, d, e, n : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*f*(p + S(1))) + (m + S(2)*p + S(3))*Int((f*x)**m*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))) - (f*x)**(m + S(1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*f*(p + S(1))))
    rubi.add(rule102)

    pattern103 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Not(RationalQ(m) & Greater(m, S(1)))), CustomConstraint(lambda n, m: IntegerQ(m) | Equal(n, S(1))), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)))
    rule103 = ReplacementRule(pattern103, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : -b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x)/(S(2)*f*(p + S(1))) + (m + S(2)*p + S(3))*Int((f*x)**m*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x)/(S(2)*d1*d2*(p + S(1))) - (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*d1*d2*f*(p + S(1))))
    rubi.add(rule103)

    pattern104 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Not(RationalQ(m) & Greater(m, S(1)))), CustomConstraint(lambda n, p, m: IntegerQ(m) | IntegerQ(p) | Equal(n, S(1))))
    rule104 = ReplacementRule(pattern104, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : -b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(S(2)*f*(p + S(1))) + (m + S(2)*p + S(3))*Int((f*x)**m*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x)/(S(2)*d1*d2*(p + S(1))) - (f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*d1*d2*f*(p + S(1))))
    rubi.add(rule104)

    pattern105 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule105 = ReplacementRule(pattern105, lambda x, b, m, c, f, a, d, e, n : -b*f*n*sqrt(c**S(2)*x**S(2) + S(1))*Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(-1)), x)/(c*m*sqrt(d + e*x**S(2))) + f*(f*x)**(m + S(-1))*(a + b*asinh(c*x))**n*sqrt(d + e*x**S(2))/(e*m) - f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*asinh(c*x))**n/sqrt(d + e*x**S(2)), x)/(c**S(2)*m))
    rubi.add(rule105)

    pattern106 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule106 = ReplacementRule(pattern106, lambda e2, e1, d1, x, b, m, c, f, a, d2, n : b*f*n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1)), x)/(c*d1*d2*m*sqrt(c*x + S(-1))*sqrt(c*x + S(1))) + f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(e1*e2*m) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)), x)/(c**S(2)*m))
    rubi.add(rule106)

    pattern107 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule107 = ReplacementRule(pattern107, lambda x, b, m, c, a, d, e, n : c**(-m + S(-1))*Subst(Int((a + b*x)**n*sinh(x)**m, x), x, asinh(c*x))/sqrt(d))
    rubi.add(rule107)

    pattern108 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda m: IntegerQ(m)))
    rule108 = ReplacementRule(pattern108, lambda e2, e1, d1, x, b, m, c, a, d2, n : c**(-m + S(-1))*Subst(Int((a + b*x)**n*Cosh(x)**m, x), x, acosh(c*x))/sqrt(-d1*d2))
    rubi.add(rule108)

    pattern109 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule109 = ReplacementRule(pattern109, lambda x, b, m, c, f, a, d, e : -b*c*(f*x)**(m + S(2))*HypergeometricPFQ(List(S(1), m/S(2) + S(1), m/S(2) + S(1)), List(m/S(2) + S(3)/2, m/S(2) + S(2)), -c**S(2)*x**S(2))/(sqrt(d)*f**S(2)*(m + S(1))*(m + S(2))) + (f*x)**(m + S(1))*(a + b*asinh(c*x))*Hypergeometric2F1(S(1)/2, m/S(2) + S(1)/2, m/S(2) + S(3)/2, -c**S(2)*x**S(2))/(sqrt(d)*f*(m + S(1))))
    rubi.add(rule109)

    pattern110 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule110 = ReplacementRule(pattern110, lambda e2, e1, d1, x, b, m, c, f, a, d2 : b*c*(f*x)**(m + S(2))*HypergeometricPFQ(List(S(1), m/S(2) + S(1), m/S(2) + S(1)), List(m/S(2) + S(3)/2, m/S(2) + S(2)), c**S(2)*x**S(2))/(f**S(2)*sqrt(-d1*d2)*(m + S(1))*(m + S(2))) + (f*x)**(m + S(1))*(a + b*acosh(c*x))*sqrt(-c**S(2)*x**S(2) + S(1))*Hypergeometric2F1(S(1)/2, m/S(2) + S(1)/2, m/S(2) + S(3)/2, c**S(2)*x**S(2))/(f*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*(m + S(1))))
    rubi.add(rule110)

    pattern111 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda d: Not(PositiveQ(d))), CustomConstraint(lambda n, m: IntegerQ(m) | Equal(n, S(1))))
    rule111 = ReplacementRule(pattern111, lambda x, b, m, c, f, a, d, e, n : sqrt(c**S(2)*x**S(2) + S(1))*Int((f*x)**m*(a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule111)

    pattern112 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda d2, d1: Not(NegativeQ(d2) & PositiveQ(d1))), CustomConstraint(lambda n, m: IntegerQ(m) | Equal(n, S(1))))
    rule112 = ReplacementRule(pattern112, lambda e2, e1, d1, x, b, m, c, f, a, d2, n : sqrt(c*x + S(-1))*sqrt(c*x + S(1))*Int((f*x)**m*(a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)))
    rubi.add(rule112)

    pattern113 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(1))), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m)))
    rule113 = ReplacementRule(pattern113, lambda p, x, b, m, c, f, a, d, e, n : -b*f*n*(-d)**p*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(c*(m + S(2)*p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(e*(m + S(2)*p + S(1))) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p, x)/(c**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule113)

    pattern114 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule114 = ReplacementRule(pattern114, lambda p, x, b, m, c, f, a, d, e, n : -b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(c*(m + S(2)*p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(e*(m + S(2)*p + S(1))) - f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p, x)/(c**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule114)

    pattern115 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(1))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)))
    rule115 = ReplacementRule(pattern115, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : -b*f*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x)/(c*(m + S(2)*p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(e1*e2*(m + S(2)*p + S(1))) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x)/(c**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule115)

    pattern116 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule116 = ReplacementRule(pattern116, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : -b*f*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x)/(c*(m + S(2)*p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(e1*e2*(m + S(2)*p + S(1))) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x)/(c**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule116)

    pattern117 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(1))), CustomConstraint(lambda p: IntegerQ(p)))
    rule117 = ReplacementRule(pattern117, lambda p, x, m, b, c, f, a, d, e, n : f*m*(-d)**p*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) + (f*x)**m*(a + b*acosh(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))))
    rubi.add(rule117)

    pattern118 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(1))))
    rule118 = ReplacementRule(pattern118, lambda p, x, m, b, c, f, a, d, e, n : -d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) + (f*x)**m*(a + b*asinh(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule118)

    pattern119 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(1))), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule119 = ReplacementRule(pattern119, lambda p, e2, e1, d1, x, m, b, c, f, a, d2, n : f*m*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) + (f*x)**m*(a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**p*(d2 + e2*x)**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))))
    rubi.add(rule119)

    pattern120 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(1))))
    rule120 = ReplacementRule(pattern120, lambda p, e2, e1, d1, x, m, b, c, f, a, d2, n : f*m*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) + (f*x)**m*(a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**p*(d2 + e2*x)**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))))
    rubi.add(rule120)

    pattern121 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda d: PositiveQ(d)))
    rule121 = ReplacementRule(pattern121, lambda x, m, b, c, f, a, d, e, n : -f*m*Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(1)), x)/(b*c*sqrt(d)*(n + S(1))) + (f*x)**m*(a + b*asinh(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule121)

    pattern122 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)))
    rule122 = ReplacementRule(pattern122, lambda e2, e1, d1, x, m, b, c, f, a, d2, n : -f*m*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1)), x)/(b*c*sqrt(-d1*d2)*(n + S(1))) + (f*x)**m*(a + b*acosh(c*x))**(n + S(1))/(b*c*sqrt(-d1*d2)*(n + S(1))))
    rubi.add(rule122)

    pattern123 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule123 = ReplacementRule(pattern123, lambda x, m, b, c, f, a, d, e, n : sqrt(c**S(2)*x**S(2) + S(1))*Int((f*x)**m*(a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule123)

    pattern124 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda d2, d1: Not(NegativeQ(d2) & PositiveQ(d1))))
    rule124 = ReplacementRule(pattern124, lambda e2, e1, d1, x, m, b, c, f, a, d2, n : sqrt(c*x + S(-1))*sqrt(c*x + S(1))*Int((f*x)**m*(a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)))
    rubi.add(rule124)

    pattern125 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(-3))), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule125 = ReplacementRule(pattern125, lambda p, x, m, b, c, f, a, d, e, n : -c*(-d)**p*(m + S(2)*p + S(1))*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(b*f*(n + S(1))) + f*m*(-d)**p*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) + (f*x)**m*(a + b*acosh(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))))
    rubi.add(rule125)

    pattern126 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(-3))), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)))
    rule126 = ReplacementRule(pattern126, lambda p, x, m, b, c, f, a, d, e, n : -c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*(m + S(2)*p + S(1))*Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*f*(n + S(1))) - d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) + (f*x)**m*(a + b*asinh(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule126)

    pattern127 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(-3))), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)))
    rule127 = ReplacementRule(pattern127, lambda p, e2, e1, d1, x, m, b, c, f, a, d2, n : -c*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*(m + S(2)*p + S(1))*Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x)/(b*f*(n + S(1))) + f*m*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) + (f*x)**m*(a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**p*(d2 + e2*x)**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))))
    rubi.add(rule127)

    pattern128 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, p: IntegerQ(p) | PositiveQ(d)))
    rule128 = ReplacementRule(pattern128, lambda p, x, b, m, c, a, d, e, n : c**(-m + S(-1))*d**p*Subst(Int((a + b*x)**n*Cosh(x)**(S(2)*p + S(1))*sinh(x)**m, x), x, asinh(c*x)))
    rubi.add(rule128)

    pattern129 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule129 = ReplacementRule(pattern129, lambda p, x, b, m, c, a, d, e, n : c**(-m + S(-1))*(-d)**p*Subst(Int((a + b*x)**n*Cosh(x)**m*sinh(x)**(S(2)*p + S(1)), x), x, acosh(c*x)))
    rubi.add(rule129)

    pattern130 = Pattern(Integral(x_**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d2, d1: NegativeQ(d2) & PositiveQ(d1)))
    rule130 = ReplacementRule(pattern130, lambda p, e2, e1, d1, x, b, m, c, a, d2, n : c**(-m + S(-1))*(-d1*d2)**p*Subst(Int((a + b*x)**n*Cosh(x)**m*sinh(x)**(S(2)*p + S(1)), x), x, acosh(c*x)))
    rubi.add(rule130)

    pattern131 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d, p: Not(IntegerQ(p) | PositiveQ(d))))
    rule131 = ReplacementRule(pattern131, lambda p, x, b, m, c, a, d, e, n : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x**m*(a + b*asinh(c*x))**n*(c**S(2)*x**S(2) + S(1))**p, x))
    rubi.add(rule131)

    pattern132 = Pattern(Integral(x_**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda d2, p, d1: Not(IntegerQ(p) | (NegativeQ(d2) & PositiveQ(d1)))))
    rule132 = ReplacementRule(pattern132, lambda p, e2, e1, d1, x, b, m, c, a, d2, n : (-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int(x**m*(a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p, x))
    rubi.add(rule132)

    pattern133 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda m: Not(PositiveIntegerQ(m/S(2) + S(1)/2))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Less(S(-3), m, S(0))))
    rule133 = ReplacementRule(pattern133, lambda p, x, b, m, c, f, a, d, e, n : Int(ExpandIntegrand((a + b*asinh(c*x))**n/sqrt(d + e*x**S(2)), (f*x)**m*(d + e*x**S(2))**(p + S(1)/2), x), x))
    rubi.add(rule133)

    pattern134 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda m: Not(PositiveIntegerQ(m/S(2) + S(1)/2))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Less(S(-3), m, S(0))))
    rule134 = ReplacementRule(pattern134, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : Int(ExpandIntegrand((a + b*acosh(c*x))**n/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)), (f*x)**m*(d1 + e1*x)**(p + S(1)/2)*(d2 + e2*x)**(p + S(1)/2), x), x))
    rubi.add(rule134)

    pattern135 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(3))))
    rule135 = ReplacementRule(pattern135, lambda x, m, b, c, f, a, d, e : -b*c*Int((f*x)**(m + S(1))*(d*(m + S(3)) + e*x**S(2)*(m + S(1)))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(f*(m + S(1))*(m + S(3))) + d*(f*x)**(m + S(1))*(a + b*acosh(c*x))/(f*(m + S(1))) + e*(f*x)**(m + S(3))*(a + b*acosh(c*x))/(f**S(3)*(m + S(3))))
    rubi.add(rule135)

    pattern136 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, d, e: NonzeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule136 = ReplacementRule(pattern136, lambda p, x, b, c, a, d, e : -b*c*Int((d + e*x**S(2))**(p + S(1))/sqrt(c**S(2)*x**S(2) + S(1)), x)/(S(2)*e*(p + S(1))) + (a + b*asinh(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule136)

    pattern137 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule137 = ReplacementRule(pattern137, lambda p, x, b, c, a, d, e : -b*c*Int((d + e*x**S(2))**(p + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(S(2)*e*(p + S(1))) + (a + b*acosh(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule137)

    pattern138 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d, e: NonzeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (LessEqual(m + p, S(0)) & PositiveIntegerQ(m/S(2) + S(-1)/2))), )
    def With138(p, x, m, b, c, f, a, d, e):
        u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asinh(c*x), u, x)
    rule138 = ReplacementRule(pattern138, lambda p, x, m, b, c, f, a, d, e : With138(p, x, m, b, c, f, a, d, e))
    rubi.add(rule138)

    pattern139 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (LessEqual(m + p, S(0)) & PositiveIntegerQ(m/S(2) + S(-1)/2))), )
    def With139(p, x, m, b, c, f, a, d, e):
        u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Dist(a + b*acosh(c*x), u, x)
    rule139 = ReplacementRule(pattern139, lambda p, x, m, b, c, f, a, d, e : With139(p, x, m, b, c, f, a, d, e))
    rubi.add(rule139)

    pattern140 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, d, e: NonzeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m)))
    rule140 = ReplacementRule(pattern140, lambda p, x, b, m, c, f, a, d, e, n : Int(ExpandIntegrand((a + b*asinh(c*x))**n, (f*x)**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule140)

    pattern141 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda c, e, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m)))
    rule141 = ReplacementRule(pattern141, lambda p, x, b, m, c, f, a, d, e, n : Int(ExpandIntegrand((a + b*acosh(c*x))**n, (f*x)**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule141)

    pattern142 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule142 = ReplacementRule(pattern142, lambda p, x, b, m, c, f, a, d, e, n : Int((f*x)**m*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule142)

    pattern143 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: IntegerQ(p)))
    rule143 = ReplacementRule(pattern143, lambda p, x, b, m, c, f, a, d, e, n : Int((f*x)**m*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule143)

    pattern144 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule144 = ReplacementRule(pattern144, lambda p, e2, e1, d1, x, b, m, c, f, a, d2, n : Int((f*x)**m*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x))
    rubi.add(rule144)

    pattern145 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda g, f, d, e: ZeroQ(d*g + e*f)), CustomConstraint(lambda g, c, f: ZeroQ(c**S(2)*f**S(2) + g**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule145 = ReplacementRule(pattern145, lambda p, h, x, b, m, g, c, f, a, d, e, n : (d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p))*Int((h*x)**m*(a + b*asinh(c*x))**n*(d*f + e*g*x**S(2))**p, x))
    rubi.add(rule145)

    pattern146 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule146 = ReplacementRule(pattern146, lambda p, x, b, m, c, f, a, d, e, n : (-d)**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((f*x)**m*(a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p, x))
    rubi.add(rule146)

    pattern147 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule147 = ReplacementRule(pattern147, lambda x, b, c, a, d, e, n : Subst(Int((a + b*x)**n*Cosh(x)/(c*d + e*sinh(x)), x), x, asinh(c*x)))
    rubi.add(rule147)

    pattern148 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule148 = ReplacementRule(pattern148, lambda x, b, c, a, d, e, n : Subst(Int((a + b*x)**n*sinh(x)/(c*d + e*Cosh(x)), x), x, acosh(c*x)))
    rubi.add(rule148)

    pattern149 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule149 = ReplacementRule(pattern149, lambda x, b, m, c, a, d, e, n : -b*c*n*Int((a + b*asinh(c*x))**(n + S(-1))*(d + e*x)**(m + S(1))/sqrt(c**S(2)*x**S(2) + S(1)), x)/(e*(m + S(1))) + (a + b*asinh(c*x))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))))
    rubi.add(rule149)

    pattern150 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule150 = ReplacementRule(pattern150, lambda x, b, m, c, a, d, e, n : -b*c*n*Int((a + b*acosh(c*x))**(n + S(-1))*(d + e*x)**(m + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x)/(e*(m + S(1))) + (a + b*acosh(c*x))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))))
    rubi.add(rule150)

    pattern151 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule151 = ReplacementRule(pattern151, lambda x, b, m, c, a, d, e, n : Int(ExpandIntegrand((a + b*asinh(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule151)

    pattern152 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule152 = ReplacementRule(pattern152, lambda x, b, m, c, a, d, e, n : Int(ExpandIntegrand((a + b*acosh(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule152)

    pattern153 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule153 = ReplacementRule(pattern153, lambda x, b, m, c, a, d, e, n : c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*d + e*sinh(x))**m*Cosh(x), x), x, asinh(c*x)))
    rubi.add(rule153)

    pattern154 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule154 = ReplacementRule(pattern154, lambda x, b, m, c, a, d, e, n : c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*d + e*Cosh(x))**m*sinh(x), x), x, acosh(c*x)))
    rubi.add(rule154)

    pattern155 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)), )
    def With155(Px, x, b, c, a):
        u = IntHide(Px, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asinh(c*x), u, x)
    rule155 = ReplacementRule(pattern155, lambda Px, x, b, c, a : With155(Px, x, b, c, a))
    rubi.add(rule155)

    pattern156 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)), )
    def With156(Px, x, b, c, a):
        u = IntHide(Px, x)
        return -b*c*sqrt(-c**S(2)*x**S(2) + S(1))*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x)/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))) + Dist(a + b*acosh(c*x), u, x)
    rule156 = ReplacementRule(pattern156, lambda Px, x, b, c, a : With156(Px, x, b, c, a))
    rubi.add(rule156)

    pattern157 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)))
    rule157 = ReplacementRule(pattern157, lambda Px, x, b, c, a, n : Int(ExpandIntegrand(Px*(a + b*asinh(c*x))**n, x), x))
    rubi.add(rule157)

    pattern158 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)))
    rule158 = ReplacementRule(pattern158, lambda Px, x, b, c, a, n : Int(ExpandIntegrand(Px*(a + b*acosh(c*x))**n, x), x))
    rubi.add(rule158)

    pattern159 = Pattern(Integral(Px_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)), )
    def With159(Px, x, b, m, c, a, d, e):
        u = IntHide(Px*(d + e*x)**m, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asinh(c*x), u, x)
    rule159 = ReplacementRule(pattern159, lambda Px, x, b, m, c, a, d, e : With159(Px, x, b, m, c, a, d, e))
    rubi.add(rule159)

    pattern160 = Pattern(Integral(Px_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)), )
    def With160(Px, x, b, m, c, a, d, e):
        u = IntHide(Px*(d + e*x)**m, x)
        return -b*c*sqrt(-c**S(2)*x**S(2) + S(1))*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x)/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))) + Dist(a + b*acosh(c*x), u, x)
    rule160 = ReplacementRule(pattern160, lambda Px, x, b, m, c, a, d, e : With160(Px, x, b, m, c, a, d, e))
    rubi.add(rule160)

    pattern161 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda p, m: Less(m + p + S(1), S(0))), )
    def With161(p, x, b, m, g, c, f, a, d, e, n):
        u = IntHide((d + e*x)**m*(f + g*x)**p, x)
        return -b*c*n*Int(SimplifyIntegrand(u*(a + b*asinh(c*x))**(n + S(-1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Dist((a + b*asinh(c*x))**n, u, x)
    rule161 = ReplacementRule(pattern161, lambda p, x, b, m, g, c, f, a, d, e, n : With161(p, x, b, m, g, c, f, a, d, e, n))
    rubi.add(rule161)

    pattern162 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda p, m: Less(m + p + S(1), S(0))), )
    def With162(p, x, b, m, g, c, f, a, d, e, n):
        u = IntHide((d + e*x)**m*(f + g*x)**p, x)
        return -b*c*n*Int(SimplifyIntegrand(u*(a + b*acosh(c*x))**(n + S(-1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Dist((a + b*acosh(c*x))**n, u, x)
    rule162 = ReplacementRule(pattern162, lambda p, x, b, m, g, c, f, a, d, e, n : With162(p, x, b, m, g, c, f, a, d, e, n))
    rubi.add(rule162)

    pattern163 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_*(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))/(d_ + x_*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)), CustomConstraint(lambda g, d, e, h: ZeroQ(-S(2)*d*h + e*g)), )
    def With163(p, h, x, b, g, c, f, a, d, e, n):
        u = IntHide((f + g*x + h*x**S(2))**p/(d + e*x)**S(2), x)
        return -b*c*n*Int(SimplifyIntegrand(u*(a + b*asinh(c*x))**(n + S(-1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Dist((a + b*asinh(c*x))**n, u, x)
    rule163 = ReplacementRule(pattern163, lambda p, h, x, b, g, c, f, a, d, e, n : With163(p, h, x, b, g, c, f, a, d, e, n))
    rubi.add(rule163)

    pattern164 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_*(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))/(d_ + x_*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)), CustomConstraint(lambda g, d, e, h: ZeroQ(-S(2)*d*h + e*g)), )
    def With164(p, h, x, b, g, c, f, a, d, e, n):
        u = IntHide((f + g*x + h*x**S(2))**p/(d + e*x)**S(2), x)
        return -b*c*n*Int(SimplifyIntegrand(u*(a + b*acosh(c*x))**(n + S(-1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Dist((a + b*acosh(c*x))**n, u, x)
    rule164 = ReplacementRule(pattern164, lambda p, h, x, b, g, c, f, a, d, e, n : With164(p, h, x, b, g, c, f, a, d, e, n))
    rubi.add(rule164)

    pattern165 = Pattern(Integral(Px_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule165 = ReplacementRule(pattern165, lambda Px, x, b, m, c, a, d, e, n : Int(ExpandIntegrand(Px*(a + b*asinh(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule165)

    pattern166 = Pattern(Integral(Px_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule166 = ReplacementRule(pattern166, lambda Px, x, b, m, c, a, d, e, n : Int(ExpandIntegrand(Px*(a + b*acosh(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule166)

    pattern167 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p, m: Greater(m, S(3)) | Less(m, -S(2)*p + S(-1))), )
    def With167(p, x, b, m, g, c, f, a, d, e):
        u = IntHide((d + e*x**S(2))**p*(f + g*x)**m, x)
        return -b*c*Int(Dist(1/sqrt(c**S(2)*x**S(2) + S(1)), u, x), x) + Dist(a + b*asinh(c*x), u, x)
    rule167 = ReplacementRule(pattern167, lambda p, x, b, m, g, c, f, a, d, e : With167(p, x, b, m, g, c, f, a, d, e))
    rubi.add(rule167)

    pattern168 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p, m: Greater(m, S(3)) | Less(m, -S(2)*p + S(-1))), )
    def With168(p, e2, e1, d1, x, b, m, g, c, f, a, d2):
        u = IntHide((d1 + e1*x)**p*(d2 + e2*x)**p*(f + g*x)**m, x)
        return -b*c*Int(Dist(S(1)/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), u, x), x) + Dist(a + b*acosh(c*x), u, x)
    rule168 = ReplacementRule(pattern168, lambda p, e2, e1, d1, x, b, m, g, c, f, a, d2 : With168(p, e2, e1, d1, x, b, m, g, c, f, a, d2))
    rubi.add(rule168)

    pattern169 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p, n, m: Equal(m, S(1)) | Greater(p, S(0)) | (Equal(m, S(2)) & Less(p, S(-2))) | (Equal(n, S(1)) & Greater(p, S(-1)))))
    rule169 = ReplacementRule(pattern169, lambda p, x, b, m, g, c, f, a, d, e, n : Int(ExpandIntegrand((a + b*asinh(c*x))**n*(d + e*x**S(2))**p, (f + g*x)**m, x), x))
    rubi.add(rule169)

    pattern170 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p, n, m: Equal(m, S(1)) | Greater(p, S(0)) | (Equal(m, S(2)) & Less(p, S(-2))) | (Equal(n, S(1)) & Greater(p, S(-1)))))
    rule170 = ReplacementRule(pattern170, lambda p, e2, e1, d1, x, b, m, g, c, f, a, d2, n : Int(ExpandIntegrand((a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, (f + g*x)**m, x), x))
    rubi.add(rule170)

    pattern171 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(x_*WC('g', S(1)) + WC('f', S(0)))**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule171 = ReplacementRule(pattern171, lambda x, b, m, g, c, f, a, d, e, n : (a + b*asinh(c*x))**(n + S(1))*(d + e*x**S(2))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))) - Int((a + b*asinh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1))*(d*g*m + S(2)*e*f*x + e*g*x**S(2)*(m + S(2))), x)/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule171)

    pattern172 = Pattern(Integral(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule172 = ReplacementRule(pattern172, lambda e2, e1, d1, x, b, m, g, c, f, a, d2, n : (a + b*acosh(c*x))**(n + S(1))*(f + g*x)**m*(d1*d2 + e1*e2*x**S(2))/(b*c*sqrt(-d1*d2)*(n + S(1))) - Int((a + b*acosh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1))*(d1*d2*g*m + S(2)*e1*e2*f*x + e1*e2*g*x**S(2)*(m + S(2))), x)/(b*c*sqrt(-d1*d2)*(n + S(1))))
    rubi.add(rule172)

    pattern173 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule173 = ReplacementRule(pattern173, lambda p, x, b, m, g, c, f, a, d, e, n : Int(ExpandIntegrand((a + b*asinh(c*x))**n*sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(-1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule173)

    pattern174 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule174 = ReplacementRule(pattern174, lambda p, e2, e1, d1, x, b, m, g, c, f, a, d2, n : Int(ExpandIntegrand((a + b*acosh(c*x))**n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x), (d1 + e1*x)**(p + S(-1)/2)*(d2 + e2*x)**(p + S(-1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule174)

    pattern175 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule175 = ReplacementRule(pattern175, lambda p, x, b, m, g, c, f, a, d, e, n : (a + b*asinh(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))) - Int(ExpandIntegrand((a + b*asinh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), (d + e*x**S(2))**(p + S(-1)/2)*(d*g*m + e*f*x*(S(2)*p + S(1)) + e*g*x**S(2)*(m + S(2)*p + S(1))), x), x)/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule175)

    pattern176 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(-1)/2)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule176 = ReplacementRule(pattern176, lambda p, e2, e1, d1, x, b, m, g, c, f, a, d2, n : (a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**(p + S(1)/2)*(d2 + e2*x)**(p + S(1)/2)*(f + g*x)**m/(b*c*sqrt(-d1*d2)*(n + S(1))) - Int(ExpandIntegrand((a + b*acosh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), (d1 + e1*x)**(p + S(-1)/2)*(d2 + e2*x)**(p + S(-1)/2)*(d1*d2*g*m + e1*e2*f*x*(S(2)*p + S(1)) + e1*e2*g*x**S(2)*(m + S(2)*p + S(1))), x), x)/(b*c*sqrt(-d1*d2)*(n + S(1))))
    rubi.add(rule176)

    pattern177 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule177 = ReplacementRule(pattern177, lambda x, b, m, g, c, f, a, d, e, n : -g*m*Int((a + b*asinh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), x)/(b*c*sqrt(d)*(n + S(1))) + (a + b*asinh(c*x))**(n + S(1))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule177)

    pattern178 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule178 = ReplacementRule(pattern178, lambda e2, e1, d1, x, b, m, g, c, f, a, d2, n : -g*m*Int((a + b*acosh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), x)/(b*c*sqrt(-d1*d2)*(n + S(1))) + (a + b*acosh(c*x))**(n + S(1))*(f + g*x)**m/(b*c*sqrt(-d1*d2)*(n + S(1))))
    rubi.add(rule178)

    pattern179 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n, m: PositiveIntegerQ(n) | Greater(m, S(0))))
    rule179 = ReplacementRule(pattern179, lambda x, b, m, g, c, f, a, d, e, n : c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*f + g*sinh(x))**m, x), x, asinh(c*x))/sqrt(d))
    rubi.add(rule179)

    pattern180 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda n, m: PositiveIntegerQ(n) | Greater(m, S(0))))
    rule180 = ReplacementRule(pattern180, lambda e2, e1, d1, x, b, m, g, c, f, a, d2, n : c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*f + g*Cosh(x))**m, x), x, acosh(c*x))/sqrt(-d1*d2))
    rubi.add(rule180)

    pattern181 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule181 = ReplacementRule(pattern181, lambda p, x, b, m, g, c, f, a, d, e, n : Int(ExpandIntegrand((a + b*asinh(c*x))**n/sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule181)

    pattern182 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule182 = ReplacementRule(pattern182, lambda p, e2, e1, d1, x, b, m, g, c, f, a, d2, n : Int(ExpandIntegrand((a + b*acosh(c*x))**n/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)), (d1 + e1*x)**(p + S(1)/2)*(d2 + e2*x)**(p + S(1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule182)

    pattern183 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule183 = ReplacementRule(pattern183, lambda p, x, b, m, g, c, f, a, d, e, n : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*asinh(c*x))**n*(f + g*x)**m*(c**S(2)*x**S(2) + S(1))**p, x))
    rubi.add(rule183)

    pattern184 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule184 = ReplacementRule(pattern184, lambda p, x, b, m, g, c, f, a, d, e, n : (-d)**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((a + b*acosh(c*x))**n*(f + g*x)**m*(c*x + S(-1))**p*(c*x + S(1))**p, x))
    rubi.add(rule184)

    pattern185 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d2, d1: Not(NegativeQ(d2) & PositiveQ(d1))))
    rule185 = ReplacementRule(pattern185, lambda p, e2, e1, d1, x, b, m, g, c, f, a, d2, n : (-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*acosh(c*x))**n*(f + g*x)**m*(c*x + S(-1))**p*(c*x + S(1))**p, x))
    rubi.add(rule185)

    pattern186 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule186 = ReplacementRule(pattern186, lambda h, x, m, b, g, f, c, a, d, e, n : -g*m*Int((a + b*asinh(c*x))**(n + S(1))/(f + g*x), x)/(b*c*sqrt(d)*(n + S(1))) + (a + b*asinh(c*x))**(n + S(1))*log(h*(f + g*x)**m)/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule186)

    pattern187 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1)))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda d1: PositiveQ(d1)), CustomConstraint(lambda d2: NegativeQ(d2)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule187 = ReplacementRule(pattern187, lambda e2, h, e1, x, d1, m, b, g, f, c, a, d2, n : -g*m*Int((a + b*acosh(c*x))**(n + S(1))/(f + g*x), x)/(b*c*sqrt(-d1*d2)*(n + S(1))) + (a + b*acosh(c*x))**(n + S(1))*log(h*(f + g*x)**m)/(b*c*sqrt(-d1*d2)*(n + S(1))))
    rubi.add(rule187)

    pattern188 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule188 = ReplacementRule(pattern188, lambda p, h, x, m, b, g, f, c, a, d, e, n : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*asinh(c*x))**n*(c**S(2)*x**S(2) + S(1))**p*log(h*(f + g*x)**m), x))
    rubi.add(rule188)

    pattern189 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule189 = ReplacementRule(pattern189, lambda p, h, x, m, b, g, f, c, a, d, e, n : (-d)**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p*log(h*(f + g*x)**m), x))
    rubi.add(rule189)

    pattern190 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d2, d1: Not(NegativeQ(d2) & PositiveQ(d1))))
    rule190 = ReplacementRule(pattern190, lambda p, e2, h, e1, x, d1, m, b, g, f, c, a, d2, n : (-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*Int((a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p*log(h*(f + g*x)**m), x))
    rubi.add(rule190)

    pattern191 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m + S(1)/2)), )
    def With191(x, b, m, g, c, f, a, d, e):
        u = IntHide((d + e*x)**m*(f + g*x)**m, x)
        return -b*c*Int(Dist(1/sqrt(c**S(2)*x**S(2) + S(1)), u, x), x) + Dist(a + b*asinh(c*x), u, x)
    rule191 = ReplacementRule(pattern191, lambda x, b, m, g, c, f, a, d, e : With191(x, b, m, g, c, f, a, d, e))
    rubi.add(rule191)

    pattern192 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m + S(1)/2)), )
    def With192(x, b, m, g, c, f, a, d, e):
        u = IntHide((d + e*x)**m*(f + g*x)**m, x)
        return -b*c*Int(Dist(S(1)/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), u, x), x) + Dist(a + b*acosh(c*x), u, x)
    rule192 = ReplacementRule(pattern192, lambda x, b, m, g, c, f, a, d, e : With192(x, b, m, g, c, f, a, d, e))
    rubi.add(rule192)

    pattern193 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule193 = ReplacementRule(pattern193, lambda x, b, m, g, c, f, a, d, e, n : Int(ExpandIntegrand((a + b*asinh(c*x))**n, (d + e*x)**m*(f + g*x)**m, x), x))
    rubi.add(rule193)

    pattern194 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule194 = ReplacementRule(pattern194, lambda x, b, m, g, c, f, a, d, e, n : Int(ExpandIntegrand((a + b*acosh(c*x))**n, (d + e*x)**m*(f + g*x)**m, x), x))
    rubi.add(rule194)

    pattern195 = Pattern(Integral(u_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda a, v, x, b, c: InverseFunctionFreeQ(v, x)))
    def With195(x, b, c, a, u):
        v = IntHide(u, x)
        return -b*c*Int(SimplifyIntegrand(v/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asinh(c*x), v, x)
    rule195 = ReplacementRule(pattern195, lambda x, b, c, a, u : With195(x, b, c, a, u))
    rubi.add(rule195)

    pattern196 = Pattern(Integral(u_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda a, v, x, b, c: InverseFunctionFreeQ(v, x)))
    def With196(x, b, c, a, u):
        v = IntHide(u, x)
        return -b*c*sqrt(-c**S(2)*x**S(2) + S(1))*Int(SimplifyIntegrand(v/sqrt(-c**S(2)*x**S(2) + S(1)), x), x)/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))) + Dist(a + b*acosh(c*x), v, x)
    rule196 = ReplacementRule(pattern196, lambda x, b, c, a, u : With196(x, b, c, a, u))
    rubi.add(rule196)

    pattern197 = Pattern(Integral(Px_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda x, u: SumQ(u)))
    def With197(Px, p, x, b, c, a, d, e, n):
        u = ExpandIntegrand(Px*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p, x)
        return Int(u, x)
    rule197 = ReplacementRule(pattern197, lambda Px, p, x, b, c, a, d, e, n : With197(Px, p, x, b, c, a, d, e, n))
    rubi.add(rule197)

    pattern198 = Pattern(Integral(Px_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda x, u: SumQ(u)))
    def With198(Px, p, e2, e1, d1, x, b, c, a, d2, n):
        u = ExpandIntegrand(Px*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x)
        return Int(u, x)
    rule198 = ReplacementRule(pattern198, lambda Px, p, e2, e1, d1, x, b, c, a, d2, n : With198(Px, p, e2, e1, d1, x, b, c, a, d2, n))
    rubi.add(rule198)

    pattern199 = Pattern(Integral((f_ + (d_ + x_**S(2)*WC('e', S(1)))**p_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))*WC('Px', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda x, u: SumQ(u)))
    def With199(Px, p, x, b, m, g, c, f, a, d, e, n):
        u = ExpandIntegrand(Px*(a + b*asinh(c*x))**n*(f + g*(d + e*x**S(2))**p)**m, x)
        return Int(u, x)
    rule199 = ReplacementRule(pattern199, lambda Px, p, x, b, m, g, c, f, a, d, e, n : With199(Px, p, x, b, m, g, c, f, a, d, e, n))
    rubi.add(rule199)

    pattern200 = Pattern(Integral((f_ + (d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*WC('Px', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda x, Px: PolynomialQ(Px, x)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda x, u: SumQ(u)))
    def With200(Px, p, e2, e1, d1, x, b, m, g, c, f, a, d2, n):
        u = ExpandIntegrand(Px*(a + b*acosh(c*x))**n*(f + g*(d1 + e1*x)**p*(d2 + e2*x)**p)**m, x)
        return Int(u, x)
    rule200 = ReplacementRule(pattern200, lambda Px, p, e2, e1, d1, x, b, m, g, c, f, a, d2, n : With200(Px, p, e2, e1, d1, x, b, m, g, c, f, a, d2, n))
    rubi.add(rule200)

    pattern201 = Pattern(Integral(RFx_*asinh(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda x, u: SumQ(u)))
    def With201(RFx, c, n, x):
        u = ExpandIntegrand(asinh(c*x)**n, RFx, x)
        return Int(u, x)
    rule201 = ReplacementRule(pattern201, lambda RFx, c, n, x : With201(RFx, c, n, x))
    rubi.add(rule201)

    pattern202 = Pattern(Integral(RFx_*acosh(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda x, u: SumQ(u)))
    def With202(RFx, c, n, x):
        u = ExpandIntegrand(acosh(c*x)**n, RFx, x)
        return Int(u, x)
    rule202 = ReplacementRule(pattern202, lambda RFx, c, n, x : With202(RFx, c, n, x))
    rubi.add(rule202)

    pattern203 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule203 = ReplacementRule(pattern203, lambda RFx, x, b, c, a, n : Int(ExpandIntegrand(RFx*(a + b*asinh(c*x))**n, x), x))
    rubi.add(rule203)

    pattern204 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule204 = ReplacementRule(pattern204, lambda RFx, x, b, c, a, n : Int(ExpandIntegrand(RFx*(a + b*acosh(c*x))**n, x), x))
    rubi.add(rule204)

    pattern205 = Pattern(Integral(RFx_*(d_ + x_**S(2)*WC('e', S(1)))**p_*asinh(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda x, u: SumQ(u)))
    def With205(RFx, p, x, c, d, e, n):
        u = ExpandIntegrand((d + e*x**S(2))**p*asinh(c*x)**n, RFx, x)
        return Int(u, x)
    rule205 = ReplacementRule(pattern205, lambda RFx, p, x, c, d, e, n : With205(RFx, p, x, c, d, e, n))
    rubi.add(rule205)

    pattern206 = Pattern(Integral(RFx_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*acosh(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda x, u: SumQ(u)))
    def With206(RFx, p, e2, e1, d1, x, c, d2, n):
        u = ExpandIntegrand((d1 + e1*x)**p*(d2 + e2*x)**p*acosh(c*x)**n, RFx, x)
        return Int(u, x)
    rule206 = ReplacementRule(pattern206, lambda RFx, p, e2, e1, d1, x, c, d2, n : With206(RFx, p, e2, e1, d1, x, c, d2, n))
    rubi.add(rule206)

    pattern207 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule207 = ReplacementRule(pattern207, lambda RFx, p, x, b, c, a, d, e, n : Int(ExpandIntegrand((d + e*x**S(2))**p, RFx*(a + b*asinh(c*x))**n, x), x))
    rubi.add(rule207)

    pattern208 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d1, x: FreeQ(d1, x)), CustomConstraint(lambda e1, x: FreeQ(e1, x)), CustomConstraint(lambda d2, x: FreeQ(d2, x)), CustomConstraint(lambda e2, x: FreeQ(e2, x)), CustomConstraint(lambda RFx, x: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda c, d1, e1: ZeroQ(-c*d1 + e1)), CustomConstraint(lambda d2, c, e2: ZeroQ(c*d2 + e2)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule208 = ReplacementRule(pattern208, lambda RFx, p, e2, e1, d1, x, b, c, a, d2, n : Int(ExpandIntegrand((d1 + e1*x)**p*(d2 + e2*x)**p, RFx*(a + b*acosh(c*x))**n, x), x))
    rubi.add(rule208)

    pattern209 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule209 = ReplacementRule(pattern209, lambda x, b, c, a, u, n : Int(u*(a + b*asinh(c*x))**n, x))
    rubi.add(rule209)

    pattern210 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule210 = ReplacementRule(pattern210, lambda x, b, c, a, u, n : Int(u*(a + b*acosh(c*x))**n, x))
    rubi.add(rule210)

    pattern211 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule211 = ReplacementRule(pattern211, lambda x, b, c, a, d, n : Subst(Int((a + b*asinh(x))**n, x), x, c + d*x)/d)
    rubi.add(rule211)

    pattern212 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule212 = ReplacementRule(pattern212, lambda x, b, c, a, d, n : Subst(Int((a + b*acosh(x))**n, x), x, c + d*x)/d)
    rubi.add(rule212)

    pattern213 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule213 = ReplacementRule(pattern213, lambda x, b, m, f, c, a, d, e, n : Subst(Int((a + b*asinh(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule213)

    pattern214 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule214 = ReplacementRule(pattern214, lambda x, b, m, f, c, a, d, e, n : Subst(Int((a + b*acosh(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule214)

    pattern215 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(-S(2)*A*c*d + B*(c**S(2) + S(1)))), CustomConstraint(lambda B, c, C, d: ZeroQ(-B*d + S(2)*C*c)))
    rule215 = ReplacementRule(pattern215, lambda B, A, C, p, x, b, c, a, d, n : Subst(Int((a + b*asinh(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule215)

    pattern216 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda B, c, C, d: ZeroQ(-B*d + S(2)*C*c)))
    rule216 = ReplacementRule(pattern216, lambda B, A, C, p, x, b, c, a, d, n : Subst(Int((a + b*acosh(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule216)

    pattern217 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(-S(2)*A*c*d + B*(c**S(2) + S(1)))), CustomConstraint(lambda B, c, C, d: ZeroQ(-B*d + S(2)*C*c)))
    rule217 = ReplacementRule(pattern217, lambda B, A, C, p, x, b, m, f, c, a, d, e, n : Subst(Int((a + b*asinh(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule217)

    pattern218 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda B, c, C, d: ZeroQ(-B*d + S(2)*C*c)))
    rule218 = ReplacementRule(pattern218, lambda B, A, C, p, x, b, m, f, c, a, d, e, n : Subst(Int((a + b*acosh(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule218)

    pattern219 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(1))))
    rule219 = ReplacementRule(pattern219, lambda x, b, c, a, d : -sqrt(Pi)*x*(-c*sinh(a/(S(2)*b)) + Cosh(a/(S(2)*b)))*FresnelC(sqrt(-c/(Pi*b))*sqrt(a + b*asinh(c + d*x**S(2))))/(sqrt(-c/b)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + Cosh(asinh(c + d*x**S(2))/S(2)))) + sqrt(Pi)*x*(c*sinh(a/(S(2)*b)) + Cosh(a/(S(2)*b)))*FresnelS(sqrt(-c/(Pi*b))*sqrt(a + b*asinh(c + d*x**S(2))))/(sqrt(-c/b)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + Cosh(asinh(c + d*x**S(2))/S(2)))) + x*sqrt(a + b*asinh(c + d*x**S(2))))
    rubi.add(rule219)

    pattern220 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule220 = ReplacementRule(pattern220, lambda x, b, c, a, d, n : S(4)*b**S(2)*n*(n + S(-1))*Int((a + b*asinh(c + d*x**S(2)))**(n + S(-2)), x) - S(2)*b*n*(a + b*asinh(c + d*x**S(2)))**(n + S(-1))*sqrt(S(2)*c*d*x**S(2) + d**S(2)*x**S(4))/(d*x) + x*(a + b*asinh(c + d*x**S(2)))**n)
    rubi.add(rule220)

    pattern221 = Pattern(Integral(1/(WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(1))))
    rule221 = ReplacementRule(pattern221, lambda x, b, c, a, d : x*(c*Cosh(a/(S(2)*b)) - sinh(a/(S(2)*b)))*CoshIntegral((a + b*asinh(c + d*x**S(2)))/(S(2)*b))/(S(2)*b*(c*sinh(asinh(c + d*x**S(2))/S(2)) + Cosh(asinh(c + d*x**S(2))/S(2)))) + x*(-c*sinh(a/(S(2)*b)) + Cosh(a/(S(2)*b)))*SinhIntegral((a + b*asinh(c + d*x**S(2)))/(S(2)*b))/(S(2)*b*(c*sinh(asinh(c + d*x**S(2))/S(2)) + Cosh(asinh(c + d*x**S(2))/S(2)))))
    rubi.add(rule221)

    pattern222 = Pattern(Integral(1/sqrt(WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(1))))
    rule222 = ReplacementRule(pattern222, lambda x, b, c, a, d : sqrt(S(2))*sqrt(Pi)*x*(c + S(-1))*(Cosh(a/(S(2)*b)) + sinh(a/(S(2)*b)))*Erf(sqrt(S(2))*sqrt(a + b*asinh(c + d*x**S(2)))/(S(2)*sqrt(b)))/(S(4)*sqrt(b)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + Cosh(asinh(c + d*x**S(2))/S(2)))) + sqrt(S(2))*sqrt(Pi)*x*(c + S(1))*(Cosh(a/(S(2)*b)) - sinh(a/(S(2)*b)))*Erfi(sqrt(S(2))*sqrt(a + b*asinh(c + d*x**S(2)))/(S(2)*sqrt(b)))/(S(4)*sqrt(b)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + Cosh(asinh(c + d*x**S(2))/S(2)))))
    rubi.add(rule222)

    pattern223 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1))))**(S(-3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(1))))
    rule223 = ReplacementRule(pattern223, lambda x, b, c, a, d : -sqrt(Pi)*x*(-c/b)**(S(3)/2)*(-c*sinh(a/(S(2)*b)) + Cosh(a/(S(2)*b)))*FresnelC(sqrt(-c/(Pi*b))*sqrt(a + b*asinh(c + d*x**S(2))))/(c*sinh(asinh(c + d*x**S(2))/S(2)) + Cosh(asinh(c + d*x**S(2))/S(2))) + sqrt(Pi)*x*(-c/b)**(S(3)/2)*(c*sinh(a/(S(2)*b)) + Cosh(a/(S(2)*b)))*FresnelS(sqrt(-c/(Pi*b))*sqrt(a + b*asinh(c + d*x**S(2))))/(c*sinh(asinh(c + d*x**S(2))/S(2)) + Cosh(asinh(c + d*x**S(2))/S(2))) - sqrt(S(2)*c*d*x**S(2) + d**S(2)*x**S(4))/(b*d*x*sqrt(a + b*asinh(c + d*x**S(2)))))
    rubi.add(rule223)

    pattern224 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1))))**(S(-2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(1))))
    rule224 = ReplacementRule(pattern224, lambda x, b, c, a, d : -sqrt(S(2)*c*d*x**S(2) + d**S(2)*x**S(4))/(S(2)*b*d*x*(a + b*asinh(c + d*x**S(2)))) + x*(c*Cosh(a/(S(2)*b)) - sinh(a/(S(2)*b)))*SinhIntegral((a + b*asinh(c + d*x**S(2)))/(S(2)*b))/(S(4)*b**S(2)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + Cosh(asinh(c + d*x**S(2))/S(2)))) + x*(-c*sinh(a/(S(2)*b)) + Cosh(a/(S(2)*b)))*CoshIntegral((a + b*asinh(c + d*x**S(2)))/(S(2)*b))/(S(4)*b**S(2)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + Cosh(asinh(c + d*x**S(2))/S(2)))))
    rubi.add(rule224)

    pattern225 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule225 = ReplacementRule(pattern225, lambda x, b, c, a, d, n : (a + b*asinh(c + d*x**S(2)))**(n + S(1))*sqrt(S(2)*c*d*x**S(2) + d**S(2)*x**S(4))/(S(2)*b*d*x*(n + S(1))) - x*(a + b*asinh(c + d*x**S(2)))**(n + S(2))/(S(4)*b**S(2)*(n + S(1))*(n + S(2))) + Int((a + b*asinh(c + d*x**S(2)))**(n + S(2)), x)/(S(4)*b**S(2)*(n + S(1))*(n + S(2))))
    rubi.add(rule225)

    pattern226 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule226 = ReplacementRule(pattern226, lambda x, a, d, b : -sqrt(S(2))*sqrt(Pi)*sqrt(b)*(Cosh(a/(S(2)*b)) - sinh(a/(S(2)*b)))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*d*x) + sqrt(S(2))*sqrt(Pi)*sqrt(b)*(Cosh(a/(S(2)*b)) + sinh(a/(S(2)*b)))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*d*x) + S(2)*sqrt(a + b*acosh(d*x**S(2) + S(1)))*sinh(acosh(d*x**S(2) + S(1))/S(2))**S(2)/(d*x))
    rubi.add(rule226)

    pattern227 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(-1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule227 = ReplacementRule(pattern227, lambda x, a, d, b : -sqrt(S(2))*sqrt(Pi)*sqrt(b)*(Cosh(a/(S(2)*b)) - sinh(a/(S(2)*b)))*Cosh(acosh(d*x**S(2) + S(-1))/S(2))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))/(S(2)*d*x) - sqrt(S(2))*sqrt(Pi)*sqrt(b)*(Cosh(a/(S(2)*b)) + sinh(a/(S(2)*b)))*Cosh(acosh(d*x**S(2) + S(-1))/S(2))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))/(S(2)*d*x) + S(2)*sqrt(a + b*acosh(d*x**S(2) + S(-1)))*Cosh(acosh(d*x**S(2) + S(-1))/S(2))**S(2)/(d*x))
    rubi.add(rule227)

    pattern228 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule228 = ReplacementRule(pattern228, lambda x, b, c, a, d, n : S(4)*b**S(2)*n*(n + S(-1))*Int((a + b*acosh(c + d*x**S(2)))**(n + S(-2)), x) - S(2)*b*n*(a + b*acosh(c + d*x**S(2)))**(n + S(-1))*(S(2)*c*d*x**S(2) + d**S(2)*x**S(4))/(d*x*sqrt(c + d*x**S(2) + S(-1))*sqrt(c + d*x**S(2) + S(1))) + x*(a + b*acosh(c + d*x**S(2)))**n)
    rubi.add(rule228)

    pattern229 = Pattern(Integral(1/(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule229 = ReplacementRule(pattern229, lambda x, a, d, b : sqrt(S(2))*x*Cosh(a/(S(2)*b))*CoshIntegral((a + b*acosh(d*x**S(2) + S(1)))/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))) - sqrt(S(2))*x*SinhIntegral((a + b*acosh(d*x**S(2) + S(1)))/(S(2)*b))*sinh(a/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))))
    rubi.add(rule229)

    pattern230 = Pattern(Integral(1/(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(-1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule230 = ReplacementRule(pattern230, lambda x, a, d, b : sqrt(S(2))*x*Cosh(a/(S(2)*b))*SinhIntegral((a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))) - sqrt(S(2))*x*CoshIntegral((a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*b))*sinh(a/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))))
    rubi.add(rule230)

    pattern231 = Pattern(Integral(1/sqrt(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule231 = ReplacementRule(pattern231, lambda x, a, d, b : sqrt(S(2))*sqrt(Pi)*(Cosh(a/(S(2)*b)) - sinh(a/(S(2)*b)))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*sqrt(b)*d*x) + sqrt(S(2))*sqrt(Pi)*(Cosh(a/(S(2)*b)) + sinh(a/(S(2)*b)))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*sqrt(b)*d*x))
    rubi.add(rule231)

    pattern232 = Pattern(Integral(1/sqrt(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(-1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule232 = ReplacementRule(pattern232, lambda x, a, d, b : sqrt(S(2))*sqrt(Pi)*(Cosh(a/(S(2)*b)) - sinh(a/(S(2)*b)))*Cosh(acosh(d*x**S(2) + S(-1))/S(2))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))/(S(2)*sqrt(b)*d*x) - sqrt(S(2))*sqrt(Pi)*(Cosh(a/(S(2)*b)) + sinh(a/(S(2)*b)))*Cosh(acosh(d*x**S(2) + S(-1))/S(2))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))/(S(2)*sqrt(b)*d*x))
    rubi.add(rule232)

    pattern233 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(1)))**(S(-3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule233 = ReplacementRule(pattern233, lambda x, a, d, b : sqrt(S(2))*sqrt(Pi)*(Cosh(a/(S(2)*b)) - sinh(a/(S(2)*b)))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*b**(S(3)/2)*d*x) - sqrt(S(2))*sqrt(Pi)*(Cosh(a/(S(2)*b)) + sinh(a/(S(2)*b)))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*b**(S(3)/2)*d*x) - sqrt(d*x**S(2))*sqrt(d*x**S(2) + S(2))/(b*d*x*sqrt(a + b*acosh(d*x**S(2) + S(1)))))
    rubi.add(rule233)

    pattern234 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(-1)))**(S(-3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule234 = ReplacementRule(pattern234, lambda x, a, d, b : sqrt(S(2))*sqrt(Pi)*(Cosh(a/(S(2)*b)) - sinh(a/(S(2)*b)))*Cosh(acosh(d*x**S(2) + S(-1))/S(2))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))/(S(2)*b**(S(3)/2)*d*x) + sqrt(S(2))*sqrt(Pi)*(Cosh(a/(S(2)*b)) + sinh(a/(S(2)*b)))*Cosh(acosh(d*x**S(2) + S(-1))/S(2))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))/(S(2)*b**(S(3)/2)*d*x) - sqrt(d*x**S(2))*sqrt(d*x**S(2) + S(-2))/(b*d*x*sqrt(a + b*acosh(d*x**S(2) + S(-1)))))
    rubi.add(rule234)

    pattern235 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(1)))**(S(-2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule235 = ReplacementRule(pattern235, lambda x, a, d, b : -sqrt(d*x**S(2))*sqrt(d*x**S(2) + S(2))/(S(2)*b*d*x*(a + b*acosh(d*x**S(2) + S(1)))) + sqrt(S(2))*x*Cosh(a/(S(2)*b))*SinhIntegral((a + b*acosh(d*x**S(2) + S(1)))/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))) - sqrt(S(2))*x*CoshIntegral((a + b*acosh(d*x**S(2) + S(1)))/(S(2)*b))*sinh(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))))
    rubi.add(rule235)

    pattern236 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(-1)))**(S(-2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule236 = ReplacementRule(pattern236, lambda x, a, d, b : -sqrt(d*x**S(2))*sqrt(d*x**S(2) + S(-2))/(S(2)*b*d*x*(a + b*acosh(d*x**S(2) + S(-1)))) + sqrt(S(2))*x*Cosh(a/(S(2)*b))*CoshIntegral((a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))) - sqrt(S(2))*x*SinhIntegral((a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*b))*sinh(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))))
    rubi.add(rule236)

    pattern237 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule237 = ReplacementRule(pattern237, lambda x, b, c, a, d, n : (a + b*acosh(c + d*x**S(2)))**(n + S(1))*(S(2)*c*x**S(2) + d*x**S(4))/(S(2)*b*x*(n + S(1))*sqrt(c + d*x**S(2) + S(-1))*sqrt(c + d*x**S(2) + S(1))) - x*(a + b*acosh(c + d*x**S(2)))**(n + S(2))/(S(4)*b**S(2)*(n + S(1))*(n + S(2))) + Int((a + b*acosh(c + d*x**S(2)))**(n + S(2)), x)/(S(4)*b**S(2)*(n + S(1))*(n + S(2))))
    rubi.add(rule237)

    pattern238 = Pattern(Integral(asinh(x_**p_*WC('a', S(1)))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule238 = ReplacementRule(pattern238, lambda x, a, p, n : Subst(Int(x**n*coth(x), x), x, asinh(a*x**p))/p)
    rubi.add(rule238)

    pattern239 = Pattern(Integral(acosh(x_**p_*WC('a', S(1)))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule239 = ReplacementRule(pattern239, lambda x, a, p, n : Subst(Int(x**n*tanh(x), x), x, acosh(a*x**p))/p)
    rubi.add(rule239)

    pattern240 = Pattern(Integral(WC('u', S(1))*asinh(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule240 = ReplacementRule(pattern240, lambda x, b, m, c, a, u, n : Int(u*acsch(a/c + b*x**n/c)**m, x))
    rubi.add(rule240)

    pattern241 = Pattern(Integral(WC('u', S(1))*acosh(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule241 = ReplacementRule(pattern241, lambda x, b, m, c, a, u, n : Int(u*asech(a/c + b*x**n/c)**m, x))
    rubi.add(rule241)

    pattern242 = Pattern(Integral(asinh(sqrt(x_**S(2)*WC('b', S(1)) + S(-1)))**WC('n', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + S(-1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule242 = ReplacementRule(pattern242, lambda x, n, b : sqrt(b*x**S(2))*Subst(Int(asinh(x)**n/sqrt(x**S(2) + S(1)), x), x, sqrt(b*x**S(2) + S(-1)))/(b*x))
    rubi.add(rule242)

    pattern243 = Pattern(Integral(acosh(sqrt(x_**S(2)*WC('b', S(1)) + S(1)))**WC('n', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule243 = ReplacementRule(pattern243, lambda x, n, b : sqrt(sqrt(b*x**S(2) + S(1)) + S(-1))*sqrt(sqrt(b*x**S(2) + S(1)) + S(1))*Subst(Int(acosh(x)**n/(sqrt(x + S(-1))*sqrt(x + S(1))), x), x, sqrt(b*x**S(2) + S(1)))/(b*x))
    rubi.add(rule243)

    pattern244 = Pattern(Integral(f_**(WC('c', S(1))*asinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule244 = ReplacementRule(pattern244, lambda x, b, c, f, a, n : Subst(Int(f**(c*x**n)*Cosh(x), x), x, asinh(a + b*x))/b)
    rubi.add(rule244)

    pattern245 = Pattern(Integral(f_**(WC('c', S(1))*acosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule245 = ReplacementRule(pattern245, lambda x, b, c, f, a, n : Subst(Int(f**(c*x**n)*sinh(x), x), x, acosh(a + b*x))/b)
    rubi.add(rule245)

    pattern246 = Pattern(Integral(f_**(WC('c', S(1))*asinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*x_**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule246 = ReplacementRule(pattern246, lambda x, b, m, c, f, a, n : Subst(Int(f**(c*x**n)*(-a/b + sinh(x)/b)**m*Cosh(x), x), x, asinh(a + b*x))/b)
    rubi.add(rule246)

    pattern247 = Pattern(Integral(f_**(WC('c', S(1))*acosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*x_**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule247 = ReplacementRule(pattern247, lambda x, b, m, c, f, a, n : Subst(Int(f**(c*x**n)*(-a/b + Cosh(x)/b)**m*sinh(x), x), x, acosh(a + b*x))/b)
    rubi.add(rule247)

    pattern248 = Pattern(Integral(asinh(u_), x_), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda x, u: Not(FunctionOfExponentialQ(u, x))))
    rule248 = ReplacementRule(pattern248, lambda x, u : x*asinh(u) - Int(SimplifyIntegrand(x*D(u, x)/sqrt(u**S(2) + S(1)), x), x))
    rubi.add(rule248)

    pattern249 = Pattern(Integral(acosh(u_), x_), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda x, u: Not(FunctionOfExponentialQ(u, x))))
    rule249 = ReplacementRule(pattern249, lambda x, u : x*acosh(u) - Int(SimplifyIntegrand(x*D(u, x)/(sqrt(u + S(-1))*sqrt(u + S(1))), x), x))
    rubi.add(rule249)

    pattern250 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(u_)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda x, m, c, d, u: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda x, u: Not(FunctionOfExponentialQ(u, x))))
    rule250 = ReplacementRule(pattern250, lambda x, b, m, c, a, d, u : -b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/sqrt(u**S(2) + S(1)), x), x)/(d*(m + S(1))) + (a + b*asinh(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule250)

    pattern251 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(u_)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda x, m, c, d, u: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda x, u: Not(FunctionOfExponentialQ(u, x))))
    rule251 = ReplacementRule(pattern251, lambda x, b, m, c, a, d, u : -b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(sqrt(u + S(-1))*sqrt(u + S(1))), x), x)/(d*(m + S(1))) + (a + b*acosh(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule251)

    pattern252 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*asinh(u_)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda a, u, x, b, w: InverseFunctionFreeQ(w, x)))
    def With252(x, b, v, a, u):
        w = IntHide(v, x)
        return -b*Int(SimplifyIntegrand(w*D(u, x)/sqrt(u**S(2) + S(1)), x), x) + Dist(a + b*asinh(u), w, x)
    rule252 = ReplacementRule(pattern252, lambda x, b, v, a, u : With252(x, b, v, a, u))
    rubi.add(rule252)

    pattern253 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acosh(u_)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda a, u, x, b, w: InverseFunctionFreeQ(w, x)))
    def With253(x, b, v, a, u):
        w = IntHide(v, x)
        return -b*Int(SimplifyIntegrand(w*D(u, x)/(sqrt(u + S(-1))*sqrt(u + S(1))), x), x) + Dist(a + b*acosh(u), w, x)
    rule253 = ReplacementRule(pattern253, lambda x, b, v, a, u : With253(x, b, v, a, u))
    rubi.add(rule253)

    pattern254 = Pattern(Integral(exp(WC('n', S(1))*asinh(u_)), x_), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda x, u: PolynomialQ(u, x)))
    rule254 = ReplacementRule(pattern254, lambda x, u, n : Int((u + sqrt(u**S(2) + S(1)))**n, x))
    rubi.add(rule254)

    pattern255 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*asinh(u_)), x_), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda x, u: PolynomialQ(u, x)))
    rule255 = ReplacementRule(pattern255, lambda x, u, n, m : Int(x**m*(u + sqrt(u**S(2) + S(1)))**n, x))
    rubi.add(rule255)

    pattern256 = Pattern(Integral(exp(WC('n', S(1))*acosh(u_)), x_), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda x, u: PolynomialQ(u, x)))
    rule256 = ReplacementRule(pattern256, lambda x, u, n : Int((u + sqrt(u + S(-1))*sqrt(u + S(1)))**n, x))
    rubi.add(rule256)

    pattern257 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*acosh(u_)), x_), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda x, u: PolynomialQ(u, x)))
    rule257 = ReplacementRule(pattern257, lambda x, u, n, m : Int(x**m*(u + sqrt(u + S(-1))*sqrt(u + S(1)))**n, x))
    rubi.add(rule257)

    pattern258 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule258 = ReplacementRule(pattern258, lambda x, b, c, a, n : -b*c*n*Int(x*(a + b*atanh(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x) + x*(a + b*atanh(c*x))**n)
    rubi.add(rule258)

    pattern259 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule259 = ReplacementRule(pattern259, lambda x, b, c, a, n : -b*c*n*Int(x*(a + b*acoth(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x) + x*(a + b*acoth(c*x))**n)
    rubi.add(rule259)

    pattern260 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule260 = ReplacementRule(pattern260, lambda x, b, c, a, n : Int((a + b*atanh(c*x))**n, x))
    rubi.add(rule260)

    pattern261 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule261 = ReplacementRule(pattern261, lambda x, b, c, a, n : Int((a + b*acoth(c*x))**n, x))
    rubi.add(rule261)

    pattern262 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d**S(2) - e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule262 = ReplacementRule(pattern262, lambda x, b, c, a, d, e, n : b*c*n*Int((a + b*atanh(c*x))**(n + S(-1))*log(S(2)*d/(d + e*x))/(-c**S(2)*x**S(2) + S(1)), x)/e - (a + b*atanh(c*x))**n*log(S(2)*d/(d + e*x))/e)
    rubi.add(rule262)

    pattern263 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d**S(2) - e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule263 = ReplacementRule(pattern263, lambda x, b, c, a, d, e, n : b*c*n*Int((a + b*acoth(c*x))**(n + S(-1))*log(S(2)*d/(d + e*x))/(-c**S(2)*x**S(2) + S(1)), x)/e - (a + b*acoth(c*x))**n*log(S(2)*d/(d + e*x))/e)
    rubi.add(rule263)

    pattern264 = Pattern(Integral(atanh(x_*WC('c', S(1)))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: PositiveQ(c*d/e + S(1))), CustomConstraint(lambda c, e, d: NegativeQ(c*d/e + S(-1))))
    rule264 = ReplacementRule(pattern264, lambda c, d, e, x : -PolyLog(S(2), Simp(c*(d + e*x)/(c*d - e), x))/(S(2)*e) + PolyLog(S(2), Simp(c*(d + e*x)/(c*d + e), x))/(S(2)*e) - log(d + e*x)*atanh(c*d/e)/e)
    rubi.add(rule264)

    pattern265 = Pattern(Integral(atanh(x_*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule265 = ReplacementRule(pattern265, lambda c, e, d, x : -Int(log(-c*x + S(1))/(d + e*x), x)/S(2) + Int(log(c*x + S(1))/(d + e*x), x)/S(2))
    rubi.add(rule265)

    pattern266 = Pattern(Integral(acoth(x_*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule266 = ReplacementRule(pattern266, lambda c, e, d, x : -Int(log(S(1) - S(1)/(c*x))/(d + e*x), x)/S(2) + Int(log(S(1) + S(1)/(c*x))/(d + e*x), x)/S(2))
    rubi.add(rule266)

    pattern267 = Pattern(Integral((a_ + WC('b', S(1))*atanh(x_*WC('c', S(1))))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule267 = ReplacementRule(pattern267, lambda x, b, c, a, d, e : a*log(RemoveContent(d + e*x, x))/e + b*Int(atanh(c*x)/(d + e*x), x))
    rubi.add(rule267)

    pattern268 = Pattern(Integral((a_ + WC('b', S(1))*acoth(x_*WC('c', S(1))))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule268 = ReplacementRule(pattern268, lambda x, b, c, a, d, e : a*log(RemoveContent(d + e*x, x))/e + b*Int(acoth(c*x)/(d + e*x), x))
    rubi.add(rule268)

    pattern269 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule269 = ReplacementRule(pattern269, lambda p, x, b, c, a, d, e : -b*c*Int((d + e*x)**(p + S(1))/(-c**S(2)*x**S(2) + S(1)), x)/(e*(p + S(1))) + (a + b*atanh(c*x))*(d + e*x)**(p + S(1))/(e*(p + S(1))))
    rubi.add(rule269)

    pattern270 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule270 = ReplacementRule(pattern270, lambda p, x, b, c, a, d, e : -b*c*Int((d + e*x)**(p + S(1))/(-c**S(2)*x**S(2) + S(1)), x)/(e*(p + S(1))) + (a + b*acoth(c*x))*(d + e*x)**(p + S(1))/(e*(p + S(1))))
    rubi.add(rule270)

    pattern271 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule271 = ReplacementRule(pattern271, lambda x, b, c, a, n : -S(2)*b*c*n*Int((a + b*atanh(c*x))**(n + S(-1))*atanh(S(1) - S(2)/(-c*x + S(1)))/(-c**S(2)*x**S(2) + S(1)), x) + S(2)*(a + b*atanh(c*x))**n*atanh(S(1) - S(2)/(-c*x + S(1))))
    rubi.add(rule271)

    pattern272 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule272 = ReplacementRule(pattern272, lambda x, b, c, a, n : -S(2)*b*c*n*Int((a + b*acoth(c*x))**(n + S(-1))*acoth(S(1) - S(2)/(-c*x + S(1)))/(-c**S(2)*x**S(2) + S(1)), x) + S(2)*(a + b*acoth(c*x))**n*acoth(S(1) - S(2)/(-c*x + S(1))))
    rubi.add(rule272)

    pattern273 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule273 = ReplacementRule(pattern273, lambda x, b, m, c, a, n : -b*c*n*Int(x**(m + S(1))*(a + b*atanh(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*atanh(c*x))**n/(m + S(1)))
    rubi.add(rule273)

    pattern274 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule274 = ReplacementRule(pattern274, lambda x, b, m, c, a, n : -b*c*n*Int(x**(m + S(1))*(a + b*acoth(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*acoth(c*x))**n/(m + S(1)))
    rubi.add(rule274)

    pattern275 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)))
    rule275 = ReplacementRule(pattern275, lambda p, x, b, c, a, d, e, n : Int(ExpandIntegrand((a + b*atanh(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule275)

    pattern276 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, n: PositiveIntegerQ(n, p)))
    rule276 = ReplacementRule(pattern276, lambda p, x, b, c, a, d, e, n : Int(ExpandIntegrand((a + b*acoth(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule276)

    pattern277 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule277 = ReplacementRule(pattern277, lambda p, x, b, c, a, d, e, n : Int((a + b*atanh(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule277)

    pattern278 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule278 = ReplacementRule(pattern278, lambda p, x, b, c, a, d, e, n : Int((a + b*acoth(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule278)

    pattern279 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d**S(2) - e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule279 = ReplacementRule(pattern279, lambda x, b, m, c, a, d, e, n : -d*Int(x**(m + S(-1))*(a + b*atanh(c*x))**n/(d + e*x), x)/e + Int(x**(m + S(-1))*(a + b*atanh(c*x))**n, x)/e)
    rubi.add(rule279)

    pattern280 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d**S(2) - e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule280 = ReplacementRule(pattern280, lambda x, b, m, c, a, d, e, n : -d*Int(x**(m + S(-1))*(a + b*acoth(c*x))**n/(d + e*x), x)/e + Int(x**(m + S(-1))*(a + b*acoth(c*x))**n, x)/e)
    rubi.add(rule280)

    pattern281 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d**S(2) - e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule281 = ReplacementRule(pattern281, lambda x, b, c, a, d, e, n : -b*c*n*Int((a + b*atanh(c*x))**(n + S(-1))*log(S(2)*e*x/(d + e*x))/(-c**S(2)*x**S(2) + S(1)), x)/d + (a + b*atanh(c*x))**n*log(S(2)*e*x/(d + e*x))/d)
    rubi.add(rule281)

    pattern282 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d**S(2) - e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule282 = ReplacementRule(pattern282, lambda x, b, c, a, d, e, n : -b*c*n*Int((a + b*acoth(c*x))**(n + S(-1))*log(S(2)*e*x/(d + e*x))/(-c**S(2)*x**S(2) + S(1)), x)/d + (a + b*acoth(c*x))**n*log(S(2)*e*x/(d + e*x))/d)
    rubi.add(rule282)

    pattern283 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d**S(2) - e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule283 = ReplacementRule(pattern283, lambda x, b, m, c, a, d, e, n : -e*Int(x**(m + S(1))*(a + b*atanh(c*x))**n/(d + e*x), x)/d + Int(x**m*(a + b*atanh(c*x))**n, x)/d)
    rubi.add(rule283)

    pattern284 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d**S(2) - e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule284 = ReplacementRule(pattern284, lambda x, b, m, c, a, d, e, n : -e*Int(x**(m + S(1))*(a + b*acoth(c*x))**n/(d + e*x), x)/d + Int(x**m*(a + b*acoth(c*x))**n, x)/d)
    rubi.add(rule284)

    pattern285 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda a, p, m: IntegerQ(m) | NonzeroQ(a) | Greater(p, S(0))))
    rule285 = ReplacementRule(pattern285, lambda p, x, b, m, c, a, d, e, n : Int(ExpandIntegrand(x**m*(a + b*atanh(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule285)

    pattern286 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda a, p, m: IntegerQ(m) | NonzeroQ(a) | Greater(p, S(0))))
    rule286 = ReplacementRule(pattern286, lambda p, x, b, m, c, a, d, e, n : Int(ExpandIntegrand(x**m*(a + b*acoth(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule286)

    pattern287 = Pattern(Integral(x_**WC('m', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule287 = ReplacementRule(pattern287, lambda p, x, b, m, c, a, d, e, n : Int(x**m*(a + b*atanh(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule287)

    pattern288 = Pattern(Integral(x_**WC('m', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule288 = ReplacementRule(pattern288, lambda p, x, b, m, c, a, d, e, n : Int(x**m*(a + b*acoth(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule288)

    pattern289 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule289 = ReplacementRule(pattern289, lambda p, x, b, c, a, d, e : b*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*atanh(c*x))*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule289)

    pattern290 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule290 = ReplacementRule(pattern290, lambda p, x, b, c, a, d, e : b*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*acoth(c*x))*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule290)

    pattern291 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule291 = ReplacementRule(pattern291, lambda p, x, b, c, a, d, e, n : -b**S(2)*d*n*(n + S(-1))*Int((a + b*atanh(c*x))**(n + S(-2))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p*(S(2)*p + S(1))) + b*n*(a + b*atanh(c*x))**(n + S(-1))*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*atanh(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule291)

    pattern292 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule292 = ReplacementRule(pattern292, lambda p, x, b, c, a, d, e, n : -b**S(2)*d*n*(n + S(-1))*Int((a + b*acoth(c*x))**(n + S(-2))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p*(S(2)*p + S(1))) + b*n*(a + b*acoth(c*x))**(n + S(-1))*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*acoth(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule292)

    pattern293 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)))
    rule293 = ReplacementRule(pattern293, lambda x, b, c, a, d, e : log(RemoveContent(a + b*atanh(c*x), x))/(b*c*d))
    rubi.add(rule293)

    pattern294 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)))
    rule294 = ReplacementRule(pattern294, lambda x, b, c, a, d, e : log(RemoveContent(a + b*acoth(c*x), x))/(b*c*d))
    rubi.add(rule294)

    pattern295 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule295 = ReplacementRule(pattern295, lambda x, b, c, a, d, e, n : (a + b*atanh(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule295)

    pattern296 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule296 = ReplacementRule(pattern296, lambda x, b, c, a, d, e, n : (a + b*acoth(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule296)

    pattern297 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule297 = ReplacementRule(pattern297, lambda x, b, c, a, d, e : -ImaginaryI*b*PolyLog(S(2), -ImaginaryI*sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)) + ImaginaryI*b*PolyLog(S(2), ImaginaryI*sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)) - S(2)*(a + b*atanh(c*x))*ArcTan(sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)))
    rubi.add(rule297)

    pattern298 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule298 = ReplacementRule(pattern298, lambda x, b, c, a, d, e : -ImaginaryI*b*PolyLog(S(2), -ImaginaryI*sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)) + ImaginaryI*b*PolyLog(S(2), ImaginaryI*sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)) - S(2)*(a + b*acoth(c*x))*ArcTan(sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)))
    rubi.add(rule298)

    pattern299 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule299 = ReplacementRule(pattern299, lambda x, b, c, a, d, e, n : Subst(Int((a + b*x)**n*sech(x), x), x, atanh(c*x))/(c*sqrt(d)))
    rubi.add(rule299)

    pattern300 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule300 = ReplacementRule(pattern300, lambda x, b, c, a, d, e, n : -x*sqrt(S(1) - S(1)/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*csch(x), x), x, acoth(c*x))/sqrt(d + e*x**S(2)))
    rubi.add(rule300)

    pattern301 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule301 = ReplacementRule(pattern301, lambda x, b, c, a, d, e, n : sqrt(-c**S(2)*x**S(2) + S(1))*Int((a + b*atanh(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule301)

    pattern302 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule302 = ReplacementRule(pattern302, lambda x, b, c, a, d, e, n : sqrt(-c**S(2)*x**S(2) + S(1))*Int((a + b*acoth(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule302)

    pattern303 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule303 = ReplacementRule(pattern303, lambda x, b, c, a, d, e, n : -b*c*n*Int(x*(a + b*atanh(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/S(2) + x*(a + b*atanh(c*x))**n/(S(2)*d*(d + e*x**S(2))) + (a + b*atanh(c*x))**(n + S(1))/(S(2)*b*c*d**S(2)*(n + S(1))))
    rubi.add(rule303)

    pattern304 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule304 = ReplacementRule(pattern304, lambda x, b, c, a, d, e, n : -b*c*n*Int(x*(a + b*acoth(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/S(2) + x*(a + b*acoth(c*x))**n/(S(2)*d*(d + e*x**S(2))) + (a + b*acoth(c*x))**(n + S(1))/(S(2)*b*c*d**S(2)*(n + S(1))))
    rubi.add(rule304)

    pattern305 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)))
    rule305 = ReplacementRule(pattern305, lambda x, b, c, a, d, e : -b/(c*d*sqrt(d + e*x**S(2))) + x*(a + b*atanh(c*x))/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule305)

    pattern306 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)))
    rule306 = ReplacementRule(pattern306, lambda x, b, c, a, d, e : -b/(c*d*sqrt(d + e*x**S(2))) + x*(a + b*acoth(c*x))/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule306)

    pattern307 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule307 = ReplacementRule(pattern307, lambda p, x, b, c, a, d, e : -b*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule307)

    pattern308 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule308 = ReplacementRule(pattern308, lambda p, x, b, c, a, d, e : -b*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule308)

    pattern309 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule309 = ReplacementRule(pattern309, lambda x, b, c, a, d, e, n : b**S(2)*n*(n + S(-1))*Int((a + b*atanh(c*x))**(n + S(-2))/(d + e*x**S(2))**(S(3)/2), x) - b*n*(a + b*atanh(c*x))**(n + S(-1))/(c*d*sqrt(d + e*x**S(2))) + x*(a + b*atanh(c*x))**n/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule309)

    pattern310 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule310 = ReplacementRule(pattern310, lambda x, b, c, a, d, e, n : b**S(2)*n*(n + S(-1))*Int((a + b*acoth(c*x))**(n + S(-2))/(d + e*x**S(2))**(S(3)/2), x) - b*n*(a + b*acoth(c*x))**(n + S(-1))/(c*d*sqrt(d + e*x**S(2))) + x*(a + b*acoth(c*x))**n/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule310)

    pattern311 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule311 = ReplacementRule(pattern311, lambda p, x, b, c, a, d, e, n : b**S(2)*n*(n + S(-1))*Int((a + b*atanh(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/(S(4)*(p + S(1))**S(2)) - b*n*(a + b*atanh(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule311)

    pattern312 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule312 = ReplacementRule(pattern312, lambda p, x, b, c, a, d, e, n : b**S(2)*n*(n + S(-1))*Int((a + b*acoth(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/(S(4)*(p + S(1))**S(2)) - b*n*(a + b*acoth(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule312)

    pattern313 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))))
    rule313 = ReplacementRule(pattern313, lambda p, x, b, c, a, d, e, n : S(2)*c*(p + S(1))*Int(x*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) + (a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule313)

    pattern314 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))))
    rule314 = ReplacementRule(pattern314, lambda p, x, b, c, a, d, e, n : S(2)*c*(p + S(1))*Int(x*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) + (a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule314)

    pattern315 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda d, p: IntegerQ(p) | PositiveQ(d)))
    rule315 = ReplacementRule(pattern315, lambda p, x, b, c, a, d, e, n : d**p*Subst(Int((a + b*x)**n*Cosh(x)**(-S(2)*p + S(-2)), x), x, atanh(c*x))/c)
    rubi.add(rule315)

    pattern316 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda d, p: Not(IntegerQ(p) | PositiveQ(d))))
    rule316 = ReplacementRule(pattern316, lambda p, x, b, c, a, d, e, n : d**(p + S(1)/2)*sqrt(-c**S(2)*x**S(2) + S(1))*Int((a + b*atanh(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x)/sqrt(d + e*x**S(2)))
    rubi.add(rule316)

    pattern317 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule317 = ReplacementRule(pattern317, lambda p, x, b, c, a, d, e, n : -(-d)**p*Subst(Int((a + b*x)**n*sinh(x)**(-S(2)*p + S(-2)), x), x, acoth(c*x))/c)
    rubi.add(rule317)

    pattern318 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule318 = ReplacementRule(pattern318, lambda p, x, b, c, a, d, e, n : -x*(-d)**(p + S(1)/2)*sqrt((c**S(2)*x**S(2) + S(-1))/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*sinh(x)**(-S(2)*p + S(-2)), x), x, acoth(c*x))/sqrt(d + e*x**S(2)))
    rubi.add(rule318)

    pattern319 = Pattern(Integral(atanh(x_*WC('c', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule319 = ReplacementRule(pattern319, lambda c, e, d, x : -Int(log(-c*x + S(1))/(d + e*x**S(2)), x)/S(2) + Int(log(c*x + S(1))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule319)

    pattern320 = Pattern(Integral(acoth(x_*WC('c', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule320 = ReplacementRule(pattern320, lambda c, e, d, x : -Int(log(S(1) - S(1)/(c*x))/(d + e*x**S(2)), x)/S(2) + Int(log(S(1) + S(1)/(c*x))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule320)

    pattern321 = Pattern(Integral((a_ + WC('b', S(1))*atanh(x_*WC('c', S(1))))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule321 = ReplacementRule(pattern321, lambda x, b, c, a, d, e : a*Int(1/(d + e*x**S(2)), x) + b*Int(atanh(c*x)/(d + e*x**S(2)), x))
    rubi.add(rule321)

    pattern322 = Pattern(Integral((a_ + WC('b', S(1))*acoth(x_*WC('c', S(1))))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule322 = ReplacementRule(pattern322, lambda x, b, c, a, d, e : a*Int(1/(d + e*x**S(2)), x) + b*Int(acoth(c*x)/(d + e*x**S(2)), x))
    rubi.add(rule322)

    pattern323 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With323(p, x, b, c, a, d, e):
        u = IntHide((d + e*x**S(2))**p, x)
        return -b*c*Int(ExpandIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*atanh(c*x), u, x)
    rule323 = ReplacementRule(pattern323, lambda p, x, b, c, a, d, e : With323(p, x, b, c, a, d, e))
    rubi.add(rule323)

    pattern324 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With324(p, x, b, c, a, d, e):
        u = IntHide((d + e*x**S(2))**p, x)
        return -b*c*Int(ExpandIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acoth(c*x), u, x)
    rule324 = ReplacementRule(pattern324, lambda p, x, b, c, a, d, e : With324(p, x, b, c, a, d, e))
    rubi.add(rule324)

    pattern325 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule325 = ReplacementRule(pattern325, lambda p, x, b, c, a, d, e, n : Int(ExpandIntegrand((a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule325)

    pattern326 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule326 = ReplacementRule(pattern326, lambda p, x, b, c, a, d, e, n : Int(ExpandIntegrand((a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule326)

    pattern327 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule327 = ReplacementRule(pattern327, lambda p, x, b, c, a, d, e, n : Int((a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule327)

    pattern328 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule328 = ReplacementRule(pattern328, lambda p, x, b, c, a, d, e, n : Int((a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule328)

    pattern329 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule329 = ReplacementRule(pattern329, lambda x, b, m, c, a, d, e, n : -d*Int(x**(m + S(-2))*(a + b*atanh(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*atanh(c*x))**n, x)/e)
    rubi.add(rule329)

    pattern330 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule330 = ReplacementRule(pattern330, lambda x, b, m, c, a, d, e, n : -d*Int(x**(m + S(-2))*(a + b*acoth(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*acoth(c*x))**n, x)/e)
    rubi.add(rule330)

    pattern331 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule331 = ReplacementRule(pattern331, lambda x, b, m, c, a, d, e, n : -e*Int(x**(m + S(2))*(a + b*atanh(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*atanh(c*x))**n, x)/d)
    rubi.add(rule331)

    pattern332 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule332 = ReplacementRule(pattern332, lambda x, b, m, c, a, d, e, n : -e*Int(x**(m + S(2))*(a + b*acoth(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*acoth(c*x))**n, x)/d)
    rubi.add(rule332)

    pattern333 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule333 = ReplacementRule(pattern333, lambda x, b, c, a, d, e, n : Int((a + b*atanh(c*x))**n/(-c*x + S(1)), x)/(c*d) + (a + b*atanh(c*x))**(n + S(1))/(b*e*(n + S(1))))
    rubi.add(rule333)

    pattern334 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule334 = ReplacementRule(pattern334, lambda x, b, c, a, d, e, n : Int((a + b*acoth(c*x))**n/(-c*x + S(1)), x)/(c*d) + (a + b*acoth(c*x))**(n + S(1))/(b*e*(n + S(1))))
    rubi.add(rule334)

    pattern335 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule335 = ReplacementRule(pattern335, lambda x, b, c, a, d, e, n : x*(a + b*atanh(c*x))**(n + S(1))/(b*c*d*(n + S(1))) - Int((a + b*atanh(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))))
    rubi.add(rule335)

    pattern336 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule336 = ReplacementRule(pattern336, lambda x, b, c, a, d, e, n : -x*(a + b*acoth(c*x))**(n + S(1))/(b*c*d*(n + S(1))) - Int((a + b*acoth(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))))
    rubi.add(rule336)

    pattern337 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule337 = ReplacementRule(pattern337, lambda x, b, m, c, a, d, e, n : -d*Int(x**(m + S(-2))*(a + b*atanh(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*atanh(c*x))**n, x)/e)
    rubi.add(rule337)

    pattern338 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule338 = ReplacementRule(pattern338, lambda x, b, m, c, a, d, e, n : -d*Int(x**(m + S(-2))*(a + b*acoth(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*acoth(c*x))**n, x)/e)
    rubi.add(rule338)

    pattern339 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule339 = ReplacementRule(pattern339, lambda x, b, c, a, d, e, n : Int((a + b*atanh(c*x))**n/(x*(c*x + S(1))), x)/d + (a + b*atanh(c*x))**(n + S(1))/(b*d*(n + S(1))))
    rubi.add(rule339)

    pattern340 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule340 = ReplacementRule(pattern340, lambda x, b, c, a, d, e, n : Int((a + b*acoth(c*x))**n/(x*(c*x + S(1))), x)/d + (a + b*acoth(c*x))**(n + S(1))/(b*d*(n + S(1))))
    rubi.add(rule340)

    pattern341 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule341 = ReplacementRule(pattern341, lambda x, b, m, c, a, d, e, n : -e*Int(x**(m + S(2))*(a + b*atanh(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*atanh(c*x))**n, x)/d)
    rubi.add(rule341)

    pattern342 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule342 = ReplacementRule(pattern342, lambda x, b, m, c, a, d, e, n : -e*Int(x**(m + S(2))*(a + b*acoth(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*acoth(c*x))**n, x)/d)
    rubi.add(rule342)

    pattern343 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule343 = ReplacementRule(pattern343, lambda x, b, m, c, a, d, e, n : -m*Int(x**(m + S(-1))*(a + b*atanh(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))) + x**m*(a + b*atanh(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule343)

    pattern344 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule344 = ReplacementRule(pattern344, lambda x, b, m, c, a, d, e, n : -m*Int(x**(m + S(-1))*(a + b*acoth(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))) + x**m*(a + b*acoth(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule344)

    pattern345 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda a, m: Not(NonzeroQ(a) & Equal(m, S(1)))))
    rule345 = ReplacementRule(pattern345, lambda x, b, m, c, a, d, e : Int(ExpandIntegrand(a + b*atanh(c*x), x**m/(d + e*x**S(2)), x), x))
    rubi.add(rule345)

    pattern346 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda a, m: Not(NonzeroQ(a) & Equal(m, S(1)))))
    rule346 = ReplacementRule(pattern346, lambda x, b, m, c, a, d, e : Int(ExpandIntegrand(a + b*acoth(c*x), x**m/(d + e*x**S(2)), x), x))
    rubi.add(rule346)

    pattern347 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule347 = ReplacementRule(pattern347, lambda p, x, b, c, a, d, e, n : b*n*Int((a + b*atanh(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(S(2)*c*(p + S(1))) + (a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule347)

    pattern348 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule348 = ReplacementRule(pattern348, lambda p, x, b, c, a, d, e, n : b*n*Int((a + b*acoth(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(S(2)*c*(p + S(1))) + (a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule348)

    pattern349 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule349 = ReplacementRule(pattern349, lambda x, b, c, a, d, e, n : x*(a + b*atanh(c*x))**(n + S(1))/(b*c*d*(d + e*x**S(2))*(n + S(1))) + S(4)*Int(x*(a + b*atanh(c*x))**(n + S(2))/(d + e*x**S(2))**S(2), x)/(b**S(2)*(n + S(1))*(n + S(2))) + (a + b*atanh(c*x))**(n + S(2))*(c**S(2)*x**S(2) + S(1))/(b**S(2)*e*(d + e*x**S(2))*(n + S(1))*(n + S(2))))
    rubi.add(rule349)

    pattern350 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule350 = ReplacementRule(pattern350, lambda x, b, c, a, d, e, n : x*(a + b*acoth(c*x))**(n + S(1))/(b*c*d*(d + e*x**S(2))*(n + S(1))) + S(4)*Int(x*(a + b*acoth(c*x))**(n + S(2))/(d + e*x**S(2))**S(2), x)/(b**S(2)*(n + S(1))*(n + S(2))) + (a + b*acoth(c*x))**(n + S(2))*(c**S(2)*x**S(2) + S(1))/(b**S(2)*e*(d + e*x**S(2))*(n + S(1))*(n + S(2))))
    rubi.add(rule350)

    pattern351 = Pattern(Integral(x_**S(2)*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-5)/2)))
    rule351 = ReplacementRule(pattern351, lambda p, x, b, c, a, d, e : -b*(d + e*x**S(2))**(p + S(1))/(S(4)*c**S(3)*d*(p + S(1))**S(2)) - x*(a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*c**S(2)*d*(p + S(1))) + Int((a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*c**S(2)*d*(p + S(1))))
    rubi.add(rule351)

    pattern352 = Pattern(Integral(x_**S(2)*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-5)/2)))
    rule352 = ReplacementRule(pattern352, lambda p, x, b, c, a, d, e : -b*(d + e*x**S(2))**(p + S(1))/(S(4)*c**S(3)*d*(p + S(1))**S(2)) - x*(a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*c**S(2)*d*(p + S(1))) + Int((a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*c**S(2)*d*(p + S(1))))
    rubi.add(rule352)

    pattern353 = Pattern(Integral(x_**S(2)*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule353 = ReplacementRule(pattern353, lambda x, b, c, a, d, e, n : -b*n*Int(x*(a + b*atanh(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/(S(2)*c) + x*(a + b*atanh(c*x))**n/(S(2)*c**S(2)*d*(d + e*x**S(2))) - (a + b*atanh(c*x))**(n + S(1))/(S(2)*b*c**S(3)*d**S(2)*(n + S(1))))
    rubi.add(rule353)

    pattern354 = Pattern(Integral(x_**S(2)*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule354 = ReplacementRule(pattern354, lambda x, b, c, a, d, e, n : -b*n*Int(x*(a + b*acoth(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/(S(2)*c) + x*(a + b*acoth(c*x))**n/(S(2)*c**S(2)*d*(d + e*x**S(2))) - (a + b*acoth(c*x))**(n + S(1))/(S(2)*b*c**S(3)*d**S(2)*(n + S(1))))
    rubi.add(rule354)

    pattern355 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule355 = ReplacementRule(pattern355, lambda p, x, b, m, c, a, d, e : -b*x**m*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) + x**(m + S(-1))*(a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) - (m + S(-1))*Int(x**(m + S(-2))*(a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule355)

    pattern356 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule356 = ReplacementRule(pattern356, lambda p, x, b, m, c, a, d, e : -b*x**m*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) + x**(m + S(-1))*(a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) - (m + S(-1))*Int(x**(m + S(-2))*(a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule356)

    pattern357 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule357 = ReplacementRule(pattern357, lambda p, x, b, m, c, a, d, e, n : b**S(2)*n*(n + S(-1))*Int(x**m*(a + b*atanh(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/m**S(2) - b*n*x**m*(a + b*atanh(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) + x**(m + S(-1))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) - (m + S(-1))*Int(x**(m + S(-2))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule357)

    pattern358 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p, n: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule358 = ReplacementRule(pattern358, lambda p, x, b, m, c, a, d, e, n : b**S(2)*n*(n + S(-1))*Int(x**m*(a + b*acoth(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/m**S(2) - b*n*x**m*(a + b*acoth(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) + x**(m + S(-1))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) - (m + S(-1))*Int(x**(m + S(-2))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule358)

    pattern359 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule359 = ReplacementRule(pattern359, lambda p, x, b, m, c, a, d, e, n : -m*Int(x**(m + S(-1))*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) + x**m*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule359)

    pattern360 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule360 = ReplacementRule(pattern360, lambda p, x, b, m, c, a, d, e, n : -m*Int(x**(m + S(-1))*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) + x**m*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule360)

    pattern361 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule361 = ReplacementRule(pattern361, lambda p, x, b, m, c, a, d, e, n : -b*c*n*Int(x**(m + S(1))*(a + b*atanh(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(m + S(1)) + x**(m + S(1))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*(m + S(1))))
    rubi.add(rule361)

    pattern362 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule362 = ReplacementRule(pattern362, lambda p, x, b, m, c, a, d, e, n : -b*c*n*Int(x**(m + S(1))*(a + b*acoth(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(m + S(1)) + x**(m + S(1))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*(m + S(1))))
    rubi.add(rule362)

    pattern363 = Pattern(Integral(x_**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: NonzeroQ(m + S(2))))
    rule363 = ReplacementRule(pattern363, lambda x, b, m, c, a, d, e : -b*c*d*Int(x**(m + S(1))/sqrt(d + e*x**S(2)), x)/(m + S(2)) + d*Int(x**m*(a + b*atanh(c*x))/sqrt(d + e*x**S(2)), x)/(m + S(2)) + x**(m + S(1))*(a + b*atanh(c*x))*sqrt(d + e*x**S(2))/(m + S(2)))
    rubi.add(rule363)

    pattern364 = Pattern(Integral(x_**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: NonzeroQ(m + S(2))))
    rule364 = ReplacementRule(pattern364, lambda x, b, m, c, a, d, e : -b*c*d*Int(x**(m + S(1))/sqrt(d + e*x**S(2)), x)/(m + S(2)) + d*Int(x**m*(a + b*acoth(c*x))/sqrt(d + e*x**S(2)), x)/(m + S(2)) + x**(m + S(1))*(a + b*acoth(c*x))*sqrt(d + e*x**S(2))/(m + S(2)))
    rubi.add(rule364)

    pattern365 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p: Greater(p, S(1))))
    rule365 = ReplacementRule(pattern365, lambda p, x, b, m, c, a, d, e, n : Int(ExpandIntegrand(x**m*(a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule365)

    pattern366 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p: Greater(p, S(1))))
    rule366 = ReplacementRule(pattern366, lambda p, x, b, m, c, a, d, e, n : Int(ExpandIntegrand(x**m*(a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule366)

    pattern367 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, p, m: RationalQ(m) | (IntegerQ(p) & EqQ(n, S(1)))))
    rule367 = ReplacementRule(pattern367, lambda p, x, b, m, c, a, d, e, n : -c**S(2)*d*Int(x**(m + S(2))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x) + d*Int(x**m*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x))
    rubi.add(rule367)

    pattern368 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, p, m: RationalQ(m) | (IntegerQ(p) & EqQ(n, S(1)))))
    rule368 = ReplacementRule(pattern368, lambda p, x, b, m, c, a, d, e, n : -c**S(2)*d*Int(x**(m + S(2))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x) + d*Int(x**m*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x))
    rubi.add(rule368)

    pattern369 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule369 = ReplacementRule(pattern369, lambda x, b, m, c, a, d, e, n : b*n*Int(x**(m + S(-1))*(a + b*atanh(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x)/(c*m) + (m + S(-1))*Int(x**(m + S(-2))*(a + b*atanh(c*x))**n/sqrt(d + e*x**S(2)), x)/(c**S(2)*m) - x**(m + S(-1))*(a + b*atanh(c*x))**n*sqrt(d + e*x**S(2))/(c**S(2)*d*m))
    rubi.add(rule369)

    pattern370 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule370 = ReplacementRule(pattern370, lambda x, b, m, c, a, d, e, n : b*n*Int(x**(m + S(-1))*(a + b*acoth(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x)/(c*m) + (m + S(-1))*Int(x**(m + S(-2))*(a + b*acoth(c*x))**n/sqrt(d + e*x**S(2)), x)/(c**S(2)*m) - x**(m + S(-1))*(a + b*acoth(c*x))**n*sqrt(d + e*x**S(2))/(c**S(2)*d*m))
    rubi.add(rule370)

    pattern371 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule371 = ReplacementRule(pattern371, lambda x, b, c, a, d, e : b*PolyLog(S(2), -sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d) - b*PolyLog(S(2), sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d) - S(2)*(a + b*atanh(c*x))*atanh(sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d))
    rubi.add(rule371)

    pattern372 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule372 = ReplacementRule(pattern372, lambda x, b, c, a, d, e : b*PolyLog(S(2), -sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d) - b*PolyLog(S(2), sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d) - S(2)*(a + b*acoth(c*x))*atanh(sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d))
    rubi.add(rule372)

    pattern373 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule373 = ReplacementRule(pattern373, lambda x, b, c, a, d, e, n : Subst(Int((a + b*x)**n*csch(x), x), x, atanh(c*x))/sqrt(d))
    rubi.add(rule373)

    pattern374 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule374 = ReplacementRule(pattern374, lambda x, b, c, a, d, e, n : -c*x*sqrt(S(1) - S(1)/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*sech(x), x), x, acoth(c*x))/sqrt(d + e*x**S(2)))
    rubi.add(rule374)

    pattern375 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule375 = ReplacementRule(pattern375, lambda x, b, c, a, d, e, n : sqrt(-c**S(2)*x**S(2) + S(1))*Int((a + b*atanh(c*x))**n/(x*sqrt(-c**S(2)*x**S(2) + S(1))), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule375)

    pattern376 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule376 = ReplacementRule(pattern376, lambda x, b, c, a, d, e, n : sqrt(-c**S(2)*x**S(2) + S(1))*Int((a + b*acoth(c*x))**n/(x*sqrt(-c**S(2)*x**S(2) + S(1))), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule376)

    pattern377 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(x_**S(2)*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule377 = ReplacementRule(pattern377, lambda x, b, c, a, d, e, n : b*c*n*Int((a + b*atanh(c*x))**(n + S(-1))/(x*sqrt(d + e*x**S(2))), x) - (a + b*atanh(c*x))**n*sqrt(d + e*x**S(2))/(d*x))
    rubi.add(rule377)

    pattern378 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(x_**S(2)*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule378 = ReplacementRule(pattern378, lambda x, b, c, a, d, e, n : b*c*n*Int((a + b*acoth(c*x))**(n + S(-1))/(x*sqrt(d + e*x**S(2))), x) - (a + b*acoth(c*x))**n*sqrt(d + e*x**S(2))/(d*x))
    rubi.add(rule378)

    pattern379 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: Unequal(m, S(-2))))
    rule379 = ReplacementRule(pattern379, lambda x, b, m, c, a, d, e, n : -b*c*n*Int(x**(m + S(1))*(a + b*atanh(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x)/(m + S(1)) + c**S(2)*(m + S(2))*Int(x**(m + S(2))*(a + b*atanh(c*x))**n/sqrt(d + e*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*atanh(c*x))**n*sqrt(d + e*x**S(2))/(d*(m + S(1))))
    rubi.add(rule379)

    pattern380 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: Unequal(m, S(-2))))
    rule380 = ReplacementRule(pattern380, lambda x, b, m, c, a, d, e, n : -b*c*n*Int(x**(m + S(1))*(a + b*acoth(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x)/(m + S(1)) + c**S(2)*(m + S(2))*Int(x**(m + S(2))*(a + b*acoth(c*x))**n/sqrt(d + e*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*acoth(c*x))**n*sqrt(d + e*x**S(2))/(d*(m + S(1))))
    rubi.add(rule380)

    pattern381 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule381 = ReplacementRule(pattern381, lambda p, x, b, m, c, a, d, e, n : -d*Int(x**(m + S(-2))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x)/e + Int(x**(m + S(-2))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/e)
    rubi.add(rule381)

    pattern382 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule382 = ReplacementRule(pattern382, lambda p, x, b, m, c, a, d, e, n : -d*Int(x**(m + S(-2))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x)/e + Int(x**(m + S(-2))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/e)
    rubi.add(rule382)

    pattern383 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Less(m, S(0))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule383 = ReplacementRule(pattern383, lambda p, x, b, m, c, a, d, e, n : -e*Int(x**(m + S(2))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x)/d + Int(x**m*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/d)
    rubi.add(rule383)

    pattern384 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Less(m, S(0))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule384 = ReplacementRule(pattern384, lambda p, x, b, m, c, a, d, e, n : -e*Int(x**(m + S(2))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x)/d + Int(x**m*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/d)
    rubi.add(rule384)

    pattern385 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(2))))
    rule385 = ReplacementRule(pattern385, lambda p, x, b, m, c, a, d, e, n : c*(m + S(2)*p + S(2))*Int(x**(m + S(1))*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) - m*Int(x**(m + S(-1))*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) + x**m*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule385)

    pattern386 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(2))))
    rule386 = ReplacementRule(pattern386, lambda p, x, b, m, c, a, d, e, n : c*(m + S(2)*p + S(2))*Int(x**(m + S(1))*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) - m*Int(x**(m + S(-1))*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) + x**m*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule386)

    pattern387 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda d, p: IntegerQ(p) | PositiveQ(d)))
    rule387 = ReplacementRule(pattern387, lambda p, x, b, m, c, a, d, e, n : c**(-m + S(-1))*d**p*Subst(Int((a + b*x)**n*Cosh(x)**(-m - S(2)*p + S(-2))*sinh(x)**m, x), x, atanh(c*x)))
    rubi.add(rule387)

    pattern388 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda d, p: Not(IntegerQ(p) | PositiveQ(d))))
    rule388 = ReplacementRule(pattern388, lambda p, x, b, m, c, a, d, e, n : d**(p + S(1)/2)*sqrt(-c**S(2)*x**S(2) + S(1))*Int(x**m*(a + b*atanh(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x)/sqrt(d + e*x**S(2)))
    rubi.add(rule388)

    pattern389 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda p: IntegerQ(p)))
    rule389 = ReplacementRule(pattern389, lambda p, x, b, m, c, a, d, e, n : -c**(-m + S(-1))*(-d)**p*Subst(Int((a + b*x)**n*Cosh(x)**m*sinh(x)**(-m - S(2)*p + S(-2)), x), x, acoth(c*x)))
    rubi.add(rule389)

    pattern390 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule390 = ReplacementRule(pattern390, lambda p, x, b, m, c, a, d, e, n : -c**(-m)*x*(-d)**(p + S(1)/2)*sqrt((c**S(2)*x**S(2) + S(-1))/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*Cosh(x)**m*sinh(x)**(-m - S(2)*p + S(-2)), x), x, acoth(c*x))/sqrt(d + e*x**S(2)))
    rubi.add(rule390)

    pattern391 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule391 = ReplacementRule(pattern391, lambda p, x, b, c, a, d, e : -b*c*Int((d + e*x**S(2))**(p + S(1))/(-c**S(2)*x**S(2) + S(1)), x)/(S(2)*e*(p + S(1))) + (a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule391)

    pattern392 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule392 = ReplacementRule(pattern392, lambda p, x, b, c, a, d, e : -b*c*Int((d + e*x**S(2))**(p + S(1))/(-c**S(2)*x**S(2) + S(1)), x)/(S(2)*e*(p + S(1))) + (a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule392)

    pattern393 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & Not(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))), )
    def With393(p, x, b, m, c, a, d, e):
        u = IntHide(x**m*(d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*atanh(c*x), u, x)
    rule393 = ReplacementRule(pattern393, lambda p, x, b, m, c, a, d, e : With393(p, x, b, m, c, a, d, e))
    rubi.add(rule393)

    pattern394 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & Not(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))), )
    def With394(p, x, b, m, c, a, d, e):
        u = IntHide(x**m*(d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acoth(c*x), u, x)
    rule394 = ReplacementRule(pattern394, lambda p, x, b, m, c, a, d, e : With394(p, x, b, m, c, a, d, e))
    rubi.add(rule394)

    pattern395 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (IntegerQ(m) & Less(p, S(-1)) & Unequal(m, S(1)))))
    rule395 = ReplacementRule(pattern395, lambda p, x, b, m, c, a, d, e, n : Int(ExpandIntegrand((a + b*atanh(c*x))**n, x**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule395)

    pattern396 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (IntegerQ(m) & Less(p, S(-1)) & Unequal(m, S(1)))))
    rule396 = ReplacementRule(pattern396, lambda p, x, b, m, c, a, d, e, n : Int(ExpandIntegrand((a + b*acoth(c*x))**n, x**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule396)

    pattern397 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*atanh(x_*WC('c', S(1))))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule397 = ReplacementRule(pattern397, lambda p, x, b, m, c, a, d, e : a*Int(x**m*(d + e*x**S(2))**p, x) + b*Int(x**m*(d + e*x**S(2))**p*atanh(c*x), x))
    rubi.add(rule397)

    pattern398 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule398 = ReplacementRule(pattern398, lambda p, x, b, m, c, a, d, e : a*Int(x**m*(d + e*x**S(2))**p, x) + b*Int(x**m*(d + e*x**S(2))**p*acoth(c*x), x))
    rubi.add(rule398)

    pattern399 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule399 = ReplacementRule(pattern399, lambda p, x, b, m, c, a, d, e, n : Int(x**m*(a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule399)

    pattern400 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule400 = ReplacementRule(pattern400, lambda p, x, b, m, c, a, d, e, n : Int(x**m*(a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule400)

    pattern401 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*atanh(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(u**S(2) - (S(1) - S(2)/(c*x + S(1)))**S(2))))
    rule401 = ReplacementRule(pattern401, lambda x, b, c, a, d, u, e, n : -Int((a + b*atanh(c*x))**n*log(-u + S(1))/(d + e*x**S(2)), x)/S(2) + Int((a + b*atanh(c*x))**n*log(u + S(1))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule401)

    pattern402 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*acoth(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(u**S(2) - (S(1) - S(2)/(c*x + S(1)))**S(2))))
    rule402 = ReplacementRule(pattern402, lambda x, b, c, a, d, u, e, n : -Int((a + b*acoth(c*x))**n*log(SimplifyIntegrand(S(1) - S(1)/u, x))/(d + e*x**S(2)), x)/S(2) + Int((a + b*acoth(c*x))**n*log(SimplifyIntegrand(S(1) + 1/u, x))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule402)

    pattern403 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*atanh(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(u**S(2) - (S(1) - S(2)/(-c*x + S(1)))**S(2))))
    rule403 = ReplacementRule(pattern403, lambda x, b, c, a, d, u, e, n : -Int((a + b*atanh(c*x))**n*log(-u + S(1))/(d + e*x**S(2)), x)/S(2) + Int((a + b*atanh(c*x))**n*log(u + S(1))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule403)

    pattern404 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*acoth(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(u**S(2) - (S(1) - S(2)/(-c*x + S(1)))**S(2))))
    rule404 = ReplacementRule(pattern404, lambda x, b, c, a, d, u, e, n : -Int((a + b*acoth(c*x))**n*log(SimplifyIntegrand(S(1) - S(1)/u, x))/(d + e*x**S(2)), x)/S(2) + Int((a + b*acoth(c*x))**n*log(SimplifyIntegrand(S(1) + 1/u, x))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule404)

    pattern405 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(-(S(1) - S(2)/(c*x + S(1)))**S(2) + (-u + S(1))**S(2))))
    rule405 = ReplacementRule(pattern405, lambda x, b, c, a, d, u, e, n : -b*n*Int((a + b*atanh(c*x))**(n + S(-1))*PolyLog(S(2), -u + S(1))/(d + e*x**S(2)), x)/S(2) + (a + b*atanh(c*x))**n*PolyLog(S(2), -u + S(1))/(S(2)*c*d))
    rubi.add(rule405)

    pattern406 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(-(S(1) - S(2)/(c*x + S(1)))**S(2) + (-u + S(1))**S(2))))
    rule406 = ReplacementRule(pattern406, lambda x, b, c, a, d, u, e, n : -b*n*Int((a + b*acoth(c*x))**(n + S(-1))*PolyLog(S(2), -u + S(1))/(d + e*x**S(2)), x)/S(2) + (a + b*acoth(c*x))**n*PolyLog(S(2), -u + S(1))/(S(2)*c*d))
    rubi.add(rule406)

    pattern407 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(-(S(1) - S(2)/(-c*x + S(1)))**S(2) + (-u + S(1))**S(2))))
    rule407 = ReplacementRule(pattern407, lambda x, b, c, a, d, u, e, n : b*n*Int((a + b*atanh(c*x))**(n + S(-1))*PolyLog(S(2), -u + S(1))/(d + e*x**S(2)), x)/S(2) - (a + b*atanh(c*x))**n*PolyLog(S(2), -u + S(1))/(S(2)*c*d))
    rubi.add(rule407)

    pattern408 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(-(S(1) - S(2)/(-c*x + S(1)))**S(2) + (-u + S(1))**S(2))))
    rule408 = ReplacementRule(pattern408, lambda x, b, c, a, d, u, e, n : b*n*Int((a + b*acoth(c*x))**(n + S(-1))*PolyLog(S(2), -u + S(1))/(d + e*x**S(2)), x)/S(2) - (a + b*acoth(c*x))**n*PolyLog(S(2), -u + S(1))/(S(2)*c*d))
    rubi.add(rule408)

    pattern409 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(u**S(2) - (S(1) - S(2)/(c*x + S(1)))**S(2))))
    rule409 = ReplacementRule(pattern409, lambda p, x, b, c, a, d, u, e, n : b*n*Int((a + b*atanh(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) - (a + b*atanh(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule409)

    pattern410 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(u**S(2) - (S(1) - S(2)/(c*x + S(1)))**S(2))))
    rule410 = ReplacementRule(pattern410, lambda p, x, b, c, a, d, u, e, n : b*n*Int((a + b*acoth(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) - (a + b*acoth(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule410)

    pattern411 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(u**S(2) - (S(1) - S(2)/(-c*x + S(1)))**S(2))))
    rule411 = ReplacementRule(pattern411, lambda p, x, b, c, a, d, u, e, n : -b*n*Int((a + b*atanh(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) + (a + b*atanh(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule411)

    pattern412 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda c, u, x: ZeroQ(u**S(2) - (S(1) - S(2)/(-c*x + S(1)))**S(2))))
    rule412 = ReplacementRule(pattern412, lambda p, x, b, c, a, d, u, e, n : -b*n*Int((a + b*acoth(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) + (a + b*acoth(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule412)

    pattern413 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)))
    rule413 = ReplacementRule(pattern413, lambda x, b, c, a, d, e : (-log(a + b*acoth(c*x)) + log(a + b*atanh(c*x)))/(b**S(2)*c*d*(acoth(c*x) - atanh(c*x))))
    rubi.add(rule413)

    pattern414 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Inequality(S(0), Less, n, LessEqual, m)))
    rule414 = ReplacementRule(pattern414, lambda x, b, m, c, a, d, e, n : -n*Int((a + b*acoth(c*x))**(m + S(1))*(a + b*atanh(c*x))**(n + S(-1))/(d + e*x**S(2)), x)/(m + S(1)) + (a + b*acoth(c*x))**(m + S(1))*(a + b*atanh(c*x))**n/(b*c*d*(m + S(1))))
    rubi.add(rule414)

    pattern415 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('m', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Less(S(0), n, m)))
    rule415 = ReplacementRule(pattern415, lambda x, b, m, c, a, d, e, n : -n*Int((a + b*acoth(c*x))**(n + S(-1))*(a + b*atanh(c*x))**(m + S(1))/(d + e*x**S(2)), x)/(m + S(1)) + (a + b*acoth(c*x))**n*(a + b*atanh(c*x))**(m + S(1))/(b*c*d*(m + S(1))))
    rubi.add(rule415)

    pattern416 = Pattern(Integral(atanh(x_*WC('a', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda c, a, d, n: Not(Equal(n, S(2)) & ZeroQ(a**S(2)*c + d))))
    rule416 = ReplacementRule(pattern416, lambda x, c, a, d, n : -Int(log(-a*x + S(1))/(c + d*x**n), x)/S(2) + Int(log(a*x + S(1))/(c + d*x**n), x)/S(2))
    rubi.add(rule416)

    pattern417 = Pattern(Integral(acoth(x_*WC('a', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda c, a, d, n: Not(Equal(n, S(2)) & ZeroQ(a**S(2)*c + d))))
    rule417 = ReplacementRule(pattern417, lambda x, c, a, d, n : -Int(log(S(1) - S(1)/(a*x))/(c + d*x**n), x)/S(2) + Int(log(S(1) + S(1)/(a*x))/(c + d*x**n), x)/S(2))
    rubi.add(rule417)

    pattern418 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)))
    rule418 = ReplacementRule(pattern418, lambda x, b, g, c, f, a, d, e : -b*c*Int(x*(d + e*log(f + g*x**S(2)))/(-c**S(2)*x**S(2) + S(1)), x) - S(2)*e*g*Int(x**S(2)*(a + b*atanh(c*x))/(f + g*x**S(2)), x) + x*(a + b*atanh(c*x))*(d + e*log(f + g*x**S(2))))
    rubi.add(rule418)

    pattern419 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)))
    rule419 = ReplacementRule(pattern419, lambda x, b, g, c, f, a, d, e : -b*c*Int(x*(d + e*log(f + g*x**S(2)))/(-c**S(2)*x**S(2) + S(1)), x) - S(2)*e*g*Int(x**S(2)*(a + b*acoth(c*x))/(f + g*x**S(2)), x) + x*(a + b*acoth(c*x))*(d + e*log(f + g*x**S(2))))
    rubi.add(rule419)

    pattern420 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2))))
    rule420 = ReplacementRule(pattern420, lambda x, b, m, g, c, f, a, d, e : -b*c*Int(x**(m + S(1))*(d + e*log(f + g*x**S(2)))/(-c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) - S(2)*e*g*Int(x**(m + S(2))*(a + b*atanh(c*x))/(f + g*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*atanh(c*x))*(d + e*log(f + g*x**S(2)))/(m + S(1)))
    rubi.add(rule420)

    pattern421 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2))))
    rule421 = ReplacementRule(pattern421, lambda x, b, m, g, c, f, a, d, e : -b*c*Int(x**(m + S(1))*(d + e*log(f + g*x**S(2)))/(-c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) - S(2)*e*g*Int(x**(m + S(2))*(a + b*acoth(c*x))/(f + g*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*acoth(c*x))*(d + e*log(f + g*x**S(2)))/(m + S(1)))
    rubi.add(rule421)

    pattern422 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m/S(2) + S(1)/2)), )
    def With422(x, b, m, g, c, f, a, d, e):
        u = IntHide(x**m*(d + e*log(f + g*x**S(2))), x)
        return -b*c*Int(ExpandIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*atanh(c*x), u, x)
    rule422 = ReplacementRule(pattern422, lambda x, b, m, g, c, f, a, d, e : With422(x, b, m, g, c, f, a, d, e))
    rubi.add(rule422)

    pattern423 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m/S(2) + S(1)/2)), )
    def With423(x, b, m, g, c, f, a, d, e):
        u = IntHide(x**m*(d + e*log(f + g*x**S(2))), x)
        return -b*c*Int(ExpandIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acoth(c*x), u, x)
    rule423 = ReplacementRule(pattern423, lambda x, b, m, g, c, f, a, d, e : With423(x, b, m, g, c, f, a, d, e))
    rubi.add(rule423)

    pattern424 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Unequal(m, S(-1))), )
    def With424(x, b, m, g, c, f, a, d, e):
        u = IntHide(x**m*(a + b*atanh(c*x)), x)
        return -S(2)*e*g*Int(ExpandIntegrand(u*x/(f + g*x**S(2)), x), x) + Dist(d + e*log(f + g*x**S(2)), u, x)
    rule424 = ReplacementRule(pattern424, lambda x, b, m, g, c, f, a, d, e : With424(x, b, m, g, c, f, a, d, e))
    rubi.add(rule424)

    pattern425 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Unequal(m, S(-1))), )
    def With425(x, b, m, g, c, f, a, d, e):
        u = IntHide(x**m*(a + b*acoth(c*x)), x)
        return -S(2)*e*g*Int(ExpandIntegrand(u*x/(f + g*x**S(2)), x), x) + Dist(d + e*log(f + g*x**S(2)), u, x)
    rule425 = ReplacementRule(pattern425, lambda x, b, m, g, c, f, a, d, e : With425(x, b, m, g, c, f, a, d, e))
    rubi.add(rule425)

    pattern426 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**S(2)*(WC('d', S(0)) + WC('e', S(1))*log(f_ + x_**S(2)*WC('g', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda g, c, f: ZeroQ(c**S(2)*f + g)))
    rule426 = ReplacementRule(pattern426, lambda x, b, g, c, f, a, d, e : b*c*e*Int(x**S(2)*(a + b*atanh(c*x))/(-c**S(2)*x**S(2) + S(1)), x) + b*Int((a + b*atanh(c*x))*(d + e*log(f + g*x**S(2))), x)/c - e*x**S(2)*(a + b*atanh(c*x))**S(2)/S(2) + (a + b*atanh(c*x))**S(2)*(d + e*log(f + g*x**S(2)))*(f + g*x**S(2))/(S(2)*g))
    rubi.add(rule426)

    pattern427 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**S(2)*(WC('d', S(0)) + WC('e', S(1))*log(f_ + x_**S(2)*WC('g', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda g, c, f: ZeroQ(c**S(2)*f + g)))
    rule427 = ReplacementRule(pattern427, lambda x, b, g, c, f, a, d, e : b*c*e*Int(x**S(2)*(a + b*acoth(c*x))/(-c**S(2)*x**S(2) + S(1)), x) + b*Int((a + b*acoth(c*x))*(d + e*log(f + g*x**S(2))), x)/c - e*x**S(2)*(a + b*acoth(c*x))**S(2)/S(2) + (a + b*acoth(c*x))**S(2)*(d + e*log(f + g*x**S(2)))*(f + g*x**S(2))/(S(2)*g))
    rubi.add(rule427)

    pattern428 = Pattern(Integral(exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: OddQ(n)))
    rule428 = ReplacementRule(pattern428, lambda x, a, n : Int((-a*x + S(1))**(-n/S(2) + S(1)/2)*(a*x + S(1))**(n/S(2) + S(1)/2)/sqrt(-a**S(2)*x**S(2) + S(1)), x))
    rubi.add(rule428)

    pattern429 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: OddQ(n)))
    rule429 = ReplacementRule(pattern429, lambda x, a, n, m : Int(x**m*(-a*x + S(1))**(-n/S(2) + S(1)/2)*(a*x + S(1))**(n/S(2) + S(1)/2)/sqrt(-a**S(2)*x**S(2) + S(1)), x))
    rubi.add(rule429)

    pattern430 = Pattern(Integral(exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(OddQ(n))))
    rule430 = ReplacementRule(pattern430, lambda x, a, n : Int((-a*x + S(1))**(-n/S(2))*(a*x + S(1))**(n/S(2)), x))
    rubi.add(rule430)

    pattern431 = Pattern(Integral(x_**WC('m', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(OddQ(n))))
    rule431 = ReplacementRule(pattern431, lambda x, a, n, m : Int(x**m*(-a*x + S(1))**(-n/S(2))*(a*x + S(1))**(n/S(2)), x))
    rubi.add(rule431)

    pattern432 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a*c + d)), CustomConstraint(lambda n: IntegerQ(n/S(2) + S(-1)/2)), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule432 = ReplacementRule(pattern432, lambda p, x, c, a, d, n : c**n*Int((c + d*x)**(-n + p)*(-a**S(2)*x**S(2) + S(1))**(n/S(2)), x))
    rubi.add(rule432)

    pattern433 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a*c + d)), CustomConstraint(lambda n: IntegerQ(n/S(2) + S(-1)/2)), CustomConstraint(lambda p, n: IntegerQ(p) | ZeroQ(-n/S(2) + p) | ZeroQ(-n/S(2) + p + S(-1))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule433 = ReplacementRule(pattern433, lambda p, x, m, f, c, a, d, e, n : c**n*Int((c + d*x)**(-n + p)*(e + f*x)**m*(-a**S(2)*x**S(2) + S(1))**(n/S(2)), x))
    rubi.add(rule433)

    pattern434 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c**S(2) - d**S(2))), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)))
    rule434 = ReplacementRule(pattern434, lambda p, x, c, a, d, u, n : c**p*Int(u*(S(1) + d*x/c)**p*(-a*x + S(1))**(-n/S(2))*(a*x + S(1))**(n/S(2)), x))
    rubi.add(rule434)

    pattern435 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c**S(2) - d**S(2))), CustomConstraint(lambda c, p: Not(IntegerQ(p) | PositiveQ(c))))
    rule435 = ReplacementRule(pattern435, lambda p, x, c, a, d, u, n : Int(u*(c + d*x)**p*(-a*x + S(1))**(-n/S(2))*(a*x + S(1))**(n/S(2)), x))
    rubi.add(rule435)

    pattern436 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(-a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule436 = ReplacementRule(pattern436, lambda p, x, c, a, d, u, n : d**p*Int(u*x**(-p)*(c*x/d + S(1))**p*exp(n*atanh(a*x)), x))
    rubi.add(rule436)

    pattern437 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(-a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(n/S(2))), CustomConstraint(lambda c: PositiveQ(c)))
    rule437 = ReplacementRule(pattern437, lambda p, x, c, a, d, u, n : (S(-1))**(n/S(2))*c**p*Int(u*(S(1) - S(1)/(a*x))**(-n/S(2))*(S(1) + S(1)/(a*x))**(n/S(2))*(S(1) + d/(c*x))**p, x))
    rubi.add(rule437)

    pattern438 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(-a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(n/S(2))), CustomConstraint(lambda c: Not(PositiveQ(c))))
    rule438 = ReplacementRule(pattern438, lambda p, x, c, a, d, u, n : Int(u*(c + d/x)**p*(-a*x + S(1))**(-n/S(2))*(a*x + S(1))**(n/S(2)), x))
    rubi.add(rule438)

    pattern439 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(-a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule439 = ReplacementRule(pattern439, lambda p, x, c, a, d, u, n : x**p*(c + d/x)**p*(c*x/d + S(1))**(-p)*Int(u*x**(-p)*(c*x/d + S(1))**p*exp(n*atanh(a*x)), x))
    rubi.add(rule439)

    pattern440 = Pattern(Integral(exp(n_*atanh(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule440 = ReplacementRule(pattern440, lambda x, c, a, d, n : (-a*x + n)*exp(n*atanh(a*x))/(a*c*sqrt(c + d*x**S(2))*(n**S(2) + S(-1))))
    rubi.add(rule440)

    pattern441 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Not(IntegerQ(n))), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) - S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule441 = ReplacementRule(pattern441, lambda p, x, c, a, d, n : -S(2)*(p + S(1))*(S(2)*p + S(3))*Int((c + d*x**S(2))**(p + S(1))*exp(n*atanh(a*x)), x)/(c*(n**S(2) - S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*atanh(a*x))/(a*c*(n**S(2) - S(4)*(p + S(1))**S(2))))
    rubi.add(rule441)

    pattern442 = Pattern(Integral(exp(WC('n', S(1))*atanh(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))))
    rule442 = ReplacementRule(pattern442, lambda x, c, a, d, n : exp(n*atanh(a*x))/(a*c*n))
    rubi.add(rule442)

    pattern443 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2) + S(1)/2)), CustomConstraint(lambda p, n: Not(IntegerQ(-n/S(2) + p))))
    rule443 = ReplacementRule(pattern443, lambda p, x, c, a, d, n : c**p*Int((a*x + S(1))**n*(-a**S(2)*x**S(2) + S(1))**(-n/S(2) + p), x))
    rubi.add(rule443)

    pattern444 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n/S(2) + S(-1)/2)), CustomConstraint(lambda p, n: Not(IntegerQ(-n/S(2) + p))))
    rule444 = ReplacementRule(pattern444, lambda p, x, c, a, d, n : c**p*Int((-a*x + S(1))**(-n)*(-a**S(2)*x**S(2) + S(1))**(n/S(2) + p), x))
    rubi.add(rule444)

    pattern445 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)))
    rule445 = ReplacementRule(pattern445, lambda p, x, c, a, d, n : c**p*Int((-a*x + S(1))**(-n/S(2) + p)*(a*x + S(1))**(n/S(2) + p), x))
    rubi.add(rule445)

    pattern446 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2))))
    rule446 = ReplacementRule(pattern446, lambda p, x, c, a, d, n : c**(n/S(2))*Int((c + d*x**S(2))**(-n/S(2) + p)*(a*x + S(1))**n, x))
    rubi.add(rule446)

    pattern447 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: NegativeIntegerQ(n/S(2))))
    rule447 = ReplacementRule(pattern447, lambda p, x, c, a, d, n : c**(-n/S(2))*Int((c + d*x**S(2))**(n/S(2) + p)*(-a*x + S(1))**(-n), x))
    rubi.add(rule447)

    pattern448 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: Not(IntegerQ(p) | PositiveQ(c))))
    rule448 = ReplacementRule(pattern448, lambda p, x, c, a, d, n : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(-a**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((-a**S(2)*x**S(2) + S(1))**p*exp(n*atanh(a*x)), x))
    rubi.add(rule448)

    pattern449 = Pattern(Integral(x_*exp(n_*atanh(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule449 = ReplacementRule(pattern449, lambda x, c, a, d, n : (-a*n*x + S(1))*exp(n*atanh(a*x))/(d*sqrt(c + d*x**S(2))*(n**S(2) + S(-1))))
    rubi.add(rule449)

    pattern450 = Pattern(Integral(x_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Not(IntegerQ(n))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule450 = ReplacementRule(pattern450, lambda p, x, c, a, d, n : -a*c*n*Int((c + d*x**S(2))**p*exp(n*atanh(a*x)), x)/(S(2)*d*(p + S(1))) + (c + d*x**S(2))**(p + S(1))*exp(n*atanh(a*x))/(S(2)*d*(p + S(1))))
    rubi.add(rule450)

    pattern451 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda p, n: ZeroQ(n**S(2) + S(2)*p + S(2))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule451 = ReplacementRule(pattern451, lambda p, x, c, a, d, n : (c + d*x**S(2))**(p + S(1))*(-a*n*x + S(1))*exp(n*atanh(a*x))/(a*d*n*(n**S(2) + S(-1))))
    rubi.add(rule451)

    pattern452 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Not(IntegerQ(n))), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) - S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule452 = ReplacementRule(pattern452, lambda p, x, c, a, d, n : (n**S(2) + S(2)*p + S(2))*Int((c + d*x**S(2))**(p + S(1))*exp(n*atanh(a*x)), x)/(d*(n**S(2) - S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(-S(2)*a*x*(p + S(1)) - n)*exp(n*atanh(a*x))/(a*d*(n**S(2) - S(4)*(p + S(1))**S(2))))
    rubi.add(rule452)

    pattern453 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2) + S(1)/2)), CustomConstraint(lambda p, n: Not(IntegerQ(-n/S(2) + p))))
    rule453 = ReplacementRule(pattern453, lambda p, x, m, c, a, d, n : c**p*Int(x**m*(a*x + S(1))**n*(-a**S(2)*x**S(2) + S(1))**(-n/S(2) + p), x))
    rubi.add(rule453)

    pattern454 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda n: NegativeIntegerQ(n/S(2) + S(-1)/2)), CustomConstraint(lambda p, n: Not(IntegerQ(-n/S(2) + p))))
    rule454 = ReplacementRule(pattern454, lambda p, x, m, c, a, d, n : c**p*Int(x**m*(-a*x + S(1))**(-n)*(-a**S(2)*x**S(2) + S(1))**(n/S(2) + p), x))
    rubi.add(rule454)

    pattern455 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)))
    rule455 = ReplacementRule(pattern455, lambda p, x, m, c, a, d, n : c**p*Int(x**m*(-a*x + S(1))**(-n/S(2) + p)*(a*x + S(1))**(n/S(2) + p), x))
    rubi.add(rule455)

    pattern456 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2))))
    rule456 = ReplacementRule(pattern456, lambda p, x, m, c, a, d, n : c**(n/S(2))*Int(x**m*(c + d*x**S(2))**(-n/S(2) + p)*(a*x + S(1))**n, x))
    rubi.add(rule456)

    pattern457 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: NegativeIntegerQ(n/S(2))))
    rule457 = ReplacementRule(pattern457, lambda p, x, m, c, a, d, n : c**(-n/S(2))*Int(x**m*(c + d*x**S(2))**(n/S(2) + p)*(-a*x + S(1))**(-n), x))
    rubi.add(rule457)

    pattern458 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))))
    rule458 = ReplacementRule(pattern458, lambda p, x, m, c, a, d, n : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(-a**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x**m*(-a**S(2)*x**S(2) + S(1))**p*exp(n*atanh(a*x)), x))
    rubi.add(rule458)

    pattern459 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)))
    rule459 = ReplacementRule(pattern459, lambda p, x, c, a, d, u, n : c**p*Int(u*(-a*x + S(1))**(-n/S(2) + p)*(a*x + S(1))**(n/S(2) + p), x))
    rubi.add(rule459)

    pattern460 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: IntegerQ(n/S(2))))
    rule460 = ReplacementRule(pattern460, lambda p, x, c, a, d, u, n : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(-a*x + S(1))**(-FracPart(p))*(a*x + S(1))**(-FracPart(p))*Int(u*(-a*x + S(1))**(-n/S(2) + p)*(a*x + S(1))**(n/S(2) + p), x))
    rubi.add(rule460)

    pattern461 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda c, p: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))))
    rule461 = ReplacementRule(pattern461, lambda p, x, c, a, d, u, n : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(-a**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(u*(-a**S(2)*x**S(2) + S(1))**p*exp(n*atanh(a*x)), x))
    rubi.add(rule461)

    pattern462 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*d + c)), CustomConstraint(lambda p: IntegerQ(p)))
    rule462 = ReplacementRule(pattern462, lambda p, x, c, a, d, u, n : d**p*Int(u*x**(-S(2)*p)*(-a**S(2)*x**S(2) + S(1))**p*exp(n*atanh(a*x)), x))
    rubi.add(rule462)

    pattern463 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*d + c)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(n/S(2))), CustomConstraint(lambda c: PositiveQ(c)))
    rule463 = ReplacementRule(pattern463, lambda p, x, c, a, d, u, n : c**p*Int(u*(S(1) - S(1)/(a*x))**p*(S(1) + S(1)/(a*x))**p*exp(n*atanh(a*x)), x))
    rubi.add(rule463)

    pattern464 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*d + c)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(n/S(2))), CustomConstraint(lambda c: Not(PositiveQ(c))))
    rule464 = ReplacementRule(pattern464, lambda p, x, c, a, d, u, n : x**(S(2)*p)*(c + d/x**S(2))**p*(-a*x + S(1))**(-p)*(a*x + S(1))**(-p)*Int(u*x**(-S(2)*p)*(-a*x + S(1))**p*(a*x + S(1))**p*exp(n*atanh(a*x)), x))
    rubi.add(rule464)

    pattern465 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*d + c)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))))
    rule465 = ReplacementRule(pattern465, lambda p, x, c, a, d, u, n : x**(S(2)*p)*(c + d/x**S(2))**p*(c*x**S(2)/d + S(1))**(-p)*Int(u*x**(-S(2)*p)*(c*x**S(2)/d + S(1))**p*exp(n*atanh(a*x)), x))
    rubi.add(rule465)

    pattern466 = Pattern(Integral(exp(WC('n', S(1))*atanh((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule466 = ReplacementRule(pattern466, lambda x, b, c, a, n : Int((-a*c - b*c*x + S(1))**(-n/S(2))*(a*c + b*c*x + S(1))**(n/S(2)), x))
    rubi.add(rule466)

    pattern467 = Pattern(Integral(x_**m_*exp(n_*atanh((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(S(-1), n, S(1))))
    rule467 = ReplacementRule(pattern467, lambda x, b, m, c, a, n : S(4)*b**(-m + S(-1))*c**(-m + S(-1))*Subst(Int(x**(S(2)/n)*(x**(S(2)/n) + S(1))**(-m + S(-2))*(-a*c + x**(S(2)/n)*(-a*c + S(1)) + S(-1))**m, x), x, (-c*(a + b*x) + S(1))**(-n/S(2))*(c*(a + b*x) + S(1))**(n/S(2)))/n)
    rubi.add(rule467)

    pattern468 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*exp(WC('n', S(1))*atanh((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule468 = ReplacementRule(pattern468, lambda x, b, m, c, a, d, e, n : Int((d + e*x)**m*(-a*c - b*c*x + S(1))**(-n/S(2))*(a*c + b*c*x + S(1))**(n/S(2)), x))
    rubi.add(rule468)

    pattern469 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(a_ + x_*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, a, d, b: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda c, a, e, b: ZeroQ(b**S(2)*c + e*(-a**S(2) + S(1)))), CustomConstraint(lambda c, a, p: IntegerQ(p) | PositiveQ(c/(-a**S(2) + S(1)))))
    rule469 = ReplacementRule(pattern469, lambda p, x, b, c, a, d, u, e, n : (c/(-a**S(2) + S(1)))**p*Int(u*(-a - b*x + S(1))**(-n/S(2) + p)*(a + b*x + S(1))**(n/S(2) + p), x))
    rubi.add(rule469)

    pattern470 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(a_ + x_*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, a, d, b: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda c, a, e, b: ZeroQ(b**S(2)*c + e*(-a**S(2) + S(1)))), CustomConstraint(lambda c, a, p: Not(IntegerQ(p) | PositiveQ(c/(-a**S(2) + S(1))))))
    rule470 = ReplacementRule(pattern470, lambda p, x, b, c, a, d, u, e, n : (c + d*x + e*x**S(2))**p*(-a**S(2) - S(2)*a*b*x - b**S(2)*x**S(2) + S(1))**(-p)*Int(u*(-a**S(2) - S(2)*a*b*x - b**S(2)*x**S(2) + S(1))**p*exp(n*atanh(a*x)), x))
    rubi.add(rule470)

    pattern471 = Pattern(Integral(WC('u', S(1))*exp(WC('n', S(1))*atanh(WC('c', S(1))/(x_*WC('b', S(1)) + WC('a', S(0))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule471 = ReplacementRule(pattern471, lambda x, b, c, a, u, n : Int(u*exp(n*acoth(a/c + b*x/c)), x))
    rubi.add(rule471)

    pattern472 = Pattern(Integral(WC('u', S(1))*exp(n_*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: IntegerQ(n/S(2))))
    rule472 = ReplacementRule(pattern472, lambda x, u, a, n : (S(-1))**(n/S(2))*Int(u*exp(n*atanh(a*x)), x))
    rubi.add(rule472)

    pattern473 = Pattern(Integral(exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: OddQ(n)))
    rule473 = ReplacementRule(pattern473, lambda x, a, n : -Subst(Int((S(1) - x/a)**(-n/S(2) + S(1)/2)*(S(1) + x/a)**(n/S(2) + S(1)/2)/(x**S(2)*sqrt(S(1) - x**S(2)/a**S(2))), x), x, 1/x))
    rubi.add(rule473)

    pattern474 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule474 = ReplacementRule(pattern474, lambda x, a, n, m : -Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2) + S(1)/2)*(S(1) + x/a)**(n/S(2) + S(1)/2)/sqrt(S(1) - x**S(2)/a**S(2)), x), x, 1/x))
    rubi.add(rule474)

    pattern475 = Pattern(Integral(exp(n_*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule475 = ReplacementRule(pattern475, lambda x, a, n : -Subst(Int((S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2))/x**S(2), x), x, 1/x))
    rubi.add(rule475)

    pattern476 = Pattern(Integral(x_**WC('m', S(1))*exp(n_*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n))), CustomConstraint(lambda m: IntegerQ(m)))
    rule476 = ReplacementRule(pattern476, lambda x, a, n, m : -Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2)), x), x, 1/x))
    rubi.add(rule476)

    pattern477 = Pattern(Integral(x_**m_*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule477 = ReplacementRule(pattern477, lambda x, a, n, m : -x**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2) + S(1)/2)*(S(1) + x/a)**(n/S(2) + S(1)/2)/sqrt(S(1) - x**S(2)/a**S(2)), x), x, 1/x))
    rubi.add(rule477)

    pattern478 = Pattern(Integral(x_**m_*exp(n_*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n))), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule478 = ReplacementRule(pattern478, lambda x, a, n, m : -x**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2)), x), x, 1/x))
    rubi.add(rule478)

    pattern479 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a*c + d)), CustomConstraint(lambda p, n: ZeroQ(-n/S(2) + p)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))))
    rule479 = ReplacementRule(pattern479, lambda p, x, c, a, d, n : (c + d*x)**p*(a*x + S(1))*exp(n*acoth(a*x))/(a*(p + S(1))))
    rubi.add(rule479)

    pattern480 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c**S(2) - d**S(2))), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda p: IntegerQ(p)))
    rule480 = ReplacementRule(pattern480, lambda p, x, c, a, d, u, n : d**p*Int(u*x**p*(c/(d*x) + S(1))**p*exp(n*acoth(a*x)), x))
    rubi.add(rule480)

    pattern481 = Pattern(Integral((c_ + x_*WC('d', S(1)))**p_*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c**S(2) - d**S(2))), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule481 = ReplacementRule(pattern481, lambda p, x, c, a, d, u, n : x**(-p)*(c + d*x)**p*(c/(d*x) + S(1))**(-p)*Int(u*x**p*(c/(d*x) + S(1))**p*exp(n*acoth(a*x)), x))
    rubi.add(rule481)

    pattern482 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a*d + c)), CustomConstraint(lambda n: IntegerQ(n/S(2) + S(-1)/2)), CustomConstraint(lambda p, n: IntegerQ(p) | ZeroQ(-n/S(2) + p) | ZeroQ(-n/S(2) + p + S(-1))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule482 = ReplacementRule(pattern482, lambda p, x, c, a, d, n : -c**n*Subst(Int((S(1) - x**S(2)/a**S(2))**(n/S(2))*(c + d*x)**(-n + p)/x**S(2), x), x, 1/x))
    rubi.add(rule482)

    pattern483 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a*d + c)), CustomConstraint(lambda n: IntegerQ(n/S(2) + S(-1)/2)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p, n, m: IntegerQ(p) | Less(S(-5), m, S(-1)) | ZeroQ(-n/S(2) + p) | ZeroQ(-n/S(2) + p + S(-1))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule483 = ReplacementRule(pattern483, lambda p, x, m, c, a, d, n : -c**n*Subst(Int(x**(-m + S(-2))*(S(1) - x**S(2)/a**S(2))**(n/S(2))*(c + d*x)**(-n + p), x), x, 1/x))
    rubi.add(rule483)

    pattern484 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(-a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)))
    rule484 = ReplacementRule(pattern484, lambda p, x, c, a, d, n : -c**p*Subst(Int((S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2))*(S(1) + d*x/c)**p/x**S(2), x), x, 1/x))
    rubi.add(rule484)

    pattern485 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(-a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda m: IntegerQ(m)))
    rule485 = ReplacementRule(pattern485, lambda p, x, m, c, a, d, n : -c**p*Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2))*(S(1) + d*x/c)**p, x), x, 1/x))
    rubi.add(rule485)

    pattern486 = Pattern(Integral(x_**m_*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(-a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule486 = ReplacementRule(pattern486, lambda p, x, m, c, a, d, n : -c**p*x**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2))*(S(1) + d*x/c)**p, x), x, 1/x))
    rubi.add(rule486)

    pattern487 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(-a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda c, p: Not(IntegerQ(p) | PositiveQ(c))))
    rule487 = ReplacementRule(pattern487, lambda p, x, c, a, d, u, n : (S(1) + d/(c*x))**(-p)*(c + d/x)**p*Int(u*(S(1) + d/(c*x))**p*exp(n*acoth(a*x)), x))
    rubi.add(rule487)

    pattern488 = Pattern(Integral(exp(WC('n', S(1))*acoth(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))))
    rule488 = ReplacementRule(pattern488, lambda x, c, a, d, n : exp(n*acoth(a*x))/(a*c*n))
    rubi.add(rule488)

    pattern489 = Pattern(Integral(exp(n_*acoth(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule489 = ReplacementRule(pattern489, lambda x, c, a, d, n : (-a*x + n)*exp(n*acoth(a*x))/(a*c*sqrt(c + d*x**S(2))*(n**S(2) + S(-1))))
    rubi.add(rule489)

    pattern490 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) - S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p, n: IntegerQ(p) | Not(IntegerQ(n))))
    rule490 = ReplacementRule(pattern490, lambda p, x, c, a, d, n : -S(2)*(p + S(1))*(S(2)*p + S(3))*Int((c + d*x**S(2))**(p + S(1))*exp(n*acoth(a*x)), x)/(c*(n**S(2) - S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*acoth(a*x))/(a*c*(n**S(2) - S(4)*(p + S(1))**S(2))))
    rubi.add(rule490)

    pattern491 = Pattern(Integral(x_*exp(n_*acoth(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule491 = ReplacementRule(pattern491, lambda x, c, a, d, n : (a*n*x + S(-1))*exp(n*acoth(a*x))/(a**S(2)*c*sqrt(c + d*x**S(2))*(n**S(2) + S(-1))))
    rubi.add(rule491)

    pattern492 = Pattern(Integral(x_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: LessEqual(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) - S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p, n: IntegerQ(p) | Not(IntegerQ(n))))
    rule492 = ReplacementRule(pattern492, lambda p, x, c, a, d, n : -n*(S(2)*p + S(3))*Int((c + d*x**S(2))**(p + S(1))*exp(n*acoth(a*x)), x)/(a*c*(n**S(2) - S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(a*n*x + S(2)*p + S(2))*exp(n*acoth(a*x))/(a**S(2)*c*(n**S(2) - S(4)*(p + S(1))**S(2))))
    rubi.add(rule492)

    pattern493 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda p, n: ZeroQ(n**S(2) + S(2)*p + S(2))), CustomConstraint(lambda n: NonzeroQ(n**S(2) + S(-1))))
    rule493 = ReplacementRule(pattern493, lambda p, x, c, a, d, n : (c + d*x**S(2))**(p + S(1))*(-S(2)*a*x*(p + S(1)) - n)*exp(n*acoth(a*x))/(a**S(3)*c*n**S(2)*(n**S(2) + S(-1))))
    rubi.add(rule493)

    pattern494 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: LessEqual(p, S(-1))), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) + S(2)*p + S(2))), CustomConstraint(lambda p, n: NonzeroQ(n**S(2) - S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p, n: IntegerQ(p) | Not(IntegerQ(n))))
    rule494 = ReplacementRule(pattern494, lambda p, x, c, a, d, n : -(n**S(2) + S(2)*p + S(2))*Int((c + d*x**S(2))**(p + S(1))*exp(n*acoth(a*x)), x)/(a**S(2)*c*(n**S(2) - S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*acoth(a*x))/(a**S(3)*c*(n**S(2) - S(4)*(p + S(1))**S(2))))
    rubi.add(rule494)

    pattern495 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p, m: LessEqual(S(3), m, -S(2)*p + S(-2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule495 = ReplacementRule(pattern495, lambda p, x, m, c, a, d, n : -a**(-m + S(-1))*(-c)**p*Subst(Int(Cosh(x)**(-S(2)*p + S(-2))*exp(n*x)*coth(x)**(m + S(2)*p + S(2)), x), x, acoth(a*x)))
    rubi.add(rule495)

    pattern496 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda p: IntegerQ(p)))
    rule496 = ReplacementRule(pattern496, lambda p, x, c, a, d, u, n : d**p*Int(u*x**(S(2)*p)*(S(1) - S(1)/(a**S(2)*x**S(2)))**p*exp(n*acoth(a*x)), x))
    rubi.add(rule496)

    pattern497 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule497 = ReplacementRule(pattern497, lambda p, x, c, a, d, u, n : x**(-S(2)*p)*(S(1) - S(1)/(a**S(2)*x**S(2)))**(-p)*(c + d*x**S(2))**p*Int(u*x**(S(2)*p)*(S(1) - S(1)/(a**S(2)*x**S(2)))**p*exp(n*acoth(a*x)), x))
    rubi.add(rule497)

    pattern498 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*d + c)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda p, n: IntegersQ(S(2)*p, n/S(2) + p)))
    rule498 = ReplacementRule(pattern498, lambda p, x, c, a, d, u, n : a**(-S(2)*p)*c**p*Int(u*x**(-S(2)*p)*(a*x + S(-1))**(-n/S(2) + p)*(a*x + S(1))**(n/S(2) + p), x))
    rubi.add(rule498)

    pattern499 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*d + c)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda p, n: Not(IntegersQ(S(2)*p, n/S(2) + p))))
    rule499 = ReplacementRule(pattern499, lambda p, x, c, a, d, n : -c**p*Subst(Int((S(1) - x/a)**(-n/S(2) + p)*(S(1) + x/a)**(n/S(2) + p)/x**S(2), x), x, 1/x))
    rubi.add(rule499)

    pattern500 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*d + c)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda p, n: Not(IntegersQ(S(2)*p, n/S(2) + p))), CustomConstraint(lambda m: IntegerQ(m)))
    rule500 = ReplacementRule(pattern500, lambda p, x, m, c, a, d, n : -c**p*Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2) + p)*(S(1) + x/a)**(n/S(2) + p), x), x, 1/x))
    rubi.add(rule500)

    pattern501 = Pattern(Integral(x_**m_*(c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*d + c)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda c, p: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda p, n: Not(IntegersQ(S(2)*p, n/S(2) + p))), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule501 = ReplacementRule(pattern501, lambda p, x, m, c, a, d, n : -c**p*x**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2) + p)*(S(1) + x/a)**(n/S(2) + p), x), x, 1/x))
    rubi.add(rule501)

    pattern502 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda c, a, d: ZeroQ(a**S(2)*d + c)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda c, p: Not(IntegerQ(p) | PositiveQ(c))))
    rule502 = ReplacementRule(pattern502, lambda p, x, c, a, d, u, n : c**IntPart(p)*(S(1) - S(1)/(a**S(2)*x**S(2)))**(-FracPart(p))*(c + d/x**S(2))**FracPart(p)*Int(u*(S(1) - S(1)/(a**S(2)*x**S(2)))**p*exp(n*acoth(a*x)), x))
    rubi.add(rule502)

    pattern503 = Pattern(Integral(WC('u', S(1))*exp(n_*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(n/S(2))))
    rule503 = ReplacementRule(pattern503, lambda x, b, c, a, u, n : (S(-1))**(n/S(2))*Int(u*exp(n*atanh(c*(a + b*x))), x))
    rubi.add(rule503)

    pattern504 = Pattern(Integral(exp(WC('n', S(1))*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))))
    rule504 = ReplacementRule(pattern504, lambda x, b, c, a, n : (c*(a + b*x))**(n/S(2))*(S(1) + S(1)/(c*(a + b*x)))**(n/S(2))*(a*c + b*c*x + S(1))**(-n/S(2))*Int((a*c + b*c*x + S(-1))**(-n/S(2))*(a*c + b*c*x + S(1))**(n/S(2)), x))
    rubi.add(rule504)

    pattern505 = Pattern(Integral(x_**m_*exp(n_*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(S(-1), n, S(1))))
    rule505 = ReplacementRule(pattern505, lambda x, b, m, c, a, n : -S(4)*b**(-m + S(-1))*c**(-m + S(-1))*Subst(Int(x**(S(2)/n)*(x**(S(2)/n) + S(-1))**(-m + S(-2))*(a*c + x**(S(2)/n)*(-a*c + S(1)) + S(1))**m, x), x, (S(1) - S(1)/(c*(a + b*x)))**(-n/S(2))*(S(1) + S(1)/(c*(a + b*x)))**(n/S(2)))/n)
    rubi.add(rule505)

    pattern506 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*exp(WC('n', S(1))*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))))
    rule506 = ReplacementRule(pattern506, lambda x, b, m, c, a, d, e, n : (c*(a + b*x))**(n/S(2))*(S(1) + S(1)/(c*(a + b*x)))**(n/S(2))*(a*c + b*c*x + S(1))**(-n/S(2))*Int((d + e*x)**m*(a*c + b*c*x + S(-1))**(-n/S(2))*(a*c + b*c*x + S(1))**(n/S(2)), x))
    rubi.add(rule506)

    pattern507 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acoth(a_ + x_*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda e, a, d, b: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda c, a, e, b: ZeroQ(b**S(2)*c + e*(-a**S(2) + S(1)))), CustomConstraint(lambda c, a, p: IntegerQ(p) | PositiveQ(c/(-a**S(2) + S(1)))))
    rule507 = ReplacementRule(pattern507, lambda p, x, b, c, a, d, u, e, n : (c/(-a**S(2) + S(1)))**p*((a + b*x + S(1))/(a + b*x))**(n/S(2))*((a + b*x)/(a + b*x + S(1)))**(n/S(2))*(-a - b*x + S(1))**(n/S(2))*(a + b*x + S(-1))**(-n/S(2))*Int(u*(-a - b*x + S(1))**(-n/S(2) + p)*(a + b*x + S(1))**(n/S(2) + p), x))
    rubi.add(rule507)

    pattern508 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acoth(a_ + x_*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: Not(IntegerQ(n/S(2)))), CustomConstraint(lambda e, a, d, b: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda c, a, e, b: ZeroQ(b**S(2)*c + e*(-a**S(2) + S(1)))), CustomConstraint(lambda c, a, p: Not(IntegerQ(p) | PositiveQ(c/(-a**S(2) + S(1))))))
    rule508 = ReplacementRule(pattern508, lambda p, x, b, c, a, d, u, e, n : (c + d*x + e*x**S(2))**p*(-a**S(2) - S(2)*a*b*x - b**S(2)*x**S(2) + S(1))**(-p)*Int(u*(-a**S(2) - S(2)*a*b*x - b**S(2)*x**S(2) + S(1))**p*exp(n*acoth(a*x)), x))
    rubi.add(rule508)

    pattern509 = Pattern(Integral(WC('u', S(1))*exp(WC('n', S(1))*acoth(WC('c', S(1))/(x_*WC('b', S(1)) + WC('a', S(0))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule509 = ReplacementRule(pattern509, lambda x, b, c, a, u, n : Int(u*exp(n*atanh(a/c + b*x/c)), x))
    rubi.add(rule509)

    pattern510 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule510 = ReplacementRule(pattern510, lambda x, b, c, a, d, n : Subst(Int((a + b*atanh(x))**n, x), x, c + d*x)/d)
    rubi.add(rule510)

    pattern511 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule511 = ReplacementRule(pattern511, lambda x, b, c, a, d, n : Subst(Int((a + b*acoth(x))**n, x), x, c + d*x)/d)
    rubi.add(rule511)

    pattern512 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule512 = ReplacementRule(pattern512, lambda x, b, c, a, d, n : Int((a + b*atanh(c + d*x))**n, x))
    rubi.add(rule512)

    pattern513 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule513 = ReplacementRule(pattern513, lambda x, b, c, a, d, n : Int((a + b*acoth(c + d*x))**n, x))
    rubi.add(rule513)

    pattern514 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule514 = ReplacementRule(pattern514, lambda x, b, m, f, c, a, d, e, n : Subst(Int((a + b*atanh(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule514)

    pattern515 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule515 = ReplacementRule(pattern515, lambda x, b, m, f, c, a, d, e, n : Subst(Int((a + b*acoth(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule515)

    pattern516 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule516 = ReplacementRule(pattern516, lambda x, b, m, f, c, a, d, e, n : Int((a + b*atanh(c + d*x))**n*(e + f*x)**m, x))
    rubi.add(rule516)

    pattern517 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule517 = ReplacementRule(pattern517, lambda x, b, m, f, c, a, d, e, n : Int((a + b*acoth(c + d*x))**n*(e + f*x)**m, x))
    rubi.add(rule517)

    pattern518 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda B, c, C, d: ZeroQ(-B*d + S(2)*C*c)))
    rule518 = ReplacementRule(pattern518, lambda B, A, C, p, x, b, c, a, d, n : Subst(Int((a + b*atanh(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule518)

    pattern519 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda B, c, C, d: ZeroQ(-B*d + S(2)*C*c)))
    rule519 = ReplacementRule(pattern519, lambda B, A, C, p, x, b, c, a, d, n : Subst(Int((a + b*acoth(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule519)

    pattern520 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda B, c, C, d: ZeroQ(-B*d + S(2)*C*c)))
    rule520 = ReplacementRule(pattern520, lambda B, A, C, p, x, b, m, f, c, a, d, e, n : Subst(Int((a + b*atanh(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule520)

    pattern521 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda B, c, C, d: ZeroQ(-B*d + S(2)*C*c)))
    rule521 = ReplacementRule(pattern521, lambda B, A, C, p, x, b, m, f, c, a, d, e, n : Subst(Int((a + b*acoth(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule521)

    pattern522 = Pattern(Integral(atanh(a_ + x_*WC('b', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: RationalQ(n)))
    rule522 = ReplacementRule(pattern522, lambda x, b, c, a, d, n : -Int(log(-a - b*x + S(1))/(c + d*x**n), x)/S(2) + Int(log(a + b*x + S(1))/(c + d*x**n), x)/S(2))
    rubi.add(rule522)

    pattern523 = Pattern(Integral(acoth(a_ + x_*WC('b', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: RationalQ(n)))
    rule523 = ReplacementRule(pattern523, lambda x, b, c, a, d, n : -Int(log((a + b*x + S(-1))/(a + b*x))/(c + d*x**n), x)/S(2) + Int(log((a + b*x + S(1))/(a + b*x))/(c + d*x**n), x)/S(2))
    rubi.add(rule523)

    pattern524 = Pattern(Integral(atanh(a_ + x_*WC('b', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(RationalQ(n))))
    rule524 = ReplacementRule(pattern524, lambda x, b, c, a, d, n : Int(atanh(a + b*x)/(c + d*x**n), x))
    rubi.add(rule524)

    pattern525 = Pattern(Integral(acoth(a_ + x_*WC('b', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(RationalQ(n))))
    rule525 = ReplacementRule(pattern525, lambda x, b, c, a, d, n : Int(acoth(a + b*x)/(c + d*x**n), x))
    rubi.add(rule525)

    pattern526 = Pattern(Integral(atanh(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule526 = ReplacementRule(pattern526, lambda x, a, n, b : -b*n*Int(x**n/(-a**S(2) - S(2)*a*b*x**n - b**S(2)*x**(S(2)*n) + S(1)), x) + x*atanh(a + b*x**n))
    rubi.add(rule526)

    pattern527 = Pattern(Integral(acoth(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule527 = ReplacementRule(pattern527, lambda x, a, n, b : -b*n*Int(x**n/(-a**S(2) - S(2)*a*b*x**n - b**S(2)*x**(S(2)*n) + S(1)), x) + x*acoth(a + b*x**n))
    rubi.add(rule527)

    pattern528 = Pattern(Integral(atanh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule528 = ReplacementRule(pattern528, lambda x, a, n, b : -Int(log(-a - b*x**n + S(1))/x, x)/S(2) + Int(log(a + b*x**n + S(1))/x, x)/S(2))
    rubi.add(rule528)

    pattern529 = Pattern(Integral(acoth(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule529 = ReplacementRule(pattern529, lambda x, a, n, b : -Int(log(S(1) - S(1)/(a + b*x**n))/x, x)/S(2) + Int(log(S(1) + 1/(a + b*x**n))/x, x)/S(2))
    rubi.add(rule529)

    pattern530 = Pattern(Integral(x_**WC('m', S(1))*atanh(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda m: Unequal(m + S(1), S(0))), CustomConstraint(lambda n, m: Unequal(m + S(1), n)))
    rule530 = ReplacementRule(pattern530, lambda x, b, m, a, n : -b*n*Int(x**(m + n)/(-a**S(2) - S(2)*a*b*x**n - b**S(2)*x**(S(2)*n) + S(1)), x)/(m + S(1)) + x**(m + S(1))*atanh(a + b*x**n)/(m + S(1)))
    rubi.add(rule530)

    pattern531 = Pattern(Integral(x_**WC('m', S(1))*acoth(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda m: Unequal(m + S(1), S(0))), CustomConstraint(lambda n, m: Unequal(m + S(1), n)))
    rule531 = ReplacementRule(pattern531, lambda x, b, m, a, n : -b*n*Int(x**(m + n)/(-a**S(2) - S(2)*a*b*x**n - b**S(2)*x**(S(2)*n) + S(1)), x)/(m + S(1)) + x**(m + S(1))*acoth(a + b*x**n)/(m + S(1)))
    rubi.add(rule531)

    pattern532 = Pattern(Integral(atanh(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule532 = ReplacementRule(pattern532, lambda x, b, c, f, a, d : -Int(log(-a - b*f**(c + d*x) + S(1)), x)/S(2) + Int(log(a + b*f**(c + d*x) + S(1)), x)/S(2))
    rubi.add(rule532)

    pattern533 = Pattern(Integral(acoth(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule533 = ReplacementRule(pattern533, lambda x, b, c, f, a, d : -Int(log(S(1) - S(1)/(a + b*f**(c + d*x))), x)/S(2) + Int(log(S(1) + 1/(a + b*f**(c + d*x))), x)/S(2))
    rubi.add(rule533)

    pattern534 = Pattern(Integral(x_**WC('m', S(1))*atanh(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule534 = ReplacementRule(pattern534, lambda x, b, m, c, f, a, d : -Int(x**m*log(-a - b*f**(c + d*x) + S(1)), x)/S(2) + Int(x**m*log(a + b*f**(c + d*x) + S(1)), x)/S(2))
    rubi.add(rule534)

    pattern535 = Pattern(Integral(x_**WC('m', S(1))*acoth(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule535 = ReplacementRule(pattern535, lambda x, b, m, c, f, a, d : -Int(x**m*log(S(1) - S(1)/(a + b*f**(c + d*x))), x)/S(2) + Int(x**m*log(S(1) + 1/(a + b*f**(c + d*x))), x)/S(2))
    rubi.add(rule535)

    pattern536 = Pattern(Integral(WC('u', S(1))*atanh(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule536 = ReplacementRule(pattern536, lambda x, b, m, c, a, u, n : Int(u*acoth(a/c + b*x**n/c)**m, x))
    rubi.add(rule536)

    pattern537 = Pattern(Integral(WC('u', S(1))*acoth(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule537 = ReplacementRule(pattern537, lambda x, b, m, c, a, u, n : Int(u*atanh(a/c + b*x**n/c)**m, x))
    rubi.add(rule537)

    pattern538 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0)))*atanh(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda c, b: ZeroQ(b - c**S(2))))
    rule538 = ReplacementRule(pattern538, lambda c, a, b, x : log(atanh(c*x/sqrt(a + b*x**S(2))))/c)
    rubi.add(rule538)

    pattern539 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0)))*acoth(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda c, b: ZeroQ(b - c**S(2))))
    rule539 = ReplacementRule(pattern539, lambda c, a, b, x : -log(acoth(c*x/sqrt(a + b*x**S(2))))/c)
    rubi.add(rule539)

    pattern540 = Pattern(Integral(atanh(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, b: ZeroQ(b - c**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule540 = ReplacementRule(pattern540, lambda x, b, m, c, a : atanh(c*x/sqrt(a + b*x**S(2)))**(m + S(1))/(c*(m + S(1))))
    rubi.add(rule540)

    pattern541 = Pattern(Integral(acoth(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, b: ZeroQ(b - c**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule541 = ReplacementRule(pattern541, lambda x, b, m, c, a : -acoth(c*x/sqrt(a + b*x**S(2)))**(m + S(1))/(c*(m + S(1))))
    rubi.add(rule541)

    pattern542 = Pattern(Integral(atanh(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, b: ZeroQ(b - c**S(2))), CustomConstraint(lambda e, a, d, b: ZeroQ(-a*e + b*d)))
    rule542 = ReplacementRule(pattern542, lambda x, b, m, c, a, d, e : sqrt(a + b*x**S(2))*Int(atanh(c*x/sqrt(a + b*x**S(2)))**m/sqrt(a + b*x**S(2)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule542)

    pattern543 = Pattern(Integral(acoth(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda c, b: ZeroQ(b - c**S(2))), CustomConstraint(lambda e, a, d, b: ZeroQ(-a*e + b*d)))
    rule543 = ReplacementRule(pattern543, lambda x, b, m, c, a, d, e : sqrt(a + b*x**S(2))*Int(acoth(c*x/sqrt(a + b*x**S(2)))**m/sqrt(a + b*x**S(2)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule543)

    pattern544 = Pattern(Integral((x_**S(2)*WC('d', S(1)) + WC('c', S(0)))**n_*atanh(x_*WC('a', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(S(2)*n)), CustomConstraint(lambda n: LessEqual(n, S(-1))), )
    def With544(x, c, a, d, n):
        u = IntHide((c + d*x**S(2))**n, x)
        return -a*Int(Dist(1/(-a**S(2)*x**S(2) + S(1)), u, x), x) + Dist(atanh(a*x), u, x)
    rule544 = ReplacementRule(pattern544, lambda x, c, a, d, n : With544(x, c, a, d, n))
    rubi.add(rule544)

    pattern545 = Pattern(Integral((x_**S(2)*WC('d', S(1)) + WC('c', S(0)))**n_*acoth(x_*WC('a', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(S(2)*n)), CustomConstraint(lambda n: LessEqual(n, S(-1))), )
    def With545(x, c, a, d, n):
        u = IntHide((c + d*x**S(2))**n, x)
        return -a*Int(Dist(1/(-a**S(2)*x**S(2) + S(1)), u, x), x) + Dist(acoth(a*x), u, x)
    rule545 = ReplacementRule(pattern545, lambda x, c, a, d, n : With545(x, c, a, d, n))
    rubi.add(rule545)

    pattern546 = Pattern(Integral(u_*v_**WC('n', S(1)), x_), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Less(n, S(0))), CustomConstraint(lambda v, x: PosQ(Discriminant(v, x))), CustomConstraint(lambda x, u: MatchQ(u, Condition(Optional(Pattern(r, Blank))*Pattern(f, Blank)**Pattern(w, Blank)))), CustomConstraint(lambda u, n, v, tmp, x, ArcTanh: Not(FalseQ(tmp)) & SameQ(Head(tmp), ArcTanh) & ZeroQ(-D(v, x)**S(2) + Discriminant(v, x)*Part(tmp, S(1))**S(2))))
    def With546(v, x, u, n):
        tmp = InverseFunctionOfLinear(u, x)
        return (-Discriminant(v, x)/(S(4)*Coefficient(v, x, S(2))))**n*Subst(Int(SimplifyIntegrand(SubstForInverseFunction(u, tmp, x)*sech(x)**(S(2)*n + S(2)), x), x), x, tmp)/Coefficient(Part(tmp, S(1)), x, S(1))
    rule546 = ReplacementRule(pattern546, lambda v, x, u, n : With546(v, x, u, n))
    rubi.add(rule546)

    pattern547 = Pattern(Integral(u_*v_**WC('n', S(1)), x_), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Less(n, S(0))), CustomConstraint(lambda v, x: PosQ(Discriminant(v, x))), CustomConstraint(lambda x, u: MatchQ(u, Condition(Optional(Pattern(r, Blank))*Pattern(f, Blank)**Pattern(w, Blank)))), CustomConstraint(lambda u, n, ArcCoth, v, tmp, x: Not(FalseQ(tmp)) & SameQ(Head(tmp), ArcCoth) & ZeroQ(-D(v, x)**S(2) + Discriminant(v, x)*Part(tmp, S(1))**S(2))))
    def With547(v, x, u, n):
        tmp = InverseFunctionOfLinear(u, x)
        return (-Discriminant(v, x)/(S(4)*Coefficient(v, x, S(2))))**n*Subst(Int(SimplifyIntegrand((-csch(x)**S(2))**(n + S(1))*SubstForInverseFunction(u, tmp, x), x), x), x, tmp)/Coefficient(Part(tmp, S(1)), x, S(1))
    rule547 = ReplacementRule(pattern547, lambda v, x, u, n : With547(v, x, u, n))
    rubi.add(rule547)

    pattern548 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(-1))))
    rule548 = ReplacementRule(pattern548, lambda x, b, c, a, d : b*Int(x/(c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*atanh(c + d*tanh(a + b*x)))
    rubi.add(rule548)

    pattern549 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(-1))))
    rule549 = ReplacementRule(pattern549, lambda x, b, c, a, d : b*Int(x/(c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*acoth(c + d*tanh(a + b*x)))
    rubi.add(rule549)

    pattern550 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(-1))))
    rule550 = ReplacementRule(pattern550, lambda x, b, c, a, d : b*Int(x/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*atanh(c + d*coth(a + b*x)))
    rubi.add(rule550)

    pattern551 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(-1))))
    rule551 = ReplacementRule(pattern551, lambda x, b, c, a, d : b*Int(x/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*acoth(c + d*coth(a + b*x)))
    rubi.add(rule551)

    pattern552 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(-1))))
    rule552 = ReplacementRule(pattern552, lambda x, b, c, a, d : b*(-c - d + S(1))*Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x) - b*(c + d + S(1))*Int(x*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x) + x*atanh(c + d*tanh(a + b*x)))
    rubi.add(rule552)

    pattern553 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(-1))))
    rule553 = ReplacementRule(pattern553, lambda x, b, c, a, d : b*(-c - d + S(1))*Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x) - b*(c + d + S(1))*Int(x*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x) + x*acoth(c + d*tanh(a + b*x)))
    rubi.add(rule553)

    pattern554 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(-1))))
    rule554 = ReplacementRule(pattern554, lambda x, b, c, a, d : -b*(-c - d + S(1))*Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x) + b*(c + d + S(1))*Int(x*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x) + x*atanh(c + d*coth(a + b*x)))
    rubi.add(rule554)

    pattern555 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(-1))))
    rule555 = ReplacementRule(pattern555, lambda x, b, c, a, d : -b*(-c - d + S(1))*Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x) + b*(c + d + S(1))*Int(x*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x) + x*acoth(c + d*coth(a + b*x)))
    rubi.add(rule555)

    pattern556 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(-1))))
    rule556 = ReplacementRule(pattern556, lambda x, b, m, c, f, a, d, e : b*Int((e + f*x)**(m + S(1))/(c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*atanh(c + d*tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule556)

    pattern557 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(-1))))
    rule557 = ReplacementRule(pattern557, lambda x, b, m, c, f, a, d, e : b*Int((e + f*x)**(m + S(1))/(c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acoth(c + d*tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule557)

    pattern558 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(-1))))
    rule558 = ReplacementRule(pattern558, lambda x, b, m, c, f, a, d, e : b*Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*atanh(c + d*coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule558)

    pattern559 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(-1))))
    rule559 = ReplacementRule(pattern559, lambda x, b, m, c, f, a, d, e : b*Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acoth(c + d*coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule559)

    pattern560 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(-1))))
    rule560 = ReplacementRule(pattern560, lambda x, b, m, c, f, a, d, e : b*(-c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x)/(f*(m + S(1))) - b*(c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*atanh(c + d*tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule560)

    pattern561 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(-1))))
    rule561 = ReplacementRule(pattern561, lambda x, b, m, c, f, a, d, e : b*(-c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x)/(f*(m + S(1))) - b*(c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acoth(c + d*tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule561)

    pattern562 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(-1))))
    rule562 = ReplacementRule(pattern562, lambda x, b, m, c, f, a, d, e : -b*(-c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x)/(f*(m + S(1))) + b*(c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*atanh(c + d*coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule562)

    pattern563 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(-1))))
    rule563 = ReplacementRule(pattern563, lambda x, b, m, c, f, a, d, e : -b*(-c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x)/(f*(m + S(1))) + b*(c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acoth(c + d*coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule563)

    pattern564 = Pattern(Integral(atanh(tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule564 = ReplacementRule(pattern564, lambda x, a, b : -b*Int(x*sec(S(2)*a + S(2)*b*x), x) + x*atanh(tan(a + b*x)))
    rubi.add(rule564)

    pattern565 = Pattern(Integral(acoth(tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule565 = ReplacementRule(pattern565, lambda x, a, b : -b*Int(x*sec(S(2)*a + S(2)*b*x), x) + x*acoth(tan(a + b*x)))
    rubi.add(rule565)

    pattern566 = Pattern(Integral(atanh(cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule566 = ReplacementRule(pattern566, lambda x, a, b : -b*Int(x*sec(S(2)*a + S(2)*b*x), x) + x*atanh(cot(a + b*x)))
    rubi.add(rule566)

    pattern567 = Pattern(Integral(acoth(cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule567 = ReplacementRule(pattern567, lambda x, a, b : -b*Int(x*sec(S(2)*a + S(2)*b*x), x) + x*acoth(cot(a + b*x)))
    rubi.add(rule567)

    pattern568 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule568 = ReplacementRule(pattern568, lambda x, b, m, f, a, e : -b*Int((e + f*x)**(m + S(1))*sec(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*atanh(tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule568)

    pattern569 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule569 = ReplacementRule(pattern569, lambda x, b, m, f, a, e : -b*Int((e + f*x)**(m + S(1))*sec(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acoth(tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule569)

    pattern570 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule570 = ReplacementRule(pattern570, lambda x, b, m, f, a, e : -b*Int((e + f*x)**(m + S(1))*sec(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*atanh(cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule570)

    pattern571 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule571 = ReplacementRule(pattern571, lambda x, b, m, f, a, e : -b*Int((e + f*x)**(m + S(1))*sec(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acoth(cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule571)

    pattern572 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((ImaginaryI*d + c)**S(2) + S(-1))))
    rule572 = ReplacementRule(pattern572, lambda x, b, c, a, d : ImaginaryI*b*Int(x/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*atanh(c + d*tan(a + b*x)))
    rubi.add(rule572)

    pattern573 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((ImaginaryI*d + c)**S(2) + S(-1))))
    rule573 = ReplacementRule(pattern573, lambda x, b, c, a, d : ImaginaryI*b*Int(x/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*acoth(c + d*tan(a + b*x)))
    rubi.add(rule573)

    pattern574 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((-ImaginaryI*d + c)**S(2) + S(-1))))
    rule574 = ReplacementRule(pattern574, lambda x, b, c, a, d : ImaginaryI*b*Int(x/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*atanh(c + d*cot(a + b*x)))
    rubi.add(rule574)

    pattern575 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((-ImaginaryI*d + c)**S(2) + S(-1))))
    rule575 = ReplacementRule(pattern575, lambda x, b, c, a, d : ImaginaryI*b*Int(x/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*acoth(c + d*cot(a + b*x)))
    rubi.add(rule575)

    pattern576 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((ImaginaryI*d + c)**S(2) + S(-1))))
    rule576 = ReplacementRule(pattern576, lambda x, b, c, a, d : -ImaginaryI*b*(-ImaginaryI*d + c + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*d + c + (-ImaginaryI*d + c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + ImaginaryI*b*(ImaginaryI*d - c + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*d - c + (ImaginaryI*d - c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*atanh(c + d*tan(a + b*x)))
    rubi.add(rule576)

    pattern577 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((ImaginaryI*d + c)**S(2) + S(-1))))
    rule577 = ReplacementRule(pattern577, lambda x, b, c, a, d : -ImaginaryI*b*(-ImaginaryI*d + c + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*d + c + (-ImaginaryI*d + c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + ImaginaryI*b*(ImaginaryI*d - c + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*d - c + (ImaginaryI*d - c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*acoth(c + d*tan(a + b*x)))
    rubi.add(rule577)

    pattern578 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((-ImaginaryI*d + c)**S(2) + S(-1))))
    rule578 = ReplacementRule(pattern578, lambda x, b, c, a, d : -ImaginaryI*b*(-ImaginaryI*d - c + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*d - c - (-ImaginaryI*d - c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + ImaginaryI*b*(ImaginaryI*d + c + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*d + c - (ImaginaryI*d + c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*atanh(c + d*cot(a + b*x)))
    rubi.add(rule578)

    pattern579 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((-ImaginaryI*d + c)**S(2) + S(-1))))
    rule579 = ReplacementRule(pattern579, lambda x, b, c, a, d : -ImaginaryI*b*(-ImaginaryI*d - c + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*d - c - (-ImaginaryI*d - c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + ImaginaryI*b*(ImaginaryI*d + c + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*d + c - (ImaginaryI*d + c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*acoth(c + d*cot(a + b*x)))
    rubi.add(rule579)

    pattern580 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((ImaginaryI*d + c)**S(2) + S(-1))))
    rule580 = ReplacementRule(pattern580, lambda x, b, m, c, f, a, d, e : ImaginaryI*b*Int((e + f*x)**(m + S(1))/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*atanh(c + d*tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule580)

    pattern581 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((ImaginaryI*d + c)**S(2) + S(-1))))
    rule581 = ReplacementRule(pattern581, lambda x, b, m, c, f, a, d, e : ImaginaryI*b*Int((e + f*x)**(m + S(1))/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acoth(c + d*tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule581)

    pattern582 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((-ImaginaryI*d + c)**S(2) + S(-1))))
    rule582 = ReplacementRule(pattern582, lambda x, b, m, c, f, a, d, e : ImaginaryI*b*Int((e + f*x)**(m + S(1))/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*atanh(c + d*cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule582)

    pattern583 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((-ImaginaryI*d + c)**S(2) + S(-1))))
    rule583 = ReplacementRule(pattern583, lambda x, b, m, c, f, a, d, e : ImaginaryI*b*Int((e + f*x)**(m + S(1))/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acoth(c + d*cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule583)

    pattern584 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((ImaginaryI*d + c)**S(2) + S(-1))))
    rule584 = ReplacementRule(pattern584, lambda x, b, m, c, f, a, d, e : -ImaginaryI*b*(-ImaginaryI*d + c + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*d + c + (-ImaginaryI*d + c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + ImaginaryI*b*(ImaginaryI*d - c + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*d - c + (ImaginaryI*d - c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*atanh(c + d*tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule584)

    pattern585 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((ImaginaryI*d + c)**S(2) + S(-1))))
    rule585 = ReplacementRule(pattern585, lambda x, b, m, c, f, a, d, e : -ImaginaryI*b*(-ImaginaryI*d + c + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*d + c + (-ImaginaryI*d + c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + ImaginaryI*b*(ImaginaryI*d - c + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*d - c + (ImaginaryI*d - c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acoth(c + d*tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule585)

    pattern586 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((-ImaginaryI*d + c)**S(2) + S(-1))))
    rule586 = ReplacementRule(pattern586, lambda x, b, m, c, f, a, d, e : -ImaginaryI*b*(-ImaginaryI*d - c + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*d - c - (-ImaginaryI*d - c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + ImaginaryI*b*(ImaginaryI*d + c + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*d + c - (ImaginaryI*d + c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*atanh(c + d*cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule586)

    pattern587 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((-ImaginaryI*d + c)**S(2) + S(-1))))
    rule587 = ReplacementRule(pattern587, lambda x, b, m, c, f, a, d, e : -ImaginaryI*b*(-ImaginaryI*d - c + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*d - c - (-ImaginaryI*d - c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + ImaginaryI*b*(ImaginaryI*d + c + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*d + c - (ImaginaryI*d + c + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acoth(c + d*cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule587)

    pattern588 = Pattern(Integral(atanh(u_), x_), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)))
    rule588 = ReplacementRule(pattern588, lambda x, u : x*atanh(u) - Int(SimplifyIntegrand(x*D(u, x)/(-u**S(2) + S(1)), x), x))
    rubi.add(rule588)

    pattern589 = Pattern(Integral(acoth(u_), x_), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)))
    rule589 = ReplacementRule(pattern589, lambda x, u : x*acoth(u) - Int(SimplifyIntegrand(x*D(u, x)/(-u**S(2) + S(1)), x), x))
    rubi.add(rule589)

    pattern590 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(u_)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda x, m, c, d, u: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda x, u, m: FalseQ(PowerVariableExpn(u, m + S(1), x))))
    rule590 = ReplacementRule(pattern590, lambda x, b, m, c, a, d, u : -b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(-u**S(2) + S(1)), x), x)/(d*(m + S(1))) + (a + b*atanh(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule590)

    pattern591 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(u_)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda x, m, c, d, u: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda x, u, m: FalseQ(PowerVariableExpn(u, m + S(1), x))))
    rule591 = ReplacementRule(pattern591, lambda x, b, m, c, a, d, u : -b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(-u**S(2) + S(1)), x), x)/(d*(m + S(1))) + (a + b*acoth(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule591)

    pattern592 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*atanh(u_)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda x, b, v, a, u: FalseQ(FunctionOfLinear(v*(a + b*atanh(u)), x))), CustomConstraint(lambda a, u, x, b, w: InverseFunctionFreeQ(w, x)))
    def With592(x, b, v, a, u):
        w = IntHide(v, x)
        return -b*Int(SimplifyIntegrand(w*D(u, x)/(-u**S(2) + S(1)), x), x) + Dist(a + b*atanh(u), w, x)
    rule592 = ReplacementRule(pattern592, lambda x, b, v, a, u : With592(x, b, v, a, u))
    rubi.add(rule592)

    pattern593 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acoth(u_)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda x, b, v, a, u: FalseQ(FunctionOfLinear(v*(a + b*acoth(u)), x))), CustomConstraint(lambda a, u, x, b, w: InverseFunctionFreeQ(w, x)))
    def With593(x, b, v, a, u):
        w = IntHide(v, x)
        return -b*Int(SimplifyIntegrand(w*D(u, x)/(-u**S(2) + S(1)), x), x) + Dist(a + b*acoth(u), w, x)
    rule593 = ReplacementRule(pattern593, lambda x, b, v, a, u : With593(x, b, v, a, u))
    rubi.add(rule593)

    pattern594 = Pattern(Integral(asech(x_*WC('c', S(1))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule594 = ReplacementRule(pattern594, lambda c, x : x*asech(c*x) + sqrt(c*x + S(1))*sqrt(1/(c*x + S(1)))*Int(S(1)/(sqrt(-c*x + S(1))*sqrt(c*x + S(1))), x))
    rubi.add(rule594)

    pattern595 = Pattern(Integral(acsch(x_*WC('c', S(1))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule595 = ReplacementRule(pattern595, lambda c, x : x*acsch(c*x) + Int(S(1)/(x*sqrt(S(1) + S(1)/(c**S(2)*x**S(2)))), x)/c)
    rubi.add(rule595)

    pattern596 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule596 = ReplacementRule(pattern596, lambda x, b, c, a, n : -Subst(Int((a + b*x)**n*tanh(x)*sech(x), x), x, asech(c*x))/c)
    rubi.add(rule596)

    pattern597 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule597 = ReplacementRule(pattern597, lambda x, b, c, a, n : -Subst(Int((a + b*x)**n*coth(x)*csch(x), x), x, acsch(c*x))/c)
    rubi.add(rule597)

    pattern598 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule598 = ReplacementRule(pattern598, lambda c, a, b, x : -Subst(Int((a + b*acosh(x/c))/x, x), x, 1/x))
    rubi.add(rule598)

    pattern599 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule599 = ReplacementRule(pattern599, lambda c, a, b, x : -Subst(Int((a + b*asinh(x/c))/x, x), x, 1/x))
    rubi.add(rule599)

    pattern600 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule600 = ReplacementRule(pattern600, lambda x, b, m, c, a : b*sqrt(c*x + S(1))*sqrt(1/(c*x + S(1)))*Int(x**m/(sqrt(-c*x + S(1))*sqrt(c*x + S(1))), x)/(m + S(1)) + x**(m + S(1))*(a + b*asech(c*x))/(m + S(1)))
    rubi.add(rule600)

    pattern601 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule601 = ReplacementRule(pattern601, lambda x, b, m, c, a : b*Int(x**(m + S(-1))/sqrt(S(1) + S(1)/(c**S(2)*x**S(2))), x)/(c*(m + S(1))) + x**(m + S(1))*(a + b*acsch(c*x))/(m + S(1)))
    rubi.add(rule601)

    pattern602 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule602 = ReplacementRule(pattern602, lambda x, b, m, c, a, n : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*tanh(x)*sech(x)**(m + S(1)), x), x, asech(c*x)))
    rubi.add(rule602)

    pattern603 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule603 = ReplacementRule(pattern603, lambda x, b, m, c, a, n : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*coth(x)*csch(x)**(m + S(1)), x), x, acsch(c*x)))
    rubi.add(rule603)

    pattern604 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule604 = ReplacementRule(pattern604, lambda x, b, m, c, a, n : Int(x**m*(a + b*asech(c*x))**n, x))
    rubi.add(rule604)

    pattern605 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule605 = ReplacementRule(pattern605, lambda x, b, m, c, a, n : Int(x**m*(a + b*acsch(c*x))**n, x))
    rubi.add(rule605)

    pattern606 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With606(p, x, b, c, a, d, e):
        u = IntHide((d + e*x**S(2))**p, x)
        return b*sqrt(c*x + S(1))*sqrt(1/(c*x + S(1)))*Int(SimplifyIntegrand(u/(x*sqrt(-c*x + S(1))*sqrt(c*x + S(1))), x), x) + Dist(a + b*asech(c*x), u, x)
    rule606 = ReplacementRule(pattern606, lambda p, x, b, c, a, d, e : With606(p, x, b, c, a, d, e))
    rubi.add(rule606)

    pattern607 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With607(p, x, b, c, a, d, e):
        u = IntHide((d + e*x**S(2))**p, x)
        return -b*c*x*Int(SimplifyIntegrand(u/(x*sqrt(-c**S(2)*x**S(2) + S(-1))), x), x)/sqrt(-c**S(2)*x**S(2)) + Dist(a + b*acsch(c*x), u, x)
    rule607 = ReplacementRule(pattern607, lambda p, x, b, c, a, d, e : With607(p, x, b, c, a, d, e))
    rubi.add(rule607)

    pattern608 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)))
    rule608 = ReplacementRule(pattern608, lambda p, x, b, c, a, d, e, n : -Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule608)

    pattern609 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)))
    rule609 = ReplacementRule(pattern609, lambda p, x, b, c, a, d, e, n : -Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule609)

    pattern610 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule610 = ReplacementRule(pattern610, lambda p, x, b, c, a, d, e, n : -sqrt(x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule610)

    pattern611 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule611 = ReplacementRule(pattern611, lambda p, x, b, c, a, d, e, n : -sqrt(x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule611)

    pattern612 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d, e: Not(Negative(d) & PositiveQ(e))))
    rule612 = ReplacementRule(pattern612, lambda p, x, b, c, a, d, e, n : -sqrt(d + e*x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*sqrt(d/x**S(2) + e)))
    rubi.add(rule612)

    pattern613 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d, e: Not(Negative(d) & PositiveQ(e))))
    rule613 = ReplacementRule(pattern613, lambda p, x, b, c, a, d, e, n : -sqrt(d + e*x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*sqrt(d/x**S(2) + e)))
    rubi.add(rule613)

    pattern614 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule614 = ReplacementRule(pattern614, lambda p, x, b, c, a, d, e, n : Int((a + b*asech(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule614)

    pattern615 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule615 = ReplacementRule(pattern615, lambda p, x, b, c, a, d, e, n : Int((a + b*acsch(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule615)

    pattern616 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule616 = ReplacementRule(pattern616, lambda p, x, b, c, a, d, e : b*sqrt(c*x + S(1))*sqrt(1/(c*x + S(1)))*Int((d + e*x**S(2))**(p + S(1))/(x*sqrt(-c*x + S(1))*sqrt(c*x + S(1))), x)/(S(2)*e*(p + S(1))) + (a + b*asech(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule616)

    pattern617 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule617 = ReplacementRule(pattern617, lambda p, x, b, c, a, d, e : -b*c*x*Int((d + e*x**S(2))**(p + S(1))/(x*sqrt(-c**S(2)*x**S(2) + S(-1))), x)/(S(2)*e*sqrt(-c**S(2)*x**S(2))*(p + S(1))) + (a + b*acsch(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule617)

    pattern618 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & Not(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))), )
    def With618(p, x, b, m, c, a, d, e):
        u = IntHide(x**m*(d + e*x**S(2))**p, x)
        return b*sqrt(c*x + S(1))*sqrt(1/(c*x + S(1)))*Int(SimplifyIntegrand(u/(x*sqrt(-c*x + S(1))*sqrt(c*x + S(1))), x), x) + Dist(a + b*asech(c*x), u, x)
    rule618 = ReplacementRule(pattern618, lambda p, x, b, m, c, a, d, e : With618(p, x, b, m, c, a, d, e))
    rubi.add(rule618)

    pattern619 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & Not(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))), )
    def With619(p, x, b, m, c, a, d, e):
        u = IntHide(x**m*(d + e*x**S(2))**p, x)
        return -b*c*x*Int(SimplifyIntegrand(u/(x*sqrt(-c**S(2)*x**S(2) + S(-1))), x), x)/sqrt(-c**S(2)*x**S(2)) + Dist(a + b*acsch(c*x), u, x)
    rule619 = ReplacementRule(pattern619, lambda p, x, b, m, c, a, d, e : With619(p, x, b, m, c, a, d, e))
    rubi.add(rule619)

    pattern620 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, m: IntegersQ(m, p)))
    rule620 = ReplacementRule(pattern620, lambda p, x, b, m, c, a, d, e, n : -Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule620)

    pattern621 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, m: IntegersQ(m, p)))
    rule621 = ReplacementRule(pattern621, lambda p, x, b, m, c, a, d, e, n : -Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule621)

    pattern622 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule622 = ReplacementRule(pattern622, lambda p, x, b, m, c, a, d, e, n : -sqrt(x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule622)

    pattern623 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule623 = ReplacementRule(pattern623, lambda p, x, b, m, c, a, d, e, n : -sqrt(x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule623)

    pattern624 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, e, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d, e: Not(Negative(d) & PositiveQ(e))))
    rule624 = ReplacementRule(pattern624, lambda p, x, b, m, c, a, d, e, n : -sqrt(d + e*x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*sqrt(d/x**S(2) + e)))
    rubi.add(rule624)

    pattern625 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda c, d, e: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d, e: Not(Negative(d) & PositiveQ(e))))
    rule625 = ReplacementRule(pattern625, lambda p, x, b, m, c, a, d, e, n : -sqrt(d + e*x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*sqrt(d/x**S(2) + e)))
    rubi.add(rule625)

    pattern626 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule626 = ReplacementRule(pattern626, lambda p, x, b, m, c, a, d, e, n : Int(x**m*(a + b*asech(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule626)

    pattern627 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule627 = ReplacementRule(pattern627, lambda p, x, b, m, c, a, d, e, n : Int(x**m*(a + b*acsch(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule627)

    pattern628 = Pattern(Integral(asech(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule628 = ReplacementRule(pattern628, lambda x, a, b : Int(sqrt((-a - b*x + S(1))/(a + b*x + S(1)))/(-a - b*x + S(1)), x) + (a + b*x)*asech(a + b*x)/b)
    rubi.add(rule628)

    pattern629 = Pattern(Integral(acsch(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule629 = ReplacementRule(pattern629, lambda x, a, b : Int(S(1)/(sqrt(S(1) + (a + b*x)**(S(-2)))*(a + b*x)), x) + (a + b*x)*acsch(a + b*x)/b)
    rubi.add(rule629)

    pattern630 = Pattern(Integral(asech(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule630 = ReplacementRule(pattern630, lambda x, a, n, b : -Subst(Int(x**n*tanh(x)*sech(x), x), x, asech(a + b*x))/b)
    rubi.add(rule630)

    pattern631 = Pattern(Integral(acsch(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule631 = ReplacementRule(pattern631, lambda x, a, n, b : -Subst(Int(x**n*coth(x)*csch(x), x), x, acsch(a + b*x))/b)
    rubi.add(rule631)

    pattern632 = Pattern(Integral(asech(a_ + x_*WC('b', S(1)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule632 = ReplacementRule(pattern632, lambda x, a, b : -PolyLog(S(2), (-sqrt(-a**S(2) + S(1)) + S(1))*exp(-asech(a + b*x))/a) - PolyLog(S(2), (sqrt(-a**S(2) + S(1)) + S(1))*exp(-asech(a + b*x))/a) + PolyLog(S(2), -exp(-S(2)*asech(a + b*x)))/S(2) + log(S(1) - (-sqrt(-a**S(2) + S(1)) + S(1))*exp(-asech(a + b*x))/a)*asech(a + b*x) + log(S(1) - (sqrt(-a**S(2) + S(1)) + S(1))*exp(-asech(a + b*x))/a)*asech(a + b*x) - log(S(1) + exp(-S(2)*asech(a + b*x)))*asech(a + b*x))
    rubi.add(rule632)

    pattern633 = Pattern(Integral(acsch(a_ + x_*WC('b', S(1)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule633 = ReplacementRule(pattern633, lambda x, a, b : PolyLog(S(2), (-sqrt(a**S(2) + S(1)) + S(-1))*exp(acsch(a + b*x))/a) + PolyLog(S(2), (sqrt(a**S(2) + S(1)) + S(-1))*exp(acsch(a + b*x))/a) + PolyLog(S(2), exp(-S(2)*acsch(a + b*x)))/S(2) + log(S(1) + (-sqrt(a**S(2) + S(1)) + S(1))*exp(acsch(a + b*x))/a)*acsch(a + b*x) + log(S(1) + (sqrt(a**S(2) + S(1)) + S(1))*exp(acsch(a + b*x))/a)*acsch(a + b*x) - log(S(1) - exp(-S(2)*acsch(a + b*x)))*acsch(a + b*x) - acsch(a + b*x)**S(2))
    rubi.add(rule633)

    pattern634 = Pattern(Integral(x_**WC('m', S(1))*asech(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule634 = ReplacementRule(pattern634, lambda x, a, b, m : b**(-m + S(-1))*(b**(m + S(1))*x**(m + S(1)) - (-a)**(m + S(1)))*asech(a + b*x)/(m + S(1)) + b**(-m + S(-1))*Subst(Int(x**(-m + S(-1))*((-a*x)**(m + S(1)) - (-a*x + S(1))**(m + S(1)))/(sqrt(x + S(-1))*sqrt(x + S(1))), x), x, 1/(a + b*x))/(m + S(1)))
    rubi.add(rule634)

    pattern635 = Pattern(Integral(x_**WC('m', S(1))*acsch(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule635 = ReplacementRule(pattern635, lambda x, a, b, m : b**(-m + S(-1))*(b**(m + S(1))*x**(m + S(1)) - (-a)**(m + S(1)))*acsch(a + b*x)/(m + S(1)) + b**(-m + S(-1))*Subst(Int(x**(-m + S(-1))*((-a*x)**(m + S(1)) - (-a*x + S(1))**(m + S(1)))/sqrt(x**S(2) + S(1)), x), x, 1/(a + b*x))/(m + S(1)))
    rubi.add(rule635)

    pattern636 = Pattern(Integral(x_**WC('m', S(1))*asech(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule636 = ReplacementRule(pattern636, lambda x, b, m, a, n : -b**(-m + S(-1))*Subst(Int(x**n*(-a + sech(x))**m*tanh(x)*sech(x), x), x, asech(a + b*x)))
    rubi.add(rule636)

    pattern637 = Pattern(Integral(x_**WC('m', S(1))*acsch(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule637 = ReplacementRule(pattern637, lambda x, b, m, a, n : -b**(-m + S(-1))*Subst(Int(x**n*(-a + csch(x))**m*coth(x)*csch(x), x), x, acsch(a + b*x)))
    rubi.add(rule637)

    pattern638 = Pattern(Integral(WC('u', S(1))*asech(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule638 = ReplacementRule(pattern638, lambda x, b, m, c, a, u, n : Int(u*acosh(a/c + b*x**n/c)**m, x))
    rubi.add(rule638)

    pattern639 = Pattern(Integral(WC('u', S(1))*acsch(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule639 = ReplacementRule(pattern639, lambda x, b, m, c, a, u, n : Int(u*asinh(a/c + b*x**n/c)**m, x))
    rubi.add(rule639)

    pattern640 = Pattern(Integral(exp(asech(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)))
    rule640 = ReplacementRule(pattern640, lambda x, a : x*exp(asech(a*x)) + Int(sqrt((-a*x + S(1))/(a*x + S(1)))/(x*(-a*x + S(1))), x)/a + log(x)/a)
    rubi.add(rule640)

    pattern641 = Pattern(Integral(exp(asech(x_**p_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule641 = ReplacementRule(pattern641, lambda x, a, p : x*exp(asech(a*x**p)) + p*sqrt(a*x**p + S(1))*sqrt(1/(a*x**p + S(1)))*Int(x**(-p)/(sqrt(-a*x**p + S(1))*sqrt(a*x**p + S(1))), x)/a + p*Int(x**(-p), x)/a)
    rubi.add(rule641)

    pattern642 = Pattern(Integral(exp(acsch(x_**WC('p', S(1))*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule642 = ReplacementRule(pattern642, lambda x, a, p : Int(sqrt(S(1) + x**(-S(2)*p)/a**S(2)), x) + Int(x**(-p), x)/a)
    rubi.add(rule642)

    pattern643 = Pattern(Integral(exp(WC('n', S(1))*asech(u_)), x_), CustomConstraint(lambda n: IntegerQ(n)))
    rule643 = ReplacementRule(pattern643, lambda x, u, n : Int((sqrt((-u + S(1))/(u + S(1))) + sqrt((-u + S(1))/(u + S(1)))/u + 1/u)**n, x))
    rubi.add(rule643)

    pattern644 = Pattern(Integral(exp(WC('n', S(1))*acsch(u_)), x_), CustomConstraint(lambda n: IntegerQ(n)))
    rule644 = ReplacementRule(pattern644, lambda x, u, n : Int((sqrt(S(1) + u**(S(-2))) + 1/u)**n, x))
    rubi.add(rule644)

    pattern645 = Pattern(Integral(exp(asech(x_**WC('p', S(1))*WC('a', S(1))))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule645 = ReplacementRule(pattern645, lambda x, a, p : sqrt(a*x**p + S(1))*sqrt(1/(a*x**p + S(1)))*Int(x**(-p + S(-1))*sqrt(-a*x**p + S(1))*sqrt(a*x**p + S(1)), x)/a - x**(-p)/(a*p))
    rubi.add(rule645)

    pattern646 = Pattern(Integral(x_**WC('m', S(1))*exp(asech(x_**WC('p', S(1))*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule646 = ReplacementRule(pattern646, lambda x, a, p, m : x**(m + S(1))*exp(asech(a*x**p))/(m + S(1)) + p*sqrt(a*x**p + S(1))*sqrt(1/(a*x**p + S(1)))*Int(x**(m - p)/(sqrt(-a*x**p + S(1))*sqrt(a*x**p + S(1))), x)/(a*(m + S(1))) + p*Int(x**(m - p), x)/(a*(m + S(1))))
    rubi.add(rule646)

    pattern647 = Pattern(Integral(x_**WC('m', S(1))*exp(acsch(x_**WC('p', S(1))*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule647 = ReplacementRule(pattern647, lambda x, a, p, m : Int(x**m*sqrt(S(1) + x**(-S(2)*p)/a**S(2)), x) + Int(x**(m - p), x)/a)
    rubi.add(rule647)

    pattern648 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*asech(u_)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule648 = ReplacementRule(pattern648, lambda x, u, n, m : Int(x**m*(sqrt((-u + S(1))/(u + S(1))) + sqrt((-u + S(1))/(u + S(1)))/u + 1/u)**n, x))
    rubi.add(rule648)

    pattern649 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*acsch(u_)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule649 = ReplacementRule(pattern649, lambda x, u, n, m : Int(x**m*(sqrt(S(1) + u**(S(-2))) + 1/u)**n, x))
    rubi.add(rule649)

    pattern650 = Pattern(Integral(asech(u_), x_), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda x, u: Not(FunctionOfExponentialQ(u, x))))
    rule650 = ReplacementRule(pattern650, lambda x, u : x*asech(u) + sqrt(-u**S(2) + S(1))*Int(SimplifyIntegrand(x*D(u, x)/(u*sqrt(-u**S(2) + S(1))), x), x)/(u*sqrt(S(-1) + 1/u)*sqrt(S(1) + 1/u)))
    rubi.add(rule650)

    pattern651 = Pattern(Integral(acsch(u_), x_), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda x, u: Not(FunctionOfExponentialQ(u, x))))
    rule651 = ReplacementRule(pattern651, lambda x, u : -u*Int(SimplifyIntegrand(x*D(u, x)/(u*sqrt(-u**S(2) + S(-1))), x), x)/sqrt(-u**S(2)) + x*acsch(u))
    rubi.add(rule651)

    pattern652 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(u_)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda x, m, c, d, u: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda x, u: Not(FunctionOfExponentialQ(u, x))))
    rule652 = ReplacementRule(pattern652, lambda x, b, m, c, a, d, u : b*sqrt(-u**S(2) + S(1))*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u*sqrt(-u**S(2) + S(1))), x), x)/(d*u*sqrt(S(-1) + 1/u)*sqrt(S(1) + 1/u)*(m + S(1))) + (a + b*asech(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule652)

    pattern653 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(u_)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda x, m, c, d, u: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda x, u: Not(FunctionOfExponentialQ(u, x))))
    rule653 = ReplacementRule(pattern653, lambda x, b, m, c, a, d, u : -b*u*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u*sqrt(-u**S(2) + S(-1))), x), x)/(d*sqrt(-u**S(2))*(m + S(1))) + (a + b*acsch(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule653)

    pattern654 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*asech(u_)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda a, u, x, b, w: InverseFunctionFreeQ(w, x)))
    def With654(x, b, v, a, u):
        w = IntHide(v, x)
        return b*sqrt(-u**S(2) + S(1))*Int(SimplifyIntegrand(w*D(u, x)/(u*sqrt(-u**S(2) + S(1))), x), x)/(u*sqrt(S(-1) + 1/u)*sqrt(S(1) + 1/u)) + Dist(a + b*asech(u), w, x)
    rule654 = ReplacementRule(pattern654, lambda x, b, v, a, u : With654(x, b, v, a, u))
    rubi.add(rule654)

    pattern655 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acsch(u_)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, u: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda a, u, x, b, w: InverseFunctionFreeQ(w, x)))
    def With655(x, b, v, a, u):
        w = IntHide(v, x)
        return -b*u*Int(SimplifyIntegrand(w*D(u, x)/(u*sqrt(-u**S(2) + S(-1))), x), x)/sqrt(-u**S(2)) + Dist(a + b*acsch(u), w, x)
    rule655 = ReplacementRule(pattern655, lambda x, b, v, a, u : With655(x, b, v, a, u))
    rubi.add(rule655)

    return rubi
