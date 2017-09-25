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

def inverse_trig(rubi):

    pattern1 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule1 = ReplacementRule(pattern1, lambda b, c, n, a, x : -b*c*n*Int(x*(a + b*asin(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x) + x*(a + b*asin(c*x))**n)
    rubi.add(rule1)

    pattern2 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule2 = ReplacementRule(pattern2, lambda b, c, n, a, x : b*c*n*Int(x*(a + b*acos(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x) + x*(a + b*acos(c*x))**n)
    rubi.add(rule2)

    pattern3 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule3 = ReplacementRule(pattern3, lambda b, c, n, a, x : c*Int(x*(a + b*asin(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*(n + S(1))) + (a + b*asin(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule3)

    pattern4 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule4 = ReplacementRule(pattern4, lambda b, c, n, a, x : -c*Int(x*(a + b*acos(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*(n + S(1))) - (a + b*acos(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule4)

    pattern5 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule5 = ReplacementRule(pattern5, lambda b, c, n, a, x : Subst(Int(x**n*cos(a/b - x/b), x), x, a + b*asin(c*x))/(b*c))
    rubi.add(rule5)

    pattern6 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule6 = ReplacementRule(pattern6, lambda b, c, n, a, x : Subst(Int(x**n*sin(a/b - x/b), x), x, a + b*acos(c*x))/(b*c))
    rubi.add(rule6)

    pattern7 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule7 = ReplacementRule(pattern7, lambda b, c, n, a, x : Subst(Int((a + b*x)**n/tan(x), x), x, asin(c*x)))
    rubi.add(rule7)

    pattern8 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule8 = ReplacementRule(pattern8, lambda b, c, n, a, x : -Subst(Int((a + b*x)**n/cot(x), x), x, acos(c*x)))
    rubi.add(rule8)

    pattern9 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule9 = ReplacementRule(pattern9, lambda b, c, m, d, n, a, x : -b*c*n*Int((d*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(d*(m + S(1))) + (d*x)**(m + S(1))*(a + b*asin(c*x))**n/(d*(m + S(1))))
    rubi.add(rule9)

    pattern10 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule10 = ReplacementRule(pattern10, lambda b, c, m, d, n, a, x : b*c*n*Int((d*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(d*(m + S(1))) + (d*x)**(m + S(1))*(a + b*acos(c*x))**n/(d*(m + S(1))))
    rubi.add(rule10)

    pattern11 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule11 = ReplacementRule(pattern11, lambda b, c, m, n, a, x : -b*c*n*Int(x**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*asin(c*x))**n/(m + S(1)))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule12 = ReplacementRule(pattern12, lambda b, c, m, n, a, x : b*c*n*Int(x**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*acos(c*x))**n/(m + S(1)))
    rubi.add(rule12)

    pattern13 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Inequality(S(-2), LessEqual, n, Less, S(-1))))
    rule13 = ReplacementRule(pattern13, lambda b, c, m, n, a, x : -c**(-m + S(-1))*Subst(Int(ExpandTrigReduce((a + b*x)**(n + S(1)), (m - (m + S(1))*sin(x)**S(2))*sin(x)**(m + S(-1)), x), x), x, asin(c*x))/(b*(n + S(1))) + x**m*(a + b*asin(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule13)

    pattern14 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Inequality(S(-2), LessEqual, n, Less, S(-1))))
    rule14 = ReplacementRule(pattern14, lambda b, c, m, n, a, x : -c**(-m + S(-1))*Subst(Int(ExpandTrigReduce((a + b*x)**(n + S(1)), (m - (m + S(1))*cos(x)**S(2))*cos(x)**(m + S(-1)), x), x), x, acos(c*x))/(b*(n + S(1))) - x**m*(a + b*acos(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule14)

    pattern15 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-2))))
    rule15 = ReplacementRule(pattern15, lambda b, c, m, n, a, x : c*(m + S(1))*Int(x**(m + S(1))*(a + b*asin(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*(n + S(1))) - m*Int(x**(m + S(-1))*(a + b*asin(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*c*(n + S(1))) + x**m*(a + b*asin(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule15)

    pattern16 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-2))))
    rule16 = ReplacementRule(pattern16, lambda b, c, m, n, a, x : -c*(m + S(1))*Int(x**(m + S(1))*(a + b*acos(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*(n + S(1))) + m*Int(x**(m + S(-1))*(a + b*acos(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(b*c*(n + S(1))) - x**m*(a + b*acos(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule16)

    pattern17 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule17 = ReplacementRule(pattern17, lambda b, c, m, n, a, x : c**(-m + S(-1))*Subst(Int((a + b*x)**n*sin(x)**m*cos(x), x), x, asin(c*x)))
    rubi.add(rule17)

    pattern18 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule18 = ReplacementRule(pattern18, lambda b, c, m, n, a, x : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*sin(x)*cos(x)**m, x), x, acos(c*x)))
    rubi.add(rule18)

    pattern19 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule19 = ReplacementRule(pattern19, lambda b, c, m, d, n, a, x : Int((d*x)**m*(a + b*asin(c*x))**n, x))
    rubi.add(rule19)

    pattern20 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule20 = ReplacementRule(pattern20, lambda b, c, m, d, n, a, x : Int((d*x)**m*(a + b*acos(c*x))**n, x))
    rubi.add(rule20)

    pattern21 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule21 = ReplacementRule(pattern21, lambda b, e, c, d, a, x : log(a + b*asin(c*x))/(b*c*sqrt(d)))
    rubi.add(rule21)

    pattern22 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule22 = ReplacementRule(pattern22, lambda b, e, c, d, a, x : -log(a + b*acos(c*x))/(b*c*sqrt(d)))
    rubi.add(rule22)

    pattern23 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule23 = ReplacementRule(pattern23, lambda b, e, c, d, n, a, x : (a + b*asin(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule23)

    pattern24 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule24 = ReplacementRule(pattern24, lambda b, e, c, d, n, a, x : -(a + b*acos(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule24)

    pattern25 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule25 = ReplacementRule(pattern25, lambda b, e, c, d, n, a, x : sqrt(-c**S(2)*x**S(2) + S(1))*Int((a + b*asin(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule25)

    pattern26 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule26 = ReplacementRule(pattern26, lambda b, e, c, d, n, a, x : sqrt(-c**S(2)*x**S(2) + S(1))*Int((a + b*acos(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule26)

    pattern27 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), )
    def With27(b, e, c, d, a, p, x):
        u = IntHide((d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asin(c*x), u, x)
    rule27 = ReplacementRule(pattern27, lambda b, e, c, d, a, p, x : With27(b, e, c, d, a, p, x))
    rubi.add(rule27)

    pattern28 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), )
    def With28(b, e, c, d, a, p, x):
        u = IntHide((d + e*x**S(2))**p, x)
        return b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acos(c*x), u, x)
    rule28 = ReplacementRule(pattern28, lambda b, e, c, d, a, p, x : With28(b, e, c, d, a, p, x))
    rubi.add(rule28)

    pattern29 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule29 = ReplacementRule(pattern29, lambda b, e, c, d, n, a, x : -b*c*n*sqrt(d + e*x**S(2))*Int(x*(a + b*asin(c*x))**(n + S(-1)), x)/(S(2)*sqrt(-c**S(2)*x**S(2) + S(1))) + x*(a + b*asin(c*x))**n*sqrt(d + e*x**S(2))/S(2) + sqrt(d + e*x**S(2))*Int((a + b*asin(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(S(2)*sqrt(-c**S(2)*x**S(2) + S(1))))
    rubi.add(rule29)

    pattern30 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule30 = ReplacementRule(pattern30, lambda b, e, c, d, n, a, x : b*c*n*sqrt(d + e*x**S(2))*Int(x*(a + b*acos(c*x))**(n + S(-1)), x)/(S(2)*sqrt(-c**S(2)*x**S(2) + S(1))) + x*(a + b*acos(c*x))**n*sqrt(d + e*x**S(2))/S(2) + sqrt(d + e*x**S(2))*Int((a + b*acos(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(S(2)*sqrt(-c**S(2)*x**S(2) + S(1))))
    rubi.add(rule30)

    pattern31 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))))
    rule31 = ReplacementRule(pattern31, lambda b, e, c, d, n, a, p, x : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(S(2)*p + S(1)) + S(2)*d*p*Int((a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*asin(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule31)

    pattern32 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))))
    rule32 = ReplacementRule(pattern32, lambda b, e, c, d, n, a, p, x : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(S(2)*p + S(1)) + S(2)*d*p*Int((a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*acos(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule32)

    pattern33 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda d: PositiveQ(d)))
    rule33 = ReplacementRule(pattern33, lambda b, e, c, d, n, a, x : -b*c*n*Int(x*(a + b*asin(c*x))**(n + S(-1))/(d + e*x**S(2)), x)/sqrt(d) + x*(a + b*asin(c*x))**n/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule33)

    pattern34 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda d: PositiveQ(d)))
    rule34 = ReplacementRule(pattern34, lambda b, e, c, d, n, a, x : b*c*n*Int(x*(a + b*acos(c*x))**(n + S(-1))/(d + e*x**S(2)), x)/sqrt(d) + x*(a + b*acos(c*x))**n/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule34)

    pattern35 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule35 = ReplacementRule(pattern35, lambda b, e, c, d, n, a, x : -b*c*n*sqrt(-c**S(2)*x**S(2) + S(1))*Int(x*(a + b*asin(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x)/(d*sqrt(d + e*x**S(2))) + x*(a + b*asin(c*x))**n/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule35)

    pattern36 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule36 = ReplacementRule(pattern36, lambda b, e, c, d, n, a, x : b*c*n*sqrt(-c**S(2)*x**S(2) + S(1))*Int(x*(a + b*acos(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x)/(d*sqrt(d + e*x**S(2))) + x*(a + b*acos(c*x))**n/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule36)

    pattern37 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule37 = ReplacementRule(pattern37, lambda b, e, c, d, n, a, p, x : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*(p + S(1))) - x*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule37)

    pattern38 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule38 = ReplacementRule(pattern38, lambda b, e, c, d, n, a, p, x : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*(p + S(1))) - x*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule38)

    pattern39 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule39 = ReplacementRule(pattern39, lambda b, e, c, d, n, a, x : Subst(Int((a + b*x)**n*sec(x), x), x, asin(c*x))/(c*d))
    rubi.add(rule39)

    pattern40 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule40 = ReplacementRule(pattern40, lambda b, e, c, d, n, a, x : -Subst(Int((a + b*x)**n*csc(x), x), x, acos(c*x))/(c*d))
    rubi.add(rule40)

    pattern41 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule41 = ReplacementRule(pattern41, lambda b, e, c, d, n, a, p, x : c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(S(2)*p + S(1))*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*asin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*(n + S(1))) + (a + b*asin(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule41)

    pattern42 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule42 = ReplacementRule(pattern42, lambda b, e, c, d, n, a, p, x : -c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(S(2)*p + S(1))*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x*(a + b*acos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*(n + S(1))) - (a + b*acos(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule42)

    pattern43 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)), CustomConstraint(lambda p, d: IntegerQ(p) | PositiveQ(d)))
    rule43 = ReplacementRule(pattern43, lambda b, e, c, d, n, a, p, x : d**p*Subst(Int((a + b*x)**n*cos(x)**(S(2)*p + S(1)), x), x, asin(c*x))/c)
    rubi.add(rule43)

    pattern44 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)), CustomConstraint(lambda p, d: IntegerQ(p) | PositiveQ(d)))
    rule44 = ReplacementRule(pattern44, lambda b, e, c, d, n, a, p, x : -d**p*Subst(Int((a + b*x)**n*sin(x)**(S(2)*p + S(1)), x), x, acos(c*x))/c)
    rubi.add(rule44)

    pattern45 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)), CustomConstraint(lambda p, d: Not(IntegerQ(p) | PositiveQ(d))))
    rule45 = ReplacementRule(pattern45, lambda b, e, c, d, n, a, p, x : d**(p + S(-1)/2)*sqrt(d + e*x**S(2))*Int((a + b*asin(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x)/sqrt(-c**S(2)*x**S(2) + S(1)))
    rubi.add(rule45)

    pattern46 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)), CustomConstraint(lambda p, d: Not(IntegerQ(p) | PositiveQ(d))))
    rule46 = ReplacementRule(pattern46, lambda b, e, c, d, n, a, p, x : d**(p + S(-1)/2)*sqrt(d + e*x**S(2))*Int((a + b*acos(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x)/sqrt(-c**S(2)*x**S(2) + S(1)))
    rubi.add(rule46)

    pattern47 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With47(b, e, c, d, a, p, x):
        u = IntHide((d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asin(c*x), u, x)
    rule47 = ReplacementRule(pattern47, lambda b, e, c, d, a, p, x : With47(b, e, c, d, a, p, x))
    rubi.add(rule47)

    pattern48 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With48(b, e, c, d, a, p, x):
        u = IntHide((d + e*x**S(2))**p, x)
        return b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acos(c*x), u, x)
    rule48 = ReplacementRule(pattern48, lambda b, e, c, d, a, p, x : With48(b, e, c, d, a, p, x))
    rubi.add(rule48)

    pattern49 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n, p: PositiveIntegerQ(n) | Greater(p, S(0))))
    rule49 = ReplacementRule(pattern49, lambda b, e, c, d, n, a, p, x : Int(ExpandIntegrand((a + b*asin(c*x))**n, (d + e*x**S(2))**p, x), x))
    rubi.add(rule49)

    pattern50 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n, p: PositiveIntegerQ(n) | Greater(p, S(0))))
    rule50 = ReplacementRule(pattern50, lambda b, e, c, d, n, a, p, x : Int(ExpandIntegrand((a + b*acos(c*x))**n, (d + e*x**S(2))**p, x), x))
    rubi.add(rule50)

    pattern51 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule51 = ReplacementRule(pattern51, lambda b, e, c, d, n, a, p, x : Int((a + b*asin(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule51)

    pattern52 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule52 = ReplacementRule(pattern52, lambda b, e, c, d, n, a, p, x : Int((a + b*acos(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule52)

    pattern53 = Pattern(Integral((d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda f, e, g, d: ZeroQ(d*g + e*f)), CustomConstraint(lambda f, g, c: ZeroQ(c**S(2)*f**S(2) - g**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule53 = ReplacementRule(pattern53, lambda b, e, c, d, n, g, a, f, p, x : (d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p))*Int((a + b*asin(c*x))**n*(d*f + e*g*x**S(2))**p, x))
    rubi.add(rule53)

    pattern54 = Pattern(Integral((d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda f, e, g, d: ZeroQ(d*g + e*f)), CustomConstraint(lambda f, g, c: ZeroQ(c**S(2)*f**S(2) - g**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule54 = ReplacementRule(pattern54, lambda b, e, c, d, n, g, a, f, p, x : (d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p))*Int((a + b*acos(c*x))**n*(d*f + e*g*x**S(2))**p, x))
    rubi.add(rule54)

    pattern55 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule55 = ReplacementRule(pattern55, lambda b, e, c, d, n, a, x : -Subst(Int((a + b*x)**n*tan(x), x), x, asin(c*x))/e)
    rubi.add(rule55)

    pattern56 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule56 = ReplacementRule(pattern56, lambda b, e, c, d, n, a, x : Subst(Int((a + b*x)**n*cot(x), x), x, acos(c*x))/e)
    rubi.add(rule56)

    pattern57 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule57 = ReplacementRule(pattern57, lambda b, e, c, d, n, a, p, x : b*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) + (a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule57)

    pattern58 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule58 = ReplacementRule(pattern58, lambda b, e, c, d, n, a, p, x : -b*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) + (a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule58)

    pattern59 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule59 = ReplacementRule(pattern59, lambda b, e, c, d, n, a, x : Subst(Int((a + b*x)**n/(sin(x)*cos(x)), x), x, asin(c*x))/d)
    rubi.add(rule59)

    pattern60 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule60 = ReplacementRule(pattern60, lambda b, e, c, d, n, a, x : -Subst(Int((a + b*x)**n/(sin(x)*cos(x)), x), x, acos(c*x))/d)
    rubi.add(rule60)

    pattern61 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule61 = ReplacementRule(pattern61, lambda b, e, c, m, d, n, a, f, p, x : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + (f*x)**(m + S(1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule61)

    pattern62 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule62 = ReplacementRule(pattern62, lambda b, e, c, m, d, n, a, f, p, x : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule62)

    pattern63 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule63 = ReplacementRule(pattern63, lambda b, e, c, d, a, p, x : -b*c*d**p*Int((-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(S(2)*p) + d*Int((a + b*asin(c*x))*(d + e*x**S(2))**(p + S(-1))/x, x) + (a + b*asin(c*x))*(d + e*x**S(2))**p/(S(2)*p))
    rubi.add(rule63)

    pattern64 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule64 = ReplacementRule(pattern64, lambda b, e, c, d, a, p, x : b*c*d**p*Int((-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(S(2)*p) + d*Int((a + b*acos(c*x))*(d + e*x**S(2))**(p + S(-1))/x, x) + (a + b*acos(c*x))*(d + e*x**S(2))**p/(S(2)*p))
    rubi.add(rule64)

    pattern65 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2) + S(1)/2)))
    rule65 = ReplacementRule(pattern65, lambda b, e, c, m, d, a, f, p, x : -b*c*d**p*Int((f*x)**(m + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*asin(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*asin(c*x))*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule65)

    pattern66 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2) + S(1)/2)))
    rule66 = ReplacementRule(pattern66, lambda b, e, c, m, d, a, f, p, x : b*c*d**p*Int((f*x)**(m + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*acos(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acos(c*x))*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule66)

    pattern67 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), )
    def With67(b, e, c, m, d, a, f, p, x):
        u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asin(c*x), u, x)
    rule67 = ReplacementRule(pattern67, lambda b, e, c, m, d, a, f, p, x : With67(b, e, c, m, d, a, f, p, x))
    rubi.add(rule67)

    pattern68 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p)), )
    def With68(b, e, c, m, d, a, f, p, x):
        u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
        return b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acos(c*x), u, x)
    rule68 = ReplacementRule(pattern68, lambda b, e, c, m, d, a, f, p, x : With68(b, e, c, m, d, a, f, p, x))
    rubi.add(rule68)

    pattern69 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda d: PositiveQ(d)), )
    def With69(b, e, c, m, d, a, p, x):
        u = IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x)
        return -b*c*d**p*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(d**p*(a + b*asin(c*x)), u, x)
    rule69 = ReplacementRule(pattern69, lambda b, e, c, m, d, a, p, x : With69(b, e, c, m, d, a, p, x))
    rubi.add(rule69)

    pattern70 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)), CustomConstraint(lambda p: Unequal(p, S(-1)/2)), CustomConstraint(lambda d: PositiveQ(d)), )
    def With70(b, e, c, m, d, a, p, x):
        u = IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x)
        return b*c*d**p*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(d**p*(a + b*acos(c*x)), u, x)
    rule70 = ReplacementRule(pattern70, lambda b, e, c, m, d, a, p, x : With70(b, e, c, m, d, a, p, x))
    rubi.add(rule70)

    pattern71 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)), )
    def With71(b, e, c, m, d, a, p, x):
        u = IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x)
        return -b*c*d**(p + S(-1)/2)*sqrt(d + e*x**S(2))*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x)/sqrt(-c**S(2)*x**S(2) + S(1)) + (a + b*asin(c*x))*Int(x**m*(d + e*x**S(2))**p, x)
    rule71 = ReplacementRule(pattern71, lambda b, e, c, m, d, a, p, x : With71(b, e, c, m, d, a, p, x))
    rubi.add(rule71)

    pattern72 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda p, m: PositiveIntegerQ(m/S(2) + S(1)/2) | NegativeIntegerQ(m/S(2) + p + S(3)/2)), )
    def With72(b, e, c, m, d, a, p, x):
        u = IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x)
        return b*c*d**(p + S(-1)/2)*sqrt(d + e*x**S(2))*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x)/sqrt(-c**S(2)*x**S(2) + S(1)) + (a + b*acos(c*x))*Int(x**m*(d + e*x**S(2))**p, x)
    rule72 = ReplacementRule(pattern72, lambda b, e, c, m, d, a, p, x : With72(b, e, c, m, d, a, p, x))
    rubi.add(rule72)

    pattern73 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule73 = ReplacementRule(pattern73, lambda b, e, c, m, d, n, a, f, x : -b*c*n*sqrt(d + e*x**S(2))*Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1)), x)/(f*(m + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))) + c**S(2)*sqrt(d + e*x**S(2))*Int((f*x)**(m + S(2))*(a + b*asin(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(f**S(2)*(m + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))) + (f*x)**(m + S(1))*(a + b*asin(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(1))))
    rubi.add(rule73)

    pattern74 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule74 = ReplacementRule(pattern74, lambda b, e, c, m, d, n, a, f, x : b*c*n*sqrt(d + e*x**S(2))*Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1)), x)/(f*(m + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))) + c**S(2)*sqrt(d + e*x**S(2))*Int((f*x)**(m + S(2))*(a + b*acos(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(f**S(2)*(m + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))) + (f*x)**(m + S(1))*(a + b*acos(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(1))))
    rubi.add(rule74)

    pattern75 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule75 = ReplacementRule(pattern75, lambda b, e, c, m, d, n, a, f, p, x : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule75)

    pattern76 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule76 = ReplacementRule(pattern76, lambda b, e, c, m, d, n, a, f, p, x : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(1))) - S(2)*e*p*Int((f*x)**(m + S(2))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(1))))
    rubi.add(rule76)

    pattern77 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Not(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule77 = ReplacementRule(pattern77, lambda b, e, c, m, d, n, a, f, x : -b*c*n*sqrt(d + e*x**S(2))*Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1)), x)/(f*(m + S(2))*sqrt(-c**S(2)*x**S(2) + S(1))) + sqrt(d + e*x**S(2))*Int((f*x)**m*(a + b*asin(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/((m + S(2))*sqrt(-c**S(2)*x**S(2) + S(1))) + (f*x)**(m + S(1))*(a + b*asin(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(2))))
    rubi.add(rule77)

    pattern78 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Not(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule78 = ReplacementRule(pattern78, lambda b, e, c, m, d, n, a, f, x : b*c*n*sqrt(d + e*x**S(2))*Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1)), x)/(f*(m + S(2))*sqrt(-c**S(2)*x**S(2) + S(1))) + sqrt(d + e*x**S(2))*Int((f*x)**m*(a + b*acos(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/((m + S(2))*sqrt(-c**S(2)*x**S(2) + S(1))) + (f*x)**(m + S(1))*(a + b*acos(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(2))))
    rubi.add(rule78)

    pattern79 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Not(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule79 = ReplacementRule(pattern79, lambda b, e, c, m, d, n, a, f, p, x : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(2)*p + S(1))) + S(2)*d*p*Int((f*x)**m*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(m + S(2)*p + S(1)) + (f*x)**(m + S(1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(2)*p + S(1))))
    rubi.add(rule79)

    pattern80 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda m: Not(RationalQ(m) & Less(m, S(-1)))), CustomConstraint(lambda n, m: RationalQ(m) | ZeroQ(n + S(-1))))
    rule80 = ReplacementRule(pattern80, lambda b, e, c, m, d, n, a, f, p, x : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(f*(m + S(2)*p + S(1))) + S(2)*d*p*Int((f*x)**m*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(m + S(2)*p + S(1)) + (f*x)**(m + S(1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(2)*p + S(1))))
    rubi.add(rule80)

    pattern81 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule81 = ReplacementRule(pattern81, lambda b, e, c, m, d, n, a, f, p, x : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + c**S(2)*(m + S(2)*p + S(3))*Int((f*x)**(m + S(2))*(a + b*asin(c*x))**n*(d + e*x**S(2))**p, x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule81)

    pattern82 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule82 = ReplacementRule(pattern82, lambda b, e, c, m, d, n, a, f, p, x : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(f*(m + S(1))) + c**S(2)*(m + S(2)*p + S(3))*Int((f*x)**(m + S(2))*(a + b*acos(c*x))**n*(d + e*x**S(2))**p, x)/(f**S(2)*(m + S(1))) + (f*x)**(m + S(1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))))
    rubi.add(rule82)

    pattern83 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule83 = ReplacementRule(pattern83, lambda b, e, c, m, d, n, a, f, p, x : b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) - f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*e*(p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule83)

    pattern84 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule84 = ReplacementRule(pattern84, lambda b, e, c, m, d, n, a, f, p, x : -b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*c*(p + S(1))) - f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*e*(p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule84)

    pattern85 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Not(RationalQ(m) & Greater(m, S(1)))), CustomConstraint(lambda n, p, m: IntegerQ(m) | IntegerQ(p) | Equal(n, S(1))))
    rule85 = ReplacementRule(pattern85, lambda b, e, c, m, d, n, a, f, p, x : b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*f*(p + S(1))) + (m + S(2)*p + S(3))*Int((f*x)**m*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))) - (f*x)**(m + S(1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*f*(p + S(1))))
    rubi.add(rule85)

    pattern86 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Not(RationalQ(m) & Greater(m, S(1)))), CustomConstraint(lambda n, p, m: IntegerQ(m) | IntegerQ(p) | Equal(n, S(1))))
    rule86 = ReplacementRule(pattern86, lambda b, e, c, m, d, n, a, f, p, x : -b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(S(2)*f*(p + S(1))) + (m + S(2)*p + S(3))*Int((f*x)**m*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))) - (f*x)**(m + S(1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*f*(p + S(1))))
    rubi.add(rule86)

    pattern87 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule87 = ReplacementRule(pattern87, lambda b, e, c, m, d, n, a, f, x : b*f*n*sqrt(-c**S(2)*x**S(2) + S(1))*Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(-1)), x)/(c*m*sqrt(d + e*x**S(2))) + f*(f*x)**(m + S(-1))*(a + b*asin(c*x))**n*sqrt(d + e*x**S(2))/(e*m) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*asin(c*x))**n/sqrt(d + e*x**S(2)), x)/(c**S(2)*m))
    rubi.add(rule87)

    pattern88 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule88 = ReplacementRule(pattern88, lambda b, e, c, m, d, n, a, f, x : -b*f*n*sqrt(-c**S(2)*x**S(2) + S(1))*Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(-1)), x)/(c*m*sqrt(d + e*x**S(2))) + f*(f*x)**(m + S(-1))*(a + b*acos(c*x))**n*sqrt(d + e*x**S(2))/(e*m) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*acos(c*x))**n/sqrt(d + e*x**S(2)), x)/(c**S(2)*m))
    rubi.add(rule88)

    pattern89 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule89 = ReplacementRule(pattern89, lambda b, e, c, m, d, n, a, x : c**(-m + S(-1))*Subst(Int((a + b*x)**n*sin(x)**m, x), x, asin(c*x))/sqrt(d))
    rubi.add(rule89)

    pattern90 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule90 = ReplacementRule(pattern90, lambda b, e, c, m, d, n, a, x : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*cos(x)**m, x), x, acos(c*x))/sqrt(d))
    rubi.add(rule90)

    pattern91 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule91 = ReplacementRule(pattern91, lambda b, e, c, m, d, a, f, x : -b*c*(f*x)**(m + S(2))*HypergeometricPFQ(List(S(1), m/S(2) + S(1), m/S(2) + S(1)), List(m/S(2) + S(3)/2, m/S(2) + S(2)), c**S(2)*x**S(2))/(sqrt(d)*f**S(2)*(m + S(1))*(m + S(2))) + (f*x)**(m + S(1))*(a + b*asin(c*x))*Hypergeometric2F1(S(1)/2, m/S(2) + S(1)/2, m/S(2) + S(3)/2, c**S(2)*x**S(2))/(sqrt(d)*f*(m + S(1))))
    rubi.add(rule91)

    pattern92 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule92 = ReplacementRule(pattern92, lambda b, e, c, m, d, a, f, x : b*c*(f*x)**(m + S(2))*HypergeometricPFQ(List(S(1), m/S(2) + S(1), m/S(2) + S(1)), List(m/S(2) + S(3)/2, m/S(2) + S(2)), c**S(2)*x**S(2))/(sqrt(d)*f**S(2)*(m + S(1))*(m + S(2))) + (f*x)**(m + S(1))*(a + b*acos(c*x))*Hypergeometric2F1(S(1)/2, m/S(2) + S(1)/2, m/S(2) + S(3)/2, c**S(2)*x**S(2))/(sqrt(d)*f*(m + S(1))))
    rubi.add(rule92)

    pattern93 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda d: Not(PositiveQ(d))), CustomConstraint(lambda n, m: IntegerQ(m) | Equal(n, S(1))))
    rule93 = ReplacementRule(pattern93, lambda b, e, c, m, d, n, a, f, x : sqrt(-c**S(2)*x**S(2) + S(1))*Int((f*x)**m*(a + b*asin(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule93)

    pattern94 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda d: Not(PositiveQ(d))), CustomConstraint(lambda n, m: IntegerQ(m) | Equal(n, S(1))))
    rule94 = ReplacementRule(pattern94, lambda b, e, c, m, d, n, a, f, x : sqrt(-c**S(2)*x**S(2) + S(1))*Int((f*x)**m*(a + b*acos(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule94)

    pattern95 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule95 = ReplacementRule(pattern95, lambda b, e, c, m, d, n, a, f, p, x : b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(c*(m + S(2)*p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(e*(m + S(2)*p + S(1))) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*asin(c*x))**n*(d + e*x**S(2))**p, x)/(c**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule95)

    pattern96 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(1))), CustomConstraint(lambda m: IntegerQ(m)))
    rule96 = ReplacementRule(pattern96, lambda b, e, c, m, d, n, a, f, p, x : -b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x)/(c*(m + S(2)*p + S(1))) + f*(f*x)**(m + S(-1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(e*(m + S(2)*p + S(1))) + f**S(2)*(m + S(-1))*Int((f*x)**(m + S(-2))*(a + b*acos(c*x))**n*(d + e*x**S(2))**p, x)/(c**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule96)

    pattern97 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(1))))
    rule97 = ReplacementRule(pattern97, lambda b, e, c, m, d, n, a, f, p, x : -d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) + (f*x)**m*(a + b*asin(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule97)

    pattern98 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(1))))
    rule98 = ReplacementRule(pattern98, lambda b, e, c, m, d, n, a, f, p, x : d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) - (f*x)**m*(a + b*acos(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule98)

    pattern99 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda d: PositiveQ(d)))
    rule99 = ReplacementRule(pattern99, lambda b, e, c, m, d, n, a, f, x : -f*m*Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(1)), x)/(b*c*sqrt(d)*(n + S(1))) + (f*x)**m*(a + b*asin(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule99)

    pattern100 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda d: PositiveQ(d)))
    rule100 = ReplacementRule(pattern100, lambda b, e, c, m, d, n, a, f, x : f*m*Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(1)), x)/(b*c*sqrt(d)*(n + S(1))) - (f*x)**m*(a + b*acos(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule100)

    pattern101 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(-3))), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)))
    rule101 = ReplacementRule(pattern101, lambda b, e, c, m, d, n, a, f, p, x : c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*(m + S(2)*p + S(1))*Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*f*(n + S(1))) - d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) + (f*x)**m*(a + b*asin(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule101)

    pattern102 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(-3))), CustomConstraint(lambda p: PositiveIntegerQ(S(2)*p)))
    rule102 = ReplacementRule(pattern102, lambda b, e, c, m, d, n, a, f, p, x : -c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*(m + S(2)*p + S(1))*Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*f*(n + S(1))) + d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x)/(b*c*(n + S(1))) - (f*x)**m*(a + b*acos(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))))
    rubi.add(rule102)

    pattern103 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, d: IntegerQ(p) | PositiveQ(d)))
    rule103 = ReplacementRule(pattern103, lambda b, e, c, m, d, n, a, p, x : c**(-m + S(-1))*d**p*Subst(Int((a + b*x)**n*sin(x)**m*cos(x)**(S(2)*p + S(1)), x), x, asin(c*x)))
    rubi.add(rule103)

    pattern104 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, d: IntegerQ(p) | PositiveQ(d)))
    rule104 = ReplacementRule(pattern104, lambda b, e, c, m, d, n, a, p, x : -c**(-m + S(-1))*d**p*Subst(Int((a + b*x)**n*sin(x)**(S(2)*p + S(1))*cos(x)**m, x), x, acos(c*x)))
    rubi.add(rule104)

    pattern105 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, d: Not(IntegerQ(p) | PositiveQ(d))))
    rule105 = ReplacementRule(pattern105, lambda b, e, c, m, d, n, a, p, x : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x**m*(a + b*asin(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x))
    rubi.add(rule105)

    pattern106 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda p: Greater(p, S(-1))), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, d: Not(IntegerQ(p) | PositiveQ(d))))
    rule106 = ReplacementRule(pattern106, lambda b, e, c, m, d, n, a, p, x : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x**m*(a + b*acos(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x))
    rubi.add(rule106)

    pattern107 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda m: Not(PositiveIntegerQ(m/S(2) + S(1)/2))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Less(S(-3), m, S(0))))
    rule107 = ReplacementRule(pattern107, lambda b, e, c, m, d, n, a, f, p, x : Int(ExpandIntegrand((a + b*asin(c*x))**n/sqrt(d + e*x**S(2)), (f*x)**m*(d + e*x**S(2))**(p + S(1)/2), x), x))
    rubi.add(rule107)

    pattern108 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda m: Not(PositiveIntegerQ(m/S(2) + S(1)/2))), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Less(S(-3), m, S(0))))
    rule108 = ReplacementRule(pattern108, lambda b, e, c, m, d, n, a, f, p, x : Int(ExpandIntegrand((a + b*acos(c*x))**n/sqrt(d + e*x**S(2)), (f*x)**m*(d + e*x**S(2))**(p + S(1)/2), x), x))
    rubi.add(rule108)

    pattern109 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule109 = ReplacementRule(pattern109, lambda b, e, c, d, a, p, x : -b*c*Int((d + e*x**S(2))**(p + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(S(2)*e*(p + S(1))) + (a + b*asin(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule109)

    pattern110 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule110 = ReplacementRule(pattern110, lambda b, e, c, d, a, p, x : b*c*Int((d + e*x**S(2))**(p + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(S(2)*e*(p + S(1))) + (a + b*acos(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule110)

    pattern111 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (LessEqual(m + p, S(0)) & PositiveIntegerQ(m/S(2) + S(-1)/2))), )
    def With111(b, e, c, m, d, a, f, p, x):
        u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asin(c*x), u, x)
    rule111 = ReplacementRule(pattern111, lambda b, e, c, m, d, a, f, p, x : With111(b, e, c, m, d, a, f, p, x))
    rubi.add(rule111)

    pattern112 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (LessEqual(m + p, S(0)) & PositiveIntegerQ(m/S(2) + S(-1)/2))), )
    def With112(b, e, c, m, d, a, f, p, x):
        u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
        return b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acos(c*x), u, x)
    rule112 = ReplacementRule(pattern112, lambda b, e, c, m, d, a, f, p, x : With112(b, e, c, m, d, a, f, p, x))
    rubi.add(rule112)

    pattern113 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m)))
    rule113 = ReplacementRule(pattern113, lambda b, e, c, m, d, n, a, f, p, x : Int(ExpandIntegrand((a + b*asin(c*x))**n, (f*x)**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule113)

    pattern114 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda e, c, d: NonzeroQ(c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m)))
    rule114 = ReplacementRule(pattern114, lambda b, e, c, m, d, n, a, f, p, x : Int(ExpandIntegrand((a + b*acos(c*x))**n, (f*x)**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule114)

    pattern115 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule115 = ReplacementRule(pattern115, lambda b, e, c, m, d, n, a, f, p, x : Int((f*x)**m*(a + b*asin(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule115)

    pattern116 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule116 = ReplacementRule(pattern116, lambda b, e, c, m, d, n, a, f, p, x : Int((f*x)**m*(a + b*acos(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule116)

    pattern117 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda f, e, g, d: ZeroQ(d*g + e*f)), CustomConstraint(lambda f, g, c: ZeroQ(c**S(2)*f**S(2) - g**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule117 = ReplacementRule(pattern117, lambda b, e, c, m, d, n, g, h, a, f, p, x : (d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p))*Int((h*x)**m*(a + b*asin(c*x))**n*(d*f + e*g*x**S(2))**p, x))
    rubi.add(rule117)

    pattern118 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda f, e, g, d: ZeroQ(d*g + e*f)), CustomConstraint(lambda f, g, c: ZeroQ(c**S(2)*f**S(2) - g**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule118 = ReplacementRule(pattern118, lambda b, e, c, m, d, n, g, h, a, f, p, x : (d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p))*Int((h*x)**m*(a + b*acos(c*x))**n*(d*f + e*g*x**S(2))**p, x))
    rubi.add(rule118)

    pattern119 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule119 = ReplacementRule(pattern119, lambda b, e, c, d, n, a, x : Subst(Int((a + b*x)**n*cos(x)/(c*d + e*sin(x)), x), x, asin(c*x)))
    rubi.add(rule119)

    pattern120 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule120 = ReplacementRule(pattern120, lambda b, e, c, d, n, a, x : -Subst(Int((a + b*x)**n*sin(x)/(c*d + e*cos(x)), x), x, acos(c*x)))
    rubi.add(rule120)

    pattern121 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule121 = ReplacementRule(pattern121, lambda b, e, c, m, d, n, a, x : -b*c*n*Int((a + b*asin(c*x))**(n + S(-1))*(d + e*x)**(m + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(e*(m + S(1))) + (a + b*asin(c*x))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))))
    rubi.add(rule121)

    pattern122 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule122 = ReplacementRule(pattern122, lambda b, e, c, m, d, n, a, x : b*c*n*Int((a + b*acos(c*x))**(n + S(-1))*(d + e*x)**(m + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x)/(e*(m + S(1))) + (a + b*acos(c*x))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))))
    rubi.add(rule122)

    pattern123 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule123 = ReplacementRule(pattern123, lambda b, e, c, m, d, n, a, x : Int(ExpandIntegrand((a + b*asin(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule123)

    pattern124 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule124 = ReplacementRule(pattern124, lambda b, e, c, m, d, n, a, x : Int(ExpandIntegrand((a + b*acos(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule124)

    pattern125 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule125 = ReplacementRule(pattern125, lambda b, e, c, m, d, n, a, x : c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*d + e*sin(x))**m*cos(x), x), x, asin(c*x)))
    rubi.add(rule125)

    pattern126 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule126 = ReplacementRule(pattern126, lambda b, e, c, m, d, n, a, x : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*d + e*cos(x))**m*sin(x), x), x, acos(c*x)))
    rubi.add(rule126)

    pattern127 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), )
    def With127(b, c, a, Px, x):
        u = IntHide(Px, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asin(c*x), u, x)
    rule127 = ReplacementRule(pattern127, lambda b, c, a, Px, x : With127(b, c, a, Px, x))
    rubi.add(rule127)

    pattern128 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), )
    def With128(b, c, a, Px, x):
        u = IntHide(Px, x)
        return b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acos(c*x), u, x)
    rule128 = ReplacementRule(pattern128, lambda b, c, a, Px, x : With128(b, c, a, Px, x))
    rubi.add(rule128)

    pattern129 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)))
    rule129 = ReplacementRule(pattern129, lambda b, c, n, a, Px, x : Int(ExpandIntegrand(Px*(a + b*asin(c*x))**n, x), x))
    rubi.add(rule129)

    pattern130 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)))
    rule130 = ReplacementRule(pattern130, lambda b, c, n, a, Px, x : Int(ExpandIntegrand(Px*(a + b*acos(c*x))**n, x), x))
    rubi.add(rule130)

    pattern131 = Pattern(Integral(Px_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), )
    def With131(b, e, c, m, d, a, Px, x):
        u = IntHide(Px*(d + e*x)**m, x)
        return -b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asin(c*x), u, x)
    rule131 = ReplacementRule(pattern131, lambda b, e, c, m, d, a, Px, x : With131(b, e, c, m, d, a, Px, x))
    rubi.add(rule131)

    pattern132 = Pattern(Integral(Px_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), )
    def With132(b, e, c, m, d, a, Px, x):
        u = IntHide(Px*(d + e*x)**m, x)
        return b*c*Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acos(c*x), u, x)
    rule132 = ReplacementRule(pattern132, lambda b, e, c, m, d, a, Px, x : With132(b, e, c, m, d, a, Px, x))
    rubi.add(rule132)

    pattern133 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(n, p)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda p, m: Less(m + p + S(1), S(0))), )
    def With133(b, e, c, m, d, n, g, a, f, p, x):
        u = IntHide((d + e*x)**m*(f + g*x)**p, x)
        return -b*c*n*Int(SimplifyIntegrand(u*(a + b*asin(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist((a + b*asin(c*x))**n, u, x)
    rule133 = ReplacementRule(pattern133, lambda b, e, c, m, d, n, g, a, f, p, x : With133(b, e, c, m, d, n, g, a, f, p, x))
    rubi.add(rule133)

    pattern134 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(n, p)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda p, m: Less(m + p + S(1), S(0))), )
    def With134(b, e, c, m, d, n, g, a, f, p, x):
        u = IntHide((d + e*x)**m*(f + g*x)**p, x)
        return b*c*n*Int(SimplifyIntegrand(u*(a + b*acos(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist((a + b*acos(c*x))**n, u, x)
    rule134 = ReplacementRule(pattern134, lambda b, e, c, m, d, n, g, a, f, p, x : With134(b, e, c, m, d, n, g, a, f, p, x))
    rubi.add(rule134)

    pattern135 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_*(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))/(d_ + x_*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(n, p)), CustomConstraint(lambda h, g, e, d: ZeroQ(-S(2)*d*h + e*g)), )
    def With135(b, e, c, d, n, g, h, a, f, p, x):
        u = IntHide((f + g*x + h*x**S(2))**p/(d + e*x)**S(2), x)
        return -b*c*n*Int(SimplifyIntegrand(u*(a + b*asin(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist((a + b*asin(c*x))**n, u, x)
    rule135 = ReplacementRule(pattern135, lambda b, e, c, d, n, g, h, a, f, p, x : With135(b, e, c, d, n, g, h, a, f, p, x))
    rubi.add(rule135)

    pattern136 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_*(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))/(d_ + x_*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(n, p)), CustomConstraint(lambda h, g, e, d: ZeroQ(-S(2)*d*h + e*g)), )
    def With136(b, e, c, d, n, g, h, a, f, p, x):
        u = IntHide((f + g*x + h*x**S(2))**p/(d + e*x)**S(2), x)
        return b*c*n*Int(SimplifyIntegrand(u*(a + b*acos(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist((a + b*acos(c*x))**n, u, x)
    rule136 = ReplacementRule(pattern136, lambda b, e, c, d, n, g, h, a, f, p, x : With136(b, e, c, d, n, g, h, a, f, p, x))
    rubi.add(rule136)

    pattern137 = Pattern(Integral(Px_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule137 = ReplacementRule(pattern137, lambda b, e, c, m, d, n, a, Px, x : Int(ExpandIntegrand(Px*(a + b*asin(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule137)

    pattern138 = Pattern(Integral(Px_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule138 = ReplacementRule(pattern138, lambda b, e, c, m, d, n, a, Px, x : Int(ExpandIntegrand(Px*(a + b*acos(c*x))**n*(d + e*x)**m, x), x))
    rubi.add(rule138)

    pattern139 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p, m: Greater(m, S(3)) | Less(m, -S(2)*p + S(-1))), )
    def With139(b, e, c, m, d, g, a, f, p, x):
        u = IntHide((d + e*x**S(2))**p*(f + g*x)**m, x)
        return -b*c*Int(Dist(1/sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x) + Dist(a + b*asin(c*x), u, x)
    rule139 = ReplacementRule(pattern139, lambda b, e, c, m, d, g, a, f, p, x : With139(b, e, c, m, d, g, a, f, p, x))
    rubi.add(rule139)

    pattern140 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda p, m: Greater(m, S(3)) | Less(m, -S(2)*p + S(-1))), )
    def With140(b, e, c, m, d, g, a, f, p, x):
        u = IntHide((d + e*x**S(2))**p*(f + g*x)**m, x)
        return b*c*Int(Dist(1/sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x) + Dist(a + b*acos(c*x), u, x)
    rule140 = ReplacementRule(pattern140, lambda b, e, c, m, d, g, a, f, p, x : With140(b, e, c, m, d, g, a, f, p, x))
    rubi.add(rule140)

    pattern141 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n, p, m: Equal(m, S(1)) | Greater(p, S(0)) | (Equal(m, S(2)) & Less(p, S(-2))) | (Equal(n, S(1)) & Greater(p, S(-1)))))
    rule141 = ReplacementRule(pattern141, lambda b, e, c, m, d, n, g, a, f, p, x : Int(ExpandIntegrand((a + b*asin(c*x))**n*(d + e*x**S(2))**p, (f + g*x)**m, x), x))
    rubi.add(rule141)

    pattern142 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n, p, m: Equal(m, S(1)) | Greater(p, S(0)) | (Equal(m, S(2)) & Less(p, S(-2))) | (Equal(n, S(1)) & Greater(p, S(-1)))))
    rule142 = ReplacementRule(pattern142, lambda b, e, c, m, d, n, g, a, f, p, x : Int(ExpandIntegrand((a + b*acos(c*x))**n*(d + e*x**S(2))**p, (f + g*x)**m, x), x))
    rubi.add(rule142)

    pattern143 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule143 = ReplacementRule(pattern143, lambda b, e, c, m, d, n, g, a, f, x : (a + b*asin(c*x))**(n + S(1))*(d + e*x**S(2))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))) - Int((a + b*asin(c*x))**(n + S(1))*(f + g*x)**(m + S(-1))*(d*g*m + S(2)*e*f*x + e*g*x**S(2)*(m + S(2))), x)/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule143)

    pattern144 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule144 = ReplacementRule(pattern144, lambda b, e, c, m, d, n, g, a, f, x : -(a + b*acos(c*x))**(n + S(1))*(d + e*x**S(2))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))) + Int((a + b*acos(c*x))**(n + S(1))*(f + g*x)**(m + S(-1))*(d*g*m + S(2)*e*f*x + e*g*x**S(2)*(m + S(2))), x)/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule144)

    pattern145 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule145 = ReplacementRule(pattern145, lambda b, e, c, m, d, n, g, a, f, p, x : Int(ExpandIntegrand((a + b*asin(c*x))**n*sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(-1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule145)

    pattern146 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule146 = ReplacementRule(pattern146, lambda b, e, c, m, d, n, g, a, f, p, x : Int(ExpandIntegrand((a + b*acos(c*x))**n*sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(-1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule146)

    pattern147 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule147 = ReplacementRule(pattern147, lambda b, e, c, m, d, n, g, a, f, p, x : (a + b*asin(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))) - Int(ExpandIntegrand((a + b*asin(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), (d + e*x**S(2))**(p + S(-1)/2)*(d*g*m + e*f*x*(S(2)*p + S(1)) + e*g*x**S(2)*(m + S(2)*p + S(1))), x), x)/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule147)

    pattern148 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: Less(m, S(0))))
    rule148 = ReplacementRule(pattern148, lambda b, e, c, m, d, n, g, a, f, p, x : -(a + b*acos(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))) + Int(ExpandIntegrand((a + b*acos(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), (d + e*x**S(2))**(p + S(-1)/2)*(d*g*m + e*f*x*(S(2)*p + S(1)) + e*g*x**S(2)*(m + S(2)*p + S(1))), x), x)/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule148)

    pattern149 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule149 = ReplacementRule(pattern149, lambda b, e, c, m, d, n, g, a, f, x : -g*m*Int((a + b*asin(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), x)/(b*c*sqrt(d)*(n + S(1))) + (a + b*asin(c*x))**(n + S(1))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule149)

    pattern150 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule150 = ReplacementRule(pattern150, lambda b, e, c, m, d, n, g, a, f, x : g*m*Int((a + b*acos(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), x)/(b*c*sqrt(d)*(n + S(1))) - (a + b*acos(c*x))**(n + S(1))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule150)

    pattern151 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n, m: PositiveIntegerQ(n) | Greater(m, S(0))))
    rule151 = ReplacementRule(pattern151, lambda b, e, c, m, d, n, g, a, f, x : c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*f + g*sin(x))**m, x), x, asin(c*x))/sqrt(d))
    rubi.add(rule151)

    pattern152 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n, m: PositiveIntegerQ(n) | Greater(m, S(0))))
    rule152 = ReplacementRule(pattern152, lambda b, e, c, m, d, n, g, a, f, x : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*(c*f + g*cos(x))**m, x), x, acos(c*x))/sqrt(d))
    rubi.add(rule152)

    pattern153 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule153 = ReplacementRule(pattern153, lambda b, e, c, m, d, n, g, a, f, p, x : Int(ExpandIntegrand((a + b*asin(c*x))**n/sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule153)

    pattern154 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule154 = ReplacementRule(pattern154, lambda b, e, c, m, d, n, g, a, f, p, x : Int(ExpandIntegrand((a + b*acos(c*x))**n/sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m, x), x))
    rubi.add(rule154)

    pattern155 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule155 = ReplacementRule(pattern155, lambda b, e, c, m, d, n, g, a, f, p, x : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*asin(c*x))**n*(f + g*x)**m*(-c**S(2)*x**S(2) + S(1))**p, x))
    rubi.add(rule155)

    pattern156 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule156 = ReplacementRule(pattern156, lambda b, e, c, m, d, n, g, a, f, p, x : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*acos(c*x))**n*(f + g*x)**m*(-c**S(2)*x**S(2) + S(1))**p, x))
    rubi.add(rule156)

    pattern157 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule157 = ReplacementRule(pattern157, lambda b, e, c, m, d, n, g, h, f, a, x : -g*m*Int((a + b*asin(c*x))**(n + S(1))/(f + g*x), x)/(b*c*sqrt(d)*(n + S(1))) + (a + b*asin(c*x))**(n + S(1))*log(h*(f + g*x)**m)/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule157)

    pattern158 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule158 = ReplacementRule(pattern158, lambda b, e, c, m, d, n, g, h, f, a, x : g*m*Int((a + b*acos(c*x))**(n + S(1))/(f + g*x), x)/(b*c*sqrt(d)*(n + S(1))) - (a + b*acos(c*x))**(n + S(1))*log(h*(f + g*x)**m)/(b*c*sqrt(d)*(n + S(1))))
    rubi.add(rule158)

    pattern159 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule159 = ReplacementRule(pattern159, lambda b, e, c, m, d, n, g, h, f, a, p, x : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*asin(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p*log(h*(f + g*x)**m), x))
    rubi.add(rule159)

    pattern160 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule160 = ReplacementRule(pattern160, lambda b, e, c, m, d, n, g, h, f, a, p, x : d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a + b*acos(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p*log(h*(f + g*x)**m), x))
    rubi.add(rule160)

    pattern161 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m + S(1)/2)), )
    def With161(b, e, c, m, d, g, a, f, x):
        u = IntHide((d + e*x)**m*(f + g*x)**m, x)
        return -b*c*Int(Dist(1/sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x) + Dist(a + b*asin(c*x), u, x)
    rule161 = ReplacementRule(pattern161, lambda b, e, c, m, d, g, a, f, x : With161(b, e, c, m, d, g, a, f, x))
    rubi.add(rule161)

    pattern162 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m + S(1)/2)), )
    def With162(b, e, c, m, d, g, a, f, x):
        u = IntHide((d + e*x)**m*(f + g*x)**m, x)
        return b*c*Int(Dist(1/sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x) + Dist(a + b*acos(c*x), u, x)
    rule162 = ReplacementRule(pattern162, lambda b, e, c, m, d, g, a, f, x : With162(b, e, c, m, d, g, a, f, x))
    rubi.add(rule162)

    pattern163 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule163 = ReplacementRule(pattern163, lambda b, e, c, m, d, n, g, a, f, x : Int(ExpandIntegrand((a + b*asin(c*x))**n*(d + e*x)**m*(f + g*x)**m, x), x))
    rubi.add(rule163)

    pattern164 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule164 = ReplacementRule(pattern164, lambda b, e, c, m, d, n, g, a, f, x : Int(ExpandIntegrand((a + b*acos(c*x))**n*(d + e*x)**m*(f + g*x)**m, x), x))
    rubi.add(rule164)

    pattern165 = Pattern(Integral(u_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda v, c, x, a, b: InverseFunctionFreeQ(v, x)))
    def With165(b, c, a, u, x):
        v = IntHide(u, x)
        return -b*c*Int(SimplifyIntegrand(v/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*asin(c*x), v, x)
    rule165 = ReplacementRule(pattern165, lambda b, c, a, u, x : With165(b, c, a, u, x))
    rubi.add(rule165)

    pattern166 = Pattern(Integral(u_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda v, c, x, a, b: InverseFunctionFreeQ(v, x)))
    def With166(b, c, a, u, x):
        v = IntHide(u, x)
        return b*c*Int(SimplifyIntegrand(v/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acos(c*x), v, x)
    rule166 = ReplacementRule(pattern166, lambda b, c, a, u, x : With166(b, c, a, u, x))
    rubi.add(rule166)

    pattern167 = Pattern(Integral(Px_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda x, u: SumQ(u)))
    def With167(b, e, c, d, n, a, Px, p, x):
        u = ExpandIntegrand(Px*(a + b*asin(c*x))**n*(d + e*x**S(2))**p, x)
        return Int(u, x)
    rule167 = ReplacementRule(pattern167, lambda b, e, c, d, n, a, Px, p, x : With167(b, e, c, d, n, a, Px, p, x))
    rubi.add(rule167)

    pattern168 = Pattern(Integral(Px_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda x, u: SumQ(u)))
    def With168(b, e, c, d, n, a, Px, p, x):
        u = ExpandIntegrand(Px*(a + b*acos(c*x))**n*(d + e*x**S(2))**p, x)
        return Int(u, x)
    rule168 = ReplacementRule(pattern168, lambda b, e, c, d, n, a, Px, p, x : With168(b, e, c, d, n, a, Px, p, x))
    rubi.add(rule168)

    pattern169 = Pattern(Integral((f_ + (d_ + x_**S(2)*WC('e', S(1)))**p_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))*WC('Px', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda x, u: SumQ(u)))
    def With169(b, e, c, m, d, n, g, a, Px, f, p, x):
        u = ExpandIntegrand(Px*(a + b*asin(c*x))**n*(f + g*(d + e*x**S(2))**p)**m, x)
        return Int(u, x)
    rule169 = ReplacementRule(pattern169, lambda b, e, c, m, d, n, g, a, Px, f, p, x : With169(b, e, c, m, d, n, g, a, Px, f, p, x))
    rubi.add(rule169)

    pattern170 = Pattern(Integral((f_ + (d_ + x_**S(2)*WC('e', S(1)))**p_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))*WC('Px', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda x, u: SumQ(u)))
    def With170(b, e, c, m, d, n, g, a, Px, f, p, x):
        u = ExpandIntegrand(Px*(a + b*acos(c*x))**n*(f + g*(d + e*x**S(2))**p)**m, x)
        return Int(u, x)
    rule170 = ReplacementRule(pattern170, lambda b, e, c, m, d, n, g, a, Px, f, p, x : With170(b, e, c, m, d, n, g, a, Px, f, p, x))
    rubi.add(rule170)

    pattern171 = Pattern(Integral(RFx_*asin(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda x, RFx: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda x, u: SumQ(u)))
    def With171(n, c, x, RFx):
        u = ExpandIntegrand(asin(c*x)**n, RFx, x)
        return Int(u, x)
    rule171 = ReplacementRule(pattern171, lambda n, c, x, RFx : With171(n, c, x, RFx))
    rubi.add(rule171)

    pattern172 = Pattern(Integral(RFx_*acos(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda x, RFx: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda x, u: SumQ(u)))
    def With172(n, c, x, RFx):
        u = ExpandIntegrand(acos(c*x)**n, RFx, x)
        return Int(u, x)
    rule172 = ReplacementRule(pattern172, lambda n, c, x, RFx : With172(n, c, x, RFx))
    rubi.add(rule172)

    pattern173 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda x, RFx: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule173 = ReplacementRule(pattern173, lambda b, c, n, RFx, a, x : Int(ExpandIntegrand(RFx*(a + b*asin(c*x))**n, x), x))
    rubi.add(rule173)

    pattern174 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda x, RFx: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule174 = ReplacementRule(pattern174, lambda b, c, n, RFx, a, x : Int(ExpandIntegrand(RFx*(a + b*acos(c*x))**n, x), x))
    rubi.add(rule174)

    pattern175 = Pattern(Integral(RFx_*(d_ + x_**S(2)*WC('e', S(1)))**p_*asin(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda x, RFx: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda x, u: SumQ(u)))
    def With175(e, c, d, n, RFx, p, x):
        u = ExpandIntegrand((d + e*x**S(2))**p*asin(c*x)**n, RFx, x)
        return Int(u, x)
    rule175 = ReplacementRule(pattern175, lambda e, c, d, n, RFx, p, x : With175(e, c, d, n, RFx, p, x))
    rubi.add(rule175)

    pattern176 = Pattern(Integral(RFx_*(d_ + x_**S(2)*WC('e', S(1)))**p_*acos(x_*WC('c', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda x, RFx: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)), CustomConstraint(lambda x, u: SumQ(u)))
    def With176(e, c, d, n, RFx, p, x):
        u = ExpandIntegrand((d + e*x**S(2))**p*acos(c*x)**n, RFx, x)
        return Int(u, x)
    rule176 = ReplacementRule(pattern176, lambda e, c, d, n, RFx, p, x : With176(e, c, d, n, RFx, p, x))
    rubi.add(rule176)

    pattern177 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda x, RFx: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule177 = ReplacementRule(pattern177, lambda b, e, c, d, n, RFx, a, p, x : Int(ExpandIntegrand((d + e*x**S(2))**p, RFx*(a + b*asin(c*x))**n, x), x))
    rubi.add(rule177)

    pattern178 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda x, RFx: RationalFunctionQ(RFx, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule178 = ReplacementRule(pattern178, lambda b, e, c, d, n, RFx, a, p, x : Int(ExpandIntegrand((d + e*x**S(2))**p, RFx*(a + b*acos(c*x))**n, x), x))
    rubi.add(rule178)

    pattern179 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule179 = ReplacementRule(pattern179, lambda b, c, n, a, u, x : Int(u*(a + b*asin(c*x))**n, x))
    rubi.add(rule179)

    pattern180 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule180 = ReplacementRule(pattern180, lambda b, c, n, a, u, x : Int(u*(a + b*acos(c*x))**n, x))
    rubi.add(rule180)

    pattern181 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule181 = ReplacementRule(pattern181, lambda b, c, d, n, a, x : Subst(Int((a + b*asin(x))**n, x), x, c + d*x)/d)
    rubi.add(rule181)

    pattern182 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule182 = ReplacementRule(pattern182, lambda b, c, d, n, a, x : Subst(Int((a + b*acos(x))**n, x), x, c + d*x)/d)
    rubi.add(rule182)

    pattern183 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule183 = ReplacementRule(pattern183, lambda b, e, c, m, d, n, a, f, x : Subst(Int((a + b*asin(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule183)

    pattern184 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule184 = ReplacementRule(pattern184, lambda b, e, c, m, d, n, a, f, x : Subst(Int((a + b*acos(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule184)

    pattern185 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda C, B, c, d: ZeroQ(-B*d + S(2)*C*c)))
    rule185 = ReplacementRule(pattern185, lambda b, C, B, c, d, n, a, p, x, A : Subst(Int((a + b*asin(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule185)

    pattern186 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda C, B, c, d: ZeroQ(-B*d + S(2)*C*c)))
    rule186 = ReplacementRule(pattern186, lambda b, C, B, c, d, n, a, p, x, A : Subst(Int((a + b*acos(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule186)

    pattern187 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda C, B, c, d: ZeroQ(-B*d + S(2)*C*c)))
    rule187 = ReplacementRule(pattern187, lambda b, C, e, B, c, m, d, n, a, f, p, x, A : Subst(Int((a + b*asin(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule187)

    pattern188 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(S(2)*A*c*d + B*(-c**S(2) + S(1)))), CustomConstraint(lambda C, B, c, d: ZeroQ(-B*d + S(2)*C*c)))
    rule188 = ReplacementRule(pattern188, lambda b, C, e, B, c, m, d, n, a, f, p, x, A : Subst(Int((a + b*acos(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule188)

    pattern189 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))))
    rule189 = ReplacementRule(pattern189, lambda b, c, d, a, x : sqrt(Pi)*x*(-c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelS(sqrt(c/(Pi*b))*sqrt(a + b*asin(c + d*x**S(2))))/(sqrt(c/b)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))) - sqrt(Pi)*x*(c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelC(sqrt(c/(Pi*b))*sqrt(a + b*asin(c + d*x**S(2))))/(sqrt(c/b)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))) + x*sqrt(a + b*asin(c + d*x**S(2))))
    rubi.add(rule189)

    pattern190 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule190 = ReplacementRule(pattern190, lambda b, a, x, d : -S(2)*sqrt(Pi)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(a/(S(2)*b))*sin(acos(d*x**S(2) + S(1))/S(2))/(d*x*sqrt(1/b)) + S(2)*sqrt(Pi)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(acos(d*x**S(2) + S(1))/S(2))*cos(a/(S(2)*b))/(d*x*sqrt(1/b)) - S(2)*sqrt(a + b*acos(d*x**S(2) + S(1)))*sin(acos(d*x**S(2) + S(1))/S(2))**S(2)/(d*x))
    rubi.add(rule190)

    pattern191 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(-1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule191 = ReplacementRule(pattern191, lambda b, a, x, d : -S(2)*sqrt(Pi)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*cos(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x*sqrt(1/b)) - S(2)*sqrt(Pi)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*sin(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x*sqrt(1/b)) + S(2)*sqrt(a + b*acos(d*x**S(2) + S(-1)))*cos(acos(d*x**S(2) + S(-1))/S(2))**S(2)/(d*x))
    rubi.add(rule191)

    pattern192 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule192 = ReplacementRule(pattern192, lambda b, c, d, n, a, x : -S(4)*b**S(2)*n*(n + S(-1))*Int((a + b*asin(c + d*x**S(2)))**(n + S(-2)), x) + S(2)*b*n*(a + b*asin(c + d*x**S(2)))**(n + S(-1))*sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(d*x) + x*(a + b*asin(c + d*x**S(2)))**n)
    rubi.add(rule192)

    pattern193 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule193 = ReplacementRule(pattern193, lambda b, c, d, n, a, x : -S(4)*b**S(2)*n*(n + S(-1))*Int((a + b*acos(c + d*x**S(2)))**(n + S(-2)), x) - S(2)*b*n*(a + b*acos(c + d*x**S(2)))**(n + S(-1))*sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(d*x) + x*(a + b*acos(c + d*x**S(2)))**n)
    rubi.add(rule193)

    pattern194 = Pattern(Integral(1/(WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))))
    rule194 = ReplacementRule(pattern194, lambda b, c, d, a, x : -x*(c*cos(a/(S(2)*b)) - sin(a/(S(2)*b)))*CosIntegral(c*(a + b*asin(c + d*x**S(2)))/(S(2)*b))/(S(2)*b*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))) - x*(c*cos(a/(S(2)*b)) + sin(a/(S(2)*b)))*SinIntegral(c*(a + b*asin(c + d*x**S(2)))/(S(2)*b))/(S(2)*b*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))))
    rubi.add(rule194)

    pattern195 = Pattern(Integral(1/(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule195 = ReplacementRule(pattern195, lambda b, a, x, d : sqrt(S(2))*x*CosIntegral((a + b*acos(d*x**S(2) + S(1)))/(S(2)*b))*cos(a/(S(2)*b))/(S(2)*b*sqrt(-d*x**S(2))) + sqrt(S(2))*x*SinIntegral((a + b*acos(d*x**S(2) + S(1)))/(S(2)*b))*sin(a/(S(2)*b))/(S(2)*b*sqrt(-d*x**S(2))))
    rubi.add(rule195)

    pattern196 = Pattern(Integral(1/(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(-1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule196 = ReplacementRule(pattern196, lambda b, a, x, d : sqrt(S(2))*x*CosIntegral((a + b*acos(d*x**S(2) + S(-1)))/(S(2)*b))*sin(a/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))) - sqrt(S(2))*x*SinIntegral((a + b*acos(d*x**S(2) + S(-1)))/(S(2)*b))*cos(a/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))))
    rubi.add(rule196)

    pattern197 = Pattern(Integral(1/sqrt(WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))))
    rule197 = ReplacementRule(pattern197, lambda b, c, d, a, x : -sqrt(Pi)*x*(-c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelC(sqrt(a + b*asin(c + d*x**S(2)))/(sqrt(Pi)*sqrt(b*c)))/(sqrt(b*c)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))) - sqrt(Pi)*x*(c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelS(sqrt(a + b*asin(c + d*x**S(2)))/(sqrt(Pi)*sqrt(b*c)))/(sqrt(b*c)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))))
    rubi.add(rule197)

    pattern198 = Pattern(Integral(1/sqrt(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule198 = ReplacementRule(pattern198, lambda b, a, x, d : -S(2)*sqrt(Pi/b)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(acos(d*x**S(2) + S(1))/S(2))*cos(a/(S(2)*b))/(d*x) - S(2)*sqrt(Pi/b)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(a/(S(2)*b))*sin(acos(d*x**S(2) + S(1))/S(2))/(d*x))
    rubi.add(rule198)

    pattern199 = Pattern(Integral(1/sqrt(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(-1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule199 = ReplacementRule(pattern199, lambda b, a, x, d : S(2)*sqrt(Pi/b)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*sin(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x) - S(2)*sqrt(Pi/b)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*cos(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x))
    rubi.add(rule199)

    pattern200 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1))))**(S(-3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))))
    rule200 = ReplacementRule(pattern200, lambda b, c, d, a, x : sqrt(Pi)*x*(c/b)**(S(3)/2)*(-c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelS(sqrt(c/(Pi*b))*sqrt(a + b*asin(c + d*x**S(2))))/(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2))) - sqrt(Pi)*x*(c/b)**(S(3)/2)*(c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelC(sqrt(c/(Pi*b))*sqrt(a + b*asin(c + d*x**S(2))))/(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2))) - sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(b*d*x*sqrt(a + b*asin(c + d*x**S(2)))))
    rubi.add(rule200)

    pattern201 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(1)))**(S(-3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule201 = ReplacementRule(pattern201, lambda b, a, x, d : -S(2)*sqrt(Pi)*(1/b)**(S(3)/2)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(a/(S(2)*b))*sin(acos(d*x**S(2) + S(1))/S(2))/(d*x) + S(2)*sqrt(Pi)*(1/b)**(S(3)/2)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(acos(d*x**S(2) + S(1))/S(2))*cos(a/(S(2)*b))/(d*x) + sqrt(-d**S(2)*x**S(4) - S(2)*d*x**S(2))/(b*d*x*sqrt(a + b*acos(d*x**S(2) + S(1)))))
    rubi.add(rule201)

    pattern202 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(-1)))**(S(-3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule202 = ReplacementRule(pattern202, lambda b, a, x, d : -S(2)*sqrt(Pi)*(1/b)**(S(3)/2)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*cos(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x) - S(2)*sqrt(Pi)*(1/b)**(S(3)/2)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*sin(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x) + sqrt(-d**S(2)*x**S(4) + S(2)*d*x**S(2))/(b*d*x*sqrt(a + b*acos(d*x**S(2) + S(-1)))))
    rubi.add(rule202)

    pattern203 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1))))**(S(-2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))))
    rule203 = ReplacementRule(pattern203, lambda b, c, d, a, x : -sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(S(2)*b*d*x*(a + b*asin(c + d*x**S(2)))) + x*(-c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*SinIntegral(c*(a + b*asin(c + d*x**S(2)))/(S(2)*b))/(S(4)*b**S(2)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))) - x*(c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*CosIntegral(c*(a + b*asin(c + d*x**S(2)))/(S(2)*b))/(S(4)*b**S(2)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))))
    rubi.add(rule203)

    pattern204 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(1)))**(S(-2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule204 = ReplacementRule(pattern204, lambda b, a, x, d : sqrt(-d**S(2)*x**S(4) - S(2)*d*x**S(2))/(S(2)*b*d*x*(a + b*acos(d*x**S(2) + S(1)))) + sqrt(S(2))*x*CosIntegral((a + b*acos(d*x**S(2) + S(1)))/(S(2)*b))*sin(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(-d*x**S(2))) - sqrt(S(2))*x*SinIntegral((a + b*acos(d*x**S(2) + S(1)))/(S(2)*b))*cos(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(-d*x**S(2))))
    rubi.add(rule204)

    pattern205 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(-1)))**(S(-2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)))
    rule205 = ReplacementRule(pattern205, lambda b, a, x, d : sqrt(-d**S(2)*x**S(4) + S(2)*d*x**S(2))/(S(2)*b*d*x*(a + b*acos(d*x**S(2) + S(-1)))) - sqrt(S(2))*x*CosIntegral((a + b*acos(d*x**S(2) + S(-1)))/(S(2)*b))*cos(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))) - sqrt(S(2))*x*SinIntegral((a + b*acos(d*x**S(2) + S(-1)))/(S(2)*b))*sin(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))))
    rubi.add(rule205)

    pattern206 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule206 = ReplacementRule(pattern206, lambda b, c, d, n, a, x : (a + b*asin(c + d*x**S(2)))**(n + S(1))*sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(S(2)*b*d*x*(n + S(1))) + x*(a + b*asin(c + d*x**S(2)))**(n + S(2))/(S(4)*b**S(2)*(n + S(1))*(n + S(2))) - Int((a + b*asin(c + d*x**S(2)))**(n + S(2)), x)/(S(4)*b**S(2)*(n + S(1))*(n + S(2))))
    rubi.add(rule206)

    pattern207 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c: ZeroQ(c**S(2) + S(-1))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule207 = ReplacementRule(pattern207, lambda b, c, d, n, a, x : -(a + b*acos(c + d*x**S(2)))**(n + S(1))*sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(S(2)*b*d*x*(n + S(1))) + x*(a + b*acos(c + d*x**S(2)))**(n + S(2))/(S(4)*b**S(2)*(n + S(1))*(n + S(2))) - Int((a + b*acos(c + d*x**S(2)))**(n + S(2)), x)/(S(4)*b**S(2)*(n + S(1))*(n + S(2))))
    rubi.add(rule207)

    pattern208 = Pattern(Integral(asin(x_**p_*WC('a', S(1)))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule208 = ReplacementRule(pattern208, lambda n, a, p, x : Subst(Int(x**n*cot(x), x), x, asin(a*x**p))/p)
    rubi.add(rule208)

    pattern209 = Pattern(Integral(acos(x_**p_*WC('a', S(1)))**WC('n', S(1))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule209 = ReplacementRule(pattern209, lambda n, a, p, x : -Subst(Int(x**n*tan(x), x), x, acos(a*x**p))/p)
    rubi.add(rule209)

    pattern210 = Pattern(Integral(WC('u', S(1))*asin(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule210 = ReplacementRule(pattern210, lambda b, c, m, n, a, u, x : Int(u*acsc(a/c + b*x**n/c)**m, x))
    rubi.add(rule210)

    pattern211 = Pattern(Integral(WC('u', S(1))*acos(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule211 = ReplacementRule(pattern211, lambda b, c, m, n, a, u, x : Int(u*asec(a/c + b*x**n/c)**m, x))
    rubi.add(rule211)

    pattern212 = Pattern(Integral(asin(sqrt(x_**S(2)*WC('b', S(1)) + S(1)))**WC('n', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule212 = ReplacementRule(pattern212, lambda b, x, n : sqrt(-b*x**S(2))*Subst(Int(asin(x)**n/sqrt(-x**S(2) + S(1)), x), x, sqrt(b*x**S(2) + S(1)))/(b*x))
    rubi.add(rule212)

    pattern213 = Pattern(Integral(acos(sqrt(x_**S(2)*WC('b', S(1)) + S(1)))**WC('n', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule213 = ReplacementRule(pattern213, lambda b, x, n : sqrt(-b*x**S(2))*Subst(Int(acos(x)**n/sqrt(-x**S(2) + S(1)), x), x, sqrt(b*x**S(2) + S(1)))/(b*x))
    rubi.add(rule213)

    pattern214 = Pattern(Integral(f_**(WC('c', S(1))*asin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule214 = ReplacementRule(pattern214, lambda b, c, n, a, f, u, x : Subst(Int(f**(c*x**n)*ReplaceAll(u, Rule(x, -a/b + sin(x)/b))*cos(x), x), x, asin(a + b*x))/b)
    rubi.add(rule214)

    pattern215 = Pattern(Integral(f_**(WC('c', S(1))*acos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule215 = ReplacementRule(pattern215, lambda b, c, n, a, f, u, x : -Subst(Int(f**(c*x**n)*ReplaceAll(u, Rule(x, -a/b + cos(x)/b))*sin(x), x), x, acos(a + b*x))/b)
    rubi.add(rule215)

    pattern216 = Pattern(Integral(asin(x_**S(2)*WC('a', S(1)) + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c: EqQ(b**S(2)*c, S(1))))
    rule216 = ReplacementRule(pattern216, lambda b, c, d, a, x : x*asin(a*x**S(2) + b*sqrt(c + d*x**S(2))) - x*sqrt(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)*Int(x*(S(2)*a*sqrt(c + d*x**S(2)) + b*d)/(sqrt(c + d*x**S(2))*sqrt(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)), x)/sqrt(-x**S(2)*(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)))
    rubi.add(rule216)

    pattern217 = Pattern(Integral(acos(x_**S(2)*WC('a', S(1)) + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, c: EqQ(b**S(2)*c, S(1))))
    rule217 = ReplacementRule(pattern217, lambda b, c, d, a, x : x*acos(a*x**S(2) + b*sqrt(c + d*x**S(2))) + x*sqrt(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)*Int(x*(S(2)*a*sqrt(c + d*x**S(2)) + b*d)/(sqrt(c + d*x**S(2))*sqrt(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)), x)/sqrt(-x**S(2)*(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)))
    rubi.add(rule217)

    pattern218 = Pattern(Integral(asin(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda u, x: Not(FunctionOfExponentialQ(u, x))))
    rule218 = ReplacementRule(pattern218, lambda u, x : x*asin(u) - Int(SimplifyIntegrand(x*D(u, x)/sqrt(-u**S(2) + S(1)), x), x))
    rubi.add(rule218)

    pattern219 = Pattern(Integral(acos(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda u, x: Not(FunctionOfExponentialQ(u, x))))
    rule219 = ReplacementRule(pattern219, lambda u, x : x*acos(u) + Int(SimplifyIntegrand(x*D(u, x)/sqrt(-u**S(2) + S(1)), x), x))
    rubi.add(rule219)

    pattern220 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(u_)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, m, d, u, x: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, x: Not(FunctionOfExponentialQ(u, x))))
    rule220 = ReplacementRule(pattern220, lambda b, c, m, d, a, u, x : -b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/sqrt(-u**S(2) + S(1)), x), x)/(d*(m + S(1))) + (a + b*asin(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule220)

    pattern221 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(u_)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, m, d, u, x: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, x: Not(FunctionOfExponentialQ(u, x))))
    rule221 = ReplacementRule(pattern221, lambda b, c, m, d, a, u, x : b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/sqrt(-u**S(2) + S(1)), x), x)/(d*(m + S(1))) + (a + b*acos(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule221)

    pattern222 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*asin(u_)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda w, u, x, a, b: InverseFunctionFreeQ(w, x)))
    def With222(b, a, u, v, x):
        w = IntHide(v, x)
        return -b*Int(SimplifyIntegrand(w*D(u, x)/sqrt(-u**S(2) + S(1)), x), x) + Dist(a + b*asin(u), w, x)
    rule222 = ReplacementRule(pattern222, lambda b, a, u, v, x : With222(b, a, u, v, x))
    rubi.add(rule222)

    pattern223 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acos(u_)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda w, u, x, a, b: InverseFunctionFreeQ(w, x)))
    def With223(b, a, u, v, x):
        w = IntHide(v, x)
        return b*Int(SimplifyIntegrand(w*D(u, x)/sqrt(-u**S(2) + S(1)), x), x) + Dist(a + b*acos(u), w, x)
    rule223 = ReplacementRule(pattern223, lambda b, a, u, v, x : With223(b, a, u, v, x))
    rubi.add(rule223)

    pattern224 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule224 = ReplacementRule(pattern224, lambda b, c, n, a, x : -b*c*n*Int(x*(a + b*ArcTan(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x) + x*(a + b*ArcTan(c*x))**n)
    rubi.add(rule224)

    pattern225 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule225 = ReplacementRule(pattern225, lambda b, c, n, a, x : b*c*n*Int(x*(a + b*acot(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x) + x*(a + b*acot(c*x))**n)
    rubi.add(rule225)

    pattern226 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule226 = ReplacementRule(pattern226, lambda b, c, n, a, x : Int((a + b*ArcTan(c*x))**n, x))
    rubi.add(rule226)

    pattern227 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule227 = ReplacementRule(pattern227, lambda b, c, n, a, x : Int((a + b*acot(c*x))**n, x))
    rubi.add(rule227)

    pattern228 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule228 = ReplacementRule(pattern228, lambda b, e, c, d, n, a, x : b*c*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*log(S(2)*d/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x)/e - (a + b*ArcTan(c*x))**n*log(S(2)*d/(d + e*x))/e)
    rubi.add(rule228)

    pattern229 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule229 = ReplacementRule(pattern229, lambda b, e, c, d, n, a, x : -b*c*n*Int((a + b*acot(c*x))**(n + S(-1))*log(S(2)*d/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x)/e - (a + b*acot(c*x))**n*log(S(2)*d/(d + e*x))/e)
    rubi.add(rule229)

    pattern230 = Pattern(Integral(ArcTan(x_*WC('c', S(1)))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: PositiveQ(ImaginaryI*c*d/e + S(1))), CustomConstraint(lambda e, c, d: NegativeQ(ImaginaryI*c*d/e + S(-1))))
    rule230 = ReplacementRule(pattern230, lambda e, c, x, d : ImaginaryI*PolyLog(S(2), Simp(ImaginaryI*c*(d + e*x)/(ImaginaryI*c*d - e), x))/(S(2)*e) - ImaginaryI*PolyLog(S(2), Simp(ImaginaryI*c*(d + e*x)/(ImaginaryI*c*d + e), x))/(S(2)*e) - ArcTan(c*d/e)*log(d + e*x)/e)
    rubi.add(rule230)

    pattern231 = Pattern(Integral(ArcTan(x_*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule231 = ReplacementRule(pattern231, lambda e, c, x, d : ImaginaryI*Int(log(-ImaginaryI*c*x + S(1))/(d + e*x), x)/S(2) - ImaginaryI*Int(log(ImaginaryI*c*x + S(1))/(d + e*x), x)/S(2))
    rubi.add(rule231)

    pattern232 = Pattern(Integral(acot(x_*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule232 = ReplacementRule(pattern232, lambda e, c, x, d : ImaginaryI*Int(log(-ImaginaryI/(c*x) + S(1))/(d + e*x), x)/S(2) - ImaginaryI*Int(log(ImaginaryI/(c*x) + S(1))/(d + e*x), x)/S(2))
    rubi.add(rule232)

    pattern233 = Pattern(Integral((a_ + ArcTan(x_*WC('c', S(1)))*WC('b', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule233 = ReplacementRule(pattern233, lambda b, e, c, d, a, x : a*log(RemoveContent(d + e*x, x))/e + b*Int(ArcTan(c*x)/(d + e*x), x))
    rubi.add(rule233)

    pattern234 = Pattern(Integral((a_ + WC('b', S(1))*acot(x_*WC('c', S(1))))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule234 = ReplacementRule(pattern234, lambda b, e, c, d, a, x : a*log(RemoveContent(d + e*x, x))/e + b*Int(acot(c*x)/(d + e*x), x))
    rubi.add(rule234)

    pattern235 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule235 = ReplacementRule(pattern235, lambda b, e, c, d, a, p, x : -b*c*Int((d + e*x)**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x)/(e*(p + S(1))) + (a + b*ArcTan(c*x))*(d + e*x)**(p + S(1))/(e*(p + S(1))))
    rubi.add(rule235)

    pattern236 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule236 = ReplacementRule(pattern236, lambda b, e, c, d, a, p, x : b*c*Int((d + e*x)**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x)/(e*(p + S(1))) + (a + b*acot(c*x))*(d + e*x)**(p + S(1))/(e*(p + S(1))))
    rubi.add(rule236)

    pattern237 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule237 = ReplacementRule(pattern237, lambda b, c, n, a, x : -S(2)*b*c*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*atanh(-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))/(c**S(2)*x**S(2) + S(1)), x) + S(2)*(a + b*ArcTan(c*x))**n*atanh(-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1)))
    rubi.add(rule237)

    pattern238 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule238 = ReplacementRule(pattern238, lambda b, c, n, a, x : S(2)*b*c*n*Int((a + b*acot(c*x))**(n + S(-1))*acoth(-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))/(c**S(2)*x**S(2) + S(1)), x) + S(2)*(a + b*acot(c*x))**n*acoth(-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1)))
    rubi.add(rule238)

    pattern239 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule239 = ReplacementRule(pattern239, lambda b, c, m, n, a, x : -b*c*n*Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcTan(c*x))**n/(m + S(1)))
    rubi.add(rule239)

    pattern240 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule240 = ReplacementRule(pattern240, lambda b, c, m, n, a, x : b*c*n*Int(x**(m + S(1))*(a + b*acot(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*acot(c*x))**n/(m + S(1)))
    rubi.add(rule240)

    pattern241 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(n, p)))
    rule241 = ReplacementRule(pattern241, lambda b, e, c, d, n, a, p, x : Int(ExpandIntegrand((a + b*ArcTan(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule241)

    pattern242 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(n, p)))
    rule242 = ReplacementRule(pattern242, lambda b, e, c, d, n, a, p, x : Int(ExpandIntegrand((a + b*acot(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule242)

    pattern243 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule243 = ReplacementRule(pattern243, lambda b, e, c, d, n, a, p, x : Int((a + b*ArcTan(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule243)

    pattern244 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule244 = ReplacementRule(pattern244, lambda b, e, c, d, n, a, p, x : Int((a + b*acot(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule244)

    pattern245 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule245 = ReplacementRule(pattern245, lambda b, e, c, m, d, n, a, x : -d*Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**n/(d + e*x), x)/e + Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**n, x)/e)
    rubi.add(rule245)

    pattern246 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule246 = ReplacementRule(pattern246, lambda b, e, c, m, d, n, a, x : -d*Int(x**(m + S(-1))*(a + b*acot(c*x))**n/(d + e*x), x)/e + Int(x**(m + S(-1))*(a + b*acot(c*x))**n, x)/e)
    rubi.add(rule246)

    pattern247 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*(d_ + x_*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule247 = ReplacementRule(pattern247, lambda b, e, c, d, n, a, x : -b*c*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*log(S(2)*e*x/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x)/d + (a + b*ArcTan(c*x))**n*log(S(2)*e*x/(d + e*x))/d)
    rubi.add(rule247)

    pattern248 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule248 = ReplacementRule(pattern248, lambda b, e, c, d, n, a, x : b*c*n*Int((a + b*acot(c*x))**(n + S(-1))*log(S(2)*e*x/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x)/d + (a + b*acot(c*x))**n*log(S(2)*e*x/(d + e*x))/d)
    rubi.add(rule248)

    pattern249 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule249 = ReplacementRule(pattern249, lambda b, e, c, m, d, n, a, x : -e*Int(x**(m + S(1))*(a + b*ArcTan(c*x))**n/(d + e*x), x)/d + Int(x**m*(a + b*ArcTan(c*x))**n, x)/d)
    rubi.add(rule249)

    pattern250 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d**S(2) + e**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule250 = ReplacementRule(pattern250, lambda b, e, c, m, d, n, a, x : -e*Int(x**(m + S(1))*(a + b*acot(c*x))**n/(d + e*x), x)/d + Int(x**m*(a + b*acot(c*x))**n, x)/d)
    rubi.add(rule250)

    pattern251 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda a, p, m: IntegerQ(m) | NonzeroQ(a) | Greater(p, S(0))))
    rule251 = ReplacementRule(pattern251, lambda b, e, c, m, d, n, a, p, x : Int(ExpandIntegrand(x**m*(a + b*ArcTan(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule251)

    pattern252 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda a, p, m: IntegerQ(m) | NonzeroQ(a) | Greater(p, S(0))))
    rule252 = ReplacementRule(pattern252, lambda b, e, c, m, d, n, a, p, x : Int(ExpandIntegrand(x**m*(a + b*acot(c*x))**n*(d + e*x)**p, x), x))
    rubi.add(rule252)

    pattern253 = Pattern(Integral(x_**WC('m', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule253 = ReplacementRule(pattern253, lambda b, e, c, m, d, n, a, p, x : Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule253)

    pattern254 = Pattern(Integral(x_**WC('m', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule254 = ReplacementRule(pattern254, lambda b, e, c, m, d, n, a, p, x : Int(x**m*(a + b*acot(c*x))**n*(d + e*x)**p, x))
    rubi.add(rule254)

    pattern255 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule255 = ReplacementRule(pattern255, lambda b, e, c, d, a, p, x : -b*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*ArcTan(c*x))*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule255)

    pattern256 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule256 = ReplacementRule(pattern256, lambda b, e, c, d, a, p, x : b*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*acot(c*x))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*acot(c*x))*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule256)

    pattern257 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule257 = ReplacementRule(pattern257, lambda b, e, c, d, n, a, p, x : b**S(2)*d*n*(n + S(-1))*Int((a + b*ArcTan(c*x))**(n + S(-2))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p*(S(2)*p + S(1))) - b*n*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule257)

    pattern258 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule258 = ReplacementRule(pattern258, lambda b, e, c, d, n, a, p, x : b**S(2)*d*n*(n + S(-1))*Int((a + b*acot(c*x))**(n + S(-2))*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p*(S(2)*p + S(1))) + b*n*(a + b*acot(c*x))**(n + S(-1))*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))) + S(2)*d*p*Int((a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x)/(S(2)*p + S(1)) + x*(a + b*acot(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)))
    rubi.add(rule258)

    pattern259 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)))
    rule259 = ReplacementRule(pattern259, lambda b, e, c, d, a, x : log(RemoveContent(a + b*ArcTan(c*x), x))/(b*c*d))
    rubi.add(rule259)

    pattern260 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)))
    rule260 = ReplacementRule(pattern260, lambda b, e, c, d, a, x : -log(RemoveContent(a + b*acot(c*x), x))/(b*c*d))
    rubi.add(rule260)

    pattern261 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule261 = ReplacementRule(pattern261, lambda b, e, c, d, n, a, x : (a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule261)

    pattern262 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule262 = ReplacementRule(pattern262, lambda b, e, c, d, n, a, x : -(a + b*acot(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule262)

    pattern263 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule263 = ReplacementRule(pattern263, lambda b, e, c, d, a, x : ImaginaryI*b*PolyLog(S(2), -ImaginaryI*sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/(c*sqrt(d)) - ImaginaryI*b*PolyLog(S(2), ImaginaryI*sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/(c*sqrt(d)) - S(2)*ImaginaryI*(a + b*ArcTan(c*x))*ArcTan(sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/(c*sqrt(d)))
    rubi.add(rule263)

    pattern264 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule264 = ReplacementRule(pattern264, lambda b, e, c, d, a, x : -ImaginaryI*b*PolyLog(S(2), -ImaginaryI*sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/(c*sqrt(d)) + ImaginaryI*b*PolyLog(S(2), ImaginaryI*sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/(c*sqrt(d)) - S(2)*ImaginaryI*(a + b*acot(c*x))*ArcTan(sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/(c*sqrt(d)))
    rubi.add(rule264)

    pattern265 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule265 = ReplacementRule(pattern265, lambda b, e, c, d, n, a, x : Subst(Int((a + b*x)**n*sec(x), x), x, ArcTan(c*x))/(c*sqrt(d)))
    rubi.add(rule265)

    pattern266 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule266 = ReplacementRule(pattern266, lambda b, e, c, d, n, a, x : -x*sqrt(S(1) + S(1)/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*csc(x), x), x, acot(c*x))/sqrt(d + e*x**S(2)))
    rubi.add(rule266)

    pattern267 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule267 = ReplacementRule(pattern267, lambda b, e, c, d, n, a, x : sqrt(c**S(2)*x**S(2) + S(1))*Int((a + b*ArcTan(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule267)

    pattern268 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule268 = ReplacementRule(pattern268, lambda b, e, c, d, n, a, x : sqrt(c**S(2)*x**S(2) + S(1))*Int((a + b*acot(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule268)

    pattern269 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule269 = ReplacementRule(pattern269, lambda b, e, c, d, n, a, x : -b*c*n*Int(x*(a + b*ArcTan(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/S(2) + x*(a + b*ArcTan(c*x))**n/(S(2)*d*(d + e*x**S(2))) + (a + b*ArcTan(c*x))**(n + S(1))/(S(2)*b*c*d**S(2)*(n + S(1))))
    rubi.add(rule269)

    pattern270 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule270 = ReplacementRule(pattern270, lambda b, e, c, d, n, a, x : b*c*n*Int(x*(a + b*acot(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/S(2) + x*(a + b*acot(c*x))**n/(S(2)*d*(d + e*x**S(2))) - (a + b*acot(c*x))**(n + S(1))/(S(2)*b*c*d**S(2)*(n + S(1))))
    rubi.add(rule270)

    pattern271 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)))
    rule271 = ReplacementRule(pattern271, lambda b, e, c, d, a, x : b/(c*d*sqrt(d + e*x**S(2))) + x*(a + b*ArcTan(c*x))/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule271)

    pattern272 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)))
    rule272 = ReplacementRule(pattern272, lambda b, e, c, d, a, x : -b/(c*d*sqrt(d + e*x**S(2))) + x*(a + b*acot(c*x))/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule272)

    pattern273 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule273 = ReplacementRule(pattern273, lambda b, e, c, d, a, p, x : b*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule273)

    pattern274 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule274 = ReplacementRule(pattern274, lambda b, e, c, d, a, p, x : -b*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule274)

    pattern275 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule275 = ReplacementRule(pattern275, lambda b, e, c, d, n, a, x : -b**S(2)*n*(n + S(-1))*Int((a + b*ArcTan(c*x))**(n + S(-2))/(d + e*x**S(2))**(S(3)/2), x) + b*n*(a + b*ArcTan(c*x))**(n + S(-1))/(c*d*sqrt(d + e*x**S(2))) + x*(a + b*ArcTan(c*x))**n/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule275)

    pattern276 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule276 = ReplacementRule(pattern276, lambda b, e, c, d, n, a, x : -b**S(2)*n*(n + S(-1))*Int((a + b*acot(c*x))**(n + S(-2))/(d + e*x**S(2))**(S(3)/2), x) - b*n*(a + b*acot(c*x))**(n + S(-1))/(c*d*sqrt(d + e*x**S(2))) + x*(a + b*acot(c*x))**n/(d*sqrt(d + e*x**S(2))))
    rubi.add(rule276)

    pattern277 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule277 = ReplacementRule(pattern277, lambda b, e, c, d, n, a, p, x : -b**S(2)*n*(n + S(-1))*Int((a + b*ArcTan(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/(S(4)*(p + S(1))**S(2)) + b*n*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule277)

    pattern278 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)))
    rule278 = ReplacementRule(pattern278, lambda b, e, c, d, n, a, p, x : -b**S(2)*n*(n + S(-1))*Int((a + b*acot(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/(S(4)*(p + S(1))**S(2)) - b*n*(a + b*acot(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)) - x*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))) + (S(2)*p + S(3))*Int((a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*d*(p + S(1))))
    rubi.add(rule278)

    pattern279 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))))
    rule279 = ReplacementRule(pattern279, lambda b, e, c, d, n, a, p, x : -S(2)*c*(p + S(1))*Int(x*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) + (a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule279)

    pattern280 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))))
    rule280 = ReplacementRule(pattern280, lambda b, e, c, d, n, a, p, x : S(2)*c*(p + S(1))*Int(x*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) - (a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule280)

    pattern281 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda p, d: IntegerQ(p) | PositiveQ(d)))
    rule281 = ReplacementRule(pattern281, lambda b, e, c, d, n, a, p, x : d**p*Subst(Int((a + b*x)**n*cos(x)**(-S(2)*p + S(-2)), x), x, ArcTan(c*x))/c)
    rubi.add(rule281)

    pattern282 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda p, d: Not(IntegerQ(p) | PositiveQ(d))))
    rule282 = ReplacementRule(pattern282, lambda b, e, c, d, n, a, p, x : d**(p + S(1)/2)*sqrt(c**S(2)*x**S(2) + S(1))*Int((a + b*ArcTan(c*x))**n*(c**S(2)*x**S(2) + S(1))**p, x)/sqrt(d + e*x**S(2)))
    rubi.add(rule282)

    pattern283 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule283 = ReplacementRule(pattern283, lambda b, e, c, d, n, a, p, x : -d**p*Subst(Int((a + b*x)**n*sin(x)**(-S(2)*p + S(-2)), x), x, acot(c*x))/c)
    rubi.add(rule283)

    pattern284 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: NegativeIntegerQ(S(2)*p + S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule284 = ReplacementRule(pattern284, lambda b, e, c, d, n, a, p, x : -d**(p + S(1)/2)*x*sqrt((c**S(2)*x**S(2) + S(1))/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*sin(x)**(-S(2)*p + S(-2)), x), x, acot(c*x))/sqrt(d + e*x**S(2)))
    rubi.add(rule284)

    pattern285 = Pattern(Integral(ArcTan(x_*WC('c', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule285 = ReplacementRule(pattern285, lambda e, c, x, d : ImaginaryI*Int(log(-ImaginaryI*c*x + S(1))/(d + e*x**S(2)), x)/S(2) - ImaginaryI*Int(log(ImaginaryI*c*x + S(1))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule285)

    pattern286 = Pattern(Integral(acot(x_*WC('c', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule286 = ReplacementRule(pattern286, lambda e, c, x, d : ImaginaryI*Int(log(-ImaginaryI/(c*x) + S(1))/(d + e*x**S(2)), x)/S(2) - ImaginaryI*Int(log(ImaginaryI/(c*x) + S(1))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule286)

    pattern287 = Pattern(Integral((a_ + ArcTan(x_*WC('c', S(1)))*WC('b', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule287 = ReplacementRule(pattern287, lambda b, e, c, d, a, x : a*Int(1/(d + e*x**S(2)), x) + b*Int(ArcTan(c*x)/(d + e*x**S(2)), x))
    rubi.add(rule287)

    pattern288 = Pattern(Integral((a_ + WC('b', S(1))*acot(x_*WC('c', S(1))))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule288 = ReplacementRule(pattern288, lambda b, e, c, d, a, x : a*Int(1/(d + e*x**S(2)), x) + b*Int(acot(c*x)/(d + e*x**S(2)), x))
    rubi.add(rule288)

    pattern289 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With289(b, e, c, d, a, p, x):
        u = IntHide((d + e*x**S(2))**p, x)
        return -b*c*Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcTan(c*x), u, x)
    rule289 = ReplacementRule(pattern289, lambda b, e, c, d, a, p, x : With289(b, e, c, d, a, p, x))
    rubi.add(rule289)

    pattern290 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With290(b, e, c, d, a, p, x):
        u = IntHide((d + e*x**S(2))**p, x)
        return b*c*Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acot(c*x), u, x)
    rule290 = ReplacementRule(pattern290, lambda b, e, c, d, a, p, x : With290(b, e, c, d, a, p, x))
    rubi.add(rule290)

    pattern291 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule291 = ReplacementRule(pattern291, lambda b, e, c, d, n, a, p, x : Int(ExpandIntegrand((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule291)

    pattern292 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule292 = ReplacementRule(pattern292, lambda b, e, c, d, n, a, p, x : Int(ExpandIntegrand((a + b*acot(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule292)

    pattern293 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule293 = ReplacementRule(pattern293, lambda b, e, c, d, n, a, p, x : Int((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule293)

    pattern294 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule294 = ReplacementRule(pattern294, lambda b, e, c, d, n, a, p, x : Int((a + b*acot(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule294)

    pattern295 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule295 = ReplacementRule(pattern295, lambda b, e, c, m, d, n, a, x : -d*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n, x)/e)
    rubi.add(rule295)

    pattern296 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule296 = ReplacementRule(pattern296, lambda b, e, c, m, d, n, a, x : -d*Int(x**(m + S(-2))*(a + b*acot(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*acot(c*x))**n, x)/e)
    rubi.add(rule296)

    pattern297 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule297 = ReplacementRule(pattern297, lambda b, e, c, m, d, n, a, x : -e*Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*ArcTan(c*x))**n, x)/d)
    rubi.add(rule297)

    pattern298 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule298 = ReplacementRule(pattern298, lambda b, e, c, m, d, n, a, x : -e*Int(x**(m + S(2))*(a + b*acot(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*acot(c*x))**n, x)/d)
    rubi.add(rule298)

    pattern299 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule299 = ReplacementRule(pattern299, lambda b, e, c, d, n, a, x : -ImaginaryI*(a + b*ArcTan(c*x))**(n + S(1))/(b*e*(n + S(1))) - Int((a + b*ArcTan(c*x))**n/(ImaginaryI - c*x), x)/(c*d))
    rubi.add(rule299)

    pattern300 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule300 = ReplacementRule(pattern300, lambda b, e, c, d, n, a, x : ImaginaryI*(a + b*acot(c*x))**(n + S(1))/(b*e*(n + S(1))) - Int((a + b*acot(c*x))**n/(ImaginaryI - c*x), x)/(c*d))
    rubi.add(rule300)

    pattern301 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule301 = ReplacementRule(pattern301, lambda b, e, c, d, n, a, x : x*(a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(n + S(1))) - Int((a + b*ArcTan(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))))
    rubi.add(rule301)

    pattern302 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule302 = ReplacementRule(pattern302, lambda b, e, c, d, n, a, x : -x*(a + b*acot(c*x))**(n + S(1))/(b*c*d*(n + S(1))) + Int((a + b*acot(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))))
    rubi.add(rule302)

    pattern303 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule303 = ReplacementRule(pattern303, lambda b, e, c, m, d, n, a, x : -d*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n, x)/e)
    rubi.add(rule303)

    pattern304 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule304 = ReplacementRule(pattern304, lambda b, e, c, m, d, n, a, x : -d*Int(x**(m + S(-2))*(a + b*acot(c*x))**n/(d + e*x**S(2)), x)/e + Int(x**(m + S(-2))*(a + b*acot(c*x))**n, x)/e)
    rubi.add(rule304)

    pattern305 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule305 = ReplacementRule(pattern305, lambda b, e, c, d, n, a, x : ImaginaryI*Int((a + b*ArcTan(c*x))**n/(x*(ImaginaryI + c*x)), x)/d - ImaginaryI*(a + b*ArcTan(c*x))**(n + S(1))/(b*d*(n + S(1))))
    rubi.add(rule305)

    pattern306 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule306 = ReplacementRule(pattern306, lambda b, e, c, d, n, a, x : ImaginaryI*Int((a + b*acot(c*x))**n/(x*(ImaginaryI + c*x)), x)/d + ImaginaryI*(a + b*acot(c*x))**(n + S(1))/(b*d*(n + S(1))))
    rubi.add(rule306)

    pattern307 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule307 = ReplacementRule(pattern307, lambda b, e, c, m, d, n, a, x : -e*Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*ArcTan(c*x))**n, x)/d)
    rubi.add(rule307)

    pattern308 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))))
    rule308 = ReplacementRule(pattern308, lambda b, e, c, m, d, n, a, x : -e*Int(x**(m + S(2))*(a + b*acot(c*x))**n/(d + e*x**S(2)), x)/d + Int(x**m*(a + b*acot(c*x))**n, x)/d)
    rubi.add(rule308)

    pattern309 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule309 = ReplacementRule(pattern309, lambda b, e, c, m, d, n, a, x : -m*Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))) + x**m*(a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule309)

    pattern310 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule310 = ReplacementRule(pattern310, lambda b, e, c, m, d, n, a, x : m*Int(x**(m + S(-1))*(a + b*acot(c*x))**(n + S(1)), x)/(b*c*d*(n + S(1))) - x**m*(a + b*acot(c*x))**(n + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule310)

    pattern311 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda a, m: Not(NonzeroQ(a) & Equal(m, S(1)))))
    rule311 = ReplacementRule(pattern311, lambda b, e, c, m, d, a, x : Int(ExpandIntegrand(a + b*ArcTan(c*x), x**m/(d + e*x**S(2)), x), x))
    rubi.add(rule311)

    pattern312 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda a, m: Not(NonzeroQ(a) & Equal(m, S(1)))))
    rule312 = ReplacementRule(pattern312, lambda b, e, c, m, d, a, x : Int(ExpandIntegrand(a + b*acot(c*x), x**m/(d + e*x**S(2)), x), x))
    rubi.add(rule312)

    pattern313 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule313 = ReplacementRule(pattern313, lambda b, e, c, d, n, a, p, x : -b*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(S(2)*c*(p + S(1))) + (a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule313)

    pattern314 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule314 = ReplacementRule(pattern314, lambda b, e, c, d, n, a, p, x : b*n*Int((a + b*acot(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(S(2)*c*(p + S(1))) + (a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule314)

    pattern315 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule315 = ReplacementRule(pattern315, lambda b, e, c, d, n, a, x : x*(a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(d + e*x**S(2))*(n + S(1))) - S(4)*Int(x*(a + b*ArcTan(c*x))**(n + S(2))/(d + e*x**S(2))**S(2), x)/(b**S(2)*(n + S(1))*(n + S(2))) - (a + b*ArcTan(c*x))**(n + S(2))*(-c**S(2)*x**S(2) + S(1))/(b**S(2)*e*(d + e*x**S(2))*(n + S(1))*(n + S(2))))
    rubi.add(rule315)

    pattern316 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule316 = ReplacementRule(pattern316, lambda b, e, c, d, n, a, x : -x*(a + b*acot(c*x))**(n + S(1))/(b*c*d*(d + e*x**S(2))*(n + S(1))) - S(4)*Int(x*(a + b*acot(c*x))**(n + S(2))/(d + e*x**S(2))**S(2), x)/(b**S(2)*(n + S(1))*(n + S(2))) - (a + b*acot(c*x))**(n + S(2))*(-c**S(2)*x**S(2) + S(1))/(b**S(2)*e*(d + e*x**S(2))*(n + S(1))*(n + S(2))))
    rubi.add(rule316)

    pattern317 = Pattern(Integral(x_**S(2)*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-5)/2)))
    rule317 = ReplacementRule(pattern317, lambda b, e, c, d, a, p, x : -b*(d + e*x**S(2))**(p + S(1))/(S(4)*c**S(3)*d*(p + S(1))**S(2)) + x*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*c**S(2)*d*(p + S(1))) - Int((a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*c**S(2)*d*(p + S(1))))
    rubi.add(rule317)

    pattern318 = Pattern(Integral(x_**S(2)*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-5)/2)))
    rule318 = ReplacementRule(pattern318, lambda b, e, c, d, a, p, x : b*(d + e*x**S(2))**(p + S(1))/(S(4)*c**S(3)*d*(p + S(1))**S(2)) + x*(a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*c**S(2)*d*(p + S(1))) - Int((a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(S(2)*c**S(2)*d*(p + S(1))))
    rubi.add(rule318)

    pattern319 = Pattern(Integral(x_**S(2)*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule319 = ReplacementRule(pattern319, lambda b, e, c, d, n, a, x : b*n*Int(x*(a + b*ArcTan(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/(S(2)*c) - x*(a + b*ArcTan(c*x))**n/(S(2)*c**S(2)*d*(d + e*x**S(2))) + (a + b*ArcTan(c*x))**(n + S(1))/(S(2)*b*c**S(3)*d**S(2)*(n + S(1))))
    rubi.add(rule319)

    pattern320 = Pattern(Integral(x_**S(2)*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule320 = ReplacementRule(pattern320, lambda b, e, c, d, n, a, x : -b*n*Int(x*(a + b*acot(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x)/(S(2)*c) - x*(a + b*acot(c*x))**n/(S(2)*c**S(2)*d*(d + e*x**S(2))) - (a + b*acot(c*x))**(n + S(1))/(S(2)*b*c**S(3)*d**S(2)*(n + S(1))))
    rubi.add(rule320)

    pattern321 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule321 = ReplacementRule(pattern321, lambda b, e, c, m, d, a, p, x : b*x**m*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) - x**(m + S(-1))*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) + (m + S(-1))*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule321)

    pattern322 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule322 = ReplacementRule(pattern322, lambda b, e, c, m, d, a, p, x : -b*x**m*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) - x**(m + S(-1))*(a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) + (m + S(-1))*Int(x**(m + S(-2))*(a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule322)

    pattern323 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule323 = ReplacementRule(pattern323, lambda b, e, c, m, d, n, a, p, x : -b**S(2)*n*(n + S(-1))*Int(x**m*(a + b*ArcTan(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/m**S(2) + b*n*x**m*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) - x**(m + S(-1))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) + (m + S(-1))*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule323)

    pattern324 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(1))))
    rule324 = ReplacementRule(pattern324, lambda b, e, c, m, d, n, a, p, x : -b**S(2)*n*(n + S(-1))*Int(x**m*(a + b*acot(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x)/m**S(2) - b*n*x**m*(a + b*acot(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)) - x**(m + S(-1))*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m) + (m + S(-1))*Int(x**(m + S(-2))*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/(c**S(2)*d*m))
    rubi.add(rule324)

    pattern325 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule325 = ReplacementRule(pattern325, lambda b, e, c, m, d, n, a, p, x : -m*Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) + x**m*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule325)

    pattern326 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule326 = ReplacementRule(pattern326, lambda b, e, c, m, d, n, a, p, x : m*Int(x**(m + S(-1))*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) - x**m*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule326)

    pattern327 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule327 = ReplacementRule(pattern327, lambda b, e, c, m, d, n, a, p, x : -b*c*n*Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*(m + S(1))))
    rubi.add(rule327)

    pattern328 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p, m: ZeroQ(m + S(2)*p + S(3))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule328 = ReplacementRule(pattern328, lambda b, e, c, m, d, n, a, p, x : b*c*n*Int(x**(m + S(1))*(a + b*acot(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x)/(m + S(1)) + x**(m + S(1))*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*(m + S(1))))
    rubi.add(rule328)

    pattern329 = Pattern(Integral(x_**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: NonzeroQ(m + S(2))))
    rule329 = ReplacementRule(pattern329, lambda b, e, c, m, d, a, x : -b*c*d*Int(x**(m + S(1))/sqrt(d + e*x**S(2)), x)/(m + S(2)) + d*Int(x**m*(a + b*ArcTan(c*x))/sqrt(d + e*x**S(2)), x)/(m + S(2)) + x**(m + S(1))*(a + b*ArcTan(c*x))*sqrt(d + e*x**S(2))/(m + S(2)))
    rubi.add(rule329)

    pattern330 = Pattern(Integral(x_**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: NonzeroQ(m + S(2))))
    rule330 = ReplacementRule(pattern330, lambda b, e, c, m, d, a, x : b*c*d*Int(x**(m + S(1))/sqrt(d + e*x**S(2)), x)/(m + S(2)) + d*Int(x**m*(a + b*acot(c*x))/sqrt(d + e*x**S(2)), x)/(m + S(2)) + x**(m + S(1))*(a + b*acot(c*x))*sqrt(d + e*x**S(2))/(m + S(2)))
    rubi.add(rule330)

    pattern331 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p: Greater(p, S(1))))
    rule331 = ReplacementRule(pattern331, lambda b, e, c, m, d, n, a, p, x : Int(ExpandIntegrand(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule331)

    pattern332 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda p: Greater(p, S(1))))
    rule332 = ReplacementRule(pattern332, lambda b, e, c, m, d, n, a, p, x : Int(ExpandIntegrand(x**m*(a + b*acot(c*x))**n*(d + e*x**S(2))**p, x), x))
    rubi.add(rule332)

    pattern333 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, p, m: RationalQ(m) | (IntegerQ(p) & EqQ(n, S(1)))))
    rule333 = ReplacementRule(pattern333, lambda b, e, c, m, d, n, a, p, x : c**S(2)*d*Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x) + d*Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x))
    rubi.add(rule333)

    pattern334 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, p, m: RationalQ(m) | (IntegerQ(p) & EqQ(n, S(1)))))
    rule334 = ReplacementRule(pattern334, lambda b, e, c, m, d, n, a, p, x : c**S(2)*d*Int(x**(m + S(2))*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x) + d*Int(x**m*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x))
    rubi.add(rule334)

    pattern335 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule335 = ReplacementRule(pattern335, lambda b, e, c, m, d, n, a, x : -b*n*Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x)/(c*m) - (m + S(-1))*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n/sqrt(d + e*x**S(2)), x)/(c**S(2)*m) + x**(m + S(-1))*(a + b*ArcTan(c*x))**n*sqrt(d + e*x**S(2))/(c**S(2)*d*m))
    rubi.add(rule335)

    pattern336 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule336 = ReplacementRule(pattern336, lambda b, e, c, m, d, n, a, x : b*n*Int(x**(m + S(-1))*(a + b*acot(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x)/(c*m) - (m + S(-1))*Int(x**(m + S(-2))*(a + b*acot(c*x))**n/sqrt(d + e*x**S(2)), x)/(c**S(2)*m) + x**(m + S(-1))*(a + b*acot(c*x))**n*sqrt(d + e*x**S(2))/(c**S(2)*d*m))
    rubi.add(rule336)

    pattern337 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule337 = ReplacementRule(pattern337, lambda b, e, c, d, a, x : ImaginaryI*b*PolyLog(S(2), -sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/sqrt(d) - ImaginaryI*b*PolyLog(S(2), sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/sqrt(d) - S(2)*(a + b*ArcTan(c*x))*atanh(sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/sqrt(d))
    rubi.add(rule337)

    pattern338 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda d: PositiveQ(d)))
    rule338 = ReplacementRule(pattern338, lambda b, e, c, d, a, x : -ImaginaryI*b*PolyLog(S(2), -sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/sqrt(d) + ImaginaryI*b*PolyLog(S(2), sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/sqrt(d) - S(2)*(a + b*acot(c*x))*atanh(sqrt(ImaginaryI*c*x + S(1))/sqrt(-ImaginaryI*c*x + S(1)))/sqrt(d))
    rubi.add(rule338)

    pattern339 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule339 = ReplacementRule(pattern339, lambda b, e, c, d, n, a, x : Subst(Int((a + b*x)**n*csc(x), x), x, ArcTan(c*x))/sqrt(d))
    rubi.add(rule339)

    pattern340 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: PositiveQ(d)))
    rule340 = ReplacementRule(pattern340, lambda b, e, c, d, n, a, x : -c*x*sqrt(S(1) + S(1)/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*sec(x), x), x, acot(c*x))/sqrt(d + e*x**S(2)))
    rubi.add(rule340)

    pattern341 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule341 = ReplacementRule(pattern341, lambda b, e, c, d, n, a, x : sqrt(c**S(2)*x**S(2) + S(1))*Int((a + b*ArcTan(c*x))**n/(x*sqrt(c**S(2)*x**S(2) + S(1))), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule341)

    pattern342 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d: Not(PositiveQ(d))))
    rule342 = ReplacementRule(pattern342, lambda b, e, c, d, n, a, x : sqrt(c**S(2)*x**S(2) + S(1))*Int((a + b*acot(c*x))**n/(x*sqrt(c**S(2)*x**S(2) + S(1))), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule342)

    pattern343 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule343 = ReplacementRule(pattern343, lambda b, e, c, d, n, a, x : b*c*n*Int((a + b*ArcTan(c*x))**(n + S(-1))/(x*sqrt(d + e*x**S(2))), x) - (a + b*ArcTan(c*x))**n*sqrt(d + e*x**S(2))/(d*x))
    rubi.add(rule343)

    pattern344 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(x_**S(2)*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))))
    rule344 = ReplacementRule(pattern344, lambda b, e, c, d, n, a, x : -b*c*n*Int((a + b*acot(c*x))**(n + S(-1))/(x*sqrt(d + e*x**S(2))), x) - (a + b*acot(c*x))**n*sqrt(d + e*x**S(2))/(d*x))
    rubi.add(rule344)

    pattern345 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: Unequal(m, S(-2))))
    rule345 = ReplacementRule(pattern345, lambda b, e, c, m, d, n, a, x : -b*c*n*Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x)/(m + S(1)) - c**S(2)*(m + S(2))*Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n/sqrt(d + e*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcTan(c*x))**n*sqrt(d + e*x**S(2))/(d*(m + S(1))))
    rubi.add(rule345)

    pattern346 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: Unequal(m, S(-2))))
    rule346 = ReplacementRule(pattern346, lambda b, e, c, m, d, n, a, x : b*c*n*Int(x**(m + S(1))*(a + b*acot(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x)/(m + S(1)) - c**S(2)*(m + S(2))*Int(x**(m + S(2))*(a + b*acot(c*x))**n/sqrt(d + e*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*acot(c*x))**n*sqrt(d + e*x**S(2))/(d*(m + S(1))))
    rubi.add(rule346)

    pattern347 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule347 = ReplacementRule(pattern347, lambda b, e, c, m, d, n, a, p, x : -d*Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x)/e + Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/e)
    rubi.add(rule347)

    pattern348 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule348 = ReplacementRule(pattern348, lambda b, e, c, m, d, n, a, p, x : -d*Int(x**(m + S(-2))*(a + b*acot(c*x))**n*(d + e*x**S(2))**p, x)/e + Int(x**(m + S(-2))*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/e)
    rubi.add(rule348)

    pattern349 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Less(m, S(0))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule349 = ReplacementRule(pattern349, lambda b, e, c, m, d, n, a, p, x : -e*Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x)/d + Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/d)
    rubi.add(rule349)

    pattern350 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p, m: IntegersQ(m, n, S(2)*p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: Less(m, S(0))), CustomConstraint(lambda n: Unequal(n, S(-1))))
    rule350 = ReplacementRule(pattern350, lambda b, e, c, m, d, n, a, p, x : -e*Int(x**(m + S(2))*(a + b*acot(c*x))**n*(d + e*x**S(2))**p, x)/d + Int(x**m*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x)/d)
    rubi.add(rule350)

    pattern351 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(2))))
    rule351 = ReplacementRule(pattern351, lambda b, e, c, m, d, n, a, p, x : -c*(m + S(2)*p + S(2))*Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) - m*Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) + x**m*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule351)

    pattern352 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, p, m: RationalQ(m, n, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda p, m: NonzeroQ(m + S(2)*p + S(2))))
    rule352 = ReplacementRule(pattern352, lambda b, e, c, m, d, n, a, p, x : c*(m + S(2)*p + S(2))*Int(x**(m + S(1))*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*(n + S(1))) + m*Int(x**(m + S(-1))*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x)/(b*c*(n + S(1))) - x**m*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))))
    rubi.add(rule352)

    pattern353 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda p, d: IntegerQ(p) | PositiveQ(d)))
    rule353 = ReplacementRule(pattern353, lambda b, e, c, m, d, n, a, p, x : c**(-m + S(-1))*d**p*Subst(Int((a + b*x)**n*sin(x)**m*cos(x)**(-m - S(2)*p + S(-2)), x), x, ArcTan(c*x)))
    rubi.add(rule353)

    pattern354 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda p, d: Not(IntegerQ(p) | PositiveQ(d))))
    rule354 = ReplacementRule(pattern354, lambda b, e, c, m, d, n, a, p, x : d**(p + S(1)/2)*sqrt(c**S(2)*x**S(2) + S(1))*Int(x**m*(a + b*ArcTan(c*x))**n*(c**S(2)*x**S(2) + S(1))**p, x)/sqrt(d + e*x**S(2)))
    rubi.add(rule354)

    pattern355 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda p: IntegerQ(p)))
    rule355 = ReplacementRule(pattern355, lambda b, e, c, m, d, n, a, p, x : -c**(-m + S(-1))*d**p*Subst(Int((a + b*x)**n*sin(x)**(-m - S(2)*p + S(-2))*cos(x)**m, x), x, acot(c*x)))
    rubi.add(rule355)

    pattern356 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p, m: NegativeIntegerQ(m + S(2)*p + S(1))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule356 = ReplacementRule(pattern356, lambda b, e, c, m, d, n, a, p, x : -c**(-m)*d**(p + S(1)/2)*x*sqrt((c**S(2)*x**S(2) + S(1))/(c**S(2)*x**S(2)))*Subst(Int((a + b*x)**n*sin(x)**(-m - S(2)*p + S(-2))*cos(x)**m, x), x, acot(c*x))/sqrt(d + e*x**S(2)))
    rubi.add(rule356)

    pattern357 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule357 = ReplacementRule(pattern357, lambda b, e, c, d, a, p, x : -b*c*Int((d + e*x**S(2))**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x)/(S(2)*e*(p + S(1))) + (a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule357)

    pattern358 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule358 = ReplacementRule(pattern358, lambda b, e, c, d, a, p, x : b*c*Int((d + e*x**S(2))**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x)/(S(2)*e*(p + S(1))) + (a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule358)

    pattern359 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & Not(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))), )
    def With359(b, e, c, m, d, a, p, x):
        u = IntHide(x**m*(d + e*x**S(2))**p, x)
        return -b*c*Int(SimplifyIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcTan(c*x), u, x)
    rule359 = ReplacementRule(pattern359, lambda b, e, c, m, d, a, p, x : With359(b, e, c, m, d, a, p, x))
    rubi.add(rule359)

    pattern360 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & Not(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))), )
    def With360(b, e, c, m, d, a, p, x):
        u = IntHide(x**m*(d + e*x**S(2))**p, x)
        return b*c*Int(SimplifyIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acot(c*x), u, x)
    rule360 = ReplacementRule(pattern360, lambda b, e, c, m, d, a, p, x : With360(b, e, c, m, d, a, p, x))
    rubi.add(rule360)

    pattern361 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (IntegerQ(m) & Less(p, S(-1)) & Unequal(m, S(1)))))
    rule361 = ReplacementRule(pattern361, lambda b, e, c, m, d, n, a, p, x : Int(ExpandIntegrand((a + b*ArcTan(c*x))**n, x**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule361)

    pattern362 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: Greater(p, S(0)) | (IntegerQ(m) & Less(p, S(-1)) & Unequal(m, S(1)))))
    rule362 = ReplacementRule(pattern362, lambda b, e, c, m, d, n, a, p, x : Int(ExpandIntegrand((a + b*acot(c*x))**n, x**m*(d + e*x**S(2))**p, x), x))
    rubi.add(rule362)

    pattern363 = Pattern(Integral(x_**WC('m', S(1))*(a_ + ArcTan(x_*WC('c', S(1)))*WC('b', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule363 = ReplacementRule(pattern363, lambda b, e, c, m, d, a, p, x : a*Int(x**m*(d + e*x**S(2))**p, x) + b*Int(x**m*(d + e*x**S(2))**p*ArcTan(c*x), x))
    rubi.add(rule363)

    pattern364 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*acot(x_*WC('c', S(1))))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule364 = ReplacementRule(pattern364, lambda b, e, c, m, d, a, p, x : a*Int(x**m*(d + e*x**S(2))**p, x) + b*Int(x**m*(d + e*x**S(2))**p*acot(c*x), x))
    rubi.add(rule364)

    pattern365 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule365 = ReplacementRule(pattern365, lambda b, e, c, m, d, n, a, p, x : Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule365)

    pattern366 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule366 = ReplacementRule(pattern366, lambda b, e, c, m, d, n, a, p, x : Int(x**m*(a + b*acot(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule366)

    pattern367 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*atanh(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule367 = ReplacementRule(pattern367, lambda b, e, c, d, n, a, u, x : -Int((a + b*ArcTan(c*x))**n*log(-u + S(1))/(d + e*x**S(2)), x)/S(2) + Int((a + b*ArcTan(c*x))**n*log(u + S(1))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule367)

    pattern368 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*acoth(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule368 = ReplacementRule(pattern368, lambda b, e, c, d, n, a, u, x : -Int((a + b*acot(c*x))**n*log(SimplifyIntegrand(S(1) - S(1)/u, x))/(d + e*x**S(2)), x)/S(2) + Int((a + b*acot(c*x))**n*log(SimplifyIntegrand(S(1) + 1/u, x))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule368)

    pattern369 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*atanh(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule369 = ReplacementRule(pattern369, lambda b, e, c, d, n, a, u, x : -Int((a + b*ArcTan(c*x))**n*log(-u + S(1))/(d + e*x**S(2)), x)/S(2) + Int((a + b*ArcTan(c*x))**n*log(u + S(1))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule369)

    pattern370 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*acoth(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule370 = ReplacementRule(pattern370, lambda b, e, c, d, n, a, u, x : -Int((a + b*acot(c*x))**n*log(SimplifyIntegrand(S(1) - S(1)/u, x))/(d + e*x**S(2)), x)/S(2) + Int((a + b*acot(c*x))**n*log(SimplifyIntegrand(S(1) + 1/u, x))/(d + e*x**S(2)), x)/S(2))
    rubi.add(rule370)

    pattern371 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ((-u + S(1))**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule371 = ReplacementRule(pattern371, lambda b, e, c, d, n, a, u, x : -ImaginaryI*b*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(S(2), -u + S(1))/(d + e*x**S(2)), x)/S(2) + ImaginaryI*(a + b*ArcTan(c*x))**n*PolyLog(S(2), -u + S(1))/(S(2)*c*d))
    rubi.add(rule371)

    pattern372 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ((-u + S(1))**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule372 = ReplacementRule(pattern372, lambda b, e, c, d, n, a, u, x : ImaginaryI*b*n*Int((a + b*acot(c*x))**(n + S(-1))*PolyLog(S(2), -u + S(1))/(d + e*x**S(2)), x)/S(2) + ImaginaryI*(a + b*acot(c*x))**n*PolyLog(S(2), -u + S(1))/(S(2)*c*d))
    rubi.add(rule372)

    pattern373 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ((-u + S(1))**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule373 = ReplacementRule(pattern373, lambda b, e, c, d, n, a, u, x : ImaginaryI*b*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(S(2), -u + S(1))/(d + e*x**S(2)), x)/S(2) - ImaginaryI*(a + b*ArcTan(c*x))**n*PolyLog(S(2), -u + S(1))/(S(2)*c*d))
    rubi.add(rule373)

    pattern374 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ((-u + S(1))**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule374 = ReplacementRule(pattern374, lambda b, e, c, d, n, a, u, x : -ImaginaryI*b*n*Int((a + b*acot(c*x))**(n + S(-1))*PolyLog(S(2), -u + S(1))/(d + e*x**S(2)), x)/S(2) - ImaginaryI*(a + b*acot(c*x))**n*PolyLog(S(2), -u + S(1))/(S(2)*c*d))
    rubi.add(rule374)

    pattern375 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule375 = ReplacementRule(pattern375, lambda b, e, c, d, n, a, p, u, x : ImaginaryI*b*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) - ImaginaryI*(a + b*ArcTan(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule375)

    pattern376 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI + c*x) + S(1))**S(2))))
    rule376 = ReplacementRule(pattern376, lambda b, e, c, d, n, a, p, u, x : -ImaginaryI*b*n*Int((a + b*acot(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) - ImaginaryI*(a + b*acot(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule376)

    pattern377 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule377 = ReplacementRule(pattern377, lambda b, e, c, d, n, a, p, u, x : -ImaginaryI*b*n*Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) + ImaginaryI*(a + b*ArcTan(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule377)

    pattern378 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda u, c, x: ZeroQ(u**S(2) - (-S(2)*ImaginaryI/(ImaginaryI - c*x) + S(1))**S(2))))
    rule378 = ReplacementRule(pattern378, lambda b, e, c, d, n, a, p, u, x : ImaginaryI*b*n*Int((a + b*acot(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x)/S(2) + ImaginaryI*(a + b*acot(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d))
    rubi.add(rule378)

    pattern379 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)))
    rule379 = ReplacementRule(pattern379, lambda b, e, c, d, a, x : (log(a + b*ArcTan(c*x)) - log(a + b*acot(c*x)))/(b*c*d*(S(2)*a + b*ArcTan(c*x) + b*acot(c*x))))
    rubi.add(rule379)

    pattern380 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('m', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Inequality(S(0), Less, n, LessEqual, m)))
    rule380 = ReplacementRule(pattern380, lambda b, e, c, m, d, n, a, x : n*Int((a + b*ArcTan(c*x))**(n + S(-1))*(a + b*acot(c*x))**(m + S(1))/(d + e*x**S(2)), x)/(m + S(1)) - (a + b*ArcTan(c*x))**n*(a + b*acot(c*x))**(m + S(1))/(b*c*d*(m + S(1))))
    rubi.add(rule380)

    pattern381 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, c, d: ZeroQ(-c**S(2)*d + e)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Less(S(0), n, m)))
    rule381 = ReplacementRule(pattern381, lambda b, e, c, m, d, n, a, x : n*Int((a + b*ArcTan(c*x))**(m + S(1))*(a + b*acot(c*x))**(n + S(-1))/(d + e*x**S(2)), x)/(m + S(1)) + (a + b*ArcTan(c*x))**(m + S(1))*(a + b*acot(c*x))**n/(b*c*d*(m + S(1))))
    rubi.add(rule381)

    pattern382 = Pattern(Integral(ArcTan(x_*WC('a', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n, a, c, d: Not(Equal(n, S(2)) & ZeroQ(-a**S(2)*c + d))))
    rule382 = ReplacementRule(pattern382, lambda c, d, n, a, x : ImaginaryI*Int(log(-ImaginaryI*a*x + S(1))/(c + d*x**n), x)/S(2) - ImaginaryI*Int(log(ImaginaryI*a*x + S(1))/(c + d*x**n), x)/S(2))
    rubi.add(rule382)

    pattern383 = Pattern(Integral(acot(x_*WC('a', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n, a, c, d: Not(Equal(n, S(2)) & ZeroQ(-a**S(2)*c + d))))
    rule383 = ReplacementRule(pattern383, lambda c, d, n, a, x : ImaginaryI*Int(log(-ImaginaryI/(a*x) + S(1))/(c + d*x**n), x)/S(2) - ImaginaryI*Int(log(ImaginaryI/(a*x) + S(1))/(c + d*x**n), x)/S(2))
    rubi.add(rule383)

    pattern384 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)))
    rule384 = ReplacementRule(pattern384, lambda b, e, c, d, g, a, f, x : -b*c*Int(x*(d + e*log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x) - S(2)*e*g*Int(x**S(2)*(a + b*ArcTan(c*x))/(f + g*x**S(2)), x) + x*(a + b*ArcTan(c*x))*(d + e*log(f + g*x**S(2))))
    rubi.add(rule384)

    pattern385 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)))
    rule385 = ReplacementRule(pattern385, lambda b, e, c, d, g, a, f, x : b*c*Int(x*(d + e*log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x) - S(2)*e*g*Int(x**S(2)*(a + b*acot(c*x))/(f + g*x**S(2)), x) + x*(a + b*acot(c*x))*(d + e*log(f + g*x**S(2))))
    rubi.add(rule385)

    pattern386 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2))))
    rule386 = ReplacementRule(pattern386, lambda b, e, c, m, d, g, a, f, x : -b*c*Int(x**(m + S(1))*(d + e*log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) - S(2)*e*g*Int(x**(m + S(2))*(a + b*ArcTan(c*x))/(f + g*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*ArcTan(c*x))*(d + e*log(f + g*x**S(2)))/(m + S(1)))
    rubi.add(rule386)

    pattern387 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: NegativeIntegerQ(m/S(2))))
    rule387 = ReplacementRule(pattern387, lambda b, e, c, m, d, g, a, f, x : b*c*Int(x**(m + S(1))*(d + e*log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x)/(m + S(1)) - S(2)*e*g*Int(x**(m + S(2))*(a + b*acot(c*x))/(f + g*x**S(2)), x)/(m + S(1)) + x**(m + S(1))*(a + b*acot(c*x))*(d + e*log(f + g*x**S(2)))/(m + S(1)))
    rubi.add(rule387)

    pattern388 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m/S(2) + S(1)/2)), )
    def With388(b, e, c, m, d, g, a, f, x):
        u = IntHide(x**m*(d + e*log(f + g*x**S(2))), x)
        return -b*c*Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*ArcTan(c*x), u, x)
    rule388 = ReplacementRule(pattern388, lambda b, e, c, m, d, g, a, f, x : With388(b, e, c, m, d, g, a, f, x))
    rubi.add(rule388)

    pattern389 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m/S(2) + S(1)/2)), )
    def With389(b, e, c, m, d, g, a, f, x):
        u = IntHide(x**m*(d + e*log(f + g*x**S(2))), x)
        return b*c*Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x) + Dist(a + b*acot(c*x), u, x)
    rule389 = ReplacementRule(pattern389, lambda b, e, c, m, d, g, a, f, x : With389(b, e, c, m, d, g, a, f, x))
    rubi.add(rule389)

    pattern390 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Unequal(m, S(-1))), )
    def With390(b, e, c, m, d, g, a, f, x):
        u = IntHide(x**m*(a + b*ArcTan(c*x)), x)
        return -S(2)*e*g*Int(ExpandIntegrand(u*x/(f + g*x**S(2)), x), x) + Dist(d + e*log(f + g*x**S(2)), u, x)
    rule390 = ReplacementRule(pattern390, lambda b, e, c, m, d, g, a, f, x : With390(b, e, c, m, d, g, a, f, x))
    rubi.add(rule390)

    pattern391 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Unequal(m, S(-1))), )
    def With391(b, e, c, m, d, g, a, f, x):
        u = IntHide(x**m*(a + b*acot(c*x)), x)
        return -S(2)*e*g*Int(ExpandIntegrand(u*x/(f + g*x**S(2)), x), x) + Dist(d + e*log(f + g*x**S(2)), u, x)
    rule391 = ReplacementRule(pattern391, lambda b, e, c, m, d, g, a, f, x : With391(b, e, c, m, d, g, a, f, x))
    rubi.add(rule391)

    pattern392 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**S(2)*(WC('d', S(0)) + WC('e', S(1))*log(f_ + x_**S(2)*WC('g', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda f, g, c: ZeroQ(-c**S(2)*f + g)))
    rule392 = ReplacementRule(pattern392, lambda b, e, c, d, g, a, f, x : b*c*e*Int(x**S(2)*(a + b*ArcTan(c*x))/(c**S(2)*x**S(2) + S(1)), x) - b*Int((a + b*ArcTan(c*x))*(d + e*log(f + g*x**S(2))), x)/c - e*x**S(2)*(a + b*ArcTan(c*x))**S(2)/S(2) + (a + b*ArcTan(c*x))**S(2)*(d + e*log(f + g*x**S(2)))*(f + g*x**S(2))/(S(2)*g))
    rubi.add(rule392)

    pattern393 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**S(2)*(WC('d', S(0)) + WC('e', S(1))*log(f_ + x_**S(2)*WC('g', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda f, g, c: ZeroQ(-c**S(2)*f + g)))
    rule393 = ReplacementRule(pattern393, lambda b, e, c, d, g, a, f, x : -b*c*e*Int(x**S(2)*(a + b*acot(c*x))/(c**S(2)*x**S(2) + S(1)), x) + b*Int((a + b*acot(c*x))*(d + e*log(f + g*x**S(2))), x)/c - e*x**S(2)*(a + b*acot(c*x))**S(2)/S(2) + (a + b*acot(c*x))**S(2)*(d + e*log(f + g*x**S(2)))*(f + g*x**S(2))/(S(2)*g))
    rubi.add(rule393)

    pattern394 = Pattern(Integral(exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: OddQ(ImaginaryI*n)))
    rule394 = ReplacementRule(pattern394, lambda n, a, x : Int((-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + S(1)/2)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + S(1)/2)/sqrt(a**S(2)*x**S(2) + S(1)), x))
    rubi.add(rule394)

    pattern395 = Pattern(Integral(x_**WC('m', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: OddQ(ImaginaryI*n)))
    rule395 = ReplacementRule(pattern395, lambda n, a, m, x : Int(x**m*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + S(1)/2)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + S(1)/2)/sqrt(a**S(2)*x**S(2) + S(1)), x))
    rubi.add(rule395)

    pattern396 = Pattern(Integral(exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(OddQ(ImaginaryI*n))))
    rule396 = ReplacementRule(pattern396, lambda n, a, x : Int((-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule396)

    pattern397 = Pattern(Integral(x_**WC('m', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(OddQ(ImaginaryI*n))))
    rule397 = ReplacementRule(pattern397, lambda n, a, m, x : Int(x**m*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule397)

    pattern398 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*c**S(2) + d**S(2))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)))
    rule398 = ReplacementRule(pattern398, lambda c, d, n, a, p, u, x : c**p*Int(u*(S(1) + d*x/c)**p*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule398)

    pattern399 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*c**S(2) + d**S(2))), CustomConstraint(lambda p, c: Not(IntegerQ(p) | PositiveQ(c))))
    rule399 = ReplacementRule(pattern399, lambda c, d, n, a, p, u, x : Int(u*(c + d*x)**p*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule399)

    pattern400 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule400 = ReplacementRule(pattern400, lambda c, d, n, a, p, u, x : d**p*Int(u*x**(-p)*(c*x/d + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule400)

    pattern401 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))), CustomConstraint(lambda c: PositiveQ(c)))
    rule401 = ReplacementRule(pattern401, lambda c, d, n, a, p, u, x : (S(-1))**(n/S(2))*c**p*Int(u*(S(1) - S(1)/(ImaginaryI*a*x))**(ImaginaryI*n/S(2))*(S(1) + S(1)/(ImaginaryI*a*x))**(-ImaginaryI*n/S(2))*(S(1) + d/(c*x))**p, x))
    rubi.add(rule401)

    pattern402 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))), CustomConstraint(lambda c: Not(PositiveQ(c))))
    rule402 = ReplacementRule(pattern402, lambda c, d, n, a, p, u, x : Int(u*(c + d/x)**p*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule402)

    pattern403 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule403 = ReplacementRule(pattern403, lambda c, d, n, a, p, u, x : x**p*(c + d/x)**p*(c*x/d + S(1))**(-p)*Int(u*x**(-p)*(c*x/d + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule403)

    pattern404 = Pattern(Integral(exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n))))
    rule404 = ReplacementRule(pattern404, lambda c, d, n, a, x : (a*x + n)*exp(n*ArcTan(a*x))/(a*c*sqrt(c + d*x**S(2))*(n**S(2) + S(1))))
    rubi.add(rule404)

    pattern405 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n))), CustomConstraint(lambda n, p: NonzeroQ(n**S(2) + S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule405 = ReplacementRule(pattern405, lambda c, d, n, a, p, x : S(2)*(p + S(1))*(S(2)*p + S(3))*Int((c + d*x**S(2))**(p + S(1))*exp(n*ArcTan(a*x)), x)/(c*(n**S(2) + S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(-S(2)*a*x*(p + S(1)) + n)*exp(n*ArcTan(a*x))/(a*c*(n**S(2) + S(4)*(p + S(1))**S(2))))
    rubi.add(rule405)

    pattern406 = Pattern(Integral(exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)))
    rule406 = ReplacementRule(pattern406, lambda c, d, n, a, x : exp(n*ArcTan(a*x))/(a*c*n))
    rubi.add(rule406)

    pattern407 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2) + S(1)/2)), CustomConstraint(lambda n, p: Not(IntegerQ(-ImaginaryI*n/S(2) + p))))
    rule407 = ReplacementRule(pattern407, lambda c, d, n, a, p, x : c**p*Int((a**S(2)*x**S(2) + S(1))**(-ImaginaryI*n/S(2) + p)*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n), x))
    rubi.add(rule407)

    pattern408 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)))
    rule408 = ReplacementRule(pattern408, lambda c, d, n, a, p, x : c**p*Int((-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + p), x))
    rubi.add(rule408)

    pattern409 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: PositiveIntegerQ(ImaginaryI*n/S(2))))
    rule409 = ReplacementRule(pattern409, lambda c, d, n, a, p, x : c**(ImaginaryI*n/S(2))*Int((c + d*x**S(2))**(-ImaginaryI*n/S(2) + p)*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n), x))
    rubi.add(rule409)

    pattern410 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: NegativeIntegerQ(ImaginaryI*n/S(2))))
    rule410 = ReplacementRule(pattern410, lambda c, d, n, a, p, x : c**(-ImaginaryI*n/S(2))*Int((c + d*x**S(2))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n), x))
    rubi.add(rule410)

    pattern411 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: Not(IntegerQ(p) | PositiveQ(c))))
    rule411 = ReplacementRule(pattern411, lambda c, d, n, a, p, x : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(a**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int((a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule411)

    pattern412 = Pattern(Integral(x_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n))))
    rule412 = ReplacementRule(pattern412, lambda c, d, n, a, x : (a*n*x + S(-1))*exp(n*ArcTan(a*x))/(d*sqrt(c + d*x**S(2))*(n**S(2) + S(1))))
    rubi.add(rule412)

    pattern413 = Pattern(Integral(x_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule413 = ReplacementRule(pattern413, lambda c, d, n, a, p, x : -a*c*n*Int((c + d*x**S(2))**p*exp(n*ArcTan(a*x)), x)/(S(2)*d*(p + S(1))) + (c + d*x**S(2))**(p + S(1))*exp(n*ArcTan(a*x))/(S(2)*d*(p + S(1))))
    rubi.add(rule413)

    pattern414 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n, p: ZeroQ(n**S(2) - S(2)*p + S(-2))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n))))
    rule414 = ReplacementRule(pattern414, lambda c, d, n, a, p, x : (c + d*x**S(2))**(p + S(1))*(a*n*x + S(-1))*exp(n*ArcTan(a*x))/(a*d*n*(n**S(2) + S(1))))
    rubi.add(rule414)

    pattern415 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n))), CustomConstraint(lambda n, p: NonzeroQ(n**S(2) + S(4)*(p + S(1))**S(2))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule415 = ReplacementRule(pattern415, lambda c, d, n, a, p, x : (n**S(2) - S(2)*p + S(-2))*Int((c + d*x**S(2))**(p + S(1))*exp(n*ArcTan(a*x)), x)/(d*(n**S(2) + S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) - n)*exp(n*ArcTan(a*x))/(a*d*(n**S(2) + S(4)*(p + S(1))**S(2))))
    rubi.add(rule415)

    pattern416 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2) + S(1)/2)), CustomConstraint(lambda n, p: Not(IntegerQ(-ImaginaryI*n/S(2) + p))))
    rule416 = ReplacementRule(pattern416, lambda c, m, d, n, a, p, x : c**p*Int(x**m*(a**S(2)*x**S(2) + S(1))**(-ImaginaryI*n/S(2) + p)*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n), x))
    rubi.add(rule416)

    pattern417 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)))
    rule417 = ReplacementRule(pattern417, lambda c, m, d, n, a, p, x : c**p*Int(x**m*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + p), x))
    rubi.add(rule417)

    pattern418 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: PositiveIntegerQ(ImaginaryI*n/S(2))))
    rule418 = ReplacementRule(pattern418, lambda c, m, d, n, a, p, x : c**(ImaginaryI*n/S(2))*Int(x**m*(c + d*x**S(2))**(-ImaginaryI*n/S(2) + p)*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n), x))
    rubi.add(rule418)

    pattern419 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: NegativeIntegerQ(ImaginaryI*n/S(2))))
    rule419 = ReplacementRule(pattern419, lambda c, m, d, n, a, p, x : c**(-ImaginaryI*n/S(2))*Int(x**m*(c + d*x**S(2))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n), x))
    rubi.add(rule419)

    pattern420 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: Not(IntegerQ(p) | PositiveQ(c))))
    rule420 = ReplacementRule(pattern420, lambda c, m, d, n, a, p, x : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(a**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(x**m*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule420)

    pattern421 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)))
    rule421 = ReplacementRule(pattern421, lambda c, d, n, a, p, u, x : c**p*Int(u*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + p), x))
    rubi.add(rule421)

    pattern422 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))))
    rule422 = ReplacementRule(pattern422, lambda c, d, n, a, p, u, x : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(-ImaginaryI*a*x + S(1))**(-FracPart(p))*(ImaginaryI*a*x + S(1))**(-FracPart(p))*Int(u*(-ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(-ImaginaryI*n/S(2) + p), x))
    rubi.add(rule422)

    pattern423 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p, c: Not(IntegerQ(p) | PositiveQ(c))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))))
    rule423 = ReplacementRule(pattern423, lambda c, d, n, a, p, u, x : c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(a**S(2)*x**S(2) + S(1))**(-FracPart(p))*Int(u*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule423)

    pattern424 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda p: IntegerQ(p)))
    rule424 = ReplacementRule(pattern424, lambda c, d, n, a, p, u, x : d**p*Int(u*x**(-S(2)*p)*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule424)

    pattern425 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))), CustomConstraint(lambda c: PositiveQ(c)))
    rule425 = ReplacementRule(pattern425, lambda c, d, n, a, p, u, x : c**p*Int(u*(-ImaginaryI/(a*x) + S(1))**p*(ImaginaryI/(a*x) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule425)

    pattern426 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))), CustomConstraint(lambda c: Not(PositiveQ(c))))
    rule426 = ReplacementRule(pattern426, lambda c, d, n, a, p, u, x : x**(S(2)*p)*(c + d/x**S(2))**p*(-ImaginaryI*a*x + S(1))**(-p)*(ImaginaryI*a*x + S(1))**(-p)*Int(u*x**(-S(2)*p)*(-ImaginaryI*a*x + S(1))**p*(ImaginaryI*a*x + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule426)

    pattern427 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))))
    rule427 = ReplacementRule(pattern427, lambda c, d, n, a, p, u, x : x**(S(2)*p)*(c + d/x**S(2))**p*(a**S(2)*x**S(2) + S(1))**(-p)*Int(u*x**(-S(2)*p)*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule427)

    pattern428 = Pattern(Integral(exp(ArcTan((a_ + x_*WC('b', S(1)))*WC('c', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule428 = ReplacementRule(pattern428, lambda b, c, n, a, x : Int((-ImaginaryI*a*c - ImaginaryI*b*c*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule428)

    pattern429 = Pattern(Integral(x_**m_*exp(n_*ArcTan((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda n: RationalQ(ImaginaryI*n)), CustomConstraint(lambda n: Less(S(-1), ImaginaryI*n, S(1))))
    rule429 = ReplacementRule(pattern429, lambda b, c, m, n, a, x : S(4)*ImaginaryI**(-m)*b**(-m + S(-1))*c**(-m + S(-1))*Subst(Int(x**(S(2)/(ImaginaryI*n))*(x**(S(2)/(ImaginaryI*n)) + S(1))**(-m + S(-2))*(-ImaginaryI*a*c - x**(S(2)/(ImaginaryI*n))*(ImaginaryI*a*c + S(1)) + S(1))**m, x), x, (-ImaginaryI*c*(a + b*x) + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*c*(a + b*x) + S(1))**(-ImaginaryI*n/S(2)))/n)
    rubi.add(rule429)

    pattern430 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*exp(ArcTan((a_ + x_*WC('b', S(1)))*WC('c', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule430 = ReplacementRule(pattern430, lambda b, e, c, m, d, n, a, x : Int((d + e*x)**m*(-ImaginaryI*a*c - ImaginaryI*b*c*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(-ImaginaryI*n/S(2)), x))
    rubi.add(rule430)

    pattern431 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(a_ + x_*WC('b', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, a, e, d: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda b, a, e, c: ZeroQ(b**S(2)*c - e*(a**S(2) + S(1)))), CustomConstraint(lambda a, p, c: IntegerQ(p) | PositiveQ(c/(a**S(2) + S(1)))))
    rule431 = ReplacementRule(pattern431, lambda b, e, c, d, n, a, p, u, x : (c/(a**S(2) + S(1)))**p*Int(u*(-ImaginaryI*a - ImaginaryI*b*x + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*a + ImaginaryI*b*x + S(1))**(-ImaginaryI*n/S(2) + p), x))
    rubi.add(rule431)

    pattern432 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(a_ + x_*WC('b', S(1)))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, a, e, d: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda b, a, e, c: ZeroQ(b**S(2)*c - e*(a**S(2) + S(1)))), CustomConstraint(lambda a, p, c: Not(IntegerQ(p) | PositiveQ(c/(a**S(2) + S(1))))))
    rule432 = ReplacementRule(pattern432, lambda b, e, c, d, n, a, p, u, x : (c + d*x + e*x**S(2))**p*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**(-p)*Int(u*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x))
    rubi.add(rule432)

    pattern433 = Pattern(Integral(WC('u', S(1))*exp(ArcTan(WC('c', S(1))/(x_*WC('b', S(1)) + WC('a', S(0))))*WC('n', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule433 = ReplacementRule(pattern433, lambda b, c, n, a, u, x : Int(u*exp(n*acot(a/c + b*x/c)), x))
    rubi.add(rule433)

    pattern434 = Pattern(Integral(WC('u', S(1))*exp(n_*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))))
    rule434 = ReplacementRule(pattern434, lambda n, a, x, u : (S(-1))**(ImaginaryI*n/S(2))*Int(u*exp(-n*ArcTan(a*x)), x))
    rubi.add(rule434)

    pattern435 = Pattern(Integral(exp(n_*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: OddQ(ImaginaryI*n)))
    rule435 = ReplacementRule(pattern435, lambda n, a, x : -Subst(Int((-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + S(1)/2)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + S(1)/2)/(x**S(2)*sqrt(S(1) + x**S(2)/a**S(2))), x), x, 1/x))
    rubi.add(rule435)

    pattern436 = Pattern(Integral(x_**WC('m', S(1))*exp(n_*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n: OddQ(ImaginaryI*n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule436 = ReplacementRule(pattern436, lambda n, a, m, x : -Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + S(1)/2)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + S(1)/2)/sqrt(S(1) + x**S(2)/a**S(2)), x), x, 1/x))
    rubi.add(rule436)

    pattern437 = Pattern(Integral(exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n))))
    rule437 = ReplacementRule(pattern437, lambda n, a, x : -Subst(Int((-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2))/x**S(2), x), x, 1/x))
    rubi.add(rule437)

    pattern438 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n))), CustomConstraint(lambda m: IntegerQ(m)))
    rule438 = ReplacementRule(pattern438, lambda n, a, m, x : -Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(n/S(2))*(ImaginaryI*x/a + S(1))**(-n/S(2)), x), x, 1/x))
    rubi.add(rule438)

    pattern439 = Pattern(Integral(x_**m_*exp(n_*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: OddQ(ImaginaryI*n)), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule439 = ReplacementRule(pattern439, lambda n, a, x, m : -x**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + S(1)/2)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + S(1)/2)/sqrt(S(1) + x**S(2)/a**S(2)), x), x, 1/x))
    rubi.add(rule439)

    pattern440 = Pattern(Integral(x_**m_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule440 = ReplacementRule(pattern440, lambda n, a, x, m : -Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(n/S(2))*(ImaginaryI*x/a + S(1))**(-n/S(2)), x), x, 1/x))
    rubi.add(rule440)

    pattern441 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*c**S(2) + d**S(2))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p: IntegerQ(p)))
    rule441 = ReplacementRule(pattern441, lambda c, d, n, a, p, u, x : d**p*Int(u*x**p*(c/(d*x) + S(1))**p*exp(n*acot(a*x)), x))
    rubi.add(rule441)

    pattern442 = Pattern(Integral((c_ + x_*WC('d', S(1)))**p_*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*c**S(2) + d**S(2))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule442 = ReplacementRule(pattern442, lambda c, d, n, a, p, u, x : x**(-p)*(c + d*x)**p*(c/(d*x) + S(1))**(-p)*Int(u*x**p*(c/(d*x) + S(1))**p*exp(n*acot(a*x)), x))
    rubi.add(rule442)

    pattern443 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)))
    rule443 = ReplacementRule(pattern443, lambda c, d, n, a, p, x : -c**p*Subst(Int((S(1) + d*x/c)**p*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2))/x**S(2), x), x, 1/x))
    rubi.add(rule443)

    pattern444 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda m: IntegerQ(m)))
    rule444 = ReplacementRule(pattern444, lambda c, m, d, n, a, p, x : -c**p*Subst(Int(x**(-m + S(-2))*(S(1) + d*x/c)**p*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2)), x), x, 1/x))
    rubi.add(rule444)

    pattern445 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: Not(IntegerQ(p) | PositiveQ(c))))
    rule445 = ReplacementRule(pattern445, lambda c, d, n, a, p, x : (S(1) + d/(c*x))**(-p)*(c + d/x)**p*Int((S(1) + d/(c*x))**p*exp(n*acot(a*x)), x))
    rubi.add(rule445)

    pattern446 = Pattern(Integral(x_**m_*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule446 = ReplacementRule(pattern446, lambda c, m, d, n, a, p, x : -c**p*x**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(S(1) + d*x/c)**p*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2)), x), x, 1/x))
    rubi.add(rule446)

    pattern447 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(a**S(2)*d**S(2) + c**S(2))), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: Not(IntegerQ(p) | PositiveQ(c))))
    rule447 = ReplacementRule(pattern447, lambda c, d, n, a, p, u, x : (S(1) + d/(c*x))**(-p)*(c + d/x)**p*Int(u*(S(1) + d/(c*x))**p*exp(n*acot(a*x)), x))
    rubi.add(rule447)

    pattern448 = Pattern(Integral(exp(WC('n', S(1))*acot(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)))
    rule448 = ReplacementRule(pattern448, lambda c, d, n, a, x : -exp(n*acot(a*x))/(a*c*n))
    rubi.add(rule448)

    pattern449 = Pattern(Integral(exp(WC('n', S(1))*acot(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: Not(OddQ(ImaginaryI*n))))
    rule449 = ReplacementRule(pattern449, lambda c, d, n, a, x : (a*x - n)*exp(n*acot(a*x))/(a*c*sqrt(c + d*x**S(2))*(n**S(2) + S(1))))
    rubi.add(rule449)

    pattern450 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)), CustomConstraint(lambda n, p: NonzeroQ(n**S(2) + S(4)*(p + S(1))**S(2))), CustomConstraint(lambda n, p: Not(IntegerQ(p) & EvenQ(ImaginaryI*n))), CustomConstraint(lambda n, p: Not(Not(IntegerQ(p)) & OddQ(ImaginaryI*n))))
    rule450 = ReplacementRule(pattern450, lambda c, d, n, a, p, x : S(2)*(p + S(1))*(S(2)*p + S(3))*Int((c + d*x**S(2))**(p + S(1))*exp(n*acot(a*x)), x)/(c*(n**S(2) + S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(-S(2)*a*x*(p + S(1)) - n)*exp(n*acot(a*x))/(a*c*(n**S(2) + S(4)*(p + S(1))**S(2))))
    rubi.add(rule450)

    pattern451 = Pattern(Integral(x_*exp(WC('n', S(1))*acot(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: Not(OddQ(ImaginaryI*n))))
    rule451 = ReplacementRule(pattern451, lambda c, d, n, a, x : (-a*n*x + S(-1))*exp(n*acot(a*x))/(a**S(2)*c*sqrt(c + d*x**S(2))*(n**S(2) + S(1))))
    rubi.add(rule451)

    pattern452 = Pattern(Integral(x_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: LessEqual(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-3)/2)), CustomConstraint(lambda n, p: NonzeroQ(n**S(2) + S(4)*(p + S(1))**S(2))), CustomConstraint(lambda n, p: Not(IntegerQ(p) & EvenQ(ImaginaryI*n))), CustomConstraint(lambda n, p: Not(Not(IntegerQ(p)) & OddQ(ImaginaryI*n))))
    rule452 = ReplacementRule(pattern452, lambda c, d, n, a, p, x : n*(S(2)*p + S(3))*Int((c + d*x**S(2))**(p + S(1))*exp(n*acot(a*x)), x)/(a*c*(n**S(2) + S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(-a*n*x + S(2)*p + S(2))*exp(n*acot(a*x))/(a**S(2)*c*(n**S(2) + S(4)*(p + S(1))**S(2))))
    rubi.add(rule452)

    pattern453 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n, p: ZeroQ(n**S(2) - S(2)*p + S(-2))), CustomConstraint(lambda n: NonzeroQ(n**S(2) + S(1))))
    rule453 = ReplacementRule(pattern453, lambda c, d, n, a, p, x : (c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*acot(a*x))/(a**S(3)*c*n**S(2)*(n**S(2) + S(1))))
    rubi.add(rule453)

    pattern454 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: LessEqual(p, S(-1))), CustomConstraint(lambda n, p: NonzeroQ(n**S(2) - S(2)*p + S(-2))), CustomConstraint(lambda n, p: NonzeroQ(n**S(2) + S(4)*(p + S(1))**S(2))), CustomConstraint(lambda n, p: Not(IntegerQ(p) & EvenQ(ImaginaryI*n))), CustomConstraint(lambda n, p: Not(Not(IntegerQ(p)) & OddQ(ImaginaryI*n))))
    rule454 = ReplacementRule(pattern454, lambda c, d, n, a, p, x : (n**S(2) - S(2)*p + S(-2))*Int((c + d*x**S(2))**(p + S(1))*exp(n*acot(a*x)), x)/(a**S(2)*c*(n**S(2) + S(4)*(p + S(1))**S(2))) + (c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*acot(a*x))/(a**S(3)*c*(n**S(2) + S(4)*(p + S(1))**S(2))))
    rubi.add(rule454)

    pattern455 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p, m: LessEqual(S(3), m, -S(2)*p + S(-2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule455 = ReplacementRule(pattern455, lambda c, m, d, n, a, p, x : -a**(-m + S(-1))*c**p*Subst(Int(exp(n*x)*cos(x)**(-S(2)*p + S(-2))*cot(x)**(m + S(2)*p + S(2)), x), x, acot(a*x)))
    rubi.add(rule455)

    pattern456 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p: IntegerQ(p)))
    rule456 = ReplacementRule(pattern456, lambda c, d, n, a, p, u, x : d**p*Int(u*x**(S(2)*p)*(S(1) + S(1)/(a**S(2)*x**S(2)))**p*exp(n*acot(a*x)), x))
    rubi.add(rule456)

    pattern457 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*c + d)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p: Not(IntegerQ(p))))
    rule457 = ReplacementRule(pattern457, lambda c, d, n, a, p, u, x : x**(-S(2)*p)*(S(1) + S(1)/(a**S(2)*x**S(2)))**(-p)*(c + d*x**S(2))**p*Int(u*x**(S(2)*p)*(S(1) + S(1)/(a**S(2)*x**S(2)))**p*exp(n*acot(a*x)), x))
    rubi.add(rule457)

    pattern458 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda n, p: IntegersQ(S(2)*p, ImaginaryI*n/S(2) + p)))
    rule458 = ReplacementRule(pattern458, lambda c, d, n, a, p, u, x : c**p*(ImaginaryI*a)**(-S(2)*p)*Int(u*x**(-S(2)*p)*(ImaginaryI*a*x + S(-1))**(-ImaginaryI*n/S(2) + p)*(ImaginaryI*a*x + S(1))**(ImaginaryI*n/S(2) + p), x))
    rubi.add(rule458)

    pattern459 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda n, p: Not(IntegersQ(S(2)*p, ImaginaryI*n/S(2) + p))))
    rule459 = ReplacementRule(pattern459, lambda c, d, n, a, p, x : -c**p*Subst(Int((-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + p)/x**S(2), x), x, 1/x))
    rubi.add(rule459)

    pattern460 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda n, p: Not(IntegersQ(S(2)*p, ImaginaryI*n/S(2) + p))), CustomConstraint(lambda m: IntegerQ(m)))
    rule460 = ReplacementRule(pattern460, lambda c, m, d, n, a, p, x : -c**p*Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + p), x), x, 1/x))
    rubi.add(rule460)

    pattern461 = Pattern(Integral(x_**m_*(c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: IntegerQ(p) | PositiveQ(c)), CustomConstraint(lambda n, p: Not(IntegersQ(S(2)*p, ImaginaryI*n/S(2) + p))), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule461 = ReplacementRule(pattern461, lambda c, m, d, n, a, p, x : -c**p*x**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(-ImaginaryI*x/a + S(1))**(ImaginaryI*n/S(2) + p)*(ImaginaryI*x/a + S(1))**(-ImaginaryI*n/S(2) + p), x), x, 1/x))
    rubi.add(rule461)

    pattern462 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda a, c, d: ZeroQ(-a**S(2)*d + c)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda p, c: Not(IntegerQ(p) | PositiveQ(c))))
    rule462 = ReplacementRule(pattern462, lambda c, d, n, a, p, u, x : (S(1) + S(1)/(a**S(2)*x**S(2)))**(-p)*(c + d/x**S(2))**p*Int(u*(S(1) + S(1)/(a**S(2)*x**S(2)))**p*exp(n*acot(a*x)), x))
    rubi.add(rule462)

    pattern463 = Pattern(Integral(WC('u', S(1))*exp(n_*acot((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(ImaginaryI*n/S(2))))
    rule463 = ReplacementRule(pattern463, lambda b, c, n, a, u, x : (S(-1))**(ImaginaryI*n/S(2))*Int(u*exp(-n*ArcTan(c*(a + b*x))), x))
    rubi.add(rule463)

    pattern464 = Pattern(Integral(exp(WC('n', S(1))*acot((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))))
    rule464 = ReplacementRule(pattern464, lambda b, c, n, a, x : (ImaginaryI*c*(a + b*x))**(ImaginaryI*n/S(2))*(S(1) + S(1)/(ImaginaryI*c*(a + b*x)))**(ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(-ImaginaryI*n/S(2))*Int((ImaginaryI*a*c + ImaginaryI*b*c*x + S(-1))**(-ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(ImaginaryI*n/S(2)), x))
    rubi.add(rule464)

    pattern465 = Pattern(Integral(x_**m_*exp(n_*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda n: RationalQ(ImaginaryI*n)), CustomConstraint(lambda n: Less(S(-1), ImaginaryI*n, S(1))))
    rule465 = ReplacementRule(pattern465, lambda b, c, m, n, a, x : S(4)*ImaginaryI**(-m)*b**(-m + S(-1))*c**(-m + S(-1))*Subst(Int(x**(S(2)/(ImaginaryI*n))*(x**(S(2)/(ImaginaryI*n)) + S(-1))**(-m + S(-2))*(ImaginaryI*a*c + x**(S(2)/(ImaginaryI*n))*(-ImaginaryI*a*c + S(1)) + S(1))**m, x), x, (S(1) - S(1)/(ImaginaryI*c*(a + b*x)))**(-ImaginaryI*n/S(2))*(S(1) + S(1)/(ImaginaryI*c*(a + b*x)))**(ImaginaryI*n/S(2)))/n)
    rubi.add(rule465)

    pattern466 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*exp(WC('n', S(1))*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))))
    rule466 = ReplacementRule(pattern466, lambda b, e, c, m, d, n, a, x : (ImaginaryI*c*(a + b*x))**(ImaginaryI*n/S(2))*(S(1) + S(1)/(ImaginaryI*c*(a + b*x)))**(ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(-ImaginaryI*n/S(2))*Int((d + e*x)**m*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(-1))**(-ImaginaryI*n/S(2))*(ImaginaryI*a*c + ImaginaryI*b*c*x + S(1))**(ImaginaryI*n/S(2)), x))
    rubi.add(rule466)

    pattern467 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acot(a_ + x_*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda b, a, e, d: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda b, a, e, c: ZeroQ(b**S(2)*c - e*(a**S(2) + S(1)))), CustomConstraint(lambda a, p, c: IntegerQ(p) | PositiveQ(c/(a**S(2) + S(1)))))
    rule467 = ReplacementRule(pattern467, lambda b, e, c, d, n, a, p, u, x : (c/(a**S(2) + S(1)))**p*((ImaginaryI*a + ImaginaryI*b*x + S(1))/(ImaginaryI*a + ImaginaryI*b*x))**(ImaginaryI*n/S(2))*((ImaginaryI*a + ImaginaryI*b*x)/(ImaginaryI*a + ImaginaryI*b*x + S(1)))**(ImaginaryI*n/S(2))*(-ImaginaryI*a - ImaginaryI*b*x + S(1))**(ImaginaryI*n/S(2))*(ImaginaryI*a + ImaginaryI*b*x + S(-1))**(-ImaginaryI*n/S(2))*Int(u*(-ImaginaryI*a - ImaginaryI*b*x + S(1))**(-ImaginaryI*n/S(2) + p)*(ImaginaryI*a + ImaginaryI*b*x + S(1))**(ImaginaryI*n/S(2) + p), x))
    rubi.add(rule467)

    pattern468 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acot(a_ + x_*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: Not(IntegerQ(ImaginaryI*n/S(2)))), CustomConstraint(lambda b, a, e, d: ZeroQ(-S(2)*a*e + b*d)), CustomConstraint(lambda b, a, e, c: ZeroQ(b**S(2)*c - e*(a**S(2) + S(1)))), CustomConstraint(lambda a, p, c: Not(IntegerQ(p) | PositiveQ(c/(a**S(2) + S(1))))))
    rule468 = ReplacementRule(pattern468, lambda b, e, c, d, n, a, p, u, x : (c + d*x + e*x**S(2))**p*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**(-p)*Int(u*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**p*exp(n*acot(a*x)), x))
    rubi.add(rule468)

    pattern469 = Pattern(Integral(WC('u', S(1))*exp(WC('n', S(1))*acot(WC('c', S(1))/(x_*WC('b', S(1)) + WC('a', S(0))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule469 = ReplacementRule(pattern469, lambda b, c, n, a, u, x : Int(u*exp(n*ArcTan(a/c + b*x/c)), x))
    rubi.add(rule469)

    pattern470 = Pattern(Integral((ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule470 = ReplacementRule(pattern470, lambda b, c, d, n, a, x : Subst(Int((a + b*ArcTan(x))**n, x), x, c + d*x)/d)
    rubi.add(rule470)

    pattern471 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule471 = ReplacementRule(pattern471, lambda b, c, d, n, a, x : Subst(Int((a + b*acot(x))**n, x), x, c + d*x)/d)
    rubi.add(rule471)

    pattern472 = Pattern(Integral((ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule472 = ReplacementRule(pattern472, lambda b, c, d, n, a, x : Int((a + b*ArcTan(c + d*x))**n, x))
    rubi.add(rule472)

    pattern473 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule473 = ReplacementRule(pattern473, lambda b, c, d, n, a, x : Int((a + b*acot(c + d*x))**n, x))
    rubi.add(rule473)

    pattern474 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule474 = ReplacementRule(pattern474, lambda b, e, c, m, d, n, a, f, x : Subst(Int((a + b*ArcTan(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule474)

    pattern475 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule475 = ReplacementRule(pattern475, lambda b, e, c, m, d, n, a, f, x : Subst(Int((a + b*acot(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule475)

    pattern476 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**m_*(ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule476 = ReplacementRule(pattern476, lambda b, e, c, m, d, n, a, f, x : Int((a + b*ArcTan(c + d*x))**n*(e + f*x)**m, x))
    rubi.add(rule476)

    pattern477 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**m_*(WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(PositiveIntegerQ(n))))
    rule477 = ReplacementRule(pattern477, lambda b, e, c, m, d, n, a, f, x : Int((a + b*acot(c + d*x))**n*(e + f*x)**m, x))
    rubi.add(rule477)

    pattern478 = Pattern(Integral((ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(-S(2)*A*c*d + B*(c**S(2) + S(1)))), CustomConstraint(lambda C, B, c, d: ZeroQ(-B*d + S(2)*C*c)))
    rule478 = ReplacementRule(pattern478, lambda b, C, B, c, d, n, a, p, x, A : Subst(Int((a + b*ArcTan(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule478)

    pattern479 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(-S(2)*A*c*d + B*(c**S(2) + S(1)))), CustomConstraint(lambda C, B, c, d: ZeroQ(-B*d + S(2)*C*c)))
    rule479 = ReplacementRule(pattern479, lambda b, C, B, c, d, n, a, p, x, A : Subst(Int((a + b*acot(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p, x), x, c + d*x)/d)
    rubi.add(rule479)

    pattern480 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(-S(2)*A*c*d + B*(c**S(2) + S(1)))), CustomConstraint(lambda C, B, c, d: ZeroQ(-B*d + S(2)*C*c)))
    rule480 = ReplacementRule(pattern480, lambda b, C, e, B, c, m, d, n, a, f, p, x, A : Subst(Int((a + b*ArcTan(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule480)

    pattern481 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda B, c, A, d: ZeroQ(-S(2)*A*c*d + B*(c**S(2) + S(1)))), CustomConstraint(lambda C, B, c, d: ZeroQ(-B*d + S(2)*C*c)))
    rule481 = ReplacementRule(pattern481, lambda b, C, e, B, c, m, d, n, a, f, p, x, A : Subst(Int((a + b*acot(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x)/d)
    rubi.add(rule481)

    pattern482 = Pattern(Integral(ArcTan(a_ + x_*WC('b', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: RationalQ(n)))
    rule482 = ReplacementRule(pattern482, lambda b, c, d, n, a, x : ImaginaryI*Int(log(-ImaginaryI*a - ImaginaryI*b*x + S(1))/(c + d*x**n), x)/S(2) - ImaginaryI*Int(log(ImaginaryI*a + ImaginaryI*b*x + S(1))/(c + d*x**n), x)/S(2))
    rubi.add(rule482)

    pattern483 = Pattern(Integral(acot(a_ + x_*WC('b', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: RationalQ(n)))
    rule483 = ReplacementRule(pattern483, lambda b, c, d, n, a, x : ImaginaryI*Int(log((-ImaginaryI + a + b*x)/(a + b*x))/(c + d*x**n), x)/S(2) - ImaginaryI*Int(log((ImaginaryI + a + b*x)/(a + b*x))/(c + d*x**n), x)/S(2))
    rubi.add(rule483)

    pattern484 = Pattern(Integral(ArcTan(a_ + x_*WC('b', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(RationalQ(n))))
    rule484 = ReplacementRule(pattern484, lambda b, c, d, n, a, x : Int(ArcTan(a + b*x)/(c + d*x**n), x))
    rubi.add(rule484)

    pattern485 = Pattern(Integral(acot(a_ + x_*WC('b', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(RationalQ(n))))
    rule485 = ReplacementRule(pattern485, lambda b, c, d, n, a, x : Int(acot(a + b*x)/(c + d*x**n), x))
    rubi.add(rule485)

    pattern486 = Pattern(Integral(ArcTan(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule486 = ReplacementRule(pattern486, lambda b, a, x, n : -b*n*Int(x**n/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x) + x*ArcTan(a + b*x**n))
    rubi.add(rule486)

    pattern487 = Pattern(Integral(acot(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule487 = ReplacementRule(pattern487, lambda b, a, x, n : b*n*Int(x**n/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x) + x*acot(a + b*x**n))
    rubi.add(rule487)

    pattern488 = Pattern(Integral(ArcTan(x_**n_*WC('b', S(1)) + WC('a', S(0)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule488 = ReplacementRule(pattern488, lambda b, a, x, n : ImaginaryI*Int(log(-ImaginaryI*a - ImaginaryI*b*x**n + S(1))/x, x)/S(2) - ImaginaryI*Int(log(ImaginaryI*a + ImaginaryI*b*x**n + S(1))/x, x)/S(2))
    rubi.add(rule488)

    pattern489 = Pattern(Integral(acot(x_**n_*WC('b', S(1)) + WC('a', S(0)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule489 = ReplacementRule(pattern489, lambda b, a, x, n : ImaginaryI*Int(log(-ImaginaryI/(a + b*x**n) + S(1))/x, x)/S(2) - ImaginaryI*Int(log(ImaginaryI/(a + b*x**n) + S(1))/x, x)/S(2))
    rubi.add(rule489)

    pattern490 = Pattern(Integral(x_**WC('m', S(1))*ArcTan(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda m: Unequal(m + S(1), S(0))), CustomConstraint(lambda n, m: Unequal(m + S(1), n)))
    rule490 = ReplacementRule(pattern490, lambda b, m, n, a, x : -b*n*Int(x**(m + n)/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x)/(m + S(1)) + x**(m + S(1))*ArcTan(a + b*x**n)/(m + S(1)))
    rubi.add(rule490)

    pattern491 = Pattern(Integral(x_**WC('m', S(1))*acot(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda m: Unequal(m + S(1), S(0))), CustomConstraint(lambda n, m: Unequal(m + S(1), n)))
    rule491 = ReplacementRule(pattern491, lambda b, m, n, a, x : b*n*Int(x**(m + n)/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x)/(m + S(1)) + x**(m + S(1))*acot(a + b*x**n)/(m + S(1)))
    rubi.add(rule491)

    pattern492 = Pattern(Integral(ArcTan(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule492 = ReplacementRule(pattern492, lambda b, c, d, a, f, x : ImaginaryI*Int(log(-ImaginaryI*a - ImaginaryI*b*f**(c + d*x) + S(1)), x)/S(2) - ImaginaryI*Int(log(ImaginaryI*a + ImaginaryI*b*f**(c + d*x) + S(1)), x)/S(2))
    rubi.add(rule492)

    pattern493 = Pattern(Integral(acot(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule493 = ReplacementRule(pattern493, lambda b, c, d, a, f, x : ImaginaryI*Int(log(-ImaginaryI/(a + b*f**(c + d*x)) + S(1)), x)/S(2) - ImaginaryI*Int(log(ImaginaryI/(a + b*f**(c + d*x)) + S(1)), x)/S(2))
    rubi.add(rule493)

    pattern494 = Pattern(Integral(x_**WC('m', S(1))*ArcTan(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule494 = ReplacementRule(pattern494, lambda b, c, m, d, a, f, x : ImaginaryI*Int(x**m*log(-ImaginaryI*a - ImaginaryI*b*f**(c + d*x) + S(1)), x)/S(2) - ImaginaryI*Int(x**m*log(ImaginaryI*a + ImaginaryI*b*f**(c + d*x) + S(1)), x)/S(2))
    rubi.add(rule494)

    pattern495 = Pattern(Integral(x_**WC('m', S(1))*acot(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule495 = ReplacementRule(pattern495, lambda b, c, m, d, a, f, x : ImaginaryI*Int(x**m*log(-ImaginaryI/(a + b*f**(c + d*x)) + S(1)), x)/S(2) - ImaginaryI*Int(x**m*log(ImaginaryI/(a + b*f**(c + d*x)) + S(1)), x)/S(2))
    rubi.add(rule495)

    pattern496 = Pattern(Integral(ArcTan(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule496 = ReplacementRule(pattern496, lambda b, c, m, n, a, u, x : Int(u*acot(a/c + b*x**n/c)**m, x))
    rubi.add(rule496)

    pattern497 = Pattern(Integral(WC('u', S(1))*acot(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule497 = ReplacementRule(pattern497, lambda b, c, m, n, a, u, x : Int(u*ArcTan(a/c + b*x**n/c)**m, x))
    rubi.add(rule497)

    pattern498 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0)))*ArcTan(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))))
    rule498 = ReplacementRule(pattern498, lambda b, a, c, x : log(ArcTan(c*x/sqrt(a + b*x**S(2))))/c)
    rubi.add(rule498)

    pattern499 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0)))*acot(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))))
    rule499 = ReplacementRule(pattern499, lambda b, a, c, x : -log(acot(c*x/sqrt(a + b*x**S(2))))/c)
    rubi.add(rule499)

    pattern500 = Pattern(Integral(ArcTan(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule500 = ReplacementRule(pattern500, lambda b, c, m, a, x : ArcTan(c*x/sqrt(a + b*x**S(2)))**(m + S(1))/(c*(m + S(1))))
    rubi.add(rule500)

    pattern501 = Pattern(Integral(acot(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule501 = ReplacementRule(pattern501, lambda b, c, m, a, x : -acot(c*x/sqrt(a + b*x**S(2)))**(m + S(1))/(c*(m + S(1))))
    rubi.add(rule501)

    pattern502 = Pattern(Integral(ArcTan(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))), CustomConstraint(lambda b, a, e, d: ZeroQ(-a*e + b*d)))
    rule502 = ReplacementRule(pattern502, lambda b, e, c, m, d, a, x : sqrt(a + b*x**S(2))*Int(ArcTan(c*x/sqrt(a + b*x**S(2)))**m/sqrt(a + b*x**S(2)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule502)

    pattern503 = Pattern(Integral(acot(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, c: ZeroQ(b + c**S(2))), CustomConstraint(lambda b, a, e, d: ZeroQ(-a*e + b*d)))
    rule503 = ReplacementRule(pattern503, lambda b, e, c, m, d, a, x : sqrt(a + b*x**S(2))*Int(acot(c*x/sqrt(a + b*x**S(2)))**m/sqrt(a + b*x**S(2)), x)/sqrt(d + e*x**S(2)))
    rubi.add(rule503)

    pattern504 = Pattern(Integral(ArcTan(v_ + sqrt(w_)*WC('s', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda s: ZeroQ(s**S(2) + S(-1))), CustomConstraint(lambda v, w: ZeroQ(-v**S(2) + w + S(-1))))
    rule504 = ReplacementRule(pattern504, lambda w, u, v, s, x : Pi*s*Int(u, x)/S(4) + Int(u*ArcTan(v), x)/S(2))
    rubi.add(rule504)

    pattern505 = Pattern(Integral(WC('u', S(1))*acot(v_ + sqrt(w_)*WC('s', S(1))), x_), CustomConstraint(lambda s: ZeroQ(s**S(2) + S(-1))), CustomConstraint(lambda v, w: ZeroQ(-v**S(2) + w + S(-1))))
    rule505 = ReplacementRule(pattern505, lambda w, u, v, s, x : Pi*s*Int(u, x)/S(4) - Int(u*ArcTan(v), x)/S(2))
    rubi.add(rule505)

    pattern506 = Pattern(Integral(u_*v_**WC('n', S(1)), x_), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Less(n, S(0))), CustomConstraint(lambda v, x: NegQ(Discriminant(v, x))), CustomConstraint(lambda u, x: MatchQ(u, Condition(Optional(Pattern(r, Blank))*Pattern(f, Blank)**Pattern(w, Blank)))), CustomConstraint(lambda v, tmp, u, n, x, ArcTan: Not(FalseQ(tmp)) & SameQ(Head(tmp), ArcTan) & ZeroQ(D(v, x)**S(2) + Discriminant(v, x)*Part(tmp, S(1))**S(2))))
    def With506(n, u, v, x):
        tmp = InverseFunctionOfLinear(u, x)
        return (-Discriminant(v, x)/(S(4)*Coefficient(v, x, S(2))))**n*Subst(Int(SimplifyIntegrand(SubstForInverseFunction(u, tmp, x)*sec(x)**(S(2)*n + S(2)), x), x), x, tmp)/Coefficient(Part(tmp, S(1)), x, S(1))
    rule506 = ReplacementRule(pattern506, lambda n, u, v, x : With506(n, u, v, x))
    rubi.add(rule506)

    pattern507 = Pattern(Integral(u_*v_**WC('n', S(1)), x_), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Less(n, S(0))), CustomConstraint(lambda v, x: NegQ(Discriminant(v, x))), CustomConstraint(lambda u, x: MatchQ(u, Condition(Optional(Pattern(r, Blank))*Pattern(f, Blank)**Pattern(w, Blank)))), CustomConstraint(lambda v, tmp, u, ArcCot, n, x: Not(FalseQ(tmp)) & SameQ(Head(tmp), ArcCot) & ZeroQ(D(v, x)**S(2) + Discriminant(v, x)*Part(tmp, S(1))**S(2))))
    def With507(n, u, v, x):
        tmp = InverseFunctionOfLinear(u, x)
        return -(-Discriminant(v, x)/(S(4)*Coefficient(v, x, S(2))))**n*Subst(Int(SimplifyIntegrand(SubstForInverseFunction(u, tmp, x)*csc(x)**(S(2)*n + S(2)), x), x), x, tmp)/Coefficient(Part(tmp, S(1)), x, S(1))
    rule507 = ReplacementRule(pattern507, lambda n, u, v, x : With507(n, u, v, x))
    rubi.add(rule507)

    pattern508 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule508 = ReplacementRule(pattern508, lambda b, c, d, a, x : -ImaginaryI*b*Int(x/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*ArcTan(c + d*tan(a + b*x)))
    rubi.add(rule508)

    pattern509 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule509 = ReplacementRule(pattern509, lambda b, c, d, a, x : ImaginaryI*b*Int(x/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*acot(c + d*tan(a + b*x)))
    rubi.add(rule509)

    pattern510 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule510 = ReplacementRule(pattern510, lambda b, c, d, a, x : -ImaginaryI*b*Int(x/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*ArcTan(c + d*cot(a + b*x)))
    rubi.add(rule510)

    pattern511 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule511 = ReplacementRule(pattern511, lambda b, c, d, a, x : ImaginaryI*b*Int(x/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x) + x*acot(c + d*cot(a + b*x)))
    rubi.add(rule511)

    pattern512 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule512 = ReplacementRule(pattern512, lambda b, c, d, a, x : b*(-ImaginaryI*c - d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c + d + (-ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) - b*(ImaginaryI*c + d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c - d + (ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*ArcTan(c + d*tan(a + b*x)))
    rubi.add(rule512)

    pattern513 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule513 = ReplacementRule(pattern513, lambda b, c, d, a, x : -b*(-ImaginaryI*c - d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c + d + (-ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + b*(ImaginaryI*c + d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c - d + (ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*acot(c + d*tan(a + b*x)))
    rubi.add(rule513)

    pattern514 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule514 = ReplacementRule(pattern514, lambda b, c, d, a, x : -b*(-ImaginaryI*c + d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c - d - (-ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + b*(ImaginaryI*c - d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c + d - (ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*ArcTan(c + d*cot(a + b*x)))
    rubi.add(rule514)

    pattern515 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule515 = ReplacementRule(pattern515, lambda b, c, d, a, x : b*(-ImaginaryI*c + d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c - d - (-ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) - b*(ImaginaryI*c - d + S(1))*Int(x*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c + d - (ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x) + x*acot(c + d*cot(a + b*x)))
    rubi.add(rule515)

    pattern516 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule516 = ReplacementRule(pattern516, lambda b, e, c, m, d, a, f, x : -ImaginaryI*b*Int((e + f*x)**(m + S(1))/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule516)

    pattern517 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule517 = ReplacementRule(pattern517, lambda b, e, c, m, d, a, f, x : ImaginaryI*b*Int((e + f*x)**(m + S(1))/(ImaginaryI*d + c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acot(c + d*tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule517)

    pattern518 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule518 = ReplacementRule(pattern518, lambda b, e, c, m, d, a, f, x : -ImaginaryI*b*Int((e + f*x)**(m + S(1))/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule518)

    pattern519 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule519 = ReplacementRule(pattern519, lambda b, e, c, m, d, a, f, x : ImaginaryI*b*Int((e + f*x)**(m + S(1))/(-ImaginaryI*d - c*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + c), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acot(c + d*cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule519)

    pattern520 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule520 = ReplacementRule(pattern520, lambda b, e, c, m, d, a, f, x : b*(-ImaginaryI*c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c + d + (-ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) - b*(ImaginaryI*c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c - d + (ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule520)

    pattern521 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((ImaginaryI*d + c)**S(2) + S(1))))
    rule521 = ReplacementRule(pattern521, lambda b, e, c, m, d, a, f, x : -b*(-ImaginaryI*c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c + d + (-ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + b*(ImaginaryI*c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c - d + (ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acot(c + d*tan(a + b*x))/(f*(m + S(1))))
    rubi.add(rule521)

    pattern522 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule522 = ReplacementRule(pattern522, lambda b, e, c, m, d, a, f, x : -b*(-ImaginaryI*c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c - d - (-ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + b*(ImaginaryI*c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c + d - (ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule522)

    pattern523 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((-ImaginaryI*d + c)**S(2) + S(1))))
    rule523 = ReplacementRule(pattern523, lambda b, e, c, m, d, a, f, x : b*(-ImaginaryI*c + d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(-ImaginaryI*c - d - (-ImaginaryI*c + d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) - b*(ImaginaryI*c - d + S(1))*Int((e + f*x)**(m + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x)/(ImaginaryI*c + d - (ImaginaryI*c - d + S(1))*exp(S(2)*ImaginaryI*a + S(2)*ImaginaryI*b*x) + S(1)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acot(c + d*cot(a + b*x))/(f*(m + S(1))))
    rubi.add(rule523)

    pattern524 = Pattern(Integral(ArcTan(tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule524 = ReplacementRule(pattern524, lambda b, a, x : -b*Int(x*sech(S(2)*a + S(2)*b*x), x) + x*ArcTan(tanh(a + b*x)))
    rubi.add(rule524)

    pattern525 = Pattern(Integral(acot(tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule525 = ReplacementRule(pattern525, lambda b, a, x : b*Int(x*sech(S(2)*a + S(2)*b*x), x) + x*acot(tanh(a + b*x)))
    rubi.add(rule525)

    pattern526 = Pattern(Integral(ArcTan(coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule526 = ReplacementRule(pattern526, lambda b, a, x : b*Int(x*sech(S(2)*a + S(2)*b*x), x) + x*ArcTan(coth(a + b*x)))
    rubi.add(rule526)

    pattern527 = Pattern(Integral(acot(coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule527 = ReplacementRule(pattern527, lambda b, a, x : -b*Int(x*sech(S(2)*a + S(2)*b*x), x) + x*acot(coth(a + b*x)))
    rubi.add(rule527)

    pattern528 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule528 = ReplacementRule(pattern528, lambda b, e, m, a, f, x : -b*Int((e + f*x)**(m + S(1))*sech(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule528)

    pattern529 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule529 = ReplacementRule(pattern529, lambda b, e, m, a, f, x : b*Int((e + f*x)**(m + S(1))*sech(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acot(tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule529)

    pattern530 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule530 = ReplacementRule(pattern530, lambda b, e, m, a, f, x : b*Int((e + f*x)**(m + S(1))*sech(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule530)

    pattern531 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule531 = ReplacementRule(pattern531, lambda b, e, m, a, f, x : -b*Int((e + f*x)**(m + S(1))*sech(S(2)*a + S(2)*b*x), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acot(coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule531)

    pattern532 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(1))))
    rule532 = ReplacementRule(pattern532, lambda b, c, d, a, x : -b*Int(x/(c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*ArcTan(c + d*tanh(a + b*x)))
    rubi.add(rule532)

    pattern533 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(1))))
    rule533 = ReplacementRule(pattern533, lambda b, c, d, a, x : b*Int(x/(c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*acot(c + d*tanh(a + b*x)))
    rubi.add(rule533)

    pattern534 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(1))))
    rule534 = ReplacementRule(pattern534, lambda b, c, d, a, x : -b*Int(x/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*ArcTan(c + d*coth(a + b*x)))
    rubi.add(rule534)

    pattern535 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(1))))
    rule535 = ReplacementRule(pattern535, lambda b, c, d, a, x : b*Int(x/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x) + x*acot(c + d*coth(a + b*x)))
    rubi.add(rule535)

    pattern536 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(1))))
    rule536 = ReplacementRule(pattern536, lambda b, c, d, a, x : ImaginaryI*b*(ImaginaryI - c - d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d + (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x) - ImaginaryI*b*(ImaginaryI + c + d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d + (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x) + x*ArcTan(c + d*tanh(a + b*x)))
    rubi.add(rule536)

    pattern537 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(1))))
    rule537 = ReplacementRule(pattern537, lambda b, c, d, a, x : -ImaginaryI*b*(ImaginaryI - c - d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d + (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x) + ImaginaryI*b*(ImaginaryI + c + d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d + (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x) + x*acot(c + d*tanh(a + b*x)))
    rubi.add(rule537)

    pattern538 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(1))))
    rule538 = ReplacementRule(pattern538, lambda b, c, d, a, x : -ImaginaryI*b*(ImaginaryI - c - d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d - (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x) + ImaginaryI*b*(ImaginaryI + c + d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d - (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x) + x*ArcTan(c + d*coth(a + b*x)))
    rubi.add(rule538)

    pattern539 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(1))))
    rule539 = ReplacementRule(pattern539, lambda b, c, d, a, x : ImaginaryI*b*(ImaginaryI - c - d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d - (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x) - ImaginaryI*b*(ImaginaryI + c + d)*Int(x*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d - (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x) + x*acot(c + d*coth(a + b*x)))
    rubi.add(rule539)

    pattern540 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(1))))
    rule540 = ReplacementRule(pattern540, lambda b, e, c, m, d, a, f, x : -b*Int((e + f*x)**(m + S(1))/(c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule540)

    pattern541 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(1))))
    rule541 = ReplacementRule(pattern541, lambda b, e, c, m, d, a, f, x : b*Int((e + f*x)**(m + S(1))/(c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acot(c + d*tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule541)

    pattern542 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(1))))
    rule542 = ReplacementRule(pattern542, lambda b, e, c, m, d, a, f, x : -b*Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule542)

    pattern543 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: ZeroQ((c - d)**S(2) + S(1))))
    rule543 = ReplacementRule(pattern543, lambda b, e, c, m, d, a, f, x : b*Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acot(c + d*coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule543)

    pattern544 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(1))))
    rule544 = ReplacementRule(pattern544, lambda b, e, c, m, d, a, f, x : ImaginaryI*b*(ImaginaryI - c - d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d + (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) - ImaginaryI*b*(ImaginaryI + c + d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d + (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule544)

    pattern545 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(1))))
    rule545 = ReplacementRule(pattern545, lambda b, e, c, m, d, a, f, x : -ImaginaryI*b*(ImaginaryI - c - d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d + (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + ImaginaryI*b*(ImaginaryI + c + d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d + (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acot(c + d*tanh(a + b*x))/(f*(m + S(1))))
    rubi.add(rule545)

    pattern546 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(1))))
    rule546 = ReplacementRule(pattern546, lambda b, e, c, m, d, a, f, x : -ImaginaryI*b*(ImaginaryI - c - d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d - (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + ImaginaryI*b*(ImaginaryI + c + d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d - (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*ArcTan(c + d*coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule546)

    pattern547 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda c, d: NonzeroQ((c - d)**S(2) + S(1))))
    rule547 = ReplacementRule(pattern547, lambda b, e, c, m, d, a, f, x : ImaginaryI*b*(ImaginaryI - c - d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI - c + d - (ImaginaryI - c - d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) - ImaginaryI*b*(ImaginaryI + c + d)*Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(ImaginaryI + c - d - (ImaginaryI + c + d)*exp(S(2)*a + S(2)*b*x)), x)/(f*(m + S(1))) + (e + f*x)**(m + S(1))*acot(c + d*coth(a + b*x))/(f*(m + S(1))))
    rubi.add(rule547)

    pattern548 = Pattern(Integral(ArcTan(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)))
    rule548 = ReplacementRule(pattern548, lambda u, x : x*ArcTan(u) - Int(SimplifyIntegrand(x*D(u, x)/(u**S(2) + S(1)), x), x))
    rubi.add(rule548)

    pattern549 = Pattern(Integral(acot(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)))
    rule549 = ReplacementRule(pattern549, lambda u, x : x*acot(u) + Int(SimplifyIntegrand(x*D(u, x)/(u**S(2) + S(1)), x), x))
    rubi.add(rule549)

    pattern550 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(ArcTan(u_)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, m, d, u, x: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, m, x: FalseQ(PowerVariableExpn(u, m + S(1), x))))
    rule550 = ReplacementRule(pattern550, lambda b, c, m, d, a, u, x : -b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u**S(2) + S(1)), x), x)/(d*(m + S(1))) + (a + b*ArcTan(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule550)

    pattern551 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(u_)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, m, d, u, x: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, m, x: FalseQ(PowerVariableExpn(u, m + S(1), x))))
    rule551 = ReplacementRule(pattern551, lambda b, c, m, d, a, u, x : b*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u**S(2) + S(1)), x), x)/(d*(m + S(1))) + (a + b*acot(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule551)

    pattern552 = Pattern(Integral(v_*(ArcTan(u_)*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda b, a, u, v, x: FalseQ(FunctionOfLinear(v*(a + b*ArcTan(u)), x))), CustomConstraint(lambda w, u, x, a, b: InverseFunctionFreeQ(w, x)))
    def With552(b, a, u, v, x):
        w = IntHide(v, x)
        return -b*Int(SimplifyIntegrand(w*D(u, x)/(u**S(2) + S(1)), x), x) + Dist(a + b*ArcTan(u), w, x)
    rule552 = ReplacementRule(pattern552, lambda b, a, u, v, x : With552(b, a, u, v, x))
    rubi.add(rule552)

    pattern553 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acot(u_)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda b, a, u, v, x: FalseQ(FunctionOfLinear(v*(a + b*acot(u)), x))), CustomConstraint(lambda w, u, x, a, b: InverseFunctionFreeQ(w, x)))
    def With553(b, a, u, v, x):
        w = IntHide(v, x)
        return b*Int(SimplifyIntegrand(w*D(u, x)/(u**S(2) + S(1)), x), x) + Dist(a + b*acot(u), w, x)
    rule553 = ReplacementRule(pattern553, lambda b, a, u, v, x : With553(b, a, u, v, x))
    rubi.add(rule553)

    pattern554 = Pattern(Integral(ArcTan(v_)*log(w_)/(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda x, w: LinearQ(w, x)), CustomConstraint(lambda b, a, v, x: ZeroQ(D(v/(a + b*x), x))), CustomConstraint(lambda b, a, x, w: ZeroQ(D(w/(a + b*x), x))))
    rule554 = ReplacementRule(pattern554, lambda b, w, a, v, x : ImaginaryI*Int(log(w)*log(-ImaginaryI*v + S(1))/(a + b*x), x)/S(2) - ImaginaryI*Int(log(w)*log(ImaginaryI*v + S(1))/(a + b*x), x)/S(2))
    rubi.add(rule554)

    pattern555 = Pattern(Integral(ArcTan(v_)*log(w_), x_), CustomConstraint(lambda v, x: InverseFunctionFreeQ(v, x)), CustomConstraint(lambda x, w: InverseFunctionFreeQ(w, x)))
    rule555 = ReplacementRule(pattern555, lambda v, x, w : x*ArcTan(v)*log(w) - Int(SimplifyIntegrand(x*ArcTan(v)*D(w, x)/w, x), x) - Int(SimplifyIntegrand(x*D(v, x)*log(w)/(v**S(2) + S(1)), x), x))
    rubi.add(rule555)

    pattern556 = Pattern(Integral(log(w_)*acot(v_), x_), CustomConstraint(lambda v, x: InverseFunctionFreeQ(v, x)), CustomConstraint(lambda x, w: InverseFunctionFreeQ(w, x)))
    rule556 = ReplacementRule(pattern556, lambda v, x, w : x*log(w)*acot(v) - Int(SimplifyIntegrand(x*D(w, x)*acot(v)/w, x), x) + Int(SimplifyIntegrand(x*D(v, x)*log(w)/(v**S(2) + S(1)), x), x))
    rubi.add(rule556)

    pattern557 = Pattern(Integral(u_*ArcTan(v_)*log(w_), x_), CustomConstraint(lambda v, x: InverseFunctionFreeQ(v, x)), CustomConstraint(lambda x, w: InverseFunctionFreeQ(w, x)), CustomConstraint(lambda z, w, v, x: InverseFunctionFreeQ(z, x)))
    def With557(u, v, x, w):
        z = IntHide(u, x)
        return Dist(ArcTan(v)*log(w), z, x) - Int(SimplifyIntegrand(z*ArcTan(v)*D(w, x)/w, x), x) - Int(SimplifyIntegrand(z*D(v, x)*log(w)/(v**S(2) + S(1)), x), x)
    rule557 = ReplacementRule(pattern557, lambda u, v, x, w : With557(u, v, x, w))
    rubi.add(rule557)

    pattern558 = Pattern(Integral(u_*log(w_)*acot(v_), x_), CustomConstraint(lambda v, x: InverseFunctionFreeQ(v, x)), CustomConstraint(lambda x, w: InverseFunctionFreeQ(w, x)), CustomConstraint(lambda z, w, v, x: InverseFunctionFreeQ(z, x)))
    def With558(u, v, x, w):
        z = IntHide(u, x)
        return Dist(log(w)*acot(v), z, x) - Int(SimplifyIntegrand(z*D(w, x)*acot(v)/w, x), x) + Int(SimplifyIntegrand(z*D(v, x)*log(w)/(v**S(2) + S(1)), x), x)
    rule558 = ReplacementRule(pattern558, lambda u, v, x, w : With558(u, v, x, w))
    rubi.add(rule558)

    pattern559 = Pattern(Integral(asec(x_*WC('c', S(1))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule559 = ReplacementRule(pattern559, lambda c, x : x*asec(c*x) - Int(S(1)/(x*sqrt(S(1) - S(1)/(c**S(2)*x**S(2)))), x)/c)
    rubi.add(rule559)

    pattern560 = Pattern(Integral(acsc(x_*WC('c', S(1))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule560 = ReplacementRule(pattern560, lambda c, x : x*acsc(c*x) + Int(S(1)/(x*sqrt(S(1) - S(1)/(c**S(2)*x**S(2)))), x)/c)
    rubi.add(rule560)

    pattern561 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule561 = ReplacementRule(pattern561, lambda b, c, n, a, x : Subst(Int((a + b*x)**n*tan(x)*sec(x), x), x, asec(c*x))/c)
    rubi.add(rule561)

    pattern562 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule562 = ReplacementRule(pattern562, lambda b, c, n, a, x : -Subst(Int((a + b*x)**n*cot(x)*csc(x), x), x, acsc(c*x))/c)
    rubi.add(rule562)

    pattern563 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule563 = ReplacementRule(pattern563, lambda b, a, c, x : -Subst(Int((a + b*acos(x/c))/x, x), x, 1/x))
    rubi.add(rule563)

    pattern564 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule564 = ReplacementRule(pattern564, lambda b, a, c, x : -Subst(Int((a + b*asin(x/c))/x, x), x, 1/x))
    rubi.add(rule564)

    pattern565 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule565 = ReplacementRule(pattern565, lambda b, c, m, a, x : -b*Int(x**(m + S(-1))/sqrt(S(1) - S(1)/(c**S(2)*x**S(2))), x)/(c*(m + S(1))) + x**(m + S(1))*(a + b*asec(c*x))/(m + S(1)))
    rubi.add(rule565)

    pattern566 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule566 = ReplacementRule(pattern566, lambda b, c, m, a, x : b*Int(x**(m + S(-1))/sqrt(S(1) - S(1)/(c**S(2)*x**S(2))), x)/(c*(m + S(1))) + x**(m + S(1))*(a + b*acsc(c*x))/(m + S(1)))
    rubi.add(rule566)

    pattern567 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule567 = ReplacementRule(pattern567, lambda b, c, m, n, a, x : c**(-m + S(-1))*Subst(Int((a + b*x)**n*tan(x)*sec(x)**(m + S(1)), x), x, asec(c*x)))
    rubi.add(rule567)

    pattern568 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule568 = ReplacementRule(pattern568, lambda b, c, m, n, a, x : -c**(-m + S(-1))*Subst(Int((a + b*x)**n*cot(x)*csc(x)**(m + S(1)), x), x, acsc(c*x)))
    rubi.add(rule568)

    pattern569 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule569 = ReplacementRule(pattern569, lambda b, c, m, n, a, x : Int(x**m*(a + b*asec(c*x))**n, x))
    rubi.add(rule569)

    pattern570 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule570 = ReplacementRule(pattern570, lambda b, c, m, n, a, x : Int(x**m*(a + b*acsc(c*x))**n, x))
    rubi.add(rule570)

    pattern571 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With571(b, e, c, d, a, p, x):
        u = IntHide((d + e*x**S(2))**p, x)
        return -b*c*x*Int(SimplifyIntegrand(u/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x), x)/sqrt(c**S(2)*x**S(2)) + Dist(a + b*asec(c*x), u, x)
    rule571 = ReplacementRule(pattern571, lambda b, e, c, d, a, p, x : With571(b, e, c, d, a, p, x))
    rubi.add(rule571)

    pattern572 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: PositiveIntegerQ(p) | NegativeIntegerQ(p + S(1)/2)), )
    def With572(b, e, c, d, a, p, x):
        u = IntHide((d + e*x**S(2))**p, x)
        return b*c*x*Int(SimplifyIntegrand(u/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x), x)/sqrt(c**S(2)*x**S(2)) + Dist(a + b*acsc(c*x), u, x)
    rule572 = ReplacementRule(pattern572, lambda b, e, c, d, a, p, x : With572(b, e, c, d, a, p, x))
    rubi.add(rule572)

    pattern573 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)))
    rule573 = ReplacementRule(pattern573, lambda b, e, c, d, n, a, p, x : -Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule573)

    pattern574 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)))
    rule574 = ReplacementRule(pattern574, lambda b, e, c, d, n, a, p, x : -Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule574)

    pattern575 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule575 = ReplacementRule(pattern575, lambda b, e, c, d, n, a, p, x : -sqrt(x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule575)

    pattern576 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule576 = ReplacementRule(pattern576, lambda b, e, c, d, n, a, p, x : -sqrt(x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule576)

    pattern577 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e, d: Not(Negative(d) & PositiveQ(e))))
    rule577 = ReplacementRule(pattern577, lambda b, e, c, d, n, a, p, x : -sqrt(d + e*x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*sqrt(d/x**S(2) + e)))
    rubi.add(rule577)

    pattern578 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e, d: Not(Negative(d) & PositiveQ(e))))
    rule578 = ReplacementRule(pattern578, lambda b, e, c, d, n, a, p, x : -sqrt(d + e*x**S(2))*Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*sqrt(d/x**S(2) + e)))
    rubi.add(rule578)

    pattern579 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule579 = ReplacementRule(pattern579, lambda b, e, c, d, n, a, p, x : Int((a + b*asec(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule579)

    pattern580 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule580 = ReplacementRule(pattern580, lambda b, e, c, d, n, a, p, x : Int((a + b*acsc(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule580)

    pattern581 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule581 = ReplacementRule(pattern581, lambda b, e, c, d, a, p, x : -b*c*x*Int((d + e*x**S(2))**(p + S(1))/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x)/(S(2)*e*sqrt(c**S(2)*x**S(2))*(p + S(1))) + (a + b*asec(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule581)

    pattern582 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule582 = ReplacementRule(pattern582, lambda b, e, c, d, a, p, x : b*c*x*Int((d + e*x**S(2))**(p + S(1))/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x)/(S(2)*e*sqrt(c**S(2)*x**S(2))*(p + S(1))) + (a + b*acsc(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))))
    rubi.add(rule582)

    pattern583 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & Not(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))), )
    def With583(b, e, c, m, d, a, p, x):
        u = IntHide(x**m*(d + e*x**S(2))**p, x)
        return -b*c*x*Int(SimplifyIntegrand(u/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x), x)/sqrt(c**S(2)*x**S(2)) + Dist(a + b*asec(c*x), u, x)
    rule583 = ReplacementRule(pattern583, lambda b, e, c, m, d, a, p, x : With583(b, e, c, m, d, a, p, x))
    rubi.add(rule583)

    pattern584 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, m: (NegativeIntegerQ(m/S(2) + p + S(1)/2) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2))) | (PositiveIntegerQ(p) & Not(NegativeIntegerQ(m/S(2) + S(-1)/2) & Greater(m + S(2)*p + S(3), S(0)))) | (PositiveIntegerQ(m/S(2) + S(1)/2) & Not(NegativeIntegerQ(p) & Greater(m + S(2)*p + S(3), S(0))))), )
    def With584(b, e, c, m, d, a, p, x):
        u = IntHide(x**m*(d + e*x**S(2))**p, x)
        return b*c*x*Int(SimplifyIntegrand(u/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x), x)/sqrt(c**S(2)*x**S(2)) + Dist(a + b*acsc(c*x), u, x)
    rule584 = ReplacementRule(pattern584, lambda b, e, c, m, d, a, p, x : With584(b, e, c, m, d, a, p, x))
    rubi.add(rule584)

    pattern585 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, m: IntegersQ(m, p)))
    rule585 = ReplacementRule(pattern585, lambda b, e, c, m, d, n, a, p, x : -Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule585)

    pattern586 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, m: IntegersQ(m, p)))
    rule586 = ReplacementRule(pattern586, lambda b, e, c, m, d, n, a, p, x : -Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x))
    rubi.add(rule586)

    pattern587 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule587 = ReplacementRule(pattern587, lambda b, e, c, m, d, n, a, p, x : -sqrt(x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule587)

    pattern588 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e: PositiveQ(e)), CustomConstraint(lambda d: Negative(d)))
    rule588 = ReplacementRule(pattern588, lambda b, e, c, m, d, n, a, p, x : -sqrt(x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/x)
    rubi.add(rule588)

    pattern589 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e, d: Not(Negative(d) & PositiveQ(e))))
    rule589 = ReplacementRule(pattern589, lambda b, e, c, m, d, n, a, p, x : -sqrt(d + e*x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*sqrt(d/x**S(2) + e)))
    rubi.add(rule589)

    pattern590 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, c, d: ZeroQ(c**S(2)*d + e)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda e, d: Not(Negative(d) & PositiveQ(e))))
    rule590 = ReplacementRule(pattern590, lambda b, e, c, m, d, n, a, p, x : -sqrt(d + e*x**S(2))*Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, 1/x)/(x*sqrt(d/x**S(2) + e)))
    rubi.add(rule590)

    pattern591 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule591 = ReplacementRule(pattern591, lambda b, e, c, m, d, n, a, p, x : Int(x**m*(a + b*asec(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule591)

    pattern592 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule592 = ReplacementRule(pattern592, lambda b, e, c, m, d, n, a, p, x : Int(x**m*(a + b*acsc(c*x))**n*(d + e*x**S(2))**p, x))
    rubi.add(rule592)

    pattern593 = Pattern(Integral(asec(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule593 = ReplacementRule(pattern593, lambda b, a, x : -Int(S(1)/(sqrt(S(1) - S(1)/(a + b*x)**S(2))*(a + b*x)), x) + (a + b*x)*asec(a + b*x)/b)
    rubi.add(rule593)

    pattern594 = Pattern(Integral(acsc(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule594 = ReplacementRule(pattern594, lambda b, a, x : Int(S(1)/(sqrt(S(1) - S(1)/(a + b*x)**S(2))*(a + b*x)), x) + (a + b*x)*acsc(a + b*x)/b)
    rubi.add(rule594)

    pattern595 = Pattern(Integral(asec(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule595 = ReplacementRule(pattern595, lambda b, a, x, n : Subst(Int(x**n*tan(x)*sec(x), x), x, asec(a + b*x))/b)
    rubi.add(rule595)

    pattern596 = Pattern(Integral(acsc(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule596 = ReplacementRule(pattern596, lambda b, a, x, n : -Subst(Int(x**n*cot(x)*csc(x), x), x, acsc(a + b*x))/b)
    rubi.add(rule596)

    pattern597 = Pattern(Integral(asec(a_ + x_*WC('b', S(1)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule597 = ReplacementRule(pattern597, lambda b, a, x : -ImaginaryI*PolyLog(S(2), (-sqrt(-a**S(2) + S(1)) + S(1))*exp(ImaginaryI*asec(a + b*x))/a) - ImaginaryI*PolyLog(S(2), (sqrt(-a**S(2) + S(1)) + S(1))*exp(ImaginaryI*asec(a + b*x))/a) + ImaginaryI*PolyLog(S(2), -exp(S(2)*ImaginaryI*asec(a + b*x)))/S(2) + log(S(1) - (-sqrt(-a**S(2) + S(1)) + S(1))*exp(ImaginaryI*asec(a + b*x))/a)*asec(a + b*x) + log(S(1) - (sqrt(-a**S(2) + S(1)) + S(1))*exp(ImaginaryI*asec(a + b*x))/a)*asec(a + b*x) - log(exp(S(2)*ImaginaryI*asec(a + b*x)) + S(1))*asec(a + b*x))
    rubi.add(rule597)

    pattern598 = Pattern(Integral(acsc(a_ + x_*WC('b', S(1)))/x_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule598 = ReplacementRule(pattern598, lambda b, a, x : ImaginaryI*PolyLog(S(2), ImaginaryI*(-sqrt(-a**S(2) + S(1)) + S(1))*exp(-ImaginaryI*acsc(a + b*x))/a) + ImaginaryI*PolyLog(S(2), ImaginaryI*(sqrt(-a**S(2) + S(1)) + S(1))*exp(-ImaginaryI*acsc(a + b*x))/a) + ImaginaryI*PolyLog(S(2), exp(S(2)*ImaginaryI*acsc(a + b*x)))/S(2) + ImaginaryI*acsc(a + b*x)**S(2) + log(-ImaginaryI*(-sqrt(-a**S(2) + S(1)) + S(1))*exp(-ImaginaryI*acsc(a + b*x))/a + S(1))*acsc(a + b*x) + log(-ImaginaryI*(sqrt(-a**S(2) + S(1)) + S(1))*exp(-ImaginaryI*acsc(a + b*x))/a + S(1))*acsc(a + b*x) - log(-exp(S(2)*ImaginaryI*acsc(a + b*x)) + S(1))*acsc(a + b*x))
    rubi.add(rule598)

    pattern599 = Pattern(Integral(x_**WC('m', S(1))*asec(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule599 = ReplacementRule(pattern599, lambda b, a, m, x : b**(-m + S(-1))*(b**(m + S(1))*x**(m + S(1)) - (-a)**(m + S(1)))*asec(a + b*x)/(m + S(1)) - b**(-m + S(-1))*Subst(Int(x**(-m + S(-1))*((-a*x)**(m + S(1)) - (-a*x + S(1))**(m + S(1)))/sqrt(-x**S(2) + S(1)), x), x, 1/(a + b*x))/(m + S(1)))
    rubi.add(rule599)

    pattern600 = Pattern(Integral(x_**WC('m', S(1))*acsc(a_ + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule600 = ReplacementRule(pattern600, lambda b, a, m, x : b**(-m + S(-1))*(b**(m + S(1))*x**(m + S(1)) - (-a)**(m + S(1)))*acsc(a + b*x)/(m + S(1)) + b**(-m + S(-1))*Subst(Int(x**(-m + S(-1))*((-a*x)**(m + S(1)) - (-a*x + S(1))**(m + S(1)))/sqrt(-x**S(2) + S(1)), x), x, 1/(a + b*x))/(m + S(1)))
    rubi.add(rule600)

    pattern601 = Pattern(Integral(x_**WC('m', S(1))*asec(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule601 = ReplacementRule(pattern601, lambda b, m, n, a, x : b**(-m + S(-1))*Subst(Int(x**n*(-a + sec(x))**m*tan(x)*sec(x), x), x, asec(a + b*x)))
    rubi.add(rule601)

    pattern602 = Pattern(Integral(x_**WC('m', S(1))*acsc(a_ + x_*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule602 = ReplacementRule(pattern602, lambda b, m, n, a, x : -b**(-m + S(-1))*Subst(Int(x**n*(-a + csc(x))**m*cot(x)*csc(x), x), x, acsc(a + b*x)))
    rubi.add(rule602)

    pattern603 = Pattern(Integral(WC('u', S(1))*asec(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule603 = ReplacementRule(pattern603, lambda b, c, m, n, a, u, x : Int(u*acos(a/c + b*x**n/c)**m, x))
    rubi.add(rule603)

    pattern604 = Pattern(Integral(WC('u', S(1))*acsc(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule604 = ReplacementRule(pattern604, lambda b, c, m, n, a, u, x : Int(u*asin(a/c + b*x**n/c)**m, x))
    rubi.add(rule604)

    pattern605 = Pattern(Integral(f_**(WC('c', S(1))*asec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule605 = ReplacementRule(pattern605, lambda b, c, n, a, f, u, x : Subst(Int(f**(c*x**n)*ReplaceAll(u, Rule(x, -a/b + sec(x)/b))*tan(x)*sec(x), x), x, asec(a + b*x))/b)
    rubi.add(rule605)

    pattern606 = Pattern(Integral(f_**(WC('c', S(1))*acsc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule606 = ReplacementRule(pattern606, lambda b, c, n, a, f, u, x : -Subst(Int(f**(c*x**n)*ReplaceAll(u, Rule(x, -a/b + csc(x)/b))*cot(x)*csc(x), x), x, acsc(a + b*x))/b)
    rubi.add(rule606)

    pattern607 = Pattern(Integral(asec(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda u, x: Not(FunctionOfExponentialQ(u, x))))
    rule607 = ReplacementRule(pattern607, lambda u, x : -u*Int(SimplifyIntegrand(x*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x)/sqrt(u**S(2)) + x*asec(u))
    rubi.add(rule607)

    pattern608 = Pattern(Integral(acsc(u_), x_), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda u, x: Not(FunctionOfExponentialQ(u, x))))
    rule608 = ReplacementRule(pattern608, lambda u, x : u*Int(SimplifyIntegrand(x*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x)/sqrt(u**S(2)) + x*acsc(u))
    rubi.add(rule608)

    pattern609 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(u_)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, m, d, u, x: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, x: Not(FunctionOfExponentialQ(u, x))))
    rule609 = ReplacementRule(pattern609, lambda b, c, m, d, a, u, x : -b*u*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x)/(d*(m + S(1))*sqrt(u**S(2))) + (a + b*asec(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule609)

    pattern610 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(u_)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda c, m, d, u, x: Not(FunctionOfQ((c + d*x)**(m + S(1)), u, x))), CustomConstraint(lambda u, x: Not(FunctionOfExponentialQ(u, x))))
    rule610 = ReplacementRule(pattern610, lambda b, c, m, d, a, u, x : b*u*Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x)/(d*(m + S(1))*sqrt(u**S(2))) + (a + b*acsc(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))))
    rubi.add(rule610)

    pattern611 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*asec(u_)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda w, u, x, a, b: InverseFunctionFreeQ(w, x)))
    def With611(b, a, u, v, x):
        w = IntHide(v, x)
        return -b*u*Int(SimplifyIntegrand(w*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x)/sqrt(u**S(2)) + Dist(a + b*asec(u), w, x)
    rule611 = ReplacementRule(pattern611, lambda b, a, u, v, x : With611(b, a, u, v, x))
    rubi.add(rule611)

    pattern612 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acsc(u_)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda u, x: InverseFunctionFreeQ(u, x)), CustomConstraint(lambda v, x: Not(MatchQ(v, Condition((x*Optional(Pattern(d, Blank)) + Optional(Pattern(c, Blank)))**Optional(Pattern(m, Blank)))))), CustomConstraint(lambda w, u, x, a, b: InverseFunctionFreeQ(w, x)))
    def With612(b, a, u, v, x):
        w = IntHide(v, x)
        return b*u*Int(SimplifyIntegrand(w*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x)/sqrt(u**S(2)) + Dist(a + b*acsc(u), w, x)
    rule612 = ReplacementRule(pattern612, lambda b, a, u, v, x : With612(b, a, u, v, x))
    rubi.add(rule612)

    return rubi
