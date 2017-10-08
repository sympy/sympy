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

def hyperbolic(rubi):

    pattern1 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule1 = ReplacementRule(pattern1, lambda m, e, x, f, d, c : -d*m*Int((c + d*x)**(m + S(-1))*Cosh(e + f*x), x)/f + (c + d*x)**m*Cosh(e + f*x)/f)
    rubi.add(rule1)

    pattern2 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cosh(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule2 = ReplacementRule(pattern2, lambda m, e, x, f, d, c : -d*m*Int((c + d*x)**(m + S(-1))*sinh(e + f*x), x)/f + (c + d*x)**m*sinh(e + f*x)/f)
    rubi.add(rule2)

    pattern3 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*sinh(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule3 = ReplacementRule(pattern3, lambda m, e, x, f, d, c : -f*Int((c + d*x)**(m + S(1))*Cosh(e + f*x), x)/(d*(m + S(1))) + (c + d*x)**(m + S(1))*sinh(e + f*x)/(d*(m + S(1))))
    rubi.add(rule3)

    pattern4 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*Cosh(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule4 = ReplacementRule(pattern4, lambda m, e, x, f, d, c : -f*Int((c + d*x)**(m + S(1))*sinh(e + f*x), x)/(d*(m + S(1))) + (c + d*x)**(m + S(1))*Cosh(e + f*x)/(d*(m + S(1))))
    rubi.add(rule4)

    pattern5 = Pattern(Integral(sinh(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, f, e, c: ZeroQ(-c*f + d*e)))
    rule5 = ReplacementRule(pattern5, lambda e, x, f, d, c : SinhIntegral(e + f*x)/d)
    rubi.add(rule5)

    pattern6 = Pattern(Integral(Cosh(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, f, e, c: ZeroQ(-c*f + d*e)))
    rule6 = ReplacementRule(pattern6, lambda e, x, f, d, c : CoshIntegral(e + f*x)/d)
    rubi.add(rule6)

    pattern7 = Pattern(Integral(sinh(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, f, e, c: NonzeroQ(-c*f + d*e)))
    rule7 = ReplacementRule(pattern7, lambda e, x, f, d, c : Cosh((-c*f + d*e)/d)*Int(sinh(c*f/d + f*x)/(c + d*x), x) + Int(Cosh(c*f/d + f*x)/(c + d*x), x)*sinh((-c*f + d*e)/d))
    rubi.add(rule7)

    pattern8 = Pattern(Integral(Cosh(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda d, f, e, c: NonzeroQ(-c*f + d*e)))
    rule8 = ReplacementRule(pattern8, lambda e, x, f, d, c : Cosh((-c*f + d*e)/d)*Int(Cosh(c*f/d + f*x)/(c + d*x), x) + Int(sinh(c*f/d + f*x)/(c + d*x), x)*sinh((-c*f + d*e)/d))
    rubi.add(rule8)

    pattern9 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule9 = ReplacementRule(pattern9, lambda m, e, x, f, d, c : -Int((c + d*x)**m*exp(-e - f*x), x)/S(2) + Int((c + d*x)**m*exp(e + f*x), x)/S(2))
    rubi.add(rule9)

    pattern10 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cosh(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule10 = ReplacementRule(pattern10, lambda m, e, x, f, d, c : Int((c + d*x)**m*exp(-e - f*x), x)/S(2) + Int((c + d*x)**m*exp(e + f*x), x)/S(2))
    rubi.add(rule10)

    pattern11 = Pattern(Integral((WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule11 = ReplacementRule(pattern11, lambda c, e, x, f, d, n, b : -b**S(2)*(n + S(-1))*Int((b*sinh(e + f*x))**(n + S(-2))*(c + d*x), x)/n + b*(b*sinh(e + f*x))**(n + S(-1))*(c + d*x)*Cosh(e + f*x)/(f*n) - d*(b*sinh(e + f*x))**n/(f**S(2)*n**S(2)))
    rubi.add(rule11)

    pattern12 = Pattern(Integral((Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule12 = ReplacementRule(pattern12, lambda e, x, f, d, n, b, c : b**S(2)*(n + S(-1))*Int((b*Cosh(e + f*x))**(n + S(-2))*(c + d*x), x)/n + b*(b*Cosh(e + f*x))**(n + S(-1))*(c + d*x)*sinh(e + f*x)/(f*n) - d*(b*Cosh(e + f*x))**n/(f**S(2)*n**S(2)))
    rubi.add(rule12)

    pattern13 = Pattern(Integral((WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule13 = ReplacementRule(pattern13, lambda m, c, e, x, f, d, n, b : -b**S(2)*(n + S(-1))*Int((b*sinh(e + f*x))**(n + S(-2))*(c + d*x)**m, x)/n + b*(b*sinh(e + f*x))**(n + S(-1))*(c + d*x)**m*Cosh(e + f*x)/(f*n) + d**S(2)*m*(m + S(-1))*Int((b*sinh(e + f*x))**n*(c + d*x)**(m + S(-2)), x)/(f**S(2)*n**S(2)) - d*m*(b*sinh(e + f*x))**n*(c + d*x)**(m + S(-1))/(f**S(2)*n**S(2)))
    rubi.add(rule13)

    pattern14 = Pattern(Integral((Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule14 = ReplacementRule(pattern14, lambda m, e, x, f, d, n, b, c : b**S(2)*(n + S(-1))*Int((b*Cosh(e + f*x))**(n + S(-2))*(c + d*x)**m, x)/n + b*(b*Cosh(e + f*x))**(n + S(-1))*(c + d*x)**m*sinh(e + f*x)/(f*n) + d**S(2)*m*(m + S(-1))*Int((b*Cosh(e + f*x))**n*(c + d*x)**(m + S(-2)), x)/(f**S(2)*n**S(2)) - d*m*(b*Cosh(e + f*x))**n*(c + d*x)**(m + S(-1))/(f**S(2)*n**S(2)))
    rubi.add(rule14)

    pattern15 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*sinh(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: Not(RationalQ(m)) | Inequality(S(-1), LessEqual, m, Less, S(1))))
    rule15 = ReplacementRule(pattern15, lambda m, e, x, f, d, n, c : Int(ExpandTrigReduce((c + d*x)**m, sinh(e + f*x)**n, x), x))
    rubi.add(rule15)

    pattern16 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*Cosh(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: Not(RationalQ(m)) | Inequality(S(-1), LessEqual, m, Less, S(1))))
    rule16 = ReplacementRule(pattern16, lambda m, e, x, f, d, n, c : Int(ExpandTrigReduce((c + d*x)**m, Cosh(e + f*x)**n, x), x))
    rubi.add(rule16)

    pattern17 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*sinh(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Inequality(S(-2), LessEqual, m, Less, S(-1))))
    rule17 = ReplacementRule(pattern17, lambda m, e, x, f, d, n, c : -f*n*Int(ExpandTrigReduce((c + d*x)**(m + S(1)), Cosh(e + f*x)*sinh(e + f*x)**(n + S(-1)), x), x)/(d*(m + S(1))) + (c + d*x)**(m + S(1))*sinh(e + f*x)**n/(d*(m + S(1))))
    rubi.add(rule17)

    pattern18 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*Cosh(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Inequality(S(-2), LessEqual, m, Less, S(-1))))
    rule18 = ReplacementRule(pattern18, lambda m, e, x, f, d, n, c : -f*n*Int(ExpandTrigReduce((c + d*x)**(m + S(1)), Cosh(e + f*x)**(n + S(-1))*sinh(e + f*x), x), x)/(d*(m + S(1))) + (c + d*x)**(m + S(1))*Cosh(e + f*x)**n/(d*(m + S(1))))
    rubi.add(rule18)

    pattern19 = Pattern(Integral((WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: Less(m, S(-2))))
    rule19 = ReplacementRule(pattern19, lambda m, c, e, x, f, d, n, b : b**S(2)*f**S(2)*n*(n + S(-1))*Int((b*sinh(e + f*x))**(n + S(-2))*(c + d*x)**(m + S(2)), x)/(d**S(2)*(m + S(1))*(m + S(2))) - b*f*n*(b*sinh(e + f*x))**(n + S(-1))*(c + d*x)**(m + S(2))*Cosh(e + f*x)/(d**S(2)*(m + S(1))*(m + S(2))) + (b*sinh(e + f*x))**n*(c + d*x)**(m + S(1))/(d*(m + S(1))) + f**S(2)*n**S(2)*Int((b*sinh(e + f*x))**n*(c + d*x)**(m + S(2)), x)/(d**S(2)*(m + S(1))*(m + S(2))))
    rubi.add(rule19)

    pattern20 = Pattern(Integral((Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: Less(m, S(-2))))
    rule20 = ReplacementRule(pattern20, lambda m, e, x, f, d, n, b, c : -b**S(2)*f**S(2)*n*(n + S(-1))*Int((b*Cosh(e + f*x))**(n + S(-2))*(c + d*x)**(m + S(2)), x)/(d**S(2)*(m + S(1))*(m + S(2))) - b*f*n*(b*Cosh(e + f*x))**(n + S(-1))*(c + d*x)**(m + S(2))*sinh(e + f*x)/(d**S(2)*(m + S(1))*(m + S(2))) + (b*Cosh(e + f*x))**n*(c + d*x)**(m + S(1))/(d*(m + S(1))) + f**S(2)*n**S(2)*Int((b*Cosh(e + f*x))**n*(c + d*x)**(m + S(2)), x)/(d**S(2)*(m + S(1))*(m + S(2))))
    rubi.add(rule20)

    pattern21 = Pattern(Integral((WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule21 = ReplacementRule(pattern21, lambda c, e, x, f, d, n, b : (b*sinh(e + f*x))**(n + S(1))*(c + d*x)*Cosh(e + f*x)/(b*f*(n + S(1))) - d*(b*sinh(e + f*x))**(n + S(2))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))) - (n + S(2))*Int((b*sinh(e + f*x))**(n + S(2))*(c + d*x), x)/(b**S(2)*(n + S(1))))
    rubi.add(rule21)

    pattern22 = Pattern(Integral((Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule22 = ReplacementRule(pattern22, lambda e, x, f, d, n, b, c : (b*Cosh(e + f*x))**(n + S(1))*(-c - d*x)*sinh(e + f*x)/(b*f*(n + S(1))) + d*(b*Cosh(e + f*x))**(n + S(2))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))) + (n + S(2))*Int((b*Cosh(e + f*x))**(n + S(2))*(c + d*x), x)/(b**S(2)*(n + S(1))))
    rubi.add(rule22)

    pattern23 = Pattern(Integral((WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule23 = ReplacementRule(pattern23, lambda m, c, e, x, f, d, n, b : (b*sinh(e + f*x))**(n + S(1))*(c + d*x)**m*Cosh(e + f*x)/(b*f*(n + S(1))) + d**S(2)*m*(m + S(-1))*Int((b*sinh(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-2)), x)/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))) - d*m*(b*sinh(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-1))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))) - (n + S(2))*Int((b*sinh(e + f*x))**(n + S(2))*(c + d*x)**m, x)/(b**S(2)*(n + S(1))))
    rubi.add(rule23)

    pattern24 = Pattern(Integral((Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule24 = ReplacementRule(pattern24, lambda m, e, x, f, d, n, b, c : -(b*Cosh(e + f*x))**(n + S(1))*(c + d*x)**m*sinh(e + f*x)/(b*f*(n + S(1))) - d**S(2)*m*(m + S(-1))*Int((b*Cosh(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-2)), x)/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))) + d*m*(b*Cosh(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-1))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))) + (n + S(2))*Int((b*Cosh(e + f*x))**(n + S(2))*(c + d*x)**m, x)/(b**S(2)*(n + S(1))))
    rubi.add(rule24)

    pattern25 = Pattern(Integral((a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, m, b, a: PositiveIntegerQ(m) | Equal(n, S(1)) | NonzeroQ(a**S(2) + b**S(2))))
    rule25 = ReplacementRule(pattern25, lambda m, c, e, x, f, d, a, n, b : Int(ExpandIntegrand((c + d*x)**m, (a + b*sinh(e + f*x))**n, x), x))
    rubi.add(rule25)

    pattern26 = Pattern(Integral((a_ + Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, m, b, a: PositiveIntegerQ(m) | Equal(n, S(1)) | NonzeroQ(a**S(2) - b**S(2))))
    rule26 = ReplacementRule(pattern26, lambda m, e, x, f, d, a, n, b, c : Int(ExpandIntegrand((c + d*x)**m, (a + b*Cosh(e + f*x))**n, x), x))
    rubi.add(rule26)

    pattern27 = Pattern(Integral((a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) + b**S(2))), CustomConstraint(lambda n: IntegerQ(n)))
    rule27 = ReplacementRule(pattern27, lambda m, c, e, x, f, d, a, n, b : (S(2)*a)**n*Int((c + d*x)**m*Cosh(-Pi*a/(S(4)*b) + e/S(2) + f*x/S(2))**(S(2)*n), x))
    rubi.add(rule27)

    pattern28 = Pattern(Integral((a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) + b**S(2))), CustomConstraint(lambda n: IntegerQ(S(2)*n)), CustomConstraint(lambda n, m: PositiveIntegerQ(m) | Greater(n, S(0))))
    rule28 = ReplacementRule(pattern28, lambda m, c, e, x, f, d, a, n, b : (S(2)*a)**IntPart(n)*(a + b*sinh(e + f*x))**FracPart(n)*Cosh(-Pi*a/(S(4)*b) + e/S(2) + f*x/S(2))**(-S(2)*FracPart(n))*Int((c + d*x)**m*Cosh(-Pi*a/(S(4)*b) + e/S(2) + f*x/S(2))**(S(2)*n), x))
    rubi.add(rule28)

    pattern29 = Pattern(Integral((a_ + Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, a: ZeroQ(a - b)), CustomConstraint(lambda n: IntegerQ(n)))
    rule29 = ReplacementRule(pattern29, lambda m, e, x, f, d, a, n, b, c : (S(2)*a)**n*Int((c + d*x)**m*Cosh(e/S(2) + f*x/S(2))**(S(2)*n), x))
    rubi.add(rule29)

    pattern30 = Pattern(Integral((a_ + Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, a: ZeroQ(a + b)), CustomConstraint(lambda n: IntegerQ(n)))
    rule30 = ReplacementRule(pattern30, lambda m, e, x, f, d, a, n, b, c : (-S(2)*a)**n*Int((c + d*x)**m*sinh(e/S(2) + f*x/S(2))**(S(2)*n), x))
    rubi.add(rule30)

    pattern31 = Pattern(Integral((a_ + Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, a: ZeroQ(a - b)), CustomConstraint(lambda n: IntegerQ(S(2)*n)), CustomConstraint(lambda n, m: PositiveIntegerQ(m) | Greater(n, S(0))))
    rule31 = ReplacementRule(pattern31, lambda m, e, x, f, d, a, n, b, c : (S(2)*a)**IntPart(n)*(a + b*Cosh(e + f*x))**FracPart(n)*Cosh(e/S(2) + f*x/S(2))**(-S(2)*FracPart(n))*Int((c + d*x)**m*Cosh(e/S(2) + f*x/S(2))**(S(2)*n), x))
    rubi.add(rule31)

    pattern32 = Pattern(Integral((a_ + Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, a: ZeroQ(a + b)), CustomConstraint(lambda n: IntegerQ(S(2)*n)), CustomConstraint(lambda n, m: PositiveIntegerQ(m) | Greater(n, S(0))))
    rule32 = ReplacementRule(pattern32, lambda m, e, x, f, d, a, n, b, c : (-S(2)*a)**IntPart(n)*(a + b*Cosh(e + f*x))**FracPart(n)*Int((c + d*x)**m*sinh(e/S(2) + f*x/S(2))**(S(2)*n), x)*sinh(e/S(2) + f*x/S(2))**(-S(2)*FracPart(n)))
    rubi.add(rule32)

    pattern33 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) + b**S(2))), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule33 = ReplacementRule(pattern33, lambda m, c, e, x, f, d, a, b : -S(2)*Int((c + d*x)**m*exp(e + f*x)/(-S(2)*a*exp(e + f*x) - b*exp(S(2)*e + S(2)*f*x) + b), x))
    rubi.add(rule33)

    pattern34 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule34 = ReplacementRule(pattern34, lambda m, e, x, f, d, a, b, c : S(2)*Int((c + d*x)**m*exp(e + f*x)/(S(2)*a*exp(e + f*x) + b*exp(S(2)*e + S(2)*f*x) + b), x))
    rubi.add(rule34)

    pattern35 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) + b**S(2))), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule35 = ReplacementRule(pattern35, lambda m, c, e, x, f, d, a, b : a*Int((c + d*x)**m/(a + b*sinh(e + f*x)), x)/(a**S(2) + b**S(2)) + b*d*m*Int((c + d*x)**(m + S(-1))*Cosh(e + f*x)/(a + b*sinh(e + f*x)), x)/(f*(a**S(2) + b**S(2))) - b*(c + d*x)**m*Cosh(e + f*x)/(f*(a + b*sinh(e + f*x))*(a**S(2) + b**S(2))))
    rubi.add(rule35)

    pattern36 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule36 = ReplacementRule(pattern36, lambda m, e, x, f, d, a, b, c : a*Int((c + d*x)**m/(a + b*Cosh(e + f*x)), x)/(a**S(2) - b**S(2)) + b*d*m*Int((c + d*x)**(m + S(-1))*sinh(e + f*x)/(a + b*Cosh(e + f*x)), x)/(f*(a**S(2) - b**S(2))) - b*(c + d*x)**m*sinh(e + f*x)/(f*(a + b*Cosh(e + f*x))*(a**S(2) - b**S(2))))
    rubi.add(rule36)

    pattern37 = Pattern(Integral((a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) + b**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n + S(2))), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule37 = ReplacementRule(pattern37, lambda m, c, e, x, f, d, a, n, b : a*Int((a + b*sinh(e + f*x))**(n + S(1))*(c + d*x)**m, x)/(a**S(2) + b**S(2)) - b*d*m*Int((a + b*sinh(e + f*x))**(n + S(1))*(c + d*x)**(m + S(-1))*Cosh(e + f*x), x)/(f*(a**S(2) + b**S(2))*(n + S(1))) - b*(n + S(2))*Int((a + b*sinh(e + f*x))**(n + S(1))*(c + d*x)**m*sinh(e + f*x), x)/((a**S(2) + b**S(2))*(n + S(1))) + b*(a + b*sinh(e + f*x))**(n + S(1))*(c + d*x)**m*Cosh(e + f*x)/(f*(a**S(2) + b**S(2))*(n + S(1))))
    rubi.add(rule37)

    pattern38 = Pattern(Integral((a_ + Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n + S(2))), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule38 = ReplacementRule(pattern38, lambda m, e, x, f, d, a, n, b, c : a*Int((a + b*Cosh(e + f*x))**(n + S(1))*(c + d*x)**m, x)/(a**S(2) - b**S(2)) - b*d*m*Int((a + b*Cosh(e + f*x))**(n + S(1))*(c + d*x)**(m + S(-1))*sinh(e + f*x), x)/(f*(a**S(2) - b**S(2))*(n + S(1))) - b*(n + S(2))*Int((a + b*Cosh(e + f*x))**(n + S(1))*(c + d*x)**m*Cosh(e + f*x), x)/((a**S(2) - b**S(2))*(n + S(1))) + b*(a + b*Cosh(e + f*x))**(n + S(1))*(c + d*x)**m*sinh(e + f*x)/(f*(a**S(2) - b**S(2))*(n + S(1))))
    rubi.add(rule38)

    pattern39 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(v_))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda v, u, x: LinearQ(List(u, v), x)), CustomConstraint(lambda v, u, x: Not(LinearMatchQ(List(u, v), x))))
    rule39 = ReplacementRule(pattern39, lambda m, x, a, n, v, u, b : Int((a + b*sinh(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x))
    rubi.add(rule39)

    pattern40 = Pattern(Integral(u_**WC('m', S(1))*(Cosh(v_)*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda v, u, x: LinearQ(List(u, v), x)), CustomConstraint(lambda v, u, x: Not(LinearMatchQ(List(u, v), x))))
    rule40 = ReplacementRule(pattern40, lambda m, x, a, n, v, u, b : Int((a + b*Cosh(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x))
    rubi.add(rule40)

    pattern41 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule41 = ReplacementRule(pattern41, lambda m, c, e, x, a, d, f, n, b : Int((a + b*sinh(e + f*x))**n*(c + d*x)**m, x))
    rubi.add(rule41)

    pattern42 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(Cosh(x_*WC('f', S(1)) + WC('e', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule42 = ReplacementRule(pattern42, lambda m, c, e, x, f, a, d, n, b : Int((a + b*Cosh(e + f*x))**n*(c + d*x)**m, x))
    rubi.add(rule42)

    pattern43 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule43 = ReplacementRule(pattern43, lambda p, x, a, d, n, b, c : Int(ExpandIntegrand(sinh(c + d*x), (a + b*x**n)**p, x), x))
    rubi.add(rule43)

    pattern44 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule44 = ReplacementRule(pattern44, lambda p, c, x, a, d, n, b : Int(ExpandIntegrand(Cosh(c + d*x), (a + b*x**n)**p, x), x))
    rubi.add(rule44)

    pattern45 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(2))))
    rule45 = ReplacementRule(pattern45, lambda p, x, a, d, n, b, c : -d*Int(x**(-n + S(1))*(a + b*x**n)**(p + S(1))*Cosh(c + d*x), x)/(b*n*(p + S(1))) + x**(-n + S(1))*(a + b*x**n)**(p + S(1))*sinh(c + d*x)/(b*n*(p + S(1))) - (-n + S(1))*Int(x**(-n)*(a + b*x**n)**(p + S(1))*sinh(c + d*x), x)/(b*n*(p + S(1))))
    rubi.add(rule45)

    pattern46 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n: Greater(n, S(2))))
    rule46 = ReplacementRule(pattern46, lambda p, c, x, a, d, n, b : -d*Int(x**(-n + S(1))*(a + b*x**n)**(p + S(1))*sinh(c + d*x), x)/(b*n*(p + S(1))) + x**(-n + S(1))*(a + b*x**n)**(p + S(1))*Cosh(c + d*x)/(b*n*(p + S(1))) - (-n + S(1))*Int(x**(-n)*(a + b*x**n)**(p + S(1))*Cosh(c + d*x), x)/(b*n*(p + S(1))))
    rubi.add(rule46)

    pattern47 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, p: Equal(n, S(2)) | Equal(p, S(-1))))
    rule47 = ReplacementRule(pattern47, lambda p, x, a, d, n, b, c : Int(ExpandIntegrand(sinh(c + d*x), (a + b*x**n)**p, x), x))
    rubi.add(rule47)

    pattern48 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, p: Equal(n, S(2)) | Equal(p, S(-1))))
    rule48 = ReplacementRule(pattern48, lambda p, c, x, a, d, n, b : Int(ExpandIntegrand(Cosh(c + d*x), (a + b*x**n)**p, x), x))
    rubi.add(rule48)

    pattern49 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule49 = ReplacementRule(pattern49, lambda p, x, a, d, n, b, c : Int(x**(n*p)*(a*x**(-n) + b)**p*sinh(c + d*x), x))
    rubi.add(rule49)

    pattern50 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule50 = ReplacementRule(pattern50, lambda p, c, x, a, d, n, b : Int(x**(n*p)*(a*x**(-n) + b)**p*Cosh(c + d*x), x))
    rubi.add(rule50)

    pattern51 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule51 = ReplacementRule(pattern51, lambda p, x, a, d, n, b, c : Int((a + b*x**n)**p*sinh(c + d*x), x))
    rubi.add(rule51)

    pattern52 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule52 = ReplacementRule(pattern52, lambda p, c, x, a, d, n, b : Int((a + b*x**n)**p*Cosh(c + d*x), x))
    rubi.add(rule52)

    pattern53 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule53 = ReplacementRule(pattern53, lambda p, m, e, x, a, d, n, b, c : Int(ExpandIntegrand(sinh(c + d*x), (e*x)**m*(a + b*x**n)**p, x), x))
    rubi.add(rule53)

    pattern54 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule54 = ReplacementRule(pattern54, lambda p, m, c, e, x, a, d, n, b : Int(ExpandIntegrand(Cosh(c + d*x), (e*x)**m*(a + b*x**n)**p, x), x))
    rubi.add(rule54)

    pattern55 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, e: IntegerQ(n) | PositiveQ(e)))
    rule55 = ReplacementRule(pattern55, lambda p, m, e, x, a, d, n, b, c : -d*e**m*Int((a + b*x**n)**(p + S(1))*Cosh(c + d*x), x)/(b*n*(p + S(1))) + e**m*(a + b*x**n)**(p + S(1))*sinh(c + d*x)/(b*n*(p + S(1))))
    rubi.add(rule55)

    pattern56 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, e: IntegerQ(n) | PositiveQ(e)))
    rule56 = ReplacementRule(pattern56, lambda p, m, c, e, x, a, d, n, b : -d*e**m*Int((a + b*x**n)**(p + S(1))*sinh(c + d*x), x)/(b*n*(p + S(1))) + e**m*(a + b*x**n)**(p + S(1))*Cosh(c + d*x)/(b*n*(p + S(1))))
    rubi.add(rule56)

    pattern57 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, m: Greater(n, S(2)) | Greater(m - n + S(1), S(0))))
    rule57 = ReplacementRule(pattern57, lambda p, m, x, a, d, n, b, c : -d*Int(x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*Cosh(c + d*x), x)/(b*n*(p + S(1))) + x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*sinh(c + d*x)/(b*n*(p + S(1))) - (m - n + S(1))*Int(x**(m - n)*(a + b*x**n)**(p + S(1))*sinh(c + d*x), x)/(b*n*(p + S(1))))
    rubi.add(rule57)

    pattern58 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, m: Greater(n, S(2)) | Greater(m - n + S(1), S(0))))
    rule58 = ReplacementRule(pattern58, lambda p, m, c, x, a, d, n, b : -d*Int(x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*sinh(c + d*x), x)/(b*n*(p + S(1))) + x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*Cosh(c + d*x)/(b*n*(p + S(1))) - (m - n + S(1))*Int(x**(m - n)*(a + b*x**n)**(p + S(1))*Cosh(c + d*x), x)/(b*n*(p + S(1))))
    rubi.add(rule58)

    pattern59 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, p: Equal(n, S(2)) | Equal(p, S(-1))))
    rule59 = ReplacementRule(pattern59, lambda p, m, x, a, d, n, b, c : Int(ExpandIntegrand(sinh(c + d*x), x**m*(a + b*x**n)**p, x), x))
    rubi.add(rule59)

    pattern60 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, p: Equal(n, S(2)) | Equal(p, S(-1))))
    rule60 = ReplacementRule(pattern60, lambda p, m, c, x, a, d, n, b : Int(ExpandIntegrand(Cosh(c + d*x), x**m*(a + b*x**n)**p, x), x))
    rubi.add(rule60)

    pattern61 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule61 = ReplacementRule(pattern61, lambda p, m, x, a, d, n, b, c : Int(x**(m + n*p)*(a*x**(-n) + b)**p*sinh(c + d*x), x))
    rubi.add(rule61)

    pattern62 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule62 = ReplacementRule(pattern62, lambda p, m, c, x, a, d, n, b : Int(x**(m + n*p)*(a*x**(-n) + b)**p*Cosh(c + d*x), x))
    rubi.add(rule62)

    pattern63 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule63 = ReplacementRule(pattern63, lambda p, m, e, x, a, d, n, b, c : Int((e*x)**m*(a + b*x**n)**p*sinh(c + d*x), x))
    rubi.add(rule63)

    pattern64 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule64 = ReplacementRule(pattern64, lambda p, m, c, e, x, a, d, n, b : Int((e*x)**m*(a + b*x**n)**p*Cosh(c + d*x), x))
    rubi.add(rule64)

    pattern65 = Pattern(Integral(sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule65 = ReplacementRule(pattern65, lambda d, n, x, c : -Int(exp(-c - d*x**n), x)/S(2) + Int(exp(c + d*x**n), x)/S(2))
    rubi.add(rule65)

    pattern66 = Pattern(Integral(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule66 = ReplacementRule(pattern66, lambda d, n, x, c : Int(exp(-c - d*x**n), x)/S(2) + Int(exp(c + d*x**n), x)/S(2))
    rubi.add(rule66)

    pattern67 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, p: IntegersQ(n, p)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda p: Greater(p, S(1))))
    rule67 = ReplacementRule(pattern67, lambda p, x, a, d, n, b, c : Int(ExpandTrigReduce((a + b*sinh(c + d*x**n))**p, x), x))
    rubi.add(rule67)

    pattern68 = Pattern(Integral((Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, p: IntegersQ(n, p)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda p: Greater(p, S(1))))
    rule68 = ReplacementRule(pattern68, lambda p, x, a, d, n, b, c : Int(ExpandTrigReduce((a + b*Cosh(c + d*x**n))**p, x), x))
    rubi.add(rule68)

    pattern69 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule69 = ReplacementRule(pattern69, lambda p, x, a, d, n, b, c : -Subst(Int((a + b*sinh(c + d*x**(-n)))**p/x**S(2), x), x, 1/x))
    rubi.add(rule69)

    pattern70 = Pattern(Integral((Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule70 = ReplacementRule(pattern70, lambda p, x, a, d, n, b, c : -Subst(Int((a + b*Cosh(c + d*x**(-n)))**p/x**S(2), x), x, 1/x))
    rubi.add(rule70)

    pattern71 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: FractionQ(n)), )
    def With71(p, x, a, d, n, b, c):
        k = Denominator(n)
        return k*Subst(Int(x**(k + S(-1))*(a + b*sinh(c + d*x**(k*n)))**p, x), x, x**(1/k))
    rule71 = ReplacementRule(pattern71, lambda p, x, a, d, n, b, c : With71(p, x, a, d, n, b, c))
    rubi.add(rule71)

    pattern72 = Pattern(Integral((Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: FractionQ(n)), )
    def With72(p, x, a, d, n, b, c):
        k = Denominator(n)
        return k*Subst(Int(x**(k + S(-1))*(a + b*Cosh(c + d*x**(k*n)))**p, x), x, x**(1/k))
    rule72 = ReplacementRule(pattern72, lambda p, x, a, d, n, b, c : With72(p, x, a, d, n, b, c))
    rubi.add(rule72)

    pattern73 = Pattern(Integral(sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule73 = ReplacementRule(pattern73, lambda d, n, x, c : -Int(exp(-c - d*x**n), x)/S(2) + Int(exp(c + d*x**n), x)/S(2))
    rubi.add(rule73)

    pattern74 = Pattern(Integral(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule74 = ReplacementRule(pattern74, lambda d, n, x, c : Int(exp(-c - d*x**n), x)/S(2) + Int(exp(c + d*x**n), x)/S(2))
    rubi.add(rule74)

    pattern75 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule75 = ReplacementRule(pattern75, lambda p, x, a, d, n, b, c : Int(ExpandTrigReduce((a + b*sinh(c + d*x**n))**p, x), x))
    rubi.add(rule75)

    pattern76 = Pattern(Integral((Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule76 = ReplacementRule(pattern76, lambda p, x, a, d, n, b, c : Int(ExpandTrigReduce((a + b*Cosh(c + d*x**n))**p, x), x))
    rubi.add(rule76)

    pattern77 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda u, x: LinearQ(u, x)), CustomConstraint(lambda u, x: NonzeroQ(u - x)))
    rule77 = ReplacementRule(pattern77, lambda p, x, a, d, n, b, u, c : Subst(Int((a + b*sinh(c + d*x**n))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule77)

    pattern78 = Pattern(Integral((Cosh(u_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda u, x: LinearQ(u, x)), CustomConstraint(lambda u, x: NonzeroQ(u - x)))
    rule78 = ReplacementRule(pattern78, lambda p, x, a, d, n, b, u, c : Subst(Int((a + b*Cosh(c + d*x**n))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule78)

    pattern79 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: LinearQ(u, x)))
    rule79 = ReplacementRule(pattern79, lambda p, x, a, d, n, u, b, c : Int((a + b*sinh(c + d*u**n))**p, x))
    rubi.add(rule79)

    pattern80 = Pattern(Integral((Cosh(u_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: LinearQ(u, x)))
    rule80 = ReplacementRule(pattern80, lambda p, x, a, d, n, u, b, c : Int((a + b*Cosh(c + d*u**n))**p, x))
    rubi.add(rule80)

    pattern81 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(u_))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule81 = ReplacementRule(pattern81, lambda p, x, a, u, b : Int((a + b*sinh(ExpandToSum(u, x)))**p, x))
    rubi.add(rule81)

    pattern82 = Pattern(Integral((Cosh(u_)*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule82 = ReplacementRule(pattern82, lambda p, x, a, u, b : Int((a + b*Cosh(ExpandToSum(u, x)))**p, x))
    rubi.add(rule82)

    pattern83 = Pattern(Integral(sinh(x_**n_*WC('d', S(1)))/x_, x_), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule83 = ReplacementRule(pattern83, lambda d, n, x : SinhIntegral(d*x**n)/n)
    rubi.add(rule83)

    pattern84 = Pattern(Integral(Cosh(x_**n_*WC('d', S(1)))/x_, x_), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule84 = ReplacementRule(pattern84, lambda d, n, x : CoshIntegral(d*x**n)/n)
    rubi.add(rule84)

    pattern85 = Pattern(Integral(sinh(c_ + x_**n_*WC('d', S(1)))/x_, x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule85 = ReplacementRule(pattern85, lambda d, n, x, c : Cosh(c)*Int(sinh(d*x**n)/x, x) + Int(Cosh(d*x**n)/x, x)*sinh(c))
    rubi.add(rule85)

    pattern86 = Pattern(Integral(Cosh(c_ + x_**n_*WC('d', S(1)))/x_, x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule86 = ReplacementRule(pattern86, lambda d, n, x, c : Cosh(c)*Int(Cosh(d*x**n)/x, x) + Int(sinh(d*x**n)/x, x)*sinh(c))
    rubi.add(rule86)

    pattern87 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda d, p, b, c, x, mn, a, n: IntegerQ(mn) & (Equal(p, S(1)) | Greater(mn, S(0)))))
    def With87(p, m, x, a, d, n, b, c):
        mn = (m + S(1))/n
        return Subst(Int(x**(mn + S(-1))*(a + b*sinh(c + d*x))**p, x), x, x**n)/n
    rule87 = ReplacementRule(pattern87, lambda p, m, x, a, d, n, b, c : With87(p, m, x, a, d, n, b, c))
    rubi.add(rule87)

    pattern88 = Pattern(Integral(x_**WC('m', S(1))*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda d, p, b, c, x, mn, a, n: IntegerQ(mn) & (Equal(p, S(1)) | Greater(mn, S(0)))))
    def With88(p, m, x, a, d, n, b, c):
        mn = (m + S(1))/n
        return Subst(Int(x**(mn + S(-1))*(a + b*Cosh(c + d*x))**p, x), x, x**n)/n
    rule88 = ReplacementRule(pattern88, lambda p, m, x, a, d, n, b, c : With88(p, m, x, a, d, n, b, c))
    rubi.add(rule88)

    pattern89 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda d, m, p, c, b, x, mn, e, a, n: IntegerQ(mn) & (Equal(p, S(1)) | Greater(mn, S(0)))))
    def With89(p, m, e, x, a, d, n, b, c):
        mn = (m + S(1))/n
        return e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*sinh(c + d*x**n))**p, x)
    rule89 = ReplacementRule(pattern89, lambda p, m, e, x, a, d, n, b, c : With89(p, m, e, x, a, d, n, b, c))
    rubi.add(rule89)

    pattern90 = Pattern(Integral((e_*x_)**m_*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda d, m, p, c, b, x, mn, e, a, n: IntegerQ(mn) & (Equal(p, S(1)) | Greater(mn, S(0)))))
    def With90(p, m, e, x, a, d, n, b, c):
        mn = (m + S(1))/n
        return e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*Cosh(c + d*x**n))**p, x)
    rule90 = ReplacementRule(pattern90, lambda p, m, e, x, a, d, n, b, c : With90(p, m, e, x, a, d, n, b, c))
    rubi.add(rule90)

    pattern91 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n, m: Less(S(0), n, m + S(1))))
    rule91 = ReplacementRule(pattern91, lambda m, e, x, d, n, c : -e**n*(m - n + S(1))*Int((e*x)**(m - n)*Cosh(c + d*x**n), x)/(d*n) + e**(n + S(-1))*(e*x)**(m - n + S(1))*Cosh(c + d*x**n)/(d*n))
    rubi.add(rule91)

    pattern92 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n, m: Less(S(0), n, m + S(1))))
    rule92 = ReplacementRule(pattern92, lambda m, e, x, d, n, c : -e**n*(m - n + S(1))*Int((e*x)**(m - n)*sinh(c + d*x**n), x)/(d*n) + e**(n + S(-1))*(e*x)**(m - n + S(1))*sinh(c + d*x**n)/(d*n))
    rubi.add(rule92)

    pattern93 = Pattern(Integral((x_*WC('e', S(1)))**m_*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule93 = ReplacementRule(pattern93, lambda m, e, x, d, n, c : -d*e**(-n)*n*Int((e*x)**(m + n)*Cosh(c + d*x**n), x)/(m + S(1)) + (e*x)**(m + S(1))*sinh(c + d*x**n)/(e*(m + S(1))))
    rubi.add(rule93)

    pattern94 = Pattern(Integral((x_*WC('e', S(1)))**m_*Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))))
    rule94 = ReplacementRule(pattern94, lambda m, e, x, d, n, c : -d*e**(-n)*n*Int((e*x)**(m + n)*sinh(c + d*x**n), x)/(m + S(1)) + (e*x)**(m + S(1))*Cosh(c + d*x**n)/(e*(m + S(1))))
    rubi.add(rule94)

    pattern95 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule95 = ReplacementRule(pattern95, lambda m, e, x, d, n, c : -Int((e*x)**m*exp(-c - d*x**n), x)/S(2) + Int((e*x)**m*exp(c + d*x**n), x)/S(2))
    rubi.add(rule95)

    pattern96 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule96 = ReplacementRule(pattern96, lambda m, e, x, d, n, c : Int((e*x)**m*exp(-c - d*x**n), x)/S(2) + Int((e*x)**m*exp(c + d*x**n), x)/S(2))
    rubi.add(rule96)

    pattern97 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, p: IntegersQ(n, p)), CustomConstraint(lambda n, m: ZeroQ(m + n)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda n: NonzeroQ(n + S(-1))))
    rule97 = ReplacementRule(pattern97, lambda p, m, x, a, n, b : b*n*p*Int(Cosh(a + b*x**n)*sinh(a + b*x**n)**(p + S(-1)), x)/(n + S(-1)) - x**(-n + S(1))*sinh(a + b*x**n)**p/(n + S(-1)))
    rubi.add(rule97)

    pattern98 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, p: IntegersQ(n, p)), CustomConstraint(lambda n, m: ZeroQ(m + n)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda n: NonzeroQ(n + S(-1))))
    rule98 = ReplacementRule(pattern98, lambda p, m, x, a, n, b : b*n*p*Int(Cosh(a + b*x**n)**(p + S(-1))*sinh(a + b*x**n), x)/(n + S(-1)) - x**(-n + S(1))*Cosh(a + b*x**n)**p/(n + S(-1)))
    rubi.add(rule98)

    pattern99 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(m - S(2)*n + S(1))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))))
    rule99 = ReplacementRule(pattern99, lambda p, m, x, a, n, b : -(p + S(-1))*Int(x**m*sinh(a + b*x**n)**(p + S(-2)), x)/p + x**n*Cosh(a + b*x**n)*sinh(a + b*x**n)**(p + S(-1))/(b*n*p) - sinh(a + b*x**n)**p/(b**S(2)*n*p**S(2)))
    rubi.add(rule99)

    pattern100 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(m - S(2)*n + S(1))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))))
    rule100 = ReplacementRule(pattern100, lambda p, m, x, a, n, b : (p + S(-1))*Int(x**m*Cosh(a + b*x**n)**(p + S(-2)), x)/p + x**n*Cosh(a + b*x**n)**(p + S(-1))*sinh(a + b*x**n)/(b*n*p) - Cosh(a + b*x**n)**p/(b**S(2)*n*p**S(2)))
    rubi.add(rule100)

    pattern101 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda n, m: Less(S(0), S(2)*n, m + S(1))))
    rule101 = ReplacementRule(pattern101, lambda p, m, x, a, n, b : -(p + S(-1))*Int(x**m*sinh(a + b*x**n)**(p + S(-2)), x)/p + x**(m - n + S(1))*Cosh(a + b*x**n)*sinh(a + b*x**n)**(p + S(-1))/(b*n*p) + x**(m - S(2)*n + S(1))*(-m + n + S(-1))*sinh(a + b*x**n)**p/(b**S(2)*n**S(2)*p**S(2)) + (m - S(2)*n + S(1))*(m - n + S(1))*Int(x**(m - S(2)*n)*sinh(a + b*x**n)**p, x)/(b**S(2)*n**S(2)*p**S(2)))
    rubi.add(rule101)

    pattern102 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda n, m: Less(S(0), S(2)*n, m + S(1))))
    rule102 = ReplacementRule(pattern102, lambda p, m, x, a, n, b : (p + S(-1))*Int(x**m*Cosh(a + b*x**n)**(p + S(-2)), x)/p + x**(m - n + S(1))*Cosh(a + b*x**n)**(p + S(-1))*sinh(a + b*x**n)/(b*n*p) + x**(m - S(2)*n + S(1))*(-m + n + S(-1))*Cosh(a + b*x**n)**p/(b**S(2)*n**S(2)*p**S(2)) + (m - S(2)*n + S(1))*(m - n + S(1))*Int(x**(m - S(2)*n)*Cosh(a + b*x**n)**p, x)/(b**S(2)*n**S(2)*p**S(2)))
    rubi.add(rule102)

    pattern103 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda n, m: Less(S(0), S(2)*n, -m + S(1))), CustomConstraint(lambda n, m: NonzeroQ(m + n + S(1))))
    rule103 = ReplacementRule(pattern103, lambda p, m, x, a, n, b : b**S(2)*n**S(2)*p**S(2)*Int(x**(m + S(2)*n)*sinh(a + b*x**n)**p, x)/((m + S(1))*(m + n + S(1))) + b**S(2)*n**S(2)*p*(p + S(-1))*Int(x**(m + S(2)*n)*sinh(a + b*x**n)**(p + S(-2)), x)/((m + S(1))*(m + n + S(1))) - b*n*p*x**(m + n + S(1))*Cosh(a + b*x**n)*sinh(a + b*x**n)**(p + S(-1))/((m + S(1))*(m + n + S(1))) + x**(m + S(1))*sinh(a + b*x**n)**p/(m + S(1)))
    rubi.add(rule103)

    pattern104 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda n, m: Less(S(0), S(2)*n, -m + S(1))), CustomConstraint(lambda n, m: NonzeroQ(m + n + S(1))))
    rule104 = ReplacementRule(pattern104, lambda p, m, x, a, n, b : b**S(2)*n**S(2)*p**S(2)*Int(x**(m + S(2)*n)*Cosh(a + b*x**n)**p, x)/((m + S(1))*(m + n + S(1))) - b**S(2)*n**S(2)*p*(p + S(-1))*Int(x**(m + S(2)*n)*Cosh(a + b*x**n)**(p + S(-2)), x)/((m + S(1))*(m + n + S(1))) - b*n*p*x**(m + n + S(1))*Cosh(a + b*x**n)**(p + S(-1))*sinh(a + b*x**n)/((m + S(1))*(m + n + S(1))) + x**(m + S(1))*Cosh(a + b*x**n)**p/(m + S(1)))
    rubi.add(rule104)

    pattern105 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)), )
    def With105(p, m, e, x, a, d, n, b, c):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*sinh(c + d*e**(-n)*x**(k*n)))**p, x), x, (e*x)**(1/k))/e
    rule105 = ReplacementRule(pattern105, lambda p, m, e, x, a, d, n, b, c : With105(p, m, e, x, a, d, n, b, c))
    rubi.add(rule105)

    pattern106 = Pattern(Integral((x_*WC('e', S(1)))**m_*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)), )
    def With106(p, m, e, x, a, d, n, b, c):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*Cosh(c + d*e**(-n)*x**(k*n)))**p, x), x, (e*x)**(1/k))/e
    rule106 = ReplacementRule(pattern106, lambda p, m, e, x, a, d, n, b, c : With106(p, m, e, x, a, d, n, b, c))
    rubi.add(rule106)

    pattern107 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: Greater(p, S(1))))
    rule107 = ReplacementRule(pattern107, lambda p, m, e, x, a, d, n, b, c : Int(ExpandTrigReduce((e*x)**m, (a + b*sinh(c + d*x**n))**p, x), x))
    rubi.add(rule107)

    pattern108 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: Greater(p, S(1))))
    rule108 = ReplacementRule(pattern108, lambda p, m, e, x, a, d, n, b, c : Int(ExpandTrigReduce((e*x)**m, (a + b*Cosh(c + d*x**n))**p, x), x))
    rubi.add(rule108)

    pattern109 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(m - S(2)*n + S(1))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))))
    rule109 = ReplacementRule(pattern109, lambda p, m, x, a, n, b : -(p + S(2))*Int(x**m*sinh(a + b*x**n)**(p + S(2)), x)/(p + S(1)) + x**n*Cosh(a + b*x**n)*sinh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))) - sinh(a + b*x**n)**(p + S(2))/(b**S(2)*n*(p + S(1))*(p + S(2))))
    rubi.add(rule109)

    pattern110 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(m - S(2)*n + S(1))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))))
    rule110 = ReplacementRule(pattern110, lambda p, m, x, a, n, b : (p + S(2))*Int(x**m*Cosh(a + b*x**n)**(p + S(2)), x)/(p + S(1)) - x**n*Cosh(a + b*x**n)**(p + S(1))*sinh(a + b*x**n)/(b*n*(p + S(1))) + Cosh(a + b*x**n)**(p + S(2))/(b**S(2)*n*(p + S(1))*(p + S(2))))
    rubi.add(rule110)

    pattern111 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))), CustomConstraint(lambda n, m: Less(S(0), S(2)*n, m + S(1))))
    rule111 = ReplacementRule(pattern111, lambda p, m, x, a, n, b : -(p + S(2))*Int(x**m*sinh(a + b*x**n)**(p + S(2)), x)/(p + S(1)) + x**(m - n + S(1))*Cosh(a + b*x**n)*sinh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))) - x**(m - S(2)*n + S(1))*(m - n + S(1))*sinh(a + b*x**n)**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))) + (m - S(2)*n + S(1))*(m - n + S(1))*Int(x**(m - S(2)*n)*sinh(a + b*x**n)**(p + S(2)), x)/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))))
    rubi.add(rule111)

    pattern112 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))), CustomConstraint(lambda n, m: Less(S(0), S(2)*n, m + S(1))))
    rule112 = ReplacementRule(pattern112, lambda p, m, x, a, n, b : (p + S(2))*Int(x**m*Cosh(a + b*x**n)**(p + S(2)), x)/(p + S(1)) - x**(m - n + S(1))*Cosh(a + b*x**n)**(p + S(1))*sinh(a + b*x**n)/(b*n*(p + S(1))) + x**(m - S(2)*n + S(1))*(m - n + S(1))*Cosh(a + b*x**n)**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))) - (m - S(2)*n + S(1))*(m - n + S(1))*Int(x**(m - S(2)*n)*Cosh(a + b*x**n)**(p + S(2)), x)/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))))
    rubi.add(rule112)

    pattern113 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule113 = ReplacementRule(pattern113, lambda p, m, x, a, d, n, b, c : -Subst(Int(x**(-m + S(-2))*(a + b*sinh(c + d*x**(-n)))**p, x), x, 1/x))
    rubi.add(rule113)

    pattern114 = Pattern(Integral(x_**WC('m', S(1))*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule114 = ReplacementRule(pattern114, lambda p, m, x, a, d, n, b, c : -Subst(Int(x**(-m + S(-2))*(a + b*Cosh(c + d*x**(-n)))**p, x), x, 1/x))
    rubi.add(rule114)

    pattern115 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)), )
    def With115(p, m, e, x, a, d, n, b, c):
        k = Denominator(m)
        return -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*sinh(c + d*e**(-n)*x**(-k*n)))**p, x), x, (e*x)**(-S(1)/k))/e
    rule115 = ReplacementRule(pattern115, lambda p, m, e, x, a, d, n, b, c : With115(p, m, e, x, a, d, n, b, c))
    rubi.add(rule115)

    pattern116 = Pattern(Integral((x_*WC('e', S(1)))**m_*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)), )
    def With116(p, m, e, x, a, d, n, b, c):
        k = Denominator(m)
        return -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*Cosh(c + d*e**(-n)*x**(-k*n)))**p, x), x, (e*x)**(-S(1)/k))/e
    rule116 = ReplacementRule(pattern116, lambda p, m, e, x, a, d, n, b, c : With116(p, m, e, x, a, d, n, b, c))
    rubi.add(rule116)

    pattern117 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: Not(RationalQ(m))))
    rule117 = ReplacementRule(pattern117, lambda p, m, e, x, a, d, n, b, c : -(e*x)**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*sinh(c + d*x**(-n)))**p, x), x, 1/x))
    rubi.add(rule117)

    pattern118 = Pattern(Integral((x_*WC('e', S(1)))**m_*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: Not(RationalQ(m))))
    rule118 = ReplacementRule(pattern118, lambda p, m, e, x, a, d, n, b, c : -(e*x)**m*(1/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*Cosh(c + d*x**(-n)))**p, x), x, 1/x))
    rubi.add(rule118)

    pattern119 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: FractionQ(n)), )
    def With119(p, m, x, a, d, n, b, c):
        k = Denominator(n)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*sinh(c + d*x**(k*n)))**p, x), x, x**(1/k))
    rule119 = ReplacementRule(pattern119, lambda p, m, x, a, d, n, b, c : With119(p, m, x, a, d, n, b, c))
    rubi.add(rule119)

    pattern120 = Pattern(Integral(x_**WC('m', S(1))*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: FractionQ(n)), )
    def With120(p, m, x, a, d, n, b, c):
        k = Denominator(n)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*Cosh(c + d*x**(k*n)))**p, x), x, x**(1/k))
    rule120 = ReplacementRule(pattern120, lambda p, m, x, a, d, n, b, c : With120(p, m, x, a, d, n, b, c))
    rubi.add(rule120)

    pattern121 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: FractionQ(n)))
    rule121 = ReplacementRule(pattern121, lambda p, m, e, x, a, d, n, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*sinh(c + d*x**n))**p, x))
    rubi.add(rule121)

    pattern122 = Pattern(Integral((e_*x_)**m_*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: FractionQ(n)))
    rule122 = ReplacementRule(pattern122, lambda p, m, e, x, a, d, n, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*Cosh(c + d*x**n))**p, x))
    rubi.add(rule122)

    pattern123 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: PositiveIntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule123 = ReplacementRule(pattern123, lambda p, m, x, a, d, n, b, c : Subst(Int((a + b*sinh(c + d*x**(n/(m + S(1)))))**p, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule123)

    pattern124 = Pattern(Integral(x_**WC('m', S(1))*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: PositiveIntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule124 = ReplacementRule(pattern124, lambda p, m, x, a, d, n, b, c : Subst(Int((a + b*Cosh(c + d*x**(n/(m + S(1)))))**p, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule124)

    pattern125 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: PositiveIntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule125 = ReplacementRule(pattern125, lambda p, m, e, x, a, d, n, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*sinh(c + d*x**n))**p, x))
    rubi.add(rule125)

    pattern126 = Pattern(Integral((e_*x_)**m_*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: PositiveIntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule126 = ReplacementRule(pattern126, lambda p, m, e, x, a, d, n, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*Cosh(c + d*x**n))**p, x))
    rubi.add(rule126)

    pattern127 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule127 = ReplacementRule(pattern127, lambda m, e, x, d, n, c : -Int((e*x)**m*exp(-c - d*x**n), x)/S(2) + Int((e*x)**m*exp(c + d*x**n), x)/S(2))
    rubi.add(rule127)

    pattern128 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule128 = ReplacementRule(pattern128, lambda m, e, x, d, n, c : Int((e*x)**m*exp(-c - d*x**n), x)/S(2) + Int((e*x)**m*exp(c + d*x**n), x)/S(2))
    rubi.add(rule128)

    pattern129 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule129 = ReplacementRule(pattern129, lambda p, m, e, x, a, d, n, b, c : Int(ExpandTrigReduce((e*x)**m, (a + b*sinh(c + d*x**n))**p, x), x))
    rubi.add(rule129)

    pattern130 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(Cosh(x_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule130 = ReplacementRule(pattern130, lambda p, m, e, x, a, d, n, b, c : Int(ExpandTrigReduce((e*x)**m, (a + b*Cosh(c + d*x**n))**p, x), x))
    rubi.add(rule130)

    pattern131 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: LinearQ(u, x)), CustomConstraint(lambda u, x: NonzeroQ(u - x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule131 = ReplacementRule(pattern131, lambda p, m, x, a, d, n, b, u, c : Coefficient(u, x, S(1))**(-m + S(-1))*Subst(Int((a + b*sinh(c + d*x**n))**p*(x - Coefficient(u, x, S(0)))**m, x), x, u))
    rubi.add(rule131)

    pattern132 = Pattern(Integral(x_**WC('m', S(1))*(Cosh(u_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: LinearQ(u, x)), CustomConstraint(lambda u, x: NonzeroQ(u - x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule132 = ReplacementRule(pattern132, lambda p, m, x, a, d, n, b, u, c : Coefficient(u, x, S(1))**(-m + S(-1))*Subst(Int((a + b*Cosh(c + d*x**n))**p*(x - Coefficient(u, x, S(0)))**m, x), x, u))
    rubi.add(rule132)

    pattern133 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: LinearQ(u, x)))
    rule133 = ReplacementRule(pattern133, lambda p, m, e, x, a, d, n, b, u, c : Int((e*x)**m*(a + b*sinh(c + d*u**n))**p, x))
    rubi.add(rule133)

    pattern134 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(Cosh(u_**n_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: LinearQ(u, x)))
    rule134 = ReplacementRule(pattern134, lambda p, m, e, x, a, d, n, b, u, c : Int((e*x)**m*(a + b*Cosh(c + d*u**n))**p, x))
    rubi.add(rule134)

    pattern135 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(u_))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule135 = ReplacementRule(pattern135, lambda p, m, e, x, a, u, b : Int((e*x)**m*(a + b*sinh(ExpandToSum(u, x)))**p, x))
    rubi.add(rule135)

    pattern136 = Pattern(Integral((e_*x_)**WC('m', S(1))*(Cosh(u_)*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule136 = ReplacementRule(pattern136, lambda p, m, e, x, a, u, b : Int((e*x)**m*(a + b*Cosh(ExpandToSum(u, x)))**p, x))
    rubi.add(rule136)

    pattern137 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule137 = ReplacementRule(pattern137, lambda p, m, x, a, n, b : sinh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))))
    rubi.add(rule137)

    pattern138 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule138 = ReplacementRule(pattern138, lambda p, m, x, a, n, b : Cosh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))))
    rubi.add(rule138)

    pattern139 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n, m: Less(S(0), n, m + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule139 = ReplacementRule(pattern139, lambda p, m, x, a, n, b : x**(m - n + S(1))*sinh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))) - (m - n + S(1))*Int(x**(m - n)*sinh(a + b*x**n)**(p + S(1)), x)/(b*n*(p + S(1))))
    rubi.add(rule139)

    pattern140 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n, m: Less(S(0), n, m + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule140 = ReplacementRule(pattern140, lambda p, m, x, a, n, b : x**(m - n + S(1))*Cosh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))) - (m - n + S(1))*Int(x**(m - n)*Cosh(a + b*x**n)**(p + S(1)), x)/(b*n*(p + S(1))))
    rubi.add(rule140)

    pattern141 = Pattern(Integral(sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule141 = ReplacementRule(pattern141, lambda c, b, x, a : -Int(exp(-a - b*x - c*x**S(2)), x)/S(2) + Int(exp(a + b*x + c*x**S(2)), x)/S(2))
    rubi.add(rule141)

    pattern142 = Pattern(Integral(Cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)))
    rule142 = ReplacementRule(pattern142, lambda c, b, x, a : Int(exp(-a - b*x - c*x**S(2)), x)/S(2) + Int(exp(a + b*x + c*x**S(2)), x)/S(2))
    rubi.add(rule142)

    pattern143 = Pattern(Integral(sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule143 = ReplacementRule(pattern143, lambda x, a, n, b, c : Int(ExpandTrigReduce(sinh(a + b*x + c*x**S(2))**n, x), x))
    rubi.add(rule143)

    pattern144 = Pattern(Integral(Cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule144 = ReplacementRule(pattern144, lambda x, a, n, b, c : Int(ExpandTrigReduce(Cosh(a + b*x + c*x**S(2))**n, x), x))
    rubi.add(rule144)

    pattern145 = Pattern(Integral(sinh(v_)**WC('n', S(1)), x_), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda v, x: Not(QuadraticMatchQ(v, x))))
    rule145 = ReplacementRule(pattern145, lambda n, v, x : Int(sinh(ExpandToSum(v, x))**n, x))
    rubi.add(rule145)

    pattern146 = Pattern(Integral(Cosh(v_)**WC('n', S(1)), x_), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda v, x: Not(QuadraticMatchQ(v, x))))
    rule146 = ReplacementRule(pattern146, lambda n, v, x : Int(Cosh(ExpandToSum(v, x))**n, x))
    rubi.add(rule146)

    pattern147 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e, b: ZeroQ(b*e - S(2)*c*d)))
    rule147 = ReplacementRule(pattern147, lambda e, x, a, d, b, c : e*Cosh(a + b*x + c*x**S(2))/(S(2)*c))
    rubi.add(rule147)

    pattern148 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*Cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e, b: ZeroQ(b*e - S(2)*c*d)))
    rule148 = ReplacementRule(pattern148, lambda e, x, a, d, b, c : e*sinh(a + b*x + c*x**S(2))/(S(2)*c))
    rubi.add(rule148)

    pattern149 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e, b: NonzeroQ(b*e - S(2)*c*d)))
    rule149 = ReplacementRule(pattern149, lambda e, x, a, d, b, c : e*Cosh(a + b*x + c*x**S(2))/(S(2)*c) - (b*e - S(2)*c*d)*Int(sinh(a + b*x + c*x**S(2)), x)/(S(2)*c))
    rubi.add(rule149)

    pattern150 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*Cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda d, c, e, b: NonzeroQ(b*e - S(2)*c*d)))
    rule150 = ReplacementRule(pattern150, lambda e, x, a, d, b, c : e*sinh(a + b*x + c*x**S(2))/(S(2)*c) - (b*e - S(2)*c*d)*Int(Cosh(a + b*x + c*x**S(2)), x)/(S(2)*c))
    rubi.add(rule150)

    pattern151 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda d, c, e, b: ZeroQ(b*e - S(2)*c*d)))
    rule151 = ReplacementRule(pattern151, lambda m, e, x, a, d, b, c : -e**S(2)*(m + S(-1))*Int((d + e*x)**(m + S(-2))*Cosh(a + b*x + c*x**S(2)), x)/(S(2)*c) + e*(d + e*x)**(m + S(-1))*Cosh(a + b*x + c*x**S(2))/(S(2)*c))
    rubi.add(rule151)

    pattern152 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*Cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda d, c, e, b: ZeroQ(b*e - S(2)*c*d)))
    rule152 = ReplacementRule(pattern152, lambda m, e, x, a, d, b, c : -e**S(2)*(m + S(-1))*Int((d + e*x)**(m + S(-2))*sinh(a + b*x + c*x**S(2)), x)/(S(2)*c) + e*(d + e*x)**(m + S(-1))*sinh(a + b*x + c*x**S(2))/(S(2)*c))
    rubi.add(rule152)

    pattern153 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda d, c, e, b: NonzeroQ(b*e - S(2)*c*d)))
    rule153 = ReplacementRule(pattern153, lambda m, e, x, a, d, b, c : -e**S(2)*(m + S(-1))*Int((d + e*x)**(m + S(-2))*Cosh(a + b*x + c*x**S(2)), x)/(S(2)*c) + e*(d + e*x)**(m + S(-1))*Cosh(a + b*x + c*x**S(2))/(S(2)*c) - (b*e - S(2)*c*d)*Int((d + e*x)**(m + S(-1))*sinh(a + b*x + c*x**S(2)), x)/(S(2)*c))
    rubi.add(rule153)

    pattern154 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*Cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(1))), CustomConstraint(lambda d, c, e, b: NonzeroQ(b*e - S(2)*c*d)))
    rule154 = ReplacementRule(pattern154, lambda m, e, x, a, d, b, c : -e**S(2)*(m + S(-1))*Int((d + e*x)**(m + S(-2))*sinh(a + b*x + c*x**S(2)), x)/(S(2)*c) + e*(d + e*x)**(m + S(-1))*sinh(a + b*x + c*x**S(2))/(S(2)*c) - (b*e - S(2)*c*d)*Int((d + e*x)**(m + S(-1))*Cosh(a + b*x + c*x**S(2)), x)/(S(2)*c))
    rubi.add(rule154)

    pattern155 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda d, c, e, b: ZeroQ(b*e - S(2)*c*d)))
    rule155 = ReplacementRule(pattern155, lambda m, e, x, a, d, b, c : -S(2)*c*Int((d + e*x)**(m + S(2))*Cosh(a + b*x + c*x**S(2)), x)/(e**S(2)*(m + S(1))) + (d + e*x)**(m + S(1))*sinh(a + b*x + c*x**S(2))/(e*(m + S(1))))
    rubi.add(rule155)

    pattern156 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*Cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda d, c, e, b: ZeroQ(b*e - S(2)*c*d)))
    rule156 = ReplacementRule(pattern156, lambda m, e, x, a, d, b, c : -S(2)*c*Int((d + e*x)**(m + S(2))*sinh(a + b*x + c*x**S(2)), x)/(e**S(2)*(m + S(1))) + (d + e*x)**(m + S(1))*Cosh(a + b*x + c*x**S(2))/(e*(m + S(1))))
    rubi.add(rule156)

    pattern157 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda d, c, e, b: NonzeroQ(b*e - S(2)*c*d)))
    rule157 = ReplacementRule(pattern157, lambda m, e, x, a, d, b, c : -S(2)*c*Int((d + e*x)**(m + S(2))*Cosh(a + b*x + c*x**S(2)), x)/(e**S(2)*(m + S(1))) + (d + e*x)**(m + S(1))*sinh(a + b*x + c*x**S(2))/(e*(m + S(1))) - (b*e - S(2)*c*d)*Int((d + e*x)**(m + S(1))*Cosh(a + b*x + c*x**S(2)), x)/(e**S(2)*(m + S(1))))
    rubi.add(rule157)

    pattern158 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*Cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda d, c, e, b: NonzeroQ(b*e - S(2)*c*d)))
    rule158 = ReplacementRule(pattern158, lambda m, e, x, a, d, b, c : -S(2)*c*Int((d + e*x)**(m + S(2))*sinh(a + b*x + c*x**S(2)), x)/(e**S(2)*(m + S(1))) + (d + e*x)**(m + S(1))*Cosh(a + b*x + c*x**S(2))/(e*(m + S(1))) - (b*e - S(2)*c*d)*Int((d + e*x)**(m + S(1))*sinh(a + b*x + c*x**S(2)), x)/(e**S(2)*(m + S(1))))
    rubi.add(rule158)

    pattern159 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule159 = ReplacementRule(pattern159, lambda m, c, e, x, a, d, b : Int((d + e*x)**m*sinh(a + b*x + c*x**S(2)), x))
    rubi.add(rule159)

    pattern160 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*Cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)))
    rule160 = ReplacementRule(pattern160, lambda m, e, x, a, d, b, c : Int((d + e*x)**m*Cosh(a + b*x + c*x**S(2)), x))
    rubi.add(rule160)

    pattern161 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule161 = ReplacementRule(pattern161, lambda m, c, e, x, a, d, n, b : Int(ExpandTrigReduce((d + e*x)**m, sinh(a + b*x + c*x**S(2))**n, x), x))
    rubi.add(rule161)

    pattern162 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*Cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule162 = ReplacementRule(pattern162, lambda m, e, x, a, d, n, b, c : Int(ExpandTrigReduce((d + e*x)**m, Cosh(a + b*x + c*x**S(2))**n, x), x))
    rubi.add(rule162)

    pattern163 = Pattern(Integral(u_**WC('m', S(1))*sinh(v_)**WC('n', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda u, x: LinearQ(u, x)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda v, u, x: Not(LinearMatchQ(u, x) & QuadraticMatchQ(v, x))))
    rule163 = ReplacementRule(pattern163, lambda m, x, n, v, u : Int(ExpandToSum(u, x)**m*sinh(ExpandToSum(v, x))**n, x))
    rubi.add(rule163)

    pattern164 = Pattern(Integral(u_**WC('m', S(1))*Cosh(v_)**WC('n', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda u, x: LinearQ(u, x)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda v, u, x: Not(LinearMatchQ(u, x) & QuadraticMatchQ(v, x))))
    rule164 = ReplacementRule(pattern164, lambda m, x, n, v, u : Int(Cosh(ExpandToSum(v, x))**n*ExpandToSum(u, x)**m, x))
    rubi.add(rule164)

    pattern165 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule165 = ReplacementRule(pattern165, lambda m, e, x, f, a, b : S(2)*Int((a + b*x)**m*exp(S(2)*e + S(2)*f*x)/(exp(S(2)*e + S(2)*f*x) + S(1)), x) - (a + b*x)**(m + S(1))/(b*(m + S(1))))
    rubi.add(rule165)

    pattern166 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule166 = ReplacementRule(pattern166, lambda m, e, x, a, f, b : -S(2)*Int((a + b*x)**m*exp(S(2)*e + S(2)*f*x)/(-exp(S(2)*e + S(2)*f*x) + S(1)), x) - (a + b*x)**(m + S(1))/(b*(m + S(1))))
    rubi.add(rule166)

    pattern167 = Pattern(Integral((WC('c', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: Greater(m, S(0))))
    rule167 = ReplacementRule(pattern167, lambda m, e, x, a, f, n, b, c : b*c*m*Int((c*tanh(e + f*x))**(n + S(-1))*(a + b*x)**(m + S(-1)), x)/(f*(n + S(-1))) + c**S(2)*Int((c*tanh(e + f*x))**(n + S(-2))*(a + b*x)**m, x) - c*(c*tanh(e + f*x))**(n + S(-1))*(a + b*x)**m/(f*(n + S(-1))))
    rubi.add(rule167)

    pattern168 = Pattern(Integral((WC('c', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda m: Greater(m, S(0))))
    rule168 = ReplacementRule(pattern168, lambda m, c, e, x, a, f, n, b : b*c*m*Int((c*coth(e + f*x))**(n + S(-1))*(a + b*x)**(m + S(-1)), x)/(f*(n + S(-1))) + c**S(2)*Int((c*coth(e + f*x))**(n + S(-2))*(a + b*x)**m, x) - c*(c*coth(e + f*x))**(n + S(-1))*(a + b*x)**m/(f*(n + S(-1))))
    rubi.add(rule168)

    pattern169 = Pattern(Integral((WC('c', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda m: Greater(m, S(0))))
    rule169 = ReplacementRule(pattern169, lambda m, e, x, a, f, n, b, c : -b*m*Int((c*tanh(e + f*x))**(n + S(1))*(a + b*x)**(m + S(-1)), x)/(c*f*(n + S(1))) + (c*tanh(e + f*x))**(n + S(1))*(a + b*x)**m/(c*f*(n + S(1))) + Int((c*tanh(e + f*x))**(n + S(2))*(a + b*x)**m, x)/c**S(2))
    rubi.add(rule169)

    pattern170 = Pattern(Integral((WC('c', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda m: Greater(m, S(0))))
    rule170 = ReplacementRule(pattern170, lambda m, c, e, x, a, f, n, b : -b*m*Int((c*coth(e + f*x))**(n + S(1))*(a + b*x)**(m + S(-1)), x)/(c*f*(n + S(1))) + (c*coth(e + f*x))**(n + S(1))*(a + b*x)**m/(c*f*(n + S(1))) + Int((c*coth(e + f*x))**(n + S(2))*(a + b*x)**m, x)/c**S(2))
    rubi.add(rule170)

    pattern171 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule171 = ReplacementRule(pattern171, lambda m, c, e, x, f, d, a, n, b : Int(ExpandIntegrand((c + d*x)**m, (a + b*tanh(e + f*x))**n, x), x))
    rubi.add(rule171)

    pattern172 = Pattern(Integral((a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule172 = ReplacementRule(pattern172, lambda m, e, x, f, d, a, n, b, c : Int(ExpandIntegrand((c + d*x)**m, (a + b*coth(e + f*x))**n, x), x))
    rubi.add(rule172)

    pattern173 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule173 = ReplacementRule(pattern173, lambda m, c, e, x, f, d, a, b : a*d*m*Int((c + d*x)**(m + S(-1))/(a + b*tanh(e + f*x)), x)/(S(2)*b*f) - a*(c + d*x)**m/(S(2)*b*f*(a + b*tanh(e + f*x))) + (c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))))
    rubi.add(rule173)

    pattern174 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule174 = ReplacementRule(pattern174, lambda m, e, x, f, d, a, b, c : a*d*m*Int((c + d*x)**(m + S(-1))/(a + b*coth(e + f*x)), x)/(S(2)*b*f) - a*(c + d*x)**m/(S(2)*b*f*(a + b*coth(e + f*x))) + (c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))))
    rubi.add(rule174)

    pattern175 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))))
    rule175 = ReplacementRule(pattern175, lambda c, e, x, f, d, a, b : -S(1)/(d*(a + b*tanh(e + f*x))*(c + d*x)) - f*Int(Cosh(S(2)*e + S(2)*f*x)/(c + d*x), x)/(b*d) + f*Int(sinh(S(2)*e + S(2)*f*x)/(c + d*x), x)/(a*d))
    rubi.add(rule175)

    pattern176 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))))
    rule176 = ReplacementRule(pattern176, lambda e, x, f, d, a, b, c : -S(1)/(d*(a + b*coth(e + f*x))*(c + d*x)) + f*Int(Cosh(S(2)*e + S(2)*f*x)/(c + d*x), x)/(b*d) - f*Int(sinh(S(2)*e + S(2)*f*x)/(c + d*x), x)/(a*d))
    rubi.add(rule176)

    pattern177 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: Unequal(m, S(-2))))
    rule177 = ReplacementRule(pattern177, lambda m, c, e, x, f, d, a, b : (c + d*x)**(m + S(1))/(d*(a + b*tanh(e + f*x))*(m + S(1))) - f*(c + d*x)**(m + S(2))/(b*d**S(2)*(m + S(1))*(m + S(2))) + S(2)*b*f*Int((c + d*x)**(m + S(1))/(a + b*tanh(e + f*x)), x)/(a*d*(m + S(1))))
    rubi.add(rule177)

    pattern178 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda m: Unequal(m, S(-2))))
    rule178 = ReplacementRule(pattern178, lambda m, e, x, f, d, a, b, c : (c + d*x)**(m + S(1))/(d*(a + b*coth(e + f*x))*(m + S(1))) - f*(c + d*x)**(m + S(2))/(b*d**S(2)*(m + S(1))*(m + S(2))) + S(2)*b*f*Int((c + d*x)**(m + S(1))/(a + b*coth(e + f*x)), x)/(a*d*(m + S(1))))
    rubi.add(rule178)

    pattern179 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))))
    rule179 = ReplacementRule(pattern179, lambda c, e, x, f, d, a, b : -Int(sinh(S(2)*e + S(2)*f*x)/(c + d*x), x)/(S(2)*b) + Int(Cosh(S(2)*e + S(2)*f*x)/(c + d*x), x)/(S(2)*a) + log(c + d*x)/(S(2)*a*d))
    rubi.add(rule179)

    pattern180 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))))
    rule180 = ReplacementRule(pattern180, lambda e, x, f, d, a, b, c : Int(sinh(S(2)*e + S(2)*f*x)/(c + d*x), x)/(S(2)*b) - Int(Cosh(S(2)*e + S(2)*f*x)/(c + d*x), x)/(S(2)*a) + log(c + d*x)/(S(2)*a*d))
    rubi.add(rule180)

    pattern181 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule181 = ReplacementRule(pattern181, lambda m, c, e, x, f, d, a, b : Int((c + d*x)**m*exp(-S(2)*a*(e + f*x)/b), x)/(S(2)*a) + (c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))))
    rubi.add(rule181)

    pattern182 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule182 = ReplacementRule(pattern182, lambda m, e, x, f, d, a, b, c : -Int((c + d*x)**m*exp(-S(2)*a*(e + f*x)/b), x)/(S(2)*a) + (c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))))
    rubi.add(rule182)

    pattern183 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda n, m: NegativeIntegerQ(m, n)))
    rule183 = ReplacementRule(pattern183, lambda m, c, e, x, f, d, a, n, b : Int(ExpandIntegrand((c + d*x)**m, (-sinh(S(2)*e + S(2)*f*x)/(S(2)*b) + Cosh(S(2)*e + S(2)*f*x)/(S(2)*a) + S(1)/(S(2)*a))**(-n), x), x))
    rubi.add(rule183)

    pattern184 = Pattern(Integral((a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda n, m: NegativeIntegerQ(m, n)))
    rule184 = ReplacementRule(pattern184, lambda m, e, x, f, d, a, n, b, c : Int(ExpandIntegrand((c + d*x)**m, (sinh(S(2)*e + S(2)*f*x)/(S(2)*b) - Cosh(S(2)*e + S(2)*f*x)/(S(2)*a) + S(1)/(S(2)*a))**(-n), x), x))
    rubi.add(rule184)

    pattern185 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule185 = ReplacementRule(pattern185, lambda m, c, e, x, f, d, a, n, b : Int(ExpandIntegrand((c + d*x)**m, (S(1)/(S(2)*a) + exp(-S(2)*a*(e + f*x)/b)/(S(2)*a))**(-n), x), x))
    rubi.add(rule185)

    pattern186 = Pattern(Integral((a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule186 = ReplacementRule(pattern186, lambda m, e, x, f, d, a, n, b, c : Int(ExpandIntegrand((c + d*x)**m, (S(1)/(S(2)*a) - exp(-S(2)*a*(e + f*x)/b)/(S(2)*a))**(-n), x), x))
    rubi.add(rule186)

    pattern187 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n + S(1))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))), )
    def With187(m, c, e, x, f, d, a, n, b):
        u = IntHide((a + b*tanh(e + f*x))**n, x)
        return -d*m*Int(Dist((c + d*x)**(m + S(-1)), u, x), x) + Dist((c + d*x)**m, u, x)
    rule187 = ReplacementRule(pattern187, lambda m, c, e, x, f, d, a, n, b : With187(m, c, e, x, f, d, a, n, b))
    rubi.add(rule187)

    pattern188 = Pattern(Integral((a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n + S(1))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))), )
    def With188(m, e, x, f, d, a, n, b, c):
        u = IntHide((a + b*coth(e + f*x))**n, x)
        return -d*m*Int(Dist((c + d*x)**(m + S(-1)), u, x), x) + Dist((c + d*x)**m, u, x)
    rule188 = ReplacementRule(pattern188, lambda m, e, x, f, d, a, n, b, c : With188(m, e, x, f, d, a, n, b, c))
    rubi.add(rule188)

    pattern189 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule189 = ReplacementRule(pattern189, lambda m, c, e, x, f, d, a, b : -S(2)*b*Int((c + d*x)**m/(a**S(2) - b**S(2) + (a - b)**S(2)*exp(-S(2)*e - S(2)*f*x)), x) + (c + d*x)**(m + S(1))/(d*(a - b)*(m + S(1))))
    rubi.add(rule189)

    pattern190 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule190 = ReplacementRule(pattern190, lambda m, e, x, f, d, a, b, c : S(2)*b*Int((c + d*x)**m/(a**S(2) - b**S(2) - (a + b)**S(2)*exp(S(2)*e + S(2)*f*x)), x) + (c + d*x)**(m + S(1))/(d*(a + b)*(m + S(1))))
    rubi.add(rule190)

    pattern191 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))))
    rule191 = ReplacementRule(pattern191, lambda c, e, x, f, d, a, b : b*(c + d*x)/(f*(a + b*tanh(e + f*x))*(a**S(2) - b**S(2))) - Int((-S(2)*a*c*f - S(2)*a*d*f*x + b*d)/(a + b*tanh(e + f*x)), x)/(f*(a**S(2) - b**S(2))) - (c + d*x)**S(2)/(S(2)*d*(a**S(2) - b**S(2))))
    rubi.add(rule191)

    pattern192 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))))
    rule192 = ReplacementRule(pattern192, lambda e, x, f, d, a, b, c : b*(c + d*x)/(f*(a + b*coth(e + f*x))*(a**S(2) - b**S(2))) - Int((-S(2)*a*c*f - S(2)*a*d*f*x + b*d)/(a + b*coth(e + f*x)), x)/(f*(a**S(2) - b**S(2))) - (c + d*x)**S(2)/(S(2)*d*(a**S(2) - b**S(2))))
    rubi.add(rule192)

    pattern193 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule193 = ReplacementRule(pattern193, lambda m, c, e, x, f, d, a, n, b : Int(ExpandIntegrand((c + d*x)**m, (-S(2)*b/(a**S(2) - b**S(2) + (a - b)**S(2)*exp(-S(2)*e - S(2)*f*x)) + 1/(a - b))**(-n), x), x))
    rubi.add(rule193)

    pattern194 = Pattern(Integral((a_ + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule194 = ReplacementRule(pattern194, lambda m, e, x, f, d, a, n, b, c : Int(ExpandIntegrand((c + d*x)**m, (S(2)*b/(a**S(2) - b**S(2) - (a + b)**S(2)*exp(S(2)*e + S(2)*f*x)) + 1/(a + b))**(-n), x), x))
    rubi.add(rule194)

    pattern195 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(v_))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda v, u, x: LinearQ(List(u, v), x)), CustomConstraint(lambda v, u, x: Not(LinearMatchQ(List(u, v), x))))
    rule195 = ReplacementRule(pattern195, lambda m, x, a, n, v, u, b : Int((a + b*tanh(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x))
    rubi.add(rule195)

    pattern196 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*coth(v_))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda v, u, x: LinearQ(List(u, v), x)), CustomConstraint(lambda v, u, x: Not(LinearMatchQ(List(u, v), x))))
    rule196 = ReplacementRule(pattern196, lambda m, x, a, n, v, u, b : Int((a + b*coth(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x))
    rubi.add(rule196)

    pattern197 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule197 = ReplacementRule(pattern197, lambda m, c, e, x, a, d, f, n, b : Int((a + b*tanh(e + f*x))**n*(c + d*x)**m, x))
    rubi.add(rule197)

    pattern198 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*coth(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule198 = ReplacementRule(pattern198, lambda m, c, e, x, f, a, d, n, b : Int((a + b*coth(e + f*x))**n*(c + d*x)**m, x))
    rubi.add(rule198)

    pattern199 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(1/n)), CustomConstraint(lambda p: IntegerQ(p)))
    rule199 = ReplacementRule(pattern199, lambda p, x, a, d, n, b, c : Subst(Int(x**(S(-1) + 1/n)*(a + b*tanh(c + d*x))**p, x), x, x**n)/n)
    rubi.add(rule199)

    pattern200 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*coth(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(1/n)), CustomConstraint(lambda p: IntegerQ(p)))
    rule200 = ReplacementRule(pattern200, lambda p, x, a, d, n, b, c : Subst(Int(x**(S(-1) + 1/n)*(a + b*coth(c + d*x))**p, x), x, x**n)/n)
    rubi.add(rule200)

    pattern201 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule201 = ReplacementRule(pattern201, lambda p, x, a, d, n, b, c : Int((a + b*tanh(c + d*x**n))**p, x))
    rubi.add(rule201)

    pattern202 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*coth(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule202 = ReplacementRule(pattern202, lambda p, x, a, d, n, b, c : Int((a + b*coth(c + d*x**n))**p, x))
    rubi.add(rule202)

    pattern203 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tanh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: LinearQ(u, x)), CustomConstraint(lambda u, x: NonzeroQ(u - x)))
    rule203 = ReplacementRule(pattern203, lambda p, x, a, d, n, b, u, c : Subst(Int((a + b*tanh(c + d*x**n))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule203)

    pattern204 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*coth(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: LinearQ(u, x)), CustomConstraint(lambda u, x: NonzeroQ(u - x)))
    rule204 = ReplacementRule(pattern204, lambda p, x, a, d, n, b, u, c : Subst(Int((a + b*coth(c + d*x**n))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule204)

    pattern205 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tanh(u_))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule205 = ReplacementRule(pattern205, lambda p, x, a, u, b : Int((a + b*tanh(ExpandToSum(u, x)))**p, x))
    rubi.add(rule205)

    pattern206 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*coth(u_))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule206 = ReplacementRule(pattern206, lambda p, x, a, u, b : Int((a + b*coth(ExpandToSum(u, x)))**p, x))
    rubi.add(rule206)

    pattern207 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: PositiveIntegerQ((m + S(1))/n)), CustomConstraint(lambda p: IntegerQ(p)))
    rule207 = ReplacementRule(pattern207, lambda p, m, x, a, d, n, b, c : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*tanh(c + d*x))**p, x), x, x**n)/n)
    rubi.add(rule207)

    pattern208 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*coth(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: PositiveIntegerQ((m + S(1))/n)), CustomConstraint(lambda p: IntegerQ(p)))
    rule208 = ReplacementRule(pattern208, lambda p, m, x, a, d, n, b, c : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*coth(c + d*x))**p, x), x, x**n)/n)
    rubi.add(rule208)

    pattern209 = Pattern(Integral(x_**WC('m', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule209 = ReplacementRule(pattern209, lambda m, x, d, n, c : Int(x**m, x) - x**(m - n + S(1))*tanh(c + d*x**n)/(d*n) + (m - n + S(1))*Int(x**(m - n)*tanh(c + d*x**n), x)/(d*n))
    rubi.add(rule209)

    pattern210 = Pattern(Integral(x_**WC('m', S(1))*coth(x_**n_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule210 = ReplacementRule(pattern210, lambda m, x, d, n, c : Int(x**m, x) - x**(m - n + S(1))*coth(c + d*x**n)/(d*n) + (m - n + S(1))*Int(x**(m - n)*coth(c + d*x**n), x)/(d*n))
    rubi.add(rule210)

    pattern211 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule211 = ReplacementRule(pattern211, lambda p, m, x, a, d, n, b, c : Int(x**m*(a + b*tanh(c + d*x**n))**p, x))
    rubi.add(rule211)

    pattern212 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*coth(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule212 = ReplacementRule(pattern212, lambda p, m, x, a, d, n, b, c : Int(x**m*(a + b*coth(c + d*x**n))**p, x))
    rubi.add(rule212)

    pattern213 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule213 = ReplacementRule(pattern213, lambda p, m, e, x, a, d, n, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*tanh(c + d*x**n))**p, x))
    rubi.add(rule213)

    pattern214 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*coth(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule214 = ReplacementRule(pattern214, lambda p, m, e, x, a, d, n, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*coth(c + d*x**n))**p, x))
    rubi.add(rule214)

    pattern215 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(u_))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule215 = ReplacementRule(pattern215, lambda p, m, e, x, a, u, b : Int((e*x)**m*(a + b*tanh(ExpandToSum(u, x)))**p, x))
    rubi.add(rule215)

    pattern216 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*coth(u_))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule216 = ReplacementRule(pattern216, lambda p, m, e, x, a, u, b : Int((e*x)**m*(a + b*coth(ExpandToSum(u, x)))**p, x))
    rubi.add(rule216)

    pattern217 = Pattern(Integral(x_**WC('m', S(1))*tanh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('q', S(1))*sech(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n, m: GreaterEqual(m - n, S(0))), CustomConstraint(lambda q: SameQ(q, S(1))))
    rule217 = ReplacementRule(pattern217, lambda p, m, x, a, n, q, b : -x**(m - n + S(1))*sech(a + b*x**n)**p/(b*n*p) + (m - n + S(1))*Int(x**(m - n)*sech(a + b*x**n)**p, x)/(b*n*p))
    rubi.add(rule217)

    pattern218 = Pattern(Integral(x_**WC('m', S(1))*coth(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('q', S(1))*csch(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n, m: GreaterEqual(m - n, S(0))), CustomConstraint(lambda q: SameQ(q, S(1))))
    rule218 = ReplacementRule(pattern218, lambda p, m, x, a, n, q, b : -x**(m - n + S(1))*csch(a + b*x**n)**p/(b*n*p) + (m - n + S(1))*Int(x**(m - n)*csch(a + b*x**n)**p, x)/(b*n*p))
    rubi.add(rule218)

    pattern219 = Pattern(Integral(tanh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule219 = ReplacementRule(pattern219, lambda x, a, n, b, c : Int(tanh(a + b*x + c*x**S(2))**n, x))
    rubi.add(rule219)

    pattern220 = Pattern(Integral(coth(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule220 = ReplacementRule(pattern220, lambda x, a, n, b, c : Int(coth(a + b*x + c*x**S(2))**n, x))
    rubi.add(rule220)

    pattern221 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*tanh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule221 = ReplacementRule(pattern221, lambda e, x, a, d, b, c : e*log(Cosh(a + b*x + c*x**S(2)))/(S(2)*c) + (-b*e + S(2)*c*d)*Int(tanh(a + b*x + c*x**S(2)), x)/(S(2)*c))
    rubi.add(rule221)

    pattern222 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*coth(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule222 = ReplacementRule(pattern222, lambda e, x, a, d, b, c : e*log(sinh(a + b*x + c*x**S(2)))/(S(2)*c) + (-b*e + S(2)*c*d)*Int(coth(a + b*x + c*x**S(2)), x)/(S(2)*c))
    rubi.add(rule222)

    pattern223 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*tanh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule223 = ReplacementRule(pattern223, lambda m, c, e, x, a, d, n, b : Int((d + e*x)**m*tanh(a + b*x + c*x**S(2))**n, x))
    rubi.add(rule223)

    pattern224 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*coth(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule224 = ReplacementRule(pattern224, lambda m, e, x, a, d, n, b, c : Int((d + e*x)**m*coth(a + b*x + c*x**S(2))**n, x))
    rubi.add(rule224)

    pattern225 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sech(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule225 = ReplacementRule(pattern225, lambda m, c, x, a, d, b : -ImaginaryI*d*m*Int((c + d*x)**(m + S(-1))*log(-ImaginaryI*exp(a + b*x) + S(1)), x)/b + ImaginaryI*d*m*Int((c + d*x)**(m + S(-1))*log(ImaginaryI*exp(a + b*x) + S(1)), x)/b + S(2)*(c + d*x)**m*ArcTan(exp(a + b*x))/b)
    rubi.add(rule225)

    pattern226 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*csch(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule226 = ReplacementRule(pattern226, lambda m, x, a, d, b, c : -d*m*Int((c + d*x)**(m + S(-1))*log(-exp(a + b*x) + S(1)), x)/b + d*m*Int((c + d*x)**(m + S(-1))*log(exp(a + b*x) + S(1)), x)/b - S(2)*(c + d*x)**m*atanh(exp(a + b*x))/b)
    rubi.add(rule226)

    pattern227 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sech(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule227 = ReplacementRule(pattern227, lambda m, c, x, a, d, b : -d*m*Int((c + d*x)**(m + S(-1))*tanh(a + b*x), x)/b + (c + d*x)**m*tanh(a + b*x)/b)
    rubi.add(rule227)

    pattern228 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*csch(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule228 = ReplacementRule(pattern228, lambda m, x, a, d, b, c : d*m*Int((c + d*x)**(m + S(-1))*coth(a + b*x), x)/b - (c + d*x)**m*coth(a + b*x)/b)
    rubi.add(rule228)

    pattern229 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))*sech(x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda n: Unequal(n, S(2))))
    rule229 = ReplacementRule(pattern229, lambda c, x, a, d, n, b : (n + S(-2))*Int((c + d*x)*sech(a + b*x)**(n + S(-2)), x)/(n + S(-1)) + (c + d*x)*tanh(a + b*x)*sech(a + b*x)**(n + S(-2))/(b*(n + S(-1))) + d*sech(a + b*x)**(n + S(-2))/(b**S(2)*(n + S(-2))*(n + S(-1))))
    rubi.add(rule229)

    pattern230 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))*csch(x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda n: Unequal(n, S(2))))
    rule230 = ReplacementRule(pattern230, lambda x, a, d, n, b, c : -(n + S(-2))*Int((c + d*x)*csch(a + b*x)**(n + S(-2)), x)/(n + S(-1)) + (-c - d*x)*coth(a + b*x)*csch(a + b*x)**(n + S(-2))/(b*(n + S(-1))) - d*csch(a + b*x)**(n + S(-2))/(b**S(2)*(n + S(-2))*(n + S(-1))))
    rubi.add(rule230)

    pattern231 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*sech(x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda n: Unequal(n, S(2))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule231 = ReplacementRule(pattern231, lambda m, c, x, a, d, n, b : (n + S(-2))*Int((c + d*x)**m*sech(a + b*x)**(n + S(-2)), x)/(n + S(-1)) + (c + d*x)**m*tanh(a + b*x)*sech(a + b*x)**(n + S(-2))/(b*(n + S(-1))) - d**S(2)*m*(m + S(-1))*Int((c + d*x)**(m + S(-2))*sech(a + b*x)**(n + S(-2)), x)/(b**S(2)*(n + S(-2))*(n + S(-1))) + d*m*(c + d*x)**(m + S(-1))*sech(a + b*x)**(n + S(-2))/(b**S(2)*(n + S(-2))*(n + S(-1))))
    rubi.add(rule231)

    pattern232 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*csch(x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda n: Unequal(n, S(2))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule232 = ReplacementRule(pattern232, lambda m, x, a, d, n, b, c : -(n + S(-2))*Int((c + d*x)**m*csch(a + b*x)**(n + S(-2)), x)/(n + S(-1)) - (c + d*x)**m*coth(a + b*x)*csch(a + b*x)**(n + S(-2))/(b*(n + S(-1))) + d**S(2)*m*(m + S(-1))*Int((c + d*x)**(m + S(-2))*csch(a + b*x)**(n + S(-2)), x)/(b**S(2)*(n + S(-2))*(n + S(-1))) - d*m*(c + d*x)**(m + S(-1))*csch(a + b*x)**(n + S(-2))/(b**S(2)*(n + S(-2))*(n + S(-1))))
    rubi.add(rule232)

    pattern233 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))*sech(x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule233 = ReplacementRule(pattern233, lambda c, x, a, d, n, b : (n + S(1))*Int((c + d*x)*sech(a + b*x)**(n + S(2)), x)/n - (c + d*x)*sinh(a + b*x)*sech(a + b*x)**(n + S(1))/(b*n) - d*sech(a + b*x)**n/(b**S(2)*n**S(2)))
    rubi.add(rule233)

    pattern234 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))*csch(x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule234 = ReplacementRule(pattern234, lambda x, a, d, n, b, c : -(n + S(1))*Int((c + d*x)*csch(a + b*x)**(n + S(2)), x)/n - (c + d*x)*Cosh(a + b*x)*csch(a + b*x)**(n + S(1))/(b*n) - d*csch(a + b*x)**n/(b**S(2)*n**S(2)))
    rubi.add(rule234)

    pattern235 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*sech(x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule235 = ReplacementRule(pattern235, lambda m, c, x, a, d, n, b : (n + S(1))*Int((c + d*x)**m*sech(a + b*x)**(n + S(2)), x)/n - (c + d*x)**m*sinh(a + b*x)*sech(a + b*x)**(n + S(1))/(b*n) + d**S(2)*m*(m + S(-1))*Int((c + d*x)**(m + S(-2))*sech(a + b*x)**n, x)/(b**S(2)*n**S(2)) - d*m*(c + d*x)**(m + S(-1))*sech(a + b*x)**n/(b**S(2)*n**S(2)))
    rubi.add(rule235)

    pattern236 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*csch(x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda m: Greater(m, S(1))))
    rule236 = ReplacementRule(pattern236, lambda m, x, a, d, n, b, c : -(n + S(1))*Int((c + d*x)**m*csch(a + b*x)**(n + S(2)), x)/n - (c + d*x)**m*Cosh(a + b*x)*csch(a + b*x)**(n + S(1))/(b*n) + d**S(2)*m*(m + S(-1))*Int((c + d*x)**(m + S(-2))*csch(a + b*x)**n, x)/(b**S(2)*n**S(2)) - d*m*(c + d*x)**(m + S(-1))*csch(a + b*x)**n/(b**S(2)*n**S(2)))
    rubi.add(rule236)

    pattern237 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sech(x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule237 = ReplacementRule(pattern237, lambda m, c, x, a, d, n, b : Cosh(a + b*x)**n*Int((c + d*x)**m*Cosh(a + b*x)**(-n), x)*sech(a + b*x)**n)
    rubi.add(rule237)

    pattern238 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*csch(x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule238 = ReplacementRule(pattern238, lambda m, x, a, d, n, b, c : Int((c + d*x)**m*sinh(a + b*x)**(-n), x)*sinh(a + b*x)**n*csch(a + b*x)**n)
    rubi.add(rule238)

    pattern239 = Pattern(Integral((a_ + WC('b', S(1))*sech(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule239 = ReplacementRule(pattern239, lambda m, c, e, x, f, d, a, n, b : Int(ExpandIntegrand((c + d*x)**m, (a + b*sech(e + f*x))**n, x), x))
    rubi.add(rule239)

    pattern240 = Pattern(Integral((a_ + WC('b', S(1))*csch(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule240 = ReplacementRule(pattern240, lambda m, e, x, f, d, a, n, b, c : Int(ExpandIntegrand((c + d*x)**m, (a + b*csch(e + f*x))**n, x), x))
    rubi.add(rule240)

    pattern241 = Pattern(Integral((a_ + WC('b', S(1))*sech(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule241 = ReplacementRule(pattern241, lambda m, c, e, x, f, d, a, n, b : Int(ExpandIntegrand((c + d*x)**m, (a*Cosh(e + f*x) + b)**n*Cosh(e + f*x)**(-n), x), x))
    rubi.add(rule241)

    pattern242 = Pattern(Integral((a_ + WC('b', S(1))*csch(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule242 = ReplacementRule(pattern242, lambda m, e, x, f, d, a, n, b, c : Int(ExpandIntegrand((c + d*x)**m, (a*sinh(e + f*x) + b)**n*sinh(e + f*x)**(-n), x), x))
    rubi.add(rule242)

    pattern243 = Pattern(Integral(u_**WC('m', S(1))*sech(v_)**WC('n', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda v, u, x: LinearQ(List(u, v), x)), CustomConstraint(lambda v, u, x: Not(LinearMatchQ(List(u, v), x))))
    rule243 = ReplacementRule(pattern243, lambda m, x, n, v, u : Int(ExpandToSum(u, x)**m*sech(ExpandToSum(v, x))**n, x))
    rubi.add(rule243)

    pattern244 = Pattern(Integral(u_**WC('m', S(1))*csch(v_)**WC('n', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda v, u, x: LinearQ(List(u, v), x)), CustomConstraint(lambda v, u, x: Not(LinearMatchQ(List(u, v), x))))
    rule244 = ReplacementRule(pattern244, lambda m, x, n, v, u : Int(ExpandToSum(u, x)**m*csch(ExpandToSum(v, x))**n, x))
    rubi.add(rule244)

    pattern245 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sech(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule245 = ReplacementRule(pattern245, lambda m, c, x, a, d, n, b : Int((c + d*x)**m*sech(a + b*x)**n, x))
    rubi.add(rule245)

    pattern246 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*csch(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule246 = ReplacementRule(pattern246, lambda m, x, a, d, n, b, c : Int((c + d*x)**m*csch(a + b*x)**n, x))
    rubi.add(rule246)

    pattern247 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sech(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(1/n)), CustomConstraint(lambda p: IntegerQ(p)))
    rule247 = ReplacementRule(pattern247, lambda p, x, a, d, n, b, c : Subst(Int(x**(S(-1) + 1/n)*(a + b*sech(c + d*x))**p, x), x, x**n)/n)
    rubi.add(rule247)

    pattern248 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*csch(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n: PositiveIntegerQ(1/n)), CustomConstraint(lambda p: IntegerQ(p)))
    rule248 = ReplacementRule(pattern248, lambda p, x, a, d, n, b, c : Subst(Int(x**(S(-1) + 1/n)*(a + b*csch(c + d*x))**p, x), x, x**n)/n)
    rubi.add(rule248)

    pattern249 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sech(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule249 = ReplacementRule(pattern249, lambda p, x, a, d, n, b, c : Int((a + b*sech(c + d*x**n))**p, x))
    rubi.add(rule249)

    pattern250 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*csch(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule250 = ReplacementRule(pattern250, lambda p, x, a, d, n, b, c : Int((a + b*csch(c + d*x**n))**p, x))
    rubi.add(rule250)

    pattern251 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sech(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: LinearQ(u, x)), CustomConstraint(lambda u, x: NonzeroQ(u - x)))
    rule251 = ReplacementRule(pattern251, lambda p, x, a, d, n, b, u, c : Subst(Int((a + b*sech(c + d*x**n))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule251)

    pattern252 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*csch(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: LinearQ(u, x)), CustomConstraint(lambda u, x: NonzeroQ(u - x)))
    rule252 = ReplacementRule(pattern252, lambda p, x, a, d, n, b, u, c : Subst(Int((a + b*csch(c + d*x**n))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule252)

    pattern253 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sech(u_))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule253 = ReplacementRule(pattern253, lambda p, x, a, u, b : Int((a + b*sech(ExpandToSum(u, x)))**p, x))
    rubi.add(rule253)

    pattern254 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*csch(u_))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule254 = ReplacementRule(pattern254, lambda p, x, a, u, b : Int((a + b*csch(ExpandToSum(u, x)))**p, x))
    rubi.add(rule254)

    pattern255 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sech(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: PositiveIntegerQ((m + S(1))/n)), CustomConstraint(lambda p: IntegerQ(p)))
    rule255 = ReplacementRule(pattern255, lambda p, m, x, a, d, n, b, c : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*sech(c + d*x))**p, x), x, x**n)/n)
    rubi.add(rule255)

    pattern256 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*csch(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: PositiveIntegerQ((m + S(1))/n)), CustomConstraint(lambda p: IntegerQ(p)))
    rule256 = ReplacementRule(pattern256, lambda p, m, x, a, d, n, b, c : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*csch(c + d*x))**p, x), x, x**n)/n)
    rubi.add(rule256)

    pattern257 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sech(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule257 = ReplacementRule(pattern257, lambda p, m, x, a, d, n, b, c : Int(x**m*(a + b*sech(c + d*x**n))**p, x))
    rubi.add(rule257)

    pattern258 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*csch(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule258 = ReplacementRule(pattern258, lambda p, m, x, a, d, n, b, c : Int(x**m*(a + b*csch(c + d*x**n))**p, x))
    rubi.add(rule258)

    pattern259 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sech(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule259 = ReplacementRule(pattern259, lambda p, m, e, x, a, d, n, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*sech(c + d*x**n))**p, x))
    rubi.add(rule259)

    pattern260 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*csch(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule260 = ReplacementRule(pattern260, lambda p, m, e, x, a, d, n, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*csch(c + d*x**n))**p, x))
    rubi.add(rule260)

    pattern261 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sech(u_))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule261 = ReplacementRule(pattern261, lambda p, m, e, x, a, u, b : Int((e*x)**m*(a + b*sech(ExpandToSum(u, x)))**p, x))
    rubi.add(rule261)

    pattern262 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*csch(u_))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda u, x: BinomialQ(u, x)), CustomConstraint(lambda u, x: Not(BinomialMatchQ(u, x))))
    rule262 = ReplacementRule(pattern262, lambda p, m, e, x, a, u, b : Int((e*x)**m*(a + b*csch(ExpandToSum(u, x)))**p, x))
    rubi.add(rule262)

    pattern263 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*sech(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n, m: GreaterEqual(m - n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(-1))))
    rule263 = ReplacementRule(pattern263, lambda p, m, x, a, n, b : -x**(m - n + S(1))*sech(a + b*x**n)**(p + S(-1))/(b*n*(p + S(-1))) + (m - n + S(1))*Int(x**(m - n)*sech(a + b*x**n)**(p + S(-1)), x)/(b*n*(p + S(-1))))
    rubi.add(rule263)

    pattern264 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*csch(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n, m: GreaterEqual(m - n, S(0))), CustomConstraint(lambda p: NonzeroQ(p + S(-1))))
    rule264 = ReplacementRule(pattern264, lambda p, m, x, a, n, b : -x**(m - n + S(1))*csch(a + b*x**n)**(p + S(-1))/(b*n*(p + S(-1))) + (m - n + S(1))*Int(x**(m - n)*csch(a + b*x**n)**(p + S(-1)), x)/(b*n*(p + S(-1))))
    rubi.add(rule264)

    pattern265 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule265 = ReplacementRule(pattern265, lambda m, x, a, d, n, b, c : -d*m*Int((c + d*x)**(m + S(-1))*sinh(a + b*x)**(n + S(1)), x)/(b*(n + S(1))) + (c + d*x)**m*sinh(a + b*x)**(n + S(1))/(b*(n + S(1))))
    rubi.add(rule265)

    pattern266 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule266 = ReplacementRule(pattern266, lambda m, x, a, d, n, b, c : -d*m*Int((c + d*x)**(m + S(-1))*Cosh(a + b*x)**(n + S(1)), x)/(b*(n + S(1))) + (c + d*x)**m*Cosh(a + b*x)**(n + S(1))/(b*(n + S(1))))
    rubi.add(rule266)

    pattern267 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(n, p)))
    rule267 = ReplacementRule(pattern267, lambda p, m, x, a, d, n, b, c : Int(ExpandTrigReduce((c + d*x)**m, Cosh(a + b*x)**p*sinh(a + b*x)**n, x), x))
    rubi.add(rule267)

    pattern268 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(n, p)))
    rule268 = ReplacementRule(pattern268, lambda p, m, c, x, a, d, n, b : Int((c + d*x)**m*sinh(a + b*x)**n*tanh(a + b*x)**(p + S(-2)), x) - Int((c + d*x)**m*sinh(a + b*x)**(n + S(-2))*tanh(a + b*x)**p, x))
    rubi.add(rule268)

    pattern269 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(n, p)))
    rule269 = ReplacementRule(pattern269, lambda p, m, x, a, d, n, b, c : Int((c + d*x)**m*Cosh(a + b*x)**n*coth(a + b*x)**(p + S(-2)), x) + Int((c + d*x)**m*Cosh(a + b*x)**(n + S(-2))*coth(a + b*x)**p, x))
    rubi.add(rule269)

    pattern270 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sech(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: SameQ(p, S(1))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule270 = ReplacementRule(pattern270, lambda p, m, c, x, a, d, n, b : d*m*Int((c + d*x)**(m + S(-1))*sech(a + b*x)**n, x)/(b*n) - (c + d*x)**m*sech(a + b*x)**n/(b*n))
    rubi.add(rule270)

    pattern271 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*csch(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: SameQ(p, S(1))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))))
    rule271 = ReplacementRule(pattern271, lambda p, m, x, a, d, n, b, c : d*m*Int((c + d*x)**(m + S(-1))*csch(a + b*x)**n, x)/(b*n) - (c + d*x)**m*csch(a + b*x)**n/(b*n))
    rubi.add(rule271)

    pattern272 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*sech(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule272 = ReplacementRule(pattern272, lambda m, c, x, a, d, n, b : -d*m*Int((c + d*x)**(m + S(-1))*tanh(a + b*x)**(n + S(1)), x)/(b*(n + S(1))) + (c + d*x)**m*tanh(a + b*x)**(n + S(1))/(b*(n + S(1))))
    rubi.add(rule272)

    pattern273 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*csch(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule273 = ReplacementRule(pattern273, lambda m, x, a, d, n, b, c : d*m*Int((c + d*x)**(m + S(-1))*coth(a + b*x)**(n + S(1)), x)/(b*(n + S(1))) - (c + d*x)**m*coth(a + b*x)**(n + S(1))/(b*(n + S(1))))
    rubi.add(rule273)

    pattern274 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**p_*sech(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: PositiveIntegerQ(p/S(2))))
    rule274 = ReplacementRule(pattern274, lambda p, m, c, x, a, d, b : Int((c + d*x)**m*tanh(a + b*x)**(p + S(-2))*sech(a + b*x), x) - Int((c + d*x)**m*tanh(a + b*x)**(p + S(-2))*sech(a + b*x)**S(3), x))
    rubi.add(rule274)

    pattern275 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**p_*sech(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p/S(2))))
    rule275 = ReplacementRule(pattern275, lambda p, m, c, x, a, d, n, b : Int((c + d*x)**m*tanh(a + b*x)**(p + S(-2))*sech(a + b*x)**n, x) - Int((c + d*x)**m*tanh(a + b*x)**(p + S(-2))*sech(a + b*x)**(n + S(2)), x))
    rubi.add(rule275)

    pattern276 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))**p_*csch(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: PositiveIntegerQ(p/S(2))))
    rule276 = ReplacementRule(pattern276, lambda p, m, x, a, d, b, c : Int((c + d*x)**m*coth(a + b*x)**(p + S(-2))*csch(a + b*x), x) + Int((c + d*x)**m*coth(a + b*x)**(p + S(-2))*csch(a + b*x)**S(3), x))
    rubi.add(rule276)

    pattern277 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))**p_*csch(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p/S(2))))
    rule277 = ReplacementRule(pattern277, lambda p, m, x, a, d, n, b, c : Int((c + d*x)**m*coth(a + b*x)**(p + S(-2))*csch(a + b*x)**n, x) + Int((c + d*x)**m*coth(a + b*x)**(p + S(-2))*csch(a + b*x)**(n + S(2)), x))
    rubi.add(rule277)

    pattern278 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sech(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n, p: EvenQ(n) | OddQ(p)), )
    def With278(p, m, c, x, a, d, n, b):
        u = IntHide(tanh(a + b*x)**p*sech(a + b*x)**n, x)
        return -d*m*Int(u*(c + d*x)**(m + S(-1)), x) + Dist((c + d*x)**m, u, x)
    rule278 = ReplacementRule(pattern278, lambda p, m, c, x, a, d, n, b : With278(p, m, c, x, a, d, n, b))
    rubi.add(rule278)

    pattern279 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*coth(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*csch(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n, p: EvenQ(n) | OddQ(p)), )
    def With279(p, m, x, a, d, n, b, c):
        u = IntHide(coth(a + b*x)**p*csch(a + b*x)**n, x)
        return -d*m*Int(u*(c + d*x)**(m + S(-1)), x) + Dist((c + d*x)**m, u, x)
    rule279 = ReplacementRule(pattern279, lambda p, m, x, a, d, n, b, c : With279(p, m, x, a, d, n, b, c))
    rubi.add(rule279)

    pattern280 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*csch(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*sech(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda n: IntegerQ(n)))
    rule280 = ReplacementRule(pattern280, lambda m, x, a, d, n, b, c : S(2)**n*Int((c + d*x)**m*csch(S(2)*a + S(2)*b*x)**n, x))
    rubi.add(rule280)

    pattern281 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*csch(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*sech(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, p: IntegersQ(n, p)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n, p: Unequal(n, p)), )
    def With281(p, m, x, a, d, n, b, c):
        u = IntHide(csch(a + b*x)**n*sech(a + b*x)**p, x)
        return -d*m*Int(u*(c + d*x)**(m + S(-1)), x) + Dist((c + d*x)**m, u, x)
    rule281 = ReplacementRule(pattern281, lambda p, m, x, a, d, n, b, c : With281(p, m, x, a, d, n, b, c))
    rubi.add(rule281)

    pattern282 = Pattern(Integral(F_**v_*G_**w_*u_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda F: HyperbolicQ(F)), CustomConstraint(lambda G: HyperbolicQ(G)), CustomConstraint(lambda v, w: ZeroQ(v - w)), CustomConstraint(lambda v, u, x, w: LinearQ(List(u, v, w), x)), CustomConstraint(lambda v, u, x, w: Not(LinearMatchQ(List(u, v, w), x))))
    rule282 = ReplacementRule(pattern282, lambda p, m, x, w, n, v, u, F, G : Int(ExpandToSum(u, x)**m*F(ExpandToSum(v, x))**n*G(ExpandToSum(v, x))**p, x))
    rubi.add(rule282)

    pattern283 = Pattern(Integral((a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule283 = ReplacementRule(pattern283, lambda m, e, x, f, d, a, n, b, c : -f*m*Int((a + b*sinh(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) + (a + b*sinh(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule283)

    pattern284 = Pattern(Integral((a_ + Cosh(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule284 = ReplacementRule(pattern284, lambda m, c, e, x, f, d, a, n, b : -f*m*Int((a + b*Cosh(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) + (a + b*Cosh(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule284)

    pattern285 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sech(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule285 = ReplacementRule(pattern285, lambda m, e, x, f, d, a, n, b, c : -f*m*Int((a + b*tanh(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) + (a + b*tanh(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule285)

    pattern286 = Pattern(Integral((a_ + WC('b', S(1))*coth(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*csch(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule286 = ReplacementRule(pattern286, lambda m, c, e, x, f, d, a, n, b : f*m*Int((a + b*coth(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) - (a + b*coth(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule286)

    pattern287 = Pattern(Integral((a_ + WC('b', S(1))*sech(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*tanh(x_*WC('d', S(1)) + WC('c', S(0)))*sech(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule287 = ReplacementRule(pattern287, lambda m, e, x, f, d, a, n, b, c : f*m*Int((a + b*sech(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) - (a + b*sech(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule287)

    pattern288 = Pattern(Integral((a_ + WC('b', S(1))*csch(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*coth(x_*WC('d', S(1)) + WC('c', S(0)))*csch(x_*WC('d', S(1)) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: NonzeroQ(n + S(1))))
    rule288 = ReplacementRule(pattern288, lambda m, c, e, x, f, d, a, n, b : f*m*Int((a + b*csch(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x)/(b*d*(n + S(1))) - (a + b*csch(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))))
    rubi.add(rule288)

    pattern289 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)), CustomConstraint(lambda m: IntegerQ(m)))
    rule289 = ReplacementRule(pattern289, lambda p, m, c, e, x, a, f, d, q, b : Int(ExpandTrigReduce((e + f*x)**m, sinh(a + b*x)**p*sinh(c + d*x)**q, x), x))
    rubi.add(rule289)

    pattern290 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*Cosh(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)), CustomConstraint(lambda m: IntegerQ(m)))
    rule290 = ReplacementRule(pattern290, lambda p, m, e, x, a, d, f, b, q, c : Int(ExpandTrigReduce((e + f*x)**m, Cosh(a + b*x)**p*Cosh(c + d*x)**q, x), x))
    rubi.add(rule290)

    pattern291 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cosh(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)))
    rule291 = ReplacementRule(pattern291, lambda p, m, e, x, f, d, a, b, q, c : Int(ExpandTrigReduce((e + f*x)**m, Cosh(c + d*x)**q*sinh(a + b*x)**p, x), x))
    rubi.add(rule291)

    pattern292 = Pattern(Integral(F_**(x_*WC('b', S(1)) + WC('a', S(0)))*G_**(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda F: MemberQ(List(Sinh, Cosh), F)), CustomConstraint(lambda G: MemberQ(List(Sech, Csch), G)), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)), CustomConstraint(lambda d, c, a, b: ZeroQ(-a*d + b*c)), CustomConstraint(lambda d, b: PositiveIntegerQ(b/d + S(-1))))
    rule292 = ReplacementRule(pattern292, lambda p, m, c, e, x, a, f, d, F, q, G, b : Int(ExpandTrigExpand((e + f*x)**m*G(c + d*x)**q, F, c + d*x, p, b/d, x), x))
    rubi.add(rule292)

    pattern293 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sinh(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda F, c, e, b: NonzeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2))))
    rule293 = ReplacementRule(pattern293, lambda c, e, x, a, d, F, b : -F**(c*(a + b*x))*b*c*log(F)*sinh(d + e*x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)) + F**(c*(a + b*x))*e*Cosh(d + e*x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)))
    rubi.add(rule293)

    pattern294 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cosh(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda F, c, e, b: NonzeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2))))
    rule294 = ReplacementRule(pattern294, lambda c, e, x, a, d, F, b : -F**(c*(a + b*x))*b*c*Cosh(d + e*x)*log(F)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)) + F**(c*(a + b*x))*e*sinh(d + e*x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)))
    rubi.add(rule294)

    pattern295 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, n, b, F, c: NonzeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule295 = ReplacementRule(pattern295, lambda c, e, x, a, d, n, F, b : -F**(c*(a + b*x))*b*c*log(F)*sinh(d + e*x)**n/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)) + F**(c*(a + b*x))*e*n*Cosh(d + e*x)*sinh(d + e*x)**(n + S(-1))/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)) - e**S(2)*n*(n + S(-1))*Int(F**(c*(a + b*x))*sinh(d + e*x)**(n + S(-2)), x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)))
    rubi.add(rule295)

    pattern296 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cosh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, n, b, F, c: NonzeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule296 = ReplacementRule(pattern296, lambda c, e, x, a, d, n, F, b : -F**(c*(a + b*x))*b*c*Cosh(d + e*x)**n*log(F)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)) + F**(c*(a + b*x))*e*n*Cosh(d + e*x)**(n + S(-1))*sinh(d + e*x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)) + e**S(2)*n*(n + S(-1))*Int(F**(c*(a + b*x))*Cosh(d + e*x)**(n + S(-2)), x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)))
    rubi.add(rule296)

    pattern297 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, n, b, F, c: ZeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))), CustomConstraint(lambda n: NonzeroQ(n + S(1))), CustomConstraint(lambda n: NonzeroQ(n + S(2))))
    rule297 = ReplacementRule(pattern297, lambda c, e, x, a, d, n, F, b : -F**(c*(a + b*x))*b*c*log(F)*sinh(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))) + F**(c*(a + b*x))*Cosh(d + e*x)*sinh(d + e*x)**(n + S(1))/(e*(n + S(1))))
    rubi.add(rule297)

    pattern298 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cosh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, n, b, F, c: ZeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))), CustomConstraint(lambda n: NonzeroQ(n + S(1))), CustomConstraint(lambda n: NonzeroQ(n + S(2))))
    rule298 = ReplacementRule(pattern298, lambda c, e, x, a, d, n, F, b : F**(c*(a + b*x))*b*c*Cosh(d + e*x)**(n + S(2))*log(F)/(e**S(2)*(n + S(1))*(n + S(2))) - F**(c*(a + b*x))*Cosh(d + e*x)**(n + S(1))*sinh(d + e*x)/(e*(n + S(1))))
    rubi.add(rule298)

    pattern299 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, n, b, F, c: NonzeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule299 = ReplacementRule(pattern299, lambda c, e, x, a, d, n, F, b : -F**(c*(a + b*x))*b*c*log(F)*sinh(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))) + F**(c*(a + b*x))*Cosh(d + e*x)*sinh(d + e*x)**(n + S(1))/(e*(n + S(1))) - (-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))*Int(F**(c*(a + b*x))*sinh(d + e*x)**(n + S(2)), x)/(e**S(2)*(n + S(1))*(n + S(2))))
    rubi.add(rule299)

    pattern300 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cosh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, n, b, F, c: NonzeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda n: Unequal(n, S(-2))))
    rule300 = ReplacementRule(pattern300, lambda c, e, x, a, d, n, F, b : F**(c*(a + b*x))*b*c*Cosh(d + e*x)**(n + S(2))*log(F)/(e**S(2)*(n + S(1))*(n + S(2))) - F**(c*(a + b*x))*Cosh(d + e*x)**(n + S(1))*sinh(d + e*x)/(e*(n + S(1))) + (-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))*Int(F**(c*(a + b*x))*Cosh(d + e*x)**(n + S(2)), x)/(e**S(2)*(n + S(1))*(n + S(2))))
    rubi.add(rule300)

    pattern301 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule301 = ReplacementRule(pattern301, lambda c, e, x, a, d, n, F, b : (exp(S(2)*d + S(2)*e*x) + S(-1))**(-n)*Int(F**(c*(a + b*x))*(exp(S(2)*d + S(2)*e*x) + S(-1))**n*exp(-n*(d + e*x)), x)*exp(n*(d + e*x))*sinh(d + e*x)**n)
    rubi.add(rule301)

    pattern302 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cosh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule302 = ReplacementRule(pattern302, lambda c, e, x, a, d, n, F, b : (exp(S(2)*d + S(2)*e*x) + S(1))**(-n)*Cosh(d + e*x)**n*Int(F**(c*(a + b*x))*(exp(S(2)*d + S(2)*e*x) + S(1))**n*exp(-n*(d + e*x)), x)*exp(n*(d + e*x)))
    rubi.add(rule302)

    pattern303 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*tanh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule303 = ReplacementRule(pattern303, lambda c, e, x, a, d, n, F, b : Int(ExpandIntegrand(F**(c*(a + b*x))*(exp(S(2)*d + S(2)*e*x) + S(-1))**n*(exp(S(2)*d + S(2)*e*x) + S(1))**(-n), x), x))
    rubi.add(rule303)

    pattern304 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*coth(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule304 = ReplacementRule(pattern304, lambda e, x, a, d, n, b, F, c : Int(ExpandIntegrand(F**(c*(a + b*x))*(exp(S(2)*d + S(2)*e*x) + S(-1))**(-n)*(exp(S(2)*d + S(2)*e*x) + S(1))**n, x), x))
    rubi.add(rule304)

    pattern305 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sech(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, n, b, F, c: NonzeroQ(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule305 = ReplacementRule(pattern305, lambda c, e, x, a, d, n, F, b : -F**(c*(a + b*x))*b*c*log(F)*sech(d + e*x)**n/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)) - F**(c*(a + b*x))*e*n*sinh(d + e*x)*sech(d + e*x)**(n + S(1))/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)) + e**S(2)*n*(n + S(1))*Int(F**(c*(a + b*x))*sech(d + e*x)**(n + S(2)), x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)))
    rubi.add(rule305)

    pattern306 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*csch(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, n, b, F, c: NonzeroQ(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))))
    rule306 = ReplacementRule(pattern306, lambda c, e, x, a, d, n, F, b : -F**(c*(a + b*x))*b*c*log(F)*csch(d + e*x)**n/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)) - F**(c*(a + b*x))*e*n*Cosh(d + e*x)*csch(d + e*x)**(n + S(1))/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)) - e**S(2)*n*(n + S(1))*Int(F**(c*(a + b*x))*csch(d + e*x)**(n + S(2)), x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)))
    rubi.add(rule306)

    pattern307 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sech(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, n, b, F, c: ZeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))), CustomConstraint(lambda n: NonzeroQ(n + S(-1))), CustomConstraint(lambda n: NonzeroQ(n + S(-2))))
    rule307 = ReplacementRule(pattern307, lambda c, e, x, a, d, n, F, b : F**(c*(a + b*x))*b*c*log(F)*sech(d + e*x)**(n + S(-2))/(e**S(2)*(n + S(-2))*(n + S(-1))) + F**(c*(a + b*x))*sinh(d + e*x)*sech(d + e*x)**(n + S(-1))/(e*(n + S(-1))))
    rubi.add(rule307)

    pattern308 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*csch(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda e, n, b, F, c: ZeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))), CustomConstraint(lambda n: NonzeroQ(n + S(-1))), CustomConstraint(lambda n: NonzeroQ(n + S(-2))))
    rule308 = ReplacementRule(pattern308, lambda c, e, x, a, d, n, F, b : -F**(c*(a + b*x))*b*c*log(F)*csch(d + e*x)**(n + S(-2))/(e**S(2)*(n + S(-2))*(n + S(-1))) - F**(c*(a + b*x))*Cosh(d + e*x)*csch(d + e*x)**(n + S(-1))/(e*(n + S(-1))))
    rubi.add(rule308)

    pattern309 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sech(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, n, b, F, c: NonzeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda n: Unequal(n, S(2))))
    rule309 = ReplacementRule(pattern309, lambda c, e, x, a, d, n, F, b : F**(c*(a + b*x))*b*c*log(F)*sech(d + e*x)**(n + S(-2))/(e**S(2)*(n + S(-2))*(n + S(-1))) + F**(c*(a + b*x))*sinh(d + e*x)*sech(d + e*x)**(n + S(-1))/(e*(n + S(-1))) + (-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))*Int(F**(c*(a + b*x))*sech(d + e*x)**(n + S(-2)), x)/(e**S(2)*(n + S(-2))*(n + S(-1))))
    rubi.add(rule309)

    pattern310 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*csch(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, n, b, F, c: NonzeroQ(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda n: Unequal(n, S(2))))
    rule310 = ReplacementRule(pattern310, lambda c, e, x, a, d, n, F, b : -F**(c*(a + b*x))*b*c*log(F)*csch(d + e*x)**(n + S(-2))/(e**S(2)*(n + S(-2))*(n + S(-1))) - F**(c*(a + b*x))*Cosh(d + e*x)*csch(d + e*x)**(n + S(-1))/(e*(n + S(-1))) - (-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))*Int(F**(c*(a + b*x))*csch(d + e*x)**(n + S(-2)), x)/(e**S(2)*(n + S(-2))*(n + S(-1))))
    rubi.add(rule310)

    pattern311 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sech(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule311 = ReplacementRule(pattern311, lambda c, e, x, a, d, n, F, b : S(2)**n*F**(c*(a + b*x))*Hypergeometric2F1(n, b*c*log(F)/(S(2)*e) + n/S(2), b*c*log(F)/(S(2)*e) + n/S(2) + S(1), -exp(S(2)*d + S(2)*e*x))*exp(n*(d + e*x))/(b*c*log(F) + e*n))
    rubi.add(rule311)

    pattern312 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*csch(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule312 = ReplacementRule(pattern312, lambda e, x, a, d, n, b, F, c : (S(-2))**n*F**(c*(a + b*x))*Hypergeometric2F1(n, b*c*log(F)/(S(2)*e) + n/S(2), b*c*log(F)/(S(2)*e) + n/S(2) + S(1), exp(S(2)*d + S(2)*e*x))*exp(n*(d + e*x))/(b*c*log(F) + e*n))
    rubi.add(rule312)

    pattern313 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sech(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule313 = ReplacementRule(pattern313, lambda c, e, x, a, d, n, F, b : (exp(S(2)*d + S(2)*e*x) + S(1))**n*Int(SimplifyIntegrand(F**(c*(a + b*x))*(exp(S(2)*d + S(2)*e*x) + S(1))**(-n)*exp(n*(d + e*x)), x), x)*exp(-n*(d + e*x))*sech(d + e*x)**n)
    rubi.add(rule313)

    pattern314 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*csch(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule314 = ReplacementRule(pattern314, lambda e, x, a, d, n, b, F, c : (-exp(-S(2)*d - S(2)*e*x) + S(1))**n*Int(SimplifyIntegrand(F**(c*(a + b*x))*(-exp(-S(2)*d - S(2)*e*x) + S(1))**(-n)*exp(-n*(d + e*x)), x), x)*exp(n*(d + e*x))*csch(d + e*x)**n)
    rubi.add(rule314)

    pattern315 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda g, f: ZeroQ(f**S(2) + g**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule315 = ReplacementRule(pattern315, lambda c, e, x, a, d, f, g, n, F, b : S(2)**n*f**n*Int(F**(c*(a + b*x))*Cosh(-Pi*f/(S(4)*g) + d/S(2) + e*x/S(2))**(S(2)*n), x))
    rubi.add(rule315)

    pattern316 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + Cosh(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda g, f: ZeroQ(f - g)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule316 = ReplacementRule(pattern316, lambda e, x, a, d, f, g, n, b, F, c : S(2)**n*g**n*Int(F**(c*(a + b*x))*Cosh(d/S(2) + e*x/S(2))**(S(2)*n), x))
    rubi.add(rule316)

    pattern317 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + Cosh(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda g, f: ZeroQ(f + g)), CustomConstraint(lambda n: NegativeIntegerQ(n)))
    rule317 = ReplacementRule(pattern317, lambda e, x, a, d, f, g, n, b, F, c : S(2)**n*g**n*Int(F**(c*(a + b*x))*sinh(d/S(2) + e*x/S(2))**(S(2)*n), x))
    rubi.add(rule317)

    pattern318 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*Cosh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda g, f: ZeroQ(f**S(2) + g**S(2))), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Equal(m + n, S(0))))
    rule318 = ReplacementRule(pattern318, lambda m, e, x, a, d, f, g, n, b, F, c : g**n*Int(F**(c*(a + b*x))*tanh(-Pi*f/(S(4)*g) + d/S(2) + e*x/S(2))**m, x))
    rubi.add(rule318)

    pattern319 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + Cosh(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1)))**WC('n', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda g, f: ZeroQ(f - g)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Equal(m + n, S(0))))
    rule319 = ReplacementRule(pattern319, lambda m, e, x, a, d, f, g, n, b, F, c : g**n*Int(F**(c*(a + b*x))*tanh(d/S(2) + e*x/S(2))**m, x))
    rubi.add(rule319)

    pattern320 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + Cosh(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1)))**WC('n', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda g, f: ZeroQ(f + g)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Equal(m + n, S(0))))
    rule320 = ReplacementRule(pattern320, lambda m, e, x, a, d, f, g, n, b, F, c : g**n*Int(F**(c*(a + b*x))*coth(d/S(2) + e*x/S(2))**m, x))
    rubi.add(rule320)

    pattern321 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(h_ + Cosh(x_*WC('e', S(1)) + WC('d', S(0)))*WC('i', S(1)))/(f_ + WC('g', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda g, f: ZeroQ(f**S(2) + g**S(2))), CustomConstraint(lambda h, i: ZeroQ(h**S(2) - i**S(2))), CustomConstraint(lambda h, g, i, f: ZeroQ(-f*i + g*h)))
    rule321 = ReplacementRule(pattern321, lambda h, c, e, x, a, d, f, g, F, i, b : S(2)*i*Int(F**(c*(a + b*x))*Cosh(d + e*x)/(f + g*sinh(d + e*x)), x) + Int(F**(c*(a + b*x))*(h - i*Cosh(d + e*x))/(f + g*sinh(d + e*x)), x))
    rubi.add(rule321)

    pattern322 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(h_ + WC('i', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0))))/(f_ + Cosh(x_*WC('e', S(1)) + WC('d', S(0)))*WC('g', S(1))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda g, f: ZeroQ(f**S(2) - g**S(2))), CustomConstraint(lambda h, i: ZeroQ(h**S(2) + i**S(2))), CustomConstraint(lambda h, g, i, f: ZeroQ(f*i + g*h)))
    rule322 = ReplacementRule(pattern322, lambda h, e, x, a, d, f, g, F, b, i, c : S(2)*i*Int(F**(c*(a + b*x))*sinh(d + e*x)/(f + g*Cosh(d + e*x)), x) + Int(F**(c*(a + b*x))*(h - i*sinh(d + e*x))/(f + g*Cosh(d + e*x)), x))
    rubi.add(rule322)

    pattern323 = Pattern(Integral(F_**(u_*WC('c', S(1)))*G_**v_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda G: HyperbolicQ(G)), CustomConstraint(lambda v, u, x: LinearQ(List(u, v), x)), CustomConstraint(lambda v, u, x: Not(LinearMatchQ(List(u, v), x))))
    rule323 = ReplacementRule(pattern323, lambda x, n, v, u, F, G, c : Int(F**(c*ExpandToSum(u, x))*G(ExpandToSum(v, x))**n, x))
    rubi.add(rule323)

    pattern324 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: PositiveIntegerQ(n)), )
    def With324(m, c, e, x, a, d, n, F, b):
        u = IntHide(F**(c*(a + b*x))*sinh(d + e*x)**n, x)
        return u*x**m - Dist(m, Int(u*x**(m + S(-1)), x))
    rule324 = ReplacementRule(pattern324, lambda m, c, e, x, a, d, n, F, b : With324(m, c, e, x, a, d, n, F, b))
    rubi.add(rule324)

    pattern325 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*Cosh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: PositiveIntegerQ(n)), )
    def With325(m, e, x, a, d, n, b, F, c):
        u = IntHide(F**(c*(a + b*x))*Cosh(d + e*x)**n, x)
        return u*x**m - Dist(m, Int(u*x**(m + S(-1)), x))
    rule325 = ReplacementRule(pattern325, lambda m, e, x, a, d, n, b, F, c : With325(m, e, x, a, d, n, b, F, c))
    rubi.add(rule325)

    pattern326 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*Cosh(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule326 = ReplacementRule(pattern326, lambda m, e, x, a, f, d, g, n, b, F, c : Int(ExpandTrigReduce(F**(c*(a + b*x)), Cosh(f + g*x)**n*sinh(d + e*x)**m, x), x))
    rubi.add(rule326)

    pattern327 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('p', S(1))*Cosh(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, m, p: PositiveIntegerQ(m, n, p)))
    rule327 = ReplacementRule(pattern327, lambda p, m, e, x, a, f, d, g, n, b, F, c : Int(ExpandTrigReduce(F**(c*(a + b*x))*x**p, Cosh(f + g*x)**n*sinh(d + e*x)**m, x), x))
    rubi.add(rule327)

    pattern328 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*G_**(x_*WC('e', S(1)) + WC('d', S(0)))*H_**(x_*WC('e', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)), CustomConstraint(lambda G: HyperbolicQ(G)), CustomConstraint(lambda H: HyperbolicQ(H)))
    rule328 = ReplacementRule(pattern328, lambda m, H, c, e, x, a, d, n, F, G, b : Int(ExpandTrigToExp(F**(c*(a + b*x)), G(d + e*x)**m*H(d + e*x)**n, x), x))
    rubi.add(rule328)

    pattern329 = Pattern(Integral(F_**u_*sinh(v_)**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda u, x: LinearQ(u, x) | PolyQ(u, x, S(2))), CustomConstraint(lambda v, x: LinearQ(v, x) | PolyQ(v, x, S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule329 = ReplacementRule(pattern329, lambda x, n, v, u, F : Int(ExpandTrigToExp(F**u, sinh(v)**n, x), x))
    rubi.add(rule329)

    pattern330 = Pattern(Integral(F_**u_*Cosh(v_)**WC('n', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda u, x: LinearQ(u, x) | PolyQ(u, x, S(2))), CustomConstraint(lambda v, x: LinearQ(v, x) | PolyQ(v, x, S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule330 = ReplacementRule(pattern330, lambda x, n, v, u, F : Int(ExpandTrigToExp(F**u, Cosh(v)**n, x), x))
    rubi.add(rule330)

    pattern331 = Pattern(Integral(F_**u_*Cosh(v_)**WC('n', S(1))*sinh(v_)**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda u, x: LinearQ(u, x) | PolyQ(u, x, S(2))), CustomConstraint(lambda v, x: LinearQ(v, x) | PolyQ(v, x, S(2))), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule331 = ReplacementRule(pattern331, lambda m, x, n, v, u, F : Int(ExpandTrigToExp(F**u, Cosh(v)**n*sinh(v)**m, x), x))
    rubi.add(rule331)

    pattern332 = Pattern(Integral(sinh(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, p, b: RationalQ(b, n, p)))
    rule332 = ReplacementRule(pattern332, lambda p, c, x, n, b : Int(((c*x**n)**b/S(2) - (c*x**n)**(-b)/S(2))**p, x))
    rubi.add(rule332)

    pattern333 = Pattern(Integral(Cosh(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, p, b: RationalQ(b, n, p)))
    rule333 = ReplacementRule(pattern333, lambda p, c, x, n, b : Int(((c*x**n)**b/S(2) + (c*x**n)**(-b)/S(2))**p, x))
    rubi.add(rule333)

    pattern334 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, p, b: ZeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + S(-1))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule334 = ReplacementRule(pattern334, lambda p, x, a, n, b, c : -x*(p + S(2))*sinh(a + b*log(c*x**n))**(p + S(2))/(p + S(1)) + x*sinh(a + b*log(c*x**n))**(p + S(2))*coth(a + b*log(c*x**n))/(b*n*(p + S(1))))
    rubi.add(rule334)

    pattern335 = Pattern(Integral(Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, p, b: ZeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + S(-1))), CustomConstraint(lambda p: NonzeroQ(p + S(1))))
    rule335 = ReplacementRule(pattern335, lambda p, x, a, n, b, c : x*(p + S(2))*Cosh(a + b*log(c*x**n))**(p + S(2))/(p + S(1)) - x*Cosh(a + b*log(c*x**n))**(p + S(2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(1))))
    rubi.add(rule335)

    pattern336 = Pattern(Integral(sqrt(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, b: ZeroQ(b*n + S(-2))))
    rule336 = ReplacementRule(pattern336, lambda x, a, n, b, c : x*Int(sqrt((c*x**n)**(S(4)/n)*exp(S(2)*a) + S(-1))/x, x)*sqrt(sinh(a + b*log(c*x**n)))/sqrt((c*x**n)**(S(4)/n)*exp(S(2)*a) + S(-1)))
    rubi.add(rule336)

    pattern337 = Pattern(Integral(sqrt(Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, b: ZeroQ(b*n + S(-2))))
    rule337 = ReplacementRule(pattern337, lambda x, a, n, b, c : x*sqrt(Cosh(a + b*log(c*x**n)))*Int(sqrt((c*x**n)**(S(4)/n)*exp(S(2)*a) + S(1))/x, x)/sqrt((c*x**n)**(S(4)/n)*exp(S(2)*a) + S(1)))
    rubi.add(rule337)

    pattern338 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda n, p, b: ZeroQ(b**S(2)*n**S(2)*p**S(2) + S(-1))))
    rule338 = ReplacementRule(pattern338, lambda p, x, a, n, b, c : Int(ExpandIntegrand(((c*x**n)**(S(1)/(n*p))*exp(a*b*n*p)/(S(2)*b*n*p) - (c*x**n)**(-S(1)/(n*p))*exp(-a*b*n*p)/(S(2)*b*n*p))**p, x), x))
    rubi.add(rule338)

    pattern339 = Pattern(Integral(Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda n, p, b: ZeroQ(b**S(2)*n**S(2)*p**S(2) + S(-1))))
    rule339 = ReplacementRule(pattern339, lambda p, x, a, n, b, c : Int(ExpandIntegrand(((c*x**n)**(S(1)/(n*p))*exp(a*b*n*p)/S(2) + (c*x**n)**(-S(1)/(n*p))*exp(-a*b*n*p)/S(2))**p, x), x))
    rubi.add(rule339)

    pattern340 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, b: NonzeroQ(b**S(2)*n**S(2) + S(-1))))
    rule340 = ReplacementRule(pattern340, lambda x, a, n, b, c : b*n*x*Cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(-1)) - x*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(-1)))
    rubi.add(rule340)

    pattern341 = Pattern(Integral(Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, b: NonzeroQ(b**S(2)*n**S(2) + S(-1))))
    rule341 = ReplacementRule(pattern341, lambda x, a, n, b, c : b*n*x*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(-1)) - x*Cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(-1)))
    rubi.add(rule341)

    pattern342 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(-1))))
    rule342 = ReplacementRule(pattern342, lambda p, x, a, n, b, c : -b**S(2)*n**S(2)*p*(p + S(-1))*Int(sinh(a + b*log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*p**S(2) + S(-1)) + b*n*p*x*Cosh(a + b*log(c*x**n))*sinh(a + b*log(c*x**n))**(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + S(-1)) - x*sinh(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(-1)))
    rubi.add(rule342)

    pattern343 = Pattern(Integral(Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(-1))))
    rule343 = ReplacementRule(pattern343, lambda p, x, a, n, b, c : b**S(2)*n**S(2)*p*(p + S(-1))*Int(Cosh(a + b*log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*p**S(2) + S(-1)) + b*n*p*x*Cosh(a + b*log(c*x**n))**(p + S(-1))*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + S(-1)) - x*Cosh(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(-1)))
    rubi.add(rule343)

    pattern344 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + S(-1))))
    rule344 = ReplacementRule(pattern344, lambda p, x, a, n, b, c : x*sinh(a + b*log(c*x**n))**(p + S(2))*coth(a + b*log(c*x**n))/(b*n*(p + S(1))) - x*sinh(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))) - (b**S(2)*n**S(2)*(p + S(2))**S(2) + S(-1))*Int(sinh(a + b*log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))))
    rubi.add(rule344)

    pattern345 = Pattern(Integral(Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) + S(-1))))
    rule345 = ReplacementRule(pattern345, lambda p, x, a, n, b, c : -x*Cosh(a + b*log(c*x**n))**(p + S(2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(1))) + x*Cosh(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))) + (b**S(2)*n**S(2)*(p + S(2))**S(2) + S(-1))*Int(Cosh(a + b*log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))))
    rubi.add(rule345)

    pattern346 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(-1))))
    rule346 = ReplacementRule(pattern346, lambda p, x, a, n, b, c : x*(S(2) - S(2)*(c*x**n)**(-S(2)*b)*exp(-S(2)*a))**(-p)*((c*x**n)**b*exp(a) - (c*x**n)**(-b)*exp(-a))**p*Hypergeometric2F1(-p, (-b*n*p + S(-1))/(S(2)*b*n), S(1) - (b*n*p + S(1))/(S(2)*b*n), (c*x**n)**(-S(2)*b)*exp(-S(2)*a))/(b*n*p + S(1)))
    rubi.add(rule346)

    pattern347 = Pattern(Integral(Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(-1))))
    rule347 = ReplacementRule(pattern347, lambda p, x, a, n, b, c : x*(S(2) + S(2)*(c*x**n)**(-S(2)*b)*exp(-S(2)*a))**(-p)*((c*x**n)**b*exp(a) + (c*x**n)**(-b)*exp(-a))**p*Hypergeometric2F1(-p, (-b*n*p + S(-1))/(S(2)*b*n), S(1) - (b*n*p + S(1))/(S(2)*b*n), -(c*x**n)**(-S(2)*b)*exp(-S(2)*a))/(b*n*p + S(1)))
    rubi.add(rule347)

    pattern348 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m, p, b: ZeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) - (m + S(1))**S(2))), CustomConstraint(lambda p: NonzeroQ(p + S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule348 = ReplacementRule(pattern348, lambda p, m, x, a, n, b, c : x**(m + S(1))*(-p + S(-2))*sinh(a + b*log(c*x**n))**(p + S(2))/((m + S(1))*(p + S(1))) + x**(m + S(1))*sinh(a + b*log(c*x**n))**(p + S(2))*coth(a + b*log(c*x**n))/(b*n*(p + S(1))))
    rubi.add(rule348)

    pattern349 = Pattern(Integral(x_**WC('m', S(1))*Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m, p, b: ZeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) - (m + S(1))**S(2))), CustomConstraint(lambda p: NonzeroQ(p + S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule349 = ReplacementRule(pattern349, lambda p, m, x, a, n, b, c : x**(m + S(1))*(p + S(2))*Cosh(a + b*log(c*x**n))**(p + S(2))/((m + S(1))*(p + S(1))) - x**(m + S(1))*Cosh(a + b*log(c*x**n))**(p + S(2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(1))))
    rubi.add(rule349)

    pattern350 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda n, m, p, b: ZeroQ(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2))))
    rule350 = ReplacementRule(pattern350, lambda p, m, x, a, n, b, c : S(2)**(-p)*Int(ExpandIntegrand(x**m*((c*x**n)**((-m + S(-1))/(n*p))*(-m + S(-1))*exp(-a*b*n*p/(m + S(1)))/(b*n*p) + (c*x**n)**((m + S(1))/(n*p))*(m + S(1))*exp(a*b*n*p/(m + S(1)))/(b*n*p))**p, x), x))
    rubi.add(rule350)

    pattern351 = Pattern(Integral(x_**WC('m', S(1))*Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda n, m, p, b: ZeroQ(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2))))
    rule351 = ReplacementRule(pattern351, lambda p, m, x, a, n, b, c : S(2)**(-p)*Int(ExpandIntegrand(x**m*((c*x**n)**((-m + S(-1))/(n*p))*exp(-a*b*n*p/(m + S(1))) + (c*x**n)**((m + S(1))/(n*p))*exp(a*b*n*p/(m + S(1))))**p, x), x))
    rubi.add(rule351)

    pattern352 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m, b: NonzeroQ(b**S(2)*n**S(2) - (m + S(1))**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule352 = ReplacementRule(pattern352, lambda m, x, a, n, b, c : b*n*x**(m + S(1))*Cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2) - (m + S(1))**S(2)) + x**(m + S(1))*(-m + S(-1))*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2) - (m + S(1))**S(2)))
    rubi.add(rule352)

    pattern353 = Pattern(Integral(x_**WC('m', S(1))*Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m, b: NonzeroQ(b**S(2)*n**S(2) - (m + S(1))**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule353 = ReplacementRule(pattern353, lambda m, x, a, n, b, c : b*n*x**(m + S(1))*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2) - (m + S(1))**S(2)) + x**(m + S(1))*(-m + S(-1))*Cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2) - (m + S(1))**S(2)))
    rubi.add(rule353)

    pattern354 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule354 = ReplacementRule(pattern354, lambda p, m, x, a, n, b, c : -b**S(2)*n**S(2)*p*(p + S(-1))*Int(x**m*sinh(a + b*log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)) + b*n*p*x**(m + S(1))*Cosh(a + b*log(c*x**n))*sinh(a + b*log(c*x**n))**(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)) + x**(m + S(1))*(-m + S(-1))*sinh(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)))
    rubi.add(rule354)

    pattern355 = Pattern(Integral(x_**WC('m', S(1))*Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule355 = ReplacementRule(pattern355, lambda p, m, x, a, n, b, c : b**S(2)*n**S(2)*p*(p + S(-1))*Int(x**m*Cosh(a + b*log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)) + b*n*p*x**(m + S(1))*Cosh(a + b*log(c*x**n))**(p + S(-1))*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)) + x**(m + S(1))*(-m + S(-1))*Cosh(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)))
    rubi.add(rule355)

    pattern356 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) - (m + S(1))**S(2))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule356 = ReplacementRule(pattern356, lambda p, m, x, a, n, b, c : x**(m + S(1))*sinh(a + b*log(c*x**n))**(p + S(2))*coth(a + b*log(c*x**n))/(b*n*(p + S(1))) - x**(m + S(1))*(m + S(1))*sinh(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))) - (b**S(2)*n**S(2)*(p + S(2))**S(2) - (m + S(1))**S(2))*Int(x**m*sinh(a + b*log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))))
    rubi.add(rule356)

    pattern357 = Pattern(Integral(x_**WC('m', S(1))*Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(2))**S(2) - (m + S(1))**S(2))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda p: Unequal(p, S(-2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule357 = ReplacementRule(pattern357, lambda p, m, x, a, n, b, c : -x**(m + S(1))*Cosh(a + b*log(c*x**n))**(p + S(2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(1))) + x**(m + S(1))*(m + S(1))*Cosh(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))) + (b**S(2)*n**S(2)*(p + S(2))**S(2) - (m + S(1))**S(2))*Int(x**m*Cosh(a + b*log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))))
    rubi.add(rule357)

    pattern358 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2))))
    rule358 = ReplacementRule(pattern358, lambda p, m, x, a, n, b, c : x**(m + S(1))*(S(2) - S(2)*(c*x**n)**(-S(2)*b)*exp(-S(2)*a))**(-p)*((c*x**n)**b*exp(a) - (c*x**n)**(-b)*exp(-a))**p*Hypergeometric2F1(-p, (-b*n*p - m + S(-1))/(S(2)*b*n), S(1) - (b*n*p + m + S(1))/(S(2)*b*n), (c*x**n)**(-S(2)*b)*exp(-S(2)*a))/(b*n*p + m + S(1)))
    rubi.add(rule358)

    pattern359 = Pattern(Integral(x_**WC('m', S(1))*Cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2))))
    rule359 = ReplacementRule(pattern359, lambda p, m, x, a, n, b, c : x**(m + S(1))*(S(2) + S(2)*(c*x**n)**(-S(2)*b)*exp(-S(2)*a))**(-p)*((c*x**n)**b*exp(a) + (c*x**n)**(-b)*exp(-a))**p*Hypergeometric2F1(-p, (-b*n*p - m + S(-1))/(S(2)*b*n), S(1) - (b*n*p + m + S(1))/(S(2)*b*n), -(c*x**n)**(-S(2)*b)*exp(-S(2)*a))/(b*n*p + m + S(1)))
    rubi.add(rule359)

    pattern360 = Pattern(Integral(sech(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, p, b: RationalQ(b, n, p)))
    rule360 = ReplacementRule(pattern360, lambda p, c, x, n, b : S(2)**p*Int(((c*x**n)**b/((c*x**n)**(S(2)*b) + S(1)))**p, x))
    rubi.add(rule360)

    pattern361 = Pattern(Integral(csch(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, p, b: RationalQ(b, n, p)))
    rule361 = ReplacementRule(pattern361, lambda p, c, x, n, b : S(2)**p*Int(((c*x**n)**b/((c*x**n)**(S(2)*b) + S(-1)))**p, x))
    rubi.add(rule361)

    pattern362 = Pattern(Integral(sech(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, b: ZeroQ(b**S(2)*n**S(2) + S(-1))))
    rule362 = ReplacementRule(pattern362, lambda x, a, n, b, c : S(2)*Int((c*x**n)**(1/n)/((c*x**n)**(S(2)/n) + exp(-S(2)*a*b*n)), x)*exp(-a*b*n))
    rubi.add(rule362)

    pattern363 = Pattern(Integral(csch(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, b: ZeroQ(b**S(2)*n**S(2) + S(-1))))
    rule363 = ReplacementRule(pattern363, lambda x, a, n, b, c : -S(2)*b*n*Int((c*x**n)**(1/n)/(-(c*x**n)**(S(2)/n) + exp(-S(2)*a*b*n)), x)*exp(-a*b*n))
    rubi.add(rule363)

    pattern364 = Pattern(Integral(sech(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, p, b: ZeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(-1))), CustomConstraint(lambda p: NonzeroQ(p + S(-1))))
    rule364 = ReplacementRule(pattern364, lambda p, x, a, n, b, c : x*(p + S(-2))*sech(a + b*log(c*x**n))**(p + S(-2))/(p + S(-1)) + x*tanh(a + b*log(c*x**n))*sech(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))))
    rubi.add(rule364)

    pattern365 = Pattern(Integral(csch(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, p, b: ZeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(-1))), CustomConstraint(lambda p: NonzeroQ(p + S(-1))))
    rule365 = ReplacementRule(pattern365, lambda p, x, a, n, b, c : x*(-p + S(2))*csch(a + b*log(c*x**n))**(p + S(-2))/(p + S(-1)) - x*coth(a + b*log(c*x**n))*csch(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))))
    rubi.add(rule365)

    pattern366 = Pattern(Integral(sech(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p: Unequal(p, S(2))), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(-1))))
    rule366 = ReplacementRule(pattern366, lambda p, x, a, n, b, c : x*tanh(a + b*log(c*x**n))*sech(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))) + x*sech(a + b*log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))) + (b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(-1))*Int(sech(a + b*log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))))
    rubi.add(rule366)

    pattern367 = Pattern(Integral(csch(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p: Unequal(p, S(2))), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(-1))))
    rule367 = ReplacementRule(pattern367, lambda p, x, a, n, b, c : -x*coth(a + b*log(c*x**n))*csch(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))) - x*csch(a + b*log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))) - (b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(-1))*Int(csch(a + b*log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))))
    rubi.add(rule367)

    pattern368 = Pattern(Integral(sech(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(-1))))
    rule368 = ReplacementRule(pattern368, lambda p, x, a, n, b, c : b**S(2)*n**S(2)*p*(p + S(1))*Int(sech(a + b*log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*p**S(2) + S(-1)) - b*n*p*x*sinh(a + b*log(c*x**n))*sech(a + b*log(c*x**n))**(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + S(-1)) - x*sech(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(-1)))
    rubi.add(rule368)

    pattern369 = Pattern(Integral(csch(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(-1))))
    rule369 = ReplacementRule(pattern369, lambda p, x, a, n, b, c : -b**S(2)*n**S(2)*p*(p + S(1))*Int(csch(a + b*log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*p**S(2) + S(-1)) - b*n*p*x*Cosh(a + b*log(c*x**n))*csch(a + b*log(c*x**n))**(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + S(-1)) - x*csch(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(-1)))
    rubi.add(rule369)

    pattern370 = Pattern(Integral(sech(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(-1))))
    rule370 = ReplacementRule(pattern370, lambda p, x, a, n, b, c : S(2)**p*x*((c*x**n)**b*exp(a)/((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(1)))**p*((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(1))**p*Hypergeometric2F1(p, (b*n*p + S(1))/(S(2)*b*n), S(1) + (b*n*p + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*b)*exp(S(2)*a))/(b*n*p + S(1)))
    rubi.add(rule370)

    pattern371 = Pattern(Integral(csch(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) + S(-1))))
    rule371 = ReplacementRule(pattern371, lambda p, x, a, n, b, c : x*((c*x**n)**b*exp(a)/((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(-1)))**p*(-S(2)*(c*x**n)**(S(2)*b)*exp(S(2)*a) + S(2))**p*Hypergeometric2F1(p, (b*n*p + S(1))/(S(2)*b*n), S(1) + (b*n*p + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*b)*exp(S(2)*a))/(b*n*p + S(1)))
    rubi.add(rule371)

    pattern372 = Pattern(Integral(x_**WC('m', S(1))*sech(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, m, p, b: RationalQ(b, m, n, p)))
    rule372 = ReplacementRule(pattern372, lambda p, m, c, x, n, b : S(2)**p*Int(x**m*((c*x**n)**b/((c*x**n)**(S(2)*b) + S(1)))**p, x))
    rubi.add(rule372)

    pattern373 = Pattern(Integral(x_**WC('m', S(1))*csch(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, m, p, b: RationalQ(b, m, n, p)))
    rule373 = ReplacementRule(pattern373, lambda p, m, c, x, n, b : S(2)**p*Int(x**m*((c*x**n)**b/((c*x**n)**(S(2)*b) + S(-1)))**p, x))
    rubi.add(rule373)

    pattern374 = Pattern(Integral(x_**WC('m', S(1))*sec(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m, b: ZeroQ(b**S(2)*n**S(2) - (m + S(1))**S(2))))
    rule374 = ReplacementRule(pattern374, lambda m, x, a, n, b, c : S(2)*Int(x**m*(c*x**n)**((m + S(1))/n)/((c*x**n)**(S(2)*(m + S(1))/n) + exp(-S(2)*a*b*n/(m + S(1)))), x)*exp(-a*b*n/(m + S(1))))
    rubi.add(rule374)

    pattern375 = Pattern(Integral(x_**WC('m', S(1))*csc(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m, b: ZeroQ(b**S(2)*n**S(2) - (m + S(1))**S(2))))
    rule375 = ReplacementRule(pattern375, lambda m, x, a, n, b, c : -S(2)*b*n*Int(x**m*(c*x**n)**((m + S(1))/n)/(-(c*x**n)**(S(2)*(m + S(1))/n) + exp(-S(2)*a*b*n/(m + S(1)))), x)*exp(-a*b*n/(m + S(1)))/(m + S(1)))
    rubi.add(rule375)

    pattern376 = Pattern(Integral(x_**WC('m', S(1))*sech(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m, p, b: ZeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(-1))))
    rule376 = ReplacementRule(pattern376, lambda p, m, x, a, n, b, c : x**(m + S(1))*(p + S(-2))*sech(a + b*log(c*x**n))**(p + S(-2))/((m + S(1))*(p + S(-1))) + x**(m + S(1))*tanh(a + b*log(c*x**n))*sech(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))))
    rubi.add(rule376)

    pattern377 = Pattern(Integral(x_**WC('m', S(1))*csch(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m, p, b: ZeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda p: NonzeroQ(p + S(-1))))
    rule377 = ReplacementRule(pattern377, lambda p, m, x, a, n, b, c : x**(m + S(1))*(-p + S(2))*csch(a + b*log(c*x**n))**(p + S(-2))/((m + S(1))*(p + S(-1))) - x**(m + S(1))*coth(a + b*log(c*x**n))*csch(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))))
    rubi.add(rule377)

    pattern378 = Pattern(Integral(x_**WC('m', S(1))*sech(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p: Unequal(p, S(2))), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) - (m + S(1))**S(2))))
    rule378 = ReplacementRule(pattern378, lambda p, m, x, a, n, b, c : x**(m + S(1))*tanh(a + b*log(c*x**n))*sech(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))) + x**(m + S(1))*(m + S(1))*sech(a + b*log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))) + (b**S(2)*n**S(2)*(p + S(-2))**S(2) - (m + S(1))**S(2))*Int(x**m*sech(a + b*log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))))
    rubi.add(rule378)

    pattern379 = Pattern(Integral(x_**WC('m', S(1))*csch(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(1))), CustomConstraint(lambda p: Unequal(p, S(2))), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*(p + S(-2))**S(2) - (m + S(1))**S(2))))
    rule379 = ReplacementRule(pattern379, lambda p, m, x, a, n, b, c : -x**(m + S(1))*coth(a + b*log(c*x**n))*csch(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))) - x**(m + S(1))*(m + S(1))*csch(a + b*log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))) - (b**S(2)*n**S(2)*(p + S(-2))**S(2) - (m + S(1))**S(2))*Int(x**m*csch(a + b*log(c*x**n))**(p + S(-2)), x)/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))))
    rubi.add(rule379)

    pattern380 = Pattern(Integral(x_**WC('m', S(1))*sech(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2))))
    rule380 = ReplacementRule(pattern380, lambda p, m, x, a, n, b, c : b**S(2)*n**S(2)*p*(p + S(1))*Int(x**m*sech(a + b*log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)) - b*n*p*x**(m + S(1))*sinh(a + b*log(c*x**n))*sech(a + b*log(c*x**n))**(p + S(1))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)) + x**(m + S(1))*(-m + S(-1))*sech(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)))
    rubi.add(rule380)

    pattern381 = Pattern(Integral(x_**WC('m', S(1))*csch(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2))))
    rule381 = ReplacementRule(pattern381, lambda p, m, x, a, n, b, c : -b**S(2)*n**S(2)*p*(p + S(1))*Int(x**m*csch(a + b*log(c*x**n))**(p + S(2)), x)/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)) - b*n*p*x**(m + S(1))*Cosh(a + b*log(c*x**n))*csch(a + b*log(c*x**n))**(p + S(1))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)) + x**(m + S(1))*(-m + S(-1))*csch(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)))
    rubi.add(rule381)

    pattern382 = Pattern(Integral(x_**WC('m', S(1))*sech(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2))))
    rule382 = ReplacementRule(pattern382, lambda p, m, x, a, n, b, c : S(2)**p*x**(m + S(1))*((c*x**n)**b*exp(a)/((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(1)))**p*((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(1))**p*Hypergeometric2F1(p, (b*n*p + m + S(1))/(S(2)*b*n), S(1) + (b*n*p + m + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*b)*exp(S(2)*a))/(b*n*p + m + S(1)))
    rubi.add(rule382)

    pattern383 = Pattern(Integral(x_**WC('m', S(1))*csch(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m, p, b: NonzeroQ(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2))))
    rule383 = ReplacementRule(pattern383, lambda p, m, x, a, n, b, c : S(2)**p*x**(m + S(1))*((c*x**n)**b*exp(a)/((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(-1)))**p*(-(c*x**n)**(S(2)*b)*exp(S(2)*a) + S(1))**p*Hypergeometric2F1(p, (b*n*p + m + S(1))/(S(2)*b*n), S(1) + (b*n*p + m + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*b)*exp(S(2)*a))/(b*n*p + m + S(1)))
    rubi.add(rule383)

    pattern384 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*sinh(x_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule384 = ReplacementRule(pattern384, lambda p, x, a, b : -p*Int(log(b*x)**(p + S(-1))*sinh(a*x*log(b*x)**p), x) + Cosh(a*x*log(b*x)**p)/a)
    rubi.add(rule384)

    pattern385 = Pattern(Integral(Cosh(x_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1)))*log(x_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule385 = ReplacementRule(pattern385, lambda p, x, a, b : -p*Int(Cosh(a*x*log(b*x)**p)*log(b*x)**(p + S(-1)), x) + sinh(a*x*log(b*x)**p)/a)
    rubi.add(rule385)

    pattern386 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*sinh(x_**n_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule386 = ReplacementRule(pattern386, lambda p, x, a, n, b : -p*Int(log(b*x)**(p + S(-1))*sinh(a*x**n*log(b*x)**p), x)/n + x**(-n + S(1))*Cosh(a*x**n*log(b*x)**p)/(a*n) + (n + S(-1))*Int(x**(-n)*Cosh(a*x**n*log(b*x)**p), x)/(a*n))
    rubi.add(rule386)

    pattern387 = Pattern(Integral(Cosh(x_**n_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1)))*log(x_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, p: RationalQ(n, p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule387 = ReplacementRule(pattern387, lambda p, x, a, n, b : -p*Int(Cosh(a*x**n*log(b*x)**p)*log(b*x)**(p + S(-1)), x)/n + x**(-n + S(1))*sinh(a*x**n*log(b*x)**p)/(a*n) + (n + S(-1))*Int(x**(-n)*sinh(a*x**n*log(b*x)**p), x)/(a*n))
    rubi.add(rule387)

    pattern388 = Pattern(Integral(x_**WC('m', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))*sinh(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule388 = ReplacementRule(pattern388, lambda p, m, x, a, n, b : -p*Int(x**(n + S(-1))*log(b*x)**(p + S(-1))*sinh(a*x**n*log(b*x)**p), x)/n - Cosh(a*x**n*log(b*x)**p)/(a*n))
    rubi.add(rule388)

    pattern389 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1)))*log(x_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))))
    rule389 = ReplacementRule(pattern389, lambda p, m, x, a, n, b : -p*Int(x**(n + S(-1))*Cosh(a*x**n*log(b*x)**p)*log(b*x)**(p + S(-1)), x)/n + sinh(a*x**n*log(b*x)**p)/(a*n))
    rubi.add(rule389)

    pattern390 = Pattern(Integral(x_**m_*log(x_*WC('b', S(1)))**WC('p', S(1))*sinh(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m, p: RationalQ(m, n, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n, m: NonzeroQ(m - n + S(1))))
    rule390 = ReplacementRule(pattern390, lambda p, m, x, a, n, b : -p*Int(x**m*log(b*x)**(p + S(-1))*sinh(a*x**n*log(b*x)**p), x)/n + x**(m - n + S(1))*Cosh(a*x**n*log(b*x)**p)/(a*n) - (m - n + S(1))*Int(x**(m - n)*Cosh(a*x**n*log(b*x)**p), x)/(a*n))
    rubi.add(rule390)

    pattern391 = Pattern(Integral(x_**m_*Cosh(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1)))*log(x_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m, p: RationalQ(m, n, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n, m: NonzeroQ(m - n + S(1))))
    rule391 = ReplacementRule(pattern391, lambda p, m, x, a, n, b : -p*Int(x**m*Cosh(a*x**n*log(b*x)**p)*log(b*x)**(p + S(-1)), x)/n + x**(m - n + S(1))*sinh(a*x**n*log(b*x)**p)/(a*n) - (m - n + S(1))*Int(x**(m - n)*sinh(a*x**n*log(b*x)**p), x)/(a*n))
    rubi.add(rule391)

    pattern392 = Pattern(Integral(sinh(WC('a', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule392 = ReplacementRule(pattern392, lambda x, a, d, n, c : -Subst(Int(sinh(a*x)**n/x**S(2), x), x, 1/(c + d*x))/d)
    rubi.add(rule392)

    pattern393 = Pattern(Integral(Cosh(WC('a', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule393 = ReplacementRule(pattern393, lambda x, a, d, n, c : -Subst(Int(Cosh(a*x)**n/x**S(2), x), x, 1/(c + d*x))/d)
    rubi.add(rule393)

    pattern394 = Pattern(Integral(sinh((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d, c, a, b: NonzeroQ(-a*d + b*c)))
    rule394 = ReplacementRule(pattern394, lambda e, x, a, d, n, b, c : -Subst(Int(sinh(b*e/d - e*x*(-a*d + b*c)/d)**n/x**S(2), x), x, 1/(c + d*x))/d)
    rubi.add(rule394)

    pattern395 = Pattern(Integral(Cosh((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d, c, a, b: NonzeroQ(-a*d + b*c)))
    rule395 = ReplacementRule(pattern395, lambda e, x, a, d, n, b, c : -Subst(Int(Cosh(b*e/d - e*x*(-a*d + b*c)/d)**n/x**S(2), x), x, 1/(c + d*x))/d)
    rubi.add(rule395)

    pattern396 = Pattern(Integral(sinh(u_)**WC('n', S(1)), x_), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda u, x: QuotientOfLinearsQ(u, x)), )
    def With396(n, u, x):
        lst = QuotientOfLinearsParts(u, x)
        return Int(sinh((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**n, x)
    rule396 = ReplacementRule(pattern396, lambda n, u, x : With396(n, u, x))
    rubi.add(rule396)

    pattern397 = Pattern(Integral(Cosh(u_)**WC('n', S(1)), x_), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda u, x: QuotientOfLinearsQ(u, x)), )
    def With397(n, u, x):
        lst = QuotientOfLinearsParts(u, x)
        return Int(Cosh((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**n, x)
    rule397 = ReplacementRule(pattern397, lambda n, u, x : With397(n, u, x))
    rubi.add(rule397)

    pattern398 = Pattern(Integral(WC('u', S(1))*sinh(v_)**WC('p', S(1))*sinh(w_)**WC('q', S(1)), x_), CustomConstraint(lambda v, w: ZeroQ(v - w)))
    rule398 = ReplacementRule(pattern398, lambda p, x, w, v, u, q : Int(u*sinh(v)**(p + q), x))
    rubi.add(rule398)

    pattern399 = Pattern(Integral(Cosh(v_)**WC('p', S(1))*Cosh(w_)**WC('q', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda v, w: ZeroQ(v - w)))
    rule399 = ReplacementRule(pattern399, lambda p, x, w, v, u, q : Int(u*Cosh(v)**(p + q), x))
    rubi.add(rule399)

    pattern400 = Pattern(Integral(sinh(v_)**WC('p', S(1))*sinh(w_)**WC('q', S(1)), x_), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)), CustomConstraint(lambda v, x, w: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(v/w, x))))
    rule400 = ReplacementRule(pattern400, lambda p, x, w, v, q : Int(ExpandTrigReduce(sinh(v)**p*sinh(w)**q, x), x))
    rubi.add(rule400)

    pattern401 = Pattern(Integral(Cosh(v_)**WC('p', S(1))*Cosh(w_)**WC('q', S(1)), x_), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)), CustomConstraint(lambda v, x, w: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(v/w, x))))
    rule401 = ReplacementRule(pattern401, lambda p, x, w, v, q : Int(ExpandTrigReduce(Cosh(v)**p*Cosh(w)**q, x), x))
    rubi.add(rule401)

    pattern402 = Pattern(Integral(x_**WC('m', S(1))*sinh(v_)**WC('p', S(1))*sinh(w_)**WC('q', S(1)), x_), CustomConstraint(lambda p, m, q: PositiveIntegerQ(m, p, q)), CustomConstraint(lambda v, x, w: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(v/w, x))))
    rule402 = ReplacementRule(pattern402, lambda p, m, x, w, v, q : Int(ExpandTrigReduce(x**m, sinh(v)**p*sinh(w)**q, x), x))
    rubi.add(rule402)

    pattern403 = Pattern(Integral(x_**WC('m', S(1))*Cosh(v_)**WC('p', S(1))*Cosh(w_)**WC('q', S(1)), x_), CustomConstraint(lambda p, m, q: PositiveIntegerQ(m, p, q)), CustomConstraint(lambda v, x, w: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(v/w, x))))
    rule403 = ReplacementRule(pattern403, lambda p, m, x, w, v, q : Int(ExpandTrigReduce(x**m, Cosh(v)**p*Cosh(w)**q, x), x))
    rubi.add(rule403)

    pattern404 = Pattern(Integral(Cosh(w_)**WC('p', S(1))*WC('u', S(1))*sinh(v_)**WC('p', S(1)), x_), CustomConstraint(lambda v, w: ZeroQ(v - w)), CustomConstraint(lambda p: IntegerQ(p)))
    rule404 = ReplacementRule(pattern404, lambda p, x, w, v, u : S(2)**(-p)*Int(u*sinh(S(2)*v)**p, x))
    rubi.add(rule404)

    pattern405 = Pattern(Integral(Cosh(w_)**WC('q', S(1))*sinh(v_)**WC('p', S(1)), x_), CustomConstraint(lambda p, q: PositiveIntegerQ(p, q)), CustomConstraint(lambda v, x, w: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(v/w, x))))
    rule405 = ReplacementRule(pattern405, lambda p, x, w, v, q : Int(ExpandTrigReduce(Cosh(w)**q*sinh(v)**p, x), x))
    rubi.add(rule405)

    pattern406 = Pattern(Integral(x_**WC('m', S(1))*Cosh(w_)**WC('q', S(1))*sinh(v_)**WC('p', S(1)), x_), CustomConstraint(lambda p, m, q: PositiveIntegerQ(m, p, q)), CustomConstraint(lambda v, x, w: (PolynomialQ(v, x) & PolynomialQ(w, x)) | (BinomialQ(List(v, w), x) & IndependentQ(v/w, x))))
    rule406 = ReplacementRule(pattern406, lambda p, m, x, w, v, q : Int(ExpandTrigReduce(x**m, Cosh(w)**q*sinh(v)**p, x), x))
    rubi.add(rule406)

    pattern407 = Pattern(Integral(sinh(v_)*tanh(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule407 = ReplacementRule(pattern407, lambda n, v, x, w : -Cosh(v - w)*Int(tanh(w)**(n + S(-1))*sech(w), x) + Int(Cosh(v)*tanh(w)**(n + S(-1)), x))
    rubi.add(rule407)

    pattern408 = Pattern(Integral(Cosh(v_)*coth(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule408 = ReplacementRule(pattern408, lambda n, v, x, w : Cosh(v - w)*Int(coth(w)**(n + S(-1))*csch(w), x) + Int(sinh(v)*coth(w)**(n + S(-1)), x))
    rubi.add(rule408)

    pattern409 = Pattern(Integral(sinh(v_)*coth(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule409 = ReplacementRule(pattern409, lambda n, v, x, w : Int(Cosh(v)*coth(w)**(n + S(-1)), x) + Int(coth(w)**(n + S(-1))*csch(w), x)*sinh(v - w))
    rubi.add(rule409)

    pattern410 = Pattern(Integral(Cosh(v_)*tanh(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule410 = ReplacementRule(pattern410, lambda n, v, x, w : Int(sinh(v)*tanh(w)**(n + S(-1)), x) - Int(tanh(w)**(n + S(-1))*sech(w), x)*sinh(v - w))
    rubi.add(rule410)

    pattern411 = Pattern(Integral(sinh(v_)*sech(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule411 = ReplacementRule(pattern411, lambda n, v, x, w : Cosh(v - w)*Int(tanh(w)*sech(w)**(n + S(-1)), x) + Int(sech(w)**(n + S(-1)), x)*sinh(v - w))
    rubi.add(rule411)

    pattern412 = Pattern(Integral(Cosh(v_)*csch(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule412 = ReplacementRule(pattern412, lambda n, v, x, w : Cosh(v - w)*Int(coth(w)*csch(w)**(n + S(-1)), x) + Int(csch(w)**(n + S(-1)), x)*sinh(v - w))
    rubi.add(rule412)

    pattern413 = Pattern(Integral(sinh(v_)*csch(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule413 = ReplacementRule(pattern413, lambda n, v, x, w : Cosh(v - w)*Int(csch(w)**(n + S(-1)), x) + Int(coth(w)*csch(w)**(n + S(-1)), x)*sinh(v - w))
    rubi.add(rule413)

    pattern414 = Pattern(Integral(Cosh(v_)*sech(w_)**WC('n', S(1)), x_), CustomConstraint(lambda v, x: FreeQ(v, x)), CustomConstraint(lambda w: FreeQ(Mul(S(-1), w), x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda v, w: NonzeroQ(v - w)))
    rule414 = ReplacementRule(pattern414, lambda n, v, x, w : Cosh(v - w)*Int(sech(w)**(n + S(-1)), x) + Int(tanh(w)*sech(w)**(n + S(-1)), x)*sinh(v - w))
    rubi.add(rule414)

    pattern415 = Pattern(Integral((a_ + Cosh(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule415 = ReplacementRule(pattern415, lambda m, e, x, f, d, a, n, b, c : Int((a + b*sinh(S(2)*c + S(2)*d*x)/S(2))**n*(e + f*x)**m, x))
    rubi.add(rule415)

    pattern416 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, a: NonzeroQ(a - b)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: Less(n, S(0))), CustomConstraint(lambda n, m: Equal(n, S(-1)) | (Equal(m, S(1)) & Equal(n, S(-2)))))
    rule416 = ReplacementRule(pattern416, lambda m, x, a, d, n, b, c : S(2)**(-n)*Int(x**m*(S(2)*a + b*Cosh(S(2)*c + S(2)*d*x) - b)**n, x))
    rubi.add(rule416)

    pattern417 = Pattern(Integral(x_**WC('m', S(1))*(a_ + Cosh(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)))**n_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, a: NonzeroQ(a + b)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda m: Greater(m, S(0))), CustomConstraint(lambda n: Less(n, S(0))), CustomConstraint(lambda n, m: Equal(n, S(-1)) | (Equal(m, S(1)) & Equal(n, S(-2)))))
    rule417 = ReplacementRule(pattern417, lambda m, c, x, a, d, n, b : S(2)**(-n)*Int(x**m*(S(2)*a + b*Cosh(S(2)*c + S(2)*d*x) + b)**n, x))
    rubi.add(rule417)

    pattern418 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p: RationalQ(p)))
    rule418 = ReplacementRule(pattern418, lambda p, m, c, e, x, a, f, d, n, b : d**(-m + S(-1))*Subst(Int((-c*f + d*e + f*x)**m*sinh(a + b*x**n)**p, x), x, c + d*x))
    rubi.add(rule418)

    pattern419 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cosh((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda p: RationalQ(p)))
    rule419 = ReplacementRule(pattern419, lambda p, m, c, e, x, a, d, f, n, b : d**(-m + S(-1))*Subst(Int((-c*f + d*e + f*x)**m*Cosh(a + b*x**n)**p, x), x, c + d*x))
    rubi.add(rule419)

    pattern420 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(Cosh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0)) + WC('c', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda b, a: NonzeroQ(a + b)), CustomConstraint(lambda c, a: NonzeroQ(a + c)))
    rule420 = ReplacementRule(pattern420, lambda m, c, e, x, a, d, f, g, b : S(2)*Int((f + g*x)**m/(S(2)*a + b - c + (b + c)*Cosh(S(2)*d + S(2)*e*x)), x))
    rubi.add(rule420)

    pattern421 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*sech(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(b_ + WC('c', S(1))*tanh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule421 = ReplacementRule(pattern421, lambda m, e, x, f, d, g, b, c : S(2)*Int((f + g*x)**m/(b - c + (b + c)*Cosh(S(2)*d + S(2)*e*x)), x))
    rubi.add(rule421)

    pattern422 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*sech(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(WC('a', S(1))*sech(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('b', S(0)) + WC('c', S(1))*tanh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda b, a: NonzeroQ(a + b)), CustomConstraint(lambda c, a: NonzeroQ(a + c)))
    rule422 = ReplacementRule(pattern422, lambda m, c, e, x, f, d, a, g, b : S(2)*Int((f + g*x)**m/(S(2)*a + b - c + (b + c)*Cosh(S(2)*d + S(2)*e*x)), x))
    rubi.add(rule422)

    pattern423 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*csch(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(c_ + WC('b', S(1))*coth(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule423 = ReplacementRule(pattern423, lambda m, c, e, x, f, d, g, b : S(2)*Int((f + g*x)**m/(b - c + (b + c)*Cosh(S(2)*d + S(2)*e*x)), x))
    rubi.add(rule423)

    pattern424 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*csch(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(WC('a', S(1))*csch(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('b', S(1))*coth(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda b, a: NonzeroQ(a + b)), CustomConstraint(lambda c, a: NonzeroQ(a + c)))
    rule424 = ReplacementRule(pattern424, lambda m, e, x, a, d, f, g, b, c : S(2)*Int((f + g*x)**m/(S(2)*a + b - c + (b + c)*Cosh(S(2)*d + S(2)*e*x)), x))
    rubi.add(rule424)

    pattern425 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cosh(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule425 = ReplacementRule(pattern425, lambda m, e, x, f, d, a, b, c : Int((e + f*x)**m*exp(c + d*x)/(a + b*exp(c + d*x) - Rt(a**S(2) + b**S(2), S(2))), x) + Int((e + f*x)**m*exp(c + d*x)/(a + b*exp(c + d*x) + Rt(a**S(2) + b**S(2), S(2))), x) - (e + f*x)**(m + S(1))/(b*f*(m + S(1))))
    rubi.add(rule425)

    pattern426 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + Cosh(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)))
    rule426 = ReplacementRule(pattern426, lambda m, c, e, x, f, d, a, b : Int((e + f*x)**m*exp(c + d*x)/(a + b*exp(c + d*x) - Rt(a**S(2) - b**S(2), S(2))), x) + Int((e + f*x)**m*exp(c + d*x)/(a + b*exp(c + d*x) + Rt(a**S(2) - b**S(2), S(2))), x) - (e + f*x)**(m + S(1))/(b*f*(m + S(1))))
    rubi.add(rule426)

    pattern427 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cosh(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda b, a: ZeroQ(a**S(2) + b**S(2))))
    rule427 = ReplacementRule(pattern427, lambda m, e, x, f, d, a, n, b, c : Int((e + f*x)**m*Cosh(c + d*x)**(n + S(-2))*sinh(c + d*x), x)/b + Int((e + f*x)**m*Cosh(c + d*x)**(n + S(-2)), x)/a)
    rubi.add(rule427)

    pattern428 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + Cosh(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))))
    rule428 = ReplacementRule(pattern428, lambda m, c, e, x, f, d, a, n, b : Int((e + f*x)**m*Cosh(c + d*x)*sinh(c + d*x)**(n + S(-2)), x)/b + Int((e + f*x)**m*sinh(c + d*x)**(n + S(-2)), x)/a)
    rubi.add(rule428)

    pattern429 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*Cosh(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) + b**S(2))))
    rule429 = ReplacementRule(pattern429, lambda m, e, x, f, d, a, n, b, c : -a*Int((e + f*x)**m*Cosh(c + d*x)**(n + S(-2)), x)/b**S(2) + Int((e + f*x)**m*Cosh(c + d*x)**(n + S(-2))*sinh(c + d*x), x)/b + (a**S(2) + b**S(2))*Int((e + f*x)**m*Cosh(c + d*x)**(n + S(-2))/(a + b*sinh(c + d*x)), x)/b**S(2))
    rubi.add(rule429)

    pattern430 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + Cosh(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m: PositiveIntegerQ(m)), CustomConstraint(lambda n: IntegerQ(n)), CustomConstraint(lambda n: Greater(n, S(1))), CustomConstraint(lambda b, a: NonzeroQ(a**S(2) - b**S(2))))
    rule430 = ReplacementRule(pattern430, lambda m, c, e, x, f, d, a, n, b : -a*Int((e + f*x)**m*sinh(c + d*x)**(n + S(-2)), x)/b**S(2) + Int((e + f*x)**m*Cosh(c + d*x)*sinh(c + d*x)**(n + S(-2)), x)/b + (a**S(2) - b**S(2))*Int((e + f*x)**m*sinh(c + d*x)**(n + S(-2))/(a + b*Cosh(c + d*x)), x)/b**S(2))
    rubi.add(rule430)

    pattern431 = Pattern(Integral((A_ + WC('B', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))))*(x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda A, B, b, a: ZeroQ(A*a + B*b)))
    rule431 = ReplacementRule(pattern431, lambda A, e, x, f, d, a, b, B, c : -B*f*Int(Cosh(c + d*x)/(a + b*sinh(c + d*x)), x)/(a*d) + B*(e + f*x)*Cosh(c + d*x)/(a*d*(a + b*sinh(c + d*x))))
    rubi.add(rule431)

    pattern432 = Pattern(Integral((A_ + Cosh(x_*WC('d', S(1)) + WC('c', S(0)))*WC('B', S(1)))*(x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + Cosh(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**S(2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda A, B, b, a: ZeroQ(A*a - B*b)))
    rule432 = ReplacementRule(pattern432, lambda A, c, e, x, f, d, a, B, b : -B*f*Int(sinh(c + d*x)/(a + b*Cosh(c + d*x)), x)/(a*d) + B*(e + f*x)*sinh(c + d*x)/(a*d*(a + b*Cosh(c + d*x))))
    rubi.add(rule432)

    pattern433 = Pattern(Integral((a_ + WC('b', S(1))*tanh(v_))**WC('n', S(1))*sech(v_)**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Equal(m + n, S(0))), CustomConstraint(lambda m: OddQ(m)))
    rule433 = ReplacementRule(pattern433, lambda m, x, a, n, v, b : Int((a*Cosh(v) + b*sinh(v))**n, x))
    rubi.add(rule433)

    pattern434 = Pattern(Integral((a_ + WC('b', S(1))*coth(v_))**WC('n', S(1))*csch(v_)**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, m: IntegersQ(m, n)), CustomConstraint(lambda n, m: Equal(m + n, S(0))), CustomConstraint(lambda m: OddQ(m)))
    rule434 = ReplacementRule(pattern434, lambda m, x, a, n, v, b : Int((a*sinh(v) + b*Cosh(v))**n, x))
    rubi.add(rule434)

    pattern435 = Pattern(Integral(WC('u', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule435 = ReplacementRule(pattern435, lambda m, x, a, d, n, u, b, c : Int(ExpandTrigReduce(u, sinh(a + b*x)**m*sinh(c + d*x)**n, x), x))
    rubi.add(rule435)

    pattern436 = Pattern(Integral(Cosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*Cosh(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, m: PositiveIntegerQ(m, n)))
    rule436 = ReplacementRule(pattern436, lambda m, x, a, d, n, b, u, c : Int(ExpandTrigReduce(u, Cosh(a + b*x)**m*Cosh(c + d*x)**n, x), x))
    rubi.add(rule436)

    pattern437 = Pattern(Integral(sech(c_ + x_*WC('d', S(1)))*sech(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, b: ZeroQ(b**S(2) - d**S(2))), CustomConstraint(lambda d, c, a, b: NonzeroQ(-a*d + b*c)))
    rule437 = ReplacementRule(pattern437, lambda x, a, d, b, c : -Int(tanh(a + b*x), x)*csch((-a*d + b*c)/d) + Int(tanh(c + d*x), x)*csch((-a*d + b*c)/b))
    rubi.add(rule437)

    pattern438 = Pattern(Integral(csch(c_ + x_*WC('d', S(1)))*csch(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, b: ZeroQ(b**S(2) - d**S(2))), CustomConstraint(lambda d, c, a, b: NonzeroQ(-a*d + b*c)))
    rule438 = ReplacementRule(pattern438, lambda x, a, d, b, c : Int(coth(a + b*x), x)*csch((-a*d + b*c)/b) - Int(coth(c + d*x), x)*csch((-a*d + b*c)/d))
    rubi.add(rule438)

    pattern439 = Pattern(Integral(tanh(c_ + x_*WC('d', S(1)))*tanh(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, b: ZeroQ(b**S(2) - d**S(2))), CustomConstraint(lambda d, c, a, b: NonzeroQ(-a*d + b*c)))
    rule439 = ReplacementRule(pattern439, lambda x, a, d, b, c : b*x/d - b*Cosh((-a*d + b*c)/d)*Int(sech(a + b*x)*sech(c + d*x), x)/d)
    rubi.add(rule439)

    pattern440 = Pattern(Integral(coth(c_ + x_*WC('d', S(1)))*coth(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda d, b: ZeroQ(b**S(2) - d**S(2))), CustomConstraint(lambda d, c, a, b: NonzeroQ(-a*d + b*c)))
    rule440 = ReplacementRule(pattern440, lambda x, a, d, b, c : b*x/d + Cosh((-a*d + b*c)/d)*Int(csch(a + b*x)*csch(c + d*x), x))
    rubi.add(rule440)

    pattern441 = Pattern(Integral((Cosh(v_)*WC('a', S(1)) + WC('b', S(1))*sinh(v_))**WC('n', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, a: ZeroQ(a**S(2) - b**S(2))))
    rule441 = ReplacementRule(pattern441, lambda x, a, n, v, u, b : Int(u*(a*exp(a*v/b))**n, x))
    rubi.add(rule441)

    return rubi
