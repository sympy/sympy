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

def miscellaneous_integration(rubi):

    pattern1 = Pattern(Integral(u_*((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, x, u: Not(MatchQ(u, Condition(x**Optional(Pattern(n1, Blank))*Optional(Pattern(v, Blank)), ZeroQ(n - n1 + S(-1)))))))
    rule1 = ReplacementRule(pattern1, lambda u, a, c, x, p, n, b : c**IntPart(p)*(c*(a + b*x)**n)**FracPart(p)*(a + b*x)**(-n*FracPart(p))*Int(u*(a + b*x)**(n*p), x))
    rubi.add(rule1)

    pattern2 = Pattern(Integral(((d_*(x_*WC('b', S(1)) + WC('a', S(0))))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda q: Not(IntegerQ(q))))
    rule2 = ReplacementRule(pattern2, lambda u, a, d, c, q, x, p, b : (c*(d*(a + b*x))**p)**q*(a + b*x)**(-p*q)*Int(u*(a + b*x)**(p*q), x))
    rubi.add(rule2)

    pattern3 = Pattern(Integral((((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('d', S(1)))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda q: Not(IntegerQ(q))))
    rule3 = ReplacementRule(pattern3, lambda u, a, d, c, q, x, p, n, b : (c*(d*(a + b*x)**n)**p)**q*(a + b*x)**(-n*p*q)*Int(u*(a + b*x)**(n*p*q), x))
    rubi.add(rule3)

    pattern4 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda g, e: ZeroQ(e + g)), CustomConstraint(lambda f, d: ZeroQ(d + f + S(-2))), CustomConstraint(lambda d, C, e, A, f: ZeroQ(A*e**S(2) + C*d*f)), CustomConstraint(lambda B, d, C, e: ZeroQ(-B*e + S(2)*C*(d + S(-1)))))
    rule4 = ReplacementRule(pattern4, lambda B, a, g, d, C, c, e, F, n, A, x, f, b : g*Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x))/C)
    rubi.add(rule4)

    pattern5 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + S(1))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda g, e: ZeroQ(e + g)), CustomConstraint(lambda A, C, e: ZeroQ(A*e**S(2) + C)))
    rule5 = ReplacementRule(pattern5, lambda a, g, c, C, e, F, A, x, n, b : g*Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1)))/C)
    rubi.add(rule5)

    pattern6 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda g, e: ZeroQ(e + g)), CustomConstraint(lambda f, d: ZeroQ(d + f + S(-2))), CustomConstraint(lambda d, C, e, A, f: ZeroQ(A*e**S(2) + C*d*f)), CustomConstraint(lambda B, d, C, e: ZeroQ(-B*e + S(2)*C*(d + S(-1)))))
    rule6 = ReplacementRule(pattern6, lambda B, a, g, d, c, e, C, F, n, A, x, f, b : g*Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x))/C)
    rubi.add(rule6)

    pattern7 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda g, e: ZeroQ(e + g)), CustomConstraint(lambda A, C, e: ZeroQ(A*e**S(2) + C)))
    rule7 = ReplacementRule(pattern7, lambda a, g, c, e, C, F, A, x, n, b : g*Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1)))/C)
    rubi.add(rule7)

    pattern8 = Pattern(Integral(u_/y_, x_), CustomConstraint(lambda q, y, x: Not(FalseQ(q))))
    def With8(y, x, u):
        q = DerivativeDivides(y, u, x)
        return q*log(RemoveContent(y, x))
    rule8 = ReplacementRule(pattern8, lambda y, x, u : With8(y, x, u))
    rubi.add(rule8)

    pattern9 = Pattern(Integral(u_/(w_*y_), x_), CustomConstraint(lambda w, x, y, q: Not(FalseQ(q))))
    def With9(w, y, x, u):
        q = DerivativeDivides(w*y, u, x)
        return q*log(RemoveContent(w*y, x))
    rule9 = ReplacementRule(pattern9, lambda w, y, x, u : With9(w, y, x, u))
    rubi.add(rule9)

    pattern10 = Pattern(Integral(u_*y_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda q, y, m: Not(FalseQ(q))))
    def With10(m, y, x, u):
        q = DerivativeDivides(y, u, x)
        return q*y**(m + S(1))/(m + S(1))
    rule10 = ReplacementRule(pattern10, lambda m, y, x, u : With10(m, y, x, u))
    rubi.add(rule10)

    pattern11 = Pattern(Integral(u_*y_**WC('m', S(1))*z_**WC('n', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda m, y, z, q: Not(FalseQ(q))))
    def With11(u, m, z, x, n, y):
        q = DerivativeDivides(y*z, u*z**(-m + n), x)
        return q*y**(m + S(1))*z**(m + S(1))/(m + S(1))
    rule11 = ReplacementRule(pattern11, lambda u, m, z, x, n, y : With11(u, m, z, x, n, y))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(u_, x_), CustomConstraint(lambda x, v, u: SimplerIntegrandQ(v, u, x)))
    def With12(x, u):
        v = SimplifyIntegrand(u, x)
        return Int(v, x)
    rule12 = ReplacementRule(pattern12, lambda x, u : With12(x, u))
    #rubi.add(rule12)

    pattern13 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda f, b, d, e: ZeroQ(b*e**S(2) - d*f**S(2))))
    rule13 = ReplacementRule(pattern13, lambda u, a, m, d, c, e, x, f, n, b : (a*e**S(2) - c*f**S(2))**m*Int(ExpandIntegrand(u*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x))
    rubi.add(rule13)

    pattern14 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda f, c, e, a: ZeroQ(a*e**S(2) - c*f**S(2))))
    rule14 = ReplacementRule(pattern14, lambda u, a, m, d, c, e, x, f, n, b : (b*e**S(2) - d*f**S(2))**m*Int(ExpandIntegrand(u*x**(m*n)*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x))
    rubi.add(rule14)

    pattern15 = Pattern(Integral(u_**WC('m', S(1))*w_*(u_**n_*WC('a', S(1)) + v_)**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda n: Not(PositiveQ(n))), CustomConstraint(lambda x, v: NFreeQ(v, x)))
    rule15 = ReplacementRule(pattern15, lambda w, u, a, m, v, x, p, n : Int(u**(m + n*p)*w*(a + u**(-n)*v)**p, x))
    rubi.add(rule15)

    pattern16 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda y, v: ZeroQ(-v + y)), CustomConstraint(lambda b, n, y, q, c, m, x, a, d: Not(FalseQ(q))))
    def With16(u, a, m, v, d, c, x, n, b, y):
        q = DerivativeDivides(y, u, x)
        return q*Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, y)
    rule16 = ReplacementRule(pattern16, lambda u, a, m, v, d, c, x, n, b, y : With16(u, a, m, v, d, c, x, n, b, y))
    rubi.add(rule16)

    pattern17 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda y, v: ZeroQ(-v + y)), CustomConstraint(lambda w, y: ZeroQ(-w + y)), CustomConstraint(lambda b, n, y, q, e, p, c, m, x, a, d, f: Not(FalseQ(q))))
    def With17(w, u, a, m, v, d, c, e, n, x, p, f, b, y):
        q = DerivativeDivides(y, u, x)
        return q*Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x, y)
    rule17 = ReplacementRule(pattern17, lambda w, u, a, m, v, d, c, e, n, x, p, f, b, y : With17(w, u, a, m, v, d, c, e, n, x, p, f, b, y))
    rubi.add(rule17)

    pattern18 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(z_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda y, v: ZeroQ(-v + y)), CustomConstraint(lambda w, y: ZeroQ(-w + y)), CustomConstraint(lambda z, y: ZeroQ(y - z)), CustomConstraint(lambda b, n, y, q, e, p, g, c, m, x, a, r, d, f, h: Not(FalseQ(r))))
    def With18(w, u, a, m, g, v, d, c, e, q, h, z, n, x, p, f, b, y):
        r = DerivativeDivides(y, u, x)
        return r*Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x, y)
    rule18 = ReplacementRule(pattern18, lambda w, u, a, m, g, v, d, c, e, q, h, z, n, x, p, f, b, y : With18(w, u, a, m, g, v, d, c, e, q, h, z, n, x, p, f, b, y))
    rubi.add(rule18)

    pattern19 = Pattern(Integral((a_ + y_**n_*WC('b', S(1)))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, x, a, y, n, u, q: Not(FalseQ(q))))
    def With19(u, a, n, x, y, b):
        q = DerivativeDivides(y, u, x)
        return a*Int(u, x) + b*q*Subst(Int(x**n, x), x, y)
    rule19 = ReplacementRule(pattern19, lambda u, a, n, x, y, b : With19(u, a, n, x, y, b))
    rubi.add(rule19)

    pattern20 = Pattern(Integral((y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p, b, x, a, n, y, q: Not(FalseQ(q))))
    def With20(u, a, n, x, p, y, b):
        q = DerivativeDivides(y, u, x)
        return q*Subst(Int((a + b*x**n)**p, x), x, y)
    rule20 = ReplacementRule(pattern20, lambda u, a, n, x, p, y, b : With20(u, a, n, x, p, y, b))
    rubi.add(rule20)

    pattern21 = Pattern(Integral(v_**WC('m', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, n, y, v, q, p, m, x, a, r, u: Not(FalseQ(Set(q, DerivativeDivides(y, u, x)))) & Not(FalseQ(Set(r, Divides(y**m, v**m, x))))))
    def With21(u, a, m, v, n, x, p, y, b):
        q = Symbol('q')
        r = Symbol('r')
        return q*r*Subst(Int(x**m*(a + b*x**n)**p, x), x, y)
    rule21 = ReplacementRule(pattern21, lambda u, a, m, v, n, x, p, y, b : With21(u, a, m, v, n, x, p, y, b))
    rubi.add(rule21)

    pattern22 = Pattern(Integral((v_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, n2: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda y, v: ZeroQ(-v + y)), CustomConstraint(lambda b, n, y, q, p, c, x, a: Not(FalseQ(q))))
    def With22(n2, u, a, v, c, n, x, p, y, b):
        q = DerivativeDivides(y, u, x)
        return q*Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, y)
    rule22 = ReplacementRule(pattern22, lambda n2, u, a, v, c, n, x, p, y, b : With22(n2, u, a, v, c, n, x, p, y, b))
    rubi.add(rule22)

    pattern23 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, n2: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda y, v: ZeroQ(-v + y)), CustomConstraint(lambda w, y: ZeroQ(-w + y)), CustomConstraint(lambda b, n, y, q, A, p, c, x, a, B: Not(FalseQ(q))))
    def With23(w, B, n2, u, a, v, c, A, x, p, n, b, y):
        q = DerivativeDivides(y, u, x)
        return q*Subst(Int((A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y)
    rule23 = ReplacementRule(pattern23, lambda w, B, n2, u, a, v, c, A, x, p, n, b, y : With23(w, B, n2, u, a, v, c, A, x, p, n, b, y))
    rubi.add(rule23)

    pattern24 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, n2: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda w, y: ZeroQ(-w + y)), CustomConstraint(lambda n, y, q, A, p, c, x, a, B: Not(FalseQ(q))))
    def With24(w, B, n2, u, a, c, n, A, x, p, y):
        q = DerivativeDivides(y, u, x)
        return q*Subst(Int((A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y)
    rule24 = ReplacementRule(pattern24, lambda w, B, n2, u, a, c, n, A, x, p, y : With24(w, B, n2, u, a, c, n, A, x, p, y))
    rubi.add(rule24)

    pattern25 = Pattern(Integral(v_**WC('m', S(1))*(w_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, n2: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda w, y: ZeroQ(-w + y)), CustomConstraint(lambda b, n, y, v, q, p, c, m, x, a, r, u: Not(FalseQ(Set(q, DerivativeDivides(y, u, x)))) & Not(FalseQ(Set(r, Divides(y**m, v**m, x))))))
    def With25(w, n2, u, a, m, v, c, n, x, p, y, b):
        q = Symbol('q')
        r = Symbol('r')
        return q*r*Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y)
    rule25 = ReplacementRule(pattern25, lambda w, n2, u, a, m, v, c, n, x, p, y, b : With25(w, n2, u, a, m, v, c, n, x, p, y, b))
    rubi.add(rule25)

    pattern26 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, n2: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda y, v: ZeroQ(-v + y)), CustomConstraint(lambda w, y: ZeroQ(-w + y)), CustomConstraint(lambda b, n, y, q, A, p, c, m, x, a, r, u, z, B: Not(FalseQ(Set(q, DerivativeDivides(y, u, x)))) & Not(FalseQ(Set(r, Divides(y**m, z**m, x))))))
    def With26(w, B, n2, u, a, m, v, c, z, A, x, p, n, b, y):
        q = Symbol('q')
        r = Symbol('r')
        return q*r*Subst(Int(x**m*(A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y)
    rule26 = ReplacementRule(pattern26, lambda w, B, n2, u, a, m, v, c, z, A, x, p, n, b, y : With26(w, B, n2, u, a, m, v, c, z, A, x, p, n, b, y))
    rubi.add(rule26)

    pattern27 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, n2: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda w, y: ZeroQ(-w + y)), CustomConstraint(lambda n, y, q, A, p, c, m, x, a, r, u, z, B: Not(FalseQ(Set(q, DerivativeDivides(y, u, x)))) & Not(FalseQ(Set(r, Divides(y**m, z**m, x))))))
    def With27(w, B, n2, u, a, m, c, z, n, A, x, p, y):
        q = Symbol('q')
        r = Symbol('r')
        return q*r*Subst(Int(x**m*(A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y)
    rule27 = ReplacementRule(pattern27, lambda w, B, n2, u, a, m, c, z, n, A, x, p, y : With27(w, B, n2, u, a, m, c, z, n, A, x, p, y))
    rubi.add(rule27)

    pattern28 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda y, v: ZeroQ(-v + y)), CustomConstraint(lambda b, n, y, q, p, c, m, x, a, d: Not(FalseQ(q))))
    def With28(u, a, m, v, d, c, x, p, n, b, y):
        q = DerivativeDivides(y, u, x)
        return q*Subst(Int((a + b*x**n)**m*(c + d*x**n)**p, x), x, y)
    rule28 = ReplacementRule(pattern28, lambda u, a, m, v, d, c, x, p, n, b, y : With28(u, a, m, v, d, c, x, p, n, b, y))
    rubi.add(rule28)

    pattern29 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('q', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda y, v: ZeroQ(-v + y)), CustomConstraint(lambda w, y: ZeroQ(-w + y)), CustomConstraint(lambda b, n, y, q, e, p, c, m, x, a, r, d, f: Not(FalseQ(r))))
    def With29(w, u, a, m, v, d, c, e, q, n, x, p, f, b, y):
        r = DerivativeDivides(y, u, x)
        return r*Subst(Int((a + b*x**n)**m*(c + d*x**n)**p*(e + f*x**n)**q, x), x, y)
    rule29 = ReplacementRule(pattern29, lambda w, u, a, m, v, d, c, e, q, n, x, p, f, b, y : With29(w, u, a, m, v, d, c, e, q, n, x, p, f, b, y))
    rubi.add(rule29)

    pattern30 = Pattern(Integral(F_**v_*u_, x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda q, F, v: Not(FalseQ(q))))
    def With30(F, u, x, v):
        q = DerivativeDivides(v, u, x)
        return F**v*q/log(F)
    rule30 = ReplacementRule(pattern30, lambda F, u, x, v : With30(F, u, x, v))
    rubi.add(rule30)

    pattern31 = Pattern(Integral(F_**v_*u_*w_**WC('m', S(1)), x_), CustomConstraint(lambda F, x: FreeQ(F, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda w, v: ZeroQ(-v + w)), CustomConstraint(lambda F, m, x, v, q: Not(FalseQ(q))))
    def With31(w, u, v, m, F, x):
        q = DerivativeDivides(v, u, x)
        return q*Subst(Int(F**x*x**m, x), x, v)
    rule31 = ReplacementRule(pattern31, lambda w, u, v, m, F, x : With31(w, u, v, m, F, x))
    rubi.add(rule31)

    pattern32 = Pattern(Integral(u_*(a_ + v_**WC('p', S(1))*w_**WC('p', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: IntegerQ(p)), )
    def With32(w, u, a, m, v, x, p, b):
        c = u/(v*D(w, x) + w*D(v, x))
        return c*Subst(Int((a + b*x**p)**m, x), x, v*w)
    rule32 = ReplacementRule(pattern32, lambda w, u, a, m, v, x, p, b : With32(w, u, a, m, v, x, p, b))
    rubi.add(rule32)

    pattern33 = Pattern(Integral(u_*v_**WC('r', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda q, r, p: ZeroQ(p - q*(r + S(1)))), CustomConstraint(lambda r: NonzeroQ(r + S(1))), CustomConstraint(lambda r, p: IntegerQ(p/(r + S(1)))), )
    def With33(w, u, v, m, a, r, q, x, p, b):
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        return c*p*Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w)/(r + S(1))
    rule33 = ReplacementRule(pattern33, lambda w, u, v, m, a, r, q, x, p, b : With33(w, u, v, m, a, r, q, x, p, b))
    rubi.add(rule33)

    pattern34 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda q, s, r, p: ZeroQ(p*(s + S(1)) - q*(r + S(1)))), CustomConstraint(lambda r: NonzeroQ(r + S(1))), CustomConstraint(lambda r, p: IntegerQ(p/(r + S(1)))), )
    def With34(w, u, v, m, a, s, r, q, x, p, b):
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        return c*p*Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w**(s + S(1)))/(r + S(1))
    rule34 = ReplacementRule(pattern34, lambda w, u, v, m, a, s, r, q, x, p, b : With34(w, u, v, m, a, s, r, q, x, p, b))
    rubi.add(rule34)

    pattern35 = Pattern(Integral(u_*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda m, q, p: ZeroQ(p + q*(m*p + S(1)))), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda m: IntegerQ(m)), )
    def With35(w, u, a, m, v, q, x, p, b):
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        return c*p*Subst(Int((a*x**p + b)**m, x), x, v*w**(m*q + S(1)))
    rule35 = ReplacementRule(pattern35, lambda w, u, a, m, v, q, x, p, b : With35(w, u, a, m, v, q, x, p, b))
    rubi.add(rule35)

    pattern36 = Pattern(Integral(u_*v_**WC('r', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda m, q, r, p: ZeroQ(p + q*(m*p + r + S(1)))), CustomConstraint(lambda q: IntegerQ(q)), CustomConstraint(lambda m: IntegerQ(m)), )
    def With36(w, u, a, m, v, r, q, x, p, b):
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        return -c*q*Subst(Int((a + b*x**q)**m, x), x, v**(m*p + r + S(1))*w)
    rule36 = ReplacementRule(pattern36, lambda w, u, a, m, v, r, q, x, p, b : With36(w, u, a, m, v, r, q, x, p, b))
    rubi.add(rule36)

    pattern37 = Pattern(Integral(u_*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda m, q, s, p: ZeroQ(p*(s + S(1)) + q*(m*p + S(1)))), CustomConstraint(lambda s: NonzeroQ(s + S(1))), CustomConstraint(lambda s, q: IntegerQ(q/(s + S(1)))), CustomConstraint(lambda m: IntegerQ(m)), )
    def With37(w, u, a, m, v, s, q, x, p, b):
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        return -c*q*Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + S(1))*w**(s + S(1)))/(s + S(1))
    rule37 = ReplacementRule(pattern37, lambda w, u, a, m, v, s, q, x, p, b : With37(w, u, a, m, v, s, q, x, p, b))
    rubi.add(rule37)

    pattern38 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda m, s, r, q, p: ZeroQ(p*(s + S(1)) + q*(m*p + r + S(1)))), CustomConstraint(lambda s: NonzeroQ(s + S(1))), CustomConstraint(lambda s, q: IntegerQ(q/(s + S(1)))), CustomConstraint(lambda m: IntegerQ(m)), )
    def With38(w, u, a, m, v, s, r, q, x, p, b):
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        return -c*q*Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + r + S(1))*w**(s + S(1)))/(s + S(1))
    rule38 = ReplacementRule(pattern38, lambda w, u, a, m, v, s, r, q, x, p, b : With38(w, u, a, m, v, s, r, q, x, p, b))
    rubi.add(rule38)

    pattern39 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda m, x, u: FunctionOfQ(x**(m + S(1)), u, x)))
    rule39 = ReplacementRule(pattern39, lambda m, x, u : Subst(Int(SubstFor(x**(m + S(1)), u, x), x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule39)

    pattern40 = Pattern(Integral(u_, x_), CustomConstraint(lambda x, u, lst: Not(FalseQ(lst)) & SubstForFractionalPowerQ(u, Part(lst, S(3)), x)))
    def With40(x, u):
        lst = SubstForFractionalPowerOfLinear(u, x)
        return Part(lst, S(2))*Part(lst, S(4))*Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(1/Part(lst, S(2))))
    rule40 = ReplacementRule(pattern40, lambda x, u : With40(x, u))
    #rubi.add(rule40)

    pattern41 = Pattern(Integral(u_, x_), CustomConstraint(lambda lst, x: Not(FalseQ(lst))))
    def With41(x, u):
        lst = SubstForFractionalPowerOfQuotientOfLinears(u, x)
        return Part(lst, S(2))*Part(lst, S(4))*Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(1/Part(lst, S(2))))
    rule41 = ReplacementRule(pattern41, lambda x, u : With41(x, u))
    #rubi.add(rule41)

    pattern42 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*z_**WC('q', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda x, v: NFreeQ(v, x)), CustomConstraint(lambda w, x: NFreeQ(w, x)), CustomConstraint(lambda z, x: NFreeQ(z, x)))
    rule42 = ReplacementRule(pattern42, lambda w, u, a, m, v, q, z, x, p, n : a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*z**(-q*FracPart(p))*(a*v**m*w**n*z**q)**FracPart(p)*Int(u*v**(m*p)*w**(n*p)*z**(p*q), x))
    rubi.add(rule42)

    pattern43 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda x, v: NFreeQ(v, x)), CustomConstraint(lambda w, x: NFreeQ(w, x)))
    rule43 = ReplacementRule(pattern43, lambda w, v, u, m, a, x, p, n : a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*(a*v**m*w**n)**FracPart(p)*Int(u*v**(m*p)*w**(n*p), x))
    rubi.add(rule43)

    pattern44 = Pattern(Integral((v_**WC('m', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda x, v: NFreeQ(v, x)), CustomConstraint(lambda m, a: Not(EqQ(a, S(1)) & EqQ(m, S(1)))), CustomConstraint(lambda m, x, v: Not(EqQ(m, S(1)) & EqQ(v, x))))
    rule44 = ReplacementRule(pattern44, lambda v, u, m, a, x, p : a**IntPart(p)*v**(-m*FracPart(p))*(a*v**m)**FracPart(p)*Int(u*v**(m*p), x))
    rubi.add(rule44)

    pattern45 = Pattern(Integral((x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda x, u: Not(RationalFunctionQ(u, x))))
    rule45 = ReplacementRule(pattern45, lambda u, a, x, p, n, b : FullSimplify(x**(-n/S(2))*sqrt(a + b*x**n)/sqrt(a*x**(-n) + b))*Int(u*x**(n*p)*(a*x**(-n) + b)**p, x))
    rubi.add(rule45)

    pattern46 = Pattern(Integral((v_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda x, v: BinomialQ(v, x)), CustomConstraint(lambda x, v: Not(LinearQ(v, x))))
    rule46 = ReplacementRule(pattern46, lambda v, u, a, x, p, n, b : v**(-n*FracPart(p))*(a + b*v**n)**FracPart(p)*(a*v**(-n) + b)**(-FracPart(p))*Int(u*v**(n*p)*(a*v**(-n) + b)**p, x))
    rubi.add(rule46)

    pattern47 = Pattern(Integral((v_**n_*x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda x, v: BinomialQ(v, x)))
    rule47 = ReplacementRule(pattern47, lambda v, u, m, a, x, p, n, b : v**(-n*FracPart(p))*(a + b*v**n*x**m)**FracPart(p)*(a*v**(-n) + b*x**m)**(-FracPart(p))*Int(u*v**(n*p)*(a*v**(-n) + b*x**m)**p, x))
    rubi.add(rule47)

    pattern48 = Pattern(Integral((x_**WC('r', S(1))*WC('a', S(1)) + x_**WC('s', S(1))*WC('b', S(1)))**m_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)), CustomConstraint(lambda m: Not(IntegerQ(m))), CustomConstraint(lambda s, r: PosQ(-r + s)), CustomConstraint(lambda b, v, s, m, x, a, r, u: Not(EqQ(v, S(1)))))
    def With48(u, a, m, s, r, x, b):
        v = x**(-r*FracPart(m))*(a + b*x**(-r + s))**(-FracPart(m))*(a*x**r + b*x**s)**FracPart(m)
        return v*Int(u*x**(m*r)*(a + b*x**(-r + s))**m, x)
    rule48 = ReplacementRule(pattern48, lambda u, a, m, s, r, x, b : With48(u, a, m, s, r, x, b))
    rubi.add(rule48)

    pattern49 = Pattern(Integral(u_/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda v, x: SumQ(v)))
    def With49(u, a, x, n, b):
        v = RationalFunctionExpand(u/(a + b*x**n), x)
        return Int(v, x)
    rule49 = ReplacementRule(pattern49, lambda u, a, x, n, b : With49(u, a, x, n, b))
    rubi.add(rule49)

    pattern50 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, n2: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda b, c, a: ZeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda x, u: Not(AlgebraicFunctionQ(u, x))))
    rule50 = ReplacementRule(pattern50, lambda n2, u, a, c, x, p, n, b : S(4)**(-p)*c**(-p)*Int(u*(b + S(2)*c*x**n)**(S(2)*p), x))
    rubi.add(rule50)

    pattern51 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, n2: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda b, c, a: ZeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda x, u: Not(AlgebraicFunctionQ(u, x))))
    rule51 = ReplacementRule(pattern51, lambda n2, u, a, c, x, p, n, b : (b + S(2)*c*x**n)**(-S(2)*p)*(a + b*x**n + c*x**(S(2)*n))**p*Int(u*(b + S(2)*c*x**n)**(S(2)*p), x))
    rubi.add(rule51)

    pattern52 = Pattern(Integral(u_/(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, n2: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda v, x: SumQ(v)))
    def With52(n2, u, a, c, x, n, b):
        v = RationalFunctionExpand(u/(a + b*x**n + c*x**(S(2)*n)), x)
        return Int(v, x)
    rule52 = ReplacementRule(pattern52, lambda n2, u, a, c, x, n, b : With52(n2, u, a, c, x, n, b))
    rubi.add(rule52)

    pattern53 = Pattern(Integral(WC('u', S(1))/(x_**WC('m', S(1))*WC('a', S(1)) + sqrt(x_**n_*WC('c', S(1)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)))
    rule53 = ReplacementRule(pattern53, lambda u, a, m, c, x, n, b : Int(u*(a*x**m - b*sqrt(c*x**n))/(a**S(2)*x**(S(2)*m) - b**S(2)*c*x**n), x))
    rubi.add(rule53)

    pattern54 = Pattern(Integral(u_, x_), CustomConstraint(lambda lst, x: Not(FalseQ(lst))))
    def With54(x, u):
        lst = FunctionOfLinear(u, x)
        return Dist(1/Part(lst, S(3)), Subst(Int(Part(lst, S(1)), x), x, x*Part(lst, S(3)) + Part(lst, S(2))), x)
    rule54 = ReplacementRule(pattern54, lambda x, u : With54(x, u))
    #rubi.add(rule54)

    pattern55 = Pattern(Integral(u_/x_, x_), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda x, u: Not(RationalFunctionQ(u, x))), CustomConstraint(lambda lst, x: Not(FalseQ(lst)) & NonzeroQ(Part(lst, S(2)))))
    def With55(x, u):
        lst = PowerVariableExpn(u, S(0), x)
        return Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2)))/Part(lst, S(2))
    rule55 = ReplacementRule(pattern55, lambda x, u : With55(x, u))
    rubi.add(rule55)

    pattern56 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda m: Unequal(m, S(-1))), CustomConstraint(lambda u: NonsumQ(u)), CustomConstraint(lambda m, x, u: Greater(m, S(0)) | Not(AlgebraicFunctionQ(u, x))), CustomConstraint(lambda lst, m, x: Not(FalseQ(lst)) & NonzeroQ(-m + Part(lst, S(2)) + S(-1))))
    def With56(m, x, u):
        lst = PowerVariableExpn(u, m + S(1), x)
        return Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2)))/Part(lst, S(2))
    rule56 = ReplacementRule(pattern56, lambda m, x, u : With56(m, x, u))
    rubi.add(rule56)

    pattern57 = Pattern(Integral(u_*x_**m_, x_), CustomConstraint(lambda m: FractionQ(m)), )
    def With57(m, x, u):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*ReplaceAll(u, Rule(x, x**k)), x), x, x**(1/k))
    rule57 = ReplacementRule(pattern57, lambda m, x, u : With57(m, x, u))
    rubi.add(rule57)

    pattern58 = Pattern(Integral(u_, x_), CustomConstraint(lambda x, u: EulerIntegrandQ(u, x)), CustomConstraint(lambda lst, x: Not(FalseQ(lst))))
    def With58(x, u):
        lst = FunctionOfSquareRootOfQuadratic(u, x)
        return S(2)*Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(2)))
    rule58 = ReplacementRule(pattern58, lambda x, u : With58(x, u))
    #rubi.add(rule58)

    pattern59 = Pattern(Integral(1/(a_ + v_**S(2)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)))
    rule59 = ReplacementRule(pattern59, lambda b, x, v, a : Int(Rt(-a/b, S(2))/(-v + Rt(-a/b, S(2))), x)/(S(2)*a) + Int(Rt(-a/b, S(2))/(v + Rt(-a/b, S(2))), x)/(S(2)*a))
    rubi.add(rule59)

    pattern60 = Pattern(Integral(1/(a_ + v_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: EvenQ(n)), CustomConstraint(lambda n: Greater(n, S(2))))
    rule60 = ReplacementRule(pattern60, lambda v, a, x, n, b : Dist(S(2)/(a*n), Sum(Int((S(-1))**(S(4)*k/n)*Rt(-a/b, n/S(2))/((S(-1))**(S(4)*k/n)*Rt(-a/b, n/S(2)) - v**S(2)), x), List(k, S(1), n/S(2))), x))
    rubi.add(rule60)

    pattern61 = Pattern(Integral(1/(a_ + v_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: OddQ(n)), CustomConstraint(lambda n: Greater(n, S(1))))
    rule61 = ReplacementRule(pattern61, lambda v, a, x, n, b : Dist(S(1)/(a*n), Sum(Int((S(-1))**(S(2)*k/n)*Rt(-a/b, n)/((S(-1))**(S(2)*k/n)*Rt(-a/b, n) - v), x), List(k, S(1), n)), x))
    rubi.add(rule61)

    pattern62 = Pattern(Integral(v_/(a_ + u_**WC('n', S(1))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda u, x, v: PolynomialInQ(v, u, x)))
    rule62 = ReplacementRule(pattern62, lambda v, u, a, x, n, b : Int(ReplaceAll(ExpandIntegrand(PolynomialInSubst(v, u, x)/(a + b*x**n), x), Rule(x, u)), x))
    rubi.add(rule62)

    pattern63 = Pattern(Integral(u_, x_), CustomConstraint(lambda u, v, x: UnsameQ(v, u)))
    def With63(x, u):
        v = NormalizeIntegrand(u, x)
        return Int(v, x)
    rule63 = ReplacementRule(pattern63, lambda x, u : With63(x, u))
    rubi.add(rule63)

    pattern64 = Pattern(Integral(u_, x_), CustomConstraint(lambda v, x: SumQ(v)))
    def With64(x, u):
        v = ExpandIntegrand(u, x)
        return Int(v, x)
    rule64 = ReplacementRule(pattern64, lambda x, u : With64(x, u))
    rubi.add(rule64)

    pattern65 = Pattern(Integral((x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda d, a: ZeroQ(a + d)), CustomConstraint(lambda b, c: ZeroQ(b + c)), CustomConstraint(lambda m, n: ZeroQ(m + n)), CustomConstraint(lambda q, p: ZeroQ(p + q)))
    rule65 = ReplacementRule(pattern65, lambda u, a, m, d, c, q, x, p, n, b : x**(-m*p)*(a + b*x**m)**p*(c + d*x**n)**q*Int(u*x**(m*p), x))
    rubi.add(rule65)

    pattern66 = Pattern(Integral(u_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, n2: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda b, c, a: ZeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: IntegerQ(p + S(-1)/2)))
    rule66 = ReplacementRule(pattern66, lambda n2, u, a, c, x, p, n, b : (S(4)*c)**(-p + S(1)/2)*sqrt(a + b*x**n + c*x**(S(2)*n))*Int(u*(b + S(2)*c*x**n)**(S(2)*p), x)/(b + S(2)*c*x**n))
    rubi.add(rule66)

    pattern67 = Pattern(Integral(u_, x_), CustomConstraint(lambda lst, x: Not(FalseQ(lst))))
    def With67(x, u):
        lst = SubstForFractionalPowerOfLinear(u, x)
        return Part(lst, S(2))*Part(lst, S(4))*Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(1/Part(lst, S(2))))
    rule67 = ReplacementRule(pattern67, lambda x, u : With67(x, u))
    #rubi.add(rule67)

    pattern68 = Pattern(Integral(u_, x_))
    rule68 = ReplacementRule(pattern68, lambda x, u : u)
    rubi.add(rule68)

    return rubi
