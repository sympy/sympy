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

def miscellaneous_algebraic(rubi):

    pattern1 = Pattern(Integral(((x_**n_*WC('c', S(1)))**q_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, q: IntegerQ(n*q)))
    rule1 = ReplacementRule(pattern1, lambda p, n, a, q, x, b, c : x*(c*x**n)**(-S(1)/n)*Subst(Int((a + b*x**(n*q))**p, x), x, (c*x**n)**(1/n)))
    rubi.add(rule1)

    pattern2 = Pattern(Integral(x_**WC('m', S(1))*((x_**n_*WC('c', S(1)))**q_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda n, q: IntegerQ(n*q)), CustomConstraint(lambda m: IntegerQ(m)))
    rule2 = ReplacementRule(pattern2, lambda p, n, a, q, x, m, b, c : x**(m + S(1))*(c*x**n)**(-(m + S(1))/n)*Subst(Int(x**m*(a + b*x**(n*q))**p, x), x, (c*x**n)**(1/n)))
    rubi.add(rule2)

    pattern3 = Pattern(Integral(x_**WC('m', S(1))*((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('r', S(1))*WC('e', S(1)))**p_*((c_ + x_**WC('n', S(1))*WC('d', S(1)))**s_*WC('f', S(1)))**q_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda s, x: FreeQ(s, x)))
    rule3 = ReplacementRule(pattern3, lambda r, p, n, a, s, q, f, x, m, b, d, e, c : (e*(a + b*x**n)**r)**p*(f*(c + d*x**n)**s)**q*(a + b*x**n)**(-p*r)*(c + d*x**n)**(-q*s)*Int(x**m*(a + b*x**n)**(p*r)*(c + d*x**n)**(q*s), x))
    rubi.add(rule3)

    pattern4 = Pattern(Integral(((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, a, d, c: ZeroQ(-a*d + b*c)))
    rule4 = ReplacementRule(pattern4, lambda p, n, a, x, u, b, d, e, c : (b*e/d)**p*Int(u, x))
    rubi.add(rule4)

    pattern5 = Pattern(Integral(((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_*WC('u', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, e, d: PositiveQ(b*d*e)), CustomConstraint(lambda b, a, d, c: PositiveQ(-a*d/b + c)))
    rule5 = ReplacementRule(pattern5, lambda p, n, a, x, u, b, d, e, c : Int(u*(e*(a + b*x**n))**p*(c + d*x**n)**(-p), x))
    rubi.add(rule5)

    pattern6 = Pattern(Integral(((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p: FractionQ(p)), CustomConstraint(lambda n: IntegerQ(1/n)), )
    def With6(e, p, n, a, x, b, d, c):
        q = Denominator(p)
        return e*q*(-a*d + b*c)*Subst(Int(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + 1/n)*(b*e - d*x**q)**(S(-1) - S(1)/n), x), x, (e*(a + b*x**n)/(c + d*x**n))**(1/q))/n
    rule6 = ReplacementRule(pattern6, lambda e, p, n, a, x, b, d, c : With6(e, p, n, a, x, b, d, c))
    rubi.add(rule6)

    pattern7 = Pattern(Integral(x_**WC('m', S(1))*((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: FractionQ(p)), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)), )
    def With7(p, n, a, x, m, b, d, e, c):
        q = Denominator(p)
        return e*q*(-a*d + b*c)*Subst(Int(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + (m + S(1))/n)*(b*e - d*x**q)**(S(-1) - (m + S(1))/n), x), x, (e*(a + b*x**n)/(c + d*x**n))**(1/q))/n
    rule7 = ReplacementRule(pattern7, lambda p, n, a, x, m, b, d, e, c : With7(p, n, a, x, m, b, d, e, c))
    rubi.add(rule7)

    pattern8 = Pattern(Integral(u_**WC('r', S(1))*((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda x, u: PolynomialQ(u, x)), CustomConstraint(lambda p: FractionQ(p)), CustomConstraint(lambda n: IntegerQ(1/n)), CustomConstraint(lambda r: IntegerQ(r)), )
    def With8(r, p, n, a, x, u, b, d, e, c):
        q = Denominator(p)
        return e*q*(-a*d + b*c)*Subst(Int(SimplifyIntegrand(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + 1/n)*(b*e - d*x**q)**(S(-1) - S(1)/n)*ReplaceAll(u, Rule(x, (-a*e + c*x**q)**(1/n)*(b*e - d*x**q)**(-S(1)/n)))**r, x), x), x, (e*(a + b*x**n)/(c + d*x**n))**(1/q))/n
    rule8 = ReplacementRule(pattern8, lambda r, p, n, a, x, u, b, d, e, c : With8(r, p, n, a, x, u, b, d, e, c))
    rubi.add(rule8)

    pattern9 = Pattern(Integral(u_**WC('r', S(1))*x_**WC('m', S(1))*((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda x, u: PolynomialQ(u, x)), CustomConstraint(lambda p: FractionQ(p)), CustomConstraint(lambda n: IntegerQ(1/n)), CustomConstraint(lambda r, m: IntegersQ(m, r)), )
    def With9(r, p, n, a, x, u, m, b, d, e, c):
        q = Denominator(p)
        return e*q*(-a*d + b*c)*Subst(Int(SimplifyIntegrand(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + (m + S(1))/n)*(b*e - d*x**q)**(S(-1) - (m + S(1))/n)*ReplaceAll(u, Rule(x, (-a*e + c*x**q)**(1/n)*(b*e - d*x**q)**(-S(1)/n)))**r, x), x), x, (e*(a + b*x**n)/(c + d*x**n))**(1/q))/n
    rule9 = ReplacementRule(pattern9, lambda r, p, n, a, x, u, m, b, d, e, c : With9(r, p, n, a, x, u, m, b, d, e, c))
    rubi.add(rule9)

    pattern10 = Pattern(Integral(((WC('c', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)))
    rule10 = ReplacementRule(pattern10, lambda p, n, a, x, b, c : -c*Subst(Int((a + b*x**n)**p/x**S(2), x), x, c/x))
    rubi.add(rule10)

    pattern11 = Pattern(Integral(x_**WC('m', S(1))*((WC('c', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: IntegerQ(m)))
    rule11 = ReplacementRule(pattern11, lambda p, n, a, x, m, b, c : -c**(m + S(1))*Subst(Int(x**(-m + S(-2))*(a + b*x**n)**p, x), x, c/x))
    rubi.add(rule11)

    pattern12 = Pattern(Integral((x_*WC('d', S(1)))**m_*((WC('c', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule12 = ReplacementRule(pattern12, lambda p, n, a, x, m, b, d, c : -c*(c/x)**m*(d*x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**n)**p, x), x, c/x))
    rubi.add(rule12)

    pattern13 = Pattern(Integral(((WC('d', S(1))/x_)**n_*WC('b', S(1)) + (WC('d', S(1))/x_)**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)))
    rule13 = ReplacementRule(pattern13, lambda n2, p, n, a, x, b, d, c : -d*Subst(Int((a + b*x**n + c*x**(S(2)*n))**p/x**S(2), x), x, d/x))
    rubi.add(rule13)

    pattern14 = Pattern(Integral(x_**WC('m', S(1))*(a_ + (WC('d', S(1))/x_)**n_*WC('b', S(1)) + (WC('d', S(1))/x_)**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda m: IntegerQ(m)))
    rule14 = ReplacementRule(pattern14, lambda n2, p, n, a, x, m, b, d, c : -d**(m + S(1))*Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*x**(S(2)*n))**p, x), x, d/x))
    rubi.add(rule14)

    pattern15 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + (WC('d', S(1))/x_)**n_*WC('b', S(1)) + (WC('d', S(1))/x_)**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda m: Not(IntegerQ(m))))
    rule15 = ReplacementRule(pattern15, lambda n2, e, p, n, a, x, m, b, d, c : -d*(d/x)**m*(e*x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*x**(S(2)*n))**p, x), x, d/x))
    rubi.add(rule15)

    pattern16 = Pattern(Integral((x_**WC('n2', S(1))*WC('c', S(1)) + (WC('d', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(S(2)*n + n2)), CustomConstraint(lambda n: IntegerQ(S(2)*n)))
    rule16 = ReplacementRule(pattern16, lambda n2, p, n, a, x, b, d, c : -d*Subst(Int((a + b*x**n + c*d**(-S(2)*n)*x**(S(2)*n))**p/x**S(2), x), x, d/x))
    rubi.add(rule16)

    pattern17 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)) + (WC('d', S(1))/x_)**n_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(S(2)*n + n2)), CustomConstraint(lambda n: IntegerQ(S(2)*n)), CustomConstraint(lambda m: IntegerQ(m)))
    rule17 = ReplacementRule(pattern17, lambda n2, p, n, a, x, m, b, d, c : -d**(m + S(1))*Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*d**(-S(2)*n)*x**(S(2)*n))**p, x), x, d/x))
    rubi.add(rule17)

    pattern18 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)) + (WC('d', S(1))/x_)**n_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(S(2)*n + n2)), CustomConstraint(lambda m: Not(IntegerQ(m))), CustomConstraint(lambda n: IntegerQ(S(2)*n)))
    rule18 = ReplacementRule(pattern18, lambda n2, e, p, n, a, x, m, b, d, c : -d*(d/x)**m*(e*x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*d**(-S(2)*n)*x**(S(2)*n))**p, x), x, d/x))
    rubi.add(rule18)

    pattern19 = Pattern(Integral(u_**m_, x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda x, u: Not(LinearMatchQ(u, x))))
    rule19 = ReplacementRule(pattern19, lambda x, u, m : Int(ExpandToSum(u, x)**m, x))
    rubi.add(rule19)

    pattern20 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda v, x, u: LinearQ(List(u, v), x)), CustomConstraint(lambda v, x, u: Not(LinearMatchQ(List(u, v), x))))
    rule20 = ReplacementRule(pattern20, lambda n, v, x, u, m : Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n, x))
    rubi.add(rule20)

    pattern21 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1))*w_**WC('p', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda w, v, x, u: LinearQ(List(u, v, w), x)), CustomConstraint(lambda w, v, x, u: Not(LinearMatchQ(List(u, v, w), x))))
    rule21 = ReplacementRule(pattern21, lambda p, n, w, v, x, u, m : Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n*ExpandToSum(w, x)**p, x))
    rubi.add(rule21)

    pattern22 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1))*w_**WC('p', S(1))*z_**WC('q', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda w, v, x, u, z: LinearQ(List(u, v, w, z), x)), CustomConstraint(lambda w, v, x, u, z: Not(LinearMatchQ(List(u, v, w, z), x))))
    rule22 = ReplacementRule(pattern22, lambda p, n, q, w, v, x, u, m, z : Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n*ExpandToSum(w, x)**p*ExpandToSum(z, x)**q, x))
    rubi.add(rule22)

    pattern23 = Pattern(Integral(u_**p_, x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: QuadraticQ(u, x)), CustomConstraint(lambda x, u: Not(QuadraticMatchQ(u, x))))
    rule23 = ReplacementRule(pattern23, lambda p, x, u : Int(ExpandToSum(u, x)**p, x))
    rubi.add(rule23)

    pattern24 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('p', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda v, x, u: Not(LinearMatchQ(u, x) & QuadraticMatchQ(v, x))))
    rule24 = ReplacementRule(pattern24, lambda p, v, x, u, m : Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**p, x))
    rubi.add(rule24)

    pattern25 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1))*w_**WC('p', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda v, x, u: LinearQ(List(u, v), x)), CustomConstraint(lambda w, x: QuadraticQ(w, x)), CustomConstraint(lambda w, v, x, u: Not(QuadraticMatchQ(w, x) & LinearMatchQ(List(u, v), x))))
    rule25 = ReplacementRule(pattern25, lambda p, n, w, v, x, u, m : Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n*ExpandToSum(w, x)**p, x))
    rubi.add(rule25)

    pattern26 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1)), x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda v, x, u: QuadraticQ(List(u, v), x)), CustomConstraint(lambda v, x, u: Not(QuadraticMatchQ(List(u, v), x))))
    rule26 = ReplacementRule(pattern26, lambda p, q, v, x, u : Int(ExpandToSum(u, x)**p*ExpandToSum(v, x)**q, x))
    rubi.add(rule26)

    pattern27 = Pattern(Integral(u_**p_, x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: BinomialQ(u, x)), CustomConstraint(lambda x, u: Not(BinomialMatchQ(u, x))))
    rule27 = ReplacementRule(pattern27, lambda p, x, u : Int(ExpandToSum(u, x)**p, x))
    rubi.add(rule27)

    pattern28 = Pattern(Integral(u_**WC('p', S(1))*(x_*WC('c', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: BinomialQ(u, x)), CustomConstraint(lambda x, u: Not(BinomialMatchQ(u, x))))
    rule28 = ReplacementRule(pattern28, lambda p, x, u, m, c : Int((c*x)**m*ExpandToSum(u, x)**p, x))
    rubi.add(rule28)

    pattern29 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1)), x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda v, x, u: BinomialQ(List(u, v), x)), CustomConstraint(lambda v, x, u: ZeroQ(BinomialDegree(u, x) - BinomialDegree(v, x))), CustomConstraint(lambda v, x, u: Not(BinomialMatchQ(List(u, v), x))))
    rule29 = ReplacementRule(pattern29, lambda p, q, v, x, u : Int(ExpandToSum(u, x)**p*ExpandToSum(v, x)**q, x))
    rubi.add(rule29)

    pattern30 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*x_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda v, x, u: BinomialQ(List(u, v), x)), CustomConstraint(lambda v, x, u: ZeroQ(BinomialDegree(u, x) - BinomialDegree(v, x))), CustomConstraint(lambda v, x, u: Not(BinomialMatchQ(List(u, v), x))))
    rule30 = ReplacementRule(pattern30, lambda p, q, v, x, u, m : Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(v, x)**q, x))
    rubi.add(rule30)

    pattern31 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('p', S(1))*w_**WC('q', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda w, v, x, u: BinomialQ(List(u, v, w), x)), CustomConstraint(lambda v, x, u: ZeroQ(BinomialDegree(u, x) - BinomialDegree(v, x))), CustomConstraint(lambda w, x, u: ZeroQ(BinomialDegree(u, x) - BinomialDegree(w, x))), CustomConstraint(lambda w, v, x, u: Not(BinomialMatchQ(List(u, v, w), x))))
    rule31 = ReplacementRule(pattern31, lambda p, q, w, v, x, u, m : Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**p*ExpandToSum(w, x)**q, x))
    rubi.add(rule31)

    pattern32 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*x_**WC('m', S(1))*z_**WC('r', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda r, x: FreeQ(r, x)), CustomConstraint(lambda z, v, x, u: BinomialQ(List(u, v, z), x)), CustomConstraint(lambda v, x, u: ZeroQ(BinomialDegree(u, x) - BinomialDegree(v, x))), CustomConstraint(lambda z, x, u: ZeroQ(BinomialDegree(u, x) - BinomialDegree(z, x))), CustomConstraint(lambda z, v, x, u: Not(BinomialMatchQ(List(u, v, z), x))))
    rule32 = ReplacementRule(pattern32, lambda r, p, q, v, x, u, m, z : Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(v, x)**q*ExpandToSum(z, x)**r, x))
    rubi.add(rule32)

    pattern33 = Pattern(Integral(u_**p_, x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: GeneralizedBinomialQ(u, x)), CustomConstraint(lambda x, u: Not(GeneralizedBinomialMatchQ(u, x))))
    rule33 = ReplacementRule(pattern33, lambda p, x, u : Int(ExpandToSum(u, x)**p, x))
    rubi.add(rule33)

    pattern34 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: GeneralizedBinomialQ(u, x)), CustomConstraint(lambda x, u: Not(GeneralizedBinomialMatchQ(u, x))))
    rule34 = ReplacementRule(pattern34, lambda p, x, u, m : Int(x**m*ExpandToSum(u, x)**p, x))
    rubi.add(rule34)

    pattern35 = Pattern(Integral(u_**p_, x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: TrinomialQ(u, x)), CustomConstraint(lambda x, u: Not(TrinomialMatchQ(u, x))))
    rule35 = ReplacementRule(pattern35, lambda p, x, u : Int(ExpandToSum(u, x)**p, x))
    rubi.add(rule35)

    pattern36 = Pattern(Integral(u_**WC('p', S(1))*(x_*WC('d', S(1)))**WC('m', S(1)), x_), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: TrinomialQ(u, x)), CustomConstraint(lambda x, u: Not(TrinomialMatchQ(u, x))))
    rule36 = ReplacementRule(pattern36, lambda p, x, u, m, d : Int((d*x)**m*ExpandToSum(u, x)**p, x))
    rubi.add(rule36)

    pattern37 = Pattern(Integral(u_**WC('q', S(1))*v_**WC('p', S(1)), x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda x, u: BinomialQ(u, x)), CustomConstraint(lambda v, x: TrinomialQ(v, x)), CustomConstraint(lambda v, x, u: Not(BinomialMatchQ(u, x) & TrinomialMatchQ(v, x))))
    rule37 = ReplacementRule(pattern37, lambda p, q, v, x, u : Int(ExpandToSum(u, x)**q*ExpandToSum(v, x)**p, x))
    rubi.add(rule37)

    pattern38 = Pattern(Integral(u_**WC('q', S(1))*v_**WC('p', S(1)), x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda x, u: BinomialQ(u, x)), CustomConstraint(lambda v, x: BinomialQ(v, x)), CustomConstraint(lambda v, x, u: Not(BinomialMatchQ(u, x) & BinomialMatchQ(v, x))))
    rule38 = ReplacementRule(pattern38, lambda p, q, v, x, u : Int(ExpandToSum(u, x)**q*ExpandToSum(v, x)**p, x))
    rubi.add(rule38)

    pattern39 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1))*z_**WC('q', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda z, x: BinomialQ(z, x)), CustomConstraint(lambda x, u: TrinomialQ(u, x)), CustomConstraint(lambda z, x, u: Not(BinomialMatchQ(z, x) & TrinomialMatchQ(u, x))))
    rule39 = ReplacementRule(pattern39, lambda p, q, x, u, m, z : Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(z, x)**q, x))
    rubi.add(rule39)

    pattern40 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1))*z_**WC('q', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda z, x: BinomialQ(z, x)), CustomConstraint(lambda x, u: BinomialQ(u, x)), CustomConstraint(lambda z, x, u: Not(BinomialMatchQ(u, x) & BinomialMatchQ(z, x))))
    rule40 = ReplacementRule(pattern40, lambda p, q, x, u, m, z : Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(z, x)**q, x))
    rubi.add(rule40)

    pattern41 = Pattern(Integral(u_**p_, x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: GeneralizedTrinomialQ(u, x)), CustomConstraint(lambda x, u: Not(GeneralizedTrinomialMatchQ(u, x))))
    rule41 = ReplacementRule(pattern41, lambda p, x, u : Int(ExpandToSum(u, x)**p, x))
    rubi.add(rule41)

    pattern42 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: GeneralizedTrinomialQ(u, x)), CustomConstraint(lambda x, u: Not(GeneralizedTrinomialMatchQ(u, x))))
    rule42 = ReplacementRule(pattern42, lambda p, x, u, m : Int(x**m*ExpandToSum(u, x)**p, x))
    rubi.add(rule42)

    pattern43 = Pattern(Integral(u_**WC('p', S(1))*z_, x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda z, x: BinomialQ(z, x)), CustomConstraint(lambda x, u: GeneralizedTrinomialQ(u, x)), CustomConstraint(lambda z, x, u: ZeroQ(BinomialDegree(z, x) - GeneralizedTrinomialDegree(u, x))), CustomConstraint(lambda z, x, u: Not(BinomialMatchQ(z, x) & GeneralizedTrinomialMatchQ(u, x))))
    rule43 = ReplacementRule(pattern43, lambda z, p, x, u : Int(ExpandToSum(u, x)**p*ExpandToSum(z, x), x))
    rubi.add(rule43)

    pattern44 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1))*z_, x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda z, x: BinomialQ(z, x)), CustomConstraint(lambda x, u: GeneralizedTrinomialQ(u, x)), CustomConstraint(lambda z, x, u: ZeroQ(BinomialDegree(z, x) - GeneralizedTrinomialDegree(u, x))), CustomConstraint(lambda z, x, u: Not(BinomialMatchQ(z, x) & GeneralizedTrinomialMatchQ(u, x))))
    rule44 = ReplacementRule(pattern44, lambda p, x, u, m, z : Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(z, x), x))
    rubi.add(rule44)

    pattern45 = Pattern(Integral(x_**WC('m', S(1))*(e_ + x_**WC('n', S(1))*WC('h', S(1)) + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)))/(a_ + x_**WC('n', S(1))*WC('c', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, q: ZeroQ(-n/S(4) + q)), CustomConstraint(lambda n, r: ZeroQ(-S(3)*n/S(4) + r)), CustomConstraint(lambda n, m: ZeroQ(S(4)*m - n + S(4))), CustomConstraint(lambda h, a, e, c: ZeroQ(a*h + c*e)))
    rule45 = ReplacementRule(pattern45, lambda r, n, a, q, g, f, x, m, h, e, c : (-S(2)*a*g - S(4)*a*h*x**(n/S(4)) + S(2)*c*f*x**(n/S(2)))/(a*c*n*sqrt(a + c*x**n)))
    rubi.add(rule45)

    pattern46 = Pattern(Integral((d_*x_)**WC('m', S(1))*(e_ + x_**WC('n', S(1))*WC('h', S(1)) + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)))/(a_ + x_**WC('n', S(1))*WC('c', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, m: ZeroQ(S(4)*m - n + S(4))), CustomConstraint(lambda n, q: ZeroQ(-n/S(4) + q)), CustomConstraint(lambda n, r: ZeroQ(-S(3)*n/S(4) + r)), CustomConstraint(lambda h, a, e, c: ZeroQ(a*h + c*e)))
    rule46 = ReplacementRule(pattern46, lambda r, e, n, a, q, g, f, x, m, h, d, c : x**(-m)*(d*x)**m*Int(x**m*(e + f*x**(n/S(4)) + g*x**(S(3)*n/S(4)) + h*x**n)/(a + c*x**n)**(S(3)/2), x))
    rubi.add(rule46)

    pattern47 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**m_*(a_ + x_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda p: FractionQ(p)), CustomConstraint(lambda m: NegativeIntegerQ(m + S(1))), )
    def With47(p, Pq, a, x, m, b, c):
        n = Denominator(p)
        return n*Subst(Int(x**(n*p + n + S(-1))*(-a*c/b + c*x**n/b)**m*ReplaceAll(Pq, Rule(x, -a/b + x**n/b)), x), x, (a + b*x)**(1/n))/b
    rule47 = ReplacementRule(pattern47, lambda p, Pq, a, x, m, b, c : With47(p, Pq, a, x, m, b, c))
    rubi.add(rule47)

    pattern48 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: PositiveIntegerQ(n/(m + S(1)))), CustomConstraint(lambda Pq, x, m: PolyQ(Pq, x**(m + S(1)))))
    rule48 = ReplacementRule(pattern48, lambda p, n, Pq, a, x, m, b : Subst(Int((a + b*x**(n/(m + S(1))))**p*SubstFor(x**(m + S(1)), Pq, x), x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule48)

    pattern49 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(n, p)), CustomConstraint(lambda n, Pq, x: NonzeroQ(Coeff(Pq, x, n + S(-1)))))
    rule49 = ReplacementRule(pattern49, lambda p, n, Pq, a, x, b : Int((a + b*x**n)**p*ExpandToSum(Pq - x**(n + S(-1))*Coeff(Pq, x, n + S(-1)), x), x) + (a + b*x**n)**(p + S(1))*Coeff(Pq, x, n + S(-1))/(b*n*(p + S(1))))
    rubi.add(rule49)

    pattern50 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(p) | ZeroQ(n + S(-1))))
    rule50 = ReplacementRule(pattern50, lambda p, n, Pq, a, x, m, b, c : Int(ExpandIntegrand(Pq*(c*x)**m*(a + b*x**n)**p, x), x))
    rubi.add(rule50)

    pattern51 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n, p: PositiveIntegerQ(p) | ZeroQ(n + S(-1))))
    rule51 = ReplacementRule(pattern51, lambda p, n, Pq, a, x, b : Int(ExpandIntegrand(Pq*(a + b*x**n)**p, x), x))
    rubi.add(rule51)

    pattern52 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)))
    rule52 = ReplacementRule(pattern52, lambda p, n, Pq, a, x, m, b : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p*SubstFor(x**n, Pq, x), x), x, x**n)/n)
    rubi.add(rule52)

    pattern53 = Pattern(Integral(Pq_*(c_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)))
    rule53 = ReplacementRule(pattern53, lambda p, n, Pq, a, x, m, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(Pq*x**m*(a + b*x**n)**p, x))
    rubi.add(rule53)

    pattern54 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule54 = ReplacementRule(pattern54, lambda p, n, Pq, a, x, m, b : Pq*(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))) - Int((a + b*x**n)**(p + S(1))*D(Pq, x), x)/(b*n*(p + S(1))))
    rubi.add(rule54)

    pattern55 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda Pq, x: ZeroQ(Coeff(Pq, x, S(0)))))
    rule55 = ReplacementRule(pattern55, lambda p, n, Pq, a, x, m, b, d : Int((d*x)**(m + S(1))*(a + b*x**n)**p*ExpandToSum(Pq/x, x), x)/d)
    rubi.add(rule55)

    pattern56 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda Pq, x: ZeroQ(Coeff(Pq, x, S(0)))), CustomConstraint(lambda Pq: SumQ(Pq)))
    rule56 = ReplacementRule(pattern56, lambda p, n, Pq, a, x, b : Int(x*(a + b*x**n)**p*ExpandToSum(Pq/x, x), x))
    rubi.add(rule56)

    pattern57 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p, m: RationalQ(m, p)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda Pq, x, m: Less(m + Expon(Pq, x) + S(1), S(0))), )
    def With57(p, n, Pq, a, x, m, b):
        u = IntHide(Pq*x**m, x)
        return -b*n*p*Int(x**(m + n)*(a + b*x**n)**(p + S(-1))*ExpandToSum(u*x**(-m + S(-1)), x), x) + u*(a + b*x**n)**p
    rule57 = ReplacementRule(pattern57, lambda p, n, Pq, a, x, m, b : With57(p, n, Pq, a, x, m, b))
    rubi.add(rule57)

    pattern58 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2) + S(-1)/2)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), )
    def With58(p, n, Pq, a, x, m, b, c):
        q = Expon(Pq, x)
        i = Symbol('i')
        return a*n*p*Int((c*x)**m*(a + b*x**n)**(p + S(-1))*Sum(x**i*Coeff(Pq, x, i)/(i + m + n*p + S(1)), List(i, S(0), q)), x) + (c*x)**m*(a + b*x**n)**p*Sum(x**(i + S(1))*Coeff(Pq, x, i)/(i + m + n*p + S(1)), List(i, S(0), q))
    rule58 = ReplacementRule(pattern58, lambda p, n, Pq, a, x, m, b, c : With58(p, n, Pq, a, x, m, b, c))
    rubi.add(rule58)

    pattern59 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2) + S(-1)/2)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Greater(p, S(0))), )
    def With59(p, n, Pq, a, x, b):
        q = Expon(Pq, x)
        i = Symbol('i')
        return a*n*p*Int((a + b*x**n)**(p + S(-1))*Sum(x**i*Coeff(Pq, x, i)/(i + n*p + S(1)), List(i, S(0), q)), x) + (a + b*x**n)**p*Sum(x**(i + S(1))*Coeff(Pq, x, i)/(i + n*p + S(1)), List(i, S(0), q))
    rule59 = ReplacementRule(pattern59, lambda p, n, Pq, a, x, b : With59(p, n, Pq, a, x, b))
    rubi.add(rule59)

    pattern60 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda b, n, Pq, a, i, x, p, q: Equal(q, n + S(-1))))
    def With60(p, n, Pq, a, x, b):
        q = Expon(Pq, x)
        i = Symbol('i')
        return Int((a + b*x**n)**(p + S(1))*Sum(x**i*(i + n*(p + S(1)) + S(1))*Coeff(Pq, x, i), List(i, S(0), q + S(-1))), x)/(a*n*(p + S(1))) + (a + b*x**n)**(p + S(1))*(a*Coeff(Pq, x, q) - b*x*ExpandToSum(Pq - x**q*Coeff(Pq, x, q), x))/(a*b*n*(p + S(1)))
    rule60 = ReplacementRule(pattern60, lambda p, n, Pq, a, x, b : With60(p, n, Pq, a, x, b))
    rubi.add(rule60)

    pattern61 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, Pq, x: Less(Expon(Pq, x), n + S(-1))))
    rule61 = ReplacementRule(pattern61, lambda p, n, Pq, a, x, b : -Pq*x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))) + Int((a + b*x**n)**(p + S(1))*ExpandToSum(Pq*n*(p + S(1)) + D(Pq*x, x), x), x)/(a*n*(p + S(1))))
    rubi.add(rule61)

    pattern62 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, a, d, g: ZeroQ(a*g + b*d)))
    rule62 = ReplacementRule(pattern62, lambda a, f, x, b, d, e, g : (-S(2)*a*f - S(4)*a*g*x + S(2)*b*e*x**S(2))/(S(4)*a*b*sqrt(a + b*x**S(4))))
    rubi.add(rule62)

    pattern63 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, a, d, g: ZeroQ(a*g + b*d)))
    rule63 = ReplacementRule(pattern63, lambda a, f, x, b, d, g : (-f - S(2)*g*x)/(S(2)*b*sqrt(a + b*x**S(4))))
    rubi.add(rule63)

    pattern64 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, a, d, g: ZeroQ(a*g + b*d)))
    rule64 = ReplacementRule(pattern64, lambda a, x, b, d, e, g : -x*(S(2)*a*g - b*e*x)/(S(2)*a*b*sqrt(a + b*x**S(4))))
    rubi.add(rule64)

    pattern65 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda b, h, a, e: ZeroQ(-S(3)*a*h + b*e)))
    rule65 = ReplacementRule(pattern65, lambda a, f, x, b, h, e : (-f + S(2)*h*x**S(3))/(S(2)*b*sqrt(a + b*x**S(4))))
    rubi.add(rule65)

    pattern66 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda b, h, a, e: ZeroQ(-S(3)*a*h + b*e)))
    rule66 = ReplacementRule(pattern66, lambda a, x, b, h, e : h*x**S(3)/(b*sqrt(a + b*x**S(4))))
    rubi.add(rule66)

    pattern67 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda b, h, a, e: ZeroQ(-S(3)*a*h + b*e)), CustomConstraint(lambda b, a, d, g: ZeroQ(a*g + b*d)))
    rule67 = ReplacementRule(pattern67, lambda a, f, x, b, d, h, e, g : (-a*f + S(2)*a*h*x**S(3) + S(2)*b*d*x)/(S(2)*a*b*sqrt(a + b*x**S(4))))
    rubi.add(rule67)

    pattern68 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(4)*WC('g', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda b, h, a, e: ZeroQ(-S(3)*a*h + b*e)), CustomConstraint(lambda b, a, d, g: ZeroQ(a*g + b*d)))
    rule68 = ReplacementRule(pattern68, lambda a, x, b, d, h, e, g : x*(a*h*x**S(2) + b*d)/(a*b*sqrt(a + b*x**S(4))))
    rubi.add(rule68)

    pattern69 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda b, n, Q, Pq, a, i, x, p, q, R: GreaterEqual(q, n)))
    def With69(p, n, Pq, a, x, b):
        q = Expon(Pq, x)
        return Module(List(Set(Q, PolynomialQuotient(Pq*b**(Floor((q + S(-1))/n) + S(1)), a + b*x**n, x)), Set(R, PolynomialRemainder(Pq*b**(Floor((q + S(-1))/n) + S(1)), a + b*x**n, x)), i), -R*b**(-Floor((q + S(-1))/n) + S(-1))*x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))) + b**(-Floor((q + S(-1))/n) + S(-1))*Int((a + b*x**n)**(p + S(1))*ExpandToSum(Q*a*n*(p + S(1)) + R*n*(p + S(1)) + D(R*x, x), x), x)/(a*n*(p + S(1))))
    rule69 = ReplacementRule(pattern69, lambda p, n, Pq, a, x, b : With69(p, n, Pq, a, x, b))
    rubi.add(rule69)

    pattern70 = Pattern(Integral(Pq_*x_**m_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: NegativeIntegerQ(m)), )
    def With70(p, n, Pq, a, x, m, b):
        q = Expon(Pq, x)
        return Module(List(Set(Q, PolynomialQuotient(Pq*a*b**(Floor((q + S(-1))/n) + S(1))*x**m, a + b*x**n, x)), Set(R, PolynomialRemainder(Pq*a*b**(Floor((q + S(-1))/n) + S(1))*x**m, a + b*x**n, x)), i), -R*b**(-Floor((q + S(-1))/n) + S(-1))*x*(a + b*x**n)**(p + S(1))/(a**S(2)*n*(p + S(1))) + b**(-Floor((q + S(-1))/n) + S(-1))*Int(x**m*(a + b*x**n)**(p + S(1))*ExpandToSum(Q*n*x**(-m)*(p + S(1)) + Sum(x**(i - m)*(i + n*(p + S(1)) + S(1))*Coeff(R, x, i)/a, List(i, S(0), n + S(-1))), x), x)/(a*n*(p + S(1))))
    rule70 = ReplacementRule(pattern70, lambda p, n, Pq, a, x, m, b : With70(p, n, Pq, a, x, m, b))
    rubi.add(rule70)

    pattern71 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda b, n, Pq, m, a, x, p, g: Unequal(g, S(1))))
    def With71(p, n, Pq, a, x, m, b):
        g = GCD(m + S(1), n)
        return Subst(Int(x**(S(-1) + (m + S(1))/g)*(a + b*x**(n/g))**p*ReplaceAll(Pq, Rule(x, x**(1/g))), x), x, x**g)/g
    rule71 = ReplacementRule(pattern71, lambda p, n, Pq, a, x, m, b : With71(p, n, Pq, a, x, m, b))
    rubi.add(rule71)

    pattern72 = Pattern(Integral((A_ + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda B, b, a, A: ZeroQ(-A**S(3)*b + B**S(3)*a)))
    rule72 = ReplacementRule(pattern72, lambda a, x, A, B, b : B**S(3)*Int(1/(A**S(2) - A*B*x + B**S(2)*x**S(2)), x)/b)
    rubi.add(rule72)

    pattern73 = Pattern(Integral((A_ + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda B, b, a, A: NonzeroQ(-A**S(3)*b + B**S(3)*a)), CustomConstraint(lambda b, a: PosQ(a/b)), )
    def With73(a, x, A, B, b):
        r = Numerator(Rt(a/b, S(3)))
        s = Denominator(Rt(a/b, S(3)))
        return -r*(-A*s + B*r)*Int(1/(r + s*x), x)/(S(3)*a*s) + r*Int((r*(S(2)*A*s + B*r) + s*x*(-A*s + B*r))/(r**S(2) - r*s*x + s**S(2)*x**S(2)), x)/(S(3)*a*s)
    rule73 = ReplacementRule(pattern73, lambda a, x, A, B, b : With73(a, x, A, B, b))
    rubi.add(rule73)

    pattern74 = Pattern(Integral((A_ + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda B, b, a, A: NonzeroQ(-A**S(3)*b + B**S(3)*a)), CustomConstraint(lambda b, a: NegQ(a/b)), )
    def With74(a, x, A, B, b):
        r = Numerator(Rt(-a/b, S(3)))
        s = Denominator(Rt(-a/b, S(3)))
        return r*(A*s + B*r)*Int(1/(r - s*x), x)/(S(3)*a*s) - r*Int((r*(-S(2)*A*s + B*r) - s*x*(A*s + B*r))/(r**S(2) + r*s*x + s**S(2)*x**S(2)), x)/(S(3)*a*s)
    rule74 = ReplacementRule(pattern74, lambda a, x, A, B, b : With74(a, x, A, B, b))
    rubi.add(rule74)

    pattern75 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, C, A: ZeroQ(-A*C + B**S(2))), CustomConstraint(lambda B, b, a, C: ZeroQ(B**S(3)*b + C**S(3)*a)))
    rule75 = ReplacementRule(pattern75, lambda a, C, x, A, B, b : -C**S(2)*Int(1/(B - C*x), x)/b)
    rubi.add(rule75)

    pattern76 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda a, C, A, B, b: ZeroQ(A*b**(S(2)/3) - B*a**(S(1)/3)*b**(S(1)/3) - S(2)*C*a**(S(2)/3))), )
    def With76(a, C, x, A, B, b):
        q = a**(S(1)/3)/b**(S(1)/3)
        return C*Int(1/(q + x), x)/b + (B + C*q)*Int(1/(q**S(2) - q*x + x**S(2)), x)/b
    rule76 = ReplacementRule(pattern76, lambda a, C, x, A, B, b : With76(a, C, x, A, B, b))
    rubi.add(rule76)

    pattern77 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, C: ZeroQ(B*a**(S(1)/3)*b**(S(1)/3) + S(2)*C*a**(S(2)/3))), )
    def With77(a, C, x, B, b):
        q = a**(S(1)/3)/b**(S(1)/3)
        return C*Int(1/(q + x), x)/b + (B + C*q)*Int(1/(q**S(2) - q*x + x**S(2)), x)/b
    rule77 = ReplacementRule(pattern77, lambda a, C, x, B, b : With77(a, C, x, B, b))
    rubi.add(rule77)

    pattern78 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a, C, A: ZeroQ(A*b**(S(2)/3) - S(2)*C*a**(S(2)/3))), )
    def With78(a, C, x, A, b):
        q = a**(S(1)/3)/b**(S(1)/3)
        return C*q*Int(1/(q**S(2) - q*x + x**S(2)), x)/b + C*Int(1/(q + x), x)/b
    rule78 = ReplacementRule(pattern78, lambda a, C, x, A, b : With78(a, C, x, A, b))
    rubi.add(rule78)

    pattern79 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda a, C, A, B, b: ZeroQ(A*(-b)**(S(2)/3) - B*(-a)**(S(1)/3)*(-b)**(S(1)/3) - S(2)*C*(-a)**(S(2)/3))), )
    def With79(a, C, x, A, B, b):
        q = (-a)**(S(1)/3)/(-b)**(S(1)/3)
        return C*Int(1/(q + x), x)/b + (B + C*q)*Int(1/(q**S(2) - q*x + x**S(2)), x)/b
    rule79 = ReplacementRule(pattern79, lambda a, C, x, A, B, b : With79(a, C, x, A, B, b))
    rubi.add(rule79)

    pattern80 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, C: ZeroQ(B*(-a)**(S(1)/3)*(-b)**(S(1)/3) + S(2)*C*(-a)**(S(2)/3))), )
    def With80(a, C, x, B, b):
        q = (-a)**(S(1)/3)/(-b)**(S(1)/3)
        return C*Int(1/(q + x), x)/b + (B + C*q)*Int(1/(q**S(2) - q*x + x**S(2)), x)/b
    rule80 = ReplacementRule(pattern80, lambda a, C, x, B, b : With80(a, C, x, B, b))
    rubi.add(rule80)

    pattern81 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a, C, A: ZeroQ(A*(-b)**(S(2)/3) - S(2)*C*(-a)**(S(2)/3))), )
    def With81(a, C, x, A, b):
        q = (-a)**(S(1)/3)/(-b)**(S(1)/3)
        return C*q*Int(1/(q**S(2) - q*x + x**S(2)), x)/b + C*Int(1/(q + x), x)/b
    rule81 = ReplacementRule(pattern81, lambda a, C, x, A, b : With81(a, C, x, A, b))
    rubi.add(rule81)

    pattern82 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda a, C, A, B, b: ZeroQ(A*b**(S(2)/3) + B*b**(S(1)/3)*(-a)**(S(1)/3) - S(2)*C*(-a)**(S(2)/3))), )
    def With82(a, C, x, A, B, b):
        q = (-a)**(S(1)/3)/b**(S(1)/3)
        return -C*Int(1/(q - x), x)/b + (B - C*q)*Int(1/(q**S(2) + q*x + x**S(2)), x)/b
    rule82 = ReplacementRule(pattern82, lambda a, C, x, A, B, b : With82(a, C, x, A, B, b))
    rubi.add(rule82)

    pattern83 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, C: ZeroQ(B*b**(S(1)/3)*(-a)**(S(1)/3) - S(2)*C*(-a)**(S(2)/3))), )
    def With83(a, C, x, B, b):
        q = (-a)**(S(1)/3)/b**(S(1)/3)
        return -C*Int(1/(q - x), x)/b + (B - C*q)*Int(1/(q**S(2) + q*x + x**S(2)), x)/b
    rule83 = ReplacementRule(pattern83, lambda a, C, x, B, b : With83(a, C, x, B, b))
    rubi.add(rule83)

    pattern84 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a, C, A: ZeroQ(A*b**(S(2)/3) - S(2)*C*(-a)**(S(2)/3))), )
    def With84(a, C, x, A, b):
        q = (-a)**(S(1)/3)/b**(S(1)/3)
        return -C*q*Int(1/(q**S(2) + q*x + x**S(2)), x)/b - C*Int(1/(q - x), x)/b
    rule84 = ReplacementRule(pattern84, lambda a, C, x, A, b : With84(a, C, x, A, b))
    rubi.add(rule84)

    pattern85 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda a, C, A, B, b: ZeroQ(A*(-b)**(S(2)/3) + B*a**(S(1)/3)*(-b)**(S(1)/3) - S(2)*C*a**(S(2)/3))), )
    def With85(a, C, x, A, B, b):
        q = a**(S(1)/3)/(-b)**(S(1)/3)
        return -C*Int(1/(q - x), x)/b + (B - C*q)*Int(1/(q**S(2) + q*x + x**S(2)), x)/b
    rule85 = ReplacementRule(pattern85, lambda a, C, x, A, B, b : With85(a, C, x, A, B, b))
    rubi.add(rule85)

    pattern86 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, C: ZeroQ(B*a**(S(1)/3)*(-b)**(S(1)/3) - S(2)*C*a**(S(2)/3))), )
    def With86(a, C, x, B, b):
        q = a**(S(1)/3)/(-b)**(S(1)/3)
        return -C*Int(1/(q - x), x)/b + (B - C*q)*Int(1/(q**S(2) + q*x + x**S(2)), x)/b
    rule86 = ReplacementRule(pattern86, lambda a, C, x, B, b : With86(a, C, x, B, b))
    rubi.add(rule86)

    pattern87 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a, C, A: ZeroQ(A*(-b)**(S(2)/3) - S(2)*C*a**(S(2)/3))), )
    def With87(a, C, x, A, b):
        q = a**(S(1)/3)/(-b)**(S(1)/3)
        return -C*q*Int(1/(q**S(2) + q*x + x**S(2)), x)/b - C*Int(1/(q - x), x)/b
    rule87 = ReplacementRule(pattern87, lambda a, C, x, A, b : With87(a, C, x, A, b))
    rubi.add(rule87)

    pattern88 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda a, C, A, B, b: ZeroQ(A - B*(a/b)**(S(1)/3) - S(2)*C*(a/b)**(S(2)/3))), )
    def With88(a, C, x, A, B, b):
        q = (a/b)**(S(1)/3)
        return C*Int(1/(q + x), x)/b + (B + C*q)*Int(1/(q**S(2) - q*x + x**S(2)), x)/b
    rule88 = ReplacementRule(pattern88, lambda a, C, x, A, B, b : With88(a, C, x, A, B, b))
    rubi.add(rule88)

    pattern89 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, C: ZeroQ(B*(a/b)**(S(1)/3) + S(2)*C*(a/b)**(S(2)/3))), )
    def With89(a, C, x, B, b):
        q = (a/b)**(S(1)/3)
        return C*Int(1/(q + x), x)/b + (B + C*q)*Int(1/(q**S(2) - q*x + x**S(2)), x)/b
    rule89 = ReplacementRule(pattern89, lambda a, C, x, B, b : With89(a, C, x, B, b))
    rubi.add(rule89)

    pattern90 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a, C, A: ZeroQ(A - S(2)*C*(a/b)**(S(2)/3))), )
    def With90(a, C, x, A, b):
        q = (a/b)**(S(1)/3)
        return C*q*Int(1/(q**S(2) - q*x + x**S(2)), x)/b + C*Int(1/(q + x), x)/b
    rule90 = ReplacementRule(pattern90, lambda a, C, x, A, b : With90(a, C, x, A, b))
    rubi.add(rule90)

    pattern91 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda a, C, A, B, b: ZeroQ(A - B*Rt(a/b, S(3)) - S(2)*C*Rt(a/b, S(3))**S(2))), )
    def With91(a, C, x, A, B, b):
        q = Rt(a/b, S(3))
        return C*Int(1/(q + x), x)/b + (B + C*q)*Int(1/(q**S(2) - q*x + x**S(2)), x)/b
    rule91 = ReplacementRule(pattern91, lambda a, C, x, A, B, b : With91(a, C, x, A, B, b))
    rubi.add(rule91)

    pattern92 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, C: ZeroQ(B*Rt(a/b, S(3)) + S(2)*C*Rt(a/b, S(3))**S(2))), )
    def With92(a, C, x, B, b):
        q = Rt(a/b, S(3))
        return C*Int(1/(q + x), x)/b + (B + C*q)*Int(1/(q**S(2) - q*x + x**S(2)), x)/b
    rule92 = ReplacementRule(pattern92, lambda a, C, x, B, b : With92(a, C, x, B, b))
    rubi.add(rule92)

    pattern93 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a, C, A: ZeroQ(A - S(2)*C*Rt(a/b, S(3))**S(2))), )
    def With93(a, C, x, A, b):
        q = Rt(a/b, S(3))
        return C*q*Int(1/(q**S(2) - q*x + x**S(2)), x)/b + C*Int(1/(q + x), x)/b
    rule93 = ReplacementRule(pattern93, lambda a, C, x, A, b : With93(a, C, x, A, b))
    rubi.add(rule93)

    pattern94 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda a, C, A, B, b: ZeroQ(A + B*(-a/b)**(S(1)/3) - S(2)*C*(-a/b)**(S(2)/3))), )
    def With94(a, C, x, A, B, b):
        q = (-a/b)**(S(1)/3)
        return -C*Int(1/(q - x), x)/b + (B - C*q)*Int(1/(q**S(2) + q*x + x**S(2)), x)/b
    rule94 = ReplacementRule(pattern94, lambda a, C, x, A, B, b : With94(a, C, x, A, B, b))
    rubi.add(rule94)

    pattern95 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, C: ZeroQ(B*(-a/b)**(S(1)/3) - S(2)*C*(-a/b)**(S(2)/3))), )
    def With95(a, C, x, B, b):
        q = (-a/b)**(S(1)/3)
        return -C*Int(1/(q - x), x)/b + (B - C*q)*Int(1/(q**S(2) + q*x + x**S(2)), x)/b
    rule95 = ReplacementRule(pattern95, lambda a, C, x, B, b : With95(a, C, x, B, b))
    rubi.add(rule95)

    pattern96 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a, C, A: ZeroQ(A - S(2)*C*(-a/b)**(S(2)/3))), )
    def With96(a, C, x, A, b):
        q = (-a/b)**(S(1)/3)
        return -C*q*Int(1/(q**S(2) + q*x + x**S(2)), x)/b - C*Int(1/(q - x), x)/b
    rule96 = ReplacementRule(pattern96, lambda a, C, x, A, b : With96(a, C, x, A, b))
    rubi.add(rule96)

    pattern97 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda a, C, A, B, b: ZeroQ(A + B*Rt(-a/b, S(3)) - S(2)*C*Rt(-a/b, S(3))**S(2))), )
    def With97(a, C, x, A, B, b):
        q = Rt(-a/b, S(3))
        return -C*Int(1/(q - x), x)/b + (B - C*q)*Int(1/(q**S(2) + q*x + x**S(2)), x)/b
    rule97 = ReplacementRule(pattern97, lambda a, C, x, A, B, b : With97(a, C, x, A, B, b))
    rubi.add(rule97)

    pattern98 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, C: ZeroQ(B*Rt(-a/b, S(3)) - S(2)*C*Rt(-a/b, S(3))**S(2))), )
    def With98(a, C, x, B, b):
        q = Rt(-a/b, S(3))
        return -C*Int(1/(q - x), x)/b + (B - C*q)*Int(1/(q**S(2) + q*x + x**S(2)), x)/b
    rule98 = ReplacementRule(pattern98, lambda a, C, x, B, b : With98(a, C, x, B, b))
    rubi.add(rule98)

    pattern99 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a, C, A: ZeroQ(A - S(2)*C*Rt(-a/b, S(3))**S(2))), )
    def With99(a, C, x, A, b):
        q = Rt(-a/b, S(3))
        return -C*q*Int(1/(q**S(2) + q*x + x**S(2)), x)/b - C*Int(1/(q - x), x)/b
    rule99 = ReplacementRule(pattern99, lambda a, C, x, A, b : With99(a, C, x, A, b))
    rubi.add(rule99)

    pattern100 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, A: Not(RationalQ(a/b)) | ZeroQ(-A**S(3)*b + B**S(3)*a)))
    rule100 = ReplacementRule(pattern100, lambda a, C, x, A, B, b : C*Int(x**S(2)/(a + b*x**S(3)), x) + Int((A + B*x)/(a + b*x**S(3)), x))
    rubi.add(rule100)

    pattern101 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a: Not(RationalQ(a/b))))
    rule101 = ReplacementRule(pattern101, lambda a, C, x, B, b : B*Int(x/(a + b*x**S(3)), x) + C*Int(x**S(2)/(a + b*x**S(3)), x))
    rubi.add(rule101)

    pattern102 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a, C, A: Not(RationalQ(a, b, A, C))))
    rule102 = ReplacementRule(pattern102, lambda a, C, x, A, b : A*Int(1/(a + b*x**S(3)), x) + C*Int(x**S(2)/(a + b*x**S(3)), x))
    rubi.add(rule102)

    pattern103 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda a, C, A, B, b: ZeroQ(A - B*(a/b)**(S(1)/3) + C*(a/b)**(S(2)/3))), )
    def With103(a, C, x, A, B, b):
        q = (a/b)**(S(1)/3)
        return q**S(2)*Int((A + C*q*x)/(q**S(2) - q*x + x**S(2)), x)/a
    rule103 = ReplacementRule(pattern103, lambda a, C, x, A, B, b : With103(a, C, x, A, B, b))
    rubi.add(rule103)

    pattern104 = Pattern(Integral(x_*(x_*WC('C', S(1)) + WC('B', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, C: ZeroQ(B*(a/b)**(S(1)/3) - C*(a/b)**(S(2)/3))), )
    def With104(a, C, x, B, b):
        q = (a/b)**(S(1)/3)
        return C*q**S(3)*Int(x/(q**S(2) - q*x + x**S(2)), x)/a
    rule104 = ReplacementRule(pattern104, lambda a, C, x, B, b : With104(a, C, x, B, b))
    rubi.add(rule104)

    pattern105 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a, C, A: ZeroQ(A + C*(a/b)**(S(2)/3))), )
    def With105(a, C, x, A, b):
        q = (a/b)**(S(1)/3)
        return q**S(2)*Int((A + C*q*x)/(q**S(2) - q*x + x**S(2)), x)/a
    rule105 = ReplacementRule(pattern105, lambda a, C, x, A, b : With105(a, C, x, A, b))
    rubi.add(rule105)

    pattern106 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda a, C, A, B, b: ZeroQ(A + B*(-a/b)**(S(1)/3) + C*(-a/b)**(S(2)/3))), )
    def With106(a, C, x, A, B, b):
        q = (-a/b)**(S(1)/3)
        return q*Int((A*q + x*(A + B*q))/(q**S(2) + q*x + x**S(2)), x)/a
    rule106 = ReplacementRule(pattern106, lambda a, C, x, A, B, b : With106(a, C, x, A, B, b))
    rubi.add(rule106)

    pattern107 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, C: ZeroQ(B*(-a/b)**(S(1)/3) + C*(-a/b)**(S(2)/3))), )
    def With107(a, C, x, B, b):
        q = (-a/b)**(S(1)/3)
        return B*q**S(2)*Int(x/(q**S(2) + q*x + x**S(2)), x)/a
    rule107 = ReplacementRule(pattern107, lambda a, C, x, B, b : With107(a, C, x, B, b))
    rubi.add(rule107)

    pattern108 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a, C, A: ZeroQ(A + C*(-a/b)**(S(2)/3))), )
    def With108(a, C, x, A, b):
        q = (-a/b)**(S(1)/3)
        return A*q*Int((q + x)/(q**S(2) + q*x + x**S(2)), x)/a
    rule108 = ReplacementRule(pattern108, lambda a, C, x, A, b : With108(a, C, x, A, b))
    rubi.add(rule108)

    pattern109 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, A: NonzeroQ(-A**S(3)*b + B**S(3)*a)), CustomConstraint(lambda b, a: RationalQ(a/b)), CustomConstraint(lambda b, a: Greater(a/b, S(0))), CustomConstraint(lambda C, a, x, q, B, A: NonzeroQ(A - B*q + C*q**S(2))))
    def With109(a, C, x, A, B, b):
        q = (a/b)**(S(1)/3)
        return q*(A - B*q + C*q**S(2))*Int(1/(q + x), x)/(S(3)*a) + q*Int((q*(S(2)*A + B*q - C*q**S(2)) - x*(A - B*q - S(2)*C*q**S(2)))/(q**S(2) - q*x + x**S(2)), x)/(S(3)*a)
    rule109 = ReplacementRule(pattern109, lambda a, C, x, A, B, b : With109(a, C, x, A, B, b))
    rubi.add(rule109)

    pattern110 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a: RationalQ(a/b)), CustomConstraint(lambda b, a: Greater(a/b, S(0))), CustomConstraint(lambda C, a, x, q, B: NonzeroQ(B*q - C*q**S(2))))
    def With110(a, C, x, B, b):
        q = (a/b)**(S(1)/3)
        return -q*(B*q - C*q**S(2))*Int(1/(q + x), x)/(S(3)*a) + q*Int((q*(B*q - C*q**S(2)) + x*(B*q + S(2)*C*q**S(2)))/(q**S(2) - q*x + x**S(2)), x)/(S(3)*a)
    rule110 = ReplacementRule(pattern110, lambda a, C, x, B, b : With110(a, C, x, B, b))
    rubi.add(rule110)

    pattern111 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a: RationalQ(a/b)), CustomConstraint(lambda b, a: Greater(a/b, S(0))), CustomConstraint(lambda C, a, x, q, A: NonzeroQ(A + C*q**S(2))))
    def With111(a, C, x, A, b):
        q = (a/b)**(S(1)/3)
        return q*(A + C*q**S(2))*Int(1/(q + x), x)/(S(3)*a) + q*Int((q*(S(2)*A - C*q**S(2)) - x*(A - S(2)*C*q**S(2)))/(q**S(2) - q*x + x**S(2)), x)/(S(3)*a)
    rule111 = ReplacementRule(pattern111, lambda a, C, x, A, b : With111(a, C, x, A, b))
    rubi.add(rule111)

    pattern112 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda B, b, a, A: NonzeroQ(-A**S(3)*b + B**S(3)*a)), CustomConstraint(lambda b, a: RationalQ(a/b)), CustomConstraint(lambda b, a: Less(a/b, S(0))), CustomConstraint(lambda C, a, x, q, B, A: NonzeroQ(A + B*q + C*q**S(2))))
    def With112(a, C, x, A, B, b):
        q = (-a/b)**(S(1)/3)
        return q*(A + B*q + C*q**S(2))*Int(1/(q - x), x)/(S(3)*a) + q*Int((q*(S(2)*A - B*q - C*q**S(2)) + x*(A + B*q - S(2)*C*q**S(2)))/(q**S(2) + q*x + x**S(2)), x)/(S(3)*a)
    rule112 = ReplacementRule(pattern112, lambda a, C, x, A, B, b : With112(a, C, x, A, B, b))
    rubi.add(rule112)

    pattern113 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a: RationalQ(a/b)), CustomConstraint(lambda b, a: Less(a/b, S(0))), CustomConstraint(lambda C, a, x, q, B: NonzeroQ(B*q + C*q**S(2))))
    def With113(a, C, x, B, b):
        q = (-a/b)**(S(1)/3)
        return q*(B*q + C*q**S(2))*Int(1/(q - x), x)/(S(3)*a) + q*Int((-q*(B*q + C*q**S(2)) + x*(B*q - S(2)*C*q**S(2)))/(q**S(2) + q*x + x**S(2)), x)/(S(3)*a)
    rule113 = ReplacementRule(pattern113, lambda a, C, x, B, b : With113(a, C, x, B, b))
    rubi.add(rule113)

    pattern114 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, a: RationalQ(a/b)), CustomConstraint(lambda b, a: Less(a/b, S(0))), CustomConstraint(lambda C, a, x, q, A: NonzeroQ(A + C*q**S(2))))
    def With114(a, C, x, A, b):
        q = (-a/b)**(S(1)/3)
        return q*(A + C*q**S(2))*Int(1/(q - x), x)/(S(3)*a) + q*Int((q*(S(2)*A - C*q**S(2)) + x*(A - S(2)*C*q**S(2)))/(q**S(2) + q*x + x**S(2)), x)/(S(3)*a)
    rule114 = ReplacementRule(pattern114, lambda a, C, x, A, b : With114(a, C, x, A, b))
    rubi.add(rule114)

    pattern115 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2))), CustomConstraint(lambda n, Pq, x: Less(Expon(Pq, x), n)), CustomConstraint(lambda v, x: SumQ(v)))
    def With115(n, Pq, a, x, m, b, c):
        v = Sum(c**(-ii)*(c*x)**(ii + m)*(x**(n/S(2))*Coeff(Pq, x, ii + n/S(2)) + Coeff(Pq, x, ii))/(a + b*x**n), List(ii, S(0), n/S(2) + S(-1)))
        return Int(v, x)
    rule115 = ReplacementRule(pattern115, lambda n, Pq, a, x, m, b, c : With115(n, Pq, a, x, m, b, c))
    rubi.add(rule115)

    pattern116 = Pattern(Integral(Pq_/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2))), CustomConstraint(lambda n, Pq, x: Less(Expon(Pq, x), n)), CustomConstraint(lambda v, x: SumQ(v)))
    def With116(n, Pq, a, x, b):
        v = Sum(x**ii*(x**(n/S(2))*Coeff(Pq, x, ii + n/S(2)) + Coeff(Pq, x, ii))/(a + b*x**n), List(ii, S(0), n/S(2) + S(-1)))
        return Int(v, x)
    rule116 = ReplacementRule(pattern116, lambda n, Pq, a, x, b : With116(n, Pq, a, x, b))
    rubi.add(rule116)

    pattern117 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda a: PosQ(a)), CustomConstraint(lambda b, a, d, c: ZeroQ(c*Rt(b/a, S(3)) - d*(-sqrt(S(3)) + S(1)))), )
    def With117(a, x, b, d, c):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return -S(3)**(S(1)/4)*d*s*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(-sqrt(S(3)) + S(2))*(r*x + s)*EllipticE(asin((r*x + s*(-sqrt(S(3)) + S(1)))/(r*x + s*(S(1) + sqrt(S(3))))), S(-7) - S(4)*sqrt(S(3)))/(r**S(2)*sqrt(s*(r*x + s)/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(a + b*x**S(3))) + S(2)*d*s**S(3)*sqrt(a + b*x**S(3))/(a*r**S(2)*(r*x + s*(S(1) + sqrt(S(3)))))
    rule117 = ReplacementRule(pattern117, lambda a, x, b, d, c : With117(a, x, b, d, c))
    rubi.add(rule117)

    pattern118 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda a: PosQ(a)), CustomConstraint(lambda b, a, d, c: NonzeroQ(c*Rt(b/a, S(3)) - d*(-sqrt(S(3)) + S(1)))), )
    def With118(a, x, b, d, c):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return d*Int((r*x + s*(-sqrt(S(3)) + S(1)))/sqrt(a + b*x**S(3)), x)/r + (c*r - d*s*(-sqrt(S(3)) + S(1)))*Int(1/sqrt(a + b*x**S(3)), x)/r
    rule118 = ReplacementRule(pattern118, lambda a, x, b, d, c : With118(a, x, b, d, c))
    rubi.add(rule118)

    pattern119 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda a: NegQ(a)), CustomConstraint(lambda b, a, d, c: ZeroQ(c*Rt(b/a, S(3)) - d*(S(1) + sqrt(S(3))))), )
    def With119(a, x, b, d, c):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return S(3)**(S(1)/4)*d*s*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(-sqrt(S(3)) + S(1)))**S(2))*sqrt(sqrt(S(3)) + S(2))*(r*x + s)*EllipticE(asin((r*x + s*(S(1) + sqrt(S(3))))/(r*x + s*(-sqrt(S(3)) + S(1)))), S(-7) + S(4)*sqrt(S(3)))/(r**S(2)*sqrt(-s*(r*x + s)/(r*x + s*(-sqrt(S(3)) + S(1)))**S(2))*sqrt(a + b*x**S(3))) + S(2)*d*s**S(3)*sqrt(a + b*x**S(3))/(a*r**S(2)*(r*x + s*(-sqrt(S(3)) + S(1))))
    rule119 = ReplacementRule(pattern119, lambda a, x, b, d, c : With119(a, x, b, d, c))
    rubi.add(rule119)

    pattern120 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda a: NegQ(a)), CustomConstraint(lambda b, a, d, c: NonzeroQ(c*Rt(b/a, S(3)) - d*(S(1) + sqrt(S(3))))), )
    def With120(a, x, b, d, c):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return d*Int((r*x + s*(S(1) + sqrt(S(3))))/sqrt(a + b*x**S(3)), x)/r + (c*r - d*s*(S(1) + sqrt(S(3))))*Int(1/sqrt(a + b*x**S(3)), x)/r
    rule120 = ReplacementRule(pattern120, lambda a, x, b, d, c : With120(a, x, b, d, c))
    rubi.add(rule120)

    pattern121 = Pattern(Integral((c_ + x_**S(4)*WC('d', S(1)))/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, a, d, c: ZeroQ(S(2)*c*Rt(b/a, S(3))**S(2) - d*(-sqrt(S(3)) + S(1)))), )
    def With121(a, x, b, d, c):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return -S(3)**(S(1)/4)*d*s*x*sqrt((r**S(2)*x**S(4) - r*s*x**S(2) + s**S(2))/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*(r*x**S(2) + s)*EllipticE(acos((r*x**S(2)*(-sqrt(S(3)) + S(1)) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)), sqrt(S(3))/S(4) + S(1)/2)/(S(2)*r**S(2)*sqrt(r*x**S(2)*(r*x**S(2) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*sqrt(a + b*x**S(6))) + d*s**S(3)*x*(S(1) + sqrt(S(3)))*sqrt(a + b*x**S(6))/(S(2)*a*r**S(2)*(r*x**S(2)*(S(1) + sqrt(S(3))) + s))
    rule121 = ReplacementRule(pattern121, lambda a, x, b, d, c : With121(a, x, b, d, c))
    rubi.add(rule121)

    pattern122 = Pattern(Integral((c_ + x_**S(4)*WC('d', S(1)))/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, a, d, c: NonzeroQ(S(2)*c*Rt(b/a, S(3))**S(2) - d*(-sqrt(S(3)) + S(1)))), )
    def With122(a, x, b, d, c):
        q = Rt(b/a, S(3))
        return d*Int((S(2)*q**S(2)*x**S(4) - sqrt(S(3)) + S(1))/sqrt(a + b*x**S(6)), x)/(S(2)*q**S(2)) + (S(2)*c*q**S(2) - d*(-sqrt(S(3)) + S(1)))*Int(1/sqrt(a + b*x**S(6)), x)/(S(2)*q**S(2))
    rule122 = ReplacementRule(pattern122, lambda a, x, b, d, c : With122(a, x, b, d, c))
    rubi.add(rule122)

    pattern123 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, a, d, c: ZeroQ(-a*d**S(4) + b*c**S(4))))
    rule123 = ReplacementRule(pattern123, lambda a, x, b, d, c : -c*d*x**S(3)*sqrt(-(c - d*x**S(2))**S(2)/(c*d*x**S(2)))*sqrt(-d**S(2)*(a + b*x**S(8))/(b*c**S(2)*x**S(4)))*EllipticF(asin(sqrt((sqrt(S(2))*c**S(2) + S(2)*c*d*x**S(2) + sqrt(S(2))*d**S(2)*x**S(4))/(c*d*x**S(2)))/S(2)), S(-2) + S(2)*sqrt(S(2)))/(sqrt(sqrt(S(2)) + S(2))*sqrt(a + b*x**S(8))*(c - d*x**S(2))))
    rubi.add(rule123)

    pattern124 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, a, d, c: NonzeroQ(-a*d**S(4) + b*c**S(4))))
    rule124 = ReplacementRule(pattern124, lambda a, x, b, d, c : -(-c*Rt(b/a, S(4)) + d)*Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/(S(2)*Rt(b/a, S(4))) + (c*Rt(b/a, S(4)) + d)*Int((x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/(S(2)*Rt(b/a, S(4))))
    rubi.add(rule124)

    pattern125 = Pattern(Integral(Pq_/(x_*sqrt(a_ + x_**n_*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda Pq, x: NonzeroQ(Coeff(Pq, x, S(0)))))
    rule125 = ReplacementRule(pattern125, lambda n, Pq, a, x, b : Coeff(Pq, x, S(0))*Int(S(1)/(x*sqrt(a + b*x**n)), x) + Int(ExpandToSum((Pq - Coeff(Pq, x, S(0)))/x, x)/sqrt(a + b*x**n), x))
    rubi.add(rule125)

    pattern126 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2))), CustomConstraint(lambda n, Pq, x: Not(PolyQ(Pq, x**(n/S(2))))), )
    def With126(p, n, Pq, a, x, m, b, c):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        return Int(Sum(c**(-j)*(c*x)**(j + m)*(a + b*x**n)**p*Sum(x**(k*n/S(2))*Coeff(Pq, x, j + k*n/S(2)), List(k, S(0), S(1) + S(2)*(-j + q)/n)), List(j, S(0), n/S(2) + S(-1))), x)
    rule126 = ReplacementRule(pattern126, lambda p, n, Pq, a, x, m, b, c : With126(p, n, Pq, a, x, m, b, c))
    rubi.add(rule126)

    pattern127 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n/S(2))), CustomConstraint(lambda n, Pq, x: Not(PolyQ(Pq, x**(n/S(2))))), )
    def With127(p, n, Pq, a, x, b):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        return Int(Sum(x**j*(a + b*x**n)**p*Sum(x**(k*n/S(2))*Coeff(Pq, x, j + k*n/S(2)), List(k, S(0), S(1) + S(2)*(-j + q)/n)), List(j, S(0), n/S(2) + S(-1))), x)
    rule127 = ReplacementRule(pattern127, lambda p, n, Pq, a, x, b : With127(p, n, Pq, a, x, b))
    rubi.add(rule127)

    pattern128 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, Pq, x: Equal(Expon(Pq, x), n + S(-1))))
    rule128 = ReplacementRule(pattern128, lambda p, n, Pq, a, x, b : Coeff(Pq, x, n + S(-1))*Int(x**(n + S(-1))*(a + b*x**n)**p, x) + Int((a + b*x**n)**p*ExpandToSum(Pq - x**(n + S(-1))*Coeff(Pq, x, n + S(-1)), x), x))
    rubi.add(rule128)

    pattern129 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule129 = ReplacementRule(pattern129, lambda n, Pq, a, x, m, b, c : Int(ExpandIntegrand(Pq*(c*x)**m/(a + b*x**n), x), x))
    rubi.add(rule129)

    pattern130 = Pattern(Integral(Pq_/(a_ + x_**n_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: IntegerQ(n)))
    rule130 = ReplacementRule(pattern130, lambda n, Pq, a, x, b : Int(ExpandIntegrand(Pq/(a + b*x**n), x), x))
    rubi.add(rule130)

    pattern131 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda n, Pq, x: LessEqual(n + S(-1), Expon(Pq, x))), CustomConstraint(lambda b, n, Pq, m, Pq0, a, x, p, c: NonzeroQ(Pq0)))
    def With131(p, n, Pq, a, x, m, b, c):
        Pq0 = Coeff(Pq, x, S(0))
        return Pq0*(c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))) + Int((c*x)**(m + S(1))*(a + b*x**n)**p*ExpandToSum(-S(2)*Pq0*b*x**(n + S(-1))*(m + n*(p + S(1)) + S(1)) + S(2)*a*(Pq - Pq0)*(m + S(1))/x, x), x)/(S(2)*a*c*(m + S(1)))
    rule131 = ReplacementRule(pattern131, lambda p, n, Pq, a, x, m, b, c : With131(p, n, Pq, a, x, m, b, c))
    rubi.add(rule131)

    pattern132 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda Pqq, b, n, Pq, m, a, x, p, q, c: GreaterEqual(-n + q, S(0)) & NonzeroQ(m + n*p + q + S(1)) & (IntegerQ(S(2)*p) | IntegerQ(p + (q + S(1))/(S(2)*n)))))
    def With132(p, n, Pq, a, x, m, b, c):
        q = Expon(Pq, x)
        return With(List(Set(Pqq, Coeff(Pq, x, q))), Pqq*c**(n - q + S(-1))*(c*x)**(m - n + q + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*p + q + S(1))) + Int((c*x)**m*(a + b*x**n)**p*ExpandToSum(-Pqq*a*x**(-n + q)*(m - n + q + S(1)) + b*(Pq - Pqq*x**q)*(m + n*p + q + S(1)), x), x)/(b*(m + n*p + q + S(1))))
    rule132 = ReplacementRule(pattern132, lambda p, n, Pq, a, x, m, b, c : With132(p, n, Pq, a, x, m, b, c))
    rubi.add(rule132)

    pattern133 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda Pqq, b, n, Pq, a, x, p, q: GreaterEqual(-n + q, S(0)) & NonzeroQ(n*p + q + S(1)) & (IntegerQ(S(2)*p) | IntegerQ(p + (q + S(1))/(S(2)*n)))))
    def With133(p, n, Pq, a, x, b):
        q = Expon(Pq, x)
        return With(List(Set(Pqq, Coeff(Pq, x, q))), Pqq*x**(-n + q + S(1))*(a + b*x**n)**(p + S(1))/(b*(n*p + q + S(1))) + Int((a + b*x**n)**p*ExpandToSum(-Pqq*a*x**(-n + q)*(-n + q + S(1)) + b*(Pq - Pqq*x**q)*(n*p + q + S(1)), x), x)/(b*(n*p + q + S(1))))
    rule133 = ReplacementRule(pattern133, lambda p, n, Pq, a, x, b : With133(p, n, Pq, a, x, b))
    rubi.add(rule133)

    pattern134 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)), )
    def With134(p, n, Pq, a, x, m, b):
        q = Expon(Pq, x)
        return -Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, 1/x)), x), x), x, 1/x)
    rule134 = ReplacementRule(pattern134, lambda p, n, Pq, a, x, m, b : With134(p, n, Pq, a, x, m, b))
    rubi.add(rule134)

    pattern135 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)), )
    def With135(p, n, Pq, a, x, m, b, c):
        g = Denominator(m)
        q = Expon(Pq, x)
        return -g*Subst(Int(x**(-g*(m + q + S(1)) + S(-1))*(a + b*c**(-n)*x**(-g*n))**p*ExpandToSum(x**(g*q)*ReplaceAll(Pq, Rule(x, x**(-g)/c)), x), x), x, (c*x)**(-S(1)/g))/c
    rule135 = ReplacementRule(pattern135, lambda p, n, Pq, a, x, m, b, c : With135(p, n, Pq, a, x, m, b, c))
    rubi.add(rule135)

    pattern136 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: Not(RationalQ(m))), )
    def With136(p, n, Pq, a, x, m, b, c):
        q = Expon(Pq, x)
        return -(c*x)**m*(1/x)**m*Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, 1/x)), x), x), x, 1/x)
    rule136 = ReplacementRule(pattern136, lambda p, n, Pq, a, x, m, b, c : With136(p, n, Pq, a, x, m, b, c))
    rubi.add(rule136)

    pattern137 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: FractionQ(n)), )
    def With137(p, n, Pq, a, x, m, b):
        g = Denominator(n)
        return g*Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + b*x**(g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(1/g))
    rule137 = ReplacementRule(pattern137, lambda p, n, Pq, a, x, m, b : With137(p, n, Pq, a, x, m, b))
    rubi.add(rule137)

    pattern138 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: FractionQ(n)), )
    def With138(p, n, Pq, a, x, b):
        g = Denominator(n)
        return g*Subst(Int(x**(g + S(-1))*(a + b*x**(g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(1/g))
    rule138 = ReplacementRule(pattern138, lambda p, n, Pq, a, x, b : With138(p, n, Pq, a, x, b))
    rubi.add(rule138)

    pattern139 = Pattern(Integral(Pq_*(c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda n: FractionQ(n)))
    rule139 = ReplacementRule(pattern139, lambda p, n, Pq, a, x, m, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(Pq*x**m*(a + b*x**n)**p, x))
    rubi.add(rule139)

    pattern140 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda n, m: IntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule140 = ReplacementRule(pattern140, lambda p, n, Pq, a, x, m, b : Subst(Int((a + b*x**(n/(m + S(1))))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(n/(m + S(1))))), x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule140)

    pattern141 = Pattern(Integral(Pq_*(c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda n, m: IntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule141 = ReplacementRule(pattern141, lambda p, n, Pq, a, x, m, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(Pq*x**m*(a + b*x**n)**p, x))
    rubi.add(rule141)

    pattern142 = Pattern(Integral((A_ + x_**WC('m', S(1))*WC('B', S(1)))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))))
    rule142 = ReplacementRule(pattern142, lambda p, n, a, m, x, A, B, b : A*Int((a + b*x**n)**p, x) + B*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule142)

    pattern143 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x) | PolyQ(Pq, x**n)))
    rule143 = ReplacementRule(pattern143, lambda p, n, Pq, a, x, m, b, c : Int(ExpandIntegrand(Pq*(c*x)**m*(a + b*x**n)**p, x), x))
    rubi.add(rule143)

    pattern144 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x) | PolyQ(Pq, x**n)))
    rule144 = ReplacementRule(pattern144, lambda p, n, Pq, a, x, b : Int(ExpandIntegrand(Pq*(a + b*x**n)**p, x), x))
    rubi.add(rule144)

    pattern145 = Pattern(Integral(Pq_*u_**WC('m', S(1))*(a_ + v_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda v, x, u: LinearPairQ(u, v, x)), CustomConstraint(lambda n, Pq, v: PolyQ(Pq, v**n)))
    rule145 = ReplacementRule(pattern145, lambda p, n, Pq, a, v, x, u, m, b : u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p*SubstFor(v, Pq, x), x), x, v)/Coeff(v, x, S(1)))
    rubi.add(rule145)

    pattern146 = Pattern(Integral(Pq_*(a_ + v_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda n, Pq, v: PolyQ(Pq, v**n)))
    rule146 = ReplacementRule(pattern146, lambda p, n, Pq, a, v, x, b : Subst(Int((a + b*x**n)**p*SubstFor(v, Pq, x), x), x, v)/Coeff(v, x, S(1)))
    rubi.add(rule146)

    pattern147 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda a2, b2, a1, b1: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda a2, a1, p: IntegerQ(p) | (PositiveQ(a1) & PositiveQ(a2))))
    rule147 = ReplacementRule(pattern147, lambda b2, p, n, Pq, a2, a1, x, b1, m, c : Int(Pq*(c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule147)

    pattern148 = Pattern(Integral(Pq_*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda a2, b2, a1, b1: ZeroQ(a1*b2 + a2*b1)), CustomConstraint(lambda a2, a1, p: IntegerQ(p) | (PositiveQ(a1) & PositiveQ(a2))))
    rule148 = ReplacementRule(pattern148, lambda b2, p, n, Pq, a2, a1, x, b1 : Int(Pq*(a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule148)

    pattern149 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda a2, b2, a1, b1: ZeroQ(a1*b2 + a2*b1)))
    rule149 = ReplacementRule(pattern149, lambda b2, p, n, Pq, a2, a1, x, b1, m, c : (a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p))*Int(Pq*(c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule149)

    pattern150 = Pattern(Integral(Pq_*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a1, x: FreeQ(a1, x)), CustomConstraint(lambda b1, x: FreeQ(b1, x)), CustomConstraint(lambda a2, x: FreeQ(a2, x)), CustomConstraint(lambda b2, x: FreeQ(b2, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda a2, b2, a1, b1: ZeroQ(a1*b2 + a2*b1)))
    rule150 = ReplacementRule(pattern150, lambda b2, p, n, Pq, a2, a1, x, b1 : (a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p))*Int(Pq*(a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule150)

    pattern151 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)) + x_**WC('n2', S(1))*WC('g', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda p, n, a, f, b, d, e, c: ZeroQ(a*c*f - e*(a*d + b*c)*(n*(p + S(1)) + S(1)))), CustomConstraint(lambda e, p, c, n, a, b, d, g: ZeroQ(a*c*g - b*d*e*(S(2)*n*(p + S(1)) + S(1)))))
    rule151 = ReplacementRule(pattern151, lambda n2, e, p, c, n, a, f, x, b, d, g : e*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c))
    rubi.add(rule151)

    pattern152 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n2', S(1))*WC('g', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, p: ZeroQ(n*(p + S(1)) + S(1))), CustomConstraint(lambda e, p, c, n, a, b, d, g: ZeroQ(a*c*g - b*d*e*(S(2)*n*(p + S(1)) + S(1)))))
    rule152 = ReplacementRule(pattern152, lambda n2, e, p, c, n, a, x, b, d, g : e*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c))
    rubi.add(rule152)

    pattern153 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)) + x_**WC('n2', S(1))*WC('g', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda p, n, a, f, m, b, d, e, c: ZeroQ(a*c*f*(m + S(1)) - e*(a*d + b*c)*(m + n*(p + S(1)) + S(1)))), CustomConstraint(lambda e, p, c, n, a, m, b, d, g: ZeroQ(a*c*g*(m + S(1)) - b*d*e*(m + S(2)*n*(p + S(1)) + S(1)))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule153 = ReplacementRule(pattern153, lambda n2, e, p, c, n, a, f, x, m, b, h, d, g : e*(h*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c*h*(m + S(1))))
    rubi.add(rule153)

    pattern154 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n2', S(1))*WC('g', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, p, m: ZeroQ(m + n*(p + S(1)) + S(1))), CustomConstraint(lambda e, p, c, n, a, m, b, d, g: ZeroQ(a*c*g*(m + S(1)) - b*d*e*(m + S(2)*n*(p + S(1)) + S(1)))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule154 = ReplacementRule(pattern154, lambda n2, e, p, c, n, a, x, m, b, h, d, g : e*(h*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c*h*(m + S(1))))
    rubi.add(rule154)

    pattern155 = Pattern(Integral((A_ + x_**WC('m', S(1))*WC('B', S(1)))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda b, a, d, c: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))))
    rule155 = ReplacementRule(pattern155, lambda p, n, a, q, x, A, m, B, b, d, c : A*Int((a + b*x**n)**p*(c + d*x**n)**q, x) + B*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule155)

    pattern156 = Pattern(Integral(Px_**WC('q', S(1))*((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Px, x: PolynomialQ(Px, x)), CustomConstraint(lambda q: IntegerQ(q)), CustomConstraint(lambda n: RationalQ(n)), )
    def With156(p, n, a, Px, q, x, b, d, c):
        k = Denominator(n)
        return k*Subst(Int(SimplifyIntegrand(x**(k + S(-1))*(a + b*x**(k*n))**p*ReplaceAll(Px, Rule(x, -c/d + x**k/d))**q, x), x), x, (c + d*x)**(1/k))/d
    rule156 = ReplacementRule(pattern156, lambda p, n, a, Px, q, x, b, d, c : With156(p, n, a, Px, q, x, b, d, c))
    rubi.add(rule156)

    pattern157 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))))
    rule157 = ReplacementRule(pattern157, lambda n2, p, n, Pq, a, x, m, b, c : Subst(Int((a + b*x + c*x**S(2))**p*SubstFor(x**n, Pq, x), x), x, x**n)/n)
    rubi.add(rule157)

    pattern158 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule158 = ReplacementRule(pattern158, lambda n2, p, n, Pq, a, x, m, b, d, c : Int(ExpandIntegrand(Pq*(d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x))
    rubi.add(rule158)

    pattern159 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)))
    rule159 = ReplacementRule(pattern159, lambda n2, p, n, Pq, a, x, b, c : Int(ExpandIntegrand(Pq*(a + b*x**n + c*x**(S(2)*n))**p, x), x))
    rubi.add(rule159)

    pattern160 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('n', S(1))*WC('e', S(1)) + x_**WC('n2', S(1))*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda e, p, n, a, b, d: ZeroQ(a*e - b*d*(n*(p + S(1)) + S(1)))), CustomConstraint(lambda p, n, a, f, d, c: ZeroQ(a*f - c*d*(S(2)*n*(p + S(1)) + S(1)))))
    rule160 = ReplacementRule(pattern160, lambda n2, p, n, a, f, x, b, d, e, c : d*x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/a)
    rubi.add(rule160)

    pattern161 = Pattern(Integral((d_ + x_**WC('n2', S(1))*WC('f', S(1)))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, p: ZeroQ(n*(p + S(1)) + S(1))), CustomConstraint(lambda f, a, d, c: ZeroQ(a*f + c*d)))
    rule161 = ReplacementRule(pattern161, lambda n2, p, n, a, f, x, b, d, c : d*x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/a)
    rubi.add(rule161)

    pattern162 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('n', S(1))*WC('e', S(1)) + x_**WC('n2', S(1))*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda e, p, n, a, m, b, d: ZeroQ(a*e*(m + S(1)) - b*d*(m + n*(p + S(1)) + S(1)))), CustomConstraint(lambda p, n, a, f, m, d, c: ZeroQ(a*f*(m + S(1)) - c*d*(m + S(2)*n*(p + S(1)) + S(1)))), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule162 = ReplacementRule(pattern162, lambda n2, p, n, a, g, f, x, m, b, d, e, c : d*(g*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*g*(m + S(1))))
    rubi.add(rule162)

    pattern163 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(d_ + x_**WC('n2', S(1))*WC('f', S(1)))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, p, m: ZeroQ(m + n*(p + S(1)) + S(1))), CustomConstraint(lambda f, a, d, c: ZeroQ(a*f + c*d)), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule163 = ReplacementRule(pattern163, lambda n2, p, n, a, g, f, x, m, b, d, c : d*(g*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*g*(m + S(1))))
    rubi.add(rule163)

    pattern164 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: ZeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: Not(IntegerQ(S(2)*p))))
    rule164 = ReplacementRule(pattern164, lambda n2, p, n, Pq, a, x, m, b, d, c : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p)*Int(Pq*(d*x)**m*(b + S(2)*c*x**n)**(S(2)*p), x))
    rubi.add(rule164)

    pattern165 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: ZeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: Not(IntegerQ(S(2)*p))))
    rule165 = ReplacementRule(pattern165, lambda n2, p, n, Pq, a, x, b, c : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p)*Int(Pq*(b + S(2)*c*x**n)**(S(2)*p), x))
    rubi.add(rule165)

    pattern166 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)))
    rule166 = ReplacementRule(pattern166, lambda n2, p, n, Pq, a, x, m, b, c : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x + c*x**S(2))**p*SubstFor(x**n, Pq, x), x), x, x**n)/n)
    rubi.add(rule166)

    pattern167 = Pattern(Integral(Pq_*(d_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)))
    rule167 = ReplacementRule(pattern167, lambda n2, p, n, Pq, a, x, m, b, d, c : x**(-m)*(d*x)**m*Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x))
    rubi.add(rule167)

    pattern168 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda Pq, x: ZeroQ(Coeff(Pq, x, S(0)))))
    rule168 = ReplacementRule(pattern168, lambda n2, p, n, Pq, a, x, m, b, d, c : Int((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p*ExpandToSum(Pq/x, x), x)/d)
    rubi.add(rule168)

    pattern169 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda Pq, x: ZeroQ(Coeff(Pq, x, S(0)))), CustomConstraint(lambda Pq: SumQ(Pq)))
    rule169 = ReplacementRule(pattern169, lambda n2, p, n, Pq, a, x, b, c : Int(x*(a + b*x**n + c*x**(S(2)*n))**p*ExpandToSum(Pq/x, x), x))
    rubi.add(rule169)

    pattern170 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)) + x_**WC('n2', S(1))*WC('f', S(1)) + x_**WC('n3', S(1))*WC('g', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, n3: ZeroQ(-S(3)*n + n3)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p, n, a, g, b, d, e, c: ZeroQ(a**S(2)*g*(n + S(1)) - c*(a*e - b*d*(n*(p + S(1)) + S(1)))*(n*(S(2)*p + S(3)) + S(1)))), CustomConstraint(lambda e, p, n, a, f, b, d, c: ZeroQ(a**S(2)*f*(n + S(1)) - a*c*d*(n + S(1))*(S(2)*n*(p + S(1)) + S(1)) - b*(a*e - b*d*(n*(p + S(1)) + S(1)))*(n*(p + S(2)) + S(1)))))
    rule170 = ReplacementRule(pattern170, lambda n2, n3, p, n, a, g, f, x, b, d, e, c : x*(S(3)*a*d - x**S(2)*(-a*e + S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)))
    rubi.add(rule170)

    pattern171 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('n2', S(1))*WC('f', S(1)) + x_**WC('n3', S(1))*WC('g', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, n3: ZeroQ(-S(3)*n + n3)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p, n, a, g, b, d, c: ZeroQ(a**S(2)*g*(n + S(1)) + b*c*d*(n*(p + S(1)) + S(1))*(n*(S(2)*p + S(3)) + S(1)))), CustomConstraint(lambda p, n, a, f, b, d, c: ZeroQ(a**S(2)*f*(n + S(1)) - a*c*d*(n + S(1))*(S(2)*n*(p + S(1)) + S(1)) + b**S(2)*d*(n*(p + S(1)) + S(1))*(n*(p + S(2)) + S(1)))))
    rule171 = ReplacementRule(pattern171, lambda n2, n3, p, n, a, g, f, x, b, d, c : x*(S(3)*a*d - x**S(2)*(S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)))
    rubi.add(rule171)

    pattern172 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)) + x_**WC('n3', S(1))*WC('g', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, n3: ZeroQ(-S(3)*n + n3)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p, n, a, g, b, d, e, c: ZeroQ(a**S(2)*g*(n + S(1)) - c*(a*e - b*d*(n*(p + S(1)) + S(1)))*(n*(S(2)*p + S(3)) + S(1)))), CustomConstraint(lambda e, p, n, a, b, d, c: ZeroQ(a*c*d*(n + S(1))*(S(2)*n*(p + S(1)) + S(1)) + b*(a*e - b*d*(n*(p + S(1)) + S(1)))*(n*(p + S(2)) + S(1)))))
    rule172 = ReplacementRule(pattern172, lambda n2, n3, p, c, n, a, x, b, d, e, g : x*(S(3)*a*d - x**S(2)*(-a*e + S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)))
    rubi.add(rule172)

    pattern173 = Pattern(Integral((d_ + x_**WC('n3', S(1))*WC('g', S(1)))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, n3: ZeroQ(-S(3)*n + n3)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p, n, a, g, b, d, c: ZeroQ(a**S(2)*g*(n + S(1)) + b*c*d*(n*(p + S(1)) + S(1))*(n*(S(2)*p + S(3)) + S(1)))), CustomConstraint(lambda p, n, a, b, d, c: ZeroQ(a*c*d*(n + S(1))*(S(2)*n*(p + S(1)) + S(1)) - b**S(2)*d*(n*(p + S(1)) + S(1))*(n*(p + S(2)) + S(1)))))
    rule173 = ReplacementRule(pattern173, lambda n2, n3, p, n, a, g, x, b, d, c : x*(S(3)*a*d - x**S(2)*(S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)))
    rubi.add(rule173)

    pattern174 = Pattern(Integral(x_**WC('m', S(1))*(e_ + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)) + x_**WC('s', S(1))*WC('h', S(1)))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, q: ZeroQ(-n/S(2) + q)), CustomConstraint(lambda n, r: ZeroQ(-S(3)*n/S(2) + r)), CustomConstraint(lambda n, s: ZeroQ(-S(2)*n + s)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n, m: ZeroQ(S(2)*m - n + S(2))), CustomConstraint(lambda h, a, e, c: ZeroQ(a*h + c*e)))
    rule174 = ReplacementRule(pattern174, lambda n2, r, n, a, q, g, f, s, x, m, b, h, e, c : (-S(2)*c*x**n*(-b*g + S(2)*c*f) - S(2)*c*(-S(2)*a*g + b*f) - S(2)*h*x**(n/S(2))*(-S(4)*a*c + b**S(2)))/(c*n*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**n + c*x**(S(2)*n))))
    rubi.add(rule174)

    pattern175 = Pattern(Integral((d_*x_)**WC('m', S(1))*(e_ + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)) + x_**WC('s', S(1))*WC('h', S(1)))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, q: ZeroQ(-n/S(2) + q)), CustomConstraint(lambda n, r: ZeroQ(-S(3)*n/S(2) + r)), CustomConstraint(lambda n, s: ZeroQ(-S(2)*n + s)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n, m: ZeroQ(S(2)*m - n + S(2))), CustomConstraint(lambda h, a, e, c: ZeroQ(a*h + c*e)))
    rule175 = ReplacementRule(pattern175, lambda n2, r, e, n, a, q, g, f, s, x, m, b, h, d, c : x**(-m)*(d*x)**m*Int(x**m*(e + f*x**(n/S(2)) + g*x**(S(3)*n/S(2)) + h*x**(S(2)*n))/(a + b*x**n + c*x**(S(2)*n))**(S(3)/2), x))
    rubi.add(rule175)

    pattern176 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda b, n, Pq, a, i, x, p, q, c: Less(q, S(2)*n)))
    def With176(n2, p, n, Pq, a, x, b, c):
        q = Expon(Pq, x)
        i = Symbol('i')
        return -x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Sum(c*x**(i + n)*(-S(2)*a*Coeff(Pq, x, i + n) + b*Coeff(Pq, x, i)) + x**i*(-a*b*Coeff(Pq, x, i + n) + (-S(2)*a*c + b**S(2))*Coeff(Pq, x, i)), List(i, S(0), n + S(-1)))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))) + Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Sum(c*x**(i + n)*(-S(2)*a*Coeff(Pq, x, i + n) + b*Coeff(Pq, x, i))*(i + n*(S(2)*p + S(3)) + S(1)) + x**i*(-a*b*(i + S(1))*Coeff(Pq, x, i + n) + (-S(2)*a*c*(i + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(i + n*(p + S(1)) + S(1)))*Coeff(Pq, x, i)), List(i, S(0), n + S(-1))), x)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2)))
    rule176 = ReplacementRule(pattern176, lambda n2, p, n, Pq, a, x, b, c : With176(n2, p, n, Pq, a, x, b, c))
    rubi.add(rule176)

    pattern177 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda g, a, d, c: ZeroQ(a*g + c*d)))
    rule177 = ReplacementRule(pattern177, lambda c, a, f, x, b, d, e, g : (-c*x**S(2)*(-b*f + S(2)*c*e) - c*(-S(2)*a*f + b*e) - g*x*(-S(4)*a*c + b**S(2)))/(c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))))
    rubi.add(rule177)

    pattern178 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda g, a, d, c: ZeroQ(a*g + c*d)))
    rule178 = ReplacementRule(pattern178, lambda c, a, f, x, b, d, g : (S(2)*a*c*f + b*c*f*x**S(2) - g*x*(-S(4)*a*c + b**S(2)))/(c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))))
    rubi.add(rule178)

    pattern179 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda g, a, d, c: ZeroQ(a*g + c*d)))
    rule179 = ReplacementRule(pattern179, lambda c, a, x, b, d, e, g : (-b*c*e - S(2)*c**S(2)*e*x**S(2) - g*x*(-S(4)*a*c + b**S(2)))/(c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))))
    rubi.add(rule179)

    pattern180 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + x_**S(2)*WC('g', S(1)) + x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda h, a, e, c: ZeroQ(-S(3)*a*h + c*e)), CustomConstraint(lambda g, h, b, c: ZeroQ(-S(2)*b*h + c*g)))
    rule180 = ReplacementRule(pattern180, lambda c, a, f, x, b, h, e, g : (S(2)*a**S(2)*c*f + a*b*c*f*x**S(2) + a*h*x**S(3)*(-S(4)*a*c + b**S(2)))/(a*c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))))
    rubi.add(rule180)

    pattern181 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + x_**S(2)*WC('g', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda h, a, e, c: ZeroQ(-S(3)*a*h + c*e)), CustomConstraint(lambda g, h, b, c: ZeroQ(-S(2)*b*h + c*g)))
    rule181 = ReplacementRule(pattern181, lambda c, a, x, b, h, e, g : h*x**S(3)/(c*sqrt(a + b*x**S(2) + c*x**S(4))))
    rubi.add(rule181)

    pattern182 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda h, a, e, c: ZeroQ(-S(3)*a*h + c*e)), CustomConstraint(lambda a, g, b, d, e, c: ZeroQ(S(3)*a*g - S(2)*b*e + S(3)*c*d)))
    rule182 = ReplacementRule(pattern182, lambda c, a, f, x, b, d, h, e, g : (S(2)*a**S(2)*c*f + a*b*c*f*x**S(2) + a*h*x**S(3)*(-S(4)*a*c + b**S(2)) + c*d*x*(-S(4)*a*c + b**S(2)))/(a*c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))))
    rubi.add(rule182)

    pattern183 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(3)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda h, a, e, c: ZeroQ(-S(3)*a*h + c*e)), CustomConstraint(lambda b, e, d, c: ZeroQ(-S(2)*b*e + S(3)*c*d)))
    rule183 = ReplacementRule(pattern183, lambda a, f, x, b, d, h, e, c : (S(2)*a**S(2)*c*f + a*b*c*f*x**S(2) + a*h*x**S(3)*(-S(4)*a*c + b**S(2)) + c*d*x*(-S(4)*a*c + b**S(2)))/(a*c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))))
    rubi.add(rule183)

    pattern184 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda b, n, Q, Pq, a, i, x, p, q, c, R: GreaterEqual(q, S(2)*n)))
    def With184(n2, p, n, Pq, a, x, b, c):
        q = Expon(Pq, x)
        return Module(List(Set(Q, PolynomialQuotient(Pq*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)), Set(R, PolynomialRemainder(Pq*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)), i), -x*(b*c)**(-Floor((q + S(-1))/n) + S(-1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Sum(c*x**(i + n)*(-S(2)*a*Coeff(R, x, i + n) + b*Coeff(R, x, i)) + x**i*(-a*b*Coeff(R, x, i + n) + (-S(2)*a*c + b**S(2))*Coeff(R, x, i)), List(i, S(0), n + S(-1)))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))) + (b*c)**(-Floor((q + S(-1))/n) + S(-1))*Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*ExpandToSum(Q*a*n*(p + S(1))*(-S(4)*a*c + b**S(2)) + Sum(c*x**(i + n)*(-S(2)*a*Coeff(R, x, i + n) + b*Coeff(R, x, i))*(i + n*(S(2)*p + S(3)) + S(1)) + x**i*(-a*b*(i + S(1))*Coeff(R, x, i + n) + (-S(2)*a*c*(i + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(i + n*(p + S(1)) + S(1)))*Coeff(R, x, i)), List(i, S(0), n + S(-1))), x), x)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rule184 = ReplacementRule(pattern184, lambda n2, p, n, Pq, a, x, b, c : With184(n2, p, n, Pq, a, x, b, c))
    rubi.add(rule184)

    pattern185 = Pattern(Integral(Pq_*x_**m_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m: NegativeIntegerQ(m)), CustomConstraint(lambda b, n, Q, Pq, m, a, i, x, p, q, c, R: GreaterEqual(q, S(2)*n)))
    def With185(n2, p, n, Pq, a, x, m, b, c):
        q = Expon(Pq, x)
        return Module(List(Set(Q, PolynomialQuotient(Pq*a*x**m*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)), Set(R, PolynomialRemainder(Pq*a*x**m*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)), i), (b*c)**(-Floor((q + S(-1))/n) + S(-1))*Int(x**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*ExpandToSum(Q*n*x**(-m)*(p + S(1))*(-S(4)*a*c + b**S(2)) + Sum(c*x**(i - m + n)*(-S(2)*Coeff(R, x, i + n) + b*Coeff(R, x, i)/a)*(i + n*(S(2)*p + S(3)) + S(1)) + x**(i - m)*(-b*(i + S(1))*Coeff(R, x, i + n) + (-S(2)*c*(i + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(i + n*(p + S(1)) + S(1))/a)*Coeff(R, x, i)), List(i, S(0), n + S(-1))), x), x)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))) - x*(b*c)**(-Floor((q + S(-1))/n) + S(-1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Sum(c*x**(i + n)*(-S(2)*a*Coeff(R, x, i + n) + b*Coeff(R, x, i)) + x**i*(-a*b*Coeff(R, x, i + n) + (-S(2)*a*c + b**S(2))*Coeff(R, x, i)), List(i, S(0), n + S(-1)))/(a**S(2)*n*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rule185 = ReplacementRule(pattern185, lambda n2, p, n, Pq, a, x, m, b, c : With185(n2, p, n, Pq, a, x, m, b, c))
    rubi.add(rule185)

    pattern186 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda b, n, Pq, m, a, x, p, g, c: Unequal(g, S(1))))
    def With186(n2, p, n, Pq, a, x, m, b, c):
        g = GCD(m + S(1), n)
        return Subst(Int(x**(S(-1) + (m + S(1))/g)*(a + b*x**(n/g) + c*x**(S(2)*n/g))**p*ReplaceAll(Pq, Rule(x, x**(1/g))), x), x, x**g)/g
    rule186 = ReplacementRule(pattern186, lambda n2, p, n, Pq, a, x, m, b, c : With186(n2, p, n, Pq, a, x, m, b, c))
    rubi.add(rule186)

    pattern187 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda b, a, c: NiceSqrtQ(-S(4)*a*c + b**S(2))))
    rule187 = ReplacementRule(pattern187, lambda n2, n, Pq, a, x, m, b, d, c : Int(ExpandIntegrand(Pq*(d*x)**m/(a + b*x**n + c*x**(S(2)*n)), x), x))
    rubi.add(rule187)

    pattern188 = Pattern(Integral(Pq_/(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, Pq, a, x, b, c: Less(Expon(Pq, x), n) | NiceSqrtQ(-S(4)*a*c + b**S(2))))
    rule188 = ReplacementRule(pattern188, lambda n2, n, Pq, a, x, b, c : Int(ExpandIntegrand(Pq/(a + b*x**n + c*x**(S(2)*n)), x), x))
    rubi.add(rule188)

    pattern189 = Pattern(Integral(Pq_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda Pqq, b, Pq, a, x, p, q, c: Equal(S(2)*p + q + S(1), S(0))))
    def With189(p, Pq, a, x, b, c):
        q = Expon(Pq, x)
        return With(List(Set(Pqq, Coeff(Pq, x, q))), Pqq*c**p*log(a + b*x + c*x**S(2))/S(2) + Int((a + b*x + c*x**S(2))**p*ExpandToSum(S(2)*Pq - Pqq*c**p*(b + S(2)*c*x)*(a + b*x + c*x**S(2))**(-p + S(-1)), x), x)/S(2))
    rule189 = ReplacementRule(pattern189, lambda p, Pq, a, x, b, c : With189(p, Pq, a, x, b, c))
    rubi.add(rule189)

    pattern190 = Pattern(Integral(Pq_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda c: PosQ(c)), CustomConstraint(lambda Pqq, b, Pq, a, x, p, q, c: Equal(S(2)*p + q + S(1), S(0))))
    def With190(p, Pq, a, x, b, c):
        q = Expon(Pq, x)
        return With(List(Set(Pqq, Coeff(Pq, x, q))), Pqq*c**p*atanh((b + S(2)*c*x)/(S(2)*sqrt(a + b*x + c*x**S(2))*Rt(c, S(2)))) + Int((a + b*x + c*x**S(2))**p*ExpandToSum(Pq - Pqq*c**(p + S(1)/2)*(a + b*x + c*x**S(2))**(-p + S(-1)/2), x), x))
    rule190 = ReplacementRule(pattern190, lambda p, Pq, a, x, b, c : With190(p, Pq, a, x, b, c))
    rubi.add(rule190)

    pattern191 = Pattern(Integral(Pq_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda c: NegQ(c)), CustomConstraint(lambda Pqq, b, Pq, a, x, p, q, c: Equal(S(2)*p + q + S(1), S(0))))
    def With191(p, Pq, a, x, b, c):
        q = Expon(Pq, x)
        return With(List(Set(Pqq, Coeff(Pq, x, q))), -Pqq*(-c)**p*atan((b + S(2)*c*x)/(S(2)*sqrt(a + b*x + c*x**S(2))*Rt(-c, S(2)))) + Int((a + b*x + c*x**S(2))**p*ExpandToSum(Pq - Pqq*(-c)**(p + S(1)/2)*(a + b*x + c*x**S(2))**(-p + S(-1)/2), x), x))
    rule191 = ReplacementRule(pattern191, lambda p, Pq, a, x, b, c : With191(p, Pq, a, x, b, c))
    rubi.add(rule191)

    pattern192 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda d, Pqq, b, n, Pq, m, a, x, p, q, c: GreaterEqual(q, S(2)*n) & Unequal(m + S(2)*n*p + q + S(1), S(0)) & (IntegerQ(S(2)*p) | (Equal(n, S(1)) & IntegerQ(S(4)*p)) | IntegerQ(p + (q + S(1))/(S(2)*n)))))
    def With192(n2, p, n, Pq, a, x, m, b, d, c):
        q = Expon(Pq, x)
        return With(List(Set(Pqq, Coeff(Pq, x, q))), Pqq*d**(S(2)*n - q + S(-1))*(d*x)**(m - S(2)*n + q + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(c*(m + S(2)*n*p + q + S(1))) + Int((d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p*ExpandToSum(Pq - Pqq*x**q - Pqq*(a*x**(-S(2)*n + q)*(m - S(2)*n + q + S(1)) + b*x**(-n + q)*(m + n*(p + S(-1)) + q + S(1)))/(c*(m + S(2)*n*p + q + S(1))), x), x))
    rule192 = ReplacementRule(pattern192, lambda n2, p, n, Pq, a, x, m, b, d, c : With192(n2, p, n, Pq, a, x, m, b, d, c))
    rubi.add(rule192)

    pattern193 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda Pqq, b, n, Pq, a, x, p, q, c: GreaterEqual(q, S(2)*n) & Unequal(S(2)*n*p + q + S(1), S(0)) & (IntegerQ(S(2)*p) | (Equal(n, S(1)) & IntegerQ(S(4)*p)) | IntegerQ(p + (q + S(1))/(S(2)*n)))))
    def With193(n2, p, n, Pq, a, x, b, c):
        q = Expon(Pq, x)
        return With(List(Set(Pqq, Coeff(Pq, x, q))), Pqq*x**(-S(2)*n + q + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(c*(S(2)*n*p + q + S(1))) + Int((a + b*x**n + c*x**(S(2)*n))**p*ExpandToSum(Pq - Pqq*x**q - Pqq*(a*x**(-S(2)*n + q)*(-S(2)*n + q + S(1)) + b*x**(-n + q)*(n*(p + S(-1)) + q + S(1)))/(c*(S(2)*n*p + q + S(1))), x), x))
    rule193 = ReplacementRule(pattern193, lambda n2, p, n, Pq, a, x, b, c : With193(n2, p, n, Pq, a, x, b, c))
    rubi.add(rule193)

    pattern194 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, Pq, x: Not(PolyQ(Pq, x**n))), )
    def With194(n2, p, n, Pq, a, x, m, b, d, c):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        return Int(Sum(d**(-j)*(d*x)**(j + m)*(a + b*x**n + c*x**(S(2)*n))**p*Sum(x**(k*n)*Coeff(Pq, x, j + k*n), List(k, S(0), S(1) + (-j + q)/n)), List(j, S(0), n + S(-1))), x)
    rule194 = ReplacementRule(pattern194, lambda n2, p, n, Pq, a, x, m, b, d, c : With194(n2, p, n, Pq, a, x, m, b, d, c))
    rubi.add(rule194)

    pattern195 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)), CustomConstraint(lambda n, Pq, x: Not(PolyQ(Pq, x**n))), )
    def With195(n2, p, n, Pq, a, x, b, c):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        return Int(Sum(x**j*(a + b*x**n + c*x**(S(2)*n))**p*Sum(x**(k*n)*Coeff(Pq, x, j + k*n), List(k, S(0), S(1) + (-j + q)/n)), List(j, S(0), n + S(-1))), x)
    rule195 = ReplacementRule(pattern195, lambda n2, p, n, Pq, a, x, b, c : With195(n2, p, n, Pq, a, x, b, c))
    rubi.add(rule195)

    pattern196 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule196 = ReplacementRule(pattern196, lambda n2, n, Pq, a, x, m, b, d, c : Int(RationalFunctionExpand(Pq*(d*x)**m/(a + b*x**n + c*x**(S(2)*n)), x), x))
    rubi.add(rule196)

    pattern197 = Pattern(Integral(Pq_/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: PositiveIntegerQ(n)))
    rule197 = ReplacementRule(pattern197, lambda n2, n, Pq, a, x, b, c : Int(RationalFunctionExpand(Pq/(a + b*x**n + c*x**(S(2)*n)), x), x))
    rubi.add(rule197)

    pattern198 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: IntegerQ(m)), )
    def With198(n2, p, n, Pq, a, x, m, b, c):
        q = Expon(Pq, x)
        return -Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, 1/x)), x), x), x, 1/x)
    rule198 = ReplacementRule(pattern198, lambda n2, p, n, Pq, a, x, m, b, c : With198(n2, p, n, Pq, a, x, m, b, c))
    rubi.add(rule198)

    pattern199 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: FractionQ(m)), )
    def With199(n2, p, n, Pq, a, x, m, b, d, c):
        g = Denominator(m)
        q = Expon(Pq, x)
        return -g*Subst(Int(x**(-g*(m + q + S(1)) + S(-1))*(a + b*d**(-n)*x**(-g*n) + c*d**(-S(2)*n)*x**(-S(2)*g*n))**p*ExpandToSum(x**(g*q)*ReplaceAll(Pq, Rule(x, x**(-g)/d)), x), x), x, (d*x)**(-S(1)/g))/d
    rule199 = ReplacementRule(pattern199, lambda n2, p, n, Pq, a, x, m, b, d, c : With199(n2, p, n, Pq, a, x, m, b, d, c))
    rubi.add(rule199)

    pattern200 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: NegativeIntegerQ(n)), CustomConstraint(lambda m: Not(RationalQ(m))), )
    def With200(n2, p, n, Pq, a, x, m, b, d, c):
        q = Expon(Pq, x)
        return -(d*x)**m*(1/x)**m*Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, 1/x)), x), x), x, 1/x)
    rule200 = ReplacementRule(pattern200, lambda n2, p, n, Pq, a, x, m, b, d, c : With200(n2, p, n, Pq, a, x, m, b, d, c))
    rubi.add(rule200)

    pattern201 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: FractionQ(n)), )
    def With201(n2, p, n, Pq, a, x, m, b, c):
        g = Denominator(n)
        return g*Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + b*x**(g*n) + c*x**(S(2)*g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(1/g))
    rule201 = ReplacementRule(pattern201, lambda n2, p, n, Pq, a, x, m, b, c : With201(n2, p, n, Pq, a, x, m, b, c))
    rubi.add(rule201)

    pattern202 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: FractionQ(n)), )
    def With202(n2, p, n, Pq, a, x, b, c):
        g = Denominator(n)
        return g*Subst(Int(x**(g + S(-1))*(a + b*x**(g*n) + c*x**(S(2)*g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(1/g))
    rule202 = ReplacementRule(pattern202, lambda n2, p, n, Pq, a, x, b, c : With202(n2, p, n, Pq, a, x, b, c))
    rubi.add(rule202)

    pattern203 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: FractionQ(n)), CustomConstraint(lambda m: PositiveIntegerQ(m + S(1)/2)))
    rule203 = ReplacementRule(pattern203, lambda n2, p, n, Pq, a, x, m, b, d, c : d**(m + S(-1)/2)*sqrt(d*x)*Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x)/sqrt(x))
    rubi.add(rule203)

    pattern204 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: FractionQ(n)), CustomConstraint(lambda m: NegativeIntegerQ(m + S(-1)/2)))
    rule204 = ReplacementRule(pattern204, lambda n2, p, n, Pq, a, x, m, b, d, c : d**(m + S(1)/2)*sqrt(x)*Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x)/sqrt(d*x))
    rubi.add(rule204)

    pattern205 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n: FractionQ(n)))
    rule205 = ReplacementRule(pattern205, lambda n2, p, n, Pq, a, x, m, b, d, c : x**(-m)*(d*x)**m*Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x))
    rubi.add(rule205)

    pattern206 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n, m: IntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule206 = ReplacementRule(pattern206, lambda n2, p, n, Pq, a, x, m, b, c : Subst(Int((a + b*x**(n/(m + S(1))) + c*x**(S(2)*n/(m + S(1))))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(n/(m + S(1))))), x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule206)

    pattern207 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda n, m: IntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule207 = ReplacementRule(pattern207, lambda n2, p, n, Pq, a, x, m, b, d, c : x**(-m)*(d*x)**m*Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x))
    rubi.add(rule207)

    pattern208 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), )
    def With208(n2, n, Pq, a, x, m, b, d, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int(Pq*(d*x)**m/(b + S(2)*c*x**n - q), x)/q - S(2)*c*Int(Pq*(d*x)**m/(b + S(2)*c*x**n + q), x)/q
    rule208 = ReplacementRule(pattern208, lambda n2, n, Pq, a, x, m, b, d, c : With208(n2, n, Pq, a, x, m, b, d, c))
    rubi.add(rule208)

    pattern209 = Pattern(Integral(Pq_/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), )
    def With209(n2, n, Pq, a, x, b, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int(Pq/(b + S(2)*c*x**n - q), x)/q - S(2)*c*Int(Pq/(b + S(2)*c*x**n + q), x)/q
    rule209 = ReplacementRule(pattern209, lambda n2, n, Pq, a, x, b, c : With209(n2, n, Pq, a, x, b, c))
    rubi.add(rule209)

    pattern210 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1))))
    rule210 = ReplacementRule(pattern210, lambda n2, p, n, Pq, a, x, m, b, d, c : Int(ExpandIntegrand(Pq*(d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x))
    rubi.add(rule210)

    pattern211 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1))))
    rule211 = ReplacementRule(pattern211, lambda n2, p, n, Pq, a, x, b, c : Int(ExpandIntegrand(Pq*(a + b*x**n + c*x**(S(2)*n))**p, x), x))
    rubi.add(rule211)

    pattern212 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x) | PolyQ(Pq, x**n)))
    rule212 = ReplacementRule(pattern212, lambda n2, p, n, Pq, a, x, m, b, d, c : Int(Pq*(d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p, x))
    rubi.add(rule212)

    pattern213 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x) | PolyQ(Pq, x**n)))
    rule213 = ReplacementRule(pattern213, lambda n2, p, n, Pq, a, x, b, c : Int(Pq*(a + b*x**n + c*x**(S(2)*n))**p, x))
    rubi.add(rule213)

    pattern214 = Pattern(Integral(Pq_*u_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda v, x, u: LinearPairQ(u, v, x)), CustomConstraint(lambda n, Pq, v: PolyQ(Pq, v**n)))
    rule214 = ReplacementRule(pattern214, lambda n2, p, n, Pq, a, v, x, u, m, b, c : u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p*SubstFor(v, Pq, x), x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule214)

    pattern215 = Pattern(Integral(Pq_*(a_ + v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda v, x: LinearQ(v, x)), CustomConstraint(lambda n, Pq, v: PolyQ(Pq, v**n)))
    rule215 = ReplacementRule(pattern215, lambda n2, p, n, Pq, a, v, x, b, c : Subst(Int((a + b*x**n + c*x**(S(2)*n))**p*SubstFor(v, Pq, x), x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule215)

    pattern216 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, p, j: ZeroQ(j*p + j - n + S(1))))
    rule216 = ReplacementRule(pattern216, lambda p, n, a, x, j, b : x**(-n + S(1))*(a*x**j + b*x**n)**(p + S(1))/(b*(-j + n)*(p + S(1))))
    rubi.add(rule216)

    pattern217 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, p, j: NegativeIntegerQ((j - n*p - n + S(-1))/(j - n))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule217 = ReplacementRule(pattern217, lambda p, n, a, x, j, b : -x**(-j + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))) + (-j + n*p + n + S(1))*Int(x**(-j)*(a*x**j + b*x**n)**(p + S(1)), x)/(a*(-j + n)*(p + S(1))))
    rubi.add(rule217)

    pattern218 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, p, j: NegativeIntegerQ((j - n*p - n + S(-1))/(j - n))), CustomConstraint(lambda p, j: NonzeroQ(j*p + S(1))))
    rule218 = ReplacementRule(pattern218, lambda p, n, a, x, j, b : -b*(-j + n*p + n + S(1))*Int(x**(-j + n)*(a*x**j + b*x**n)**p, x)/(a*(j*p + S(1))) + x**(-j + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(j*p + S(1))))
    rubi.add(rule218)

    pattern219 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, p, j: RationalQ(j, n, p)), CustomConstraint(lambda n, j: Less(S(0), j, n)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, j: Less(j*p + S(1), S(0))))
    rule219 = ReplacementRule(pattern219, lambda p, n, a, x, j, b : -b*p*(-j + n)*Int(x**n*(a*x**j + b*x**n)**(p + S(-1)), x)/(j*p + S(1)) + x*(a*x**j + b*x**n)**p/(j*p + S(1)))
    rubi.add(rule219)

    pattern220 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, p, j: RationalQ(j, n, p)), CustomConstraint(lambda n, j: Less(S(0), j, n)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n, p: NonzeroQ(n*p + S(1))))
    rule220 = ReplacementRule(pattern220, lambda p, n, a, x, j, b : a*p*(-j + n)*Int(x**j*(a*x**j + b*x**n)**(p + S(-1)), x)/(n*p + S(1)) + x*(a*x**j + b*x**n)**p/(n*p + S(1)))
    rubi.add(rule220)

    pattern221 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, p, j: RationalQ(j, n, p)), CustomConstraint(lambda n, j: Less(S(0), j, n)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, p, j: Greater(j*p + S(1), -j + n)))
    rule221 = ReplacementRule(pattern221, lambda p, n, a, x, j, b : x**(-n + S(1))*(a*x**j + b*x**n)**(p + S(1))/(b*(-j + n)*(p + S(1))) - (j*p + j - n + S(1))*Int(x**(-n)*(a*x**j + b*x**n)**(p + S(1)), x)/(b*(-j + n)*(p + S(1))))
    rubi.add(rule221)

    pattern222 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, p, j: RationalQ(j, n, p)), CustomConstraint(lambda n, j: Less(S(0), j, n)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule222 = ReplacementRule(pattern222, lambda p, n, a, x, j, b : -x**(-j + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))) + (-j + n*p + n + S(1))*Int(x**(-j)*(a*x**j + b*x**n)**(p + S(1)), x)/(a*(-j + n)*(p + S(1))))
    rubi.add(rule222)

    pattern223 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda p, j: ZeroQ(j*p + S(1))))
    rule223 = ReplacementRule(pattern223, lambda p, n, a, x, j, b : a*Int(x**j*(a*x**j + b*x**n)**(p + S(-1)), x) + x*(a*x**j + b*x**n)**p/(p*(-j + n)))
    rubi.add(rule223)

    pattern224 = Pattern(Integral(1/sqrt(x_**S(2)*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n: NonzeroQ(n + S(-2))))
    rule224 = ReplacementRule(pattern224, lambda n, b, a, x : S(2)*Subst(Int(1/(-a*x**S(2) + S(1)), x), x, x/sqrt(a*x**S(2) + b*x**n))/(-n + S(2)))
    rubi.add(rule224)

    pattern225 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda p, j: ZeroQ(j*p + S(1))))
    rule225 = ReplacementRule(pattern225, lambda p, n, a, x, j, b : -x**(-j + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))) + (-j + n*p + n + S(1))*Int(x**(-j)*(a*x**j + b*x**n)**(p + S(1)), x)/(a*(-j + n)*(p + S(1))))
    rubi.add(rule225)

    pattern226 = Pattern(Integral(1/sqrt(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda n, j: RationalQ(j, n)), CustomConstraint(lambda n, j: Less(S(2)*n + S(-2), j, n)))
    rule226 = ReplacementRule(pattern226, lambda n, a, x, j, b : -a*(-j + S(2)*n + S(-2))*Int(x**(j - n)/sqrt(a*x**j + b*x**n), x)/(b*(n + S(-2))) - S(2)*x**(-n + S(1))*sqrt(a*x**j + b*x**n)/(b*(n + S(-2))))
    rubi.add(rule226)

    pattern227 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: PosQ(-j + n)))
    rule227 = ReplacementRule(pattern227, lambda p, n, a, x, j, b : x**(-j*FracPart(p))*(a + b*x**(-j + n))**(-FracPart(p))*(a*x**j + b*x**n)**FracPart(p)*Int(x**(j*p)*(a + b*x**(-j + n))**p, x))
    rubi.add(rule227)

    pattern228 = Pattern(Integral((u_**WC('j', S(1))*WC('a', S(1)) + u_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda x, u: NonzeroQ(u - x)))
    rule228 = ReplacementRule(pattern228, lambda p, n, a, x, u, j, b : Subst(Int((a*x**j + b*x**n)**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule228)

    pattern229 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, m: ZeroQ(m - n + S(1))))
    rule229 = ReplacementRule(pattern229, lambda p, n, a, m, x, j, b : Subst(Int((a*x**(j/n) + b*x)**p, x), x, x**n)/n)
    rubi.add(rule229)

    pattern230 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, p, j, m: ZeroQ(-j + m + n*p + n + S(1))), CustomConstraint(lambda c, j: IntegerQ(j) | PositiveQ(c)))
    rule230 = ReplacementRule(pattern230, lambda p, j, n, a, x, m, b, c : -c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))))
    rubi.add(rule230)

    pattern231 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, p, j, m: NegativeIntegerQ((j - m - n*p - n + S(-1))/(j - n))), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda c, j: IntegerQ(j) | PositiveQ(c)))
    rule231 = ReplacementRule(pattern231, lambda p, j, n, a, x, m, b, c : c**j*(-j + m + n*p + n + S(1))*Int((c*x)**(-j + m)*(a*x**j + b*x**n)**(p + S(1)), x)/(a*(-j + n)*(p + S(1))) - c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))))
    rubi.add(rule231)

    pattern232 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, p, j, m: NegativeIntegerQ((j - m - n*p - n + S(-1))/(j - n))), CustomConstraint(lambda p, j, m: NonzeroQ(j*p + m + S(1))), CustomConstraint(lambda n, c, j: PositiveQ(c) | IntegersQ(j, n)))
    rule232 = ReplacementRule(pattern232, lambda p, j, n, a, x, m, b, c : -b*c**(j - n)*(-j + m + n*p + n + S(1))*Int((c*x)**(-j + m + n)*(a*x**j + b*x**n)**p, x)/(a*(j*p + m + S(1))) + c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(j*p + m + S(1))))
    rubi.add(rule232)

    pattern233 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, p, j, m: NegativeIntegerQ((j - m - n*p - n + S(-1))/(j - n))))
    rule233 = ReplacementRule(pattern233, lambda p, n, a, m, x, j, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a*x**j + b*x**n)**p, x))
    rubi.add(rule233)

    pattern234 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)), CustomConstraint(lambda n: NonzeroQ(n**S(2) + S(-1))))
    rule234 = ReplacementRule(pattern234, lambda p, n, a, m, x, j, b : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a*x**(j/n) + b*x)**p, x), x, x**n)/n)
    rubi.add(rule234)

    pattern235 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)), CustomConstraint(lambda n: NonzeroQ(n**S(2) + S(-1))))
    rule235 = ReplacementRule(pattern235, lambda p, n, a, m, x, j, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a*x**j + b*x**n)**p, x))
    rubi.add(rule235)

    pattern236 = Pattern(Integral((x_*WC('c', S(1)))**m_*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, m, p, j: RationalQ(j, m, n, p)), CustomConstraint(lambda n, j: Less(S(0), j, n)), CustomConstraint(lambda n, c, j: PositiveQ(c) | IntegersQ(j, n)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda p, j, m: Less(j*p + m + S(1), S(0))))
    rule236 = ReplacementRule(pattern236, lambda p, n, a, m, x, j, b, c : -b*c**(-n)*p*(-j + n)*Int((c*x)**(m + n)*(a*x**j + b*x**n)**(p + S(-1)), x)/(j*p + m + S(1)) + (c*x)**(m + S(1))*(a*x**j + b*x**n)**p/(c*(j*p + m + S(1))))
    rubi.add(rule236)

    pattern237 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, p, j: RationalQ(j, n, p)), CustomConstraint(lambda n, j: Less(S(0), j, n)), CustomConstraint(lambda n, c, j: PositiveQ(c) | IntegersQ(j, n)), CustomConstraint(lambda p: Greater(p, S(0))), CustomConstraint(lambda n, p, m: NonzeroQ(m + n*p + S(1))))
    rule237 = ReplacementRule(pattern237, lambda p, j, n, a, x, m, b, c : a*c**(-j)*p*(-j + n)*Int((c*x)**(j + m)*(a*x**j + b*x**n)**(p + S(-1)), x)/(m + n*p + S(1)) + (c*x)**(m + S(1))*(a*x**j + b*x**n)**p/(c*(m + n*p + S(1))))
    rubi.add(rule237)

    pattern238 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, m, p, j: RationalQ(j, m, n, p)), CustomConstraint(lambda n, j: Less(S(0), j, n)), CustomConstraint(lambda n, c, j: PositiveQ(c) | IntegersQ(j, n)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda n, p, j, m: Greater(j*p + m + S(1), -j + n)))
    rule238 = ReplacementRule(pattern238, lambda p, j, n, a, x, m, b, c : -c**n*(j*p + j + m - n + S(1))*Int((c*x)**(m - n)*(a*x**j + b*x**n)**(p + S(1)), x)/(b*(-j + n)*(p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a*x**j + b*x**n)**(p + S(1))/(b*(-j + n)*(p + S(1))))
    rubi.add(rule238)

    pattern239 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, p, j: RationalQ(j, n, p)), CustomConstraint(lambda n, j: Less(S(0), j, n)), CustomConstraint(lambda n, c, j: PositiveQ(c) | IntegersQ(j, n)), CustomConstraint(lambda p: Less(p, S(-1))))
    rule239 = ReplacementRule(pattern239, lambda p, j, n, a, x, m, b, c : c**j*(-j + m + n*p + n + S(1))*Int((c*x)**(-j + m)*(a*x**j + b*x**n)**(p + S(1)), x)/(a*(-j + n)*(p + S(1))) - c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))))
    rubi.add(rule239)

    pattern240 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: RationalQ(j, n)), CustomConstraint(lambda n, j: Less(S(0), j, n)), CustomConstraint(lambda n, c, j: PositiveQ(c) | IntegersQ(j, n)), CustomConstraint(lambda n, p, j, m: PositiveQ(j*p + j + m - n + S(1))), CustomConstraint(lambda n, p, m: NonzeroQ(m + n*p + S(1))))
    rule240 = ReplacementRule(pattern240, lambda p, j, n, a, x, m, b, c : -a*c**(-j + n)*(j*p + j + m - n + S(1))*Int((c*x)**(j + m - n)*(a*x**j + b*x**n)**p, x)/(b*(m + n*p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a*x**j + b*x**n)**(p + S(1))/(b*(m + n*p + S(1))))
    rubi.add(rule240)

    pattern241 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: RationalQ(j, n)), CustomConstraint(lambda n, j: Less(S(0), j, n)), CustomConstraint(lambda n, c, j: PositiveQ(c) | IntegersQ(j, n)), CustomConstraint(lambda p, j, m: NegativeQ(j*p + m + S(1))))
    rule241 = ReplacementRule(pattern241, lambda p, j, n, a, x, m, b, c : -b*c**(j - n)*(-j + m + n*p + n + S(1))*Int((c*x)**(-j + m + n)*(a*x**j + b*x**n)**p, x)/(a*(j*p + m + S(1))) + c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(j*p + m + S(1))))
    rubi.add(rule241)

    pattern242 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: IntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule242 = ReplacementRule(pattern242, lambda p, n, a, m, x, j, b : Subst(Int((a*x**(j/(m + S(1))) + b*x**(n/(m + S(1))))**p, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule242)

    pattern243 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: IntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule243 = ReplacementRule(pattern243, lambda p, n, a, m, x, j, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a*x**j + b*x**n)**p, x))
    rubi.add(rule243)

    pattern244 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: PositiveIntegerQ(p + S(1)/2)), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda p, j, m: ZeroQ(j*p + m + S(1))), CustomConstraint(lambda c, j: IntegerQ(j) | PositiveQ(c)))
    rule244 = ReplacementRule(pattern244, lambda p, j, n, a, x, m, b, c : a*c**(-j)*Int((c*x)**(j + m)*(a*x**j + b*x**n)**(p + S(-1)), x) + (c*x)**(m + S(1))*(a*x**j + b*x**n)**p/(c*p*(-j + n)))
    rubi.add(rule244)

    pattern245 = Pattern(Integral(x_**WC('m', S(1))/sqrt(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda j, m: ZeroQ(-j/S(2) + m + S(1))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)))
    rule245 = ReplacementRule(pattern245, lambda n, a, m, x, j, b : -S(2)*Subst(Int(1/(-a*x**S(2) + S(1)), x), x, x**(j/S(2))/sqrt(a*x**j + b*x**n))/(-j + n))
    rubi.add(rule245)

    pattern246 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1)/2)), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda p, j, m: ZeroQ(j*p + m + S(1))), CustomConstraint(lambda c, j: IntegerQ(j) | PositiveQ(c)))
    rule246 = ReplacementRule(pattern246, lambda p, j, n, a, x, m, b, c : c**j*(-j + m + n*p + n + S(1))*Int((c*x)**(-j + m)*(a*x**j + b*x**n)**(p + S(1)), x)/(a*(-j + n)*(p + S(1))) - c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))))
    rubi.add(rule246)

    pattern247 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: IntegerQ(p + S(1)/2)), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda p, j, m: ZeroQ(j*p + m + S(1))))
    rule247 = ReplacementRule(pattern247, lambda p, n, a, m, x, j, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a*x**j + b*x**n)**p, x))
    rubi.add(rule247)

    pattern248 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: PosQ(-j + n)))
    rule248 = ReplacementRule(pattern248, lambda p, j, n, a, x, m, b, c : c**IntPart(m)*x**(-j*FracPart(p) - FracPart(m))*(c*x)**FracPart(m)*(a + b*x**(-j + n))**(-FracPart(p))*(a*x**j + b*x**n)**FracPart(p)*Int(x**(j*p + m)*(a + b*x**(-j + n))**p, x))
    rubi.add(rule248)

    pattern249 = Pattern(Integral(u_**WC('m', S(1))*(v_**WC('j', S(1))*WC('a', S(1)) + v_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda v, x, u: LinearPairQ(u, v, x)))
    rule249 = ReplacementRule(pattern249, lambda p, n, a, m, v, x, u, j, b : u**m*v**(-m)*Subst(Int(x**m*(a*x**j + b*x**n)**p, x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule249)

    pattern250 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda k, x: FreeQ(k, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda k, j: NonzeroQ(-j + k)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, k: IntegerQ(k/n)), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)), CustomConstraint(lambda n: NonzeroQ(n**S(2) + S(-1))))
    rule250 = ReplacementRule(pattern250, lambda p, j, n, a, k, q, x, m, b, d, c : Subst(Int(x**(S(-1) + (m + S(1))/n)*(c + d*x)**q*(a*x**(j/n) + b*x**(k/n))**p, x), x, x**n)/n)
    rubi.add(rule250)

    pattern251 = Pattern(Integral((e_*x_)**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda k, x: FreeQ(k, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda k, j: NonzeroQ(-j + k)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, k: IntegerQ(k/n)), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)), CustomConstraint(lambda n: NonzeroQ(n**S(2) + S(-1))))
    rule251 = ReplacementRule(pattern251, lambda e, p, j, n, k, a, q, x, m, b, d, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(c + d*x**n)**q*(a*x**j + b*x**k)**p, x))
    rubi.add(rule251)

    pattern252 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, jn, j: ZeroQ(jn - j - n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, a, d, c: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda p, n, a, m, j, b, d, c: ZeroQ(a*d*(j*p + m + S(1)) - b*c*(m + n + p*(j + n) + S(1)))), CustomConstraint(lambda e, j: IntegersQ(j) | PositiveQ(e)), CustomConstraint(lambda p, j, m: NonzeroQ(j*p + m + S(1))))
    rule252 = ReplacementRule(pattern252, lambda p, n, a, m, jn, x, j, b, d, e, c : c*e**(j + S(-1))*(e*x)**(-j + m + S(1))*(a*x**j + b*x**(j + n))**(p + S(1))/(a*(j*p + m + S(1))))
    rubi.add(rule252)

    pattern253 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, jn, j: ZeroQ(jn - j - n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, a, d, c: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda m, p, j: RationalQ(j, m, p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda m, j: Inequality(S(0), Less, j, LessEqual, m)), CustomConstraint(lambda e, j: IntegerQ(j) | PositiveQ(e)))
    rule253 = ReplacementRule(pattern253, lambda p, n, a, m, jn, x, j, b, d, e, c : -e**j*(a*d*(j*p + m + S(1)) - b*c*(m + n + p*(j + n) + S(1)))*Int((e*x)**(-j + m)*(a*x**j + b*x**(j + n))**(p + S(1)), x)/(a*b*n*(p + S(1))) - e**(j + S(-1))*(e*x)**(-j + m + S(1))*(-a*d + b*c)*(a*x**j + b*x**(j + n))**(p + S(1))/(a*b*n*(p + S(1))))
    rubi.add(rule253)

    pattern254 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, jn, j: ZeroQ(jn - j - n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, a, d, c: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, m: RationalQ(m, n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda n, p, j, m: Less(j*p + m, S(-1)) | (Less(p, S(0)) & IntegersQ(m + S(-1)/2, p + S(-1)/2) & Less(m, -n*p + S(-1)))), CustomConstraint(lambda n, e, j: PositiveQ(e) | IntegersQ(j, n)), CustomConstraint(lambda p, j, m: NonzeroQ(j*p + m + S(1))), CustomConstraint(lambda n, p, j, m: NonzeroQ(j*p + m - n + S(1))))
    rule254 = ReplacementRule(pattern254, lambda p, n, a, m, jn, x, j, b, d, e, c : c*e**(j + S(-1))*(e*x)**(-j + m + S(1))*(a*x**j + b*x**(j + n))**(p + S(1))/(a*(j*p + m + S(1))) + e**(-n)*(a*d*(j*p + m + S(1)) - b*c*(m + n + p*(j + n) + S(1)))*Int((e*x)**(m + n)*(a*x**j + b*x**(j + n))**p, x)/(a*(j*p + m + S(1))))
    rubi.add(rule254)

    pattern255 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, jn, j: ZeroQ(jn - j - n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, a, d, c: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, p, j, m: NonzeroQ(m + n + p*(j + n) + S(1))), CustomConstraint(lambda e, j: IntegerQ(j) | PositiveQ(e)))
    rule255 = ReplacementRule(pattern255, lambda p, n, a, m, jn, x, j, b, d, e, c : d*e**(j + S(-1))*(e*x)**(-j + m + S(1))*(a*x**j + b*x**(j + n))**(p + S(1))/(b*(m + n + p*(j + n) + S(1))) - (a*d*(j*p + m + S(1)) - b*c*(m + n + p*(j + n) + S(1)))*Int((e*x)**m*(a*x**j + b*x**(j + n))**p, x)/(b*(m + n + p*(j + n) + S(1))))
    rubi.add(rule255)

    pattern256 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda k, x: FreeQ(k, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda k, j: NonzeroQ(-j + k)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, k: IntegerQ(k/n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: IntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule256 = ReplacementRule(pattern256, lambda p, j, n, k, a, q, x, m, b, d, c : Subst(Int((c + d*x**(n/(m + S(1))))**q*(a*x**(j/(m + S(1))) + b*x**(k/(m + S(1))))**p, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule256)

    pattern257 = Pattern(Integral((e_*x_)**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda k, x: FreeQ(k, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda k, j: NonzeroQ(-j + k)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, k: IntegerQ(k/n)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: IntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule257 = ReplacementRule(pattern257, lambda e, p, j, n, k, a, q, x, m, b, d, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(c + d*x**n)**q*(a*x**j + b*x**k)**p, x))
    rubi.add(rule257)

    pattern258 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda q, x: FreeQ(q, x)), CustomConstraint(lambda n, jn, j: ZeroQ(jn - j - n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, a, d, c: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda n, j: Not(ZeroQ(j + S(-1)) & ZeroQ(n + S(-1)))))
    rule258 = ReplacementRule(pattern258, lambda p, n, a, q, m, jn, x, j, b, d, e, c : e**IntPart(m)*x**(-j*FracPart(p) - FracPart(m))*(e*x)**FracPart(m)*(a + b*x**n)**(-FracPart(p))*(a*x**j + b*x**(j + n))**FracPart(p)*Int(x**(j*p + m)*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule258)

    pattern259 = Pattern(Integral(Pq_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: RationalQ(j, n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n: Less(S(-1), n, S(1))), )
    def With259(p, n, Pq, a, x, j, b):
        d = Denominator(n)
        return d*Subst(Int(x**(d + S(-1))*(a*x**(d*j) + b*x**(d*n))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(d*n))), x), x, x**(1/d))
    rule259 = ReplacementRule(pattern259, lambda p, n, Pq, a, x, j, b : With259(p, n, Pq, a, x, j, b))
    rubi.add(rule259)

    pattern260 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)))
    rule260 = ReplacementRule(pattern260, lambda p, n, Pq, a, m, x, j, b : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a*x**(j/n) + b*x)**p*SubstFor(x**n, Pq, x), x), x, x**n)/n)
    rubi.add(rule260)

    pattern261 = Pattern(Integral(Pq_*(c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m**S(2), S(1))))
    rule261 = ReplacementRule(pattern261, lambda p, n, Pq, a, m, x, j, b, c : c**(Quotient(m, sign(m))*sign(m))*x**(-Mod(m, sign(m)))*(c*x)**Mod(m, sign(m))*Int(Pq*x**m*(a*x**j + b*x**n)**p, x))
    rubi.add(rule261)

    pattern262 = Pattern(Integral(Pq_*(c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)))
    rule262 = ReplacementRule(pattern262, lambda p, n, Pq, a, m, x, j, b, c : x**(-m)*(c*x)**m*Int(Pq*x**m*(a*x**j + b*x**n)**p, x))
    rubi.add(rule262)

    pattern263 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: PositiveIntegerQ(j, n, j/n)), CustomConstraint(lambda m: IntegerQ(m)), CustomConstraint(lambda b, n, Pq, m, j, a, x, p, g: Unequal(g, S(1))))
    def With263(p, n, Pq, a, m, x, j, b):
        g = GCD(m + S(1), n)
        return Subst(Int(x**(S(-1) + (m + S(1))/g)*(a*x**(j/g) + b*x**(n/g))**p*ReplaceAll(Pq, Rule(x, x**(1/g))), x), x, x**g)/g
    rule263 = ReplacementRule(pattern263, lambda p, n, Pq, a, m, x, j, b : With263(p, n, Pq, a, m, x, j, b))
    rubi.add(rule263)

    pattern264 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda Pq, x: PolyQ(Pq, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: PositiveIntegerQ(j, n)), CustomConstraint(lambda n, j: Less(j, n)), CustomConstraint(lambda Pqq, b, n, Pq, m, j, a, x, p, q, c: Greater(q, n + S(-1)) & Unequal(m + n*p + q + S(1), S(0)) & (IntegerQ(S(2)*p) | IntegerQ(p + (q + S(1))/(S(2)*n)))))
    def With264(p, j, n, Pq, a, x, m, b, c):
        q = Expon(Pq, x)
        return With(List(Set(Pqq, Coeff(Pq, x, q))), Pqq*c**(n - q + S(-1))*(c*x)**(m - n + q + S(1))*(a*x**j + b*x**n)**(p + S(1))/(b*(m + n*p + q + S(1))) + Int((c*x)**m*(a*x**j + b*x**n)**p*ExpandToSum(Pq - Pqq*a*x**(-n + q)*(m - n + q + S(1))/(b*(m + n*p + q + S(1))) - Pqq*x**q, x), x))
    rule264 = ReplacementRule(pattern264, lambda p, j, n, Pq, a, x, m, b, c : With264(p, j, n, Pq, a, x, m, b, c))
    rubi.add(rule264)

    pattern265 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, m: IntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule265 = ReplacementRule(pattern265, lambda p, n, Pq, a, m, x, j, b : Subst(Int((a*x**(j/(m + S(1))) + b*x**(n/(m + S(1))))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(n/(m + S(1))))), x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule265)

    pattern266 = Pattern(Integral(Pq_*(c_*x_)**m_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, m: IntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Greater(m**S(2), S(1))))
    rule266 = ReplacementRule(pattern266, lambda p, n, Pq, a, m, x, j, b, c : c**(Quotient(m, sign(m))*sign(m))*x**(-Mod(m, sign(m)))*(c*x)**Mod(m, sign(m))*Int(Pq*x**m*(a*x**j + b*x**n)**p, x))
    rubi.add(rule266)

    pattern267 = Pattern(Integral(Pq_*(c_*x_)**m_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x**n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)), CustomConstraint(lambda n, j: IntegerQ(j/n)), CustomConstraint(lambda n, m: IntegerQ(n/(m + S(1)))), CustomConstraint(lambda n: Not(IntegerQ(n))))
    rule267 = ReplacementRule(pattern267, lambda p, n, Pq, a, m, x, j, b, c : x**(-m)*(c*x)**m*Int(Pq*x**m*(a*x**j + b*x**n)**p, x))
    rubi.add(rule267)

    pattern268 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x) | PolyQ(Pq, x**n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)))
    rule268 = ReplacementRule(pattern268, lambda p, j, n, Pq, a, x, m, b, c : Int(ExpandIntegrand(Pq*(c*x)**m*(a*x**j + b*x**n)**p, x), x))
    rubi.add(rule268)

    pattern269 = Pattern(Integral(Pq_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda n, Pq, x: PolyQ(Pq, x) | PolyQ(Pq, x**n)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda n, j: NonzeroQ(-j + n)))
    rule269 = ReplacementRule(pattern269, lambda p, n, Pq, a, x, j, b : Int(ExpandIntegrand(Pq*(a*x**j + b*x**n)**p, x), x))
    rubi.add(rule269)

    pattern270 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda b, a, d: ZeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))))
    rule270 = ReplacementRule(pattern270, lambda p, a, x, b, d : S(3)**(-S(3)*p)*a**(-S(2)*p)*Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p), x))
    rubi.add(rule270)

    pattern271 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda b, a, d: NonzeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))))
    rule271 = ReplacementRule(pattern271, lambda p, a, x, b, d : Int(ExpandToSum((a + b*x + d*x**S(3))**p, x), x))
    rubi.add(rule271)

    pattern272 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda b, a, d: NonzeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))), CustomConstraint(lambda x, p, u: ProductQ(NonfreeFactors(u, x))))
    def With272(p, a, x, b, d):
        u = Factor(a + b*x + d*x**S(3))
        return FreeFactors(u, x)**p*Int(DistributeDegree(NonfreeFactors(u, x), p), x)
    rule272 = ReplacementRule(pattern272, lambda p, a, x, b, d : With272(p, a, x, b, d))
    rubi.add(rule272)

    pattern273 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda b, a, d: NonzeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))), )
    def With273(p, a, x, b, d):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        return S(3)**(-S(3)*p)*d**(-S(2)*p)*Int((-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(-sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x)
    rule273 = ReplacementRule(pattern273, lambda p, a, x, b, d : With273(p, a, x, b, d))
    rubi.add(rule273)

    pattern274 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, a, d: ZeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))))
    rule274 = ReplacementRule(pattern274, lambda p, a, x, b, d : (S(3)*a - b*x)**(-p)*(S(3)*a + S(2)*b*x)**(-S(2)*p)*(a + b*x + d*x**S(3))**p*Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p), x))
    rubi.add(rule274)

    pattern275 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, a, d: NonzeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))), CustomConstraint(lambda d, a, b, x, p, u: ProductQ(u)))
    def With275(p, a, x, b, d):
        u = NonfreeFactors(Factor(a + b*x + d*x**S(3)), x)
        return (a + b*x + d*x**S(3))**p*Int(DistributeDegree(u, p), x)/DistributeDegree(u, p)
    rule275 = ReplacementRule(pattern275, lambda p, a, x, b, d : With275(p, a, x, b, d))
    rubi.add(rule275)

    pattern276 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, a, d: NonzeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))), )
    def With276(p, a, x, b, d):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        return (-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(-sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(a + b*x + d*x**S(3))**p*Int((-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(-sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x)
    rule276 = ReplacementRule(pattern276, lambda p, a, x, b, d : With276(p, a, x, b, d))
    rubi.add(rule276)

    pattern277 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda b, a, d: ZeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))))
    rule277 = ReplacementRule(pattern277, lambda p, a, f, x, m, b, d, e : S(3)**(-S(3)*p)*a**(-S(2)*p)*Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p)*(e + f*x)**m, x))
    rubi.add(rule277)

    pattern278 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda b, a, d: NonzeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))))
    rule278 = ReplacementRule(pattern278, lambda p, a, f, x, m, b, d, e : Int(ExpandIntegrand((e + f*x)**m*(a + b*x + d*x**S(3))**p, x), x))
    rubi.add(rule278)

    pattern279 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda b, a, d: NonzeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))), CustomConstraint(lambda f, e, x, p, u, m: ProductQ(NonfreeFactors(u, x))))
    def With279(p, a, f, x, m, b, d, e):
        u = Factor(a + b*x + d*x**S(3))
        return FreeFactors(u, x)**p*Int((e + f*x)**m*DistributeDegree(NonfreeFactors(u, x), p), x)
    rule279 = ReplacementRule(pattern279, lambda p, a, f, x, m, b, d, e : With279(p, a, f, x, m, b, d, e))
    rubi.add(rule279)

    pattern280 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda b, a, d: NonzeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))), )
    def With280(p, a, f, x, m, b, d, e):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        return S(3)**(-S(3)*p)*d**(-S(2)*p)*Int((e + f*x)**m*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(-sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x)
    rule280 = ReplacementRule(pattern280, lambda p, a, f, x, m, b, d, e : With280(p, a, f, x, m, b, d, e))
    rubi.add(rule280)

    pattern281 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, a, d: ZeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))))
    rule281 = ReplacementRule(pattern281, lambda p, a, f, x, m, b, d, e : (S(3)*a - b*x)**(-p)*(S(3)*a + S(2)*b*x)**(-S(2)*p)*(a + b*x + d*x**S(3))**p*Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p)*(e + f*x)**m, x))
    rubi.add(rule281)

    pattern282 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, a, d: NonzeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))), CustomConstraint(lambda f, d, b, m, e, a, x, p, u: ProductQ(u)))
    def With282(p, a, f, x, m, b, d, e):
        u = NonfreeFactors(Factor(a + b*x + d*x**S(3)), x)
        return (a + b*x + d*x**S(3))**p*Int((e + f*x)**m*DistributeDegree(u, p), x)/DistributeDegree(u, p)
    rule282 = ReplacementRule(pattern282, lambda p, a, f, x, m, b, d, e : With282(p, a, f, x, m, b, d, e))
    rubi.add(rule282)

    pattern283 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, a, d: NonzeroQ(S(27)*a**S(2)*d + S(4)*b**S(3))), )
    def With283(p, a, f, x, m, b, d, e):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        return (-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(-sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(a + b*x + d*x**S(3))**p*Int((e + f*x)**m*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(-sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(sqrt(S(3))*ImaginaryI + S(1)) - S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x)
    rule283 = ReplacementRule(pattern283, lambda p, a, f, x, m, b, d, e : With283(p, a, f, x, m, b, d, e))
    rubi.add(rule283)

    pattern284 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda a, d, c: ZeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))))
    rule284 = ReplacementRule(pattern284, lambda p, a, x, d, c : -S(3)**(-S(3)*p)*d**(-S(2)*p)*Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p), x))
    rubi.add(rule284)

    pattern285 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda a, d, c: NonzeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))))
    rule285 = ReplacementRule(pattern285, lambda p, a, x, d, c : Int(ExpandToSum((a + c*x**S(2) + d*x**S(3))**p, x), x))
    rubi.add(rule285)

    pattern286 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda a, d, c: NonzeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))), CustomConstraint(lambda x, p, u: ProductQ(NonfreeFactors(u, x))))
    def With286(p, a, x, d, c):
        u = Factor(a + c*x**S(2) + d*x**S(3))
        return FreeFactors(u, x)**p*Int(DistributeDegree(NonfreeFactors(u, x), p), x)
    rule286 = ReplacementRule(pattern286, lambda p, a, x, d, c : With286(p, a, x, d, c))
    rubi.add(rule286)

    pattern287 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda a, d, c: NonzeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))), )
    def With287(p, a, x, d, c):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        return S(3)**(-S(3)*p)*d**(-S(2)*p)*Int((c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p, x)
    rule287 = ReplacementRule(pattern287, lambda p, a, x, d, c : With287(p, a, x, d, c))
    rubi.add(rule287)

    pattern288 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda a, d, c: ZeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))))
    rule288 = ReplacementRule(pattern288, lambda p, a, x, d, c : (c - S(3)*d*x)**(-p)*(S(2)*c + S(3)*d*x)**(-S(2)*p)*(a + c*x**S(2) + d*x**S(3))**p*Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p), x))
    rubi.add(rule288)

    pattern289 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda a, d, c: NonzeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))), CustomConstraint(lambda d, a, x, p, c, u: ProductQ(u)))
    def With289(p, a, x, d, c):
        u = NonfreeFactors(Factor(a + c*x**S(2) + d*x**S(3)), x)
        return (a + c*x**S(2) + d*x**S(3))**p*Int(DistributeDegree(u, p), x)/DistributeDegree(u, p)
    rule289 = ReplacementRule(pattern289, lambda p, a, x, d, c : With289(p, a, x, d, c))
    rubi.add(rule289)

    pattern290 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda a, d, c: NonzeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))), )
    def With290(p, a, x, d, c):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        return (a + c*x**S(2) + d*x**S(3))**p*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*Int((c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p, x)
    rule290 = ReplacementRule(pattern290, lambda p, a, x, d, c : With290(p, a, x, d, c))
    rubi.add(rule290)

    pattern291 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda a, d, c: ZeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))))
    rule291 = ReplacementRule(pattern291, lambda p, a, f, x, m, d, e, c : -S(3)**(-S(3)*p)*d**(-S(2)*p)*Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p)*(e + f*x)**m, x))
    rubi.add(rule291)

    pattern292 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda a, d, c: NonzeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))))
    rule292 = ReplacementRule(pattern292, lambda p, a, f, x, m, d, e, c : Int(ExpandIntegrand((e + f*x)**m*(a + c*x**S(2) + d*x**S(3))**p, x), x))
    rubi.add(rule292)

    pattern293 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda a, d, c: NonzeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))), CustomConstraint(lambda f, e, x, p, u, m: ProductQ(NonfreeFactors(u, x))))
    def With293(p, a, f, x, m, d, e, c):
        u = Factor(a + c*x**S(2) + d*x**S(3))
        return FreeFactors(u, x)**p*Int((e + f*x)**m*DistributeDegree(NonfreeFactors(u, x), p), x)
    rule293 = ReplacementRule(pattern293, lambda p, a, f, x, m, d, e, c : With293(p, a, f, x, m, d, e, c))
    rubi.add(rule293)

    pattern294 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda a, d, c: NonzeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))), )
    def With294(p, a, f, x, m, d, e, c):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        return S(3)**(-S(3)*p)*d**(-S(2)*p)*Int((e + f*x)**m*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p, x)
    rule294 = ReplacementRule(pattern294, lambda p, a, f, x, m, d, e, c : With294(p, a, f, x, m, d, e, c))
    rubi.add(rule294)

    pattern295 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda a, d, c: ZeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))))
    rule295 = ReplacementRule(pattern295, lambda p, a, f, x, m, d, e, c : (c - S(3)*d*x)**(-p)*(S(2)*c + S(3)*d*x)**(-S(2)*p)*(a + c*x**S(2) + d*x**S(3))**p*Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p)*(e + f*x)**m, x))
    rubi.add(rule295)

    pattern296 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda a, d, c: NonzeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))), CustomConstraint(lambda f, d, m, e, a, x, p, c, u: ProductQ(u)))
    def With296(p, a, f, x, m, d, e, c):
        u = NonfreeFactors(Factor(a + c*x**S(2) + d*x**S(3)), x)
        return (a + c*x**S(2) + d*x**S(3))**p*Int((e + f*x)**m*DistributeDegree(u, p), x)/DistributeDegree(u, p)
    rule296 = ReplacementRule(pattern296, lambda p, a, f, x, m, d, e, c : With296(p, a, f, x, m, d, e, c))
    rubi.add(rule296)

    pattern297 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda a, d, c: NonzeroQ(S(27)*a*d**S(2) + S(4)*c**S(3))), )
    def With297(p, a, f, x, m, d, e, c):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        return (a + c*x**S(2) + d*x**S(3))**p*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*Int((e + f*x)**m*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(sqrt(S(3))*ImaginaryI + S(1)) + S(2)**(S(1)/3)*r**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p, x)
    rule297 = ReplacementRule(pattern297, lambda p, a, f, x, m, d, e, c : With297(p, a, f, x, m, d, e, c))
    rubi.add(rule297)

    pattern298 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda b, d, c: ZeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: ZeroQ(-S(3)*a*c + b**S(2))))
    rule298 = ReplacementRule(pattern298, lambda p, a, x, b, d, c : S(3)**(-p)*b**(-p)*c**(-p)*Int((b + c*x)**(S(3)*p), x))
    rubi.add(rule298)

    pattern299 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda b, d, c: ZeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))))
    rule299 = ReplacementRule(pattern299, lambda p, a, x, b, d, c : S(3)**(-p)*b**(-p)*c**(-p)*Subst(Int((S(3)*a*b*c - b**S(3) + c**S(3)*x**S(3))**p, x), x, c/(S(3)*d) + x))
    rubi.add(rule299)

    pattern300 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: ZeroQ(-S(3)*a*c + b**S(2))), )
    def With300(p, a, x, b, d, c):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        return S(3)**(-p)*b**(-p)*c**(-p)*Int((b + x*(c - r))**p*(b + x*(c + r*(-sqrt(S(3))*ImaginaryI + S(1))/S(2)))**p*(b + x*(c + r*(sqrt(S(3))*ImaginaryI + S(1))/S(2)))**p, x)
    rule300 = ReplacementRule(pattern300, lambda p, a, x, b, d, c : With300(p, a, x, b, d, c))
    rubi.add(rule300)

    pattern301 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))))
    rule301 = ReplacementRule(pattern301, lambda p, a, x, b, d, c : Int(ExpandToSum((a + b*x + c*x**S(2) + d*x**S(3))**p, x), x))
    rubi.add(rule301)

    pattern302 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))), CustomConstraint(lambda x, p, u: ProductQ(NonfreeFactors(u, x))))
    def With302(p, a, x, b, d, c):
        u = Factor(a + b*x + c*x**S(2) + d*x**S(3))
        return FreeFactors(u, x)**p*Int(DistributeDegree(NonfreeFactors(u, x), p), x)
    rule302 = ReplacementRule(pattern302, lambda p, a, x, b, d, c : With302(p, a, x, b, d, c))
    rubi.add(rule302)

    pattern303 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))))
    rule303 = ReplacementRule(pattern303, lambda p, a, x, b, d, c : S(3)**(-S(3)*p)*d**(-S(2)*p)*Subst(Int((S(27)*a*d**S(2) - S(9)*b*c*d + S(2)*c**S(3) + S(27)*d**S(3)*x**S(3) - S(9)*d*x*(-S(3)*b*d + c**S(2)))**p, x), x, c/(S(3)*d) + x))
    rubi.add(rule303)

    pattern304 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, d, c: ZeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: ZeroQ(-S(3)*a*c + b**S(2))))
    rule304 = ReplacementRule(pattern304, lambda p, a, x, b, d, c : (b + c*x)**(-S(3)*p)*(a + b*x + c*x**S(2) + d*x**S(3))**p*Int((b + c*x)**(S(3)*p), x))
    rubi.add(rule304)

    pattern305 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, d, c: ZeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))), )
    def With305(p, a, x, b, d, c):
        r = Rt(-S(3)*a*b*c + b**S(3), S(3))
        return (b + c*x - r)**(-p)*(b + c*x + r*(-sqrt(S(3))*ImaginaryI + S(1))/S(2))**(-p)*(b + c*x + r*(sqrt(S(3))*ImaginaryI + S(1))/S(2))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p*Int((b + c*x - r)**p*(b + c*x + r*(-sqrt(S(3))*ImaginaryI + S(1))/S(2))**p*(b + c*x + r*(sqrt(S(3))*ImaginaryI + S(1))/S(2))**p, x)
    rule305 = ReplacementRule(pattern305, lambda p, a, x, b, d, c : With305(p, a, x, b, d, c))
    rubi.add(rule305)

    pattern306 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: ZeroQ(-S(3)*a*c + b**S(2))), )
    def With306(p, a, x, b, d, c):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        return (b + x*(c - r))**(-p)*(b + x*(c + r*(-sqrt(S(3))*ImaginaryI + S(1))/S(2)))**(-p)*(b + x*(c + r*(sqrt(S(3))*ImaginaryI + S(1))/S(2)))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p*Int((b + x*(c - r))**p*(b + x*(c + r*(-sqrt(S(3))*ImaginaryI + S(1))/S(2)))**p*(b + x*(c + r*(sqrt(S(3))*ImaginaryI + S(1))/S(2)))**p, x)
    rule306 = ReplacementRule(pattern306, lambda p, a, x, b, d, c : With306(p, a, x, b, d, c))
    rubi.add(rule306)

    pattern307 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))), CustomConstraint(lambda d, a, b, x, p, c, u: ProductQ(u)))
    def With307(p, a, x, b, d, c):
        u = NonfreeFactors(Factor(a + b*x + c*x**S(2) + d*x**S(3)), x)
        return (a + b*x + c*x**S(2) + d*x**S(3))**p*Int(DistributeDegree(u, p), x)/DistributeDegree(u, p)
    rule307 = ReplacementRule(pattern307, lambda p, a, x, b, d, c : With307(p, a, x, b, d, c))
    rubi.add(rule307)

    pattern308 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))), )
    def With308(p, a, x, b, d, c):
        r = Rt(-S(27)*a*d**S(2) + S(9)*b*c*d - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) - S(18)*a*b*c*d + S(4)*a*c**S(3) + S(4)*b**S(3)*d - b**S(2)*c**S(2)), S(3))
        return (c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)**(S(1)/3)*ImaginaryI*r**S(2)*(-ImaginaryI + sqrt(S(3))) - S(6)*b*d*(-sqrt(S(3))*ImaginaryI + S(1)) + S(2)*c**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(2)**(S(1)/3)*ImaginaryI*r**S(2)*(ImaginaryI + sqrt(S(3))) - S(6)*b*d*(sqrt(S(3))*ImaginaryI + S(1)) + S(2)*c**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p*Int((c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)**(S(1)/3)*ImaginaryI*r**S(2)*(-ImaginaryI + sqrt(S(3))) - S(6)*b*d*(-sqrt(S(3))*ImaginaryI + S(1)) + S(2)*c**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(2)**(S(1)/3)*ImaginaryI*r**S(2)*(ImaginaryI + sqrt(S(3))) - S(6)*b*d*(sqrt(S(3))*ImaginaryI + S(1)) + S(2)*c**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p, x)
    rule308 = ReplacementRule(pattern308, lambda p, a, x, b, d, c : With308(p, a, x, b, d, c))
    rubi.add(rule308)

    pattern309 = Pattern(Integral(u_**p_, x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: PolyQ(u, x, S(3))), CustomConstraint(lambda x, u: Not(CubicMatchQ(u, x))))
    rule309 = ReplacementRule(pattern309, lambda p, x, u : Int(ExpandToSum(u, x)**p, x))
    rubi.add(rule309)

    pattern310 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda b, d, c: ZeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: ZeroQ(-S(3)*a*c + b**S(2))))
    rule310 = ReplacementRule(pattern310, lambda p, a, f, x, m, b, d, e, c : S(3)**(-p)*b**(-p)*c**(-p)*Int((b + c*x)**(S(3)*p)*(e + f*x)**m, x))
    rubi.add(rule310)

    pattern311 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda b, d, c: ZeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))), )
    def With311(p, a, f, x, m, b, d, e, c):
        r = Rt(-S(3)*a*b*c + b**S(3), S(3))
        return S(3)**(-p)*b**(-p)*c**(-p)*Int((e + f*x)**m*(b + c*x - r)**p*(b + c*x + r*(-sqrt(S(3))*ImaginaryI + S(1))/S(2))**p*(b + c*x + r*(sqrt(S(3))*ImaginaryI + S(1))/S(2))**p, x)
    rule311 = ReplacementRule(pattern311, lambda p, a, f, x, m, b, d, e, c : With311(p, a, f, x, m, b, d, e, c))
    rubi.add(rule311)

    pattern312 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: IntegerQ(p)), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: ZeroQ(-S(3)*a*c + b**S(2))), )
    def With312(p, a, f, x, m, b, d, e, c):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        return S(3)**(-p)*b**(-p)*c**(-p)*Int((b + x*(c - r))**p*(b + x*(c + r*(-sqrt(S(3))*ImaginaryI + S(1))/S(2)))**p*(b + x*(c + r*(sqrt(S(3))*ImaginaryI + S(1))/S(2)))**p*(e + f*x)**m, x)
    rule312 = ReplacementRule(pattern312, lambda p, a, f, x, m, b, d, e, c : With312(p, a, f, x, m, b, d, e, c))
    rubi.add(rule312)

    pattern313 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: PositiveIntegerQ(p)), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))))
    rule313 = ReplacementRule(pattern313, lambda p, a, f, x, m, b, d, e, c : Int(ExpandIntegrand((e + f*x)**m*(a + b*x + c*x**S(2) + d*x**S(3))**p, x), x))
    rubi.add(rule313)

    pattern314 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))), CustomConstraint(lambda f, e, x, p, u, m: ProductQ(NonfreeFactors(u, x))))
    def With314(p, a, f, x, m, b, d, e, c):
        u = Factor(a + b*x + c*x**S(2) + d*x**S(3))
        return FreeFactors(u, x)**p*Int((e + f*x)**m*DistributeDegree(NonfreeFactors(u, x), p), x)
    rule314 = ReplacementRule(pattern314, lambda p, a, f, x, m, b, d, e, c : With314(p, a, f, x, m, b, d, e, c))
    rubi.add(rule314)

    pattern315 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p: NegativeIntegerQ(p)), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))))
    rule315 = ReplacementRule(pattern315, lambda p, a, f, x, m, b, d, e, c : S(3)**(-S(3)*p)*d**(-S(2)*p)*Subst(Int((S(27)*a*d**S(2) - S(9)*b*c*d + S(2)*c**S(3) + S(27)*d**S(3)*x**S(3) - S(9)*d*x*(-S(3)*b*d + c**S(2)))**p, x), x, c/(S(3)*d) + x))
    rubi.add(rule315)

    pattern316 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, d, c: ZeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: ZeroQ(-S(3)*a*c + b**S(2))))
    rule316 = ReplacementRule(pattern316, lambda p, a, f, x, m, b, d, e, c : (b + c*x)**(-S(3)*p)*(a + b*x + c*x**S(2) + d*x**S(3))**p*Int((b + c*x)**(S(3)*p)*(e + f*x)**m, x))
    rubi.add(rule316)

    pattern317 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, d, c: ZeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))), )
    def With317(p, a, f, x, m, b, d, e, c):
        r = Rt(-S(3)*a*b*c + b**S(3), S(3))
        return (b + c*x - r)**(-p)*(b + c*x + r*(-sqrt(S(3))*ImaginaryI + S(1))/S(2))**(-p)*(b + c*x + r*(sqrt(S(3))*ImaginaryI + S(1))/S(2))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p*Int((e + f*x)**m*(b + c*x - r)**p*(b + c*x + r*(-sqrt(S(3))*ImaginaryI + S(1))/S(2))**p*(b + c*x + r*(sqrt(S(3))*ImaginaryI + S(1))/S(2))**p, x)
    rule317 = ReplacementRule(pattern317, lambda p, a, f, x, m, b, d, e, c : With317(p, a, f, x, m, b, d, e, c))
    rubi.add(rule317)

    pattern318 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: ZeroQ(-S(3)*a*c + b**S(2))), )
    def With318(p, a, f, x, m, b, d, e, c):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        return (b + x*(c - r))**(-p)*(b + x*(c + r*(-sqrt(S(3))*ImaginaryI + S(1))/S(2)))**(-p)*(b + x*(c + r*(sqrt(S(3))*ImaginaryI + S(1))/S(2)))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p*Int((b + x*(c - r))**p*(b + x*(c + r*(-sqrt(S(3))*ImaginaryI + S(1))/S(2)))**p*(b + x*(c + r*(sqrt(S(3))*ImaginaryI + S(1))/S(2)))**p*(e + f*x)**m, x)
    rule318 = ReplacementRule(pattern318, lambda p, a, f, x, m, b, d, e, c : With318(p, a, f, x, m, b, d, e, c))
    rubi.add(rule318)

    pattern319 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))), CustomConstraint(lambda f, d, b, m, e, a, x, p, c, u: ProductQ(u)))
    def With319(p, a, f, x, m, b, d, e, c):
        u = NonfreeFactors(Factor(a + b*x + c*x**S(2) + d*x**S(3)), x)
        return (a + b*x + c*x**S(2) + d*x**S(3))**p*Int((e + f*x)**m*DistributeDegree(u, p), x)/DistributeDegree(u, p)
    rule319 = ReplacementRule(pattern319, lambda p, a, f, x, m, b, d, e, c : With319(p, a, f, x, m, b, d, e, c))
    rubi.add(rule319)

    pattern320 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda p: Not(IntegerQ(p))), CustomConstraint(lambda b, d, c: NonzeroQ(-S(3)*b*d + c**S(2))), CustomConstraint(lambda b, a, c: NonzeroQ(-S(3)*a*c + b**S(2))), )
    def With320(p, a, f, x, m, b, d, e, c):
        r = Rt(-S(27)*a*d**S(2) + S(9)*b*c*d - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) - S(18)*a*b*c*d + S(4)*a*c**S(3) + S(4)*b**S(3)*d - b**S(2)*c**S(2)), S(3))
        return (c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)**(S(1)/3)*ImaginaryI*r**S(2)*(-ImaginaryI + sqrt(S(3))) - S(6)*b*d*(-sqrt(S(3))*ImaginaryI + S(1)) + S(2)*c**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(2)**(S(1)/3)*ImaginaryI*r**S(2)*(ImaginaryI + sqrt(S(3))) - S(6)*b*d*(sqrt(S(3))*ImaginaryI + S(1)) + S(2)*c**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p*Int((e + f*x)**m*(c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)**(S(1)/3)*ImaginaryI*r**S(2)*(-ImaginaryI + sqrt(S(3))) - S(6)*b*d*(-sqrt(S(3))*ImaginaryI + S(1)) + S(2)*c**S(2)*(-sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(2)**(S(1)/3)*ImaginaryI*r**S(2)*(ImaginaryI + sqrt(S(3))) - S(6)*b*d*(sqrt(S(3))*ImaginaryI + S(1)) + S(2)*c**S(2)*(sqrt(S(3))*ImaginaryI + S(1)))/(S(4)*r))**p, x)
    rule320 = ReplacementRule(pattern320, lambda p, a, f, x, m, b, d, e, c : With320(p, a, f, x, m, b, d, e, c))
    rubi.add(rule320)

    pattern321 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('p', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda v, x: PolyQ(v, x, S(3))), CustomConstraint(lambda v, x, u: Not(CubicMatchQ(v, x) & LinearMatchQ(u, x))))
    rule321 = ReplacementRule(pattern321, lambda p, v, x, u, m : Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**p, x))
    rubi.add(rule321)

    pattern322 = Pattern(Integral((f_ + x_**S(2)*WC('g', S(1)))/((d_ + x_**S(2)*WC('d', S(1)) + x_*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('a', S(1)) + x_**S(3)*WC('b', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, e, a, d: ZeroQ(-a*e + b*d)), CustomConstraint(lambda f, g: ZeroQ(f + g)), CustomConstraint(lambda a, c: PosQ(a**S(2)*(S(2)*a - c))))
    rule322 = ReplacementRule(pattern322, lambda e, a, g, f, x, b, d, c : a*f*atan((a*b*x**S(2) + a*b + x*(S(4)*a**S(2) - S(2)*a*c + b**S(2)))/(S(2)*sqrt(a*x**S(4) + a + b*x**S(3) + b*x + c*x**S(2))*Rt(a**S(2)*(S(2)*a - c), S(2))))/(d*Rt(a**S(2)*(S(2)*a - c), S(2))))
    rubi.add(rule322)

    pattern323 = Pattern(Integral((f_ + x_**S(2)*WC('g', S(1)))/((d_ + x_**S(2)*WC('d', S(1)) + x_*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('a', S(1)) + x_**S(3)*WC('b', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda b, e, a, d: ZeroQ(-a*e + b*d)), CustomConstraint(lambda f, g: ZeroQ(f + g)), CustomConstraint(lambda a, c: NegQ(a**S(2)*(S(2)*a - c))))
    rule323 = ReplacementRule(pattern323, lambda e, a, g, f, x, b, d, c : -a*f*atanh((a*b*x**S(2) + a*b + x*(S(4)*a**S(2) - S(2)*a*c + b**S(2)))/(S(2)*sqrt(a*x**S(4) + a + b*x**S(3) + b*x + c*x**S(2))*Rt(-a**S(2)*(S(2)*a - c), S(2))))/(d*Rt(-a**S(2)*(S(2)*a - c), S(2))))
    rubi.add(rule323)

    pattern324 = Pattern(Integral((x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda b, e, d, c: ZeroQ(S(8)*b*e**S(2) - S(4)*c*d*e + d**S(3))), CustomConstraint(lambda p: UnsameQ(p, S(2))), CustomConstraint(lambda p: UnsameQ(p, S(3))))
    rule324 = ReplacementRule(pattern324, lambda p, a, x, b, d, e, c : Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p, x), x), x, d/(S(4)*e) + x))
    rubi.add(rule324)

    pattern325 = Pattern(Integral(v_**p_, x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda v, x: PolynomialQ(v, x)), CustomConstraint(lambda v, x: Equal(Exponent(v, x), S(4))), CustomConstraint(lambda p: UnsameQ(p, S(2))), CustomConstraint(lambda p: UnsameQ(p, S(3))), CustomConstraint(lambda e, a, d, b, x, p, c: NonzeroQ(d) & ZeroQ(S(8)*b*e**S(2) - S(4)*c*d*e + d**S(3))))
    def With325(v, p, x):
        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        return Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p, x), x), x, d/(S(4)*e) + x)
    rule325 = ReplacementRule(pattern325, lambda v, p, x : With325(v, p, x))
    rubi.add(rule325)

    pattern326 = Pattern(Integral(u_*(x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: PolynomialQ(u, x)), CustomConstraint(lambda b, e, d, c: ZeroQ(S(8)*b*e**S(2) - S(4)*c*d*e + d**S(3))), CustomConstraint(lambda p: Not(PositiveIntegerQ(p))))
    rule326 = ReplacementRule(pattern326, lambda p, a, x, u, b, d, e, c : Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p*ReplaceAll(u, Rule(x, -d/(S(4)*e) + x)), x), x), x, d/(S(4)*e) + x))
    rubi.add(rule326)

    pattern327 = Pattern(Integral(u_*v_**p_, x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda x, u: PolynomialQ(u, x)), CustomConstraint(lambda v, x: PolynomialQ(v, x)), CustomConstraint(lambda v, x: Equal(Exponent(v, x), S(4))), CustomConstraint(lambda p: Not(PositiveIntegerQ(p))), CustomConstraint(lambda d, b, e, a, x, p, c, u: NonzeroQ(d) & ZeroQ(S(8)*b*e**S(2) - S(4)*c*d*e + d**S(3))))
    def With327(v, p, x, u):
        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        return Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p*ReplaceAll(u, Rule(x, -d/(S(4)*e) + x)), x), x), x, d/(S(4)*e) + x)
    rule327 = ReplacementRule(pattern327, lambda v, p, x, u : With327(v, p, x, u))
    rubi.add(rule327)

    pattern328 = Pattern(Integral((a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda b, a, d, c: ZeroQ(S(8)*a**S(2)*d - S(4)*a*b*c + b**S(3))), CustomConstraint(lambda p: IntegerQ(S(2)*p)))
    rule328 = ReplacementRule(pattern328, lambda p, a, x, b, d, e, c : -S(16)*a**S(2)*Subst(Int((a*(S(256)*a**S(4)*x**S(4) + S(256)*a**S(3)*e - S(64)*a**S(2)*b*d - S(32)*a**S(2)*x**S(2)*(-S(8)*a*c + S(3)*b**S(2)) + S(16)*a*b**S(2)*c - S(3)*b**S(4))/(-S(4)*a*x + b)**S(4))**p/(-S(4)*a*x + b)**S(2), x), x, 1/x + b/(S(4)*a)))
    rubi.add(rule328)

    pattern329 = Pattern(Integral(v_**p_, x_), CustomConstraint(lambda p, x: FreeQ(p, x)), CustomConstraint(lambda v, x: PolynomialQ(v, x)), CustomConstraint(lambda v, x: Equal(Exponent(v, x), S(4))), CustomConstraint(lambda p: IntegerQ(S(2)*p)), CustomConstraint(lambda d, e, a, b, x, p, c: NonzeroQ(a) & NonzeroQ(b) & ZeroQ(S(8)*a**S(2)*d - S(4)*a*b*c + b**S(3))))
    def With329(v, p, x):
        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        return -S(16)*a**S(2)*Subst(Int((a*(S(256)*a**S(4)*x**S(4) + S(256)*a**S(3)*e - S(64)*a**S(2)*b*d - S(32)*a**S(2)*x**S(2)*(-S(8)*a*c + S(3)*b**S(2)) + S(16)*a*b**S(2)*c - S(3)*b**S(4))/(-S(4)*a*x + b)**S(4))**p/(-S(4)*a*x + b)**S(2), x), x, 1/x + b/(S(4)*a))
    rule329 = ReplacementRule(pattern329, lambda v, p, x : With329(v, p, x))
    rubi.add(rule329)

    pattern330 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda D, x: FreeQ(D, x)), CustomConstraint(lambda b, d: ZeroQ(-b + d)), CustomConstraint(lambda a, e: ZeroQ(-a + e)), CustomConstraint(lambda b, a, x, c: SumQ(Factor(a*x**S(4) + a + b*x**S(3) + b*x + c*x**S(2)))), )
    def With330(e, a, C, x, A, B, D, d, b, c):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        return -Int((A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x)/q + Int((A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x)/q
    rule330 = ReplacementRule(pattern330, lambda e, a, C, x, A, B, D, d, b, c : With330(e, a, C, x, A, B, D, d, b, c))
    rubi.add(rule330)

    pattern331 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda D, x: FreeQ(D, x)), CustomConstraint(lambda b, d: ZeroQ(-b + d)), CustomConstraint(lambda a, e: ZeroQ(-a + e)), CustomConstraint(lambda b, a, x, c: SumQ(Factor(a*x**S(4) + a + b*x**S(3) + b*x + c*x**S(2)))), )
    def With331(e, a, x, A, B, D, d, b, c):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        return -Int((A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x)/q + Int((A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x)/q
    rule331 = ReplacementRule(pattern331, lambda e, a, x, A, B, D, d, b, c : With331(e, a, x, A, B, D, d, b, c))
    rubi.add(rule331)

    pattern332 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(3)*WC('D', S(1)) + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda D, x: FreeQ(D, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, d: ZeroQ(-b + d)), CustomConstraint(lambda a, e: ZeroQ(-a + e)), CustomConstraint(lambda b, a, x, c: SumQ(Factor(a*x**S(4) + a + b*x**S(3) + b*x + c*x**S(2)))), )
    def With332(e, a, C, m, x, A, B, D, d, b, c):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        return -Int(x**m*(A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x)/q + Int(x**m*(A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x)/q
    rule332 = ReplacementRule(pattern332, lambda e, a, C, m, x, A, B, D, d, b, c : With332(e, a, C, m, x, A, B, D, d, b, c))
    rubi.add(rule332)

    pattern333 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(3)*WC('D', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda D, x: FreeQ(D, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, d: ZeroQ(-b + d)), CustomConstraint(lambda a, e: ZeroQ(-a + e)), CustomConstraint(lambda b, a, x, c: SumQ(Factor(a*x**S(4) + a + b*x**S(3) + b*x + c*x**S(2)))), )
    def With333(e, a, m, x, A, B, D, d, b, c):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        return -Int(x**m*(A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x)/q + Int(x**m*(A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x)/q
    rule333 = ReplacementRule(pattern333, lambda e, a, m, x, A, B, D, d, b, c : With333(e, a, m, x, A, B, D, d, b, c))
    rubi.add(rule333)

    pattern334 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda e, C, A, B, b, d, c: ZeroQ(B**S(2)*d - S(2)*B*(S(2)*A*e + C*c) + S(2)*C*(A*d + C*b))), CustomConstraint(lambda e, a, C, A, B, d, c: ZeroQ(-S(4)*A*B*C*d + S(4)*A*e*(S(2)*A*C + B**S(2)) - B**S(3)*d + S(2)*B**S(2)*C*c - S(8)*C**S(3)*a)), CustomConstraint(lambda e, C, A, B, d, c: PosQ(C*(C*(-S(4)*c*e + d**S(2)) + S(2)*e*(-S(4)*A*e + B*d)))), )
    def With334(e, a, C, x, A, B, b, d, c):
        q = Rt(C*(C*(-S(4)*c*e + d**S(2)) + S(2)*e*(-S(4)*A*e + B*d)), S(2))
        return -S(2)*C**S(2)*atanh((-B*e + C*d + S(2)*C*e*x)/q)/q + S(2)*C**S(2)*atanh(C*(S(12)*A*B*e - S(4)*A*C*d - S(3)*B**S(2)*d + S(4)*B*C*c + S(8)*C**S(2)*e*x**S(3) + S(4)*C*x**S(2)*(-B*e + S(2)*C*d) + S(4)*C*x*(S(2)*A*e - B*d + S(2)*C*c))/(q*(-S(4)*A*C + B**S(2))))/q
    rule334 = ReplacementRule(pattern334, lambda e, a, C, x, A, B, b, d, c : With334(e, a, C, x, A, B, b, d, c))
    rubi.add(rule334)

    pattern335 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, d, A, C: ZeroQ(A*d + C*b)), CustomConstraint(lambda a, e, A, C: ZeroQ(-A**S(2)*e + C**S(2)*a)), CustomConstraint(lambda e, C, A, d, c: PosQ(C*(-S(8)*A*e**S(2) + C*(-S(4)*c*e + d**S(2))))), )
    def With335(e, a, C, x, A, b, d, c):
        q = Rt(C*(-S(8)*A*e**S(2) + C*(-S(4)*c*e + d**S(2))), S(2))
        return -S(2)*C**S(2)*atanh(C*(d + S(2)*e*x)/q)/q + S(2)*C**S(2)*atanh(C*(A*d - S(2)*C*d*x**S(2) - S(2)*C*e*x**S(3) - S(2)*x*(A*e + C*c))/(A*q))/q
    rule335 = ReplacementRule(pattern335, lambda e, a, C, x, A, b, d, c : With335(e, a, C, x, A, b, d, c))
    rubi.add(rule335)

    pattern336 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda e, C, A, B, b, d, c: ZeroQ(B**S(2)*d - S(2)*B*(S(2)*A*e + C*c) + S(2)*C*(A*d + C*b))), CustomConstraint(lambda e, a, C, A, B, d, c: ZeroQ(-S(4)*A*B*C*d + S(4)*A*e*(S(2)*A*C + B**S(2)) - B**S(3)*d + S(2)*B**S(2)*C*c - S(8)*C**S(3)*a)), CustomConstraint(lambda e, C, A, B, d, c: NegQ(C*(C*(-S(4)*c*e + d**S(2)) + S(2)*e*(-S(4)*A*e + B*d)))), )
    def With336(e, a, C, x, A, B, b, d, c):
        q = Rt(-C*(C*(-S(4)*c*e + d**S(2)) + S(2)*e*(-S(4)*A*e + B*d)), S(2))
        return S(2)*C**S(2)*atan((-B*e + C*d + S(2)*C*e*x)/q)/q - S(2)*C**S(2)*atan(C*(S(12)*A*B*e - S(4)*A*C*d - S(3)*B**S(2)*d + S(4)*B*C*c + S(8)*C**S(2)*e*x**S(3) + S(4)*C*x**S(2)*(-B*e + S(2)*C*d) + S(4)*C*x*(S(2)*A*e - B*d + S(2)*C*c))/(q*(-S(4)*A*C + B**S(2))))/q
    rule336 = ReplacementRule(pattern336, lambda e, a, C, x, A, B, b, d, c : With336(e, a, C, x, A, B, b, d, c))
    rubi.add(rule336)

    pattern337 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda b, d, A, C: ZeroQ(A*d + C*b)), CustomConstraint(lambda a, e, A, C: ZeroQ(-A**S(2)*e + C**S(2)*a)), CustomConstraint(lambda e, C, A, d, c: NegQ(C*(-S(8)*A*e**S(2) + C*(-S(4)*c*e + d**S(2))))), )
    def With337(e, a, C, x, A, b, d, c):
        q = Rt(-C*(-S(8)*A*e**S(2) + C*(-S(4)*c*e + d**S(2))), S(2))
        return S(2)*C**S(2)*atan((C*d + S(2)*C*e*x)/q)/q + S(2)*C**S(2)*atan(C*(-A*d + S(2)*C*d*x**S(2) + S(2)*C*e*x**S(3) + S(2)*x*(A*e + C*c))/(A*q))/q
    rule337 = ReplacementRule(pattern337, lambda e, a, C, x, A, b, d, c : With337(e, a, C, x, A, b, d, c))
    rubi.add(rule337)

    pattern338 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda C, x: FreeQ(C, x)), CustomConstraint(lambda D, x: FreeQ(D, x)), CustomConstraint(lambda b, C, A, B, D, d, e, c: ZeroQ(S(4)*d*(-S(2)*B*e + D*c)**S(2) - S(4)*(-S(2)*B*e + D*c)*(-S(8)*A*e**S(2) - S(4)*C*c*e + S(2)*D*b*e + S(3)*D*c*d) + S(8)*(-S(4)*C*e + S(3)*D*d)*(-A*d*e - C*b*e + D*b*d))), CustomConstraint(lambda b, a, C, A, B, D, d, e, c: ZeroQ(S(8)*a*(-S(4)*C*e + S(3)*D*d)**S(3) - S(8)*c*(-S(2)*B*e + D*c)**S(2)*(-S(4)*C*e + S(3)*D*d) + S(8)*d*(-S(4)*A*e + D*b)*(-S(2)*B*e + D*c)*(-S(4)*C*e + S(3)*D*d) + S(8)*d*(-S(2)*B*e + D*c)**S(3) - S(4)*e*(-S(4)*A*e + D*b)*(S(2)*(-S(4)*A*e + D*b)*(-S(4)*C*e + S(3)*D*d) + S(4)*(-S(2)*B*e + D*c)**S(2)))))
    rule338 = ReplacementRule(pattern338, lambda e, a, C, x, A, B, D, d, b, c : D*log(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4))/(S(4)*e) - Int((-S(4)*A*e + D*b + x**S(2)*(-S(4)*C*e + S(3)*D*d) + S(2)*x*(-S(2)*B*e + D*c))/(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4)), x)/(S(4)*e))
    rubi.add(rule338)

    pattern339 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda D, x: FreeQ(D, x)), CustomConstraint(lambda e, A, D, d, b, c: ZeroQ(D**S(2)*c**S(2)*d - D*c*(-S(8)*A*e**S(2) - S(4)*C*c*e + S(2)*D*b*e + S(3)*D*c*d) + S(2)*(-S(4)*C*e + S(3)*D*d)*(-A*d*e - C*b*e + D*b*d))), CustomConstraint(lambda e, b, a, A, B, D, d, c: ZeroQ(S(54)*D**S(3)*a*d**S(3) - S(6)*D*c*d*(-S(2)*B*e + D*c)**S(2) + S(6)*D*d**S(2)*(-S(4)*A*e + D*b)*(-S(2)*B*e + D*c) + S(2)*d*(-S(2)*B*e + D*c)**S(3) - e*(-S(4)*A*e + D*b)*(S(6)*D*d*(-S(4)*A*e + D*b) + S(4)*(-S(2)*B*e + D*c)**S(2)))))
    rule339 = ReplacementRule(pattern339, lambda e, a, x, A, B, D, d, b, c : D*log(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4))/(S(4)*e) - Int((-S(4)*A*e + D*b + S(3)*D*d*x**S(2) + S(2)*x*(-S(2)*B*e + D*c))/(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4)), x)/(S(4)*e))
    rubi.add(rule339)

    pattern340 = Pattern(Integral(u_/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a, d, c: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda f, a, e, c: ZeroQ(a*e**S(2) - c*f**S(2))))
    rule340 = ReplacementRule(pattern340, lambda a, f, x, u, b, d, e, c : -a*Int(u*sqrt(c + d*x)/x, x)/(f*(-a*d + b*c)) + c*Int(u*sqrt(a + b*x)/x, x)/(e*(-a*d + b*c)))
    rubi.add(rule340)

    pattern341 = Pattern(Integral(u_/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a, d, c: NonzeroQ(-a*d + b*c)), CustomConstraint(lambda b, d, f, e: ZeroQ(b*e**S(2) - d*f**S(2))))
    rule341 = ReplacementRule(pattern341, lambda a, f, x, u, b, d, e, c : b*Int(u*sqrt(c + d*x), x)/(f*(-a*d + b*c)) - d*Int(u*sqrt(a + b*x), x)/(e*(-a*d + b*c)))
    rubi.add(rule341)

    pattern342 = Pattern(Integral(u_/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda f, a, e, c: NonzeroQ(a*e**S(2) - c*f**S(2))), CustomConstraint(lambda b, d, f, e: NonzeroQ(b*e**S(2) - d*f**S(2))))
    rule342 = ReplacementRule(pattern342, lambda a, f, x, u, b, d, e, c : e*Int(u*sqrt(a + b*x)/(a*e**S(2) - c*f**S(2) + x*(b*e**S(2) - d*f**S(2))), x) - f*Int(u*sqrt(c + d*x)/(a*e**S(2) - c*f**S(2) + x*(b*e**S(2) - d*f**S(2))), x))
    rubi.add(rule342)

    pattern343 = Pattern(Integral(WC('u', S(1))/(x_**WC('n', S(1))*WC('d', S(1)) + sqrt(x_**WC('p', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, p: ZeroQ(-S(2)*n + p)), CustomConstraint(lambda b, d, c: ZeroQ(b*c**S(2) - d**S(2))))
    rule343 = ReplacementRule(pattern343, lambda p, n, a, x, u, b, d, c : -b*Int(u*x**n, x)/(a*d) + Int(u*sqrt(a + b*x**(S(2)*n)), x)/(a*c))
    rubi.add(rule343)

    pattern344 = Pattern(Integral(x_**WC('m', S(1))/(x_**WC('n', S(1))*WC('d', S(1)) + sqrt(x_**WC('p', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, p: ZeroQ(-S(2)*n + p)), CustomConstraint(lambda b, d, c: NonzeroQ(b*c**S(2) - d**S(2))))
    rule344 = ReplacementRule(pattern344, lambda p, n, a, x, m, b, d, c : c*Int(x**m*sqrt(a + b*x**(S(2)*n))/(a*c**S(2) + x**(S(2)*n)*(b*c**S(2) - d**S(2))), x) - d*Int(x**(m + n)/(a*c**S(2) + x**(S(2)*n)*(b*c**S(2) - d**S(2))), x))
    rubi.add(rule344)

    pattern345 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: PosQ(a/b)), )
    def With345(a, f, x, b, d, e):
        r = Numerator(Rt(a/b, S(3)))
        s = Denominator(Rt(a/b, S(3)))
        return r*Int(S(1)/((r + s*x)*sqrt(d + e*x + f*x**S(2))), x)/(S(3)*a) + r*Int((S(2)*r - s*x)/(sqrt(d + e*x + f*x**S(2))*(r**S(2) - r*s*x + s**S(2)*x**S(2))), x)/(S(3)*a)
    rule345 = ReplacementRule(pattern345, lambda a, f, x, b, d, e : With345(a, f, x, b, d, e))
    rubi.add(rule345)

    pattern346 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: PosQ(a/b)), )
    def With346(a, f, x, b, d):
        r = Numerator(Rt(a/b, S(3)))
        s = Denominator(Rt(a/b, S(3)))
        return r*Int(S(1)/(sqrt(d + f*x**S(2))*(r + s*x)), x)/(S(3)*a) + r*Int((S(2)*r - s*x)/(sqrt(d + f*x**S(2))*(r**S(2) - r*s*x + s**S(2)*x**S(2))), x)/(S(3)*a)
    rule346 = ReplacementRule(pattern346, lambda a, f, x, b, d : With346(a, f, x, b, d))
    rubi.add(rule346)

    pattern347 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NegQ(a/b)), )
    def With347(a, f, x, b, d, e):
        r = Numerator(Rt(-a/b, S(3)))
        s = Denominator(Rt(-a/b, S(3)))
        return r*Int(S(1)/((r - s*x)*sqrt(d + e*x + f*x**S(2))), x)/(S(3)*a) + r*Int((S(2)*r + s*x)/(sqrt(d + e*x + f*x**S(2))*(r**S(2) + r*s*x + s**S(2)*x**S(2))), x)/(S(3)*a)
    rule347 = ReplacementRule(pattern347, lambda a, f, x, b, d, e : With347(a, f, x, b, d, e))
    rubi.add(rule347)

    pattern348 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda b, a: NegQ(a/b)), )
    def With348(a, f, x, b, d):
        r = Numerator(Rt(-a/b, S(3)))
        s = Denominator(Rt(-a/b, S(3)))
        return r*Int(S(1)/(sqrt(d + f*x**S(2))*(r - s*x)), x)/(S(3)*a) + r*Int((S(2)*r + s*x)/(sqrt(d + f*x**S(2))*(r**S(2) + r*s*x + s**S(2)*x**S(2))), x)/(S(3)*a)
    rule348 = ReplacementRule(pattern348, lambda a, f, x, b, d : With348(a, f, x, b, d))
    rubi.add(rule348)

    pattern349 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule349 = ReplacementRule(pattern349, lambda e, a, x, b, d, c : d*Int(S(1)/((d**S(2) - e**S(2)*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x) - e*Int(x/((d**S(2) - e**S(2)*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x))
    rubi.add(rule349)

    pattern350 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)))
    rule350 = ReplacementRule(pattern350, lambda e, a, x, d, c : d*Int(S(1)/(sqrt(a + c*x**S(4))*(d**S(2) - e**S(2)*x**S(2))), x) - e*Int(x/(sqrt(a + c*x**S(4))*(d**S(2) - e**S(2)*x**S(2))), x))
    rubi.add(rule350)

    pattern351 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))**S(2)*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda a, b, d, e, c: NonzeroQ(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4))), CustomConstraint(lambda b, e, d, c: ZeroQ(b*d*e**S(2) + S(2)*c*d**S(3))))
    rule351 = ReplacementRule(pattern351, lambda e, a, x, b, d, c : -c*Int((d**S(2) - e**S(2)*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x)/(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4)) - e**S(3)*sqrt(a + b*x**S(2) + c*x**S(4))/((d + e*x)*(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4))))
    rubi.add(rule351)

    pattern352 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))**S(2)*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda a, b, d, e, c: NonzeroQ(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4))), CustomConstraint(lambda b, e, d, c: NonzeroQ(b*d*e**S(2) + S(2)*c*d**S(3))))
    rule352 = ReplacementRule(pattern352, lambda e, a, x, b, d, c : -c*Int((d**S(2) - e**S(2)*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x)/(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4)) - e**S(3)*sqrt(a + b*x**S(2) + c*x**S(4))/((d + e*x)*(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4))) + (b*d*e**S(2) + S(2)*c*d**S(3))*Int(S(1)/((d + e*x)*sqrt(a + b*x**S(2) + c*x**S(4))), x)/(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4)))
    rubi.add(rule352)

    pattern353 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))**S(2)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, a, d, c: NonzeroQ(a*e**S(4) + c*d**S(4))))
    rule353 = ReplacementRule(pattern353, lambda e, a, x, d, c : S(2)*c*d**S(3)*Int(S(1)/(sqrt(a + c*x**S(4))*(d + e*x)), x)/(a*e**S(4) + c*d**S(4)) - c*Int((d**S(2) - e**S(2)*x**S(2))/sqrt(a + c*x**S(4)), x)/(a*e**S(4) + c*d**S(4)) - e**S(3)*sqrt(a + c*x**S(4))/((d + e*x)*(a*e**S(4) + c*d**S(4))))
    rubi.add(rule353)

    pattern354 = Pattern(Integral((A_ + x_**S(2)*WC('B', S(1)))/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda B, e, d, A: ZeroQ(A*e + B*d)), CustomConstraint(lambda e, a, d, c: ZeroQ(-a*e**S(2) + c*d**S(2))))
    rule354 = ReplacementRule(pattern354, lambda a, x, A, B, b, d, e, c : A*Subst(Int(1/(d - x**S(2)*(-S(2)*a*e + b*d)), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4))))
    rubi.add(rule354)

    pattern355 = Pattern(Integral((A_ + x_**S(2)*WC('B', S(1)))/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda B, e, d, A: ZeroQ(A*e + B*d)), CustomConstraint(lambda e, a, d, c: ZeroQ(-a*e**S(2) + c*d**S(2))))
    rule355 = ReplacementRule(pattern355, lambda a, x, A, B, d, e, c : A*Subst(Int(1/(S(2)*a*e*x**S(2) + d), x), x, x/sqrt(a + c*x**S(4))))
    rubi.add(rule355)

    pattern356 = Pattern(Integral((A_ + x_**S(4)*WC('B', S(1)))/(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))*(d_ + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda f, a, d, c: ZeroQ(-a*f + c*d)), CustomConstraint(lambda B, a, c, A: ZeroQ(A*c + B*a)))
    rule356 = ReplacementRule(pattern356, lambda a, f, x, A, B, b, d, e, c : A*Subst(Int(1/(d - x**S(2)*(-a*e + b*d)), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4))))
    rubi.add(rule356)

    pattern357 = Pattern(Integral((A_ + x_**S(4)*WC('B', S(1)))/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda f, a, d, c: ZeroQ(-a*f + c*d)), CustomConstraint(lambda B, a, c, A: ZeroQ(A*c + B*a)))
    rule357 = ReplacementRule(pattern357, lambda a, f, x, A, B, d, e, c : A*Subst(Int(1/(a*e*x**S(2) + d), x), x, x/sqrt(a + c*x**S(4))))
    rubi.add(rule357)

    pattern358 = Pattern(Integral((A_ + x_**S(4)*WC('B', S(1)))/((d_ + x_**S(4)*WC('f', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda f, a, d, c: ZeroQ(-a*f + c*d)), CustomConstraint(lambda B, a, c, A: ZeroQ(A*c + B*a)))
    rule358 = ReplacementRule(pattern358, lambda a, f, x, A, B, b, d, c : A*Subst(Int(1/(-b*d*x**S(2) + d), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4))))
    rubi.add(rule358)

    pattern359 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/(d_ + x_**S(4)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, a, d, c: ZeroQ(a*e + c*d)), CustomConstraint(lambda a, c: PosQ(a*c)))
    rule359 = ReplacementRule(pattern359, lambda a, x, b, d, e, c : a*Subst(Int(1/(-S(2)*b*x**S(2) + x**S(4)*(-S(4)*a*c + b**S(2)) + S(1)), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4)))/d)
    rubi.add(rule359)

    pattern360 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/(d_ + x_**S(4)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda e, a, d, c: ZeroQ(a*e + c*d)), CustomConstraint(lambda a, c: NegQ(a*c)), )
    def With360(a, x, b, d, e, c):
        q = sqrt(-S(4)*a*c + b**S(2))
        return sqrt(S(2))*a*sqrt(-b + q)*atanh(sqrt(S(2))*x*sqrt(-b + q)*(b + S(2)*c*x**S(2) + q)/(S(4)*sqrt(a + b*x**S(2) + c*x**S(4))*Rt(-a*c, S(2))))/(S(4)*d*Rt(-a*c, S(2))) - sqrt(S(2))*a*sqrt(b + q)*atan(sqrt(S(2))*x*sqrt(b + q)*(b + S(2)*c*x**S(2) - q)/(S(4)*sqrt(a + b*x**S(2) + c*x**S(4))*Rt(-a*c, S(2))))/(S(4)*d*Rt(-a*c, S(2)))
    rule360 = ReplacementRule(pattern360, lambda a, x, b, d, e, c : With360(a, x, b, d, e, c))
    rubi.add(rule360)

    pattern361 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)))
    rule361 = ReplacementRule(pattern361, lambda e, a, f, x, b, d, c : a*Int(S(1)/((a**S(2) - b**S(2)*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x) - b*Int(x/((a**S(2) - b**S(2)*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x))
    rubi.add(rule361)

    pattern362 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*sqrt(x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda e, c, a, f, b, h, d, g: ZeroQ(-f**S(2)*(a*h**S(2) - b*g*h + c*g**S(2)) + (-d*h + e*g)**S(2))), CustomConstraint(lambda e, c, f, b, h, d, g: ZeroQ(-S(2)*d*e*h + S(2)*e**S(2)*g - f**S(2)*(-b*h + S(2)*c*g))))
    rule362 = ReplacementRule(pattern362, lambda c, a, f, x, b, d, h, e, g : S(2)*sqrt(d + e*x + f*sqrt(a + b*x + c*x**S(2)))*(S(9)*c**S(2)*f*g*h*x**S(2) + S(3)*c**S(2)*f*h**S(2)*x**S(3) + c*f*x*(a*h**S(2) - b*g*h + S(10)*c*g**S(2)) + f*(S(2)*a*b*h**S(2) - S(3)*a*c*g*h - S(2)*b**S(2)*g*h + S(5)*b*c*g**S(2)) - (-d*h + e*g)*sqrt(a + b*x + c*x**S(2))*(-S(2)*b*h + S(5)*c*g + c*h*x))/(S(15)*c**S(2)*f*(g + h*x)))
    rubi.add(rule362)

    pattern363 = Pattern(Integral((u_ + (sqrt(v_)*WC('k', S(1)) + WC('j', S(0)))*WC('f', S(1)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda k, x: FreeQ(k, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda f, v, x, u, j: Not(LinearMatchQ(u, x) & QuadraticMatchQ(v, x) & (ZeroQ(j) | ZeroQ(f + S(-1))))), CustomConstraint(lambda k, f, v, x, u, j, h, g: ZeroQ(-f**S(2)*k**S(2)*(g**S(2)*Coefficient(v, x, S(2)) - g*h*Coefficient(v, x, S(1)) + h**S(2)*Coefficient(v, x, S(0))) + (g*Coefficient(u, x, S(1)) - h*(f*j + Coefficient(u, x, S(0))))**S(2))))
    rule363 = ReplacementRule(pattern363, lambda n, k, m, f, v, x, u, j, h, g : Int((g + h*x)**m*(f*k*sqrt(ExpandToSum(v, x)) + ExpandToSum(f*j + u, x))**n, x))
    rubi.add(rule363)

    pattern364 = Pattern(Integral(((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**n_*WC('h', S(1)) + WC('g', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, e, c: ZeroQ(-c*f**S(2) + e**S(2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule364 = ReplacementRule(pattern364, lambda p, c, n, a, f, x, b, d, h, e, g : S(2)*Subst(Int((g + h*x**n)**p*(d**S(2)*e + e*x**S(2) - f**S(2)*(-a*e + b*d) - x*(-b*f**S(2) + S(2)*d*e))/(b*f**S(2) - S(2)*d*e + S(2)*e*x)**S(2), x), x, d + e*x + f*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule364)

    pattern365 = Pattern(Integral(((x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**n_*WC('h', S(1)) + WC('g', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, e, c: ZeroQ(-c*f**S(2) + e**S(2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule365 = ReplacementRule(pattern365, lambda p, c, n, a, f, x, d, h, e, g : Subst(Int((g + h*x**n)**p*(a*f**S(2) + d**S(2) - S(2)*d*x + x**S(2))/(d - x)**S(2), x), x, d + e*x + f*sqrt(a + c*x**S(2)))/(S(2)*e))
    rubi.add(rule365)

    pattern366 = Pattern(Integral(((u_ + sqrt(v_)*WC('f', S(1)))**n_*WC('h', S(1)) + WC('g', S(0)))**WC('p', S(1)), x_), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda v, x: QuadraticQ(v, x)), CustomConstraint(lambda v, x, u: Not(LinearMatchQ(u, x) & QuadraticMatchQ(v, x))), CustomConstraint(lambda f, v, x, u: ZeroQ(-f**S(2)*Coefficient(v, x, S(2)) + Coefficient(u, x, S(1))**S(2))), CustomConstraint(lambda p: IntegerQ(p)))
    rule366 = ReplacementRule(pattern366, lambda p, n, f, v, x, u, h, g : Int((g + h*(f*sqrt(ExpandToSum(v, x)) + ExpandToSum(u, x))**n)**p, x))
    rubi.add(rule366)

    pattern367 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*WC('f', S(1)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, e, c: ZeroQ(-c*f**S(2) + e**S(2))), CustomConstraint(lambda m: IntegerQ(m)))
    rule367 = ReplacementRule(pattern367, lambda c, n, a, f, x, m, h, e, g : S(2)**(-m + S(-1))*e**(-m + S(-1))*Subst(Int(x**(-m + n + S(-2))*(a*f**S(2) + x**S(2))*(-a*f**S(2)*h + S(2)*e*g*x + h*x**S(2))**m, x), x, e*x + f*sqrt(a + c*x**S(2))))
    rubi.add(rule367)

    pattern368 = Pattern(Integral(x_**WC('p', S(1))*(g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, e, c: ZeroQ(-c*f**S(2) + e**S(2))), CustomConstraint(lambda g, i, a, c: ZeroQ(-a*i + c*g)), CustomConstraint(lambda p, m: IntegersQ(p, S(2)*m)), CustomConstraint(lambda i, c, m: IntegerQ(m) | PositiveQ(i/c)))
    rule368 = ReplacementRule(pattern368, lambda i, p, n, a, g, f, x, m, e, c : S(2)**(-S(2)*m - p + S(-1))*e**(-p + S(-1))*f**(-S(2)*m)*(i/c)**m*Subst(Int(x**(-S(2)*m + n - p + S(-2))*(-a*f**S(2) + x**S(2))**p*(a*f**S(2) + x**S(2))**(S(2)*m + S(1)), x), x, e*x + f*sqrt(a + c*x**S(2))))
    rubi.add(rule368)

    pattern369 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1))*(x_**S(2)*WC('i', S(1)) + x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, e, c: ZeroQ(-c*f**S(2) + e**S(2))), CustomConstraint(lambda g, i, a, c: ZeroQ(-a*i + c*g)), CustomConstraint(lambda b, i, h, c: ZeroQ(-b*i + c*h)), CustomConstraint(lambda m: IntegerQ(S(2)*m)), CustomConstraint(lambda i, c, m: IntegerQ(m) | PositiveQ(i/c)))
    rule369 = ReplacementRule(pattern369, lambda i, e, c, n, a, f, x, m, b, h, d, g : S(2)*f**(-S(2)*m)*(i/c)**m*Subst(Int(x**n*(b*f**S(2) - S(2)*d*e + S(2)*e*x)**(-S(2)*m + S(-2))*(d**S(2)*e + e*x**S(2) - f**S(2)*(-a*e + b*d) - x*(-b*f**S(2) + S(2)*d*e))**(S(2)*m + S(1)), x), x, d + e*x + f*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule369)

    pattern370 = Pattern(Integral((g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, e, c: ZeroQ(-c*f**S(2) + e**S(2))), CustomConstraint(lambda g, i, a, c: ZeroQ(-a*i + c*g)), CustomConstraint(lambda m: IntegerQ(S(2)*m)), CustomConstraint(lambda i, c, m: IntegerQ(m) | PositiveQ(i/c)))
    rule370 = ReplacementRule(pattern370, lambda i, n, a, g, f, x, m, d, e, c : S(2)**(-S(2)*m + S(-1))*f**(-S(2)*m)*(i/c)**m*Subst(Int(x**n*(-d + x)**(-S(2)*m + S(-2))*(a*f**S(2) + d**S(2) - S(2)*d*x + x**S(2))**(S(2)*m + S(1)), x), x, d + e*x + f*sqrt(a + c*x**S(2)))/e)
    rubi.add(rule370)

    pattern371 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1))*(x_**S(2)*WC('i', S(1)) + x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, e, c: ZeroQ(-c*f**S(2) + e**S(2))), CustomConstraint(lambda g, i, a, c: ZeroQ(-a*i + c*g)), CustomConstraint(lambda b, i, h, c: ZeroQ(-b*i + c*h)), CustomConstraint(lambda m: PositiveIntegerQ(m + S(1)/2)), CustomConstraint(lambda i, c: Not(PositiveQ(i/c))))
    rule371 = ReplacementRule(pattern371, lambda i, e, c, n, a, f, x, m, b, h, d, g : (i/c)**(m + S(-1)/2)*sqrt(g + h*x + i*x**S(2))*Int((a + b*x + c*x**S(2))**m*(d + e*x + f*sqrt(a + b*x + c*x**S(2)))**n, x)/sqrt(a + b*x + c*x**S(2)))
    rubi.add(rule371)

    pattern372 = Pattern(Integral((g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, e, c: ZeroQ(-c*f**S(2) + e**S(2))), CustomConstraint(lambda g, i, a, c: ZeroQ(-a*i + c*g)), CustomConstraint(lambda m: PositiveIntegerQ(m + S(1)/2)), CustomConstraint(lambda i, c: Not(PositiveQ(i/c))))
    rule372 = ReplacementRule(pattern372, lambda i, n, a, g, f, x, m, d, e, c : (i/c)**(m + S(-1)/2)*sqrt(g + i*x**S(2))*Int((a + c*x**S(2))**m*(d + e*x + f*sqrt(a + c*x**S(2)))**n, x)/sqrt(a + c*x**S(2)))
    rubi.add(rule372)

    pattern373 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1))*(x_**S(2)*WC('i', S(1)) + x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda h, x: FreeQ(h, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, e, c: ZeroQ(-c*f**S(2) + e**S(2))), CustomConstraint(lambda g, i, a, c: ZeroQ(-a*i + c*g)), CustomConstraint(lambda b, i, h, c: ZeroQ(-b*i + c*h)), CustomConstraint(lambda m: NegativeIntegerQ(m + S(-1)/2)), CustomConstraint(lambda i, c: Not(PositiveQ(i/c))))
    rule373 = ReplacementRule(pattern373, lambda i, e, c, n, a, f, x, m, b, h, d, g : (i/c)**(m + S(1)/2)*sqrt(a + b*x + c*x**S(2))*Int((a + b*x + c*x**S(2))**m*(d + e*x + f*sqrt(a + b*x + c*x**S(2)))**n, x)/sqrt(g + h*x + i*x**S(2)))
    rubi.add(rule373)

    pattern374 = Pattern(Integral((g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda i, x: FreeQ(i, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda f, e, c: ZeroQ(-c*f**S(2) + e**S(2))), CustomConstraint(lambda g, i, a, c: ZeroQ(-a*i + c*g)), CustomConstraint(lambda m: NegativeIntegerQ(m + S(-1)/2)), CustomConstraint(lambda i, c: Not(PositiveQ(i/c))))
    rule374 = ReplacementRule(pattern374, lambda i, n, a, g, f, x, m, d, e, c : (i/c)**(m + S(1)/2)*sqrt(a + c*x**S(2))*Int((a + c*x**S(2))**m*(d + e*x + f*sqrt(a + c*x**S(2)))**n, x)/sqrt(g + i*x**S(2)))
    rubi.add(rule374)

    pattern375 = Pattern(Integral(w_**WC('m', S(1))*(u_ + (sqrt(v_)*WC('k', S(1)) + WC('j', S(0)))*WC('f', S(1)))**WC('n', S(1)), x_), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda j, x: FreeQ(j, x)), CustomConstraint(lambda k, x: FreeQ(k, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, u: LinearQ(u, x)), CustomConstraint(lambda w, v, x: QuadraticQ(List(v, w), x)), CustomConstraint(lambda w, v, f, x, u, j: Not(LinearMatchQ(u, x) & QuadraticMatchQ(List(v, w), x) & (ZeroQ(j) | ZeroQ(f + S(-1))))), CustomConstraint(lambda k, f, v, x, u: ZeroQ(-f**S(2)*k**S(2)*Coefficient(v, x, S(2)) + Coefficient(u, x, S(1))**S(2))))
    rule375 = ReplacementRule(pattern375, lambda j, n, k, w, f, v, x, u, m : Int((f*k*sqrt(ExpandToSum(v, x)) + ExpandToSum(f*j + u, x))**n*ExpandToSum(w, x)**m, x))
    rubi.add(rule375)

    pattern376 = Pattern(Integral(S(1)/((a_ + x_**WC('n', S(1))*WC('b', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + (a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n, p: ZeroQ(p - S(2)/n)))
    rule376 = ReplacementRule(pattern376, lambda p, n, a, x, b, d, c : Subst(Int(1/(-c*x**S(2) + S(1)), x), x, x/sqrt(c*x**S(2) + d*(a + b*x**n)**(S(2)/n)))/a)
    rubi.add(rule376)

    pattern377 = Pattern(Integral(sqrt(a_ + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, a, c: ZeroQ(a**S(2) - b**S(2)*c)))
    rule377 = ReplacementRule(pattern377, lambda a, x, b, d, c : S(2)*a*x/sqrt(a + b*sqrt(c + d*x**S(2))) + S(2)*b**S(2)*d*x**S(3)/(S(3)*(a + b*sqrt(c + d*x**S(2)))**(S(3)/2)))
    rubi.add(rule377)

    pattern378 = Pattern(Integral(sqrt(x_**S(2)*WC('a', S(1)) + x_*sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)))/(x_*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, a, d: ZeroQ(a**S(2) - b**S(2)*d)), CustomConstraint(lambda b, a, c: ZeroQ(a + b**S(2)*c)))
    rule378 = ReplacementRule(pattern378, lambda a, x, b, d, c : sqrt(S(2))*b*Subst(Int(1/sqrt(S(1) + x**S(2)/a), x), x, a*x + b*sqrt(c + d*x**S(2)))/a)
    rubi.add(rule378)

    pattern379 = Pattern(Integral(sqrt(x_*(x_*WC('a', S(1)) + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)))*WC('e', S(1)))/(x_*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda b, a, d: ZeroQ(a**S(2) - b**S(2)*d)), CustomConstraint(lambda b, a, e, c: ZeroQ(a + b**S(2)*c*e)))
    rule379 = ReplacementRule(pattern379, lambda a, x, b, d, e, c : Int(sqrt(a*e*x**S(2) + b*e*x*sqrt(c + d*x**S(2)))/(x*sqrt(c + d*x**S(2))), x))
    rubi.add(rule379)

    pattern380 = Pattern(Integral(sqrt(x_**S(2)*WC('c', S(1)) + sqrt(a_ + x_**S(4)*WC('b', S(1)))*WC('d', S(1)))/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda b, d, c: ZeroQ(-b*d**S(2) + c**S(2))))
    rule380 = ReplacementRule(pattern380, lambda a, x, b, d, c : d*Subst(Int(1/(-S(2)*c*x**S(2) + S(1)), x), x, x/sqrt(c*x**S(2) + d*sqrt(a + b*x**S(4)))))
    rubi.add(rule380)

    pattern381 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sqrt(x_**S(2)*WC('b', S(1)) + sqrt(a_ + x_**S(4)*WC('e', S(1))))/sqrt(a_ + x_**S(4)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda b, e: ZeroQ(-b**S(2) + e)), CustomConstraint(lambda a: PositiveQ(a)))
    rule381 = ReplacementRule(pattern381, lambda a, x, m, b, d, e, c : (-ImaginaryI/S(2) + S(1)/2)*Int((c + d*x)**m/sqrt(-ImaginaryI*b*x**S(2) + sqrt(a)), x) + (ImaginaryI/S(2) + S(1)/2)*Int((c + d*x)**m/sqrt(ImaginaryI*b*x**S(2) + sqrt(a)), x))
    rubi.add(rule381)

    pattern382 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda a: PosQ(a)), CustomConstraint(lambda b: PosQ(b)), )
    def With382(a, x, b, d, c):
        q = Rt(b/a, S(3))
        return d*Int((q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x)/(-c*q + d*(S(1) + sqrt(S(3)))) - q*Int(1/sqrt(a + b*x**S(3)), x)/(-c*q + d*(S(1) + sqrt(S(3))))
    rule382 = ReplacementRule(pattern382, lambda a, x, b, d, c : With382(a, x, b, d, c))
    rubi.add(rule382)

    pattern383 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda a: PosQ(a)), CustomConstraint(lambda b: NegQ(b)), )
    def With383(a, x, b, d, c):
        q = Rt(-b/a, S(3))
        return d*Int((-q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x)/(c*q + d*(S(1) + sqrt(S(3)))) + q*Int(1/sqrt(a + b*x**S(3)), x)/(c*q + d*(S(1) + sqrt(S(3))))
    rule383 = ReplacementRule(pattern383, lambda a, x, b, d, c : With383(a, x, b, d, c))
    rubi.add(rule383)

    pattern384 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda a: NegQ(a)), CustomConstraint(lambda b: PosQ(b)), )
    def With384(a, x, b, d, c):
        q = Rt(-b/a, S(3))
        return d*Int((-q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x)/(c*q + d*(-sqrt(S(3)) + S(1))) + q*Int(1/sqrt(a + b*x**S(3)), x)/(c*q + d*(-sqrt(S(3)) + S(1)))
    rule384 = ReplacementRule(pattern384, lambda a, x, b, d, c : With384(a, x, b, d, c))
    rubi.add(rule384)

    pattern385 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda a: NegQ(a)), CustomConstraint(lambda b: NegQ(b)), )
    def With385(a, x, b, d, c):
        q = Rt(b/a, S(3))
        return d*Int((q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x)/(-c*q + d*(-sqrt(S(3)) + S(1))) - q*Int(1/sqrt(a + b*x**S(3)), x)/(-c*q + d*(-sqrt(S(3)) + S(1)))
    rule385 = ReplacementRule(pattern385, lambda a, x, b, d, c : With385(a, x, b, d, c))
    rubi.add(rule385)

    pattern386 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda a: PosQ(a)), CustomConstraint(lambda b: PosQ(b)), CustomConstraint(lambda f, d, a, e, b, x, q, c: ZeroQ(-e*q + f*(S(1) + sqrt(S(3))))))
    def With386(e, a, f, x, b, d, c):
        q = Rt(b/a, S(3))
        return S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) - q*x + S(1))/(q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(-sqrt(S(3)) + S(2))*(q*x + S(1))*Subst(Int(S(1)/(sqrt(-x**S(2) + S(1))*sqrt(x**S(2) - S(4)*sqrt(S(3)) + S(7))*(-c*q + d*(-sqrt(S(3)) + S(1)) + x*(-c*q + d*(S(1) + sqrt(S(3)))))), x), x, (-q*x + S(-1) + sqrt(S(3)))/(q*x + S(1) + sqrt(S(3))))/(q*sqrt((q*x + S(1))/(q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(a + b*x**S(3)))
    rule386 = ReplacementRule(pattern386, lambda e, a, f, x, b, d, c : With386(e, a, f, x, b, d, c))
    rubi.add(rule386)

    pattern387 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda a: PosQ(a)), CustomConstraint(lambda b: PosQ(b)), CustomConstraint(lambda f, d, b, e, a, x, q, c: NonzeroQ(-e*q + f*(S(1) + sqrt(S(3))))))
    def With387(a, f, x, b, d, e, c):
        q = Rt(b/a, S(3))
        return (-c*f + d*e)*Int((q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x)/(-c*q + d*(S(1) + sqrt(S(3)))) + (-e*q + f*(S(1) + sqrt(S(3))))*Int(1/sqrt(a + b*x**S(3)), x)/(-c*q + d*(S(1) + sqrt(S(3))))
    rule387 = ReplacementRule(pattern387, lambda a, f, x, b, d, e, c : With387(a, f, x, b, d, e, c))
    rubi.add(rule387)

    pattern388 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda a: PosQ(a)), CustomConstraint(lambda b: NegQ(b)), CustomConstraint(lambda f, d, a, e, b, x, q, c: ZeroQ(e*q + f*(S(1) + sqrt(S(3))))))
    def With388(e, a, f, x, b, d, c):
        q = Rt(-b/a, S(3))
        return -S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) + q*x + S(1))/(-q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(-sqrt(S(3)) + S(2))*(-q*x + S(1))*Subst(Int(S(1)/(sqrt(-x**S(2) + S(1))*sqrt(x**S(2) - S(4)*sqrt(S(3)) + S(7))*(c*q + d*(-sqrt(S(3)) + S(1)) + x*(c*q + d*(S(1) + sqrt(S(3)))))), x), x, (q*x + S(-1) + sqrt(S(3)))/(-q*x + S(1) + sqrt(S(3))))/(q*sqrt((-q*x + S(1))/(-q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(a + b*x**S(3)))
    rule388 = ReplacementRule(pattern388, lambda e, a, f, x, b, d, c : With388(e, a, f, x, b, d, c))
    rubi.add(rule388)

    pattern389 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda a: PosQ(a)), CustomConstraint(lambda b: NegQ(b)), CustomConstraint(lambda f, d, b, e, a, x, q, c: NonzeroQ(e*q + f*(S(1) + sqrt(S(3))))))
    def With389(a, f, x, b, d, e, c):
        q = Rt(-b/a, S(3))
        return (-c*f + d*e)*Int((-q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x)/(c*q + d*(S(1) + sqrt(S(3)))) + (e*q + f*(S(1) + sqrt(S(3))))*Int(1/sqrt(a + b*x**S(3)), x)/(c*q + d*(S(1) + sqrt(S(3))))
    rule389 = ReplacementRule(pattern389, lambda a, f, x, b, d, e, c : With389(a, f, x, b, d, e, c))
    rubi.add(rule389)

    pattern390 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda a: NegQ(a)), CustomConstraint(lambda b: PosQ(b)), CustomConstraint(lambda f, d, a, e, b, x, q, c: ZeroQ(e*q + f*(-sqrt(S(3)) + S(1)))))
    def With390(e, a, f, x, b, d, c):
        q = Rt(-b/a, S(3))
        return S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) + q*x + S(1))/(-q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(sqrt(S(3)) + S(2))*(-q*x + S(1))*Subst(Int(S(1)/(sqrt(-x**S(2) + S(1))*sqrt(x**S(2) + S(4)*sqrt(S(3)) + S(7))*(c*q + d*(S(1) + sqrt(S(3))) + x*(c*q + d*(-sqrt(S(3)) + S(1))))), x), x, (-q*x + S(1) + sqrt(S(3)))/(q*x + S(-1) + sqrt(S(3))))/(q*sqrt((q*x + S(-1))/(-q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(a + b*x**S(3)))
    rule390 = ReplacementRule(pattern390, lambda e, a, f, x, b, d, c : With390(e, a, f, x, b, d, c))
    rubi.add(rule390)

    pattern391 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda a: NegQ(a)), CustomConstraint(lambda b: PosQ(b)), CustomConstraint(lambda f, d, b, e, a, x, q, c: NonzeroQ(e*q + f*(-sqrt(S(3)) + S(1)))))
    def With391(a, f, x, b, d, e, c):
        q = Rt(-b/a, S(3))
        return (-c*f + d*e)*Int((-q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x)/(c*q + d*(-sqrt(S(3)) + S(1))) + (e*q + f*(-sqrt(S(3)) + S(1)))*Int(1/sqrt(a + b*x**S(3)), x)/(c*q + d*(-sqrt(S(3)) + S(1)))
    rule391 = ReplacementRule(pattern391, lambda a, f, x, b, d, e, c : With391(a, f, x, b, d, e, c))
    rubi.add(rule391)

    pattern392 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda a: NegQ(a)), CustomConstraint(lambda b: NegQ(b)), CustomConstraint(lambda f, d, a, e, b, x, q, c: ZeroQ(-e*q + f*(-sqrt(S(3)) + S(1)))))
    def With392(e, a, f, x, b, d, c):
        q = Rt(b/a, S(3))
        return -S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) - q*x + S(1))/(q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(sqrt(S(3)) + S(2))*(q*x + S(1))*Subst(Int(S(1)/(sqrt(-x**S(2) + S(1))*sqrt(x**S(2) + S(4)*sqrt(S(3)) + S(7))*(-c*q + d*(S(1) + sqrt(S(3))) + x*(-c*q + d*(-sqrt(S(3)) + S(1))))), x), x, (q*x + S(1) + sqrt(S(3)))/(-q*x + S(-1) + sqrt(S(3))))/(q*sqrt((-q*x + S(-1))/(q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(a + b*x**S(3)))
    rule392 = ReplacementRule(pattern392, lambda e, a, f, x, b, d, c : With392(e, a, f, x, b, d, c))
    rubi.add(rule392)

    pattern393 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda a: NegQ(a)), CustomConstraint(lambda b: NegQ(b)), CustomConstraint(lambda f, d, b, e, a, x, q, c: NonzeroQ(-e*q + f*(-sqrt(S(3)) + S(1)))))
    def With393(a, f, x, b, d, e, c):
        q = Rt(b/a, S(3))
        return (-c*f + d*e)*Int((q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x)/(-c*q + d*(-sqrt(S(3)) + S(1))) + (-e*q + f*(-sqrt(S(3)) + S(1)))*Int(1/sqrt(a + b*x**S(3)), x)/(-c*q + d*(-sqrt(S(3)) + S(1)))
    rule393 = ReplacementRule(pattern393, lambda a, f, x, b, d, e, c : With393(a, f, x, b, d, e, c))
    rubi.add(rule393)

    pattern394 = Pattern(Integral(x_**WC('m', S(1))/(c_ + x_**n_*WC('d', S(1)) + sqrt(a_ + x_**n_*WC('b', S(1)))*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, a, d, c: ZeroQ(-a*d + b*c, S(0))), CustomConstraint(lambda n, m: IntegerQ((m + S(1))/n)))
    rule394 = ReplacementRule(pattern394, lambda n, a, x, m, b, d, e, c : Subst(Int(x**(S(-1) + (m + S(1))/n)/(c + d*x + e*sqrt(a + b*x)), x), x, x**n)/n)
    rubi.add(rule394)

    pattern395 = Pattern(Integral(WC('u', S(1))/(c_ + x_**n_*WC('d', S(1)) + sqrt(a_ + x_**n_*WC('b', S(1)))*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda b, a, d, c: ZeroQ(-a*d + b*c, S(0))))
    rule395 = ReplacementRule(pattern395, lambda n, a, x, u, b, d, e, c : -a*e*Int(u/(sqrt(a + b*x**n)*(-a*e**S(2) + c**S(2) + c*d*x**n)), x) + c*Int(u/(-a*e**S(2) + c**S(2) + c*d*x**n), x))
    rubi.add(rule395)

    pattern396 = Pattern(Integral((A_ + x_**n_*WC('B', S(1)))/(a_ + x_**S(2)*WC('b', S(1)) + x_**n2_*WC('d', S(1)) + x_**n_*WC('c', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n: NonzeroQ(n + S(-2))), CustomConstraint(lambda n, a, A, B, d: ZeroQ(-A**S(2)*d*(n + S(-1))**S(2) + B**S(2)*a)), CustomConstraint(lambda n, A, B, d, c: ZeroQ(S(2)*A*d*(n + S(-1)) + B*c)))
    rule396 = ReplacementRule(pattern396, lambda n2, n, a, x, A, B, b, d, c : A**S(2)*(n + S(-1))*Subst(Int(1/(A**S(2)*b*x**S(2)*(n + S(-1))**S(2) + a), x), x, x/(A*(n + S(-1)) - B*x**n)))
    rubi.add(rule396)

    pattern397 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('n', S(1))*WC('B', S(1)))/(a_ + x_**n2_*WC('d', S(1)) + x_**WC('k', S(1))*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda A, x: FreeQ(A, x)), CustomConstraint(lambda B, x: FreeQ(B, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda k, m: ZeroQ(k - S(2)*m + S(-2))), CustomConstraint(lambda n, a, m, A, B, d: ZeroQ(-A**S(2)*d*(m - n + S(1))**S(2) + B**S(2)*a*(m + S(1))**S(2))), CustomConstraint(lambda n, A, m, B, d, c: ZeroQ(-S(2)*A*d*(m - n + S(1)) + B*c*(m + S(1)))))
    rule397 = ReplacementRule(pattern397, lambda n2, n, k, a, x, A, m, B, b, d, c : A**S(2)*(m - n + S(1))*Subst(Int(1/(A**S(2)*b*x**S(2)*(m - n + S(1))**S(2) + a), x), x, x**(m + S(1))/(A*(m - n + S(1)) + B*x**n*(m + S(1))))/(m + S(1)))
    rubi.add(rule397)

    pattern398 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, n3: ZeroQ(-S(3)*n + n3)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1))))
    rule398 = ReplacementRule(pattern398, lambda n2, n3, p, n, a, g, f, x, b, d, e, c : -x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*(a*g + c*e) - S(2)*a*c*(-a*f + c*d) + b**S(2)*c*d + x**n*(-a*b**S(2)*g - S(2)*a*c*(-a*g + c*e) + b*c*(a*f + c*d)))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a*b*(a*g + c*e) - S(2)*a*c*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(a*b**S(2)*g*(n*(p + S(2)) + S(1)) - S(2)*a*c*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))) - b*c*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1))), x), x)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule398)

    pattern399 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1))))
    rule399 = ReplacementRule(pattern399, lambda n2, p, n, a, f, x, b, d, e, c : -x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d + x**n*(-S(2)*a*c*e + b*(a*f + c*d)))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a*b*e - S(2)*a*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) - b**S(2)*d*(n*p + n + S(1)) - x**n*(-S(2)*a*c*e*(n*(S(2)*p + S(3)) + S(1)) + b*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1))), x), x)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule399)

    pattern400 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, n3: ZeroQ(-S(3)*n + n3)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1))))
    rule400 = ReplacementRule(pattern400, lambda n2, n3, p, n, a, g, x, b, d, e, c : -x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*(a*g + c*e) - S(2)*a*c**S(2)*d + b**S(2)*c*d + x**n*(-a*b**S(2)*g - S(2)*a*c*(-a*g + c*e) + b*c**S(2)*d))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a*b*(a*g + c*e) + S(2)*a*c**S(2)*d*(S(2)*n*(p + S(1)) + S(1)) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(a*b**S(2)*g*(n*(p + S(2)) + S(1)) - S(2)*a*c*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))) - b*c**S(2)*d*(n*(S(2)*p + S(3)) + S(1))), x), x)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule400)

    pattern401 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n2_*WC('f', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, n3: ZeroQ(-S(3)*n + n3)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1))))
    rule401 = ReplacementRule(pattern401, lambda n2, n3, p, n, a, g, f, x, b, d, c : -x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a**S(2)*b*g - S(2)*a*c*(-a*f + c*d) + b**S(2)*c*d + x**n*(S(2)*a**S(2)*c*g - a*b**S(2)*g + b*c*(a*f + c*d)))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a**S(2)*b*g - S(2)*a*c*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(-S(2)*a**S(2)*c*g*(n + S(1)) + a*b**S(2)*g*(n*(p + S(2)) + S(1)) - b*c*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1))), x), x)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule401)

    pattern402 = Pattern(Integral((d_ + x_**n2_*WC('f', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1))))
    rule402 = ReplacementRule(pattern402, lambda n2, p, n, a, f, x, b, d, c : -x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*(-a*f + c*d) + b**S(2)*d + b*x**n*(a*f + c*d))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))) + Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(S(2)*a*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) + b**S(2)*d*(n*p + n + S(1)) + b*x**n*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1)), x), x)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule402)

    pattern403 = Pattern(Integral((d_ + g_*x_**n3_)*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, n3: ZeroQ(-S(3)*n + n3)), CustomConstraint(lambda b, a, c: NonzeroQ(-S(4)*a*c + b**S(2))), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1))))
    rule403 = ReplacementRule(pattern403, lambda n2, n3, p, n, a, g, x, b, d, c : -x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a**S(2)*b*g - S(2)*a*c**S(2)*d + b**S(2)*c*d + x**n*(S(2)*a**S(2)*c*g - a*b**S(2)*g + b*c**S(2)*d))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a**S(2)*b*g + S(2)*a*c**S(2)*d*(S(2)*n*(p + S(1)) + S(1)) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(-S(2)*a**S(2)*c*g*(n + S(1)) + a*b**S(2)*g*(n*(p + S(2)) + S(1)) - b*c**S(2)*d*(n*(S(2)*p + S(3)) + S(1))), x), x)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule403)

    pattern404 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, n3: ZeroQ(-S(3)*n + n3)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1))))
    rule404 = ReplacementRule(pattern404, lambda n2, n3, p, n, a, g, f, x, d, e, c : x*(a + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c*x**n*(-a*g + c*e) - S(2)*a*c*(-a*f + c*d))/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))) + Int((a + c*x**(S(2)*n))**(p + S(1))*Simp(-S(2)*a*c*x**n*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))) - S(2)*a*c*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))), x), x)/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))))
    rubi.add(rule404)

    pattern405 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1))))
    rule405 = ReplacementRule(pattern405, lambda n2, p, n, a, f, x, d, e, c : x*(a + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c*e*x**n - S(2)*a*(-a*f + c*d))/(S(4)*a**S(2)*c*n*(p + S(1))) + Int((a + c*x**(S(2)*n))**(p + S(1))*Simp(S(2)*a*c*e*x**n*(n*(S(2)*p + S(3)) + S(1)) - S(2)*a*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))), x), x)/(S(4)*a**S(2)*c*n*(p + S(1))))
    rubi.add(rule405)

    pattern406 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n_*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda n2, n: ZeroQ(-S(2)*n + n2)), CustomConstraint(lambda n, n3: ZeroQ(-S(3)*n + n3)), CustomConstraint(lambda p: NegativeIntegerQ(p + S(1))))
    rule406 = ReplacementRule(pattern406, lambda n2, n3, p, n, a, g, x, d, e, c : x*(a + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c**S(2)*d - S(2)*a*c*x**n*(-a*g + c*e))/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))) + Int((a + c*x**(S(2)*n))**(p + S(1))*Simp(S(2)*a*c**S(2)*d*(S(2)*n*(p + S(1)) + S(1)) - S(2)*a*c*x**n*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))), x), x)/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))))
    rubi.add(rule406)

    pattern407 = Pattern(Integral((x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)) + WC('a', S(0)))/(d_ + x_**S(6)*WC('g', S(1)) + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda a, g, f, b, d, c: ZeroQ(-S(12)*a**S(3)*g**S(2) + a**S(2)*c*f**S(2) + S(2)*a*b*g*(a*f + S(3)*c*d) + S(9)*c**S(3)*d**S(2) - c*d*f*(S(6)*a*c + b**S(2)))), CustomConstraint(lambda a, g, f, b, d, e, c: ZeroQ(a**S(3)*c*f**S(2)*g + S(2)*a**S(3)*g**S(2)*(-S(6)*a*g + b*f) - S(3)*a**S(2)*c**S(2)*d*f*g + S(3)*c**S(4)*d**S(2)*e - c**S(3)*d*(-S(12)*a*d*g + a*e*f + S(2)*b*d*f))), CustomConstraint(lambda f, a, d, c: NonzeroQ(-a*f + S(3)*c*d)), CustomConstraint(lambda c, a, b, d, g: NonzeroQ(-S(2)*a**S(2)*g + b*c*d)), CustomConstraint(lambda a, g, f, b, d, c: NonzeroQ(S(4)*a**S(2)*g - a*b*f + b*c*d)), CustomConstraint(lambda c, a, f, b, d, g: PosQ((S(12)*a**S(2)*g**S(2) - a*c*f**S(2) + f*(-S(2)*a*b*g + S(3)*c**S(2)*d))/(c*g*(-a*f + S(3)*c*d)))), )
    def With407(a, g, f, x, b, d, e, c):
        q = Rt((S(12)*a**S(2)*g**S(2) - a*c*f**S(2) + f*(-S(2)*a*b*g + S(3)*c**S(2)*d))/(c*g*(-a*f + S(3)*c*d)), S(2))
        r = Rt((a*c*f**S(2) - f*(S(2)*a*b*g + S(3)*c**S(2)*d) + S(4)*g*(a**S(2)*g + b*c*d))/(c*g*(-a*f + S(3)*c*d)), S(2))
        return -c*atan((r - S(2)*x)/q)/(g*q) + c*atan((r + S(2)*x)/q)/(g*q) - c*atan(x*(-a*f + S(3)*c*d)*(S(6)*a**S(2)*b*g**S(2) - S(2)*a**S(2)*c*f*g - a*b**S(2)*f*g + b*c**S(2)*d*f + c**S(2)*g*x**S(4)*(-a*f + S(3)*c*d) + c*x**S(2)*(S(2)*a**S(2)*g**S(2) - a*c*f**S(2) - b*c*d*g + S(3)*c**S(2)*d*f))/(g*q*(-S(2)*a**S(2)*g + b*c*d)*(S(4)*a**S(2)*g - a*b*f + b*c*d)))/(g*q)
    rule407 = ReplacementRule(pattern407, lambda a, g, f, x, b, d, e, c : With407(a, g, f, x, b, d, e, c))
    rubi.add(rule407)

    pattern408 = Pattern(Integral((x_**S(4)*WC('c', S(1)) + WC('a', S(0)))/(d_ + x_**S(6)*WC('g', S(1)) + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda c, x: FreeQ(c, x)), CustomConstraint(lambda d, x: FreeQ(d, x)), CustomConstraint(lambda e, x: FreeQ(e, x)), CustomConstraint(lambda f, x: FreeQ(f, x)), CustomConstraint(lambda g, x: FreeQ(g, x)), CustomConstraint(lambda c, a, f, d, g: ZeroQ(-S(12)*a**S(3)*g**S(2) + a**S(2)*c*f**S(2) - S(6)*a*c**S(2)*d*f + S(9)*c**S(3)*d**S(2))), CustomConstraint(lambda a, g, f, d, e, c: ZeroQ(-S(12)*a**S(4)*g**S(3) + a**S(3)*c*f**S(2)*g - S(3)*a**S(2)*c**S(2)*d*f*g - a*c**S(3)*d*(-S(12)*d*g + e*f) + S(3)*c**S(4)*d**S(2)*e)), CustomConstraint(lambda f, a, d, c: NonzeroQ(-a*f + S(3)*c*d)), CustomConstraint(lambda c, a, f, d, g: PosQ((S(12)*a**S(2)*g**S(2) - a*c*f**S(2) + S(3)*c**S(2)*d*f)/(c*g*(-a*f + S(3)*c*d)))), )
    def With408(c, a, f, x, d, e, g):
        q = Rt((S(12)*a**S(2)*g**S(2) - a*c*f**S(2) + S(3)*c**S(2)*d*f)/(c*g*(-a*f + S(3)*c*d)), S(2))
        r = Rt((S(4)*a**S(2)*g**S(2) + a*c*f**S(2) - S(3)*c**S(2)*d*f)/(c*g*(-a*f + S(3)*c*d)), S(2))
        return -c*atan((r - S(2)*x)/q)/(g*q) + c*atan((r + S(2)*x)/q)/(g*q) - c*atan(c*x*(-a*f + S(3)*c*d)*(S(2)*a**S(2)*f*g - c*g*x**S(4)*(-a*f + S(3)*c*d) - x**S(2)*(S(2)*a**S(2)*g**S(2) - a*c*f**S(2) + S(3)*c**S(2)*d*f))/(S(8)*a**S(4)*g**S(3)*q))/(g*q)
    rule408 = ReplacementRule(pattern408, lambda c, a, f, x, d, e, g : With408(c, a, f, x, d, e, g))
    rubi.add(rule408)

    pattern409 = Pattern(Integral(u_*v_**p_, x_), CustomConstraint(lambda p: RationalQ(p)), CustomConstraint(lambda p: Less(p, S(-1))), CustomConstraint(lambda x, u: PolynomialQ(u, x)), CustomConstraint(lambda v, x: PolynomialQ(v, x)), CustomConstraint(lambda v: SumQ(v)), CustomConstraint(lambda v, x, u: Not(BinomialQ(v, x) & MonomialQ(u, x))), CustomConstraint(lambda v, x, u: Not(ZeroQ(Coefficient(u, x, S(0))) & ZeroQ(Coefficient(v, x, S(0))))), CustomConstraint(lambda v, n, w, LessEqual, Less, m, x, p, c, u: FalseQ(DerivativeDivides(v, u, x)) & Less(m + n*p, S(-1)) & Inequality(S(1), Less, n, LessEqual, m + S(1))))
    def With409(v, p, x, u):
        m = Exponent(u, x)
        n = Exponent(v, x)
        return Module(List(Set(c, Coefficient(u, x, m)/((m + n*p + S(1))*Coefficient(v, x, n))), w), CompoundExpression(Set(c, Coefficient(u, x, m)/((m + n*p + S(1))*Coefficient(v, x, n))), Set(w, Apart(-c*x**(m - n)*(v*(m - n + S(1)) + x*(p + S(1))*D(v, x)) + u, x)), If(ZeroQ(w), c*v**(p + S(1))*x**(m - n + S(1)), c*v**(p + S(1))*x**(m - n + S(1)) + Int(v**p*w, x))))
    rule409 = ReplacementRule(pattern409, lambda v, p, x, u : With409(v, p, x, u))
    rubi.add(rule409)

    return rubi
