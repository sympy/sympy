from sympy.core import Basic
from sympy.utilities.pytest import XFAIL

def _test_args(obj):
    return all(isinstance(arg, Basic) for arg in obj.args)

@XFAIL
def test_sympy__assumptions__assume__AppliedPredicate():
    from sympy.assumptions.assume import AppliedPredicate
    assert _test_args(AppliedPredicate())

@XFAIL
def test_sympy__assumptions__assume__Predicate():
    from sympy.assumptions.assume import Predicate
    assert _test_args(Predicate())

@XFAIL
def test_sympy__combinatorics__permutations__Permutation():
    from sympy.combinatorics.permutations import Permutation
    assert _test_args(Permutation())

@XFAIL
def test_sympy__combinatorics__prufer__Prufer():
    from sympy.combinatorics.prufer import Prufer
    assert _test_args(Prufer())

@XFAIL
def test_sympy__concrete__products__Product():
    from sympy.concrete.products import Product
    assert _test_args(Product())

@XFAIL
def test_sympy__concrete__summations__Sum():
    from sympy.concrete.summations import Sum
    assert _test_args(Sum())

@XFAIL
def test_sympy__core__add__Add():
    from sympy.core.add import Add
    assert _test_args(Add())

@XFAIL
def test_sympy__core__basic__Atom():
    from sympy.core.basic import Atom
    assert _test_args(Atom())

@XFAIL
def test_sympy__core__basic__Basic():
    from sympy.core.basic import Basic
    assert _test_args(Basic())

@XFAIL
def test_sympy__core__containers__Dict():
    from sympy.core.containers import Dict
    assert _test_args(Dict())

@XFAIL
def test_sympy__core__containers__Tuple():
    from sympy.core.containers import Tuple
    assert _test_args(Tuple())

@XFAIL
def test_sympy__core__expr__AtomicExpr():
    from sympy.core.expr import AtomicExpr
    assert _test_args(AtomicExpr())

@XFAIL
def test_sympy__core__expr__Expr():
    from sympy.core.expr import Expr
    assert _test_args(Expr())

@XFAIL
def test_sympy__core__function__Application():
    from sympy.core.function import Application
    assert _test_args(Application())

@XFAIL
def test_sympy__core__function__AppliedUndef():
    from sympy.core.function import AppliedUndef
    assert _test_args(AppliedUndef())

@XFAIL
def test_sympy__core__function__Derivative():
    from sympy.core.function import Derivative
    assert _test_args(Derivative())

@XFAIL
def test_sympy__core__function__Function():
    from sympy.core.function import Function
    assert _test_args(Function())

@XFAIL
def test_sympy__core__function__Lambda():
    from sympy.core.function import Lambda
    assert _test_args(Lambda())

@XFAIL
def test_sympy__core__function__Subs():
    from sympy.core.function import Subs
    assert _test_args(Subs())

@XFAIL
def test_sympy__core__function__WildFunction():
    from sympy.core.function import WildFunction
    assert _test_args(WildFunction())

@XFAIL
def test_sympy__core__mul__Mul():
    from sympy.core.mul import Mul
    assert _test_args(Mul())

@XFAIL
def test_sympy__core__numbers__Catalan():
    from sympy.core.numbers import Catalan
    assert _test_args(Catalan())

@XFAIL
def test_sympy__core__numbers__ComplexInfinity():
    from sympy.core.numbers import ComplexInfinity
    assert _test_args(ComplexInfinity())

@XFAIL
def test_sympy__core__numbers__EulerGamma():
    from sympy.core.numbers import EulerGamma
    assert _test_args(EulerGamma())

@XFAIL
def test_sympy__core__numbers__Exp1():
    from sympy.core.numbers import Exp1
    assert _test_args(Exp1())

@XFAIL
def test_sympy__core__numbers__Float():
    from sympy.core.numbers import Float
    assert _test_args(Float())

@XFAIL
def test_sympy__core__numbers__GoldenRatio():
    from sympy.core.numbers import GoldenRatio
    assert _test_args(GoldenRatio())

@XFAIL
def test_sympy__core__numbers__Half():
    from sympy.core.numbers import Half
    assert _test_args(Half())

@XFAIL
def test_sympy__core__numbers__ImaginaryUnit():
    from sympy.core.numbers import ImaginaryUnit
    assert _test_args(ImaginaryUnit())

@XFAIL
def test_sympy__core__numbers__Infinity():
    from sympy.core.numbers import Infinity
    assert _test_args(Infinity())

@XFAIL
def test_sympy__core__numbers__Integer():
    from sympy.core.numbers import Integer
    assert _test_args(Integer())

@XFAIL
def test_sympy__core__numbers__IntegerConstant():
    from sympy.core.numbers import IntegerConstant
    assert _test_args(IntegerConstant())

@XFAIL
def test_sympy__core__numbers__NaN():
    from sympy.core.numbers import NaN
    assert _test_args(NaN())

@XFAIL
def test_sympy__core__numbers__NegativeInfinity():
    from sympy.core.numbers import NegativeInfinity
    assert _test_args(NegativeInfinity())

@XFAIL
def test_sympy__core__numbers__NegativeOne():
    from sympy.core.numbers import NegativeOne
    assert _test_args(NegativeOne())

@XFAIL
def test_sympy__core__numbers__Number():
    from sympy.core.numbers import Number
    assert _test_args(Number())

@XFAIL
def test_sympy__core__numbers__NumberSymbol():
    from sympy.core.numbers import NumberSymbol
    assert _test_args(NumberSymbol())

@XFAIL
def test_sympy__core__numbers__One():
    from sympy.core.numbers import One
    assert _test_args(One())

@XFAIL
def test_sympy__core__numbers__Pi():
    from sympy.core.numbers import Pi
    assert _test_args(Pi())

@XFAIL
def test_sympy__core__numbers__Rational():
    from sympy.core.numbers import Rational
    assert _test_args(Rational())

@XFAIL
def test_sympy__core__numbers__RationalConstant():
    from sympy.core.numbers import RationalConstant
    assert _test_args(RationalConstant())

@XFAIL
def test_sympy__core__numbers__Zero():
    from sympy.core.numbers import Zero
    assert _test_args(Zero())

@XFAIL
def test_sympy__core__operations__AssocOp():
    from sympy.core.operations import AssocOp
    assert _test_args(AssocOp())

@XFAIL
def test_sympy__core__operations__LatticeOp():
    from sympy.core.operations import LatticeOp
    assert _test_args(LatticeOp())

@XFAIL
def test_sympy__core__power__Pow():
    from sympy.core.power import Pow
    assert _test_args(Pow())

@XFAIL
def test_sympy__core__relational__Equality():
    from sympy.core.relational import Equality
    assert _test_args(Equality())

@XFAIL
def test_sympy__core__relational__Inequality():
    from sympy.core.relational import Inequality
    assert _test_args(Inequality())

@XFAIL
def test_sympy__core__relational__Relational():
    from sympy.core.relational import Relational
    assert _test_args(Relational())

@XFAIL
def test_sympy__core__relational__StrictInequality():
    from sympy.core.relational import StrictInequality
    assert _test_args(StrictInequality())

@XFAIL
def test_sympy__core__relational__Unequality():
    from sympy.core.relational import Unequality
    assert _test_args(Unequality())

@XFAIL
def test_sympy__core__sets__CountableSet():
    from sympy.core.sets import CountableSet
    assert _test_args(CountableSet())

@XFAIL
def test_sympy__core__sets__EmptySet():
    from sympy.core.sets import EmptySet
    assert _test_args(EmptySet())

@XFAIL
def test_sympy__core__sets__FiniteSet():
    from sympy.core.sets import FiniteSet
    assert _test_args(FiniteSet())

@XFAIL
def test_sympy__core__sets__Interval():
    from sympy.core.sets import Interval
    assert _test_args(Interval())

@XFAIL
def test_sympy__core__sets__ProductSet():
    from sympy.core.sets import ProductSet
    assert _test_args(ProductSet())

@XFAIL
def test_sympy__core__sets__RealFiniteSet():
    from sympy.core.sets import RealFiniteSet
    assert _test_args(RealFiniteSet())

@XFAIL
def test_sympy__core__sets__RealSet():
    from sympy.core.sets import RealSet
    assert _test_args(RealSet())

@XFAIL
def test_sympy__core__sets__RealUnion():
    from sympy.core.sets import RealUnion
    assert _test_args(RealUnion())

@XFAIL
def test_sympy__core__sets__Set():
    from sympy.core.sets import Set
    assert _test_args(Set())

@XFAIL
def test_sympy__core__sets__Union():
    from sympy.core.sets import Union
    assert _test_args(Union())

@XFAIL
def test_sympy__core__symbol__Dummy():
    from sympy.core.symbol import Dummy
    assert _test_args(Dummy())

@XFAIL
def test_sympy__core__symbol__Symbol():
    from sympy.core.symbol import Symbol
    assert _test_args(Symbol())

@XFAIL
def test_sympy__core__symbol__Wild():
    from sympy.core.symbol import Wild
    assert _test_args(Wild())

@XFAIL
def test_sympy__functions__combinatorial__factorials__CombinatorialFunction():
    from sympy.functions.combinatorial.factorials import CombinatorialFunction
    assert _test_args(CombinatorialFunction())

@XFAIL
def test_sympy__functions__combinatorial__factorials__FallingFactorial():
    from sympy.functions.combinatorial.factorials import FallingFactorial
    assert _test_args(FallingFactorial())

@XFAIL
def test_sympy__functions__combinatorial__factorials__MultiFactorial():
    from sympy.functions.combinatorial.factorials import MultiFactorial
    assert _test_args(MultiFactorial())

@XFAIL
def test_sympy__functions__combinatorial__factorials__RisingFactorial():
    from sympy.functions.combinatorial.factorials import RisingFactorial
    assert _test_args(RisingFactorial())

@XFAIL
def test_sympy__functions__combinatorial__factorials__binomial():
    from sympy.functions.combinatorial.factorials import binomial
    assert _test_args(binomial())

@XFAIL
def test_sympy__functions__combinatorial__factorials__factorial():
    from sympy.functions.combinatorial.factorials import factorial
    assert _test_args(factorial())

@XFAIL
def test_sympy__functions__combinatorial__factorials__factorial2():
    from sympy.functions.combinatorial.factorials import factorial2
    assert _test_args(factorial2())

@XFAIL
def test_sympy__functions__combinatorial__numbers__bell():
    from sympy.functions.combinatorial.numbers import bell
    assert _test_args(bell())

@XFAIL
def test_sympy__functions__combinatorial__numbers__bernoulli():
    from sympy.functions.combinatorial.numbers import bernoulli
    assert _test_args(bernoulli())

@XFAIL
def test_sympy__functions__combinatorial__numbers__catalan():
    from sympy.functions.combinatorial.numbers import catalan
    assert _test_args(catalan())

@XFAIL
def test_sympy__functions__combinatorial__numbers__euler():
    from sympy.functions.combinatorial.numbers import euler
    assert _test_args(euler())

@XFAIL
def test_sympy__functions__combinatorial__numbers__fibonacci():
    from sympy.functions.combinatorial.numbers import fibonacci
    assert _test_args(fibonacci())

@XFAIL
def test_sympy__functions__combinatorial__numbers__harmonic():
    from sympy.functions.combinatorial.numbers import harmonic
    assert _test_args(harmonic())

@XFAIL
def test_sympy__functions__combinatorial__numbers__lucas():
    from sympy.functions.combinatorial.numbers import lucas
    assert _test_args(lucas())

@XFAIL
def test_sympy__functions__elementary__complexes__Abs():
    from sympy.functions.elementary.complexes import Abs
    assert _test_args(Abs())

@XFAIL
def test_sympy__functions__elementary__complexes__arg():
    from sympy.functions.elementary.complexes import arg
    assert _test_args(arg())

@XFAIL
def test_sympy__functions__elementary__complexes__conjugate():
    from sympy.functions.elementary.complexes import conjugate
    assert _test_args(conjugate())

@XFAIL
def test_sympy__functions__elementary__complexes__im():
    from sympy.functions.elementary.complexes import im
    assert _test_args(im())

@XFAIL
def test_sympy__functions__elementary__complexes__re():
    from sympy.functions.elementary.complexes import re
    assert _test_args(re())

@XFAIL
def test_sympy__functions__elementary__complexes__sign():
    from sympy.functions.elementary.complexes import sign
    assert _test_args(sign())

@XFAIL
def test_sympy__functions__elementary__exponential__LambertW():
    from sympy.functions.elementary.exponential import LambertW
    assert _test_args(LambertW())

@XFAIL
def test_sympy__functions__elementary__exponential__exp():
    from sympy.functions.elementary.exponential import exp
    assert _test_args(exp())

@XFAIL
def test_sympy__functions__elementary__exponential__log():
    from sympy.functions.elementary.exponential import log
    assert _test_args(log())

@XFAIL
def test_sympy__functions__elementary__hyperbolic__HyperbolicFunction():
    from sympy.functions.elementary.hyperbolic import HyperbolicFunction
    assert _test_args(HyperbolicFunction())

@XFAIL
def test_sympy__functions__elementary__hyperbolic__acosh():
    from sympy.functions.elementary.hyperbolic import acosh
    assert _test_args(acosh())

@XFAIL
def test_sympy__functions__elementary__hyperbolic__acoth():
    from sympy.functions.elementary.hyperbolic import acoth
    assert _test_args(acoth())

@XFAIL
def test_sympy__functions__elementary__hyperbolic__asinh():
    from sympy.functions.elementary.hyperbolic import asinh
    assert _test_args(asinh())

@XFAIL
def test_sympy__functions__elementary__hyperbolic__atanh():
    from sympy.functions.elementary.hyperbolic import atanh
    assert _test_args(atanh())

@XFAIL
def test_sympy__functions__elementary__hyperbolic__cosh():
    from sympy.functions.elementary.hyperbolic import cosh
    assert _test_args(cosh())

@XFAIL
def test_sympy__functions__elementary__hyperbolic__coth():
    from sympy.functions.elementary.hyperbolic import coth
    assert _test_args(coth())

@XFAIL
def test_sympy__functions__elementary__hyperbolic__sinh():
    from sympy.functions.elementary.hyperbolic import sinh
    assert _test_args(sinh())

@XFAIL
def test_sympy__functions__elementary__hyperbolic__tanh():
    from sympy.functions.elementary.hyperbolic import tanh
    assert _test_args(tanh())

@XFAIL
def test_sympy__functions__elementary__integers__RoundFunction():
    from sympy.functions.elementary.integers import RoundFunction
    assert _test_args(RoundFunction())

@XFAIL
def test_sympy__functions__elementary__integers__ceiling():
    from sympy.functions.elementary.integers import ceiling
    assert _test_args(ceiling())

@XFAIL
def test_sympy__functions__elementary__integers__floor():
    from sympy.functions.elementary.integers import floor
    assert _test_args(floor())

@XFAIL
def test_sympy__functions__elementary__miscellaneous__IdentityFunction():
    from sympy.functions.elementary.miscellaneous import IdentityFunction
    assert _test_args(IdentityFunction())

@XFAIL
def test_sympy__functions__elementary__miscellaneous__Max():
    from sympy.functions.elementary.miscellaneous import Max
    assert _test_args(Max())

@XFAIL
def test_sympy__functions__elementary__miscellaneous__Min():
    from sympy.functions.elementary.miscellaneous import Min
    assert _test_args(Min())

@XFAIL
def test_sympy__functions__elementary__miscellaneous__MinMaxBase():
    from sympy.functions.elementary.miscellaneous import MinMaxBase
    assert _test_args(MinMaxBase())

@XFAIL
def test_sympy__functions__elementary__piecewise__ExprCondPair():
    from sympy.functions.elementary.piecewise import ExprCondPair
    assert _test_args(ExprCondPair())

@XFAIL
def test_sympy__functions__elementary__piecewise__Piecewise():
    from sympy.functions.elementary.piecewise import Piecewise
    assert _test_args(Piecewise())

@XFAIL
def test_sympy__functions__elementary__trigonometric__TrigonometricFunction():
    from sympy.functions.elementary.trigonometric import TrigonometricFunction
    assert _test_args(TrigonometricFunction())

@XFAIL
def test_sympy__functions__elementary__trigonometric__acos():
    from sympy.functions.elementary.trigonometric import acos
    assert _test_args(acos())

@XFAIL
def test_sympy__functions__elementary__trigonometric__acot():
    from sympy.functions.elementary.trigonometric import acot
    assert _test_args(acot())

@XFAIL
def test_sympy__functions__elementary__trigonometric__asin():
    from sympy.functions.elementary.trigonometric import asin
    assert _test_args(asin())

@XFAIL
def test_sympy__functions__elementary__trigonometric__atan():
    from sympy.functions.elementary.trigonometric import atan
    assert _test_args(atan())

@XFAIL
def test_sympy__functions__elementary__trigonometric__atan2():
    from sympy.functions.elementary.trigonometric import atan2
    assert _test_args(atan2())

@XFAIL
def test_sympy__functions__elementary__trigonometric__cos():
    from sympy.functions.elementary.trigonometric import cos
    assert _test_args(cos())

@XFAIL
def test_sympy__functions__elementary__trigonometric__cot():
    from sympy.functions.elementary.trigonometric import cot
    assert _test_args(cot())

@XFAIL
def test_sympy__functions__elementary__trigonometric__sin():
    from sympy.functions.elementary.trigonometric import sin
    assert _test_args(sin())

@XFAIL
def test_sympy__functions__elementary__trigonometric__tan():
    from sympy.functions.elementary.trigonometric import tan
    assert _test_args(tan())

@XFAIL
def test_sympy__functions__special__bessel__BesselBase():
    from sympy.functions.special.bessel import BesselBase
    assert _test_args(BesselBase())

@XFAIL
def test_sympy__functions__special__bessel__SphericalBesselBase():
    from sympy.functions.special.bessel import SphericalBesselBase
    assert _test_args(SphericalBesselBase())

@XFAIL
def test_sympy__functions__special__bessel__besseli():
    from sympy.functions.special.bessel import besseli
    assert _test_args(besseli())

@XFAIL
def test_sympy__functions__special__bessel__besselj():
    from sympy.functions.special.bessel import besselj
    assert _test_args(besselj())

@XFAIL
def test_sympy__functions__special__bessel__besselk():
    from sympy.functions.special.bessel import besselk
    assert _test_args(besselk())

@XFAIL
def test_sympy__functions__special__bessel__bessely():
    from sympy.functions.special.bessel import bessely
    assert _test_args(bessely())

@XFAIL
def test_sympy__functions__special__bessel__hankel1():
    from sympy.functions.special.bessel import hankel1
    assert _test_args(hankel1())

@XFAIL
def test_sympy__functions__special__bessel__hankel2():
    from sympy.functions.special.bessel import hankel2
    assert _test_args(hankel2())

@XFAIL
def test_sympy__functions__special__bessel__jn():
    from sympy.functions.special.bessel import jn
    assert _test_args(jn())

@XFAIL
def test_sympy__functions__special__bessel__yn():
    from sympy.functions.special.bessel import yn
    assert _test_args(yn())

@XFAIL
def test_sympy__functions__special__delta_functions__DiracDelta():
    from sympy.functions.special.delta_functions import DiracDelta
    assert _test_args(DiracDelta())

@XFAIL
def test_sympy__functions__special__delta_functions__Heaviside():
    from sympy.functions.special.delta_functions import Heaviside
    assert _test_args(Heaviside())

@XFAIL
def test_sympy__functions__special__error_functions__erf():
    from sympy.functions.special.error_functions import erf
    assert _test_args(erf())

@XFAIL
def test_sympy__functions__special__gamma_functions__gamma():
    from sympy.functions.special.gamma_functions import gamma
    assert _test_args(gamma())

@XFAIL
def test_sympy__functions__special__gamma_functions__loggamma():
    from sympy.functions.special.gamma_functions import loggamma
    assert _test_args(loggamma())

@XFAIL
def test_sympy__functions__special__gamma_functions__lowergamma():
    from sympy.functions.special.gamma_functions import lowergamma
    assert _test_args(lowergamma())

@XFAIL
def test_sympy__functions__special__gamma_functions__polygamma():
    from sympy.functions.special.gamma_functions import polygamma
    assert _test_args(polygamma())

@XFAIL
def test_sympy__functions__special__gamma_functions__uppergamma():
    from sympy.functions.special.gamma_functions import uppergamma
    assert _test_args(uppergamma())

@XFAIL
def test_sympy__functions__special__hyper__TupleParametersBase():
    from sympy.functions.special.hyper import TupleParametersBase
    assert _test_args(TupleParametersBase())

@XFAIL
def test_sympy__functions__special__hyper__hyper():
    from sympy.functions.special.hyper import hyper
    assert _test_args(hyper())

@XFAIL
def test_sympy__functions__special__hyper__meijerg():
    from sympy.functions.special.hyper import meijerg
    assert _test_args(meijerg())

@XFAIL
def test_sympy__functions__special__polynomials__PolynomialSequence():
    from sympy.functions.special.polynomials import PolynomialSequence
    assert _test_args(PolynomialSequence())

@XFAIL
def test_sympy__functions__special__polynomials__assoc_legendre():
    from sympy.functions.special.polynomials import assoc_legendre
    assert _test_args(assoc_legendre())

@XFAIL
def test_sympy__functions__special__polynomials__chebyshevt():
    from sympy.functions.special.polynomials import chebyshevt
    assert _test_args(chebyshevt())

@XFAIL
def test_sympy__functions__special__polynomials__chebyshevt_root():
    from sympy.functions.special.polynomials import chebyshevt_root
    assert _test_args(chebyshevt_root())

@XFAIL
def test_sympy__functions__special__polynomials__chebyshevu():
    from sympy.functions.special.polynomials import chebyshevu
    assert _test_args(chebyshevu())

@XFAIL
def test_sympy__functions__special__polynomials__chebyshevu_root():
    from sympy.functions.special.polynomials import chebyshevu_root
    assert _test_args(chebyshevu_root())

@XFAIL
def test_sympy__functions__special__polynomials__hermite():
    from sympy.functions.special.polynomials import hermite
    assert _test_args(hermite())

@XFAIL
def test_sympy__functions__special__polynomials__legendre():
    from sympy.functions.special.polynomials import legendre
    assert _test_args(legendre())

@XFAIL
def test_sympy__functions__special__tensor_functions__Dij():
    from sympy.functions.special.tensor_functions import Dij
    assert _test_args(Dij())

@XFAIL
def test_sympy__functions__special__tensor_functions__LeviCivita():
    from sympy.functions.special.tensor_functions import LeviCivita
    assert _test_args(LeviCivita())

@XFAIL
def test_sympy__functions__special__zeta_functions__dirichlet_eta():
    from sympy.functions.special.zeta_functions import dirichlet_eta
    assert _test_args(dirichlet_eta())

@XFAIL
def test_sympy__functions__special__zeta_functions__zeta():
    from sympy.functions.special.zeta_functions import zeta
    assert _test_args(zeta())

@XFAIL
def test_sympy__integrals__integrals__Integral():
    from sympy.integrals.integrals import Integral
    assert _test_args(Integral())

@XFAIL
def test_sympy__logic__boolalg__And():
    from sympy.logic.boolalg import And
    assert _test_args(And())

@XFAIL
def test_sympy__logic__boolalg__Boolean():
    from sympy.logic.boolalg import Boolean
    assert _test_args(Boolean())

@XFAIL
def test_sympy__logic__boolalg__BooleanFunction():
    from sympy.logic.boolalg import BooleanFunction
    assert _test_args(BooleanFunction())

@XFAIL
def test_sympy__logic__boolalg__Equivalent():
    from sympy.logic.boolalg import Equivalent
    assert _test_args(Equivalent())

@XFAIL
def test_sympy__logic__boolalg__ITE():
    from sympy.logic.boolalg import ITE
    assert _test_args(ITE())

@XFAIL
def test_sympy__logic__boolalg__Implies():
    from sympy.logic.boolalg import Implies
    assert _test_args(Implies())

@XFAIL
def test_sympy__logic__boolalg__Nand():
    from sympy.logic.boolalg import Nand
    assert _test_args(Nand())

@XFAIL
def test_sympy__logic__boolalg__Nor():
    from sympy.logic.boolalg import Nor
    assert _test_args(Nor())

@XFAIL
def test_sympy__logic__boolalg__Not():
    from sympy.logic.boolalg import Not
    assert _test_args(Not())

@XFAIL
def test_sympy__logic__boolalg__Or():
    from sympy.logic.boolalg import Or
    assert _test_args(Or())

@XFAIL
def test_sympy__logic__boolalg__Xor():
    from sympy.logic.boolalg import Xor
    assert _test_args(Xor())

@XFAIL
def test_sympy__matrices__expressions__blockmatrix__BlockDiagMatrix():
    from sympy.matrices.expressions.blockmatrix import BlockDiagMatrix
    assert _test_args(BlockDiagMatrix())

@XFAIL
def test_sympy__matrices__expressions__blockmatrix__BlockMatrix():
    from sympy.matrices.expressions.blockmatrix import BlockMatrix
    assert _test_args(BlockMatrix())

@XFAIL
def test_sympy__matrices__expressions__inverse__Inverse():
    from sympy.matrices.expressions.inverse import Inverse
    assert _test_args(Inverse())

@XFAIL
def test_sympy__matrices__expressions__matadd__MatAdd():
    from sympy.matrices.expressions.matadd import MatAdd
    assert _test_args(MatAdd())

@XFAIL
def test_sympy__matrices__expressions__matexpr__Identity():
    from sympy.matrices.expressions.matexpr import Identity
    assert _test_args(Identity())

@XFAIL
def test_sympy__matrices__expressions__matexpr__MatrixExpr():
    from sympy.matrices.expressions.matexpr import MatrixExpr
    assert _test_args(MatrixExpr())

@XFAIL
def test_sympy__matrices__expressions__matexpr__MatrixSymbol():
    from sympy.matrices.expressions.matexpr import MatrixSymbol
    assert _test_args(MatrixSymbol())

@XFAIL
def test_sympy__matrices__expressions__matexpr__ZeroMatrix():
    from sympy.matrices.expressions.matexpr import ZeroMatrix
    assert _test_args(ZeroMatrix())

@XFAIL
def test_sympy__matrices__expressions__matmul__MatMul():
    from sympy.matrices.expressions.matmul import MatMul
    assert _test_args(MatMul())

@XFAIL
def test_sympy__matrices__expressions__matpow__MatPow():
    from sympy.matrices.expressions.matpow import MatPow
    assert _test_args(MatPow())

@XFAIL
def test_sympy__matrices__expressions__transpose__Transpose():
    from sympy.matrices.expressions.transpose import Transpose
    assert _test_args(Transpose())

@XFAIL
def test_sympy__physics__gaussopt__BeamParameter():
    from sympy.physics.gaussopt import BeamParameter
    assert _test_args(BeamParameter())

@XFAIL
def test_sympy__physics__paulialgebra__Pauli():
    from sympy.physics.paulialgebra import Pauli
    assert _test_args(Pauli())

@XFAIL
def test_sympy__physics__quantum__anticommutator__AntiCommutator():
    from sympy.physics.quantum.anticommutator import AntiCommutator
    assert _test_args(AntiCommutator())

@XFAIL
def test_sympy__physics__quantum__cartesian__PositionBra3D():
    from sympy.physics.quantum.cartesian import PositionBra3D
    assert _test_args(PositionBra3D())

@XFAIL
def test_sympy__physics__quantum__cartesian__PositionKet3D():
    from sympy.physics.quantum.cartesian import PositionKet3D
    assert _test_args(PositionKet3D())

@XFAIL
def test_sympy__physics__quantum__cartesian__PositionState3D():
    from sympy.physics.quantum.cartesian import PositionState3D
    assert _test_args(PositionState3D())

@XFAIL
def test_sympy__physics__quantum__cartesian__PxBra():
    from sympy.physics.quantum.cartesian import PxBra
    assert _test_args(PxBra())

@XFAIL
def test_sympy__physics__quantum__cartesian__PxKet():
    from sympy.physics.quantum.cartesian import PxKet
    assert _test_args(PxKet())

@XFAIL
def test_sympy__physics__quantum__cartesian__PxOp():
    from sympy.physics.quantum.cartesian import PxOp
    assert _test_args(PxOp())

@XFAIL
def test_sympy__physics__quantum__cartesian__XBra():
    from sympy.physics.quantum.cartesian import XBra
    assert _test_args(XBra())

@XFAIL
def test_sympy__physics__quantum__cartesian__XKet():
    from sympy.physics.quantum.cartesian import XKet
    assert _test_args(XKet())

@XFAIL
def test_sympy__physics__quantum__cartesian__XOp():
    from sympy.physics.quantum.cartesian import XOp
    assert _test_args(XOp())

@XFAIL
def test_sympy__physics__quantum__cartesian__YOp():
    from sympy.physics.quantum.cartesian import YOp
    assert _test_args(YOp())

@XFAIL
def test_sympy__physics__quantum__cartesian__ZOp():
    from sympy.physics.quantum.cartesian import ZOp
    assert _test_args(ZOp())

@XFAIL
def test_sympy__physics__quantum__cg__CG():
    from sympy.physics.quantum.cg import CG
    assert _test_args(CG())

@XFAIL
def test_sympy__physics__quantum__cg__Wigner3j():
    from sympy.physics.quantum.cg import Wigner3j
    assert _test_args(Wigner3j())

@XFAIL
def test_sympy__physics__quantum__constants__HBar():
    from sympy.physics.quantum.constants import HBar
    assert _test_args(HBar())

@XFAIL
def test_sympy__physics__quantum__gate__CGate():
    from sympy.physics.quantum.gate import CGate
    assert _test_args(CGate())

@XFAIL
def test_sympy__physics__quantum__gate__CNotGate():
    from sympy.physics.quantum.gate import CNotGate
    assert _test_args(CNotGate())

@XFAIL
def test_sympy__physics__quantum__gate__Gate():
    from sympy.physics.quantum.gate import Gate
    assert _test_args(Gate())

@XFAIL
def test_sympy__physics__quantum__gate__HadamardGate():
    from sympy.physics.quantum.gate import HadamardGate
    assert _test_args(HadamardGate())

@XFAIL
def test_sympy__physics__quantum__gate__IdentityGate():
    from sympy.physics.quantum.gate import IdentityGate
    assert _test_args(IdentityGate())

@XFAIL
def test_sympy__physics__quantum__gate__OneQubitGate():
    from sympy.physics.quantum.gate import OneQubitGate
    assert _test_args(OneQubitGate())

@XFAIL
def test_sympy__physics__quantum__gate__PhaseGate():
    from sympy.physics.quantum.gate import PhaseGate
    assert _test_args(PhaseGate())

@XFAIL
def test_sympy__physics__quantum__gate__SwapGate():
    from sympy.physics.quantum.gate import SwapGate
    assert _test_args(SwapGate())

@XFAIL
def test_sympy__physics__quantum__gate__TGate():
    from sympy.physics.quantum.gate import TGate
    assert _test_args(TGate())

@XFAIL
def test_sympy__physics__quantum__gate__TwoQubitGate():
    from sympy.physics.quantum.gate import TwoQubitGate
    assert _test_args(TwoQubitGate())

@XFAIL
def test_sympy__physics__quantum__gate__UGate():
    from sympy.physics.quantum.gate import UGate
    assert _test_args(UGate())

@XFAIL
def test_sympy__physics__quantum__gate__XGate():
    from sympy.physics.quantum.gate import XGate
    assert _test_args(XGate())

@XFAIL
def test_sympy__physics__quantum__gate__YGate():
    from sympy.physics.quantum.gate import YGate
    assert _test_args(YGate())

@XFAIL
def test_sympy__physics__quantum__gate__ZGate():
    from sympy.physics.quantum.gate import ZGate
    assert _test_args(ZGate())

@XFAIL
def test_sympy__physics__quantum__grover__OracleGate():
    from sympy.physics.quantum.grover import OracleGate
    assert _test_args(OracleGate())

@XFAIL
def test_sympy__physics__quantum__grover__WGate():
    from sympy.physics.quantum.grover import WGate
    assert _test_args(WGate())

@XFAIL
def test_sympy__physics__quantum__hilbert__ComplexSpace():
    from sympy.physics.quantum.hilbert import ComplexSpace
    assert _test_args(ComplexSpace())

@XFAIL
def test_sympy__physics__quantum__hilbert__DirectSumHilbertSpace():
    from sympy.physics.quantum.hilbert import DirectSumHilbertSpace
    assert _test_args(DirectSumHilbertSpace())

@XFAIL
def test_sympy__physics__quantum__hilbert__FockSpace():
    from sympy.physics.quantum.hilbert import FockSpace
    assert _test_args(FockSpace())

@XFAIL
def test_sympy__physics__quantum__hilbert__HilbertSpace():
    from sympy.physics.quantum.hilbert import HilbertSpace
    assert _test_args(HilbertSpace())

@XFAIL
def test_sympy__physics__quantum__hilbert__L2():
    from sympy.physics.quantum.hilbert import L2
    assert _test_args(L2())

@XFAIL
def test_sympy__physics__quantum__hilbert__TensorPowerHilbertSpace():
    from sympy.physics.quantum.hilbert import TensorPowerHilbertSpace
    assert _test_args(TensorPowerHilbertSpace())

@XFAIL
def test_sympy__physics__quantum__hilbert__TensorProductHilbertSpace():
    from sympy.physics.quantum.hilbert import TensorProductHilbertSpace
    assert _test_args(TensorProductHilbertSpace())

@XFAIL
def test_sympy__physics__quantum__operator__DifferentialOperator():
    from sympy.physics.quantum.operator import DifferentialOperator
    assert _test_args(DifferentialOperator())

@XFAIL
def test_sympy__physics__quantum__operator__HermitianOperator():
    from sympy.physics.quantum.operator import HermitianOperator
    assert _test_args(HermitianOperator())

@XFAIL
def test_sympy__physics__quantum__operator__OuterProduct():
    from sympy.physics.quantum.operator import OuterProduct
    assert _test_args(OuterProduct())

@XFAIL
def test_sympy__physics__quantum__operator__UnitaryOperator():
    from sympy.physics.quantum.operator import UnitaryOperator
    assert _test_args(UnitaryOperator())

@XFAIL
def test_sympy__physics__quantum__piab__PIABBra():
    from sympy.physics.quantum.piab import PIABBra
    assert _test_args(PIABBra())

@XFAIL
def test_sympy__physics__quantum__piab__PIABHamiltonian():
    from sympy.physics.quantum.piab import PIABHamiltonian
    assert _test_args(PIABHamiltonian())

@XFAIL
def test_sympy__physics__quantum__piab__PIABKet():
    from sympy.physics.quantum.piab import PIABKet
    assert _test_args(PIABKet())

@XFAIL
def test_sympy__physics__quantum__qexpr__QExpr():
    from sympy.physics.quantum.qexpr import QExpr
    assert _test_args(QExpr())

@XFAIL
def test_sympy__physics__quantum__qft__Fourier():
    from sympy.physics.quantum.qft import Fourier
    assert _test_args(Fourier())

@XFAIL
def test_sympy__physics__quantum__qft__IQFT():
    from sympy.physics.quantum.qft import IQFT
    assert _test_args(IQFT())

@XFAIL
def test_sympy__physics__quantum__qft__QFT():
    from sympy.physics.quantum.qft import QFT
    assert _test_args(QFT())

@XFAIL
def test_sympy__physics__quantum__qft__RkGate():
    from sympy.physics.quantum.qft import RkGate
    assert _test_args(RkGate())

@XFAIL
def test_sympy__physics__quantum__qubit__IntQubit():
    from sympy.physics.quantum.qubit import IntQubit
    assert _test_args(IntQubit())

@XFAIL
def test_sympy__physics__quantum__qubit__IntQubitBra():
    from sympy.physics.quantum.qubit import IntQubitBra
    assert _test_args(IntQubitBra())

@XFAIL
def test_sympy__physics__quantum__qubit__IntQubitState():
    from sympy.physics.quantum.qubit import IntQubitState
    assert _test_args(IntQubitState())

@XFAIL
def test_sympy__physics__quantum__qubit__Qubit():
    from sympy.physics.quantum.qubit import Qubit
    assert _test_args(Qubit())

@XFAIL
def test_sympy__physics__quantum__qubit__QubitBra():
    from sympy.physics.quantum.qubit import QubitBra
    assert _test_args(QubitBra())

@XFAIL
def test_sympy__physics__quantum__qubit__QubitState():
    from sympy.physics.quantum.qubit import QubitState
    assert _test_args(QubitState())

@XFAIL
def test_sympy__physics__quantum__shor__CMod():
    from sympy.physics.quantum.shor import CMod
    assert _test_args(CMod())

@XFAIL
def test_sympy__physics__quantum__spin__CoupledSpinState():
    from sympy.physics.quantum.spin import CoupledSpinState
    assert _test_args(CoupledSpinState())

@XFAIL
def test_sympy__physics__quantum__spin__J2Op():
    from sympy.physics.quantum.spin import J2Op
    assert _test_args(J2Op())

@XFAIL
def test_sympy__physics__quantum__spin__JminusOp():
    from sympy.physics.quantum.spin import JminusOp
    assert _test_args(JminusOp())

@XFAIL
def test_sympy__physics__quantum__spin__JplusOp():
    from sympy.physics.quantum.spin import JplusOp
    assert _test_args(JplusOp())

@XFAIL
def test_sympy__physics__quantum__spin__JxBra():
    from sympy.physics.quantum.spin import JxBra
    assert _test_args(JxBra())

@XFAIL
def test_sympy__physics__quantum__spin__JxBraCoupled():
    from sympy.physics.quantum.spin import JxBraCoupled
    assert _test_args(JxBraCoupled())

@XFAIL
def test_sympy__physics__quantum__spin__JxKet():
    from sympy.physics.quantum.spin import JxKet
    assert _test_args(JxKet())

@XFAIL
def test_sympy__physics__quantum__spin__JxKetCoupled():
    from sympy.physics.quantum.spin import JxKetCoupled
    assert _test_args(JxKetCoupled())

@XFAIL
def test_sympy__physics__quantum__spin__JxOp():
    from sympy.physics.quantum.spin import JxOp
    assert _test_args(JxOp())

@XFAIL
def test_sympy__physics__quantum__spin__JyBra():
    from sympy.physics.quantum.spin import JyBra
    assert _test_args(JyBra())

@XFAIL
def test_sympy__physics__quantum__spin__JyBraCoupled():
    from sympy.physics.quantum.spin import JyBraCoupled
    assert _test_args(JyBraCoupled())

@XFAIL
def test_sympy__physics__quantum__spin__JyKet():
    from sympy.physics.quantum.spin import JyKet
    assert _test_args(JyKet())

@XFAIL
def test_sympy__physics__quantum__spin__JyKetCoupled():
    from sympy.physics.quantum.spin import JyKetCoupled
    assert _test_args(JyKetCoupled())

@XFAIL
def test_sympy__physics__quantum__spin__JyOp():
    from sympy.physics.quantum.spin import JyOp
    assert _test_args(JyOp())

@XFAIL
def test_sympy__physics__quantum__spin__JzBra():
    from sympy.physics.quantum.spin import JzBra
    assert _test_args(JzBra())

@XFAIL
def test_sympy__physics__quantum__spin__JzBraCoupled():
    from sympy.physics.quantum.spin import JzBraCoupled
    assert _test_args(JzBraCoupled())

@XFAIL
def test_sympy__physics__quantum__spin__JzKet():
    from sympy.physics.quantum.spin import JzKet
    assert _test_args(JzKet())

@XFAIL
def test_sympy__physics__quantum__spin__JzKetCoupled():
    from sympy.physics.quantum.spin import JzKetCoupled
    assert _test_args(JzKetCoupled())

@XFAIL
def test_sympy__physics__quantum__spin__JzOp():
    from sympy.physics.quantum.spin import JzOp
    assert _test_args(JzOp())

@XFAIL
def test_sympy__physics__quantum__spin__Rotation():
    from sympy.physics.quantum.spin import Rotation
    assert _test_args(Rotation())

@XFAIL
def test_sympy__physics__quantum__spin__SpinState():
    from sympy.physics.quantum.spin import SpinState
    assert _test_args(SpinState())

@XFAIL
def test_sympy__physics__quantum__spin__WignerD():
    from sympy.physics.quantum.spin import WignerD
    assert _test_args(WignerD())

@XFAIL
def test_sympy__physics__quantum__state__Bra():
    from sympy.physics.quantum.state import Bra
    assert _test_args(Bra())

@XFAIL
def test_sympy__physics__quantum__state__BraBase():
    from sympy.physics.quantum.state import BraBase
    assert _test_args(BraBase())

@XFAIL
def test_sympy__physics__quantum__state__Ket():
    from sympy.physics.quantum.state import Ket
    assert _test_args(Ket())

@XFAIL
def test_sympy__physics__quantum__state__KetBase():
    from sympy.physics.quantum.state import KetBase
    assert _test_args(KetBase())

@XFAIL
def test_sympy__physics__quantum__state__State():
    from sympy.physics.quantum.state import State
    assert _test_args(State())

@XFAIL
def test_sympy__physics__quantum__state__StateBase():
    from sympy.physics.quantum.state import StateBase
    assert _test_args(StateBase())

@XFAIL
def test_sympy__physics__quantum__state__TimeDepBra():
    from sympy.physics.quantum.state import TimeDepBra
    assert _test_args(TimeDepBra())

@XFAIL
def test_sympy__physics__quantum__state__TimeDepKet():
    from sympy.physics.quantum.state import TimeDepKet
    assert _test_args(TimeDepKet())

@XFAIL
def test_sympy__physics__quantum__state__TimeDepState():
    from sympy.physics.quantum.state import TimeDepState
    assert _test_args(TimeDepState())

@XFAIL
def test_sympy__physics__quantum__state__Wavefunction():
    from sympy.physics.quantum.state import Wavefunction
    assert _test_args(Wavefunction())

@XFAIL
def test_sympy__physics__quantum__tensorproduct__TensorProduct():
    from sympy.physics.quantum.tensorproduct import TensorProduct
    assert _test_args(TensorProduct())

@XFAIL
def test_sympy__physics__secondquant__AnnihilateBoson():
    from sympy.physics.secondquant import AnnihilateBoson
    assert _test_args(AnnihilateBoson())

@XFAIL
def test_sympy__physics__secondquant__AnnihilateFermion():
    from sympy.physics.secondquant import AnnihilateFermion
    assert _test_args(AnnihilateFermion())

@XFAIL
def test_sympy__physics__secondquant__Annihilator():
    from sympy.physics.secondquant import Annihilator
    assert _test_args(Annihilator())

@XFAIL
def test_sympy__physics__secondquant__AntiSymmetricTensor():
    from sympy.physics.secondquant import AntiSymmetricTensor
    assert _test_args(AntiSymmetricTensor())

@XFAIL
def test_sympy__physics__secondquant__BosonState():
    from sympy.physics.secondquant import BosonState
    assert _test_args(BosonState())

@XFAIL
def test_sympy__physics__secondquant__BosonicOperator():
    from sympy.physics.secondquant import BosonicOperator
    assert _test_args(BosonicOperator())

@XFAIL
def test_sympy__physics__secondquant__Commutator():
    from sympy.physics.secondquant import Commutator
    assert _test_args(Commutator())

@XFAIL
def test_sympy__physics__secondquant__CreateBoson():
    from sympy.physics.secondquant import CreateBoson
    assert _test_args(CreateBoson())

@XFAIL
def test_sympy__physics__secondquant__CreateFermion():
    from sympy.physics.secondquant import CreateFermion
    assert _test_args(CreateFermion())

@XFAIL
def test_sympy__physics__secondquant__Creator():
    from sympy.physics.secondquant import Creator
    assert _test_args(Creator())

@XFAIL
def test_sympy__physics__secondquant__Dagger():
    from sympy.physics.secondquant import Dagger
    assert _test_args(Dagger())

@XFAIL
def test_sympy__physics__secondquant__FermionState():
    from sympy.physics.secondquant import FermionState
    assert _test_args(FermionState())

@XFAIL
def test_sympy__physics__secondquant__FermionicOperator():
    from sympy.physics.secondquant import FermionicOperator
    assert _test_args(FermionicOperator())

@XFAIL
def test_sympy__physics__secondquant__FockState():
    from sympy.physics.secondquant import FockState
    assert _test_args(FockState())

@XFAIL
def test_sympy__physics__secondquant__FockStateBosonBra():
    from sympy.physics.secondquant import FockStateBosonBra
    assert _test_args(FockStateBosonBra())

@XFAIL
def test_sympy__physics__secondquant__FockStateBosonKet():
    from sympy.physics.secondquant import FockStateBosonKet
    assert _test_args(FockStateBosonKet())

@XFAIL
def test_sympy__physics__secondquant__FockStateBra():
    from sympy.physics.secondquant import FockStateBra
    assert _test_args(FockStateBra())

@XFAIL
def test_sympy__physics__secondquant__FockStateFermionBra():
    from sympy.physics.secondquant import FockStateFermionBra
    assert _test_args(FockStateFermionBra())

@XFAIL
def test_sympy__physics__secondquant__FockStateFermionKet():
    from sympy.physics.secondquant import FockStateFermionKet
    assert _test_args(FockStateFermionKet())

@XFAIL
def test_sympy__physics__secondquant__FockStateKet():
    from sympy.physics.secondquant import FockStateKet
    assert _test_args(FockStateKet())

@XFAIL
def test_sympy__physics__secondquant__InnerProduct():
    from sympy.physics.secondquant import InnerProduct
    assert _test_args(InnerProduct())

@XFAIL
def test_sympy__physics__secondquant__KroneckerDelta():
    from sympy.physics.secondquant import KroneckerDelta
    assert _test_args(KroneckerDelta())

@XFAIL
def test_sympy__physics__secondquant__NO():
    from sympy.physics.secondquant import NO
    assert _test_args(NO())

@XFAIL
def test_sympy__physics__secondquant__PermutationOperator():
    from sympy.physics.secondquant import PermutationOperator
    assert _test_args(PermutationOperator())

@XFAIL
def test_sympy__physics__secondquant__SqOperator():
    from sympy.physics.secondquant import SqOperator
    assert _test_args(SqOperator())

@XFAIL
def test_sympy__physics__secondquant__TensorSymbol():
    from sympy.physics.secondquant import TensorSymbol
    assert _test_args(TensorSymbol())

@XFAIL
def test_sympy__physics__units__Unit():
    from sympy.physics.units import Unit
    assert _test_args(Unit())

@XFAIL
def test_sympy__polys__numberfields__AlgebraicNumber():
    from sympy.polys.numberfields import AlgebraicNumber
    assert _test_args(AlgebraicNumber())

@XFAIL
def test_sympy__polys__polytools__GroebnerBasis():
    from sympy.polys.polytools import GroebnerBasis
    assert _test_args(GroebnerBasis())

@XFAIL
def test_sympy__polys__polytools__Poly():
    from sympy.polys.polytools import Poly
    assert _test_args(Poly())

@XFAIL
def test_sympy__polys__polytools__PurePoly():
    from sympy.polys.polytools import PurePoly
    assert _test_args(PurePoly())

@XFAIL
def test_sympy__polys__rootoftools__RootOf():
    from sympy.polys.rootoftools import RootOf
    assert _test_args(RootOf())

@XFAIL
def test_sympy__polys__rootoftools__RootSum():
    from sympy.polys.rootoftools import RootSum
    assert _test_args(RootSum())

@XFAIL
def test_sympy__series__limits__Limit():
    from sympy.series.limits import Limit
    assert _test_args(Limit())

@XFAIL
def test_sympy__series__order__Order():
    from sympy.series.order import Order
    assert _test_args(Order())

@XFAIL
def test_sympy__simplify__cse_opts__Sub():
    from sympy.simplify.cse_opts import Sub
    assert _test_args(Sub())

@XFAIL
def test_sympy__tensor__indexed__Idx():
    from sympy.tensor.indexed import Idx
    assert _test_args(Idx())

@XFAIL
def test_sympy__tensor__indexed__Indexed():
    from sympy.tensor.indexed import Indexed
    assert _test_args(Indexed())

@XFAIL
def test_sympy__tensor__indexed__IndexedBase():
    from sympy.tensor.indexed import IndexedBase
    assert _test_args(IndexedBase())
