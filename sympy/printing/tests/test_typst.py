from sympy import MatAdd, MatMul, Array
from sympy.algebras.quaternion import Quaternion
from sympy.calculus.accumulationbounds import AccumBounds
from sympy.combinatorics.permutations import Cycle, Permutation, AppliedPermutation
from sympy.concrete.products import Product
from sympy.concrete.summations import Sum
from sympy.core.containers import Tuple, Dict
from sympy.core.expr import UnevaluatedExpr
from sympy.core.function import (Derivative, Function, Lambda, Subs, diff)
from sympy.core.mod import Mod
from sympy.core.mul import Mul
from sympy.core.numbers import (AlgebraicNumber, Float, I, Integer, Rational, oo, pi)
from sympy.core.parameters import evaluate
from sympy.core.power import Pow
from sympy.core.relational import Eq, Ne
from sympy.core.singleton import S
from sympy.core.symbol import (Symbol, Wild, symbols)
from sympy.functions.combinatorial.factorials import (FallingFactorial, RisingFactorial, binomial, factorial, factorial2, subfactorial)
from sympy.functions.combinatorial.numbers import (bernoulli, bell, catalan, euler, genocchi,
                                                   lucas, fibonacci, tribonacci, divisor_sigma, udivisor_sigma,
                                                   mobius, primenu, primeomega,
                                                   totient, reduced_totient)
from sympy.functions.elementary.complexes import (Abs, arg, conjugate, im, polar_lift, re)
from sympy.functions.elementary.exponential import (LambertW, exp, log)
from sympy.functions.elementary.hyperbolic import (asinh, coth)
from sympy.functions.elementary.integers import (ceiling, floor, frac)
from sympy.functions.elementary.miscellaneous import (Max, Min, root, sqrt)
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import (acsc, asin, cos, cot, sin, tan)
from sympy.functions.special.beta_functions import beta
from sympy.functions.special.delta_functions import (DiracDelta, Heaviside)
from sympy.functions.special.elliptic_integrals import (elliptic_e, elliptic_f, elliptic_k, elliptic_pi)
from sympy.functions.special.error_functions import (Chi, Ci, Ei, Shi, Si, expint)
from sympy.functions.special.gamma_functions import (gamma, uppergamma)
from sympy.functions.special.hyper import (hyper, meijerg)
from sympy.functions.special.mathieu_functions import (mathieuc, mathieucprime, mathieus, mathieusprime)
from sympy.functions.special.polynomials import (assoc_laguerre, assoc_legendre, chebyshevt, chebyshevu, gegenbauer, hermite, jacobi, laguerre, legendre)
from sympy.functions.special.singularity_functions import SingularityFunction
from sympy.functions.special.spherical_harmonics import (Ynm, Znm)
from sympy.functions.special.tensor_functions import (KroneckerDelta, LeviCivita)
from sympy.functions.special.zeta_functions import (dirichlet_eta, lerchphi, polylog, stieltjes, zeta)
from sympy.integrals.integrals import Integral
from sympy.integrals.transforms import (CosineTransform, FourierTransform, InverseCosineTransform, InverseFourierTransform, InverseLaplaceTransform, InverseMellinTransform, InverseSineTransform, LaplaceTransform, MellinTransform, SineTransform)
from sympy.logic import Implies
from sympy.logic.boolalg import (And, Or, Xor, Equivalent, false, Not, true)
from sympy.matrices.dense import Matrix
from sympy.matrices.expressions.kronecker import KroneckerProduct
from sympy.matrices.expressions.matexpr import MatrixSymbol
from sympy.matrices.expressions.permutation import PermutationMatrix
from sympy.matrices.expressions.slice import MatrixSlice
from sympy.matrices.expressions.dotproduct import DotProduct
from sympy.physics.control.lti import TransferFunction, Series, Parallel, Feedback, TransferFunctionMatrix, MIMOSeries, MIMOParallel, MIMOFeedback
from sympy.physics.quantum import Commutator, Operator
from sympy.physics.quantum.trace import Tr
from sympy.physics.units import meter, gibibyte, gram, microgram, second, milli, micro
from sympy.polys.domains.integerring import ZZ
from sympy.polys.fields import field
from sympy.polys.polytools import Poly
from sympy.polys.rings import ring
from sympy.polys.rootoftools import (RootSum, rootof)
from sympy.series.formal import fps
from sympy.series.fourier import fourier_series
from sympy.series.limits import Limit
from sympy.series.order import Order
from sympy.series.sequences import (SeqAdd, SeqFormula, SeqMul, SeqPer)
from sympy.sets.conditionset import ConditionSet
from sympy.sets.contains import Contains
from sympy.sets.fancysets import (ComplexRegion, ImageSet, Range)
from sympy.sets.ordinals import Ordinal, OrdinalOmega, OmegaPower
from sympy.sets.powerset import PowerSet
from sympy.sets.sets import (FiniteSet, Interval, Union, Intersection, Complement, SymmetricDifference, ProductSet)
from sympy.sets.setexpr import SetExpr
from sympy.stats.crv_types import Normal
from sympy.stats.symbolic_probability import (Covariance, Expectation,
                                              Probability, Variance)
from sympy.tensor.array import (ImmutableDenseNDimArray,
                                ImmutableSparseNDimArray,
                                MutableSparseNDimArray,
                                MutableDenseNDimArray,
                                tensorproduct)
from sympy.tensor.array.expressions.array_expressions import ArraySymbol, ArrayElement
from sympy.tensor.indexed import (Idx, Indexed, IndexedBase)
from sympy.tensor.toperators import PartialDerivative
from sympy.vector import CoordSys3D, Cross, Curl, Dot, Divergence, Gradient, Laplacian


from sympy.testing.pytest import (XFAIL, raises, _both_exp_pow,
                                  warns_deprecated_sympy)
from sympy.printing.typst import (typst, translate, greek_letters_set,
                                  typst_greek_dictionary)

import sympy as sym

from sympy.abc import mu, tau


class lowergamma(sym.lowergamma):
    pass   # testing notation inheritance by a subclass with same name

x, y, z, t, w, a, b, c, s, p = symbols('x y z t w a b c s p')
k, m, n = symbols('k m n', integer=True)

def test_printmethod():
    class R(Abs):
        def _typst(self, printer):
            return "foo(%s)" % printer._print(self.args[0])
    assert typst(R(x)) == r"foo(x)"

    class R(Abs):
        def _typst(self, printer):
            return "foo"
    assert typst(R(x)) == r"foo"


def test_typst_basic():
    assert typst(1 + x) == r"x + 1"
    assert typst(x**2) == r"x^(2)"
    assert typst(x**(1 + x)) == r"x^(x + 1)"
    assert typst(x**3 + x + 1 + x**2) == r"x^(3) + x^(2) + x + 1"

    assert typst(2*x*y) == r"2 x y"
    assert typst(2*x*y, mul_symbol='dot') == r"2 dot x dot y"
    assert typst(3*x**2*y, mul_symbol='#h(1cm)') == r"3#h(1cm)x^(2)#h(1cm)y"
    assert typst(1.5*3**x, mul_symbol='#h(1cm)') == r"1.5#h(1cm)3^(x)"

    assert typst(x**S.Half**5) == r"root(32, x)"
    assert typst(Mul(S.Half, x**2, -5, evaluate=False)) == r"1/2 x^(2) (-5)"
    assert typst(Mul(S.Half, x**2, 5, evaluate=False)) == r"1/2 x^(2) 5"
    assert typst(Mul(-5, -5, evaluate=False)) == r"(-5) (-5)"
    assert typst(Mul(5, -5, evaluate=False)) == r"5 (-5)"
    assert typst(Mul(S.Half, -5, S.Half, evaluate=False)) == r"1/2 (-5) 1/2"
    assert typst(Mul(5, I, 5, evaluate=False)) == r"5 i 5"
    assert typst(Mul(5, I, -5, evaluate=False)) == r"5 i (-5)"
    assert typst(Mul(Pow(x, 2), S.Half*x + 1)) == r"x^(2) (x/2 + 1)"
    assert typst(Mul(Pow(x, 3), Rational(2, 3)*x + 1)) == r"x^(3) ((2 x)/3 + 1)"
    assert typst(Mul(Pow(x, 11), 2*x + 1)) == r"x^(11) (2 x + 1)"

    assert typst(Mul(0, 1, evaluate=False)) == r'0 dot 1'
    assert typst(Mul(1, 0, evaluate=False)) == r'1 dot 0'
    assert typst(Mul(1, 1, evaluate=False)) == r'1 dot 1'
    assert typst(Mul(-1, 1, evaluate=False)) == r'(-1) 1'
    assert typst(Mul(1, 1, 1, evaluate=False)) == r'1 dot 1 dot 1'
    assert typst(Mul(1, 2, evaluate=False)) == r'1 dot 2'
    assert typst(Mul(1, S.Half, evaluate=False)) == r'1 dot 1/2'
    assert typst(Mul(1, 1, S.Half, evaluate=False)) == \
        r'1 dot 1 dot 1/2'
    assert typst(Mul(1, 1, 2, 3, x, evaluate=False)) == \
        r'1 dot 1 dot 2 dot 3 x'
    assert typst(Mul(1, -1, evaluate=False)) == r'1 (-1)'
    assert typst(Mul(4, 3, 2, 1, 0, y, x, evaluate=False)) == \
        r'4 dot 3 dot 2 dot 1 dot 0 y x'
    assert typst(Mul(4, 3, 2, 1+z, 0, y, x, evaluate=False)) == \
        r'4 dot 3 dot 2 (z + 1) 0 y x'
    assert typst(Mul(Rational(2, 3), Rational(5, 7), evaluate=False)) == \
        r'2/3 dot 5/7'

    assert typst(1/x) == r"1/x"
    assert typst(1/x, fold_short_frac=True) == r"1/x"
    assert typst(-S(3)/2) == r"- 3/2"
    assert typst(-S(3)/2, fold_short_frac=True) == r"- 3/2"
    assert typst(1/x**2) == r"1/x^(2)"
    assert typst(1/(x + y)/2) == r"1/(2 (x + y))"
    assert typst(x/2) == r"x/2"
    assert typst(x/2, fold_short_frac=True) == r"x/2"
    assert typst((x + y)/(2*x)) == r"(x + y)/(2 x)"
    assert typst((x + y)/(2*x), fold_short_frac=True) == \
        r"(x + y)/(2 x)"
    assert typst((x + y)/(2*x), long_frac_ratio=0) == \
        r"1/(2 x) (x + y)"
    assert typst((x + y)/x) == r"(x + y)/x"
    assert typst((x + y)/x, long_frac_ratio=3) == r"(x + y)/x"
    assert typst((2*sqrt(2)*x)/3) == r"(2 sqrt(2) x)/3"
    assert typst((2*sqrt(2)*x)/3, long_frac_ratio=2) == \
        r"(2 x)/3 sqrt(2)"
    assert typst(binomial(x, y)) == r'binom(x, y)'

    x_star = Symbol('x^*')
    f = Function('f')
    assert typst(x_star**2) == r"(x^(*))^(2)"
    assert typst(x_star**2, parenthesize_super=False) == r"x^*^(2)"
    assert typst(Derivative(f(x_star), x_star,2)) == r'd^(2) / (d (x^(*))^(2)) f(x^(*))'
    assert typst(Derivative(f(x_star), x_star,2), parenthesize_super=False) == r'd^(2) / (d x^*^(2)) f(x^*)'

    assert typst(2*Integral(x, x)/3) == r"(2 integral x d x)/3"

    assert typst(sqrt(x)) == r"sqrt(x)"
    assert typst(x**Rational(1, 3)) == r"root(3, x)"
    assert typst(x**Rational(1, 3), root_notation=False) == r"x^(1/3)"
    assert typst(sqrt(x)**3) == r"x^(3/2)"
    assert typst(sqrt(x)) == r"sqrt(x)"
    assert typst(x**Rational(1, 3)) == r"root(3, x)"
    assert typst(sqrt(x)**3) == r"x^(3/2)"
    assert typst(x**Rational(3, 4)) == r"x^(3/4)"
    assert typst(x**Rational(3, 4), fold_frac_powers=True) == r"x^(3/4)"
    assert typst((x + 1)**Rational(3, 4)) == \
        r"(x + 1)^(3/4)"
    assert typst((x + 1)**Rational(3, 4), fold_frac_powers=True) == \
        r"(x + 1)^(3/4)"
    assert typst(AlgebraicNumber(sqrt(2))) == r"sqrt(2)"
    assert typst(AlgebraicNumber(sqrt(2), [3, -7])) == r"-7 + 3 sqrt(2)"
    assert typst(AlgebraicNumber(sqrt(2), alias='alpha')) == r"alpha"
    assert typst(AlgebraicNumber(sqrt(2), [3, -7], alias='alpha')) == \
        r"3 alpha - 7"
    assert typst(AlgebraicNumber(2**(S(1)/3), [1, 3, -7], alias='beta')) == \
        r"beta^(2) + 3 beta - 7"

    k = ZZ.cyclotomic_field(5)
    assert typst(k.ext.field_element([1, 2, 3, 4])) == \
        r"zeta^(3) + 2 zeta^(2) + 3 zeta + 4"
    assert typst(k.ext.field_element([1, 2, 3, 4]), order='old') == \
        r"4 + 3 zeta + 2 zeta^(2) + zeta^(3)"
    assert typst(k.primes_above(19)[0]) == \
        r"(19, zeta^(2) + 5 zeta + 1)"
    assert typst(k.primes_above(19)[0], order='old') == \
           r"(19, 1 + 5 zeta + zeta^(2))"
    assert typst(k.primes_above(7)[0]) == r"(7)"

    assert typst(1.5e20*x) == r"1.5 dot 10^(20) x"
    assert typst(1.5e20*x, mul_symbol='dot') == r"1.5 dot 10^(20) dot x"
    assert typst(1.5e20*x, mul_symbol='times') == \
        r"1.5 times 10^(20) times x"

    assert typst(1/sin(x)) == r"1/(sin(x))"
    assert typst(sin(x)**-1) == r"1/(sin(x))"
    assert typst(sin(x)**Rational(3, 2)) == \
        r"sin^(3/2)(x)"
    assert typst(sin(x)**Rational(3, 2), fold_frac_powers=True) == \
        r"sin^(3/2)(x)"

    assert typst(~x) == r"not x"
    assert typst(x & y) == r"x and y"
    assert typst(x & y & z) == r"x and y and z"
    assert typst(x | y) == r"x or y"
    assert typst(x | y | z) == r"x or y or z"
    assert typst((x & y) | z) == r"z or (x and y)"
    assert typst(Implies(x, y)) == r"x => y"
    assert typst(~(x >> ~y)) == r"x arrow.r.double.not not y"
    assert typst(Implies(Or(x,y), z)) == r"(x or y) => z"
    assert typst(Implies(z, Or(x,y))) == r"z => (x or y)"
    assert typst(~(x & y)) == r"not (x and y)"

    assert typst(~x, symbol_names={x: "x_i"}) == r"not x_i"
    assert typst(x & y, symbol_names={x: "x_i", y: "y_i"}) == \
        r"x_i and y_i"
    assert typst(x & y & z, symbol_names={x: "x_i", y: "y_i", z: "z_i"}) == \
        r"x_i and y_i and z_i"
    assert typst(x | y, symbol_names={x: "x_i", y: "y_i"}) == r"x_i or y_i"
    assert typst(x | y | z, symbol_names={x: "x_i", y: "y_i", z: "z_i"}) == \
        r"x_i or y_i or z_i"
    assert typst((x & y) | z, symbol_names={x: "x_i", y: "y_i", z: "z_i"}) == \
        r"z_i or (x_i and y_i)"
    assert typst(Implies(x, y), symbol_names={x: "x_i", y: "y_i"}) == \
        r"x_i => y_i"
    assert typst(Pow(Rational(1, 3), -1, evaluate=False)) == r"1/(1/3)"
    assert typst(Pow(Rational(1, 3), -2, evaluate=False)) == r"1/(1/3)^(2)"
    assert typst(Pow(Integer(1)/100, -1, evaluate=False)) == r"1/(1/100)"

    p = Symbol('p', positive=True)
    assert typst(exp(-p)*log(p)) == r"e^(- p) log(p)"

    assert typst(Pow(Rational(2, 3), -1, evaluate=False)) == r'1/(2/3)'
    assert typst(Pow(Rational(4, 3), -1, evaluate=False)) == r'1/(4/3)'
    assert typst(Pow(Rational(-3, 4), -1, evaluate=False)) == r'1/(- 3/4)'
    assert typst(Pow(Rational(-4, 4), -1, evaluate=False)) == r'1/(-1)'
    assert typst(Pow(Rational(1, 3), -1, evaluate=False)) == r'1/(1/3)'
    assert typst(Pow(Rational(-1, 3), -1, evaluate=False)) == r'1/(- 1/3)'

def test_typst_builtins():
    assert typst(True) == r'upright("True")'
    assert typst(False) == r'upright("False")'
    assert typst(None) == r'upright("None")'
    assert typst(true) == r'upright("True")'
    assert typst(false) == r'upright("False")'


def test_typst_SingularityFunction():
    assert typst(SingularityFunction(x, 4, 5)) == \
        r"angle.l x - 4 angle.r ^ (5)"
    assert typst(SingularityFunction(x, -3, 4)) == \
        r"angle.l x + 3 angle.r ^ (4)"
    assert typst(SingularityFunction(x, 0, 4)) == \
        r"angle.l x angle.r ^ (4)"
    assert typst(SingularityFunction(x, a, n)) == \
        r"angle.l - a + x angle.r ^ (n)"
    assert typst(SingularityFunction(x, 4, -2)) == \
        r"angle.l x - 4 angle.r ^ (-2)"
    assert typst(SingularityFunction(x, 4, -1)) == \
        r"angle.l x - 4 angle.r ^ (-1)"

    assert typst(SingularityFunction(x, 4, 5)**3) == \
        r"(angle.l x - 4 angle.r ^ (5))^(3)"
    assert typst(SingularityFunction(x, -3, 4)**3) == \
        r"(angle.l x + 3 angle.r ^ (4))^(3)"
    assert typst(SingularityFunction(x, 0, 4)**3) == \
        r"(angle.l x angle.r ^ (4))^(3)"
    assert typst(SingularityFunction(x, a, n)**3) == \
        r"(angle.l - a + x angle.r ^ (n))^(3)"
    assert typst(SingularityFunction(x, 4, -2)**3) == \
        r"(angle.l x - 4 angle.r ^ (-2))^(3)"
    assert typst((SingularityFunction(x, 4, -1)**3)**3) == \
        r"(angle.l x - 4 angle.r ^ (-1))^(9)"

def test_typst_cycle():
    assert typst(Cycle(1, 2, 4)) == r"(1 space 2 space 4)"
    assert typst(Cycle(1, 2)(4, 5, 6)) == \
        r"(1 space 2)(4 space 5 space 6)"
    assert typst(Cycle()) == r"()"


def test_typst_permutation():
    assert typst(Permutation(1, 2, 4)) == r"(1 space 2 space 4)"
    assert typst(Permutation(1, 2)(4, 5, 6)) == \
        r"(1 space 2)(4 space 5 space 6)"
    assert typst(Permutation()) == r"()"
    assert typst(Permutation(2, 4)*Permutation(5)) == \
        r"(2 space 4)(5)"
    assert typst(Permutation(5)) == r"(5)"

    assert typst(Permutation(0, 1), perm_cyclic=False) == \
        r"mat(0, 1; 1, 0)"
    assert typst(Permutation(0, 1)(2, 3), perm_cyclic=False) == \
        r"mat(0, 1, 2, 3; 1, 0, 3, 2)"
    assert typst(Permutation(), perm_cyclic=False) == \
        r"()"

    with warns_deprecated_sympy():
        old_print_cyclic = Permutation.print_cyclic
        Permutation.print_cyclic = False
        assert typst(Permutation(0, 1)(2, 3)) == \
            r"mat(0, 1, 2, 3; 1, 0, 3, 2)"
        Permutation.print_cyclic = old_print_cyclic

def test_typst_Float():
    assert typst(Float(1.0e100)) == r"1.0 dot 10^(100)"
    assert typst(Float(1.0e-100)) == r"1.0 dot 10^(-100)"
    assert typst(Float(1.0e-100), mul_symbol="times") == \
        r"1.0 times 10^(-100)"
    assert typst(Float('10000.0'), full_prec=False, min=-2, max=2) == \
        r"1.0 dot 10^(4)"
    assert typst(Float('10000.0'), full_prec=False, min=-2, max=4) == \
        r"1.0 dot 10^(4)"
    assert typst(Float('10000.0'), full_prec=False, min=-2, max=5) == \
        r"10000.0"
    assert typst(Float('0.099999'), full_prec=True,  min=-2, max=5) == \
        r"9.99990000000000 dot 10^(-2)"

def test_typst_vector_expressions():
    A = CoordSys3D('A')

    assert typst(Cross(A.i, A.j*A.x*3+A.k)) == \
        r"bold(hat(i)_(A)) times ((3 bold(x_(A)))bold(hat(j)_(A)) + bold(hat(k)_(A)))"
    assert typst(Cross(A.i, A.j)) == \
        r"bold(hat(i)_(A)) times bold(hat(j)_(A))"
    assert typst(x*Cross(A.i, A.j)) == \
        r"x (bold(hat(i)_(A)) times bold(hat(j)_(A)))"
    assert typst(Cross(x*A.i, A.j)) == \
        r'- bold(hat(j)_(A)) times ((x)bold(hat(i)_(A)))'

    assert typst(Curl(3*A.x*A.j)) == \
        r"nabla times ((3 bold(x_(A)))bold(hat(j)_(A)))"
    assert typst(Curl(3*A.x*A.j+A.i)) == \
        r"nabla times (bold(hat(i)_(A)) + (3 bold(x_(A)))bold(hat(j)_(A)))"
    assert typst(Curl(3*x*A.x*A.j)) == \
        r"nabla times ((3 bold(x_(A)) x)bold(hat(j)_(A)))"
    assert typst(x*Curl(3*A.x*A.j)) == \
        r"x (nabla times ((3 bold(x_(A)))bold(hat(j)_(A))))"

    assert typst(Divergence(3*A.x*A.j+A.i)) == \
        r"nabla dot (bold(hat(i)_(A)) + (3 bold(x_(A)))bold(hat(j)_(A)))"
    assert typst(Divergence(3*A.x*A.j)) == \
        r"nabla dot ((3 bold(x_(A)))bold(hat(j)_(A)))"
    assert typst(x*Divergence(3*A.x*A.j)) == \
        r"x (nabla dot ((3 bold(x_(A)))bold(hat(j)_(A))))"

    assert typst(Dot(A.i, A.j*A.x*3+A.k)) == \
        r"bold(hat(i)_(A)) dot ((3 bold(x_(A)))bold(hat(j)_(A)) + bold(hat(k)_(A)))"
    assert typst(Dot(A.i, A.j)) == \
        r"bold(hat(i)_(A)) dot bold(hat(j)_(A))"
    assert typst(Dot(x*A.i, A.j)) == \
        r"bold(hat(j)_(A)) dot ((x)bold(hat(i)_(A)))"
    assert typst(x*Dot(A.i, A.j)) == \
        r"x (bold(hat(i)_(A)) dot bold(hat(j)_(A)))"

    assert typst(Gradient(A.x)) == r"nabla bold(x_(A))"
    assert typst(Gradient(A.x + 3*A.y)) == \
        r"nabla (bold(x_(A)) + 3 bold(y_(A)))"
    assert typst(x*Gradient(A.x)) == r"x (nabla bold(x_(A)))"
    assert typst(Gradient(x*A.x)) == r"nabla (bold(x_(A)) x)"

    assert typst(Laplacian(A.x)) == r"Delta bold(x_(A))"
    assert typst(Laplacian(A.x + 3*A.y)) == \
        r"Delta (bold(x_(A)) + 3 bold(y_(A)))"
    assert typst(x*Laplacian(A.x)) == r"x (Delta bold(x_(A)))"
    assert typst(Laplacian(x*A.x)) == r"Delta (bold(x_(A)) x)"

def test_AppliedPermutation():
    p = Permutation(0, 1, 2)
    x = Symbol('x')
    assert typst(AppliedPermutation(p, x)) == \
        r'sigma_((0 space 1 space 2))(x)'

def test_imaginary_unit():
    assert typst(1 + I) == r'1 + i'
    assert typst(1 + I, imaginary_unit='i') == r'1 + i'
    assert typst(1 + I, imaginary_unit='j') == r'1 + j'
    assert typst(1 + I, imaginary_unit='foo') == r'1 + foo'
    assert typst(I, imaginary_unit="ti") == r'text(i)'
    assert typst(I, imaginary_unit="tj") == r'text(j)'


def test_text_re_im():
    assert typst(im(x), gothic_re_im=True) == r'Im(x)'
    assert typst(im(x), gothic_re_im=False) == r'upright("im")(x)'
    assert typst(re(x), gothic_re_im=True) == r'Re(x)'
    assert typst(re(x), gothic_re_im=False) == r'upright("re")(x)'

def test_typst_limits():
    assert typst(Limit(x, x, oo)) == r"lim_(x -> infinity) x"

    # issue 8175
    f = Function('f')
    assert typst(Limit(f(x), x, 0)) == r"lim_(x -> 0^+) f(x)"
    assert typst(Limit(f(x), x, 0, "-")) == \
        r"lim_(x -> 0^-) f(x)"

    # issue #10806
    assert typst(Limit(f(x), x, 0)**2) == \
        r"(lim_(x -> 0^+) f(x))^(2)"
    # bi-directional limit
    assert typst(Limit(f(x), x, 0, dir='+-')) == \
        r"lim_(x -> 0) f(x)"


def test_typst_log():
    assert typst(log(x)) == r"log(x)"
    assert typst(log(x), ln_notation=True) == r"ln(x)"
    assert typst(log(x) + log(y)) == \
        r"log(x) + log(y)"
    assert typst(log(x) + log(y), ln_notation=True) == \
        r"ln(x) + ln(y)"
    assert typst(pow(log(x), x)) == r"log(x)^(x)"
    assert typst(pow(log(x), x), ln_notation=True) == \
        r"ln(x)^(x)"


def test_typst_sum():
    assert typst(Sum(x*y**2, (x, -2, 2), (y, -5, 5))) == \
        r"sum_(-2 <= x <= 2 \ -5 <= y <= 5) x y^(2)"
    assert typst(Sum(x**2, (x, -2, 2))) == \
        r"sum_(x=-2)^(2) x^(2)"
    assert typst(Sum(x**2 + y, (x, -2, 2))) == \
        r"sum_(x=-2)^(2) (x^(2) + y)"
    assert typst(Sum(x**2 + y, (x, -2, 2))**2) == \
        r"(sum_(x=-2)^(2) (x^(2) + y))^(2)"

def test_typst_product():
    assert typst(Product(x*y**2, (x, -2, 2), (y, -5, 5))) == \
        r"product_(-2 <= x <= 2 \ -5 <= y <= 5) x y^(2)"
    assert typst(Product(x**2, (x, -2, 2))) == \
        r"product_(x=-2)^(2) x^(2)"
    assert typst(Product(x**2 + y, (x, -2, 2))) == \
        r"product_(x=-2)^(2) (x^(2) + y)"

    assert typst(Product(x, (x, -2, 2))**2) == \
        r"(product_(x=-2)^(2) x)^(2)"

def test_typst_symbols():
    Gamma, lmbda, rho = symbols('Gamma, lambda, rho')
    tau, Tau, TAU, taU = symbols('tau, Tau, TAU, taU')
    assert typst(tau) == r"tau"
    assert typst(Tau) == r"Tau"
    assert typst(TAU) == r"tau"
    assert typst(taU) == r"tau"
    # Check that all capitalized greek letters are handled explicitly
    capitalized_letters = {l.capitalize() for l in greek_letters_set}
    assert len(capitalized_letters - set(typst_greek_dictionary.keys())) == 0
    assert typst(Gamma + lmbda) == r"Gamma + lambda"
    assert typst(Gamma * lmbda) == r"Gamma lambda"
    assert typst(Symbol('q1')) == r"q_(1)"
    assert typst(Symbol('q21')) == r"q_(21)"
    assert typst(Symbol('epsilon0')) == r"epsilon.alt_(0)"
    assert typst(Symbol('omega1')) == r"omega_(1)"
    assert typst(Symbol('91')) == r'"91"'
    assert typst(Symbol('alpha_new')) == r"alpha_(new)"
    assert typst(Symbol('C^orig')) == r"C^(orig)"
    assert typst(Symbol('x^alpha')) == r"x^(alpha)"
    assert typst(Symbol('beta^alpha')) == r"beta^(alpha)"
    assert typst(Symbol('e^Alpha')) == r"e^(Alpha)"
    assert typst(Symbol('omega_alpha^beta')) == r"omega^(beta)_(alpha)"
    assert typst(Symbol('omega') ** Symbol('beta')) == r"omega^(beta)"

@_both_exp_pow
def test_typst_functions():
    assert typst(exp(x)) == r"e^(x)"
    assert typst(exp(1) + exp(2)) == r"e + e^(2)"

    f = Function('f')
    assert typst(f(x)) == r'f(x)'
    assert typst(f) == r'f'

    g = Function('g')
    assert typst(g(x, y)) == r'g(x, y)'
    assert typst(g) == r'g'

    h = Function('h')
    assert typst(h(x, y, z)) == r'h(x, y, z)'
    assert typst(h) == r'h'

    Li = Function('Li')
    assert typst(Li) == r'upright("Li")'
    assert typst(Li(x)) == r'upright("Li")(x)'

    mybeta = Function('beta')
    # not to be confused with the beta function
    assert typst(mybeta(x, y, z)) == r'beta(x, y, z)'
    assert typst(beta(x, y)) == r'Beta(x, y)'
    assert typst(beta(x, evaluate=False)) == r'Beta(x, x)'
    assert typst(beta(x, y)**2) == r'Beta^(2)(x, y)'
    assert typst(mybeta(x)) == r'beta(x)'
    assert typst(mybeta) == r'beta'

    g = Function('gamma')
    # not to be confused with the gamma function
    assert typst(g(x, y, z)) == r'gamma(x, y, z)'
    assert typst(g(x)) == r'gamma(x)'
    assert typst(g) == r'gamma'

    a_1 = Function('a_1')
    assert typst(a_1) == r'a_(1)'
    assert typst(a_1(x)) == r'a_(1)(x)'
    assert typst(Function('a_1')) == r'a_(1)'

    assert typst(Abs(x)) == r"abs(x)"
    assert typst(Abs(x)**2) == r"abs(x)^(2)"

    assert typst(floor(x)) == r"floor(x)"
    assert typst(ceiling(x)) == r"ceil(x)"
    assert typst(floor(x)**2) == r"floor(x)^(2)"
    assert typst(ceiling(x)**2) == r"ceil(x)^(2)"

    assert typst(conjugate(x)) == r"overline(x)"
    assert typst(conjugate(x)**2) == r"overline(x)^(2)"
    assert typst(conjugate(x**2)) == r"overline(x)^(2)"

    assert typst(arg(x)) == r'arg(x)'

    # Test latex printing of function names with "_"
    assert typst(polar_lift(0)) == \
        r'upright("polar_lift")(0)'
    assert typst(polar_lift(0)**3) == \
        r'upright("polar_lift")^(3)(0)'
    
    assert typst(Function('ab')) == r'upright("ab")'
    assert typst(Function('ab1')) == r'upright("ab")_(1)'
    assert typst(Function('ab12')) == r'upright("ab")_(12)'
    assert typst(Function('ab_1')) == r'upright("ab")_(1)'
    assert typst(Function('ab_12')) == r'upright("ab")_(12)'
    assert typst(Function('ab_c')) == r'upright("ab")_(c)'
    assert typst(Function('ab_cd')) == r'upright("ab")_(cd)'
    # > with argument
    assert typst(Function('ab')(Symbol('x'))) == r'upright("ab")(x)'
    assert typst(Function('ab1')(Symbol('x'))) == r'upright("ab")_(1)(x)'
    assert typst(Function('ab12')(Symbol('x'))) == r'upright("ab")_(12)(x)'
    assert typst(Function('ab_1')(Symbol('x'))) == r'upright("ab")_(1)(x)'
    assert typst(Function('ab_c')(Symbol('x'))) == r'upright("ab")_(c)(x)'
    assert typst(Function('ab_cd')(Symbol('x'))) == r'upright("ab")_(cd)(x)'

    # > with power
    #   does not work on functions without brackets

    # > with argument and power combined
    assert typst(Function('a')()**2) == r'a^(2)()'
    assert typst(Function('a1')()**2) == r'a_(1)^(2)()'
    assert typst(Function('a12')()**2) == r'a_(12)^(2)()'
    assert typst(Function('a_1')()**2) == r'a_(1)^(2)()'
    assert typst(Function('a_12')()**2) == r'a_(12)^(2)()'
    assert typst(Function('a')(Symbol('x'))**2) == r'a^(2)(x)'
    assert typst(Function('a1')(Symbol('x'))**2) == r'a_(1)^(2)(x)'
    assert typst(Function('a12')(Symbol('x'))**2) == r'a_(12)^(2)(x)'
    assert typst(Function('a_1')(Symbol('x'))**2) == r'a_(1)^(2)(x)'
    assert typst(Function('a_12')(Symbol('x'))**2) == r'a_(12)^(2)(x)'

    assert typst(Function('a')()**32) == r'a^(32)()'
    assert typst(Function('a1')()**32) == r'a_(1)^(32)()'
    assert typst(Function('a12')()**32) == r'a_(12)^(32)()'
    assert typst(Function('a_1')()**32) == r'a_(1)^(32)()'
    assert typst(Function('a_12')()**32) == r'a_(12)^(32)()'
    assert typst(Function('a')(Symbol('x'))**32) == r'a^(32)(x)'
    assert typst(Function('a1')(Symbol('x'))**32) == r'a_(1)^(32)(x)'
    assert typst(Function('a12')(Symbol('x'))**32) == r'a_(12)^(32)(x)'
    assert typst(Function('a_1')(Symbol('x'))**32) == r'a_(1)^(32)(x)'
    assert typst(Function('a_12')(Symbol('x'))**32) == r'a_(12)^(32)(x)'

    assert typst(Function('a')()**a) == r'a^(a)()'
    assert typst(Function('a1')()**a) == r'a_(1)^(a)()'
    assert typst(Function('a12')()**a) == r'a_(12)^(a)()'
    assert typst(Function('a_1')()**a) == r'a_(1)^(a)()'
    assert typst(Function('a_12')()**a) == r'a_(12)^(a)()'
    assert typst(Function('a')(Symbol('x'))**a) == r'a^(a)(x)'
    assert typst(Function('a1')(Symbol('x'))**a) == r'a_(1)^(a)(x)'
    assert typst(Function('a12')(Symbol('x'))**a) == r'a_(12)^(a)(x)'
    assert typst(Function('a_1')(Symbol('x'))**a) == r'a_(1)^(a)(x)'
    assert typst(Function('a_12')(Symbol('x'))**a) == r'a_(12)^(a)(x)'

    ab = Symbol('ab')
    assert typst(Function('a')()**ab) == r'a^("ab")()'
    assert typst(Function('a1')()**ab) == r'a_(1)^("ab")()'
    assert typst(Function('a12')()**ab) == r'a_(12)^("ab")()'
    assert typst(Function('a_1')()**ab) == r'a_(1)^("ab")()'
    assert typst(Function('a_12')()**ab) == r'a_(12)^("ab")()'
    assert typst(Function('a')(Symbol('x'))**ab) == r'a^("ab")(x)'
    assert typst(Function('a1')(Symbol('x'))**ab) == r'a_(1)^("ab")(x)'
    assert typst(Function('a12')(Symbol('x'))**ab) == r'a_(12)^("ab")(x)'
    assert typst(Function('a_1')(Symbol('x'))**ab) == r'a_(1)^("ab")(x)'
    assert typst(Function('a_12')(Symbol('x'))**ab) == r'a_(12)^("ab")(x)'

    assert typst(Function('a^12')(x)) == r'a^(12)(x)'
    assert typst(Function('a^12')(x) ** ab) == r'(a^(12))^("ab")(x)'
    assert typst(Function('a__12')(x)) == r'a^(12)(x)'
    assert typst(Function('a__12')(x) ** ab) == r'(a^(12))^("ab")(x)'
    assert typst(Function('a_1__1_2')(x)) == r'a^(1)_(1 2)(x)'

    # issue 5868
    omega1 = Function('omega1')
    assert typst(omega1) == r"omega_(1)"
    assert typst(omega1(x)) == r"omega_(1)(x)"

    assert typst(sin(x)) == r"sin(x)"
    assert typst(sin(x), fold_func_brackets=True) == r"sin x"
    assert typst(sin(2*x**2), fold_func_brackets=True) == \
        r"sin 2 x^(2)"
    assert typst(sin(x**2), fold_func_brackets=True) == \
        r"sin x^(2)"

    assert typst(asin(x)**2) == r'upright("asin")^(2)(x)'
    assert typst(asin(x)**2, inv_trig_style="full") == \
        r"arcsin^(2)(x)"
    assert typst(asin(x)**2, inv_trig_style="power") == \
        r"sin^(-1)(x)^(2)"
    assert typst(asin(x**2), inv_trig_style="power",
                 fold_func_brackets=True) == \
        r"sin^(-1) x^(2)"
    assert typst(acsc(x), inv_trig_style="full") == \
        r'upright("arccsc")(x)'
    assert typst(asinh(x), inv_trig_style="full") == \
        r'upright("arsinh")(x)'

    assert typst(factorial(k)) == r"k!"
    assert typst(factorial(-k)) == r"(- k)!"
    assert typst(factorial(k)**2) == r"k!^(2)"

    assert typst(subfactorial(k)) == r"!k"
    assert typst(subfactorial(-k)) == r"!(- k)"
    assert typst(subfactorial(k)**2) == r"(!k)^(2)"

    assert typst(factorial2(k)) == r"k!!"
    assert typst(factorial2(-k)) == r"(- k)!!"
    assert typst(factorial2(k)**2) == r"k!!^(2)"

    assert typst(binomial(2, k)) == r"binom(2, k)"
    assert typst(binomial(2, k)**2) == r"binom(2, k)^(2)"

    assert typst(FallingFactorial(3, k)) == r"(3)_(k)"
    assert typst(RisingFactorial(3, k)) == r"3^((k))"

    assert typst(floor(x)) == r"floor(x)"
    assert typst(ceiling(x)) == r"ceil(x)"
    assert typst(frac(x)) == r'upright("frac")(x)'
    assert typst(floor(x)**2) == r"floor(x)^(2)"
    assert typst(ceiling(x)**2) == r"ceil(x)^(2)"
    assert typst(frac(x)**2) == r'upright("frac")^(2)(x)'

    assert typst(Min(x, 2, x**3)) == r"min(2, x, x^(3))"
    assert typst(Min(x, y)**2) == r"min(x, y)^(2)"
    assert typst(Max(x, 2, x**3)) == r"max(2, x, x^(3))"
    assert typst(Max(x, y)**2) == r"max(x, y)^(2)"
    assert typst(Abs(x)) == r"abs(x)"
    assert typst(Abs(x)**2) == r"abs(x)^(2)"
    assert typst(re(x)) == r'upright("re")(x)'
    assert typst(re(x + y)) == \
        r'upright("re")(x) + upright("re")(y)'
    assert typst(im(x)) == r'upright("im")(x)'
    assert typst(conjugate(x)) == r"overline(x)"
    assert typst(conjugate(x)**2) == r"overline(x)^(2)"
    assert typst(conjugate(x**2)) == r"overline(x)^(2)"
    assert typst(gamma(x)) == r"gamma(x)"
    w = Wild('w')
    assert typst(gamma(w)) == r"gamma(w)"
    assert typst(Order(x)) == r"O(x)"
    assert typst(Order(x, x)) == r"O(x)"
    assert typst(Order(x, (x, 0))) == r"O(x)"
    assert typst(Order(x, (x, oo))) == r"O(x; x->infinity)"
    assert typst(Order(x - y, (x, y))) == \
        r"O(x - y; x->y)"
    assert typst(Order(x, x, y)) == \
        r"O(x; (x, #h(0.5em)y)->(0, #h(0.5em)0))"
    assert typst(Order(x, x, y)) == \
        r"O(x; (x, #h(0.5em)y)->(0, #h(0.5em)0))"
    assert typst(Order(x, (x, oo), (y, oo))) == \
        r"O(x; (x, #h(0.5em)y)->(infinity, #h(0.5em)infinity))"
    assert typst(lowergamma(x, y)) == r'gamma(x, y)'
    assert typst(lowergamma(x, y)**2) == r'gamma^(2)(x, y)'
    assert typst(uppergamma(x, y)) == r'Gamma(x, y)'
    assert typst(uppergamma(x, y)**2) == r'Gamma^(2)(x, y)'

    assert typst(cot(x)) == r'cot(x)'
    assert typst(coth(x)) == r'coth(x)'
    assert typst(re(x)) == r'upright("re")(x)'
    assert typst(im(x)) == r'upright("im")(x)'
    assert typst(root(x, y)) == r'x^(1/y)'
    assert typst(arg(x)) == r'arg(x)'

    # assert typst(zeta(x)) == r'zeta(x)'
    # assert typst(zeta(x)**2) == r'zeta^(2)(x)'
    # assert typst(zeta(x, y)) == r'zeta(x, y)'
    # assert typst(zeta(x, y)**2) == r'zeta^(2)(x, y)'
    # assert typst(dirichlet_eta(x)) == r'upright("dirichlet_eta")(x)'
    # assert typst(dirichlet_eta(x)**2) == r'upright("dirichlet_eta")^(2)(x)'
    # assert typst(polylog(x, y)) == r'polylog(x, y)'
    # assert typst(polylog(x, y)**2) == r'polylog^(2)(x, y)'
    # assert typst(lerchphi(x, y, n)) == r'lerchphi(x, y, n)'
    # assert typst(lerchphi(x, y, n)**2) == r'lerchphi^(2)(x, y, n)'
    # assert typst(stieltjes(x)) == r'stieltjes(x)'
    # assert typst(stieltjes(x)**2) == r'stieltjes^(2)(x)'
    # assert typst(stieltjes(x, y)) == r'stieltjes(x, y)'
    # assert typst(stieltjes(x, y)**2) == r'stieltjes^(2)(x, y)'

    # assert typst(elliptic_k(z)) == r'elliptic_k(z)'
    # assert typst(elliptic_k(z)**2) == r'elliptic_k^(2)(z)'
    # assert typst(elliptic_f(x, y)) == r'elliptic_f(x, y)'
    # assert typst(elliptic_f(x, y)**2) == r'elliptic_f^(2)(x, y)'
    # assert typst(elliptic_e(x, y)) == r'elliptic_e(x, y)'
    # assert typst(elliptic_e(x, y)**2) == r'elliptic_e^(2)(x, y)'
    # assert typst(elliptic_e(z)) == r'elliptic_e(z)'
    # assert typst(elliptic_e(z)**2) == r'elliptic_e^(2)(z)'
    # assert typst(elliptic_pi(x, y, z)) == r'elliptic_pi(x, y, z)'
    # assert typst(elliptic_pi(x, y, z)**2) == r'elliptic_pi^(2)(x, y, z)'
    # assert typst(elliptic_pi(x, y)) == r'elliptic_pi(x, y)'
    # assert typst(elliptic_pi(x, y)**2) == r'elliptic_pi^(2)(x, y)'

    # assert typst(Ei(x)) == r'upright("Ei")(x)'
    # assert typst(Ei(x)**2) == r'upright("Ei")^(2)(x)'
    # assert typst(expint(x, y)) == r'expint(x, y)'
    # assert typst(expint(x, y)**2) == r'expint^(2)(x, y)'
    # assert typst(Shi(x)**2) == r'upright("Shi")^(2)(x)'
    # assert typst(Si(x)**2) == r'upright("Si")^(2)(x)'
    # assert typst(Ci(x)**2) == r'upright("Ci")^(2)(x)'
    # assert typst(Chi(x)**2) == r'Chi^(2)(x)'
    # assert typst(Chi(x)) == r'Chi(x)'
    # assert typst(jacobi(n, a, b, x)) == r'jacobi(n, a, b, x)'
    # assert typst(jacobi(n, a, b, x)**2) == r'jacobi^(2)(n, a, b, x)'
    # assert typst(gegenbauer(n, a, x)) == r'gegenbauer(n, a, x)'
    # assert typst(gegenbauer(n, a, x)**2) == r'gegenbauer^(2)(n, a, x)'
    # assert typst(chebyshevt(n, x)) == r'chebyshevt(n, x)'
    # assert typst(chebyshevt(n, x)**2) == r'chebyshevt^(2)(n, x)'
    # assert typst(chebyshevu(n, x)) == r'chebyshevu(n, x)'
    # assert typst(chebyshevu(n, x)**2) == r'chebyshevu^(2)(n, x)'
    # assert typst(legendre(n, x)) == r'legendre(n, x)'
    # assert typst(legendre(n, x)**2) == r'legendre^(2)(n, x)'
    # assert typst(assoc_legendre(n, a, x)) == r'assoc_legendre(n, a, x)'
    # assert typst(assoc_legendre(n, a, x)**2) == r'assoc_legendre^(2)(n, a, x)'
    # assert typst(laguerre(n, x)) == r'laguerre(n, x)'
    # assert typst(laguerre(n, x)**2) == r'laguerre^(2)(n, x)'
    # assert typst(assoc_laguerre(n, a, x)) == r'assoc_laguerre(n, a, x)'
    # assert typst(assoc_laguerre(n, a, x)**2) == r'assoc_laguerre^(2)(n, a, x)'
    # assert typst(hermite(n, x)) == r'hermite(n, x)'
    # assert typst(hermite(n, x)**2) == r'hermite^(2)(n, x)'

    # theta = Symbol("theta", real=True)
    # phi = Symbol("phi", real=True)
    # assert typst(Ynm(n, m, theta, phi)) == r'Ynm(n, m, theta, phi)'
    # assert typst(Ynm(n, m, theta, phi)**3) == r'Ynm^(3)(n, m, theta, phi)'
    # assert typst(Znm(n, m, theta, phi)) == r'Znm(n, m, theta, phi)'
    # assert typst(Znm(n, m, theta, phi)**3) == r'Znm^(3)(n, m, theta, phi)'

    # # Test typst printing of function names with "_"
    # assert typst(polar_lift(0)) == r'upright("polar_lift")(0)'
    # assert typst(polar_lift(0)**3) == r'upright("polar_lift")^(3)(0)'

    # assert typst(totient(n)) == r'totient(n)'
    # assert typst(totient(n) ** 2) == r'totient^(2)(n)'

    # assert typst(reduced_totient(n)) == r'reduced_totient(n)'
    # assert typst(reduced_totient(n) ** 2) == r'reduced_totient^(2)(n)'

    # assert typst(divisor_sigma(x)) == r'divisor_sigma(x)'
    # assert typst(divisor_sigma(x)**2) == r'divisor_sigma^(2)(x)'
    # assert typst(divisor_sigma(x, y)) == r'divisor_sigma(x, y)'
    # assert typst(divisor_sigma(x, y)**2) == r'divisor_sigma^(2)(x, y)'

    # assert typst(udivisor_sigma(x)) == r'udivisor_sigma(x)'
    # assert typst(udivisor_sigma(x)**2) == r'udivisor_sigma^(2)(x)'
    # assert typst(udivisor_sigma(x, y)) == r'udivisor_sigma(x, y)'
    # assert typst(udivisor_sigma(x, y)**2) == r'udivisor_sigma^(2)(x, y)'

    # assert typst(primenu(n)) == r'primenu(n)'
    # assert typst(primenu(n) ** 2) == r'primenu^(2)(n)'

    # assert typst(primeomega(n)) == r'primeomega(n)'
    # assert typst(primeomega(n) ** 2) == r'primeomega^(2)(n)'

    # assert typst(LambertW(n)) == r'LambertW(n)'
    # assert typst(LambertW(n, -1)) == r'LambertW(-1, n)'
    # assert typst(LambertW(n, k)) == r'LambertW(k, n)'
    # assert typst(LambertW(n) * LambertW(n)) == r'LambertW(n) LambertW(n)'
    # assert typst(Pow(LambertW(n), 2)) == r'LambertW^(2)(n)'
    # assert typst(LambertW(n)**k) == r'LambertW^(k)(n)'
    # assert typst(LambertW(n, k)**p) == r'LambertW^(p)(k, n)'

    # assert typst(Mod(x, 7)) == r'Mod(x, 7)'
    # assert typst(Mod(x + 1, 7)) == r'Mod(x + 1, 7)'
    # assert typst(Mod(7, x + 1)) == r'Mod(7, x + 1)'
    # assert typst(Mod(2 * x, 7)) == r'Mod(2 x, 7)'
    # assert typst(Mod(7, 2 * x)) == r'Mod(7, 2 x)'
    # assert typst(Mod(x, 7) + 1) == r'Mod(x, 7) + 1'
    # assert typst(2 * Mod(x, 7)) == r'2 Mod(x, 7)'
    # assert typst(Mod(7, 2 * x)**n) == r'Mod(7, 2 x)^(n)'

    # # some unknown function name should get rendered with upright
    # fjlkd = Function('fjlkd')
    # assert typst(fjlkd(x)) == r'upright("fjlkd")(x)'
    # # even when it is referred to without an argument
    # assert typst(fjlkd) == r'upright("fjlkd")'

# test that notation passes to subclasses of the same name only
def test_function_subclass_different_name():
    class mygamma(gamma):
        pass
    assert typst(mygamma) == r'upright("mygamma")'
    assert typst(mygamma(x)) == r'upright("mygamma")(x)'


def test_boolean_args_order():
    syms = symbols('a:f')

    expr = And(*syms)
    assert typst(expr) == r'a and b and c and d and e and f'

    expr = Or(*syms)
    assert typst(expr) == r'a or b or c or d or e or f'

    expr = Equivalent(*syms)
    assert typst(expr) == \
        r'a <=> b <=> c <=> d <=> e <=> f'

    expr = Xor(*syms)
    assert typst(expr) == \
        r'a \u{22BB} b \u{22BB} c \u{22BB} d \u{22BB} e \u{22BB} f'

def test_issue_7180():
    assert typst(Equivalent(x, y)) == r"x <=> y"
    assert typst(Not(Equivalent(x, y))) == r"x arrow.l.r.double.not y"


def test_translate():
    s = 'Alpha'
    assert translate(s) == r'Alpha'
    s = 'Beta'
    assert translate(s) == r'Beta'
    s = 'Eta'
    assert translate(s) == r'Eta'
    s = 'omicron'
    assert translate(s) == r'omicron'
    s = 'Pi'
    assert translate(s) == r'Pi'
    s = 'pi'
    assert translate(s) == r'pi'
    s = 'LamdaHatDOT'
    assert translate(s) == r'accent(accent(Lambda, hat), dot)'


def test_other_symbols():
    from sympy.printing.typst import other_symbols
    for s in other_symbols:
        assert typst(symbols(s)) == "" + s


def test_modifiers():
    # Test each modifier individually in the simplest case
    # (with funny capitalizations)
    assert typst(symbols("xMathring")) == r"accent(x, circle)"
    assert typst(symbols("xCheck")) == r"accent(x, caron)"
    assert typst(symbols("xBreve")) == r"accent(x, breve)"
    assert typst(symbols("xAcute")) == r"accent(x, acute)"
    assert typst(symbols("xGrave")) == r"accent(x, grave)"
    assert typst(symbols("xTilde")) == r"accent(x, tilde)"
    assert typst(symbols("xPrime")) == r"(x)'"
    assert typst(symbols("xddDDot")) == r"accent(x, dot.quad)"
    assert typst(symbols("xDdDot")) == r"accent(x, dot.triple)"
    assert typst(symbols("xDDot")) == r"accent(x, dot.double)"
    assert typst(symbols("xBold")) == r"bold(x)"
    assert typst(symbols("xnOrM")) == r"norm(x)"
    assert typst(symbols("xAVG")) == r"lr(angle.l x angle.r)"
    assert typst(symbols("xHat")) == r"accent(x, hat)"
    assert typst(symbols("xDot")) == r"accent(x, dot)"
    assert typst(symbols("xBar")) == r"accent(x, macron)"
    assert typst(symbols("xVec")) == r"accent(x, arrow)"
    assert typst(symbols("xAbs")) == r"abs(x)"
    assert typst(symbols("xMag")) == r"abs(x)"
    assert typst(symbols("xPrM")) == r"(x)'"
    assert typst(symbols("xBM")) == r"bold(x)"
    # Test strings that are *only* the names of modifiers
    assert typst(symbols("Mathring")) == r'"Mathring"'
    assert typst(symbols("Check")) == r'"Check"'
    assert typst(symbols("Breve")) == r'"Breve"'
    assert typst(symbols("Acute")) == r'"Acute"'
    assert typst(symbols("Grave")) == r'"Grave"'
    assert typst(symbols("Tilde")) == r'"Tilde"'
    assert typst(symbols("Prime")) == r'"Prime"'
    assert typst(symbols("DDot")) == r"accent(D, dot)"
    assert typst(symbols("Bold")) == r'"Bold"'
    assert typst(symbols("NORm")) == r'"NORm"'
    assert typst(symbols("AVG")) == r'"AVG"'
    assert typst(symbols("Hat")) == r'"Hat"'
    assert typst(symbols("Dot")) == r'"Dot"'
    assert typst(symbols("Bar")) == r'"Bar"'
    assert typst(symbols("Vec")) == r'"Vec"'
    assert typst(symbols("Abs")) == r'"Abs"'
    assert typst(symbols("Mag")) == r'"Mag"'
    assert typst(symbols("PrM")) == r'"PrM"'
    assert typst(symbols("BM")) == r'"BM"'
    assert typst(symbols("hbar")) == r"plank.reduce"
    # Check a few combinations
    assert typst(symbols("xvecdot")) == r"accent(accent(x, arrow), dot)"
    assert typst(symbols("xDotVec")) == r"accent(accent(x, dot), arrow)"
    assert typst(symbols("xHATNorm")) == r"norm(accent(x, hat))"
    # Check a couple big, ugly combinations
    assert typst(symbols('xMathringBm_yCheckPRM__zbreveAbs')) == \
        r"bold(accent(x, circle))^(abs(accent(z, breve)))_((accent(y, caron))')"
    assert typst(symbols('alphadothat_nVECDOT__tTildePrime')) == \
        r"accent(accent(alpha, dot), hat)^((accent(t, tilde))')_(accent(accent(n, arrow), dot))"


def test_greek_symbols():
    assert typst(Symbol('alpha'))   == 'alpha'
    assert typst(Symbol('beta'))    == 'beta'
    assert typst(Symbol('gamma'))   == 'gamma'
    assert typst(Symbol('delta'))   == 'delta'
    assert typst(Symbol('epsilon')) == 'epsilon.alt'
    assert typst(Symbol('zeta'))    == 'zeta'
    assert typst(Symbol('eta'))     == 'eta'
    assert typst(Symbol('theta'))   == 'theta'
    assert typst(Symbol('iota'))    == 'iota'
    assert typst(Symbol('kappa'))   == 'kappa'
    assert typst(Symbol('lambda'))  == 'lambda'
    assert typst(Symbol('mu'))      == 'mu'
    assert typst(Symbol('nu'))      == 'nu'
    assert typst(Symbol('xi'))      == 'xi'
    assert typst(Symbol('omicron')) == 'omicron'
    assert typst(Symbol('pi'))      == 'pi'
    assert typst(Symbol('rho'))     == 'rho'
    assert typst(Symbol('sigma'))   == 'sigma'
    assert typst(Symbol('tau'))     == 'tau'
    assert typst(Symbol('upsilon')) == 'upsilon'
    assert typst(Symbol('phi'))     == 'phi.alt'
    assert typst(Symbol('chi'))     == 'chi'
    assert typst(Symbol('psi'))     == 'psi'
    assert typst(Symbol('omega'))   == 'omega'

    assert typst(Symbol('Alpha'))   == 'Alpha'
    assert typst(Symbol('Beta'))    == 'Beta'
    assert typst(Symbol('Gamma'))   == 'Gamma'
    assert typst(Symbol('Delta'))   == 'Delta'
    assert typst(Symbol('Epsilon')) == 'Epsilon'
    assert typst(Symbol('Zeta'))    == 'Zeta'
    assert typst(Symbol('Eta'))     == 'Eta'
    assert typst(Symbol('Theta'))   == 'Theta'
    assert typst(Symbol('Iota'))    == 'Iota'
    assert typst(Symbol('Kappa'))   == 'Kappa'
    assert typst(Symbol('Lambda'))  == 'Lambda'
    assert typst(Symbol('Mu'))      == 'Mu'
    assert typst(Symbol('Nu'))      == 'Nu'
    assert typst(Symbol('Xi'))      == 'Xi'
    assert typst(Symbol('Omicron')) == 'Omicron'
    assert typst(Symbol('Pi'))      == 'Pi'
    assert typst(Symbol('Rho'))     == 'Rho'
    assert typst(Symbol('Sigma'))   == 'Sigma'
    assert typst(Symbol('Tau'))     == 'Tau'
    assert typst(Symbol('Upsilon')) == 'Upsilon'
    assert typst(Symbol('Phi'))     == 'Phi'
    assert typst(Symbol('Chi'))     == 'Chi'
    assert typst(Symbol('Psi'))     == 'Psi'
    assert typst(Symbol('Omega'))   == 'Omega'

    assert typst(Symbol('varepsilon')) == 'epsilon'
    assert typst(Symbol('varkappa')) == 'kappa'
    assert typst(Symbol('varphi')) == 'phi'
    assert typst(Symbol('varpi')) == 'pi.alt'
    assert typst(Symbol('varrho')) == 'rho.alt'
    assert typst(Symbol('varsigma')) == 'sigma.alt'
    assert typst(Symbol('vartheta')) == 'theta.alt'

def test_fancyset_symbols():
    assert typst(S.Rationals) == r'QQ'
    assert typst(S.Naturals) == r'NN'
    assert typst(S.Naturals0) == r'NN_(0)'
    assert typst(S.Integers) == r'ZZ'
    assert typst(S.Reals) == r'RR'
    assert typst(S.Complexes) == r'CC'

def test_builtin_no_args():
        assert typst(Chi) == r'Chi'
        assert typst(beta) == r'Beta'
        assert typst(gamma) == r'Gamma'
        assert typst(KroneckerDelta) == r'delta'
        assert typst(DiracDelta) == r'delta'
        assert typst(lowergamma) == r'gamma'


def test_issue_6853():
    p = Function('Pi')
    assert typst(p(x)) == r'Pi(x)'


def test_Mul():
    e = Mul(-2, x + 1, evaluate=False)
    assert typst(e) == r'- 2 (x + 1)'
    e = Mul(2, x + 1, evaluate=False)
    assert typst(e) == r'2 (x + 1)'
    e = Mul(S.Half, x + 1, evaluate=False)
    assert typst(e) == r'(x + 1)/2'
    e = Mul(y, x + 1, evaluate=False)
    assert typst(e) == r'y (x + 1)'
    e = Mul(-y, x + 1, evaluate=False)
    assert typst(e) == r'- y (x + 1)'
    e = Mul(-2, x + 1)
    assert typst(e) == r'- 2 x - 2'
    e = Mul(2, x + 1)
    assert typst(e) == r'2 x + 2'


def test_Pow():
    e = Pow(2, 2, evaluate=False)
    assert typst(e) == r'2^(2)'
    assert typst(x**(Rational(-1, 3))) == r'1/root(3, x)'
    x2 = Symbol(r'x^2')
    assert typst(x2**2) == r'(x^(2))^(2)'
    # Issue 11011
    assert typst(S('1.453e4500')**x) == r'(1.453 dot 10^(4500))^(x)'


def test_issue_7180():
    assert typst(Equivalent(x, y)) == r"x <=> y"
    assert typst(Not(Equivalent(x, y))) == r"x arrow.l.r.double.not y"


def test_issue_8409():
    assert typst(S.Half**n) == r"(1/2)^(n)"


def test_issue_8470():
    from sympy.parsing.sympy_parser import parse_expr
    e = parse_expr("-B*A", evaluate=False)
    assert typst(e) == r"A (- B)"


def test_issue_15439():
    x = MatrixSymbol('x', 2, 2)
    y = MatrixSymbol('y', 2, 2)
    assert typst((x * y).subs(y, -y)) == r"x (- y)"
    assert typst((x * y).subs(y, -2*y)) == r"x (- 2 y)"
    assert typst((x * y).subs(x, -x)) == r"(- x) y"


# def test_issue_2934():
#     assert typst(Symbol(r'\frac{a_1}{b_1}')) == r'"\\frac{a_1}{b_1}"'


# def test_issue_10489():
#     latexSymbolWithBrace = r'C_{x_{0}}'
#     s = Symbol(latexSymbolWithBrace)
#     assert typst(s) == r'"C_{x_{0}}"'
#     assert typst(cos(s)) == r'cos("C_{x_{0}}")'


def test_issue_12886():
    m__1, l__1 = symbols('m__1, l__1')
    assert typst(m__1**2 + l__1**2) == \
        r'(l^(1))^(2) + (m^(1))^(2)'


def test_issue_13559():
    from sympy.parsing.sympy_parser import parse_expr
    expr = parse_expr('5/1', evaluate=False)
    assert typst(expr) == r"5/1"


def test_issue_13651():
    expr = c + Mul(-1, a + b, evaluate=False)
    assert typst(expr) == r"c - (a + b)"


def test_typst_UnevaluatedExpr():
    x = symbols("x")
    he = UnevaluatedExpr(1/x)
    assert typst(he) == typst(1/x) == r"1/x"
    assert typst(he**2) == r"(1/x)^(2)"
    assert typst(he + 1) == r"1 + 1/x"
    assert typst(x*he) == r"x 1/x"

def test_MatrixElement_printing():
    # test cases for issue #11821
    A = MatrixSymbol("A", 1, 3)
    B = MatrixSymbol("B", 1, 3)
    C = MatrixSymbol("C", 1, 3)

    assert typst(A[0, 0]) == r"A_(0,0)"
    assert typst(3 * A[0, 0]) == r"3 A_(0,0)"

    F = C[0, 0].subs(C, A - B)
    assert typst(F) == r"(A - B)_(0,0)"

    i, j, k = symbols("i j k")
    M = MatrixSymbol("M", k, k)
    N = MatrixSymbol("N", k, k)
    assert typst((M*N)[i, j]) == \
        r"sum_(i_(1)=0)^(k - 1) M_(i,i_(1)) N_(i_(1),j)"

    X_a = MatrixSymbol('X_a', 3, 3)
    assert typst(X_a[0, 0]) == r"X_(a)_(0,0)"


def test_MatrixSymbol_printing():
    # test cases for issue #14237
    A = MatrixSymbol("A", 3, 3)
    B = MatrixSymbol("B", 3, 3)
    C = MatrixSymbol("C", 3, 3)

    assert typst(-A) == r"- A"
    assert typst(A - A*B - B) == r"A - A B - B"
    assert typst(-A*B - A*B*C - B) == r"- A B - A B C - B"


def test_DotProduct_printing():
    X = MatrixSymbol('X', 3, 1)
    Y = MatrixSymbol('Y', 3, 1)
    a = Symbol('a')
    assert typst(DotProduct(X, Y)) == r"X dot Y"
    assert typst(DotProduct(a * X, Y)) == r"a X dot Y"
    assert typst(a * DotProduct(X, Y)) == r"a (X dot Y)"


def test_KroneckerProduct_printing():
    A = MatrixSymbol('A', 3, 3)
    B = MatrixSymbol('B', 2, 2)
    assert typst(KroneckerProduct(A, B)) == r"A times.circle B"


def test_Series_printing():
    tf1 = TransferFunction(x*y**2 - z, y**3 - t**3, y)
    tf2 = TransferFunction(x - y, x + y, y)
    tf3 = TransferFunction(t*x**2 - t**w*x + w, t - y, y)
    assert typst(Series(tf1, tf2)) == \
        r"((x y^(2) - z)/(- t^(3) + y^(3))) ((x - y)/(x + y))"
    assert typst(Series(tf1, tf2, tf3)) == \
        r"((x y^(2) - z)/(- t^(3) + y^(3))) ((x - y)/(x + y)) ((t x^(2) - t^(w) x + w)/(t - y))"
    assert typst(Series(-tf2, tf1)) == \
        r"((- x + y)/(x + y)) ((x y^(2) - z)/(- t^(3) + y^(3)))"

    M_1 = Matrix([[5/s], [5/(2*s)]])
    T_1 = TransferFunctionMatrix.from_Matrix(M_1, s)
    M_2 = Matrix([[5, 6*s**3]])
    T_2 = TransferFunctionMatrix.from_Matrix(M_2, s)
    # Brackets
    assert typst(T_1*(T_2 + T_2)) == \
        r'mat(delim: "[", 5/s; 5/(2 s))_tau dot (mat(delim: "[", 5/1, (6 s^(3))/1)_tau +' \
        r' mat(delim: "[", 5/1, (6 s^(3))/1)_tau)' \
            == typst(MIMOSeries(MIMOParallel(T_2, T_2), T_1))
    # No Brackets
    M_3 = Matrix([[5, 6], [6, 5/s]])
    T_3 = TransferFunctionMatrix.from_Matrix(M_3, s)
    assert typst(T_1*T_2 + T_3) == r'mat(delim: "[", 5/s; 5/(2 s))_tau dot ' \
        r'mat(delim: "[", 5/1, (6 s^(3))/1)_tau + mat(delim: "[", 5/1, 6/1; 6/1, 5/s)_tau' \
            == typst(MIMOParallel(MIMOSeries(T_2, T_1), T_3))


def test_TransferFunction_printing():
    tf1 = TransferFunction(x - 1, x + 1, x)
    assert typst(tf1) == r"(x - 1)/(x + 1)"
    tf2 = TransferFunction(x + 1, 2 - y, x)
    assert typst(tf2) == r"(x + 1)/(2 - y)"
    tf3 = TransferFunction(y, y**2 + 2*y + 3, y)
    assert typst(tf3) == r"y/(y^(2) + 2 y + 3)"

def test_Parallel_printing():
    tf1 = TransferFunction(x*y**2 - z, y**3 - t**3, y)
    tf2 = TransferFunction(x - y, x + y, y)
    assert typst(Parallel(tf1, tf2)) == \
        r'(x y^(2) - z)/(- t^(3) + y^(3)) + (x - y)/(x + y)'
    assert typst(Parallel(-tf2, tf1)) == \
        r'(- x + y)/(x + y) + (x y^(2) - z)/(- t^(3) + y^(3))'

    M_1 = Matrix([[5, 6], [6, 5/s]])
    T_1 = TransferFunctionMatrix.from_Matrix(M_1, s)
    M_2 = Matrix([[5/s, 6], [6, 5/(s - 1)]])
    T_2 = TransferFunctionMatrix.from_Matrix(M_2, s)
    M_3 = Matrix([[6, 5/(s*(s - 1))], [5, 6]])
    T_3 = TransferFunctionMatrix.from_Matrix(M_3, s)
    assert typst(T_1 + T_2 + T_3) == \
        r'mat(delim: "[", 5/1, 6/1; 6/1, 5/s)_tau + mat(delim: "[", 5/s, 6/1; 6/1, 5/(s - 1))_tau + mat(delim: "[", 6/1, 5/(s (s - 1)); 5/1, 6/1)_tau' \
            == typst(MIMOParallel(T_1, T_2, T_3)) == typst(MIMOParallel(T_1, MIMOParallel(T_2, T_3))) == typst(MIMOParallel(MIMOParallel(T_1, T_2), T_3))


def test_TransferFunctionMatrix_printing():
    tf1 = TransferFunction(p, p + x, p)
    tf2 = TransferFunction(-s + p, p + s, p)
    tf3 = TransferFunction(p, y**2 + 2*y + 3, p)
    assert typst(TransferFunctionMatrix([[tf1], [tf2]])) == \
        r'mat(delim: "[", p/(p + x); (p - s)/(p + s))_tau'
    assert typst(TransferFunctionMatrix([[tf1, tf2], [tf3, -tf1]])) == \
        r'mat(delim: "[", p/(p + x), (p - s)/(p + s); p/(y^(2) + 2 y + 3), ((-1) p)/(p + x))_tau'


def test_Feedback_printing():
    tf1 = TransferFunction(p, p + x, p)
    tf2 = TransferFunction(-s + p, p + s, p)
    # Negative Feedback (Default)
    assert typst(Feedback(tf1, tf2)) == \
        r'(p/(p + x))/(1/1 + (p/(p + x)) ((p - s)/(p + s)))'
    assert typst(Feedback(tf1*tf2, TransferFunction(1, 1, p))) == \
        r'((p/(p + x)) ((p - s)/(p + s)))/(1/1 + (p/(p + x)) ((p - s)/(p + s)))'
    # Positive Feedback
    assert typst(Feedback(tf1, tf2, 1)) == \
        r'(p/(p + x))/(1/1 - (p/(p + x)) ((p - s)/(p + s)))'
    assert typst(Feedback(tf1*tf2, sign=1)) == \
        r'((p/(p + x)) ((p - s)/(p + s)))/(1/1 - (p/(p + x)) ((p - s)/(p + s)))'


def test_MIMOFeedback_printing():
    tf1 = TransferFunction(1, s, s)
    tf2 = TransferFunction(s, s**2 - 1, s)
    tf3 = TransferFunction(s, s - 1, s)
    tf4 = TransferFunction(s**2, s**2 - 1, s)

    tfm_1 = TransferFunctionMatrix([[tf1, tf2], [tf3, tf4]])
    tfm_2 = TransferFunctionMatrix([[tf4, tf3], [tf2, tf1]])

    # Negative Feedback (Default)
    assert typst(MIMOFeedback(tfm_1, tfm_2)) == \
        r'(I_(tau) + mat(delim: "[", 1/s, s/(s^(2) - 1); s/(s - 1), s^(2)/(s^(2) - 1))_tau dot mat(delim: "[", s^(2)/(s^(2) - 1), s/(s - 1); s/(s^(2) - 1), 1/s)_tau)^(-1) dot mat(delim: "[", 1/s, s/(s^(2) - 1); s/(s - 1), s^(2)/(s^(2) - 1))_tau'

    # Positive Feedback
    assert typst(MIMOFeedback(tfm_1*tfm_2, tfm_1, 1)) == \
    r'(I_(tau) - mat(delim: "[", 1/s, s/(s^(2) - 1); s/(s - 1), s^(2)/(s^(2) - 1))_tau dot mat(delim: "[", s^(2)/(s^(2) - 1), s/(s - 1); s/(s^(2) - 1), 1/s)_tau dot mat(delim: "[", 1/s, s/(s^(2) - 1); s/(s - 1), s^(2)/(s^(2) - 1))_tau)^(-1) dot mat(delim: "[", 1/s, s/(s^(2) - 1); s/(s - 1), s^(2)/(s^(2) - 1))_tau dot mat(delim: "[", s^(2)/(s^(2) - 1), s/(s - 1); s/(s^(2) - 1), 1/s)_tau'

def test_Quaternion_typst_printing():
    q = Quaternion(x, y, z, t)
    assert typst(q) == r"x + y i + z j + t k"
    q = Quaternion(x, y, z, x*t)
    assert typst(q) == r"x + y i + z j + t x k"
    q = Quaternion(x, y, z, x + t)
    assert typst(q) == r"x + y i + z j + (t + x) k"


def test_TensorProduct_printing():
    from sympy.tensor.functions import TensorProduct
    A = MatrixSymbol("A", 3, 3)
    B = MatrixSymbol("B", 3, 3)
    assert typst(TensorProduct(A, B)) == r"A times.circle B"


def test_WedgeProduct_printing():
    from sympy.diffgeom.rn import R2
    from sympy.diffgeom import WedgeProduct
    wp = WedgeProduct(R2.dx, R2.dy)
    assert typst(wp) == r'upright(d)x and upright(d)y'


def test_issue_9216():
    expr_1 = Pow(1, -1, evaluate=False)
    assert typst(expr_1) == r"1^(-1)"

    expr_2 = Pow(1, Pow(1, -1, evaluate=False), evaluate=False)
    assert typst(expr_2) == r"1^(1^(-1))"

    expr_3 = Pow(3, -2, evaluate=False)
    assert typst(expr_3) == r"1/9"

    expr_4 = Pow(1, -2, evaluate=False)
    assert typst(expr_4) == r"1^(-2)"

def test_typst_printer_tensor():
    from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorHead, tensor_heads
    L = TensorIndexType("L")
    i, j, k, l = tensor_indices("i j k l", L)
    i0 = tensor_indices("i_0", L)
    A, B, C, D = tensor_heads("A B C D", [L])
    H = TensorHead("H", [L, L])
    K = TensorHead("K", [L, L, L, L])

    assert typst(i) == r'scripts("")^(i)'
    assert typst(-i) == r'scripts("")_(i)'

    expr = A(i)
    assert typst(expr) == r'A scripts("")^(i)'

    expr = A(i0)
    assert typst(expr) == r'A scripts("")^(i_(0))'

    expr = A(-i)
    assert typst(expr) == r'A scripts("")_(i)'

    expr = -3*A(i)
    assert typst(expr) == r'-3A scripts("")^(i)'

    expr = K(i, j, -k, -i0)
    assert typst(expr) == r'K scripts("")^(i) scripts("")^(j) scripts("")_(k) scripts("")_(i_(0))'

    expr = K(i, -j, -k, i0)
    assert typst(expr) == r'K scripts("")^(i) scripts("")_(j) scripts("")_(k) scripts("")^(i_(0))'

    expr = K(i, -j, k, -i0)
    assert typst(expr) == r'K scripts("")^(i) scripts("")_(j) scripts("")^(k) scripts("")_(i_(0))'

    expr = H(i, -j)
    assert typst(expr) == r'H scripts("")^(i) scripts("")_(j)'

    expr = H(i, j)
    assert typst(expr) == r'H scripts("")^(i) scripts("")^(j)'

    expr = H(-i, -j)
    assert typst(expr) == r'H scripts("")_(i) scripts("")_(j)'

    expr = (1+x)*A(i)
    assert typst(expr) == r'(x + 1)A scripts("")^(i)'

    expr = H(i, -i)
    assert typst(expr) == r'H scripts("")^(L_(0)) scripts("")_(L_(0))'

    expr = H(i, -j)*A(j)*B(k)
    assert typst(expr) == r'H scripts("")^(i) scripts("")_(L_(0))A scripts("")^(L_(0))B scripts("")^(k)'

    expr = A(i) + 3*B(i)
    assert typst(expr) == r'3B scripts("")^(i) + A scripts("")^(i)'

    # Test ``TensorElement``:
    from sympy.tensor.tensor import TensorElement

    expr = TensorElement(K(i, j, k, l), {i: 3, k: 2})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")^(,j) scripts("")^(,k=2) scripts("")^(,l)'

    expr = TensorElement(K(i, j, k, l), {i: 3})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")^(,j) scripts("")^(k) scripts("")^(l)'

    expr = TensorElement(K(i, -j, k, l), {i: 3, k: 2})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")_(j) scripts("")^(k=2) scripts("")^(,l)'

    expr = TensorElement(K(i, -j, k, -l), {i: 3, k: 2})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")_(j) scripts("")^(k=2) scripts("")_(l)'

    expr = TensorElement(K(i, j, -k, -l), {i: 3, -k: 2})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")^(,j) scripts("")_(k=2) scripts("")_(,l)'

    expr = TensorElement(K(i, j, -k, -l), {i: 3})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")^(,j) scripts("")_(k) scripts("")_(l)'

    expr = PartialDerivative(A(i), A(i))
    assert typst(expr) == r'partial / (partial A scripts("")^(L_(0))) A scripts("")^(L_(0))'

    expr = PartialDerivative(A(-i), A(-j))
    assert typst(expr) == r'partial / (partial A scripts("")_(j)) A scripts("")_(i)'

    expr = PartialDerivative(K(i, j, -k, -l), A(m), A(-n))
    assert typst(expr) == r'partial^(2) / (partial A scripts("")^(m) partial A scripts("")_(n)) '\
                          r'K scripts("")^(i) scripts("")^(j) scripts("")_(k) scripts("")_(l)'

    expr = PartialDerivative(B(-i) + A(-i), A(-j), A(-n))
    assert typst(expr) == r'partial^(2) / (partial A scripts("")_(j) partial A scripts("")_(n)) '\
                          r'(A scripts("")_(i) + B scripts("")_(i))'

    expr = PartialDerivative(3*A(-i), A(-j), A(-n))
    assert typst(expr) == r'partial^(2) / (partial A scripts("")_(j) partial A scripts("")_(n)) '\
                          r'(3A scripts("")_(i))'


def test_multiline_typst():
    # TODO
    pass

def test_issue_15353():
    a, x = symbols('a x')
    # Obtained from nonlinsolve([(sin(a*x)),cos(a*x)],[x,a])
    sol = ConditionSet(
        Tuple(x, a), Eq(sin(a*x), 0) & Eq(cos(a*x), 0), S.Complexes**2)
    assert typst(sol) == \
        r'{(x, #h(0.5em)a) | (x, #h(0.5em)a) in CC^(2) and sin(a x) = 0 and cos(a x) = 0}'

def test_typst_symbolic_probability():
    mu = symbols("mu")
    sigma = symbols("sigma", positive=True)
    X = Normal("X", mu, sigma)
    assert typst(Expectation(X)) == r'upright("E")[X]'
    assert typst(Variance(X)) == r'upright("Var")(X)'
    assert typst(Probability(X > 0)) == r'upright("P")(X > 0)'
    Y = Normal("Y", mu, sigma)
    assert typst(Covariance(X, Y)) == r'upright("Cov")(X, Y)'


def test_trace():
    # Issue 15303
    from sympy.matrices.expressions.trace import trace
    A = MatrixSymbol("A", 2, 2)
    assert typst(trace(A)) == r'upright("tr")(A )'
    assert typst(trace(A**2)) == r'upright("tr")(A^(2) )'


def test_print_basic():
    # Issue 15303
    from sympy.core.basic import Basic
    from sympy.core.expr import Expr

    # dummy class for testing printing where the function is not
    # implemented in typst.py
    class UnimplementedExpr(Expr):
        def __new__(cls, e):
            return Basic.__new__(cls, e)

    # dummy function for testing
    def unimplemented_expr(expr):
        return UnimplementedExpr(expr).doit()

    # override class name to use superscript / subscript
    def unimplemented_expr_sup_sub(expr):
        result = UnimplementedExpr(expr)
        result.__class__.__name__ = 'UnimplementedExpr_x^1'
        return result

    assert typst(unimplemented_expr(x)) == r'upright("UnimplementedExpr")(x)'
    assert typst(unimplemented_expr(x**2)) == \
        r'upright("UnimplementedExpr")(x^(2))'
    assert typst(unimplemented_expr_sup_sub(x)) == \
        r'upright("UnimplementedExpr"^(1)_(x))(x)'


def test_MatrixSymbol_bold():
    # Issue #15871
    from sympy.matrices.expressions.trace import trace
    A = MatrixSymbol("A", 2, 2)
    assert typst(trace(A), mat_symbol_style='bold') == \
        r'upright("tr")(bold(A) )'
    assert typst(trace(A), mat_symbol_style='plain') == \
        r'upright("tr")(A )'

    A = MatrixSymbol("A", 3, 3)
    B = MatrixSymbol("B", 3, 3)
    C = MatrixSymbol("C", 3, 3)

    assert typst(-A, mat_symbol_style='bold') == r'- bold(A)'
    assert typst(A - A*B - B, mat_symbol_style='bold') == \
        r'bold(A) - bold(A) bold(B) - bold(B)'
    assert typst(-A*B - A*B*C - B, mat_symbol_style='bold') == \
        r'- bold(A) bold(B) - bold(A) bold(B) bold(C) - bold(B)'

    A_k = MatrixSymbol("A_k", 3, 3)
    assert typst(A_k, mat_symbol_style='bold') == r'bold(A)_(k)'

    A = MatrixSymbol(r"\nabla_k", 3, 3)
    assert typst(A, mat_symbol_style='bold') == r'bold(\nabla)_(k)'

def test_AppliedPermutation():
    p = Permutation(0, 1, 2)
    x = Symbol('x')
    assert typst(AppliedPermutation(p, x)) == \
        r'sigma_((0 space 1 space 2))(x)'


def test_PermutationMatrix():
    p = Permutation(0, 1, 2)
    assert typst(PermutationMatrix(p)) == r'P_((0 space 1 space 2))'
    p = Permutation(0, 3)(1, 2)
    assert typst(PermutationMatrix(p)) == \
        r'P_((0 space 3)(1 space 2))'

def test_issue_21758():
    from sympy.functions.elementary.piecewise import piecewise_fold
    from sympy.series.fourier import FourierSeries
    x = Symbol('x')
    k, n = symbols('k n')
    fo = FourierSeries(x, (x, -pi, pi), (0, SeqFormula(0, (k, 1, oo)), SeqFormula(
        Piecewise((-2*pi*cos(n*pi)/n + 2*sin(n*pi)/n**2, (n > -oo) & (n < oo) & Ne(n, 0)),
                  (0, True))*sin(n*x)/pi, (n, 1, oo))))
    assert typst(piecewise_fold(fo)) == r'cases( 2 sin(x) - sin(2 x) + (2 sin(3 x))/3 + ... '\
        r'& upright("for") n > -infinity and n < infinity and n != 0, 0 & upright("otherwise") )'
    assert typst(FourierSeries(x, (x, -pi, pi), (0, SeqFormula(0, (k, 1, oo)),
                                                 SeqFormula(0, (n, 1, oo))))) == '0'
    
def test_typst_diffgeom():
    from sympy.diffgeom import Manifold, Patch, CoordSystem, BaseScalarField, Differential
    from sympy.diffgeom.rn import R2
    x, y = symbols('x y', real=True)
    m = Manifold('M', 2)
    assert typst(m) == r'upright("M")'
    p = Patch('P', m)
    assert typst(p) == r'upright("P")_(upright("M"))'
    rect = CoordSystem('rect', p, [x, y])
    assert typst(rect) == r'upright("rect")^(upright("P"))_(upright("M"))'
    b = BaseScalarField(rect, 0)
    assert typst(b) == r'bold(x)'

    g = Function('g')
    s_field = g(R2.x, R2.y)
    assert typst(Differential(s_field)) == \
        r'upright(d)(g(bold(x), bold(y)))'


def test_typst_unit_printing():
    assert typst(5*meter) == r'5 upright("m")'
    assert typst(3*gibibyte) == r'3 upright("gibibyte")'
    assert typst(4*microgram/second) == r'(4 mu upright("g"))/(upright("s"))'
    assert typst(4*micro*gram/second) == r'(4 mu upright("g"))/(upright("s"))'
    assert typst(5*milli*meter) == r'5 upright("m") upright("m")'
    assert typst(milli) == r'upright("m")'


def test_issue_17092_typst():
    x_star = Symbol('x^*')
    assert typst(Derivative(x_star, x_star,2)) == r'd^(2) / (d (x^(*))^(2)) x^(*)'


def test_typst_decimal_separator():

    x, y, z, t = symbols('x y z t')
    k, m, n = symbols('k m n', integer=True)
    f, g, h = symbols('f g h', cls=Function)

    # comma decimal_separator
    assert(typst([1, 2.3, 4.5], decimal_separator='comma') == r'[1; #h(0.5em)2,3; #h(0.5em)4,5]')
    assert(typst(FiniteSet(1, 2.3, 4.5), decimal_separator='comma') == r'{1; 2,3; 4,5}')
    assert(typst((1, 2.3, 4.6), decimal_separator = 'comma') == r'(1; #h(0.5em)2,3; #h(0.5em)4,6)')
    assert(typst((1,), decimal_separator='comma') == r'(1;)')

    # period decimal_separator
    assert(typst([1, 2.3, 4.5], decimal_separator='period') == r'[1, #h(0.5em)2.3, #h(0.5em)4.5]' )
    assert(typst(FiniteSet(1, 2.3, 4.5), decimal_separator='period') == r'{1, 2.3, 4.5}')
    assert(typst((1, 2.3, 4.6), decimal_separator = 'period') == r'(1, #h(0.5em)2.3, #h(0.5em)4.6)')
    assert(typst((1,), decimal_separator='period') == r'(1,)')

    # default decimal_separator
    assert(typst([1, 2.3, 4.5]) == r'[1, #h(0.5em)2.3, #h(0.5em)4.5]')
    assert(typst(FiniteSet(1, 2.3, 4.5)) == r'{1, 2.3, 4.5}')
    assert(typst((1, 2.3, 4.6)) == r'(1, #h(0.5em)2.3, #h(0.5em)4.6)')
    assert(typst((1,)) == r'(1,)')

    assert(typst(Mul(3.4,5.3), decimal_separator = 'comma') == r'18,02')
    assert(typst(3.4*5.3, decimal_separator = 'comma') == r'18,02')
    x = symbols('x')
    y = symbols('y')
    z = symbols('z')
    assert(typst(x*5.3 + 2**y**3.4 + 4.5 + z, decimal_separator = 'comma') == r'2^(y^(3,4)) + 5,3 x + z + 4,5')

    assert(typst(0.987, decimal_separator='comma') == r'0,987')
    assert(typst(S(0.987), decimal_separator='comma') == r'0,987')
    assert(typst(.3, decimal_separator='comma') == r'0,3')
    assert(typst(S(.3), decimal_separator='comma') == r'0,3')


    assert(typst(5.8*10**(-7), decimal_separator='comma') == r'5,8 dot 10^(-7)')
    assert(typst(S(5.7)*10**(-7), decimal_separator='comma') == r'5,7 dot 10^(-7)')
    assert(typst(S(5.7*10**(-7)), decimal_separator='comma') == r'5,7 dot 10^(-7)')

    x = symbols('x')
    assert(typst(1.2*x+3.4, decimal_separator='comma') == r'1,2 x + 3,4')
    assert(typst(FiniteSet(1, 2.3, 4.5), decimal_separator='period') == r'{1, 2.3, 4.5}')

    # Error Handling tests
    raises(ValueError, lambda: typst([1,2.3,4.5], decimal_separator='non_existing_decimal_separator_in_list'))
    raises(ValueError, lambda: typst(FiniteSet(1,2.3,4.5), decimal_separator='non_existing_decimal_separator_in_set'))
    raises(ValueError, lambda: typst((1,2.3,4.5), decimal_separator='non_existing_decimal_separator_in_tuple'))

def test_Str():
    from sympy.core.symbol import Str
    assert str(Str('x')) == r'x'

def test_emptyPrinter():
    class MyObject:
        def __repr__(self):
            return "<MyObject with {...}>"

    # unknown objects are monospaced
    assert typst(MyObject()) == r'"<MyObject with {...}>"'

    # even if they are nested within other objects
    assert typst((MyObject(),)) == r'("<MyObject with {...}>",)'

def test_global_settings():
    import inspect

    # settings should be visible in the signature of `typst`
    assert inspect.signature(typst).parameters['imaginary_unit'].default == r'i'
    assert typst(I) == r'i'
    try:
        # but changing the defaults...
        from sympy.printing.typst import TypstPrinter
        TypstPrinter.set_global_settings(imaginary_unit='j')
        # ... should change the signature
        assert inspect.signature(typst).parameters['imaginary_unit'].default == r'j'
        assert typst(I) == r'j'
    finally:
        # there's no public API to undo this, but we need to make sure we do
        # so as not to impact other tests
        del TypstPrinter._global_settings['imaginary_unit']

    # check we really did undo it
    assert inspect.signature(typst).parameters['imaginary_unit'].default == r'i'
    assert typst(I) == r'i'

def test_pickleable():
    # this tests that the _PrintFunction instance is pickleable
    import pickle
    assert pickle.loads(pickle.dumps(typst)) is typst

def test_printing_typst_array_expressions():
    assert typst(ArraySymbol("A", (2, 3, 4))) == "A"
    assert typst(ArrayElement("A", (2, 1/(1-x), 0))) == r"A_(2, 1/(1 - x), 0)"
    M = MatrixSymbol("M", 3, 3)
    N = MatrixSymbol("N", 3, 3)
    assert typst(ArrayElement(M*N, [x, 0])) == r"(M N)_(x, 0)"

def test_Array():
    arr = Array(range(10))
    assert typst(arr) == r'mat(delim: "[", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 )'

    arr = Array(range(11))
    assert typst(arr) == r'mat(delim: "[", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 )'

def test_typst_with_unevaluated():
    with evaluate(False):
        assert typst(a * a) == r"a a"

def test_typst_indexed():
    Psi_symbol = Symbol('Psi_0', complex=True, real=False)
    Psi_indexed = IndexedBase(Symbol('Psi', complex=True, real=False))
    symbol_typst = typst(Psi_symbol * conjugate(Psi_symbol))
    indexed_typst = typst(Psi_indexed[0] * conjugate(Psi_indexed[0]))
    # overline(Psi_(0)) Psi_(0)  vs.  Psi_(0) overline(Psi_(0))
    assert symbol_typst == r'Psi_(0) overline(Psi_(0))'
    assert indexed_typst == r'overline(Psi_(0)) Psi_(0)'

    interval = r'.. '
    assert typst(Indexed('x1', Symbol('i'))) == r'x_1_(i)'
    assert typst(Indexed('x2', Idx('i'))) == r'x_2_(i)'
    assert typst(Indexed('x3', Idx('i', Symbol('N')))) == r'x_3_(i_(0.. N - 1))'
    assert typst(Indexed('x3', Idx('i', Symbol('N')+1))) == r'x_3_(i_(0.. N))'
    assert typst(Indexed('x4', Idx('i', (Symbol('a'),Symbol('b'))))) == r'x_4_(i_(a.. b))'
    assert typst(IndexedBase('gamma')) == r'gamma'
    assert typst(IndexedBase('a b')) == r'"a b"'
    assert typst(IndexedBase('a_b')) == r'a_(b)'

def test_typst_derivatives():
    # regular "d" for ordinary derivatives
    assert typst(diff(x**3, x, evaluate=False)) == \
        r"d / (d x) x^(3)"
    assert typst(diff(sin(x) + x**2, x, evaluate=False)) == \
        r"d / (d x) (x^(2) + sin(x))"
    assert typst(diff(diff(sin(x) + x**2, x, evaluate=False), evaluate=False))\
        == \
        r"d^(2) / (d x^(2)) (x^(2) + sin(x))"
    assert typst(diff(diff(diff(sin(x) + x**2, x, evaluate=False), evaluate=False), evaluate=False)) == \
        r"d^(3) / (d x^(3)) (x^(2) + sin(x))"

    # partial for partial derivatives
    assert typst(diff(sin(x * y), x, evaluate=False)) == \
        r"partial / (partial x) sin(x y)"
    assert typst(diff(sin(x * y) + x**2, x, evaluate=False)) == \
        r"partial / (partial x) (x^(2) + sin(x y))"
    assert typst(diff(diff(sin(x*y) + x**2, x, evaluate=False), x, evaluate=False)) == \
        r"partial^(2) / (partial x^(2)) (x^(2) + sin(x y))"
    assert typst(diff(diff(diff(sin(x*y) + x**2, x, evaluate=False), x, evaluate=False), x, evaluate=False)) == \
        r"partial^(3) / (partial x^(3)) (x^(2) + sin(x y))"

    # mixed partial derivatives
    f = Function("f")
    assert typst(diff(diff(f(x, y), x, evaluate=False), y, evaluate=False)) == \
        r"partial^(2) / (partial y partial x) " + typst(f(x, y))

    assert typst(diff(diff(diff(f(x, y), x, evaluate=False), x, evaluate=False), y, evaluate=False)) == \
        r"partial^(3) / (partial y partial x^(2)) " + typst(f(x, y))

    # for negative nested Derivative
    assert typst(diff(-diff(y**2,x,evaluate=False),x,evaluate=False)) == r'd / (d x) (- d / (d x) y^(2))'
    assert typst(diff(diff(-diff(diff(y,x,evaluate=False),x,evaluate=False),x,evaluate=False),x,evaluate=False)) == \
        r'd^(2) / (d x^(2)) (- d^(2) / (d x^(2)) y)'

    # # use ordinary d when one of the variables has been integrated out
    # assert typst(diff(Integral(exp(-x*y), (x, 0, oo)), y, evaluate=False)) == \
    #     r"(d / d y) int(limits_(0)^(infty) e^(- x y), dx"

    # Derivative wrapped in power:
    assert typst(diff(x, x, evaluate=False)**2) == \
        r"(d / (d x) x)^(2)"

    assert typst(diff(f(x), x)**2) == \
        r"(d / (d x) f(x))^(2)"

    assert typst(diff(f(x), (x, n))) == \
        r"d^(n) / (d x^(n)) f(x)"

    x1 = Symbol('x1')
    x2 = Symbol('x2')
    assert typst(diff(f(x1, x2), x1)) == r'partial / (partial x_(1)) f(x_(1), x_(2))'

    n1 = Symbol('n1')
    assert typst(diff(f(x), (x, n1))) == r'd^(n_(1)) / (d x^(n_(1))) f(x)'

    n2 = Symbol('n2')
    assert typst(diff(f(x), (x, Max(n1, n2)))) == \
        r'd^(max(n_(1), n_(2))) / (d x^(max(n_(1), n_(2)))) f(x)'

    # set diff operator
    assert typst(diff(f(x), x), diff_operator="rd") == r'bold(d) / (bold(d) x) f(x)'

def test_issue_17092():
    x_star = Symbol('x^*')
    assert typst(Derivative(x_star, x_star,2)) == r'd^(2) / (d (x^(*))^(2)) x^(*)'

def test_typst_Piecewise():
    p = Piecewise((x, x < 1), (x**2, True))
    assert typst(p) == r'cases( x & upright("for") x < 1, x^(2) &' \
                       r' upright("otherwise") )'
    p = Piecewise((x, x < 0), (0, x >= 0))
    assert typst(p) == r'cases( x & upright("for") x < 0, 0 &' \
                       r' upright("otherwise") )'
    A, B = symbols("A B", commutative=False)
    p = Piecewise((A**2, Eq(A, B)), (A*B, True))
    s = r'cases( A^(2) & upright("for") A = B, A B & upright("otherwise") )'
    assert typst(p) == s
    assert typst(A*p) == r'A (%s)' % s
    assert typst(p*A) == r'(%s) A' % s
    assert typst(Piecewise((x, x < 1), (x**2, x < 2))) == \
        r'cases( x & ' \
        r'upright("for") x < 1, x^(2) & upright("for") x < 2 )'

def test_typst_Matrix():
    M = Matrix([[1 + x, y], [y, x - 1]])
    assert typst(M) == \
        r'mat(delim: "[", x + 1, y; y, x - 1)'
    assert typst(M, mat_delim=None) == \
        r'mat(delim: #none, x + 1, y; y, x - 1)'

    M2 = Matrix(1, 11, range(11))
    assert typst(M2) == \
        r'mat(delim: "[", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)'


def test_typst_matrix_with_functions():
    t = symbols('t')
    theta1 = symbols('theta1', cls=Function)

    M = Matrix([[sin(theta1(t)), cos(theta1(t))],
                [cos(theta1(t).diff(t)), sin(theta1(t).diff(t))]])

    expected = (r'mat(delim: "[", sin(theta_(1)(t)), '\
                r'cos(theta_(1)(t)); '\
                r'cos(d / (d t) theta_(1)(t)), '\
                r'sin(d / (d t) theta_(1)(t)))')

    assert typst(M) == expected

def test_typst_NDimArray():
    x, y, z, w = symbols("x y z w")

    for ArrayType in (ImmutableDenseNDimArray, ImmutableSparseNDimArray,
                      MutableDenseNDimArray, MutableSparseNDimArray):
        # Basic: scalar array
        M = ArrayType(x)

        assert typst(M) == r"x"

        M = ArrayType([[1 / x, y], [z, w]])
        M1 = ArrayType([1 / x, y, z])

        M2 = tensorproduct(M1, M)
        M3 = tensorproduct(M, M)

        assert typst(M) == \
            r'mat(delim: "[", 1/x, y; z, w )'
        assert typst(M1) == \
            r'mat(delim: "[", 1/x, y, z )'
        assert typst(M2) == \
            r'mat(delim: "[", mat(delim: "[", 1/x^(2), y/x; z/x, w/x ), '\
            r'mat(delim: "[", y/x, y^(2); y z, w y ), '\
            r'mat(delim: "[", z/x, y z; z^(2), w z ) )'
        assert typst(M3) == \
            r'mat(delim: "[", mat(delim: "[", 1/x^(2), y/x; z/x, w/x ), '\
            r'mat(delim: "[", y/x, y^(2); y z, w y ); '\
            r'mat(delim: "[", z/x, y z; z^(2), w z ), '\
            r'mat(delim: "[", w/x, w y; w z, w^(2) ) )'

        Mrow = ArrayType([[x, y, 1/z]])
        Mcolumn = ArrayType([[x], [y], [1/z]])
        Mcol2 = ArrayType([Mcolumn.tolist()])

        assert typst(Mrow) == \
            r'[mat(delim: "[", x, y, 1/z )]'
        assert typst(Mcolumn) == \
            r'mat(delim: "[", x; y; 1/z )'
        assert typst(Mcol2) == \
            r'mat(delim: "[", mat(delim: "[", x; y; 1/z ) )'

def test_matAdd():
    C = MatrixSymbol('C', 5, 5)
    B = MatrixSymbol('B', 5, 5)

    n = symbols("n")
    h = MatrixSymbol("h", 1, 1)

    assert typst(C - 2*B) in [r'- 2 B + C', r'C -2 B']
    assert typst(C + 2*B) in [r'2 B + C', r'C + 2 B']
    assert typst(B - 2*C) in [r'B - 2 C', r'- 2 C + B']
    assert typst(B + 2*C) in [r'B + 2 C', r'2 C + B']

    assert typst(n * h - (-h + h.T) * (h + h.T)) == r'n h - (- h + h^(T)) (h + h^(T))'
    assert typst(MatAdd(MatAdd(h, h), MatAdd(h, h))) == r'(h + h) + (h + h)'
    assert typst(MatMul(MatMul(h, h), MatMul(h, h))) == r'(h h) (h h)'


def test_matMul():
    A = MatrixSymbol('A', 5, 5)
    B = MatrixSymbol('B', 5, 5)
    x = Symbol('x')
    assert typst(2*A) == r'2 A'
    assert typst(2*x*A) == r'2 x A'
    assert typst(-2*A) == r'- 2 A'
    assert typst(1.5*A) == r'1.5 A'
    assert typst(sqrt(2)*A) == r'sqrt(2) A'
    assert typst(-sqrt(2)*A) == r'- sqrt(2) A'
    assert typst(2*sqrt(2)*x*A) == r'2 sqrt(2) x A'
    assert typst(-2*A*(A + 2*B)) in [r'- 2 A (A + 2 B)',
                                        r'- 2 A (2 B + A)']


def test_typst_MatrixSlice():
    n = Symbol('n', integer=True)
    x, y, z, w, t, = symbols('x y z w t')
    X = MatrixSymbol('X', n, n)
    Y = MatrixSymbol('Y', 10, 10)
    Z = MatrixSymbol('Z', 10, 10)

    assert typst(MatrixSlice(X, (None, None, None), (None, None, None))) == r'X[:, :]'
    assert typst(X[x:x + 1, y:y + 1]) == r'X[x:x + 1, y:y + 1]'
    assert typst(X[x:x + 1:2, y:y + 1:2]) == r'X[x:x + 1:2, y:y + 1:2]'
    assert typst(X[:x, y:]) == r'X[:x, y:]'
    assert typst(X[:x, y:]) == r'X[:x, y:]'
    assert typst(X[x:, :y]) == r'X[x:, :y]'
    assert typst(X[x:y, z:w]) == r'X[x:y, z:w]'
    assert typst(X[x:y:t, w:t:x]) == r'X[x:y:t, w:t:x]'
    assert typst(X[x::y, t::w]) == r'X[x::y, t::w]'
    assert typst(X[:x:y, :t:w]) == r'X[:x:y, :t:w]'
    assert typst(X[::x, ::y]) == r'X[::x, ::y]'
    assert typst(MatrixSlice(X, (0, None, None), (0, None, None))) == r'X[:, :]'
    assert typst(MatrixSlice(X, (None, n, None), (None, n, None))) == r'X[:, :]'
    assert typst(MatrixSlice(X, (0, n, None), (0, n, None))) == r'X[:, :]'
    assert typst(MatrixSlice(X, (0, n, 2), (0, n, 2))) == r'X[::2, ::2]'
    assert typst(X[1:2:3, 4:5:6]) == r'X[1:2:3, 4:5:6]'
    assert typst(X[1:3:5, 4:6:8]) == r'X[1:3:5, 4:6:8]'
    assert typst(X[1:10:2]) == r'X[1:10:2, :]'
    assert typst(Y[:5, 1:9:2]) == r'Y[:5, 1:9:2]'
    assert typst(Y[:5, 1:10:2]) == r'Y[:5, 1::2]'
    assert typst(Y[5, :5:2]) == r'Y[5:6, :5:2]'
    assert typst(X[0:1, 0:1]) == r'X[:1, :1]'
    assert typst(X[0:1:2, 0:1:2]) == r'X[:1:2, :1:2]'
    assert typst((Y + Z)[2:, 2:]) == r'(Y + Z)[2:, 2:]'

def test_typst_printer_tensor():
    from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorHead, tensor_heads
    L = TensorIndexType("L")
    i, j, k, l = tensor_indices("i j k l", L)
    i0 = tensor_indices("i_0", L)
    A, B, C, D = tensor_heads("A B C D", [L])
    H = TensorHead("H", [L, L])
    K = TensorHead("K", [L, L, L, L])

    assert typst(i) == r'scripts("")^(i)'
    assert typst(-i) == r'scripts("")_(i)'

    expr = A(i)
    assert typst(expr) == r'A scripts("")^(i)'

    expr = A(i0)
    assert typst(expr) == r'A scripts("")^(i_(0))'

    expr = A(-i)
    assert typst(expr) == r'A scripts("")_(i)'

    expr = -3*A(i)
    assert typst(expr) == r'-3A scripts("")^(i)'

    expr = K(i, j, -k, -i0)
    assert typst(expr) == r'K scripts("")^(i) scripts("")^(j) scripts("")_(k) scripts("")_(i_(0))'

    expr = K(i, -j, -k, i0)
    assert typst(expr) == r'K scripts("")^(i) scripts("")_(j) scripts("")_(k) scripts("")^(i_(0))'

    expr = K(i, -j, k, -i0)
    assert typst(expr) == r'K scripts("")^(i) scripts("")_(j) scripts("")^(k) scripts("")_(i_(0))'

    expr = H(i, -j)
    assert typst(expr) == r'H scripts("")^(i) scripts("")_(j)'

    expr = H(i, j)
    assert typst(expr) == r'H scripts("")^(i) scripts("")^(j)'

    expr = H(-i, -j)
    assert typst(expr) == r'H scripts("")_(i) scripts("")_(j)'

    expr = (1+x)*A(i)
    assert typst(expr) == r'(x + 1)A scripts("")^(i)'

    expr = H(i, -i)
    assert typst(expr) == r'H scripts("")^(L_(0)) scripts("")_(L_(0))'

    expr = H(i, -j)*A(j)*B(k)
    assert typst(expr) == r'H scripts("")^(i) scripts("")_(L_(0))A scripts("")^(L_(0))B scripts("")^(k)'

    expr = A(i) + 3*B(i)
    assert typst(expr) == r'3B scripts("")^(i) + A scripts("")^(i)'

    # Test ``TensorElement``:
    from sympy.tensor.tensor import TensorElement

    expr = TensorElement(K(i, j, k, l), {i: 3, k: 2})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")^(,j) scripts("")^(,k=2) scripts("")^(,l)'

    expr = TensorElement(K(i, j, k, l), {i: 3})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")^(,j) scripts("")^(k) scripts("")^(l)'

    expr = TensorElement(K(i, -j, k, l), {i: 3, k: 2})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")_(j) scripts("")^(k=2) scripts("")^(,l)'

    expr = TensorElement(K(i, -j, k, -l), {i: 3, k: 2})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")_(j) scripts("")^(k=2) scripts("")_(l)'

    expr = TensorElement(K(i, j, -k, -l), {i: 3, -k: 2})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")^(,j) scripts("")_(k=2) scripts("")_(,l)'

    expr = TensorElement(K(i, j, -k, -l), {i: 3})
    assert typst(expr) == r'K scripts("")^(i=3) scripts("")^(,j) scripts("")_(k) scripts("")_(l)'

    expr = PartialDerivative(A(i), A(i))
    assert typst(expr) == r'partial / (partial A scripts("")^(L_(0))) A scripts("")^(L_(0))'

    expr = PartialDerivative(A(-i), A(-j))
    assert typst(expr) == r'partial / (partial A scripts("")_(j)) A scripts("")_(i)'

    expr = PartialDerivative(K(i, j, -k, -l), A(m), A(-n))
    assert typst(expr) == r'partial^(2) / (partial A scripts("")^(m) partial A scripts("")_(n)) '\
                          r'K scripts("")^(i) scripts("")^(j) scripts("")_(k) scripts("")_(l)'

    expr = PartialDerivative(B(-i) + A(-i), A(-j), A(-n))
    assert typst(expr) == r'partial^(2) / (partial A scripts("")_(j) partial A scripts("")_(n)) '\
                          r'(A scripts("")_(i) + B scripts("")_(i))'

    expr = PartialDerivative(3*A(-i), A(-j), A(-n))
    assert typst(expr) == r'partial^(2) / (partial A scripts("")_(j) partial A scripts("")_(n)) '\
                          r'(3A scripts("")_(i))'

def test_trace():
    # Issue 15303
    from sympy.matrices.expressions.trace import trace
    A = MatrixSymbol("A", 2, 2)
    assert typst(trace(A)) == r'upright("tr")(A )'
    assert typst(trace(A**2)) == r'upright("tr")(A^(2) )'

def test_printing_typst_array_expressions():
    assert typst(ArraySymbol("A", (2, 3, 4))) == "A"
    assert typst(ArrayElement("A", (2, 1/(1-x), 0))) == r"A_(2, 1/(1 - x), 0)"
    M = MatrixSymbol("M", 3, 3)
    N = MatrixSymbol("N", 3, 3)
    assert typst(ArrayElement(M*N, [x, 0])) == r"(M N)_(x, 0)"


def test_typst_subs():
    assert typst(Subs(x*y, (x, y), (1, 2))) == r'x y|_(x=1 \ y=2)'