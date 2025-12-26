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


from sympy.testing.pytest import (raises, _both_exp_pow,
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
    assert typst(frac(x)**2) == r'upright("frac")(x)^(2)'

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
    assert typst(gamma(x)) == r"Gamma(x)"
    w = Wild('w')
    assert typst(gamma(w)) == r"Gamma(w)"
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

    assert typst(zeta(x)) == r'zeta(x)'
    assert typst(zeta(x)**2) == r'zeta^(2)(x)'
    assert typst(zeta(x, y)) == r'zeta(x, y)'
    assert typst(zeta(x, y)**2) == r'zeta^(2)(x, y)'
    assert typst(dirichlet_eta(x)) == r'eta(x)'
    assert typst(dirichlet_eta(x)**2) == r'eta^(2)(x)'
    assert typst(polylog(x, y)) == r'upright("Li")_(x)(y)'
    assert typst(polylog(x, y)**2) == r'upright("Li")_(x)^(2)(y)'
    assert typst(lerchphi(x, y, n)) == r'Phi(x, y, n)'
    assert typst(lerchphi(x, y, n)**2) == r'Phi^(2)(x, y, n)'
    assert typst(stieltjes(x)) == r'gamma_(x)'
    assert typst(stieltjes(x)**2) == r'gamma_(x)^(2)'
    assert typst(stieltjes(x, y)) == r'gamma_(x)(y)'
    assert typst(stieltjes(x, y)**2) == r'gamma_(x)(y)^(2)'

    assert typst(elliptic_k(z)) == r'K(z)'
    assert typst(elliptic_k(z)**2) == r'K^(2)(z)'
    assert typst(elliptic_f(x, y)) == r'F(x mid(|) y)'
    assert typst(elliptic_f(x, y)**2) == r'F^(2)(x mid(|) y)'
    assert typst(elliptic_e(x, y)) == r'E(x mid(|) y)'
    assert typst(elliptic_e(x, y)**2) == r'E^(2)(x mid(|) y)'
    assert typst(elliptic_e(z)) == r'E(z)'
    assert typst(elliptic_e(z)**2) == r'E^(2)(z)'
    assert typst(elliptic_pi(x, y, z)) == r'Pi(x; y mid(|) z)'
    assert typst(elliptic_pi(x, y, z)**2) == r'Pi^(2)(x; y mid(|) z)'
    assert typst(elliptic_pi(x, y)) == r'Pi(x mid(|) y)'
    assert typst(elliptic_pi(x, y)**2) == r'Pi^(2)(x mid(|) y)'

    assert typst(Ei(x)) == r'upright("Ei")(x)'
    assert typst(Ei(x)**2) == r'upright("Ei")^(2)(x)'
    assert typst(expint(x, y)) == r'upright(E)_(x)(y)'
    assert typst(expint(x, y)**2) == r'upright(E)_(x)^(2)(y)'
    assert typst(Shi(x)**2) == r'upright("Shi")^(2)(x)'
    assert typst(Si(x)**2) == r'upright("Si")^(2)(x)'
    assert typst(Ci(x)**2) == r'upright("Ci")^(2)(x)'
    assert typst(Chi(x)**2) == r'upright("Chi")^(2)(x)'
    assert typst(Chi(x)) == r'upright("Chi")(x)'
    assert typst(jacobi(n, a, b, x)) == r'P_(n)^(a,b)(x)'
    assert typst(jacobi(n, a, b, x)**2) == r'(P_(n)^(a,b)(x))^(2)'
    assert typst(gegenbauer(n, a, x)) == r'C_(n)^((a))(x)'
    assert typst(gegenbauer(n, a, x)**2) == r'(C_(n)^((a))(x))^(2)'
    assert typst(chebyshevt(n, x)) == r'T_(n)(x)'
    assert typst(chebyshevt(n, x)**2) == r'(T_(n)(x))^(2)'
    assert typst(chebyshevu(n, x)) == r'U_(n)(x)'
    assert typst(chebyshevu(n, x)**2) == r'(U_(n)(x))^(2)'
    assert typst(legendre(n, x)) == r'P_(n)(x)'
    assert typst(legendre(n, x)**2) == r'(P_(n)(x))^(2)'
    assert typst(assoc_legendre(n, a, x)) == r'P_(n)^(a)(x)'
    assert typst(assoc_legendre(n, a, x)**2) == r'(P_(n)^(a)(x))^(2)'
    assert typst(laguerre(n, x)) == r'L_(n)(x)'
    assert typst(laguerre(n, x)**2) == r'(L_(n)(x))^(2)'
    assert typst(assoc_laguerre(n, a, x)) == r'L_(n)^(a)(x)'
    assert typst(assoc_laguerre(n, a, x)**2) == r'(L_(n)^(a)(x))^(2)'
    assert typst(hermite(n, x)) == r'H_(n)(x)'
    assert typst(hermite(n, x)**2) == r'(H_(n)(x))^(2)'

    theta = Symbol("theta", real=True)
    phi = Symbol("phi", real=True)
    assert typst(Ynm(n, m, theta, phi)) == r'Y_(n)^(m)(theta, phi.alt)'
    assert typst(Ynm(n, m, theta, phi)**3) == r'(Y_(n)^(m)(theta, phi.alt))^(3)'
    assert typst(Znm(n, m, theta, phi)) == r'Z_(n)^(m)(theta, phi.alt)'
    assert typst(Znm(n, m, theta, phi)**3) == r'(Z_(n)^(m)(theta, phi.alt))^(3)'

    # Test typst printing of function names with "_"
    assert typst(polar_lift(0)) == r'upright("polar_lift")(0)'
    assert typst(polar_lift(0)**3) == r'upright("polar_lift")^(3)(0)'

    assert typst(totient(n)) == r'phi(n)'
    assert typst(totient(n) ** 2) == r'(phi(n))^(2)'

    assert typst(reduced_totient(n)) == r'lambda(n)'
    assert typst(reduced_totient(n) ** 2) == r'(lambda(n))^(2)'

    assert typst(divisor_sigma(x)) == r'sigma(x)'
    assert typst(divisor_sigma(x)**2) == r'sigma^(2)(x)'
    assert typst(divisor_sigma(x, y)) == r'sigma_(y)(x)'
    assert typst(divisor_sigma(x, y)**2) == r'sigma^(2)_(y)(x)'

    assert typst(udivisor_sigma(x)) == r'sigma^*(x)'
    assert typst(udivisor_sigma(x)**2) == r'sigma^*^(2)(x)'
    assert typst(udivisor_sigma(x, y)) == r'sigma^*_(y)(x)'
    assert typst(udivisor_sigma(x, y)**2) == r'sigma^*^(2)_(y)(x)'

    assert typst(primenu(n)) == r'nu(n)'
    assert typst(primenu(n) ** 2) == r'(nu(n))^(2)'

    assert typst(primeomega(n)) == r'Omega(n)'
    assert typst(primeomega(n) ** 2) == r'(Omega(n))^(2)'

    assert typst(LambertW(n)) == r'W(n)'
    assert typst(LambertW(n, -1)) == r'W_(-1)(n)'
    assert typst(LambertW(n, k)) == r'W_(k)(n)'
    assert typst(LambertW(n) * LambertW(n)) == r'W^(2)(n)'
    assert typst(Pow(LambertW(n), 2)) == r'W^(2)(n)'
    assert typst(LambertW(n)**k) == r'W^(k)(n)'
    assert typst(LambertW(n, k)**p) == r'W^(p)_(k)(n)'

    assert typst(Mod(x, 7)) == r'x mod 7'
    assert typst(Mod(x + 1, 7)) == r'(x + 1) mod 7'
    assert typst(Mod(7, x + 1)) == r'7 mod (x + 1)'
    assert typst(Mod(2 * x, 7)) == r'2 x mod 7'
    assert typst(Mod(7, 2 * x)) == r'7 mod 2 x'
    assert typst(Mod(x, 7) + 1) == r'(x mod 7) + 1'
    assert typst(2 * Mod(x, 7)) == r'2 (x mod 7)'
    assert typst(Mod(7, 2 * x)**n) == r'(7 mod 2 x)^(n)'

    # some unknown function name should get rendered with upright
    fjlkd = Function('fjlkd')
    assert typst(fjlkd(x)) == r'upright("fjlkd")(x)'
    # even when it is referred to without an argument
    assert typst(fjlkd) == r'upright("fjlkd")'

# test that notation passes to subclasses of the same name only
def test_function_subclass_different_name():
    class mygamma(gamma):
        pass
    assert typst(mygamma) == r'upright("mygamma")'
    assert typst(mygamma(x)) == r'upright("mygamma")(x)'

def test_hyper_printing():
    from sympy.abc import x, z

    assert typst(meijerg(Tuple(pi, pi, x), Tuple(1),
                         (0, 1), Tuple(1, 2, 3/pi), z)) == \
        r'G_(4, 5)^(2, 3)(mat(delim:#none, pi"," pi"," x, 1;' \
        r'0"," 1 , 1"," 2"," 3/(pi)) mid(|) z)'
    assert typst(meijerg(Tuple(), Tuple(1), (0,), Tuple(), z)) == \
        r'G_(1, 1)^(1, 0)(mat(delim:#none, , 1;0 , ) mid(|) z)'
    assert typst(hyper((x, 2), (3,), z)) == \
        r'scripts("")_(2)F_(1)(mat(delim: "[", 2"," x; 3) mid(|) z)'
    assert typst(hyper(Tuple(), Tuple(1), z)) == \
        r'scripts("")_(0)F_(1)(mat(delim: "[", ; 1) mid(|) z)'


def test_typst_bessel():
    from sympy.functions.special.bessel import (besselj, bessely, besseli,
                                                besselk, hankel1, hankel2,
                                                jn, yn, hn1, hn2)
    from sympy.abc import z
    assert typst(besselj(n, z**2)**k) == r'J^(k)_(n)(z^(2))'
    assert typst(bessely(n, z)) == r'Y_(n)(z)'
    assert typst(besseli(n, z)) == r'I_(n)(z)'
    assert typst(besselk(n, z)) == r'K_(n)(z)'
    assert typst(hankel1(n, z**2)**2) == \
        r'(H^((1))_(n)(z^(2)))^(2)'
    assert typst(hankel2(n, z)) == r'H^((2))_(n)(z)'
    assert typst(jn(n, z)) == r'j_(n)(z)'
    assert typst(yn(n, z)) == r'y_(n)(z)'
    assert typst(hn1(n, z)) == r'h^((1))_(n)(z)'
    assert typst(hn2(n, z)) == r'h^((2))_(n)(z)'


def test_typst_fresnel():
    from sympy.functions.special.error_functions import (fresnels, fresnelc)
    from sympy.abc import z
    assert typst(fresnels(z)) == r'S(z)'
    assert typst(fresnelc(z)) == r'C(z)'
    assert typst(fresnels(z)**2) == r'S^(2)(z)'
    assert typst(fresnelc(z)**2) == r'C^(2)(z)'


def test_typst_brackets():
    assert typst((-1)**x) == r"(-1)^(x)"


def test_typst_indexed():
    Psi_symbol = Symbol('Psi_0', complex=True, real=False)
    Psi_indexed = IndexedBase(Symbol('Psi', complex=True, real=False))
    symbol_typst = typst(Psi_symbol * conjugate(Psi_symbol))
    indexed_typst = typst(Psi_indexed[0] * conjugate(Psi_indexed[0]))
    # overline(Psi_(0)) Psi_(0)  vs.  Psi_(0) overline(Psi_(0))
    assert symbol_typst == r'Psi_(0) overline(Psi_(0))'
    assert indexed_typst == r'overline(Psi_(0)) Psi_(0)'

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

def test_typst_subs():
    assert typst(Subs(x*y, (x, y), (1, 2))) == r'x y|_(x=1 \ y=2)'

def test_typst_integrals():
    assert typst(Integral(log(x), x)) == r"integral log(x) d x"
    assert typst(Integral(x**2, (x, 0, 1))) == r"integral_(0)^(1) x^(2) d x"
    assert typst(Integral(x**2, (x, 10, 20))) == r"integral_(10)^(20) x^(2) d x"
    assert typst(Integral(y*x**2, (x, 0, 1), y)) == r"integral integral_(0)^(1) x^(2) y d x d y"
    assert typst(Integral(y*x**2, (x, 0, 1), y), mode='equation*') == r"integral integral_(0)^(1) x^(2) y d x d y"
    assert typst(Integral(y*x**2, (x, 0, 1), y), mode='equation*') == r"integral integral_(0)^(1) x^(2) y d x d y"
    assert typst(Integral(x, (x, 0))) == r"integral^(0) x d x"
    assert typst(Integral(x*y, x, y)) == r"integral.double x y d x d y"
    assert typst(Integral(x*y*z, x, y, z)) == r"integral.triple x y z d x d y d z"
    assert typst(Integral(x*y*z*t, x, y, z, t)) == r"integral.quad t x y z d x d y d z d t"
    assert typst(Integral(x, x, x, x, x, x, x)) == r"integral integral integral integral integral integral x d x d x d x d x d x d x"
    assert typst(Integral(x, x, y, (z, 0, 1))) == r"integral_(0)^(1) integral integral x d x d y d z"

    # for negative nested Integral
    assert typst(Integral(-Integral(y**2,x),x)) == r"integral (- integral y^(2) d x) d x"
    assert typst(Integral(-Integral(-Integral(y,x),x),x)) == r"integral (- integral (- integral y d x) d x) d x"

    # fix issue #10806
    assert typst(Integral(z, z)**2) == r"(integral z d z)^(2)"
    assert typst(Integral(x + z, z)) == r"integral (x + z) d z"
    assert typst(Integral(x+z/2, z)) == r"integral (x + z/2) d z"
    assert typst(Integral(x**y, z)) == r"integral x^(y) d z"

    # set diff operator
    assert typst(Integral(x, x), diff_operator="rd") == r"integral x bold(d) x"
    assert typst(Integral(x, (x, 0, 1)), diff_operator="rd") == r"integral_(0)^(1) x bold(d) x"


def test_typst_sets():
    for s in (frozenset, set):
        assert typst(s([x*y, x**2])) == r"{x^(2), x y}"
        assert typst(s(range(1, 6))) == r"{1, 2, 3, 4, 5}"
        assert typst(s(range(1, 13))) == r"{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}"

    s = FiniteSet
    assert typst(s(*[x*y, x**2])) == r"{x^(2), x y}"
    assert typst(s(*range(1, 6))) == r"{1, 2, 3, 4, 5}"
    assert typst(s(*range(1, 13))) == r"{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}"


def test_typst_SetExpr():
    iv = Interval(1, 3)
    se = SetExpr(iv)
    assert typst(se) == r'upright("SetExpr")([1, 3])'


def test_typst_Range():
    assert typst(Range(1, 51)) == r"{1, 2, ..., 50}"
    assert typst(Range(1, 4)) == r"{1, 2, 3}"
    assert typst(Range(0, 3, 1)) == r"{0, 1, 2}"
    assert typst(Range(0, 30, 1)) == r"{0, 1, ..., 29}"
    assert typst(Range(30, 1, -1)) == r"{30, 29, ..., 2}"
    assert typst(Range(0, oo, 2)) == r"{0, 2, ...}"
    assert typst(Range(oo, -2, -2)) == r"{..., 2, 0}"
    assert typst(Range(-2, -oo, -1)) == r"{-2, -3, ...}"
    assert typst(Range(-oo, oo)) == r"{..., -1, 0, 1, ...}"
    assert typst(Range(oo, -oo, -1)) == r"{..., 1, 0, -1, ...}"

    a, b, c = symbols('a:c')
    assert typst(Range(a, b, c)) == r'upright("Range")(a, b, c)'
    assert typst(Range(a, 10, 1)) == r'upright("Range")(a, 10)'
    assert typst(Range(0, b, 1)) == r'upright("Range")(b)'
    assert typst(Range(0, 10, c)) == r'upright("Range")(0, 10, c)'

    i = Symbol('i', integer=True)
    n = Symbol('n', negative=True, integer=True)
    p = Symbol('p', positive=True, integer=True)

    assert typst(Range(i, i + 3)) == r"{i, i + 1, i + 2}"
    assert typst(Range(-oo, n, 2)) == r"{..., n - 4, n - 2}"
    assert typst(Range(p, oo)) == r"{p, p + 1, ...}"
    # The following will work if __iter__ is improved
    # assert typst(Range(-3, p + 7)) == r"{-3, -2, .., p + 6}"
    # Must have integer assumptions
    assert typst(Range(a, a + 3)) == r'upright("Range")(a, a + 3)'


def test_typst_sequences():
    s1 = SeqFormula(a**2, (0, oo))
    s2 = SeqPer((1, 2))

    typst_str = r'[0, 1, 4, 9, ...]'
    assert typst(s1) == typst_str

    typst_str = r'[1, 2, 1, 2, ...]'
    assert typst(s2) == typst_str

    s3 = SeqFormula(a**2, (0, 2))
    s4 = SeqPer((1, 2), (0, 2))

    typst_str = r'[0, 1, 4]'
    assert typst(s3) == typst_str

    typst_str = r'[1, 2, 1]'
    assert typst(s4) == typst_str

    s5 = SeqFormula(a**2, (-oo, 0))
    s6 = SeqPer((1, 2), (-oo, 0))

    typst_str = r'[..., 9, 4, 1, 0]'
    assert typst(s5) == typst_str

    typst_str = r'[..., 2, 1, 2, 1]'
    assert typst(s6) == typst_str

    typst_str = r'[1, 3, 5, 11, ...]'
    assert typst(SeqAdd(s1, s2)) == typst_str

    typst_str = r'[1, 3, 5]'
    assert typst(SeqAdd(s3, s4)) == typst_str

    typst_str = r'[..., 11, 5, 3, 1]'
    assert typst(SeqAdd(s5, s6)) == typst_str

    typst_str = r'[0, 2, 4, 18, ...]'
    assert typst(SeqMul(s1, s2)) == typst_str

    typst_str = r'[0, 2, 4]'
    assert typst(SeqMul(s3, s4)) == typst_str

    typst_str = r'[..., 18, 4, 2, 0]'
    assert typst(SeqMul(s5, s6)) == typst_str

    # Sequences with symbolic limits, issue 12629
    s7 = SeqFormula(a**2, (a, 0, x))
    typst_str = r'{a^(2)}_(a=0)^(x)'
    assert typst(s7) == typst_str

    b = Symbol('b')
    s8 = SeqFormula(b*a**2, (a, 0, 2))
    typst_str = r'[0, b, 4 b]'
    assert typst(s8) == typst_str


def test_typst_FourierSeries():
    typst_str = \
        r'2 sin(x) - sin(2 x) + (2 sin(3 x))/3 + ...'
    assert typst(fourier_series(x, (x, -pi, pi))) == typst_str


def test_typst_FormalPowerSeries():
    typst_str = r'sum_(k=1)^(infinity) - ((-1)^(- k) x^(k))/k'
    assert typst(fps(log(1 + x))) == typst_str


def test_typst_intervals():
    a = Symbol('a', real=True)
    assert typst(Interval(0, 0)) == r"{0}"
    assert typst(Interval(0, a)) == r"[0, a]"
    assert typst(Interval(0, a, False, False)) == r"[0, a]"
    assert typst(Interval(0, a, True, False)) == r"(0, a]"
    assert typst(Interval(0, a, False, True)) == r"[0, a)"
    assert typst(Interval(0, a, True, True)) == r"(0, a)"


def test_typst_AccumuBounds():
    a = Symbol('a', real=True)
    assert typst(AccumBounds(0, 1)) == r"angle.l 0, 1 angle.r"
    assert typst(AccumBounds(0, a)) == r"angle.l 0, a angle.r"
    assert typst(AccumBounds(a + 1, a + 2)) == \
        r"angle.l a + 1, a + 2 angle.r"


def test_typst_emptyset():
    assert typst(S.EmptySet) == r"nothing"


def test_typst_universalset():
    assert typst(S.UniversalSet) == r"UU"


def test_typst_commutator():
    A = Operator('A')
    B = Operator('B')
    comm = Commutator(B, A)
    assert typst(comm.doit()) == r"- (A B - B A)"


def test_typst_union():
    assert typst(Union(Interval(0, 1), Interval(2, 3))) == \
        r"[0, 1] union [2, 3]"
    assert typst(Union(Interval(1, 1), Interval(2, 2), Interval(3, 4))) == \
        r"{1, 2} union [3, 4]"


def test_typst_intersection():
    assert typst(Intersection(Interval(0, 1), Interval(x, y))) == \
        r"[0, 1] inter [x, y]"


def test_typst_symmetric_difference():
    assert typst(SymmetricDifference(Interval(2, 5), Interval(4, 7),
                                     evaluate=False)) == \
        r'[2, 5] triangle [4, 7]'


def test_typst_Complement():
    assert typst(Complement(S.Reals, S.Naturals)) == \
        r'RR without NN'


def test_typst_productset():
    line = Interval(0, 1)
    bigline = Interval(0, 10)
    fset = FiniteSet(1, 2, 3)
    assert typst(line**2) == r'[0, 1]^(2)'
    assert typst(line**10) == r'[0, 1]^(10)'
    assert typst((line * bigline * fset).flatten()) == r'[0, 1] times [0, 10] times {1, 2, 3}'


def test_typst_powerset():
    fset = FiniteSet(1, 2, 3)
    assert typst(PowerSet(fset)) == r'cal(P)({1, 2, 3})'


def test_typst_ordinals():
    w = OrdinalOmega()
    assert typst(w) == r'omega'
    wp = OmegaPower(2, 3)
    assert typst(wp) == r'3 omega^(2)'
    assert typst(Ordinal(wp, OmegaPower(1, 1))) == r'3 omega^(2) + omega'
    assert typst(Ordinal(OmegaPower(2, 1), OmegaPower(1, 2))) == r'omega^(2) + 2 omega'


def test_set_operators_parenthesis_typst():
    a, b, c, d = symbols('a:d')
    A = FiniteSet(a)
    B = FiniteSet(b)
    C = FiniteSet(c)
    D = FiniteSet(d)

    U1 = Union(A, B, evaluate=False)
    U2 = Union(C, D, evaluate=False)
    I1 = Intersection(A, B, evaluate=False)
    I2 = Intersection(C, D, evaluate=False)
    C1 = Complement(A, B, evaluate=False)
    C2 = Complement(C, D, evaluate=False)
    D1 = SymmetricDifference(A, B, evaluate=False)
    D2 = SymmetricDifference(C, D, evaluate=False)
    # XXX ProductSet does not support evaluate keyword
    P1 = ProductSet(A, B)
    P2 = ProductSet(C, D)

    assert typst(Intersection(A, U2, evaluate=False)) == \
        r'{a} inter ({c} union {d})'
    assert typst(Intersection(U1, U2, evaluate=False)) == \
        r'({a} union {b}) inter ({c} union {d})'
    assert typst(Intersection(C1, C2, evaluate=False)) == \
        r'({a} without {b}) inter ({c} without {d})'
    assert typst(Intersection(D1, D2, evaluate=False)) == \
        r'({a} triangle {b}) inter ({c} triangle {d})'
    assert typst(Intersection(P1, P2, evaluate=False)) == \
        r'({a} times {b}) inter ({c} times {d})'

    assert typst(Union(A, I2, evaluate=False)) == \
        r'{a} union ({c} inter {d})'
    assert typst(Union(I1, I2, evaluate=False)) == \
        r'({a} inter {b}) union ({c} inter {d})'
    assert typst(Union(C1, C2, evaluate=False)) == \
        r'({a} without {b}) union ({c} without {d})'
    assert typst(Union(D1, D2, evaluate=False)) == \
        r'({a} triangle {b}) union ({c} triangle {d})'
    assert typst(Union(P1, P2, evaluate=False)) == \
        r'({a} times {b}) union ({c} times {d})'

    assert typst(Complement(A, C2, evaluate=False)) == \
        r'{a} without ({c} without {d})'
    assert typst(Complement(U1, U2, evaluate=False)) == \
        r'({a} union {b}) without ({c} union {d})'
    assert typst(Complement(I1, I2, evaluate=False)) == \
        r'({a} inter {b}) without ({c} inter {d})'
    assert typst(Complement(D1, D2, evaluate=False)) == \
        r'({a} triangle {b}) without ({c} triangle {d})'
    assert typst(Complement(P1, P2, evaluate=False)) == \
        r'({a} times {b}) without ({c} times {d})'

    assert typst(SymmetricDifference(A, D2, evaluate=False)) == \
        r'{a} triangle ({c} triangle {d})'
    assert typst(SymmetricDifference(U1, U2, evaluate=False)) == \
        r'({a} union {b}) triangle ({c} union {d})'
    assert typst(SymmetricDifference(I1, I2, evaluate=False)) == \
        r'({a} inter {b}) triangle ({c} inter {d})'
    assert typst(SymmetricDifference(C1, C2, evaluate=False)) == \
        r'({a} without {b}) triangle ({c} without {d})'
    assert typst(SymmetricDifference(P1, P2, evaluate=False)) == \
        r'({a} times {b}) triangle ({c} times {d})'

    # XXX This can be incorrect since cartesian product is not associative
    assert typst(ProductSet(A, P2).flatten()) == \
        r'{a} times {c} times {d}'
    assert typst(ProductSet(U1, U2)) == \
        r'({a} union {b}) times ({c} union {d})'
    assert typst(ProductSet(I1, I2)) == \
        r'({a} inter {b}) times ({c} inter {d})'
    assert typst(ProductSet(C1, C2)) == \
        r'({a} without {b}) times ({c} without {d})'
    assert typst(ProductSet(D1, D2)) == \
        r'({a} triangle {b}) times ({c} triangle {d})'


def test_typst_Complexes():
    assert typst(S.Complexes) == r'CC'


def test_typst_Naturals():
    assert typst(S.Naturals) == r'NN'


def test_typst_Naturals0():
    assert typst(S.Naturals0) == r'NN_(0)'


def test_typst_Integers():
    assert typst(S.Integers) == r'ZZ'


def test_typst_ImageSet():
    x = Symbol('x')
    assert typst(ImageSet(Lambda(x, x**2), S.Naturals)) == \
        r'{x^(2) mid(|) x in NN}'

    y = Symbol('y')
    imgset = ImageSet(Lambda((x, y), x + y), {1, 2, 3}, {3, 4})
    assert typst(imgset) == \
        r'{x + y mid(|) x in {1, 2, 3}, y in {3, 4}}'

    imgset = ImageSet(Lambda(((x, y),), x + y), ProductSet({1, 2, 3}, {3, 4}))
    assert typst(imgset) == \
        r'{x + y mid(|) (x, #h(0.5em)y) in {1, 2, 3} times {3, 4}}'


def test_typst_ConditionSet():
    x = Symbol('x')
    assert typst(ConditionSet(x, Eq(x**2, 1), S.Reals)) == \
        r'{x mid(|) x in RR and x^(2) = 1}'
    assert typst(ConditionSet(x, Eq(x**2, 1), S.UniversalSet)) == \
        r'{x mid(|) x^(2) = 1}'


def test_typst_ComplexRegion():
    assert typst(ComplexRegion(Interval(3, 5)*Interval(4, 6))) == \
        r'{x + y i mid(|) x, y in [3, 5] times [4, 6]}'
    assert typst(ComplexRegion(Interval(0, 1)*Interval(0, 2*pi), polar=True)) == \
        r'{r (i sin(theta) + cos(theta)) mid(|) r, theta in [0, 1] times [0, 2 pi)}'


def test_typst_Contains():
    x = Symbol('x')
    assert typst(Contains(x, S.Naturals)) == r'x in NN'


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


def test_issue_3568():
    beta = Symbol(r'\beta')
    y = beta + x
    assert typst(y) in [r'"\beta" + x', r'x + "\beta"']

    beta = Symbol(r'beta')
    y = beta + x
    assert typst(y) in [r'beta + x', r'x + beta']


def test_typst():
    assert typst((2*tau)**Rational(7, 2)) == r"8 sqrt(2) tau^(7/2)"
    assert typst((2*mu)**Rational(7, 2), mode='equation*') == \
        r"8 sqrt(2) mu^(7/2)"
    assert typst([2/x, y]) == r"[2/x, #h(0.5em)y]"


def test_typst_dict():
    d = {Rational(1): 1, x**2: 2, x: 3, x**3: 4}
    assert typst(d) == \
        r'{1 : 1, #h(0.5em)x : 3, #h(0.5em)x^(2) : 2, #h(0.5em)x^(3) : 4}'
    D = Dict(d)
    assert typst(D) == \
        r'{1 : 1, #h(0.5em)x : 3, #h(0.5em)x^(2) : 2, #h(0.5em)x^(3) : 4}'


def test_typst_list():
    ll = [Symbol('omega1'), Symbol('a'), Symbol('alpha')]
    assert typst(ll) == r'[omega_(1), #h(0.5em)a, #h(0.5em)alpha]'


def test_typst_NumberSymbols():
    assert typst(S.Catalan) == "G"
    assert typst(S.EulerGamma) == r"gamma"
    assert typst(S.Exp1) == "e"
    assert typst(S.GoldenRatio) == r"phi"
    assert typst(S.Pi) == r"pi"
    assert typst(S.TribonacciConstant) == r'upright("TribonacciConstant")'


def test_typst_rational():
    # tests issue 3973
    assert typst(-Rational(1, 2)) == r"- 1/2"
    assert typst(Rational(-1, 2)) == r"- 1/2"
    assert typst(Rational(1, -2)) == r"- 1/2"
    assert typst(-Rational(-1, 2)) == r"1/2"
    assert typst(-Rational(1, 2)*x) == r"- x/2"
    assert typst(-Rational(1, 2)*x + Rational(-2, 3)*y) == \
        r"- x/2 - (2 y)/3"


def test_typst_inverse():
    # tests issue 4129
    assert typst(1/x) == r"1/x"
    assert typst(1/(x + y)) == r"1/(x + y)"


def test_typst_DiracDelta():
    assert typst(DiracDelta(x)) == r"delta(x)"
    assert typst(DiracDelta(x)**2) == r"(delta(x))^(2)"
    assert typst(DiracDelta(x, 0)) == r"delta(x)"
    assert typst(DiracDelta(x, 5)) == \
        r"delta^((5))(x)"
    assert typst(DiracDelta(x, 5)**2) == \
        r"(delta^((5))(x))^(2)"


def test_typst_Heaviside():
    assert typst(Heaviside(x)) == r"theta(x)"
    assert typst(Heaviside(x)**2) == r"(theta(x))^(2)"


def test_typst_KroneckerDelta():
    assert typst(KroneckerDelta(x, y)) == r"delta_(x y)"
    assert typst(KroneckerDelta(x, y + 1)) == r"delta_(x, y + 1)"
    # issue 6578
    assert typst(KroneckerDelta(x + 1, y)) == r"delta_(y, x + 1)"
    assert typst(Pow(KroneckerDelta(x, y), 2, evaluate=False)) == \
        r"(delta_(x y))^(2)"


def test_typst_LeviCivita():
    assert typst(LeviCivita(x, y, z)) == r"epsilon_(x y z)"
    assert typst(LeviCivita(x, y, z)**2) == \
        r"(epsilon_(x y z))^(2)"
    assert typst(LeviCivita(x, y, z + 1)) == r"epsilon_(x, y, z + 1)"
    assert typst(LeviCivita(x, y + 1, z)) == r"epsilon_(x, y + 1, z)"
    assert typst(LeviCivita(x + 1, y, z)) == r"epsilon_(x + 1, y, z)"



def test_typst_mathieu():
    assert typst(mathieuc(x, y, z)) == r'C(x, y, z)'
    assert typst(mathieus(x, y, z)) == r'S(x, y, z)'
    assert typst(mathieuc(x, y, z)**2) == r'C(x, y, z)^(2)'
    assert typst(mathieus(x, y, z)**2) == r'S(x, y, z)^(2)'
    assert typst(mathieucprime(x, y, z)) == r"C^(prime)(x, y, z)"
    assert typst(mathieusprime(x, y, z)) == r"S^(prime)(x, y, z)"
    assert typst(mathieucprime(x, y, z)**2) == r"C^(prime)(x, y, z)^(2)"
    assert typst(mathieusprime(x, y, z)**2) == r"S^(prime)(x, y, z)^(2)"


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

def test_typst_mul_symbol():
    assert typst(4*4**x, mul_symbol='times') == r"4 times 4^(x)"
    assert typst(4*4**x, mul_symbol='dot') == r"4 dot 4^(x)"
    # assert typst(4*4**x, mul_symbol='ldot') == r"4 #h(0.5em).#h(0.5em) 4^(x)"

    assert typst(4*x, mul_symbol='times') == r"4 times x"
    assert typst(4*x, mul_symbol='dot') == r"4 dot x"
    # assert typst(4*x, mul_symbol='ldot') == r"4 #h(0.5em).#h(0.5em) x"


def test_typst_issue_4381():
    y = 4*4**log(2)
    assert typst(y) == r'4 dot 4^(log(2))'
    assert typst(1/y) == r'1/(4 dot 4^(log(2)))'


def test_typst_issue_4576():
    assert typst(Symbol("beta_13_2")) == r"beta_(13 2)"
    assert typst(Symbol("beta_132_20")) == r"beta_(132 20)"
    assert typst(Symbol("beta_13")) == r"beta_(13)"
    assert typst(Symbol("x_a_b")) == r"x_(a b)"
    assert typst(Symbol("x_1_2_3")) == r"x_(1 2 3)"
    assert typst(Symbol("x_a_b1")) == r"x_(a b1)"
    assert typst(Symbol("x_a_1")) == r"x_(a 1)"
    assert typst(Symbol("x_1_a")) == r"x_(1 a)"
    assert typst(Symbol("x_1^aa")) == r"x^(aa)_(1)"
    assert typst(Symbol("x_1__aa")) == r"x^(aa)_(1)"
    assert typst(Symbol("x_11^a")) == r"x^(a)_(11)"
    assert typst(Symbol("x_11__a")) == r"x^(a)_(11)"
    assert typst(Symbol("x_a_a_a_a")) == r"x_(a a a a)"
    assert typst(Symbol("x_a_a^a^a")) == r"x^(a a)_(a a)"
    assert typst(Symbol("x_a_a__a__a")) == r"x^(a a)_(a a)"
    assert typst(Symbol("alpha_11")) == r"alpha_(11)"
    assert typst(Symbol("alpha_11_11")) == r"alpha_(11 11)"
    assert typst(Symbol("alpha_alpha")) == r"alpha_(alpha)"
    assert typst(Symbol("alpha^aleph")) == r"alpha^(aleph)"
    assert typst(Symbol("alpha__aleph")) == r"alpha^(aleph)"


def test_typst_pow_fraction():
    x = Symbol('x')
    # Testing exp
    assert r'e^(-x)' in typst(exp(-x)/2).replace(' ', '')  # Remove Whitespace

    # Testing e^{-x} in case future changes alter behavior of muls or fracs
    # In particular current output is 1/2 e^(- x) but perhaps this will
    # change to 1/(e^(-x))/2

    # Testing general, non-exp, power
    assert r'3^(-x)' in typst(3**-x/2).replace(' ', '')


def test_typst_noncommutative():
    A, B, C = symbols('A,B,C', commutative=False)

    assert typst(A*B*C**-1) == r"A B C^(-1)"
    assert typst(C**-1*A*B) == r"C^(-1) A B"
    assert typst(A*C**-1*B) == r"A C^(-1) B"


def test_typst_order():
    expr = x**3 + x**2*y + y**4 + 3*x*y**3

    assert typst(expr, order='lex') == r"x^(3) + x^(2) y + 3 x y^(3) + y^(4)"
    assert typst(
        expr, order='rev-lex') == r"y^(4) + 3 x y^(3) + x^(2) y + x^(3)"
    assert typst(expr, order='none') == r"x^(3) + y^(4) + y x^(2) + 3 x y^(3)"


def test_typst_Lambda():
    assert typst(Lambda(x, x + 1)) == r"(x arrow.r.bar x + 1)"
    assert typst(Lambda((x, y), x + 1)) == r"((x, #h(0.5em)y) arrow.r.bar x + 1)"
    assert typst(Lambda(x, x)) == r"(x arrow.r.bar x)"

def test_typst_PolyElement():
    Ruv, u, v = ring("u,v", ZZ)
    Rxyz, x, y, z = ring("x,y,z", Ruv)

    assert typst(x - x) == r"0"
    assert typst(x - 1) == r"x - 1"
    assert typst(x + 1) == r"x + 1"

    assert typst((u**2 + 3*u*v + 1)*x**2*y + u + 1) == \
        r"(u^(2) + 3 u v + 1) x^(2) y + u + 1"
    assert typst((u**2 + 3*u*v + 1)*x**2*y + (u + 1)*x) == \
        r"(u^(2) + 3 u v + 1) x^(2) y + (u + 1) x"
    assert typst((u**2 + 3*u*v + 1)*x**2*y + (u + 1)*x + 1) == \
        r"(u^(2) + 3 u v + 1) x^(2) y + (u + 1) x + 1"
    assert typst((-u**2 + 3*u*v - 1)*x**2*y - (u + 1)*x - 1) == \
        r"-(u^(2) - 3 u v + 1) x^(2) y - (u + 1) x - 1"

    assert typst(-(v**2 + v + 1)*x + 3*u*v + 1) == \
        r"-(v^(2) + v + 1) x + 3 u v + 1"
    assert typst(-(v**2 + v + 1)*x - 3*u*v + 1) == \
        r"-(v^(2) + v + 1) x - 3 u v + 1"


def test_typst_FracElement():
    Fuv, u, v = field("u,v", ZZ)
    Fxyzt, x, y, z, t = field("x,y,z,t", Fuv)

    assert typst(x - x) == r"0"
    assert typst(x - 1) == r"x - 1"
    assert typst(x + 1) == r"x + 1"

    assert typst(x/3) == r"frac(x, 3)"
    assert typst(x/z) == r"frac(x, z)"
    assert typst(x*y/z) == r"frac(x y, z)"
    assert typst(x/(z*t)) == r"frac(x, z t)"
    assert typst(x*y/(z*t)) == r"frac(x y, z t)"

    assert typst((x - 1)/y) == r"frac(x - 1, y)"
    assert typst((x + 1)/y) == r"frac(x + 1, y)"
    assert typst((-x - 1)/y) == r"frac(-x - 1, y)"
    assert typst((x + 1)/(y*z)) == r"frac(x + 1, y z)"
    assert typst(-y/(x + 1)) == r"frac(-y, x + 1)"
    assert typst(y*z/(x + 1)) == r"frac(y z, x + 1)"

    assert typst(((u + 1)*x*y + 1)/((v - 1)*z - 1)) == \
        r"frac((u + 1) x y + 1, (v - 1) z - 1)"
    assert typst(((u + 1)*x*y + 1)/((v - 1)*z - t*u*v - 1)) == \
        r"frac((u + 1) x y + 1, (v - 1) z - u v t - 1)"


def test_typst_Poly():
    assert typst(Poly(x**2 + 2 * x, x)) == \
        r'upright("Poly")(x^(2) + 2 x, x, domain=ZZ)'
    assert typst(Poly(x/y, x)) == \
        r'upright("Poly")(1/y x, x, domain=ZZ(y))'
    assert typst(Poly(2.0*x + y)) == \
        r'upright("Poly")(2.0 x + 1.0 y, x, y, domain=RR)'


def test_typst_Poly_order():
    assert typst(Poly([a, 1, b, 2, c, 3], x)) == \
        r'upright("Poly")(a x^(5) + x^(4) + b x^(3) + 2 x^(2) + c x + 3, x, domain=ZZ[a, b, c])'
    assert typst(Poly([a, 1, b+c, 2, 3], x)) == \
        r'upright("Poly")(a x^(4) + x^(3) + (b + c) x^(2) + 2 x + 3, x, domain=ZZ[a, b, c])'
    assert typst(Poly(a*x**3 + x**2*y - x*y - c*y**3 - b*x*y**2 + y - a*x + b,
                      (x, y))) == \
        r'upright("Poly")(a x^(3) + x^(2)y -  b xy^(2) - xy -  a x -  c y^(3) + y + b, x, y, domain=ZZ[a, b, c])'


def test_typst_ComplexRootOf():
    assert typst(rootof(x**5 + x + 3, 0)) == \
        r'upright("CRootOf")(x^(5) + x + 3, 0)'


def test_typst_RootSum():
    assert typst(RootSum(x**5 + x + 3, sin)) == \
        r'upright("RootSum")(x^(5) + x + 3, (x arrow.r.bar sin(x)))'


def test_settings():
    raises(TypeError, lambda: typst(x*y, method="garbage"))


def test_typst_numbers():
    assert typst(catalan(n)) == r'G_(n)'
    assert typst(catalan(n)**2) == r'G_(n)^(2)'
    assert typst(bernoulli(n)) == r'B_(n)'
    assert typst(bernoulli(n, x)) == r'B_(n)(x)'
    assert typst(bernoulli(n)**2) == r'B_(n)^(2)'
    assert typst(bernoulli(n, x)**2) == r'B_(n)^(2)(x)'
    assert typst(genocchi(n)) == r'G_(n)'
    assert typst(genocchi(n, x)) == r'G_(n)(x)'
    assert typst(genocchi(n)**2) == r'G_(n)^(2)'
    assert typst(genocchi(n, x)**2) == r'G_(n)^(2)(x)'
    assert typst(bell(n)) == r'B_(n)'
    assert typst(bell(n, x)) == r'B_(n)(x)'
    assert typst(bell(n, m, (x, y))) == r'B_(n, m)(x, y)'
    assert typst(bell(n)**2) == r'B_(n)^(2)'
    assert typst(bell(n, x)**2) == r'B_(n)^(2)(x)'
    assert typst(bell(n, m, (x, y))**2) == r'B_(n, m)^(2)(x, y)'
    assert typst(fibonacci(n)) == r'F_(n)'
    assert typst(fibonacci(n, x)) == r'F_(n)(x)'
    assert typst(fibonacci(n)**2) == r'F_(n)^(2)'
    assert typst(fibonacci(n, x)**2) == r'F_(n)^(2)(x)'
    assert typst(lucas(n)) == r'L_(n)'
    assert typst(lucas(n)**2) == r'L_(n)^(2)'
    assert typst(tribonacci(n)) == r'T_(n)'
    assert typst(tribonacci(n, x)) == r'T_(n)(x)'
    assert typst(tribonacci(n)**2) == r'T_(n)^(2)'
    assert typst(tribonacci(n, x)**2) == r'T_(n)^(2)(x)'
    assert typst(mobius(n)) == r'mu(n)'
    assert typst(mobius(n)**2) == r'mu^(2)(n)'


def test_typst_euler():
    assert typst(euler(n)) == r'E_(n)'
    assert typst(euler(n, x)) == r'E_(n)(x)'
    assert typst(euler(n, x)**2) == r'E_(n)^(2)(x)'


def test_lamda():
    assert typst(Symbol('lamda')) == r'lambda'
    assert typst(Symbol('Lamda')) == r'Lambda'


def test_custom_symbol_names():
    x = Symbol('x')
    y = Symbol('y')
    assert typst(x) == r'x'
    assert typst(x, symbol_names={x: "x_i"}) == r'x_i'
    assert typst(x + y, symbol_names={x: "x_i"}) == r'x_i + y'
    assert typst(x**2, symbol_names={x: "x_i"}) == r'x_i^(2)'
    assert typst(x + y, symbol_names={x: "x_i", y: "y_j"}) == r'x_i + y_j'

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


def test_typst_RandomDomain():
    from sympy.stats import Normal, Die, Exponential, pspace, where
    from sympy.stats.rv import RandomDomain

    X = Normal('x1', 0, 1)
    assert typst(where(X > 0)) == r'upright("Domain: ")0 < x_(1) and x_(1) < infinity'

    D = Die('d1', 6)
    assert typst(where(D > 4)) == r'upright("Domain: ")d_(1) = 5 or d_(1) = 6'

    A = Exponential('a', 1)
    B = Exponential('b', 1)
    assert typst(
        pspace(Tuple(A, B)).domain) == \
        r'upright("Domain: ")0 <= a and 0 <= b and a < infinity and b < infinity'

    assert typst(RandomDomain(FiniteSet(x), FiniteSet(1, 2))) == \
        r'upright("Domain: "){x} in {1, 2}'

def test_typst_PrettyPoly():
    from sympy.polys.domains import QQ
    F = QQ.frac_field(x, y)
    R = QQ[x, y]

    assert typst(F.convert(x/(x + y))) == r'frac(x, x + y)'
    assert typst(R.convert(x + y)) == typst(x + y)


def test_typst_integral_transforms():
    x = Symbol("x")
    k = Symbol("k")
    f = Function("f")
    a = Symbol("a")
    b = Symbol("b")

    assert typst(MellinTransform(f(x), x, k)) == \
        r'cal(M)_(x)[f(x)](k)'
    assert typst(InverseMellinTransform(f(k), k, x, a, b)) == \
        r'cal(M)^(-1)_(k)[f(k)](x)'

    assert typst(LaplaceTransform(f(x), x, k)) == \
        r'cal(L)_(x)[f(x)](k)'
    assert typst(InverseLaplaceTransform(f(k), k, x, (a, b))) == \
        r'cal(L)^(-1)_(k)[f(k)](x)'

    assert typst(FourierTransform(f(x), x, k)) == \
        r'cal(F)_(x)[f(x)](k)'
    assert typst(InverseFourierTransform(f(k), k, x)) == \
        r'cal(F)^(-1)_(k)[f(k)](x)'

    assert typst(CosineTransform(f(x), x, k)) == \
        r'cal(COS)_(x)[f(x)](k)'
    assert typst(InverseCosineTransform(f(k), k, x)) == \
        r'cal(COS)^(-1)_(k)[f(k)](x)'

    assert typst(SineTransform(f(x), x, k)) == \
        r'cal(SIN)_(x)[f(x)](k)'
    assert typst(InverseSineTransform(f(k), k, x)) == \
        r'cal(SIN)^(-1)_(k)[f(k)](x)'


def test_typst_PolynomialRingBase():
    from sympy.polys.domains import QQ
    assert typst(QQ.old_poly_ring(x, y)) == r'QQ[x, y]'
    assert typst(QQ.old_poly_ring(x, y, order="ilex")) == \
        r'S_<^(-1)QQ[x, y]'


def test_typst_categories():
    from sympy.categories import (Object, IdentityMorphism,
                                  NamedMorphism, Category, Diagram,
                                  DiagramGrid)

    A1 = Object("A1")
    A2 = Object("A2")
    A3 = Object("A3")

    f1 = NamedMorphism(A1, A2, "f1")
    f2 = NamedMorphism(A2, A3, "f2")
    id_A1 = IdentityMorphism(A1)

    K1 = Category("K1")

    assert typst(A1) == r'A_(1)'
    assert typst(f1) == r'f_(1):A_(1) arrow.r.bar A_(2)'
    assert typst(id_A1) == r'"id":A_(1) arrow.r.bar A_(1)'
    assert typst(f2*f1) == r'f_(2)circle.small f_(1):A_(1) arrow.r.bar A_(3)'

    assert typst(K1) == r'upright(bold(K_(1)))'

    d = Diagram()
    assert typst(d) == r'nothing'

    d = Diagram({f1: "unique", f2: S.EmptySet})
    assert typst(d) == r'{f_(2)circle.small f_(1):A_(1) arrow.r.bar A_(3) : ' \
        r'nothing, #h(0.5em)"id":A_(1) arrow.r.bar A_(1) : nothing, #h(0.5em)' \
        r'"id":A_(2) arrow.r.bar A_(2) : nothing, #h(0.5em)"id":A_(3) arrow.r.bar ' \
        r'A_(3) : nothing, #h(0.5em)f_(1):A_(1) arrow.r.bar A_(2) : {"unique"}, #h(0.5em)f_(2):A_(2) arrow.r.bar A_(3) : nothing}'

    d = Diagram({f1: "unique", f2: S.EmptySet}, {f2 * f1: "unique"})
    assert typst(d) == r'{f_(2)circle.small f_(1):A_(1) arrow.r.bar A_(3) : ' \
        r'nothing, #h(0.5em)"id":A_(1) arrow.r.bar A_(1) : nothing, #h(0.5em)' \
        r'"id":A_(2) arrow.r.bar A_(2) : nothing, #h(0.5em)"id":A_(3) arrow.r.bar A_(3) ' \
        r': nothing, #h(0.5em)f_(1):A_(1) arrow.r.bar A_(2) : {"unique"}, #h(0.5em)' \
        r'f_(2):A_(2) arrow.r.bar A_(3) : nothing}arrow.r.long {f_(2)circle.small f_(1):A_(1) arrow.r.bar A_(3) : {"unique"}}'

    # A linear diagram.
    A = Object("A")
    B = Object("B")
    C = Object("C")
    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, C, "g")
    d = Diagram([f, g])
    grid = DiagramGrid(d)

    assert typst(grid) == r'mat(delim:#none, A, B;, C)'


def test_typst_Modules():
    from sympy.polys.domains import QQ
    from sympy.polys.agca import homomorphism

    R = QQ.old_poly_ring(x, y)
    F = R.free_module(2)
    M = F.submodule([x, y], [1, x**2])

    assert typst(F) == r'(QQ[x, y])^(2)'
    assert typst(M) == \
        r'angle.l [x, y], [1, x^(2)] angle.r'

    I = R.ideal(x**2, y)
    assert typst(I) == r'angle.l x^(2), y angle.r'

    Q = F / M
    assert typst(Q) == \
        r'frac((QQ[x, y])^(2), angle.l [x, y], [1, x^(2)] angle.r)'
    assert typst(Q.submodule([1, x**3/2], [2, y])) == \
        r'angle.l [1,x^(3)/2] + angle.l [x, y], [1, x^(2)] angle.r, [2,y] + angle.l [x, y], [1, x^(2)] angle.r angle.r'

    h = homomorphism(QQ.old_poly_ring(x).free_module(2),
                     QQ.old_poly_ring(x).free_module(2), [0, 0])

    assert typst(h) == \
        r'mat(delim: "[", 0, 0; 0, 0) : (QQ[x])^(2) arrow.r.bar (QQ[x])^(2)'


def test_typst_QuotientRing():
    from sympy.polys.domains import QQ
    R = QQ.old_poly_ring(x)/[x**2 + 1]

    assert typst(R) == \
        r'frac(QQ[x], angle.l x^(2) + 1 angle.r)'
    assert typst(R.one) == r'1 + angle.l x^(2) + 1 angle.r'


def test_typst_Tr():
    #TODO: Handle indices
    A, B = symbols('A B', commutative=False)
    t = Tr(A*B)
    assert typst(t) == r'upright("tr")(A B)'


def test_typst_Determinant():
    from sympy.matrices import Determinant, Inverse, BlockMatrix, OneMatrix, ZeroMatrix
    m = Matrix(((1, 2), (3, 4)))
    assert typst(Determinant(m)) == r'mat(delim: "|", 1, 2; 3, 4)'
    assert typst(Determinant(Inverse(m))) == \
        r'mat(delim: "|", mat(delim: "[", 1, 2; 3, 4)^(-1))'
    X = MatrixSymbol('X', 2, 2)
    assert typst(Determinant(X)) == r'mat(delim: "|", X)'
    assert typst(Determinant(X + m)) == \
        r'mat(delim: "|", mat(delim: "[", 1, 2; 3, 4) + X)'
    assert typst(Determinant(BlockMatrix(((OneMatrix(2, 2), X),
                                          (m, ZeroMatrix(2, 2)))))) == \
        r'mat(delim: "|", 1, X; mat(delim: "[", 1, 2; 3, 4), 0)'

def test_Adjoint():
    from sympy.matrices import Adjoint, Inverse, Transpose
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    assert typst(Adjoint(X)) == r'X^(dagger)'
    assert typst(Adjoint(X + Y)) == r'(X + Y)^(dagger)'
    assert typst(Adjoint(X) + Adjoint(Y)) == r'X^(dagger) + Y^(dagger)'
    assert typst(Adjoint(X*Y)) == r'(X Y)^(dagger)'
    assert typst(Adjoint(Y)*Adjoint(X)) == r'Y^(dagger) X^(dagger)'
    assert typst(Adjoint(X**2)) == r'(X^(2))^(dagger)'
    assert typst(Adjoint(X)**2) == r'(X^(dagger))^(2)'
    assert typst(Adjoint(Inverse(X))) == r'(X^(-1))^(dagger)'
    assert typst(Inverse(Adjoint(X))) == r'(X^(dagger))^(-1)'
    assert typst(Adjoint(Transpose(X))) == r'(X^(T))^(dagger)'
    assert typst(Transpose(Adjoint(X))) == r'(X^(dagger))^(T)'
    assert typst(Transpose(Adjoint(X) + Y)) == r'(X^(dagger) + Y)^(T)'
    m = Matrix(((1, 2), (3, 4)))
    assert typst(Adjoint(m)) == r'mat(delim: "[", 1, 2; 3, 4)^(dagger)'
    assert typst(Adjoint(m+X)) == r'(mat(delim: "[", 1, 2; 3, 4) + X)^(dagger)'
    from sympy.matrices import BlockMatrix, OneMatrix, ZeroMatrix
    assert typst(Adjoint(BlockMatrix(((OneMatrix(2, 2), X),
                                      (m, ZeroMatrix(2, 2)))))) == \
        r'mat(delim: "[", 1, X; mat(delim: "[", 1, 2; 3, 4), 0)^(dagger)'
    # Issue 20959
    Mx = MatrixSymbol('M^x', 2, 2)
    assert typst(Adjoint(Mx)) == r'(M^(x))^(dagger)'

    # adjoint style
    assert typst(Adjoint(X), adjoint_style="star") == r'X^(star)'
    assert typst(Adjoint(X + Y), adjoint_style="hermitian") == r'(X + Y)^(sans(upright(H)))'
    assert typst(Adjoint(X) + Adjoint(Y), adjoint_style="dagger") == r'X^(dagger) + Y^(dagger)'
    assert typst(Adjoint(Y)*Adjoint(X)) == r'Y^(dagger) X^(dagger)'
    assert typst(Adjoint(X**2), adjoint_style="star") == r'(X^(2))^(star)'
    assert typst(Adjoint(X)**2, adjoint_style="hermitian") == r'(X^(sans(upright(H))))^(2)'

def test_Transpose():
    from sympy.matrices import Transpose, MatPow, HadamardPower
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    assert typst(Transpose(X)) == r'X^(T)'
    assert typst(Transpose(X + Y)) == r'(X + Y)^(T)'

    assert typst(Transpose(HadamardPower(X, 2))) == r'(X^(circle.small 2))^(T)'
    assert typst(HadamardPower(Transpose(X), 2)) == r'(X^(T))^(circle.small 2)'
    assert typst(Transpose(MatPow(X, 2))) == r'(X^(2))^(T)'
    assert typst(MatPow(Transpose(X), 2)) == r'(X^(T))^(2)'
    m = Matrix(((1, 2), (3, 4)))
    assert typst(Transpose(m)) == r'mat(delim: "[", 1, 2; 3, 4)^(T)'
    assert typst(Transpose(m+X)) == r'(mat(delim: "[", 1, 2; 3, 4) + X)^(T)'
    from sympy.matrices import BlockMatrix, OneMatrix, ZeroMatrix
    assert typst(Transpose(BlockMatrix(((OneMatrix(2, 2), X),
                                        (m, ZeroMatrix(2, 2)))))) == \
        r'mat(delim: "[", 1, X; mat(delim: "[", 1, 2; 3, 4), 0)^(T)'
    # Issue 20959
    Mx = MatrixSymbol('M^x', 2, 2)
    assert typst(Transpose(Mx)) == r'(M^(x))^(T)'


def test_Hadamard():
    from sympy.matrices import HadamardProduct, HadamardPower
    from sympy.matrices.expressions import MatAdd, MatMul, MatPow
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    assert typst(HadamardProduct(X, Y*Y)) == r'X circle.small Y^(2)'
    assert typst(HadamardProduct(X, Y)*Y) == r'(X circle.small Y) Y'

    assert typst(HadamardPower(X, 2)) == r'X^(circle.small 2)'
    assert typst(HadamardPower(X, -1)) == r'X^(circle.small (-1))'
    assert typst(HadamardPower(MatAdd(X, Y), 2)) == r'(X + Y)^(circle.small 2)'
    assert typst(HadamardPower(MatMul(X, Y), 2)) == r'(X Y)^(circle.small 2)'

    assert typst(HadamardPower(MatPow(X, -1), -1)) == r'(X^(-1))^(circle.small (-1))'
    assert typst(MatPow(HadamardPower(X, -1), -1)) == r'(X^(circle.small (-1)))^(-1)'

    assert typst(HadamardPower(X, n+1)) == r'X^(circle.small (n + 1))'


def test_MatPow():
    from sympy.matrices.expressions import MatPow
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    assert typst(MatPow(X, 2)) == r'X^(2)'
    assert typst(MatPow(X*X, 2)) == r'(X^(2))^(2)'
    assert typst(MatPow(X*Y, 2)) == r'(X Y)^(2)'
    assert typst(MatPow(X + Y, 2)) == r'(X + Y)^(2)'
    assert typst(MatPow(X + X, 2)) == r'(2 X)^(2)'
    # Issue 20959
    Mx = MatrixSymbol('M^x', 2, 2)
    assert typst(MatPow(Mx, 2)) == r'(M^(x))^(2)'


def test_ElementwiseApplyFunction():
    X = MatrixSymbol('X', 2, 2)
    expr = (X.T*X).applyfunc(sin)
    assert typst(expr) == r'(d arrow.r.bar sin(d))_circle.stroked.small (X^(T) X)'
    expr = X.applyfunc(Lambda(x, 1/x))
    assert typst(expr) == r'(x arrow.r.bar 1/x)_circle.stroked.small (X)'


def test_ZeroMatrix():
    from sympy.matrices.expressions.special import ZeroMatrix
    assert typst(ZeroMatrix(1, 1), mat_symbol_style='plain') == r"0"
    assert typst(ZeroMatrix(1, 1), mat_symbol_style='bold') == r"bold(0)"


def test_OneMatrix():
    from sympy.matrices.expressions.special import OneMatrix
    assert typst(OneMatrix(3, 4), mat_symbol_style='plain') == r"1"
    assert typst(OneMatrix(3, 4), mat_symbol_style='bold') == r"bold(1)"


def test_Identity():
    from sympy.matrices.expressions.special import Identity
    assert typst(Identity(1), mat_symbol_style='plain') == r"II"
    assert typst(Identity(1), mat_symbol_style='bold') == r"bold(I)"


def test_typst_DFT_IDFT():
    from sympy.matrices.expressions.fourier import DFT, IDFT
    assert typst(DFT(13)) == r'upright("DFT")_(13)'
    assert typst(IDFT(x)) == r'upright("IDFT")_(x)'

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

def test_imaginary():
    i = sqrt(-1)
    assert typst(i) == r'i'


def test_builtins_without_args():
    assert typst(sin) == r'sin'
    assert typst(cos) == r'cos'
    assert typst(tan) == r'tan'
    assert typst(log) == r'log'
    assert typst(Ei) == r'upright("Ei")'
    assert typst(zeta) == r'zeta'


def test_typst_greek_functions():
    # bug because capital greeks that have roman equivalents should not use
    # Alpha, Beta, Eta, etc.
    s = Function('Alpha')
    assert typst(s) == r'Alpha'
    assert typst(s(x)) == r'Alpha(x)'
    s = Function('Beta')
    assert typst(s) == r'Beta'
    s = Function('Eta')
    assert typst(s) == r'Eta'
    assert typst(s(x)) == r'Eta(x)'

    # bug because sympy.core.numbers.Pi is special
    p = Function('Pi')
    # assert typst(p(x)) == r'Pi(x)'
    assert typst(p) == r'Pi'

    # bug because not all greeks are included
    c = Function('chi')
    assert typst(c(x)) == r'chi(x)'
    assert typst(c) == r'chi'




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
        r'{(x, #h(0.5em)a) mid(|) (x, #h(0.5em)a) in CC^(2) and sin(a x) = 0 and cos(a x) = 0}'

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
