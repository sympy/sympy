from sympy.core.add import Add
from sympy.core.function import Function
from sympy.core.mul import Mul
from sympy.core.numbers import (I, pi, Rational, oo, Float)
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.core.numbers import Number
from sympy.functions.elementary.exponential import (exp, log)
from sympy.functions.special.delta_functions import Heaviside
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import atan
from sympy.functions.elementary.complexes import re, im
from sympy.matrices.dense import eye
from sympy.physics.control.lti import SISOLinearTimeInvariant
from sympy.polys.polytools import factor
from sympy.polys.rootoftools import CRootOf
from sympy.simplify.simplify import simplify
from sympy.core.containers import Tuple
from sympy.matrices import ImmutableMatrix, Matrix, ShapeError
from sympy.functions.elementary.trigonometric import sin, cos
from sympy.physics.control.lti import (
    create_transfer_function, TransferFunctionBase, TransferFunction,
    DiscreteTransferFunction, PIDController,Series, Parallel, Feedback,
    TransferFunctionMatrix, MIMOSeries, MIMOParallel, MIMOFeedback, StateSpace,
    DiscreteStateSpace, create_state_space, gbt, bilinear, forward_diff,
    backward_diff, phase_margin, gain_margin)
from sympy.testing.pytest import raises
from sympy.logic.boolalg import false, true

from math import isclose

a, x, b, c, s, g, d, p, k, tau, zeta, wn, T, z = symbols('a, x, b, c, s, g, d,\
    p, k, tau, zeta, wn, T, z')
a0, a1, a2, a3, b0, b1, b2, b3, b4, c0, c1, c2, c3, d0, d1, d2, d3 = \
    symbols('a0:4, b0:5, c0:4, d0:4')

TF1 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
TF2 = TransferFunction(k, 1, s)
TF3 = TransferFunction(a2*p - s, a2*s + p, s)


def test_create_transfer_function():
    cont_tf1 = create_transfer_function(s+1, s**2 + 2, s)
    assert isinstance(cont_tf1, TransferFunction)

    cont_tf2 = create_transfer_function(s+1, s**2 + 2, s, sampling_time = 0)
    assert isinstance(cont_tf2, TransferFunction)

    disc_tf1 = create_transfer_function(z, z + 1, z, 0.1)
    assert isinstance(disc_tf1, DiscreteTransferFunction)

    disc_tf2 = create_transfer_function(z, z + 1, z, T)
    assert isinstance(disc_tf2, DiscreteTransferFunction)

def test_TransferFunctionBase():
    raises(NotImplementedError,
           lambda: TransferFunctionBase.from_rational_expression(2 / s))
    raises(NotImplementedError,
        lambda: TransferFunctionBase.from_coeff_lists([1], [1], s))
    raises(NotImplementedError,
        lambda: TransferFunctionBase.from_zpk([1], [1], 1, s))

def test_TransferFunction_construction():
    tf = TransferFunction(s + 1, s**2 + s + 1, s)
    assert tf.num == (s + 1)
    assert tf.den == (s**2 + s + 1)
    assert tf.args == (s + 1, s**2 + s + 1, s)

    tf1 = TransferFunction(s + 4, s - 5, s)
    assert tf1.num == (s + 4)
    assert tf1.den == (s - 5)
    assert tf1.args == (s + 4, s - 5, s)

    # using different polynomial variables.
    tf2 = TransferFunction(p + 3, p**2 - 9, p)
    assert tf2.num == (p + 3)
    assert tf2.den == (p**2 - 9)
    assert tf2.args == (p + 3, p**2 - 9, p)

    tf3 = TransferFunction(p**3 + 5*p**2 + 4, p**4 + 3*p + 1, p)
    assert tf3.args == (p**3 + 5*p**2 + 4, p**4 + 3*p + 1, p)

    # no pole-zero cancellation on its own.
    tf4 = TransferFunction((s + 3)*(s - 1), (s - 1)*(s + 5), s)
    assert tf4.den == (s - 1)*(s + 5)
    assert tf4.args == ((s + 3)*(s - 1), (s - 1)*(s + 5), s)

    tf4_ = TransferFunction(p + 2, p + 2, p)
    assert tf4_.args == (p + 2, p + 2, p)

    tf5 = TransferFunction(s - 1, 4 - p, s)
    assert tf5.args == (s - 1, 4 - p, s)

    tf5_ = TransferFunction(s - 1, s - 1, s)
    assert tf5_.args == (s - 1, s - 1, s)

    tf6 = TransferFunction(5, 6, s)
    assert tf6.num == 5
    assert tf6.den == 6
    assert tf6.args == (5, 6, s)

    tf6_ = TransferFunction(1/2, 4, s)
    assert tf6_.num == 0.5
    assert tf6_.den == 4
    assert tf6_.args == (0.500000000000000, 4, s)

    tf7 = TransferFunction(3*s**2 + 2*p + 4*s, 8*p**2 + 7*s, s)
    tf8 = TransferFunction(3*s**2 + 2*p + 4*s, 8*p**2 + 7*s, p)
    assert not tf7 == tf8

    tf7_ = TransferFunction(a0*s + a1*s**2 + a2*s**3, b0*p - b1*s, s)
    tf8_ = TransferFunction(a0*s + a1*s**2 + a2*s**3, b0*p - b1*s, s)
    assert tf7_ == tf8_
    assert -(-tf7_) == tf7_ == -(-(-(-tf7_)))

    tf9 = TransferFunction(a*s**3 + b*s**2 + g*s + d, d*p + g*p**2 + g*s, s)
    assert tf9.args == (a*s**3 + b*s**2 + d + g*s, d*p + g*p**2 + g*s, s)

    tf10 = TransferFunction(p**3 + d, g*s**2 + d*s + a, p)
    tf10_ = TransferFunction(p**3 + d, g*s**2 + d*s + a, p)
    assert tf10.args == (d + p**3, a + d*s + g*s**2, p)
    assert tf10_ == tf10

    tf11 = TransferFunction(a1*s + a0, b2*s**2 + b1*s + b0, s)
    assert tf11.num == (a0 + a1*s)
    assert tf11.den == (b0 + b1*s + b2*s**2)
    assert tf11.args == (a0 + a1*s, b0 + b1*s + b2*s**2, s)

    # when just the numerator is 0, leave the denominator alone.
    tf12 = TransferFunction(0, p**2 - p + 1, p)
    assert tf12.args == (0, p**2 - p + 1, p)

    tf13 = TransferFunction(0, 1, s)
    assert tf13.args == (0, 1, s)

    # float exponents
    tf14 = TransferFunction(a0*s**0.5 + a2*s**0.6 - a1, a1*p**(-8.7), s)
    assert tf14.args == (a0*s**0.5 - a1 + a2*s**0.6, a1*p**(-8.7), s)

    tf15 = TransferFunction(a2**2*p**(1/4) + a1*s**(-4/5), a0*s - p, p)
    assert tf15.args == (a1*s**(-0.8) + a2**2*p**0.25, a0*s - p, p)

    omega_o, k_p, k_o, k_i = symbols('omega_o, k_p, k_o, k_i')
    tf18 = TransferFunction((k_p + k_o*s + k_i/s),
                            s**2 + 2*omega_o*s + omega_o**2, s)
    assert tf18.num == k_i/s + k_o*s + k_p
    assert tf18.args == (k_i/s + k_o*s + k_p,
                         omega_o**2 + 2*omega_o*s + s**2, s)

    # ValueError when denominator is zero.
    raises(ValueError, lambda: TransferFunction(4, 0, s))
    raises(ValueError, lambda: TransferFunction(s, 0, s))
    raises(ValueError, lambda: TransferFunction(0, 0, s))

    raises(TypeError, lambda: TransferFunction(Matrix([1, 2, 3]), s, s))

    raises(TypeError, lambda: TransferFunction(s**2 + 2*s - 1, s + 3, 3))
    raises(TypeError, lambda: TransferFunction(p + 1, 5 - p, 4))
    raises(TypeError, lambda: TransferFunction(3, 4, 8))


def test_TransferFunction_functions():
    # classmethod from_rational_expression
    expr_1 = Mul(0, Pow(s, -1, evaluate=False), evaluate=False)
    expr_2 = s/0
    expr_3 = (p*s**2 + 5*s)/(s + 1)**3
    expr_4 = 6
    expr_5 = ((2 + 3*s)*(5 + 2*s))/((9 + 3*s)*(5 + 2*s**2))
    expr_6 = (9*s**4 + 4*s**2 + 8)/((s + 1)*(s + 9))
    tf = TransferFunction(s + 1, s**2 + 2, s)
    delay = exp(-s/tau)
    expr_7 = delay*tf.to_expr()
    H1 = TransferFunction.from_rational_expression(expr_7, s)
    H2 = TransferFunction(s + 1, (s**2 + 2)*exp(s/tau), s)
    expr_8 = Add(2,  3*s/(s**2 + 1), evaluate=False)

    assert TransferFunction.from_rational_expression(expr_1) == TransferFunction(0, s, s)
    raises(ZeroDivisionError, lambda: TransferFunction.from_rational_expression(expr_2))
    raises(ValueError, lambda: TransferFunction.from_rational_expression(expr_3))
    assert TransferFunction.from_rational_expression(expr_3, s) == TransferFunction((p*s**2 + 5*s), (s + 1)**3, s)
    assert TransferFunction.from_rational_expression(expr_3, p) == TransferFunction((p*s**2 + 5*s), (s + 1)**3, p)
    raises(ValueError, lambda: TransferFunction.from_rational_expression(expr_4))
    assert TransferFunction.from_rational_expression(expr_4, s) == TransferFunction(6, 1, s)
    assert TransferFunction.from_rational_expression(expr_5, s) == \
        TransferFunction((2 + 3*s)*(5 + 2*s), (9 + 3*s)*(5 + 2*s**2), s)
    assert TransferFunction.from_rational_expression(expr_6, s) == \
        TransferFunction((9*s**4 + 4*s**2 + 8), (s + 1)*(s + 9), s)
    assert H1 == H2
    assert TransferFunction.from_rational_expression(expr_8, s) == \
        TransferFunction(2*s**2 + 3*s + 2, s**2 + 1, s)

    # classmethod from_coeff_lists
    tf1 = TransferFunction.from_coeff_lists([1, 2], [3, 4, 5], s)
    num2 = [p**2, 2*p]
    den2 = [p**3, p + 1, 4]
    tf2 = TransferFunction.from_coeff_lists(num2, den2, s)
    num3 = [1, 2, 3]
    den3 = [0, 0]

    assert tf1 == TransferFunction(s + 2, 3*s**2 + 4*s + 5, s)
    assert tf2 == TransferFunction(p**2*s + 2*p, p**3*s**2 + s*(p + 1) + 4, s)
    raises(ZeroDivisionError, lambda: TransferFunction.from_coeff_lists(num3, den3, s))

    # classmethod from_zpk
    zeros = [4]
    poles = [-1+2j, -1-2j]
    gain = 3
    tf1 = TransferFunction.from_zpk(zeros, poles, gain, s)

    assert tf1 == TransferFunction(3*s - 12, (s + 1.0 - 2.0*I)*(s + 1.0 + 2.0*I), s)

    # explicitly cancel poles and zeros.
    tf0 = TransferFunction(s**5 + s**3 + s, s - s**2, s)
    a = TransferFunction(-(s**4 + s**2 + 1), s - 1, s)
    assert tf0.simplify() == simplify(tf0) == a

    tf1 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
    b = TransferFunction(p + 3, p + 5, p)
    assert tf1.simplify() == simplify(tf1) == b

    # expand the numerator and the denominator.
    G1 = TransferFunction((1 - s)**2, (s**2 + 1)**2, s)
    G2 = TransferFunction(1, -3, p)
    c = (a2*s**p + a1*s**s + a0*p**p)*(p**s + s**p)
    d = (b0*s**s + b1*p**s)*(b2*s*p + p**p)
    e = a0*p**p*p**s + a0*p**p*s**p + a1*p**s*s**s + a1*s**p*s**s + a2*p**s*s**p + a2*s**(2*p)
    f = b0*b2*p*s*s**s + b0*p**p*s**s + b1*b2*p*p**s*s + b1*p**p*p**s
    g = a1*a2*s*s**p + a1*p*s + a2*b1*p*s*s**p + b1*p**2*s
    G3 = TransferFunction(c, d, s)
    G4 = TransferFunction(a0*s**s - b0*p**p, (a1*s + b1*s*p)*(a2*s**p + p), p)

    assert G1.expand() == TransferFunction(s**2 - 2*s + 1, s**4 + 2*s**2 + 1, s)
    assert tf1.expand() == TransferFunction(p**2 + 2*p - 3, p**2 + 4*p - 5, p)
    assert G2.expand() == G2
    assert G3.expand() == TransferFunction(e, f, s)
    assert G4.expand() == TransferFunction(a0*s**s - b0*p**p, g, p)

    # testing that subs works.
    p1 = a1*s + a0
    p2 = b2*s**2 + b1*s + b0
    SP1 = TransferFunction(p1, p2, s)
    expect1 = TransferFunction(2.0*s + 1.0, 5.0*s**2 + 4.0*s + 3.0, s)
    expect1_ = TransferFunction(2*s + 1, 5*s**2 + 4*s + 3, s)
    assert SP1.subs({a0: 1, a1: 2, b0: 3, b1: 4, b2: 5}) == expect1_
    assert SP1.subs({a0: 1, a1: 2, b0: 3, b1: 4, b2: 5}).evalf() == expect1
    assert expect1_.evalf() == expect1

    c1, d0, d1, d2 = symbols('c1, d0:3')
    p3, p4 = c1*p, d2*p**3 + d1*p**2 - d0
    SP2 = TransferFunction(p3, p4, p)
    expect2 = TransferFunction(2.0*p, 5.0*p**3 + 2.0*p**2 - 3.0, p)
    expect2_ = TransferFunction(2*p, 5*p**3 + 2*p**2 - 3, p)
    assert SP2.subs({c1: 2, d0: 3, d1: 2, d2: 5}) == expect2_
    assert SP2.subs({c1: 2, d0: 3, d1: 2, d2: 5}).evalf() == expect2
    assert expect2_.evalf() == expect2

    SP3 = TransferFunction(a0*p**3 + a1*s**2 - b0*s + b1, a1*s + p, s)
    expect3 = TransferFunction(2.0*p**3 + 4.0*s**2 - s + 5.0, p + 4.0*s, s)
    expect3_ = TransferFunction(2*p**3 + 4*s**2 - s + 5, p + 4*s, s)
    assert SP3.subs({a0: 2, a1: 4, b0: 1, b1: 5}) == expect3_
    assert SP3.subs({a0: 2, a1: 4, b0: 1, b1: 5}).evalf() == expect3
    assert expect3_.evalf() == expect3

    SP4 = TransferFunction(s - a1*p**3, a0*s + p, p)
    expect4 = TransferFunction(7.0*p**3 + s, p - s, p)
    expect4_ = TransferFunction(7*p**3 + s, p - s, p)
    assert SP4.subs({a0: -1, a1: -7}) == expect4_
    assert SP4.subs({a0: -1, a1: -7}).evalf() == expect4
    assert expect4_.evalf() == expect4

    # evaluate the transfer function at particular frequencies.
    assert simplify(tf1.eval_frequency(wn) - (wn**2/(wn**2 + 4*wn - 5) +
                                              2*wn/(wn**2 + 4*wn - 5) -
                                              3/(wn**2 + 4*wn - 5))) == 0
    assert G1.eval_frequency(1 + I) == S(3)/25 + S(4)*I/25
    assert G4.eval_frequency(S(5)/3) == \
        a0*s**s/(a1*a2*s**(S(8)/3) + S(5)*a1*s/3 + 5*a2*b1*s**(S(8)/3)/3 + S(25)*b1*s/9) - 5*3**(S(1)/3)*5**(S(2)/3)*b0/(9*a1*a2*s**(S(8)/3) + 15*a1*s + 15*a2*b1*s**(S(8)/3) + 25*b1*s)

    # Low-frequency (or DC) gain.
    assert tf0.dc_gain() == 1
    assert tf1.dc_gain() == Rational(3, 5)
    assert SP2.dc_gain() == 0
    assert expect4.dc_gain() == -1
    assert expect2_.dc_gain() == 0
    assert TransferFunction(1, s, s).dc_gain() == oo

    # Poles of a transfer function.
    tf_ = TransferFunction(x**3 - k, k, x)
    _tf = TransferFunction(k, x**4 - k, x)
    TF_ = TransferFunction(x**2, x**10 + x + x**2, x)
    _TF = TransferFunction(x**10 + x + x**2, x**2, x)
    assert G1.poles() == [I, I, -I, -I]
    assert G2.poles() == []
    assert tf1.poles() == [-5, 1]
    assert expect4_.poles() == [s]
    assert SP4.poles() == [-a0*s]
    assert expect3.poles() == [-0.25*p]
    assert str(expect2.poles()) == str([0.729001428685125, -0.564500714342563 - 0.710198984796332*I, -0.564500714342563 + 0.710198984796332*I])
    assert str(expect1.poles()) == str([-0.4 - 0.66332495807108*I, -0.4 + 0.66332495807108*I])
    assert _tf.poles() == [k**(Rational(1, 4)), -k**(Rational(1, 4)), I*k**(Rational(1, 4)), -I*k**(Rational(1, 4))]
    assert TF_.poles() == [CRootOf(x**9 + x + 1, 0), 0, CRootOf(x**9 + x + 1, 1), CRootOf(x**9 + x + 1, 2),
        CRootOf(x**9 + x + 1, 3), CRootOf(x**9 + x + 1, 4), CRootOf(x**9 + x + 1, 5), CRootOf(x**9 + x + 1, 6),
        CRootOf(x**9 + x + 1, 7), CRootOf(x**9 + x + 1, 8)]
    raises(NotImplementedError, lambda: TransferFunction(x**2, a0*x**10 + x + x**2, x).poles())

    # Stability of a transfer function.
    q, r = symbols('q, r', negative=True)
    t = symbols('t', positive=True)
    TF_ = TransferFunction(s**2 + a0 - a1*p, q*s - r, s)
    stable_tf = TransferFunction(s**2 + a0 - a1*p, q*s - 1, s)
    stable_tf_ = TransferFunction(s**2 + a0 - a1*p, q*s - t, s)

    assert G1.is_stable() is False
    assert G2.is_stable() is True
    assert tf1.is_stable() is False
    assert tf1.is_stable(True) is True
    assert expect2.is_stable() is False
    assert expect1.is_stable() is True
    assert stable_tf.is_stable() is True
    assert stable_tf_.is_stable() is True
    assert TF_.is_stable() is False
    assert expect4_.is_stable() is None   # no assumption provided for the only pole 's'.
    assert SP4.is_stable() is None

    generic_den = b4 * s**4 + b3 * s**3 + b2 * s**2 + b1 * s + b0

    stab_cond = TransferFunction(1, generic_den, s).get_asymptotic_stability_conditions()
    assert stab_cond == [
        b3*b4 > 0, -b1*b4 + b2*b3 > 0,
        -b0*b3**2*b4 -b1**2*b4**2 + b1*b2*b3*b4 > 0,
        b0*b4 > 0]
    assert TransferFunction(1, (s+1)*(s+2*I)*(s-2*I), s).get_asymptotic_stability_conditions() == [false]
    assert TransferFunction(1, (s+1)*(s+2)*(s+1/2), s).get_asymptotic_stability_conditions() == [true, true, true]
    assert stable_tf.get_asymptotic_stability_conditions() == [True]

    # Zeros of a transfer function.
    assert G1.zeros() == [1, 1]
    assert G2.zeros() == []
    assert tf1.zeros() == [-3, 1]
    assert expect4_.zeros() == [
        -7**(S(2)/3)*s**(S(1)/3)/7,
        7**(S(2)/3)*s**(S(1)/3)/14 - sqrt(3)*7**(S(2)/3)*I*s**(S(1)/3)/14,
        7**(S(2)/3)*s**(S(1)/3)/14 + sqrt(3)*7**(S(2)/3)*I*s**(S(1)/3)/14
    ]
    assert SP4.zeros() == [
        s**(S(1)/3)/a1**(S(1)/3),
        -s**(S(1)/3)/(2*a1**(S(1)/3)) - sqrt(3)*I*s**(S(1)/3)/(2*a1**(S(1)/3)),
        -s**(S(1)/3)/(2*a1**(S(1)/3)) + sqrt(3)*I*s**(S(1)/3)/(2*a1**(S(1)/3))
    ]
    assert str(expect3.zeros()) == str([0.125 - 1.11102430216445*sqrt(-0.405063291139241*p**3 - 1.0),
        1.11102430216445*sqrt(-0.405063291139241*p**3 - 1.0) + 0.125])
    assert tf_.zeros() == [k**(Rational(1, 3)), -k**(Rational(1, 3))/2 - sqrt(3)*I*k**(Rational(1, 3))/2,
        -k**(Rational(1, 3))/2 + sqrt(3)*I*k**(Rational(1, 3))/2]
    assert _TF.zeros() == [CRootOf(x**9 + x + 1, 0), 0, CRootOf(x**9 + x + 1, 1), CRootOf(x**9 + x + 1, 2),
        CRootOf(x**9 + x + 1, 3), CRootOf(x**9 + x + 1, 4), CRootOf(x**9 + x + 1, 5), CRootOf(x**9 + x + 1, 6),
        CRootOf(x**9 + x + 1, 7), CRootOf(x**9 + x + 1, 8)]
    raises(NotImplementedError, lambda: TransferFunction(a0*x**10 + x + x**2, x**2, x).zeros())

    # negation of TF.
    tf2 = TransferFunction(s + 3, s**2 - s**3 + 9, s)
    tf3 = TransferFunction(-3*p + 3, 1 - p, p)
    assert -tf2 == TransferFunction(-s - 3, s**2 - s**3 + 9, s)
    assert -tf3 == TransferFunction(3*p - 3, 1 - p, p)

    # taking power of a TF.
    tf4 = TransferFunction(p + 4, p - 3, p)
    tf5 = TransferFunction(s**2 + 1, 1 - s, s)
    expect2 = TransferFunction((s**2 + 1)**3, (1 - s)**3, s)
    expect1 = TransferFunction((p + 4)**2, (p - 3)**2, p)
    assert (tf4*tf4).doit() == tf4**2 == pow(tf4, 2) == expect1
    assert (tf5*tf5*tf5).doit() == tf5**3 == pow(tf5, 3) == expect2
    assert tf5**0 == pow(tf5, 0) == TransferFunction(1, 1, s)
    assert Series(tf4).doit()**-1 == tf4**-1 == pow(tf4, -1) == TransferFunction(p - 3, p + 4, p)
    assert (tf5*tf5).doit()**-1 == tf5**-2 == pow(tf5, -2) == TransferFunction((1 - s)**2, (s**2 + 1)**2, s)

    raises(ValueError, lambda: tf4**(s**2 + s - 1))
    raises(ValueError, lambda: tf5**s)
    raises(ValueError, lambda: tf4**tf5)

    # SymPy's own functions.
    tf = TransferFunction(s - 1, s**2 - 2*s + 1, s)
    tf6 = TransferFunction(s + p, p**2 - 5, s)
    assert factor(tf) == TransferFunction(s - 1, (s - 1)**2, s)
    assert tf.num.subs(s, 2) == tf.den.subs(s, 2) == 1
    # subs & xreplace
    assert tf.subs(s, 2) == TransferFunction(s - 1, s**2 - 2*s + 1, s)
    assert tf6.subs(p, 3) == TransferFunction(s + 3, 4, s)
    assert tf3.xreplace({p: s}) == TransferFunction(-3*s + 3, 1 - s, s)
    raises(TypeError, lambda: tf3.xreplace({p: exp(2)}))
    assert tf3.subs(p, exp(2)) == tf3

    tf7 = TransferFunction(a0*s**p + a1*p**s, a2*p - s, s)
    assert tf7.xreplace({s: k}) == TransferFunction(a0*k**p + a1*p**k, a2*p - k, k)
    assert tf7.subs(s, k) == TransferFunction(a0*s**p + a1*p**s, a2*p - s, s)

    # Conversion to Expr with to_expr()
    tf8 = TransferFunction(a0*s**5 + 5*s**2 + 3, s**6 - 3, s)
    tf9 = TransferFunction((5 + s), (5 + s)*(6 + s), s)
    tf10 = TransferFunction(0, 1, s)
    tf11 = TransferFunction(1, 1, s)
    assert tf8.to_expr() == Mul((a0*s**5 + 5*s**2 + 3), Pow((s**6 - 3), -1, evaluate=False), evaluate=False)
    assert tf9.to_expr() == Mul((s + 5), Pow((5 + s)*(6 + s), -1, evaluate=False), evaluate=False)
    assert tf10.to_expr() == Mul(S(0), Pow(1, -1, evaluate=False), evaluate=False)
    assert tf11.to_expr() == Pow(1, -1, evaluate=False)


def test_TransferFunction_addition_and_subtraction():
    tf1 = TransferFunction(s + 6, s - 5, s)
    tf2 = TransferFunction(s + 3, s + 1, s)
    tf3 = TransferFunction(s + 1, s**2 + s + 1, s)
    tf4 = TransferFunction(p, 2 - p, p)

    # addition
    assert tf1 + tf2 == Parallel(tf1, tf2)
    assert tf3 + tf1 == Parallel(tf3, tf1)
    assert -tf1 + tf2 + tf3 == Parallel(-tf1, tf2, tf3)
    assert tf1 + (tf2 + tf3) == Parallel(tf1, tf2, tf3)

    c = symbols("c", commutative=False)
    raises(ValueError, lambda: tf1 + Matrix([1, 2, 3]))
    raises(ValueError, lambda: tf2 + c)
    raises(ValueError, lambda: tf3 + tf4)
    raises(ValueError, lambda: tf1 + (s - 1))
    raises(ValueError, lambda: tf1 + 8)
    raises(ValueError, lambda: (1 - p**3) + tf1)

    # subtraction
    assert tf1 - tf2 == Parallel(tf1, -tf2)
    assert tf3 - tf2 == Parallel(tf3, -tf2)
    assert -tf1 - tf3 == Parallel(-tf1, -tf3)
    assert tf1 - tf2 + tf3 == Parallel(tf1, -tf2, tf3)

    raises(ValueError, lambda: tf1 - Matrix([1, 2, 3]))
    raises(ValueError, lambda: tf3 - tf4)
    raises(ValueError, lambda: tf1 - (s - 1))
    raises(ValueError, lambda: tf1 - 8)
    raises(ValueError, lambda: (s + 5) - tf2)
    raises(ValueError, lambda: (1 + p**4) - tf1)

    dtf1 = DiscreteTransferFunction(s, s+1, s, 0.01)
    dtf2 = DiscreteTransferFunction(s + 1, s + 2, s, 0.01)
    # addition and subtraction with discrete transfer functions raises TypeError
    raises(TypeError, lambda: dtf1 + tf1)
    raises(TypeError, lambda: dtf1 - tf1)
    raises(TypeError, lambda: (dtf2 + dtf1) - tf1)
    raises(TypeError, lambda: (dtf1 - dtf2) + tf1)

def test_TransferFunction_multiplication_and_division():
    G1 = TransferFunction(s + 3, -s**3 + 9, s)
    G2 = TransferFunction(s + 1, s - 5, s)
    G3 = TransferFunction(p, p**4 - 6, p)
    G4 = TransferFunction(p + 4, p - 5, p)
    G5 = TransferFunction(s + 6, s - 5, s)
    G6 = TransferFunction(s + 3, s + 1, s)
    G7 = TransferFunction(1, 1, s)

    # multiplication
    assert G1*G2 == Series(G1, G2)
    assert -G1*G5 == Series(-G1, G5)
    assert -G2*G5*-G6 == Series(-G2, G5, -G6)
    assert -G1*-G2*-G5*-G6 == Series(-G1, -G2, -G5, -G6)
    assert G3*G4 == Series(G3, G4)
    assert (G1*G2)*-(G5*G6) == \
        Series(G1, G2, TransferFunction(-1, 1, s), Series(G5, G6))
    assert G1*G2*(G5 + G6) == Series(G1, G2, Parallel(G5, G6))

    # division - See ``test_Feedback_functions()`` for division by Parallel objects.
    assert G5/G6 == Series(G5, pow(G6, -1))
    assert -G3/G4 == Series(-G3, pow(G4, -1))
    assert (G5*G6)/G7 == Series(G5, G6, pow(G7, -1))

    c = symbols("c", commutative=False)
    raises(ValueError, lambda: G3 * Matrix([1, 2, 3]))
    raises(ValueError, lambda: G1 * c)
    raises(ValueError, lambda: G3 * G5)
    raises(ValueError, lambda: G5 * (s - 1))
    raises(ValueError, lambda: 9 * G5)

    raises(ValueError, lambda: G3 / Matrix([1, 2, 3]))
    raises(ValueError, lambda: G6 / 0)
    raises(ValueError, lambda: G3 / G5)
    raises(ValueError, lambda: G5 / 2)
    raises(ValueError, lambda: G5 / s**2)
    raises(ValueError, lambda: (s - 4*s**2) / G2)
    raises(ValueError, lambda: 0 / G4)
    raises(ValueError, lambda: G7 / (1 + G6))
    raises(ValueError, lambda: G7 / (G5 * G6))
    raises(ValueError, lambda: G7 / (G7 + (G5 + G6)))

    dtf1 = DiscreteTransferFunction(s, s+1, s, 0.01)
    dtf2 = DiscreteTransferFunction(s + 1, s + 2, s, 0.01)
    # multiplication and division with discrete transfer functions raises TypeError
    raises(TypeError, lambda: dtf1 * G1)
    raises(TypeError, lambda: dtf1 / G1)
    raises(TypeError, lambda: (dtf1 * dtf2) / G2)
    raises(TypeError, lambda: (dtf1 / dtf2) * G2)


def test_TransferFunction_is_proper():
    omega_o, zeta, tau = symbols('omega_o, zeta, tau')
    G1 = TransferFunction(omega_o**2, s**2 + p*omega_o*zeta*s + omega_o**2, omega_o)
    G2 = TransferFunction(tau - s**3, tau + p**4, tau)
    G3 = TransferFunction(a*b*s**3 + s**2 - a*p + s, b - s*p**2, p)
    G4 = TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
    assert G1.is_proper
    assert G2.is_proper
    assert G3.is_proper
    assert not G4.is_proper


def test_TransferFunction_is_strictly_proper():
    omega_o, zeta, tau = symbols('omega_o, zeta, tau')
    tf1 = TransferFunction(omega_o**2, s**2 + p*omega_o*zeta*s + omega_o**2, omega_o)
    tf2 = TransferFunction(tau - s**3, tau + p**4, tau)
    tf3 = TransferFunction(a*b*s**3 + s**2 - a*p + s, b - s*p**2, p)
    tf4 = TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
    assert not tf1.is_strictly_proper
    assert not tf2.is_strictly_proper
    assert tf3.is_strictly_proper
    assert not tf4.is_strictly_proper


def test_TransferFunction_is_biproper():
    tau, omega_o, zeta = symbols('tau, omega_o, zeta')
    tf1 = TransferFunction(omega_o**2, s**2 + p*omega_o*zeta*s + omega_o**2, omega_o)
    tf2 = TransferFunction(tau - s**3, tau + p**4, tau)
    tf3 = TransferFunction(a*b*s**3 + s**2 - a*p + s, b - s*p**2, p)
    tf4 = TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
    assert tf1.is_biproper
    assert tf2.is_biproper
    assert not tf3.is_biproper
    assert not tf4.is_biproper


def test_DiscreteTransferFunction_construction():
    tf = DiscreteTransferFunction(z + 1, z**2 + z + 1, z, 10)
    assert tf.num == (z + 1)
    assert tf.den == (z**2 + z + 1)
    assert tf.args == (z + 1, z**2 + z + 1, z, 10)
    assert tf.sampling_time == 10

    tf1 = DiscreteTransferFunction(z + 4, z - 5, z, 0.1)
    assert tf1.num == (z + 4)
    assert tf1.den == (z - 5)
    assert tf1.sampling_time == 0.1
    assert isinstance(tf1.sampling_time, Number) is True # ensure it is a sympy object, not just a python float
    assert tf1.args == (z + 4, z - 5, z, 0.1)

    # using different polynomial variables.
    tf2 = DiscreteTransferFunction(p + 3, p**2 - 9, p, 12)
    assert tf2.num == (p + 3)
    assert tf2.den == (p**2 - 9)
    assert tf2.sampling_time == 12
    assert tf2.args == (p + 3, p**2 - 9, p, 12)

    tf3 = DiscreteTransferFunction(p**3 + 5*p**2 + 4, p**4 + 3*p + 1, p)
    assert tf3.args == (p**3 + 5*p**2 + 4, p**4 + 3*p + 1, p, 1)

    # no pole-zero cancellation on its own.
    tf4 = DiscreteTransferFunction((s + 3)*(s - 1), (s - 1)*(s + 5), s, 6.4)
    assert tf4.den == (s - 1)*(s + 5)
    assert tf4.args == ((s + 3)*(s - 1), (s - 1)*(s + 5), s, 6.4)

    tf4_ = DiscreteTransferFunction(p + 2, p + 2, p, 1.2)
    assert tf4_.args == (p + 2, p + 2, p, 1.2)

    tf5 = DiscreteTransferFunction(s - 1, 4 - p, s)
    assert tf5.args == (s - 1, 4 - p, s, 1)

    tf5_ = DiscreteTransferFunction(s - 1, s - 1, s)
    assert tf5_.args == (s - 1, s - 1, s, 1)

    tf6 = DiscreteTransferFunction(5, 6, s)
    assert tf6.num == 5
    assert tf6.den == 6
    assert tf6.args == (5, 6, s, 1)

    tf6_ = DiscreteTransferFunction(0.5, 4, s)
    assert tf6_.num == 0.5
    assert tf6_.den == 4
    assert tf6_.args == (0.500000000000000, 4, s, 1)

    tf7 = DiscreteTransferFunction(3*s**2 + 2*p + 4*s, 8*p**2 + 7*s, s)
    tf8 = DiscreteTransferFunction(3*s**2 + 2*p + 4*s, 8*p**2 + 7*s, p)
    assert not tf7 == tf8

    tf7_ = DiscreteTransferFunction(a0*s + a1*s**2 + a2*s**3, b0*p - b1*s, s)
    tf8_ = DiscreteTransferFunction(a0*s + a1*s**2 + a2*s**3, b0*p - b1*s, s)
    assert tf7_ == tf8_
    assert -(-tf7_) == tf7_ == -(-(-(-tf7_)))

    tf7_2 = DiscreteTransferFunction(a0*s + a1*s**2 + a2*s**3, b0*p - b1*s, s, 2)
    tf8_2 = DiscreteTransferFunction(a0*s + a1*s**2 + a2*s**3, b0*p - b1*s, s)
    assert not tf7_2 == tf8_2

    tf9 = DiscreteTransferFunction(a*s**3 + b*s**2 + g*s + d, d*p + g*p**2 + g*s, s)
    assert tf9.args == (a*s**3 + b*s**2 + d + g*s, d*p + g*p**2 + g*s, s, 1)

    tf10 = DiscreteTransferFunction(p**3 + d, g*s**2 + d*s + a, p, 2)
    tf10_ = DiscreteTransferFunction(p**3 + d, g*s**2 + d*s + a, p, 2)
    assert tf10.args == (d + p**3, a + d*s + g*s**2, p, 2)
    assert tf10_ == tf10

    tf11 = DiscreteTransferFunction(a1*s + a0, b2*s**2 + b1*s + b0, s)
    assert tf11.num == (a0 + a1*s)
    assert tf11.den == (b0 + b1*s + b2*s**2)
    assert tf11.args == (a0 + a1*s, b0 + b1*s + b2*s**2, s, 1)

    # when just the numerator is 0, leave the denominator alone.
    tf12 = DiscreteTransferFunction(0, p**2 - p + 1, p)
    assert tf12.args == (0, p**2 - p + 1, p, 1)

    tf13 = DiscreteTransferFunction(0, 1, s)
    assert tf13.args == (0, 1, s, 1)

    # float exponents
    tf14 = DiscreteTransferFunction(a0*s**0.5 + a2*s**0.6 - a1, a1*p**(-8.7), s, 6)
    assert tf14.args == (a0*s**0.5 - a1 + a2*s**0.6, a1*p**(-8.7), s, 6)

    tf15 = DiscreteTransferFunction(a2**2*p**(0.25) + a1*s**(-0.8), a0*s - p, p)
    assert tf15.args == (a1*s**(-0.8) + a2**2*p**0.25, a0*s - p, p, 1)

    omega_o, k_p, k_o, k_i = symbols('omega_o, k_p, k_o, k_i')
    tf18 = DiscreteTransferFunction((k_p + k_o*s + k_i/s),
                              s**2 + 2*omega_o*s + omega_o**2, s)
    assert tf18.num == k_i/s + k_o*s + k_p
    assert tf18.args == (k_i/s + k_o*s + k_p,
                         omega_o**2 + 2*omega_o*s + s**2, s, 1)

    # ValueError when denominator is zero.
    raises(ValueError, lambda: DiscreteTransferFunction(4, 0, s))
    raises(ValueError, lambda: DiscreteTransferFunction(s, 0, s))
    raises(ValueError, lambda: DiscreteTransferFunction(0, 0, s))

    raises(TypeError, lambda: DiscreteTransferFunction(Matrix([1, 2, 3]), s, s))

    raises(TypeError, lambda: DiscreteTransferFunction(s**2 + 2*s - 1, s + 3, 3))
    raises(TypeError, lambda: DiscreteTransferFunction(p + 1, 5 - p, 4))
    raises(TypeError, lambda: DiscreteTransferFunction(3, 4, 8))

    raises(ValueError, lambda: DiscreteTransferFunction(s + 1, s**2 + 2, s, 0)) # sampling time cannot be zero

def _are_floats_equals(n1, n2):
    """
    Helper function to compare two floats.
    Returns True if they are approximately equal, otherwise False.

    """
    return Float(n1).is_same(Float(n2), isclose)

def _are_complex_equals(n1, n2):
    """
    Helper function to compare two complex numbers.
    Returns True if they are approximately equal, otherwise False.

    """
    return (_are_floats_equals(re(n1), re(n2)) and
            _are_floats_equals(im(n1), im(n2)))

def test_DiscreteTransferFunction_functions():
    # classmethod from_rational_expression
    expr_1 = Mul(0, Pow(s, -1, evaluate=False), evaluate=False)
    expr_2 = s/0
    expr_3 = (p*s**2 + 5*s)/(s + 1)**3
    expr_4 = 6
    expr_5 = ((2 + 3*s)*(5 + 2*s))/((9 + 3*s)*(5 + 2*s**2))
    expr_6 = (9*s**4 + 4*s**2 + 8)/((s + 1)*(s + 9))
    tf = DiscreteTransferFunction(s + 1, s**2 + 2, s)
    delay = s**(-1/tau)
    expr_7 = delay*tf.to_expr()
    H1 = DiscreteTransferFunction.from_rational_expression(expr_7, s)
    H2 = DiscreteTransferFunction(s + 1, s**(1/tau)*(s**2 + 2), s)
    expr_8 = Add(2,  3*s/(s**2 + 1), evaluate=False)

    assert DiscreteTransferFunction.from_rational_expression(expr_1,
                                                       sampling_time = 12) == \
        DiscreteTransferFunction(0, s, s, 12)
    raises(ZeroDivisionError, lambda:
           DiscreteTransferFunction.from_rational_expression(expr_2))
    raises(ValueError, lambda:
           DiscreteTransferFunction.from_rational_expression(expr_3))
    assert DiscreteTransferFunction.from_rational_expression(expr_3, s) == \
        DiscreteTransferFunction((p*s**2 + 5*s), (s + 1)**3, s)
    assert DiscreteTransferFunction.from_rational_expression(expr_3, p,
                                                       sampling_time = 1.4) == \
        DiscreteTransferFunction((p*s**2 + 5*s), (s + 1)**3, p, 1.4)
    raises(ValueError, lambda:
           DiscreteTransferFunction.from_rational_expression(expr_4))
    assert DiscreteTransferFunction.from_rational_expression(expr_4, s) == \
        DiscreteTransferFunction(6, 1, s)
    assert DiscreteTransferFunction.from_rational_expression(expr_5, s) == \
        DiscreteTransferFunction((2 + 3*s)*(5 + 2*s), (9 + 3*s)*(5 + 2*s**2), s)
    assert DiscreteTransferFunction.from_rational_expression(expr_6, s) == \
        DiscreteTransferFunction((9*s**4 + 4*s**2 + 8), (s + 1)*(s + 9), s)
    assert H1 == H2
    assert DiscreteTransferFunction.from_rational_expression(expr_8, s,
                                                       sampling_time = 0.3) == \
        DiscreteTransferFunction(2*s**2 + 3*s + 2, s**2 + 1, s, 0.3)

    # classmethod from_coeff_lists
    tf1 = DiscreteTransferFunction.from_coeff_lists([1, 2], [3, 4, 5], s)
    num2 = [p**2, 2*p]
    den2 = [p**3, p + 1, 4]
    tf2 = DiscreteTransferFunction.from_coeff_lists(num2, den2, s)
    num3 = [1, 2, 3]
    den3 = [0, 0]

    assert tf1 == DiscreteTransferFunction(s + 2, 3*s**2 + 4*s + 5, s)
    assert tf2 == DiscreteTransferFunction(p**2*s + 2*p, p**3*s**2 + s*(p + 1) + 4, s)
    raises(ZeroDivisionError, lambda:
           DiscreteTransferFunction.from_coeff_lists(num3, den3, s))

    # classmethod from_zpk
    zeros = [4]
    poles = [-1+2j, -1-2j]
    gain = 3
    tf1 = DiscreteTransferFunction.from_zpk(zeros, poles, gain, s, 3)

    assert tf1 == DiscreteTransferFunction(3*s - 12,
                                     (s + 1.0 - 2.0*I)*(s + 1.0 + 2.0*I), s, 3)

    # explicitly cancel poles and zeros.
    tf0 = DiscreteTransferFunction(s**5 + s**3 + s, s - s**2, s)
    a = DiscreteTransferFunction(-(s**4 + s**2 + 1), s - 1, s)
    assert tf0.simplify() == simplify(tf0) == a

    tf1 = DiscreteTransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
    b = DiscreteTransferFunction(p + 3, p + 5, p)
    assert tf1.simplify() == simplify(tf1) == b

    # expand the numerator and the denominator.
    G1 = DiscreteTransferFunction((1 - s)**2, (s**2 + 1)**2, s)
    G2 = DiscreteTransferFunction(1, -3, p)
    c = (a2*s**p + a1*s**s + a0*p**p)*(p**s + s**p)
    d = (b0*s**s + b1*p**s)*(b2*s*p + p**p)
    e = a0*p**p*p**s + a0*p**p*s**p + a1*p**s*s**s + a1*s**p*s**s + \
        a2*p**s*s**p + a2*s**(2*p)
    f = b0*b2*p*s*s**s + b0*p**p*s**s + b1*b2*p*p**s*s + b1*p**p*p**s
    g = a1*a2*s*s**p + a1*p*s + a2*b1*p*s*s**p + b1*p**2*s
    G3 = DiscreteTransferFunction(c, d, s)
    G4 = DiscreteTransferFunction(a0*s**s - b0*p**p, (a1*s + b1*s*p)*(a2*s**p + p), p,
                             2.4)

    assert G1.expand() == DiscreteTransferFunction(s**2 - 2*s + 1, s**4 + 2*s**2 + 1,
                                             s)
    assert tf1.expand() == DiscreteTransferFunction(p**2 + 2*p - 3, p**2 + 4*p - 5, p)
    assert G2.expand() == G2
    assert G3.expand() == DiscreteTransferFunction(e, f, s)
    assert G4.expand() == DiscreteTransferFunction(a0*s**s - b0*p**p, g, p, 2.4)
    assert G4.expand() != DiscreteTransferFunction(a0*s**s - b0*p**p, g, p)

    # testing that subs works.
    p1 = a1*s + a0
    p2 = b2*s**2 + b1*s + b0
    SP1 = DiscreteTransferFunction(p1, p2, s)
    expect1 = DiscreteTransferFunction(2.0*s + 1.0, 5.0*s**2 + 4.0*s + 3.0, s)
    expect1_ = DiscreteTransferFunction(2*s + 1, 5*s**2 + 4*s + 3, s)
    assert SP1.subs({a0: 1, a1: 2, b0: 3, b1: 4, b2: 5}) == expect1_
    assert SP1.subs({a0: 1, a1: 2, b0: 3, b1: 4, b2: 5}).evalf() == expect1
    assert expect1_.evalf() == expect1

    c1, d0, d1, d2 = symbols('c1, d0:3')
    p3, p4 = c1*p, d2*p**3 + d1*p**2 - d0
    SP2 = DiscreteTransferFunction(p3, p4, p)
    expect2 = DiscreteTransferFunction(2.0*p, 5.0*p**3 + 2.0*p**2 - 3.0, p)
    expect2_ = DiscreteTransferFunction(2*p, 5*p**3 + 2*p**2 - 3, p)
    assert SP2.subs({c1: 2, d0: 3, d1: 2, d2: 5}) == expect2_
    assert SP2.subs({c1: 2, d0: 3, d1: 2, d2: 5}).evalf() == expect2
    assert expect2_.evalf() == expect2

    SP3 = DiscreteTransferFunction(a0*p**3 + a1*s**2 - b0*s + b1, a1*s + p, s)
    expect3 = DiscreteTransferFunction(2.0*p**3 + 4.0*s**2 - s + 5.0, p + 4.0*s, s)
    expect3_ = DiscreteTransferFunction(2*p**3 + 4*s**2 - s + 5, p + 4*s, s)
    assert SP3.subs({a0: 2, a1: 4, b0: 1, b1: 5}) == expect3_
    assert SP3.subs({a0: 2, a1: 4, b0: 1, b1: 5}).evalf() == expect3
    assert expect3_.evalf() == expect3

    SP4 = DiscreteTransferFunction(s - a1*p**3, a0*s + p, p)
    expect4 = DiscreteTransferFunction(7.0*p**3 + s, p - s, p)
    expect4_ = DiscreteTransferFunction(7*p**3 + s, p - s, p)
    assert SP4.subs({a0: -1, a1: -7}) == expect4_
    assert SP4.subs({a0: -1, a1: -7}).evalf() == expect4
    assert expect4_.evalf() == expect4

    # evaluate the transfer function at particular frequencies.
    assert tf1.eval_frequency(wn) == wn/(wn + 5) + 3/(wn + 5)
    dtf_fr = DiscreteTransferFunction(z - 1, z**2 + z + 1, z)
    assert dtf_fr.eval_frequency(1 + I) == Rational(3, 13) + Rational(2, 13) * I
    assert G4.eval_frequency(S(5)/3) == (
        a0*s**s/(a1*a2*s**(S(8)/3) + S(5)*a1*s/3 +
                 5*a2*b1*s**(S(8)/3)/3 + S(25)*b1*s/9) -
                    5*3**(S(1)/3)*5**(S(2)/3)*
                        b0/(9*a1*a2*s**(S(8)/3) + 15*a1*s +
                            15*a2*b1*s**(S(8)/3) + 25*b1*s))

    # Low-frequency (or DC) gain.
    assert DiscreteTransferFunction(z, z - 1, z).dc_gain() == oo
    assert DiscreteTransferFunction(z - 1, z**2 + 3*z + 1, z).dc_gain() == 0
    assert DiscreteTransferFunction(z, z + 2, z, 0.2).dc_gain() == Rational(1, 3)

    # Poles of a transfer function.
    tf_ = DiscreteTransferFunction(x**3 - k, k, x, 3)
    _tf = DiscreteTransferFunction(k, x**4 - k, x, 1.6)
    TF_ = DiscreteTransferFunction(x**2, x**10 + x + x**2, x)
    _TF = DiscreteTransferFunction(x**10 + x + x**2, x**2, x)
    assert G1.poles() == [I, I, -I, -I]
    assert G2.poles() == []
    assert tf1.poles() == [-5, 1]
    assert expect4_.poles() == [s]
    assert SP4.poles() == [-a0*s]
    assert expect3.poles() == [-0.25*p]

    exp_poles1 = [-0.4 - 0.66332495807108*I, -0.4 + 0.66332495807108*I]
    for i, pol in enumerate(expect1.poles()):
        assert _are_complex_equals(pol, exp_poles1[i])

    exp_poles2 = [0.729001428685125, -0.564500714342563 - 0.710198984796332*I,
             -0.564500714342563 + 0.710198984796332*I]
    for i, pol in enumerate(expect2.poles()):
        assert _are_complex_equals(pol, exp_poles2[i])

    assert _tf.poles() == [k**(Rational(1, 4)), -k**(Rational(1, 4)),
                           I*k**(Rational(1, 4)), -I*k**(Rational(1, 4))]
    assert TF_.poles() == [CRootOf(x**9 + x + 1, 0), 0,
                           CRootOf(x**9 + x + 1, 1), CRootOf(x**9 + x + 1, 2),
                           CRootOf(x**9 + x + 1, 3), CRootOf(x**9 + x + 1, 4),
                           CRootOf(x**9 + x + 1, 5), CRootOf(x**9 + x + 1, 6),
                           CRootOf(x**9 + x + 1, 7), CRootOf(x**9 + x + 1, 8)]
    # test generic poles denominator
    assert DiscreteTransferFunction(1, b0*x**2 + b1*x + b2, x).poles() == [
        -b1/(2*b0) - sqrt(-4*b0*b2 + b1**2)/(2*b0),
        -b1/(2*b0) + sqrt(-4*b0*b2 + b1**2)/(2*b0)]
    raises(NotImplementedError, lambda:
           DiscreteTransferFunction(x**2, a0*x**10 + x + x**2, x).poles())

    # Stability of a transfer function.
    stable_tf = DiscreteTransferFunction(
        z,
        (z - Rational(1,2))*(z + Rational(9,10)), z, 0.1)
    marginally_stable_tf = DiscreteTransferFunction(
        z, (z - Rational(1,4))*(z + 1), z)
    unstable_tf = DiscreteTransferFunction(
        z, (z - Rational(1,4))*(z + 2), z, T)

    assert stable_tf.is_stable() == True
    assert marginally_stable_tf.is_stable() == False
    assert unstable_tf.is_stable() == False
    assert DiscreteTransferFunction(
        z, (z+Rational(1,8))*(z+k), z).is_stable() == None

    generic_den = b4 * s**4 + b3 * s**3 + b2 * s**2 + b1 * s + b0

    stab_cond = DiscreteTransferFunction(1, generic_den, s).\
                get_asymptotic_stability_conditions()
    assert stab_cond == [
        (-4*b0**2 + 6*b0*b1 - 4*b0*b2 + 2*b0*b3 - 2*b1**2 + 2*b1*b2 - 2*b1*b4 -
         2*b2*b3 + 4*b2*b4 + 2*b3**2 - 6*b3*b4 + 4*b4**2) > 0,
        (-20*b0**2 + 10*b0*b1 + 12*b0*b2 - 18*b0*b3 - 2*b1**2 - 2*b1*b2 +
         18*b1*b4 + 2*b2*b3 - 12*b2*b4 + 2*b3**2 - 10*b3*b4 + 20*b4**2) > 0,
        (64*b0**4 - 64*b0**3*b1 - 64*b0**3*b3 + 64*b0**2*b1*b2 +
         64*b0**2*b1*b3 + 64*b0**2*b1*b4 - 64*b0**2*b2**2 + 64*b0**2*b2*b3 -
         64*b0**2*b3**2 + 64*b0**2*b3*b4 - 128*b0**2*b4**2 - 64*b0*b1**2*b3 -
         64*b0*b1**2*b4 + 64*b0*b1*b2*b3 - 128*b0*b1*b2*b4 + 128*b0*b1*b3*b4 +
         64*b0*b1*b4**2 + 128*b0*b2**2*b4 - 64*b0*b2*b3**2 - 128*b0*b2*b3*b4 +
         64*b0*b3**3 - 64*b0*b3**2*b4 + 64*b0*b3*b4**2 + 64*b1**3*b4 -
         64*b1**2*b2*b4 - 64*b1**2*b4**2 + 64*b1*b2*b3*b4 + 64*b1*b2*b4**2 -
         64*b1*b3**2*b4 + 64*b1*b3*b4**2 - 64*b1*b4**3 - 64*b2**2*b4**2 +
         64*b2*b3*b4**2 - 64*b3*b4**3 + 64*b4**4) > 0,
        (b0**2 + 2*b0*b2 + 2*b0*b4 - b1**2 - 2*b1*b3 + b2**2 + 2*b2*b4 -
         b3**2 + b4**2) > 0]

    # Zeros of a transfer function.
    assert G1.zeros() == [1, 1]
    assert G2.zeros() == []
    assert tf1.zeros() == [-3, 1]
    assert expect4_.zeros() == [
        -7**(S(2)/3)*s**(S(1)/3)/7,
        7**(S(2)/3)*s**(S(1)/3)/14 - sqrt(3)*7**(S(2)/3)*I*s**(S(1)/3)/14,
        7**(S(2)/3)*s**(S(1)/3)/14 + sqrt(3)*7**(S(2)/3)*I*s**(S(1)/3)/14
    ]
    assert SP4.zeros() == [
        s**(S(1)/3)/a1**(S(1)/3),
        -s**(S(1)/3)/(2*a1**(S(1)/3)) - sqrt(3)*I*s**(S(1)/3)/(2*a1**(S(1)/3)),
        -s**(S(1)/3)/(2*a1**(S(1)/3)) + sqrt(3)*I*s**(S(1)/3)/(2*a1**(S(1)/3))
    ]

    assert str(expect3.zeros()) == \
        str([0.125 - 1.11102430216445*sqrt(-0.405063291139241*p**3 - 1.0),
            1.11102430216445*sqrt(-0.405063291139241*p**3 - 1.0) + 0.125])

    assert tf_.zeros() == [
        k**(Rational(1, 3)),
        -k**(Rational(1, 3))/2 - sqrt(3)*I*k**(Rational(1, 3))/2,
        -k**(Rational(1, 3))/2 + sqrt(3)*I*k**(Rational(1, 3))/2]
    assert _TF.zeros() == [
        CRootOf(x**9 + x + 1, 0), 0,
        CRootOf(x**9 + x + 1, 1), CRootOf(x**9 + x + 1, 2),
        CRootOf(x**9 + x + 1, 3), CRootOf(x**9 + x + 1, 4),
        CRootOf(x**9 + x + 1, 5), CRootOf(x**9 + x + 1, 6),
        CRootOf(x**9 + x + 1, 7), CRootOf(x**9 + x + 1, 8)]
    raises(NotImplementedError, lambda:
           DiscreteTransferFunction(a0*x**10 + x + x**2, x**2, x).zeros())

    # negation of TF.
    tf2 = DiscreteTransferFunction(s + 3, s**2 - s**3 + 9, s)
    tf3 = DiscreteTransferFunction(-3*p + 3, 1 - p, p)
    assert -tf2 == DiscreteTransferFunction(-s - 3, s**2 - s**3 + 9, s)
    assert -tf3 == DiscreteTransferFunction(3*p - 3, 1 - p, p)

    # taking power of a TF.
    tf4 = DiscreteTransferFunction(p + 4, p - 3, p)
    tf5 = DiscreteTransferFunction(s**2 + 1, 1 - s, s)
    expect2 = DiscreteTransferFunction((s**2 + 1)**3, (1 - s)**3, s)
    expect1 = DiscreteTransferFunction((p + 4)**2, (p - 3)**2, p)
    assert (tf4*tf4).doit() == tf4**2 == pow(tf4, 2) == expect1
    assert (tf5*tf5*tf5).doit() == tf5**3 == pow(tf5, 3) == expect2
    assert tf5**0 == pow(tf5, 0) == DiscreteTransferFunction(1, 1, s)
    assert Series(tf4).doit()**-1 == tf4**-1 == pow(tf4, -1) == \
        DiscreteTransferFunction(p - 3, p + 4, p)
    assert (tf5*tf5).doit()**-1 == tf5**-2 == pow(tf5, -2) == \
        DiscreteTransferFunction((1 - s)**2, (s**2 + 1)**2, s)

    raises(ValueError, lambda: tf4**(s**2 + s - 1))
    raises(ValueError, lambda: tf5**s)
    raises(ValueError, lambda: tf4**tf5)

    # SymPy's own functions.
    tf = DiscreteTransferFunction(s - 1, s**2 - 2*s + 1, s)
    tf6 = DiscreteTransferFunction(s + p, p**2 - 5, s)
    assert factor(tf) == DiscreteTransferFunction(s - 1, (s - 1)**2, s)
    assert tf.num.subs(s, 2) == tf.den.subs(s, 2) == 1
    # subs & xreplace
    assert tf.subs(s, 2) == DiscreteTransferFunction(s - 1, s**2 - 2*s + 1, s)
    assert tf6.subs(p, 3) == DiscreteTransferFunction(s + 3, 4, s)
    assert tf3.xreplace({p: s}) == DiscreteTransferFunction(-3*s + 3, 1 - s, s)
    raises(TypeError, lambda: tf3.xreplace({p: exp(2)}))
    assert tf3.subs(p, exp(2)) == tf3

    tf7 = DiscreteTransferFunction(a0*s**p + a1*p**s, a2*p - s, s)
    assert tf7.xreplace({s: k}) == DiscreteTransferFunction(a0*k**p + a1*p**k,
                                                      a2*p - k, k)
    assert tf7.subs(s, k) == DiscreteTransferFunction(a0*s**p + a1*p**s, a2*p - s, s)

    # Conversion to Expr with to_expr()
    tf8 = DiscreteTransferFunction(a0*s**5 + 5*s**2 + 3, s**6 - 3, s)
    tf9 = DiscreteTransferFunction((5 + s), (5 + s)*(6 + s), s)
    tf10 = DiscreteTransferFunction(0, 1, s)
    tf11 = DiscreteTransferFunction(1, 1, s)
    assert tf8.to_expr() == Mul((a0*s**5 + 5*s**2 + 3),
                                Pow((s**6 - 3), -1, evaluate=False),
                                evaluate=False)
    assert tf9.to_expr() == Mul((s + 5),
                                Pow((5 + s)*(6 + s), -1, evaluate=False),
                                evaluate=False)
    assert tf10.to_expr() == Mul(S(0), Pow(1, -1, evaluate=False),
                                 evaluate=False)
    assert tf11.to_expr() == Pow(1, -1, evaluate=False)


def test_DiscreteTransferFunction_addition_and_subtraction():
    tf1 = DiscreteTransferFunction(s + 6, s - 5, s)
    tf2 = DiscreteTransferFunction(s + 3, s + 1, s)
    tf3 = DiscreteTransferFunction(s + 1, s**2 + s + 1, s)
    tf4 = DiscreteTransferFunction(p, 2 - p, p)
    tf5 = DiscreteTransferFunction(p, p + 5, p, 20)
    tf6 = DiscreteTransferFunction(1, p - 1, p, 20)

    # addition
    assert tf1 + tf2 == Parallel(tf1, tf2)
    assert tf3 + tf1 == Parallel(tf3, tf1)
    assert -tf1 + tf2 + tf3 == Parallel(-tf1, tf2, tf3)
    assert tf1 + (tf2 + tf3) == Parallel(tf1, tf2, tf3)
    assert tf5 + tf6 == Parallel(tf5, tf6)

    c = symbols("c", commutative=False)
    raises(ValueError, lambda: tf1 + Matrix([1, 2, 3]))
    raises(ValueError, lambda: tf2 + c)
    raises(ValueError, lambda: tf3 + tf4)
    raises(ValueError, lambda: tf1 + (s - 1))
    raises(ValueError, lambda: tf1 + 8)
    raises(ValueError, lambda: (1 - p**3) + tf1)

    # subtraction
    assert tf1 - tf2 == Parallel(tf1, -tf2)
    assert tf3 - tf2 == Parallel(tf3, -tf2)
    assert -tf1 - tf3 == Parallel(-tf1, -tf3)
    assert tf1 - tf2 + tf3 == Parallel(tf1, -tf2, tf3)
    assert tf5 - tf6 == Parallel(tf5, -tf6)

    raises(ValueError, lambda: tf1 - Matrix([1, 2, 3]))
    raises(ValueError, lambda: tf3 - tf4)
    raises(ValueError, lambda: tf1 - (s - 1))
    raises(ValueError, lambda: tf1 - 8)
    raises(ValueError, lambda: (s + 5) - tf2)
    raises(ValueError, lambda: (1 + p**4) - tf1)
    # can't add systems with different sampling time
    raises(TypeError, lambda: tf5 + tf1)
    raises(TypeError, lambda: tf5 - tf1)
    raises(TypeError, lambda: (tf5 + tf6) + tf1)
    raises(TypeError, lambda: (tf5 - tf6) - tf1)

def test_DiscreteTransferFunction_multiplication_and_division():
    G1 = DiscreteTransferFunction(s + 3, -s**3 + 9, s)
    G2 = DiscreteTransferFunction(s + 1, s - 5, s)
    G3 = DiscreteTransferFunction(p, p**4 - 6, p)
    G4 = DiscreteTransferFunction(p + 4, p - 5, p)
    G5 = DiscreteTransferFunction(s + 6, s - 5, s)
    G6 = DiscreteTransferFunction(s + 3, s + 1, s)
    G7 = DiscreteTransferFunction(1, 1, s)
    G8 = DiscreteTransferFunction(p, p + 1, p, 12)
    G9 = DiscreteTransferFunction(1, p, p, 12)

    # multiplication
    assert G1*G2 == Series(G1, G2)
    assert -G1*G5 == Series(-G1, G5)
    assert -G2*G5*-G6 == Series(-G2, G5, -G6)
    assert -G1*-G2*-G5*-G6 == Series(-G1, -G2, -G5, -G6)
    assert G3*G4 == Series(G3, G4)
    assert (G1*G2)*-(G5*G6) == \
        Series(G1, G2, DiscreteTransferFunction(-1, 1, s), Series(G5, G6))
    assert G1*G2*(G5 + G6) == Series(G1, G2, Parallel(G5, G6))
    assert G8*G9 == Series(G8, G9)

    # division - See ``test_Feedback_functions()`` for division by Parallel objects.
    assert G5/G6 == Series(G5, pow(G6, -1))
    assert -G3/G4 == Series(-G3, pow(G4, -1))
    assert (G5*G6)/G7 == Series(G5, G6, pow(G7, -1))
    assert G8/G9 == Series(G8, pow(G9, -1))

    c = symbols("c", commutative=False)
    raises(ValueError, lambda: G3 * Matrix([1, 2, 3]))
    raises(ValueError, lambda: G1 * c)
    raises(ValueError, lambda: G3 * G5)
    raises(ValueError, lambda: G5 * (s - 1))
    raises(ValueError, lambda: 9 * G5)

    raises(ValueError, lambda: G3 / Matrix([1, 2, 3]))
    raises(ValueError, lambda: G6 / 0)
    raises(ValueError, lambda: G3 / G5)
    raises(ValueError, lambda: G5 / 2)
    raises(ValueError, lambda: G5 / s**2)
    raises(ValueError, lambda: (s - 4*s**2) / G2)
    raises(ValueError, lambda: 0 / G4)
    raises(ValueError, lambda: G7 / (1 + G6))
    raises(ValueError, lambda: G7 / (G5 * G6))
    raises(ValueError, lambda: G7 / (G7 + (G5 + G6)))

    # can't add systems with different sampling time
    raises(TypeError, lambda: G8 * G2)
    raises(TypeError, lambda: (G8 * G9) * G2)
    raises(TypeError, lambda: (G8 / G9) * G2)


def test_DiscreteTransferFunction_is_proper():
    omega_o, zeta, tau = symbols('omega_o, zeta, tau')
    G1 = DiscreteTransferFunction(omega_o**2, s**2 + p*omega_o*zeta*s + omega_o**2,
                          omega_o)
    G2 = DiscreteTransferFunction(tau - s**3, tau + p**4, tau)
    G3 = DiscreteTransferFunction(a*b*s**3 + s**2 - a*p + s, b - s*p**2, p)
    G4 = DiscreteTransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
    assert G1.is_proper
    assert G2.is_proper
    assert G3.is_proper
    assert not G4.is_proper


def test_DiscreteTransferFunction_is_strictly_proper():
    omega_o, zeta, tau = symbols('omega_o, zeta, tau')
    tf1 = DiscreteTransferFunction(omega_o**2, s**2 + p*omega_o*zeta*s + omega_o**2,
                             omega_o)
    tf2 = DiscreteTransferFunction(tau - s**3, tau + p**4, tau)
    tf3 = DiscreteTransferFunction(a*b*s**3 + s**2 - a*p + s, b - s*p**2, p)
    tf4 = DiscreteTransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
    assert not tf1.is_strictly_proper
    assert not tf2.is_strictly_proper
    assert tf3.is_strictly_proper
    assert not tf4.is_strictly_proper


def test_DiscreteTransferFunction_is_biproper():
    tau, omega_o, zeta = symbols('tau, omega_o, zeta')
    tf1 = DiscreteTransferFunction(omega_o**2, s**2 + p*omega_o*zeta*s + omega_o**2, omega_o)
    tf2 = DiscreteTransferFunction(tau - s**3, tau + p**4, tau)
    tf3 = DiscreteTransferFunction(a*b*s**3 + s**2 - a*p + s, b - s*p**2, p)
    tf4 = DiscreteTransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
    assert tf1.is_biproper
    assert tf2.is_biproper
    assert not tf3.is_biproper
    assert not tf4.is_biproper


def test_PIDController():
    kp, ki, kd, tf = symbols("kp ki kd tf")
    p1 = PIDController(kp, ki, kd, tf)
    p2 = PIDController()

    # Type Checking
    assert isinstance(p1, PIDController)
    assert isinstance(p1, TransferFunction)

    # Properties checking
    assert p1 == PIDController(kp, ki, kd, tf, s)
    assert p2 == PIDController(kp, ki, kd, 0, s)
    assert p1.num == kd*s**2 + ki*s*tf + ki + kp*s**2*tf + kp*s
    assert p1.den == s**2*tf + s
    assert p1.var == s
    assert p1.kp == kp
    assert p1.ki == ki
    assert p1.kd == kd
    assert p1.tf == tf

    # Functionality checking
    assert p1.doit() == TransferFunction(kd*s**2 + ki*s*tf + ki + kp*s**2*tf + kp*s, s**2*tf + s, s)
    assert p1.is_proper == True
    assert p1.is_biproper == True
    assert p1.is_strictly_proper == False
    assert p2.doit() == TransferFunction(kd*s**2 + ki + kp*s, s, s)

    # Using PIDController with TransferFunction
    tf1 = TransferFunction(s, s + 1, s)
    par1 = Parallel(p1, tf1)
    ser1 = Series(p1, tf1)
    fed1 = Feedback(p1, tf1)
    assert par1 == Parallel(PIDController(kp, ki, kd, tf, s), TransferFunction(s, s + 1, s))
    assert ser1 == Series(PIDController(kp, ki, kd, tf, s), TransferFunction(s, s + 1, s))
    assert fed1 == Feedback(PIDController(kp, ki, kd, tf, s), TransferFunction(s, s + 1, s))
    assert par1.doit() == TransferFunction(s*(s**2*tf + s) + (s + 1)*(kd*s**2 + ki*s*tf + ki + kp*s**2*tf + kp*s),
                                           (s + 1)*(s**2*tf + s), s)
    assert ser1.doit() == TransferFunction(s*(kd*s**2 + ki*s*tf + ki + kp*s**2*tf + kp*s),
                                           (s + 1)*(s**2*tf + s), s)
    assert fed1.doit() == TransferFunction((s + 1)*(s**2*tf + s)*(kd*s**2 + ki*s*tf + ki + kp*s**2*tf + kp*s),
                                           (s*(kd*s**2 + ki*s*tf + ki + kp*s**2*tf + kp*s) + (s + 1)*(s**2*tf + s))*(s**2*tf + s), s)


def test_Series_construction():
    tf = TransferFunction(a0*s**3 + a1*s**2 - a2*s, b0*p**4 + b1*p**3 - b2*s*p, s)
    tf2 = TransferFunction(a2*p - s, a2*s + p, s)
    tf3 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf4 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    inp = Function('X_d')(s)
    out = Function('X')(s)

    dtf = DiscreteTransferFunction(a0*s**3 + a1*s**2 - a2*s,
                             b0*p**4 + b1*p**3 - b2*s*p, s, 0.1)
    dtf2 = DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.1)
    dtf3 = DiscreteTransferFunction(a0*p + p**a1 - s, p, p, 0.1)
    dtf4 = DiscreteTransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s, 0.1)

    assert tf.is_continuous is True
    assert dtf.is_continuous is False

    # continuous-time tests
    s0 = Series(tf, tf2)
    assert s0.args == (tf, tf2)
    assert s0.var == s

    s1 = Series(Parallel(tf, -tf2), tf2)
    assert s1.args == (Parallel(tf, -tf2), tf2)
    assert s1.var == s

    tf3_ = TransferFunction(inp, 1, s)
    tf4_ = TransferFunction(-out, 1, s)
    s2 = Series(tf, Parallel(tf3_, tf4_), tf2)
    assert s2.args == (tf, Parallel(tf3_, tf4_), tf2)

    s3 = Series(tf, tf2, tf4)
    assert s3.args == (tf, tf2, tf4)

    s4 = Series(tf3_, tf4_)
    assert s4.args == (tf3_, tf4_)
    assert s4.var == s

    s6 = Series(tf2, tf4, Parallel(tf2, -tf), tf4)
    assert s6.args == (tf2, tf4, Parallel(tf2, -tf), tf4)

    s7 = Series(tf, tf2)
    assert s0 == s7
    assert not s0 == s2

    raises(ValueError, lambda: Series(tf, tf3))
    raises(ValueError, lambda: Series(tf, tf2, tf3, tf4))
    raises(ValueError, lambda: Series(-tf3, tf2))
    raises(TypeError, lambda: Series(2, tf, tf4))
    raises(TypeError, lambda: Series(s**2 + p*s, tf3, tf2))
    raises(TypeError, lambda: Series(tf3, Matrix([1, 2, 3, 4])))

    # discrete-time tests
    ds0 = Series(dtf, dtf2)
    assert ds0.args == (dtf, dtf2)
    assert ds0.var == s
    assert ds0.sampling_time == 0.1

    ds1 = Series(Parallel(dtf, -dtf2), dtf2)
    assert ds1.args == (Parallel(dtf, -dtf2), dtf2)
    assert ds1.var == s
    assert ds1.sampling_time == 0.1

    dtf3_ = DiscreteTransferFunction(inp, 1, s, 0.1)
    dtf4_ = DiscreteTransferFunction(-out, 1, s, 0.1)
    ds2 = Series(dtf, Parallel(dtf3_, dtf4_), dtf2)
    assert ds2.args == (dtf, Parallel(dtf3_, dtf4_), dtf2)
    assert ds2.sampling_time == 0.1

    ds3 = Series(dtf, dtf2, dtf4)
    assert ds3.args == (dtf, dtf2, dtf4)
    assert ds3.sampling_time == 0.1

    ds4 = Series(dtf3_, dtf4_)
    assert ds4.args == (dtf3_, dtf4_)
    assert ds4.var == s
    assert ds4.sampling_time == 0.1

    ds6 = Series(dtf2, dtf4, Parallel(dtf2, -dtf), dtf4)
    assert ds6.args == (dtf2, dtf4, Parallel(dtf2, -dtf), dtf4)
    assert ds6.sampling_time == 0.1

    ds7 = Series(dtf, dtf2)
    assert ds0 == ds7
    assert not ds0 == ds2
    assert ds7.sampling_time == 0.1

    dtf5 = DiscreteTransferFunction(s, s-1, s, 1)
    raises(ValueError, lambda: Series(dtf, dtf3))
    raises(ValueError, lambda: Series(dtf, dtf2, dtf3, dtf4))
    raises(ValueError, lambda: Series(-dtf3, dtf2))
    raises(TypeError, lambda: Series(dtf5, dtf)) # can't do Series with different sampling time systems
    raises(TypeError, lambda: Series(dtf5, tf)) # can't do Series with continuous and discrete time systems
    raises(TypeError, lambda: Series(2, dtf, dtf4))
    raises(TypeError, lambda: Series(s**2 + p*s, dtf3, dtf2))
    raises(TypeError, lambda: Series(dtf3, Matrix([1, 2, 3, 4])))


def test_MIMOSeries_construction():
    tf_1 = TransferFunction(a0*s**3 + a1*s**2 - a2*s, b0*p**4 + b1*p**3 - b2*s*p, s)
    tf_2 = TransferFunction(a2*p - s, a2*s + p, s)
    tf_3 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)

    dtf_1 = DiscreteTransferFunction(a0*s**3 + a1*s**2 - a2*s,
                               b0*p**4 + b1*p**3 - b2*s*p, s, 2)
    dtf_2 = DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2)
    dtf_3 = DiscreteTransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s, 2)
    dTF1 = DiscreteTransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s, 2)
    dTF2 = DiscreteTransferFunction(k, 1, s, 2)
    dTF3 = DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2)

    tfm_1 = TransferFunctionMatrix([[tf_1, tf_2, tf_3], [-tf_3, -tf_2, tf_1]])
    tfm_2 = TransferFunctionMatrix([[-tf_2], [-tf_2], [-tf_3]])
    tfm_3 = TransferFunctionMatrix([[-tf_3]])
    tfm_4 = TransferFunctionMatrix([[TF3], [TF2], [-TF1]])
    tfm_5 = TransferFunctionMatrix.from_Matrix(Matrix([1/p]), p)

    dtfm_1 = TransferFunctionMatrix([[dtf_1, dtf_2, dtf_3],
                                     [-dtf_3, -dtf_2, dtf_1]])
    dtfm_2 = TransferFunctionMatrix([[-dtf_2], [-dtf_2], [-dtf_3]])
    dtfm_3 = TransferFunctionMatrix([[-dtf_3]])
    dtfm_4 = TransferFunctionMatrix([[dTF3], [dTF2], [-dTF1]])

    # continuous-time tests
    s8 = MIMOSeries(tfm_2, tfm_1)
    assert s8.args == (tfm_2, tfm_1)
    assert s8.var == s
    assert s8.shape == (s8.num_outputs, s8.num_inputs) == (2, 1)
    assert s8.is_continuous is True

    s9 = MIMOSeries(tfm_3, tfm_2, tfm_1)
    assert s9.args == (tfm_3, tfm_2, tfm_1)
    assert s9.var == s
    assert s9.shape == (s9.num_outputs, s9.num_inputs) == (2, 1)

    s11 = MIMOSeries(tfm_3, MIMOParallel(-tfm_2, -tfm_4), tfm_1)
    assert s11.args == (tfm_3, MIMOParallel(-tfm_2, -tfm_4), tfm_1)
    assert s11.shape == (s11.num_outputs, s11.num_inputs) == (2, 1)

    # discrete-time tests
    ds8 = MIMOSeries(dtfm_2, dtfm_1)
    assert ds8.args == (dtfm_2, dtfm_1)
    assert ds8.var == s
    assert ds8.shape == (ds8.num_outputs, ds8.num_inputs) == (2, 1)
    assert ds8.is_continuous is False
    assert ds8.sampling_time == 2

    ds9 = MIMOSeries(dtfm_3, dtfm_2, dtfm_1)
    assert ds9.args == (dtfm_3, dtfm_2, dtfm_1)
    assert ds9.var == s
    assert ds9.shape == (ds9.num_outputs, ds9.num_inputs) == (2, 1)

    ds11 = MIMOSeries(dtfm_3, MIMOParallel(-dtfm_2, -dtfm_4), dtfm_1)
    assert ds11.args == (dtfm_3, MIMOParallel(-dtfm_2, -dtfm_4), dtfm_1)
    assert ds11.shape == (ds11.num_outputs, ds11.num_inputs) == (2, 1)
    assert ds11.sampling_time == 2

    # arg cannot be empty tuple.
    raises(ValueError, lambda: MIMOSeries())

    # arg cannot contain SISO as well as MIMO systems.
    raises(TypeError, lambda: MIMOSeries(tfm_1, tf_1))

    # for all the adjacent transfer function matrices:
    # no. of inputs of first TFM must be equal to the no. of outputs of the second TFM.
    raises(ValueError, lambda: MIMOSeries(tfm_1, tfm_2, -tfm_1))

    # all the TFMs must use the same complex variable.
    raises(ValueError, lambda: MIMOSeries(tfm_3, tfm_5))

    # Number or expression not allowed in the arguments.
    raises(TypeError, lambda: MIMOSeries(2, tfm_2, tfm_3))
    raises(TypeError, lambda: MIMOSeries(s**2 + p*s, -tfm_2, tfm_3))
    raises(TypeError, lambda: MIMOSeries(Matrix([1/p]), tfm_3))

    # compatibility checks
    dtfm_6 = TransferFunctionMatrix.from_Matrix(Matrix(eye(3)) * s, s, 12)
    raises(TypeError, lambda: MIMOSeries(dtfm_6, dtfm_1))
    raises(TypeError, lambda: MIMOSeries(dtfm_2, tfm_1))


def test_Series_functions():
    tf1 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    tf2 = TransferFunction(k, 1, s)
    tf3 = TransferFunction(a2*p - s, a2*s + p, s)
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)

    dtf1 = DiscreteTransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s, 0.2)
    dtf2 = DiscreteTransferFunction(k, 1, s, 0.2)
    dtf3 = DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.2)
    dtf4 = DiscreteTransferFunction(a0*p + p**a1 - s, p, p, 0.2)
    dtf5 = DiscreteTransferFunction(a1*s**2 + a2*s - a0, s + a0, s, 0.2)

    # test continuous-time Transfer Functions
    assert tf1*tf2*tf3 == Series(tf1, tf2, tf3) == Series(Series(tf1, tf2), tf3) == Series(tf1, Series(tf2, tf3))
    assert tf1*(tf2 + tf3) == Series(tf1, Parallel(tf2, tf3))
    assert tf1*tf2 + tf5 == Parallel(Series(tf1, tf2), tf5)
    assert tf1*tf2 - tf5 == Parallel(Series(tf1, tf2), -tf5)
    assert tf1*tf2 + tf3 + tf5 == Parallel(Series(tf1, tf2), tf3, tf5)
    assert tf1*tf2 - tf3 - tf5 == Parallel(Series(tf1, tf2), -tf3, -tf5)
    assert tf1*tf2 - tf3 + tf5 == Parallel(Series(tf1, tf2), -tf3, tf5)
    assert tf1*tf2 + tf3*tf5 == Parallel(Series(tf1, tf2), Series(tf3, tf5))
    assert tf1*tf2 - tf3*tf5 == Parallel(Series(tf1, tf2), Series(TransferFunction(-1, 1, s), Series(tf3, tf5)))
    assert tf2*tf3*(tf2 - tf1)*tf3 == Series(tf2, tf3, Parallel(tf2, -tf1), tf3)
    assert -tf1*tf2 == Series(-tf1, tf2)
    assert -(tf1*tf2) == Series(TransferFunction(-1, 1, s), Series(tf1, tf2))
    raises(ValueError, lambda: tf1*tf2*tf4)
    raises(ValueError, lambda: tf1*(tf2 - tf4))
    raises(ValueError, lambda: tf3*Matrix([1, 2, 3]))

    # evaluate=True -> doit()
    assert Series(tf1, tf2, evaluate=True) == Series(tf1, tf2).doit() == \
        TransferFunction(k, s**2 + 2*s*wn*zeta + wn**2, s)
    assert Series(tf1, tf2, Parallel(tf1, -tf3), evaluate=True) == Series(tf1, tf2, Parallel(tf1, -tf3)).doit() == \
        TransferFunction(k*(a2*s + p + (-a2*p + s)*(s**2 + 2*s*wn*zeta + wn**2)), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2)**2, s)
    assert Series(tf2, tf1, -tf3, evaluate=True) == Series(tf2, tf1, -tf3).doit() == \
        TransferFunction(k*(-a2*p + s), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert not Series(tf1, -tf2, evaluate=False) == Series(tf1, -tf2).doit()

    assert Series(Parallel(tf1, tf2), Parallel(tf2, -tf3)).doit() == \
        TransferFunction((k*(s**2 + 2*s*wn*zeta + wn**2) + 1)*(-a2*p + k*(a2*s + p) + s), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Series(-tf1, -tf2, -tf3).doit() == \
        TransferFunction(k*(-a2*p + s), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert -Series(tf1, tf2, tf3).doit() == \
        TransferFunction(-k*(a2*p - s), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Series(tf2, tf3, Parallel(tf2, -tf1), tf3).doit() == \
        TransferFunction(k*(a2*p - s)**2*(k*(s**2 + 2*s*wn*zeta + wn**2) - 1), (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2), s)

    assert Series(tf1, tf2).rewrite(TransferFunction) == TransferFunction(k, s**2 + 2*s*wn*zeta + wn**2, s)
    assert Series(tf2, tf1, -tf3).rewrite(TransferFunction) == \
        TransferFunction(k*(-a2*p + s), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)

    S1 = Series(Parallel(tf1, tf2), Parallel(tf2, -tf3))
    assert S1.is_proper
    assert not S1.is_strictly_proper
    assert S1.is_biproper

    S2 = Series(tf1, tf2, tf3)
    assert S2.is_proper
    assert S2.is_strictly_proper
    assert not S2.is_biproper

    S3 = Series(tf1, -tf2, Parallel(tf1, -tf3))
    assert S3.is_proper
    assert S3.is_strictly_proper
    assert not S3.is_biproper

    # test discrete-time Transfer Functions
    assert dtf1*dtf2*dtf3 == Series(dtf1, dtf2, dtf3) == \
        Series(Series(dtf1, dtf2), dtf3) == Series(dtf1, Series(dtf2, dtf3))
    assert dtf1*(dtf2 + dtf3) == Series(dtf1, Parallel(dtf2, dtf3))
    assert dtf1*dtf2 + dtf5 == Parallel(Series(dtf1, dtf2), dtf5)
    assert dtf1*dtf2 - dtf5 == Parallel(Series(dtf1, dtf2), -dtf5)
    assert dtf1*dtf2 + dtf3 + dtf5 == Parallel(Series(dtf1, dtf2), dtf3, dtf5)
    assert dtf1*dtf2 - dtf3 - dtf5 == Parallel(Series(dtf1, dtf2), -dtf3, -dtf5)
    assert dtf1*dtf2 - dtf3 + dtf5 == Parallel(Series(dtf1, dtf2), -dtf3, dtf5)
    assert dtf1*dtf2 + dtf3*dtf5 == Parallel(Series(dtf1, dtf2),
                                             Series(dtf3, dtf5))
    assert dtf1*dtf2 - dtf3*dtf5 == \
        Parallel(Series(dtf1, dtf2),
        Series(DiscreteTransferFunction(-1, 1, s, 0.2), Series(dtf3, dtf5)))
    assert dtf2*dtf3*(dtf2 - dtf1)*dtf3 == Series(dtf2, dtf3,
                                                  Parallel(dtf2, -dtf1), dtf3)
    assert -dtf1*dtf2 == Series(-dtf1, dtf2)
    assert -(dtf1*dtf2) == Series(DiscreteTransferFunction(-1, 1, s, 0.2),
                                  Series(dtf1, dtf2))

    raises(ValueError, lambda: dtf1*dtf2*dtf4)
    raises(ValueError, lambda: dtf1*(dtf2 - dtf4))
    raises(ValueError, lambda: dtf3*Matrix([1, 2, 3]))

    # evaluate=True -> doit()
    assert Series(dtf1, dtf2, evaluate=True) == Series(dtf1, dtf2).doit() == \
        DiscreteTransferFunction(k, s**2 + 2*s*wn*zeta + wn**2, s, 0.2)
    assert Series(dtf1, dtf2, Parallel(dtf1, -dtf3), evaluate=True) == \
        Series(dtf1,dtf2, Parallel(dtf1, -dtf3)).doit() == \
        DiscreteTransferFunction(k*(a2*s + p + (-a2*p + s)*\
                              (s**2 + 2*s*wn*zeta + wn**2)),
                           (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2)**2, s, 0.2)
    assert Series(dtf2, dtf1, -dtf3, evaluate=True) == \
        Series(dtf2, dtf1, -dtf3).doit() == \
        DiscreteTransferFunction(k*(-a2*p + s), (a2*s + p)*\
                           (s**2 + 2*s*wn*zeta + wn**2), s, 0.2)
    assert not Series(dtf1, -dtf2, evaluate=False) == Series(dtf1, -dtf2).doit()

    assert Series(Parallel(dtf1, dtf2), Parallel(dtf2, -dtf3)).doit() == \
        DiscreteTransferFunction((k*(s**2 + 2*s*wn*zeta + wn**2) + 1)*\
                           (-a2*p +k*(a2*s + p) + s), (a2*s + p)*\
                            (s**2 + 2*s*wn*zeta + wn**2), s, 0.2)
    assert Series(-dtf1, -dtf2, -dtf3).doit() == \
        DiscreteTransferFunction(k*(-a2*p + s), (a2*s + p)*\
                           (s**2 + 2*s*wn*zeta + wn**2), s, 0.2)
    assert -Series(dtf1, dtf2, dtf3).doit() == \
        DiscreteTransferFunction(-k*(a2*p - s), (a2*s + p)*\
                           (s**2 + 2*s*wn*zeta + wn**2), s, 0.2)
    assert Series(dtf2, dtf3, Parallel(dtf2, -dtf1), dtf3).doit() == \
        DiscreteTransferFunction(k*(a2*p - s)**2*(k*(s**2 + 2*s*wn*zeta + wn**2) - 1),
                           (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2), s, 0.2)

    assert Series(dtf1, dtf2).rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction(k, s**2 + 2*s*wn*zeta + wn**2, s, 0.2)
    assert Series(dtf2, dtf1, -dtf3).rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction(k*(-a2*p + s),
                           (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2),
                           s, 0.2)

    S1 = Series(Parallel(dtf1, dtf2), Parallel(dtf2, -dtf3))
    assert S1.is_proper
    assert not S1.is_strictly_proper
    assert S1.is_biproper

    S2 = Series(dtf1, dtf2, dtf3)
    assert S2.is_proper
    assert S2.is_strictly_proper
    assert not S2.is_biproper

    S3 = Series(dtf1, -dtf2, Parallel(dtf1, -dtf3))
    assert S3.is_proper
    assert S3.is_strictly_proper
    assert not S3.is_biproper

    # test compatibility
    raises(TypeError, lambda:
           Series(dtf1, dtf2).rewrite(TransferFunction))
    raises(TypeError, lambda:
           Series(tf1, tf2).rewrite(DiscreteTransferFunction))

def test_MIMOSeries_functions():
    tfm1 = TransferFunctionMatrix([[TF1, TF2, TF3], [-TF3, -TF2, TF1]])
    tfm2 = TransferFunctionMatrix([[-TF1], [-TF2], [-TF3]])
    tfm3 = TransferFunctionMatrix([[-TF1]])
    tfm4 = TransferFunctionMatrix([[-TF2, -TF3], [-TF1, TF2]])
    tfm5 = TransferFunctionMatrix([[TF2, -TF2], [-TF3, -TF2]])
    tfm6 = TransferFunctionMatrix([[-TF3], [TF1]])
    tfm7 = TransferFunctionMatrix([[TF1], [-TF2]])

    dTF1 = DiscreteTransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s, 0.001)
    dTF2 = DiscreteTransferFunction(k, 1, s, 0.001)
    dTF3 = DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.001)

    dtfm1 = TransferFunctionMatrix([[dTF1, dTF2, dTF3], [-dTF3, -dTF2, dTF1]])
    dtfm2 = TransferFunctionMatrix([[-dTF1], [-dTF2], [-dTF3]])
    dtfm3 = TransferFunctionMatrix([[-dTF1]])
    dtfm4 = TransferFunctionMatrix([[-dTF2, -dTF3], [-dTF1, dTF2]])
    dtfm5 = TransferFunctionMatrix([[dTF2, -dTF2], [-dTF3, -dTF2]])
    dtfm6 = TransferFunctionMatrix([[-dTF3], [dTF1]])
    dtfm7 = TransferFunctionMatrix([[dTF1], [-dTF2]])

    # continuous-time tests
    assert tfm1*tfm2 + tfm6 == MIMOParallel(MIMOSeries(tfm2, tfm1), tfm6)
    assert tfm1*tfm2 + tfm7 + tfm6 == MIMOParallel(MIMOSeries(tfm2, tfm1), tfm7, tfm6)
    assert tfm1*tfm2 - tfm6 - tfm7 == MIMOParallel(MIMOSeries(tfm2, tfm1), -tfm6, -tfm7)
    assert tfm4*tfm5 + (tfm4 - tfm5) == MIMOParallel(MIMOSeries(tfm5, tfm4), tfm4, -tfm5)
    assert tfm4*-tfm6 + (-tfm4*tfm6) == MIMOParallel(MIMOSeries(-tfm6, tfm4), MIMOSeries(tfm6, -tfm4))

    raises(ValueError, lambda: tfm1*tfm2 + TF1)
    raises(TypeError, lambda: tfm1*tfm2 + a0)
    raises(TypeError, lambda: tfm4*tfm6 - (s - 1))
    raises(TypeError, lambda: tfm4*-tfm6 - 8)
    raises(TypeError, lambda: (-1 + p**5) + tfm1*tfm2)

    # Shape criteria.

    raises(TypeError, lambda: -tfm1*tfm2 + tfm4)
    raises(TypeError, lambda: tfm1*tfm2 - tfm4 + tfm5)
    raises(TypeError, lambda: tfm1*tfm2 - tfm4*tfm5)

    assert tfm1*tfm2*-tfm3 == MIMOSeries(-tfm3, tfm2, tfm1)
    assert (tfm1*-tfm2)*tfm3 == MIMOSeries(tfm3, -tfm2, tfm1)

    # Multiplication of a Series object with a SISO TF not allowed.

    raises(ValueError, lambda: tfm4*tfm5*TF1)
    raises(TypeError, lambda: tfm4*tfm5*a1)
    raises(TypeError, lambda: tfm4*-tfm5*(s - 2))
    raises(TypeError, lambda: tfm5*tfm4*9)
    raises(TypeError, lambda: (-p**3 + 1)*tfm5*tfm4)

    # Transfer function matrix in the arguments.
    assert (MIMOSeries(tfm2, tfm1, evaluate=True) == MIMOSeries(tfm2, tfm1).doit()
        == TransferFunctionMatrix(((TransferFunction(-k**2*(a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**2 + (-a2*p + s)*(a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2)**2 - (a2*s + p)**2,
        (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**2, s),),
        (TransferFunction(k**2*(a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**2 + (-a2*p + s)*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*p - s)*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2),
        (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**2, s),))))

    # doit() should not cancel poles and zeros.
    mat_1 = Matrix([[1/(1+s), (1+s)/(1+s**2+2*s)**3]])
    mat_2 = Matrix([[(1+s)], [(1+s**2+2*s)**3/(1+s)]])
    tm_1, tm_2 = TransferFunctionMatrix.from_Matrix(mat_1, s), TransferFunctionMatrix.from_Matrix(mat_2, s)
    assert (MIMOSeries(tm_2, tm_1).doit()
        == TransferFunctionMatrix(((TransferFunction(2*(s + 1)**2*(s**2 + 2*s + 1)**3, (s + 1)**2*(s**2 + 2*s + 1)**3, s),),)))
    assert MIMOSeries(tm_2, tm_1).doit().simplify() == TransferFunctionMatrix(((TransferFunction(2, 1, s),),))

    # calling doit() will expand the internal Series and Parallel objects.
    assert (MIMOSeries(-tfm3, -tfm2, tfm1, evaluate=True)
        == MIMOSeries(-tfm3, -tfm2, tfm1).doit()
        == TransferFunctionMatrix(((TransferFunction(k**2*(a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**2 + (a2*p - s)**2*(s**2 + 2*s*wn*zeta + wn**2)**2 + (a2*s + p)**2,
        (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**3, s),),
        (TransferFunction(-k**2*(a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**2 + (-a2*p + s)*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*p - s)*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2),
        (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**3, s),))))
    assert (MIMOSeries(MIMOParallel(tfm4, tfm5), tfm5, evaluate=True)
        == MIMOSeries(MIMOParallel(tfm4, tfm5), tfm5).doit()
        == TransferFunctionMatrix((
            (TransferFunction(-k*(-a2*s - p + (-a2*p + s)*\
                                  (s**2 + 2*s*wn*zeta + wn**2)),
                                  (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s),
                                  TransferFunction(k*(-a2*p - k*(a2*s + p) + s),
                                                    a2*s + p, s)),
            (TransferFunction(-k*(-a2*s - p + (-a2*p + s)*\
                                  (s**2 + 2*s*wn*zeta + wn**2)),
                             (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s),
             TransferFunction((-a2*p + s)*(-a2*p - k*(a2*s + p) + s),
                              (a2*s + p)**2, s))))
        == MIMOSeries(MIMOParallel(tfm4, tfm5), tfm5).\
            rewrite(TransferFunctionMatrix))

    # discrete-time tests
    assert dtfm1*dtfm2 + dtfm6 == MIMOParallel(MIMOSeries(dtfm2, dtfm1), dtfm6)
    assert dtfm1*dtfm2 + dtfm7 + dtfm6 == MIMOParallel(
        MIMOSeries(dtfm2, dtfm1), dtfm7, dtfm6)
    assert dtfm1*dtfm2 - dtfm6 - dtfm7 == MIMOParallel(
        MIMOSeries(dtfm2, dtfm1), -dtfm6, -dtfm7)
    assert dtfm4*dtfm5 + (dtfm4 - dtfm5) == MIMOParallel(
        MIMOSeries(dtfm5, dtfm4), dtfm4, -dtfm5)
    assert dtfm4*-dtfm6 + (-dtfm4*dtfm6) == MIMOParallel(
        MIMOSeries(-dtfm6, dtfm4), MIMOSeries(dtfm6, -dtfm4))

    raises(ValueError, lambda: dtfm1*dtfm2 + dTF1)
    raises(TypeError, lambda: dtfm1*dtfm2 + a0)
    raises(TypeError, lambda: dtfm4*dtfm6 - (s - 1))
    raises(TypeError, lambda: dtfm4*-dtfm6 - 8)
    raises(TypeError, lambda: (-1 + p**5) + dtfm1*dtfm2)

    # Shape criteria.

    raises(TypeError, lambda: -dtfm1*dtfm2 + dtfm4)
    raises(TypeError, lambda: dtfm1*dtfm2 - dtfm4 + dtfm5)
    raises(TypeError, lambda: dtfm1*dtfm2 - dtfm4*dtfm5)

    assert dtfm1*dtfm2*-dtfm3 == MIMOSeries(-dtfm3, dtfm2, dtfm1)
    assert (dtfm1*-dtfm2)*dtfm3 == MIMOSeries(dtfm3, -dtfm2, dtfm1)

    # Multiplication of a Series object with a SISO TF not allowed.

    raises(ValueError, lambda: dtfm4*dtfm5*dTF1)
    raises(TypeError, lambda: dtfm4*dtfm5*a1)
    raises(TypeError, lambda: dtfm4*-dtfm5*(s - 2))
    raises(TypeError, lambda: dtfm5*dtfm4*9)
    raises(TypeError, lambda: (-p**3 + 1)*dtfm5*dtfm4)

    # Transfer function matrix in the arguments.
    assert (MIMOSeries(dtfm2, dtfm1, evaluate=True) == \
            MIMOSeries(dtfm2, dtfm1).doit() == \
            TransferFunctionMatrix((
                (DiscreteTransferFunction(-k**2*(a2*s + p)**2*\
                                  (s**2 + 2*s*wn*zeta + wn**2)**2 + \
                                    (-a2*p + s)*(a2*p - s)*\
                                        (s**2 + 2*s*wn*zeta + wn**2)**2 - \
                                            (a2*s + p)**2,
        (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**2, s, 0.001),),
        (DiscreteTransferFunction(k**2*(a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**2 + \
                          (-a2*p + s)*(a2*s + p)*\
                            (s**2 + 2*s*wn*zeta + wn**2) + \
                                (a2*p - s)*(a2*s + p)*\
                                    (s**2 + 2*s*wn*zeta + wn**2),
        (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**2, s, 0.001),))))

    # doit() should not cancel poles and zeros.
    dmat_1 = Matrix([[1/(1+s), (1+s)/(1+s**2+2*s)**3]])
    dmat_2 = Matrix([[(1+s)], [(1+s**2+2*s)**3/(1+s)]])
    dtm_1 = TransferFunctionMatrix.from_Matrix(dmat_1, s, 0.001)
    dtm_2 = TransferFunctionMatrix.from_Matrix(dmat_2, s, 0.001)
    assert (MIMOSeries(dtm_2, dtm_1).doit()
        == TransferFunctionMatrix((
            (DiscreteTransferFunction(2*(s + 1)**2*(s**2 + 2*s + 1)**3,
                              (s + 1)**2*(s**2 + 2*s + 1)**3, s, 0.001),),)))
    assert MIMOSeries(dtm_2, dtm_1).doit().simplify() == \
        TransferFunctionMatrix(((DiscreteTransferFunction(2, 1, s, 0.001),),))

    # calling doit() will expand the internal Series and Parallel objects.
    assert (MIMOSeries(-dtfm3, -dtfm2, dtfm1, evaluate=True)
        == MIMOSeries(-dtfm3, -dtfm2, dtfm1).doit()
        == TransferFunctionMatrix((
            (DiscreteTransferFunction(k**2*(a2*s + p)**2*\
                              (s**2 + 2*s*wn*zeta + wn**2)**2 + \
                                (a2*p - s)**2*\
                                    (s**2 + 2*s*wn*zeta + wn**2)**2 + \
                                        (a2*s + p)**2,
        (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**3, s, 0.001),),
        (DiscreteTransferFunction(-k**2*(a2*s + p)**2*\
                          (s**2 + 2*s*wn*zeta + wn**2)**2 + (-a2*p + s)*\
                            (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) + \
                                (a2*p - s)*(a2*s + p)*\
                                    (s**2 + 2*s*wn*zeta + wn**2),
        (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2)**3, s, 0.001),))))
    assert (MIMOSeries(MIMOParallel(dtfm4, dtfm5), dtfm5, evaluate=True)
        == MIMOSeries(MIMOParallel(dtfm4, dtfm5), dtfm5).doit()
        == TransferFunctionMatrix((
            (DiscreteTransferFunction(-k*(-a2*s - p + (-a2*p + s)*\
                                  (s**2 + 2*s*wn*zeta + wn**2)),
                                  (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s, 0.001),
                                  DiscreteTransferFunction(k*(-a2*p - k*(a2*s + p) + s),
                                                    a2*s + p, s, 0.001)),
            (DiscreteTransferFunction(-k*(-a2*s - p + (-a2*p + s)*\
                                  (s**2 + 2*s*wn*zeta + wn**2)),
                             (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s, 0.001),
             DiscreteTransferFunction((-a2*p + s)*(-a2*p - k*(a2*s + p) + s),
                              (a2*s + p)**2, s, 0.001))))
        == MIMOSeries(MIMOParallel(dtfm4, dtfm5), dtfm5).\
            rewrite(TransferFunctionMatrix))


def test_Parallel_construction():
    tf = TransferFunction(a0*s**3 + a1*s**2 - a2*s, b0*p**4 + b1*p**3 - b2*s*p, s)
    tf2 = TransferFunction(a2*p - s, a2*s + p, s)
    tf3 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf4 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    inp = Function('X_d')(s)
    out = Function('X')(s)

    dtf = DiscreteTransferFunction(a0*s**3 + a1*s**2 - a2*s,
                             b0*p**4 + b1*p**3 - b2*s*p, s, 0.5)
    dtf2 = DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5)
    dtf3 = DiscreteTransferFunction(a0*p + p**a1 - s, p, p, 0.5)
    dtf4 = DiscreteTransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s, 0.5)

    assert tf.is_continuous is True
    assert dtf.is_continuous is False

    # continuous-time tests
    p0 = Parallel(tf, tf2)
    assert p0.args == (tf, tf2)
    assert p0.var == s

    p1 = Parallel(Series(tf, -tf2), tf2)
    assert p1.args == (Series(tf, -tf2), tf2)
    assert p1.var == s

    tf3_ = TransferFunction(inp, 1, s)
    tf4_ = TransferFunction(-out, 1, s)
    p2 = Parallel(tf, Series(tf3_, -tf4_), tf2)
    assert p2.args == (tf, Series(tf3_, -tf4_), tf2)

    p3 = Parallel(tf, tf2, tf4)
    assert p3.args == (tf, tf2, tf4)

    p4 = Parallel(tf3_, tf4_)
    assert p4.args == (tf3_, tf4_)
    assert p4.var == s

    p5 = Parallel(tf, tf2)
    assert p0 == p5
    assert not p0 == p1

    p6 = Parallel(tf2, tf4, Series(tf2, -tf4))
    assert p6.args == (tf2, tf4, Series(tf2, -tf4))

    p7 = Parallel(tf2, tf4, Series(tf2, -tf), tf4)
    assert p7.args == (tf2, tf4, Series(tf2, -tf), tf4)

    raises(ValueError, lambda: Parallel(tf, tf3))
    raises(ValueError, lambda: Parallel(tf, tf2, tf3, tf4))
    raises(ValueError, lambda: Parallel(-tf3, tf4))
    raises(TypeError, lambda: Parallel(2, tf, tf4))
    raises(TypeError, lambda: Parallel(s**2 + p*s, tf3, tf2))
    raises(TypeError, lambda: Parallel(tf3, Matrix([1, 2, 3, 4])))

    # discrete-time tests
    dp0 = Parallel(dtf, dtf2)
    assert dp0.args == (dtf, dtf2)
    assert dp0.var == s
    assert dp0.sampling_time == 0.5

    dp1 = Parallel(Series(dtf, -dtf2), dtf2)
    assert dp1.args == (Series(dtf, -dtf2), dtf2)
    assert dp1.var == s
    assert dp1.sampling_time == 0.5


    dtf3_ = DiscreteTransferFunction(inp, 1, s, 0.5)
    dtf4_ = DiscreteTransferFunction(-out, 1, s, 0.5)
    dp2 = Parallel(dtf, Series(dtf3_, -dtf4_), dtf2)
    assert dp2.args == (dtf, Series(dtf3_, -dtf4_), dtf2)
    assert dp2.sampling_time == 0.5

    dp3 = Parallel(dtf, dtf2, dtf4)
    assert dp3.args == (dtf, dtf2, dtf4)
    assert dp3.sampling_time == 0.5

    dp4 = Parallel(dtf3_, dtf4_)
    assert dp4.args == (dtf3_, dtf4_)
    assert dp4.var == s
    assert dp4.sampling_time == 0.5

    dp5 = Parallel(dtf, dtf2)
    assert dp0 == dp5
    assert not dp0 == dp1
    assert dp5.sampling_time == 0.5

    dp6 = Parallel(dtf2, dtf4, Series(dtf2, -dtf4))
    assert dp6.args == (dtf2, dtf4, Series(dtf2, -dtf4))
    assert dp6.sampling_time == 0.5

    dp7 = Parallel(dtf2, dtf4, Series(dtf2, -dtf), dtf4)
    assert dp7.args == (dtf2, dtf4, Series(dtf2, -dtf), dtf4)
    assert dp7.sampling_time == 0.5

    dtf5 = DiscreteTransferFunction(s, s-1, s, 1)
    raises(ValueError, lambda: Parallel(dtf, dtf3))
    raises(ValueError, lambda: Parallel(dtf, dtf2, dtf3, dtf4))
    raises(ValueError, lambda: Parallel(-dtf3, dtf4))
    raises(TypeError, lambda: Parallel(dtf5, dtf)) # can't do Parallel with different sampling time systems
    raises(TypeError, lambda: Parallel(dtf5, tf)) # can't do Parallel with continuous and discrete time systems
    raises(TypeError, lambda: Parallel(2, dtf, dtf4))
    raises(TypeError, lambda: Parallel(s**2 + p*s, dtf3, dtf2))
    raises(TypeError, lambda: Parallel(dtf3, Matrix([1, 2, 3, 4])))

def test_MIMOParallel_construction():
    tfm1 = TransferFunctionMatrix([[TF1], [TF2], [TF3]])
    tfm2 = TransferFunctionMatrix([[-TF3], [TF2], [TF1]])
    tfm3 = TransferFunctionMatrix([[TF1]])
    tfm4 = TransferFunctionMatrix([[TF2], [TF1], [TF3]])
    tfm5 = TransferFunctionMatrix([[TF1, TF2], [TF2, TF1]])
    tfm6 = TransferFunctionMatrix([[TF2, TF1], [TF1, TF2]])
    tfm7 = TransferFunctionMatrix.from_Matrix(Matrix([[1/p]]), p)

    dTF1 = DiscreteTransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s, 0.5)
    dTF2 = DiscreteTransferFunction(k, 1, s, 0.5)
    dTF3 = DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.5)

    dtfm1 = TransferFunctionMatrix([[dTF1], [dTF2], [dTF3]])
    dtfm2 = TransferFunctionMatrix([[-dTF3], [dTF2], [dTF1]])
    dtfm3 = TransferFunctionMatrix([[dTF1]])
    dtfm4 = TransferFunctionMatrix([[dTF2], [dTF1], [dTF3]])
    dtfm5 = TransferFunctionMatrix([[dTF1, dTF2], [dTF2, dTF1]])
    dtfm6 = TransferFunctionMatrix([[dTF2, dTF1], [dTF1, dTF2]])

    # continuous-time tests
    p8 = MIMOParallel(tfm1, tfm2)
    assert p8.args == (tfm1, tfm2)
    assert p8.var == s
    assert p8.shape == (p8.num_outputs, p8.num_inputs) == (3, 1)
    assert p8.is_continuous is True

    p9 = MIMOParallel(MIMOSeries(tfm3, tfm1), tfm2)
    assert p9.args == (MIMOSeries(tfm3, tfm1), tfm2)
    assert p9.var == s
    assert p9.shape == (p9.num_outputs, p9.num_inputs) == (3, 1)

    p10 = MIMOParallel(tfm1, MIMOSeries(tfm3, tfm4), tfm2)
    assert p10.args == (tfm1, MIMOSeries(tfm3, tfm4), tfm2)
    assert p10.var == s
    assert p10.shape == (p10.num_outputs, p10.num_inputs) == (3, 1)

    p11 = MIMOParallel(tfm2, tfm1, tfm4)
    assert p11.args == (tfm2, tfm1, tfm4)
    assert p11.shape == (p11.num_outputs, p11.num_inputs) == (3, 1)

    p12 = MIMOParallel(tfm6, tfm5)
    assert p12.args == (tfm6, tfm5)
    assert p12.shape == (p12.num_outputs, p12.num_inputs) == (2, 2)

    p13 = MIMOParallel(tfm2, tfm4, MIMOSeries(-tfm3, tfm4), -tfm4)
    assert p13.args == (tfm2, tfm4, MIMOSeries(-tfm3, tfm4), -tfm4)
    assert p13.shape == (p13.num_outputs, p13.num_inputs) == (3, 1)

    # discrete-time tests
    dp8 = MIMOParallel(dtfm1, dtfm2)
    assert dp8.args == (dtfm1, dtfm2)
    assert dp8.var == s
    assert dp8.shape == (dp8.num_outputs, dp8.num_inputs) == (3, 1)
    assert dp8.sampling_time == 0.5
    assert dp8.is_continuous is False

    dp9 = MIMOParallel(MIMOSeries(dtfm3, dtfm1), dtfm2)
    assert dp9.args == (MIMOSeries(dtfm3, dtfm1), dtfm2)
    assert dp9.var == s
    assert dp9.shape == (dp9.num_outputs, dp9.num_inputs) == (3, 1)

    dp10 = MIMOParallel(dtfm1, MIMOSeries(dtfm3, dtfm4), dtfm2)
    assert dp10.args == (dtfm1, MIMOSeries(dtfm3, dtfm4), dtfm2)
    assert dp10.var == s
    assert dp10.shape == (dp10.num_outputs, dp10.num_inputs) == (3, 1)

    dp11 = MIMOParallel(dtfm2, dtfm1, dtfm4)
    assert dp11.args == (dtfm2, dtfm1, dtfm4)
    assert dp11.shape == (dp11.num_outputs, dp11.num_inputs) == (3, 1)

    dp12 = MIMOParallel(dtfm6, dtfm5)
    assert dp12.args == (dtfm6, dtfm5)
    assert dp12.shape == (dp12.num_outputs, dp12.num_inputs) == (2, 2)

    dp13 = MIMOParallel(dtfm2, dtfm4, MIMOSeries(-dtfm3, dtfm4), -dtfm4)
    assert dp13.args == (dtfm2, dtfm4, MIMOSeries(-dtfm3, dtfm4), -dtfm4)
    assert dp13.shape == (dp13.num_outputs, dp13.num_inputs) == (3, 1)

    # arg cannot be empty tuple.
    raises(TypeError, lambda: MIMOParallel(()))

    # arg cannot contain SISO as well as MIMO systems.
    raises(TypeError, lambda: MIMOParallel(tfm1, tfm2, TF1))

    # all TFMs must have same shapes.
    raises(TypeError, lambda: MIMOParallel(tfm1, tfm3, tfm4))

    # all TFMs must be using the same complex variable.
    raises(ValueError, lambda: MIMOParallel(tfm3, tfm7))

    # Number or expression not allowed in the arguments.
    raises(TypeError, lambda: MIMOParallel(2, tfm1, tfm4))
    raises(TypeError, lambda: MIMOParallel(s**2 + p*s, -tfm4, tfm2))

    # compatibility checks
    dtfm_6 = TransferFunctionMatrix.from_Matrix(Matrix(eye(2)) * s, s, T)
    raises(TypeError, lambda: MIMOParallel(dtfm_6, dtfm1))
    raises(TypeError, lambda: MIMOParallel(dtfm1, tfm2))

def test_Parallel_functions():
    tf1 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    tf2 = TransferFunction(k, 1, s)
    tf3 = TransferFunction(a2*p - s, a2*s + p, s)
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)

    dtf1 = DiscreteTransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s, 0.01)
    dtf2 = DiscreteTransferFunction(k, 1, s, 0.01)
    dtf3 = DiscreteTransferFunction(a2*p - s, a2*s + p, s, 0.01)
    dtf4 = DiscreteTransferFunction(a0*p + p**a1 - s, p, p, 0.01)
    dtf5 = DiscreteTransferFunction(a1*s**2 + a2*s - a0, s + a0, s, 0.01)

    # continuous-time tests
    assert tf1 + tf2 + tf3 == Parallel(tf1, tf2, tf3)
    assert tf1 + tf2 + tf3 + tf5 == Parallel(tf1, tf2, tf3, tf5)
    assert tf1 + tf2 - tf3 - tf5 == Parallel(tf1, tf2, -tf3, -tf5)
    assert tf1 + tf2*tf3 == Parallel(tf1, Series(tf2, tf3))
    assert tf1 - tf2*tf3 == Parallel(tf1, -Series(tf2,tf3))
    assert -tf1 - tf2 == Parallel(-tf1, -tf2)
    assert -(tf1 + tf2) == Series(TransferFunction(-1, 1, s), Parallel(tf1, tf2))
    assert (tf2 + tf3)*tf1 == Series(Parallel(tf2, tf3), tf1)
    assert (tf1 + tf2)*(tf3*tf5) == Series(Parallel(tf1, tf2), tf3, tf5)
    assert -(tf2 + tf3)*-tf5 == Series(TransferFunction(-1, 1, s), Parallel(tf2, tf3), -tf5)
    assert tf2 + tf3 + tf2*tf1 + tf5 == Parallel(tf2, tf3, Series(tf2, tf1), tf5)
    assert tf2 + tf3 + tf2*tf1 - tf3 == Parallel(tf2, tf3, Series(tf2, tf1), -tf3)
    assert (tf1 + tf2 + tf5)*(tf3 + tf5) == Series(Parallel(tf1, tf2, tf5), Parallel(tf3, tf5))
    raises(ValueError, lambda: tf1 + tf2 + tf4)
    raises(ValueError, lambda: tf1 - tf2*tf4)
    raises(ValueError, lambda: tf3 + Matrix([1, 2, 3]))

    # evaluate=True -> doit()
    assert Parallel(tf1, tf2, evaluate=True) == Parallel(tf1, tf2).doit() == \
        TransferFunction(k*(s**2 + 2*s*wn*zeta + wn**2) + 1, s**2 + 2*s*wn*zeta + wn**2, s)
    assert Parallel(tf1, tf2, Series(-tf1, tf3), evaluate=True) == \
        Parallel(tf1, tf2, Series(-tf1, tf3)).doit() == TransferFunction(k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2)**2 + \
            (-a2*p + s)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), (a2*s + p)*(s**2 + \
                2*s*wn*zeta + wn**2)**2, s)
    assert Parallel(tf2, tf1, -tf3, evaluate=True) == Parallel(tf2, tf1, -tf3).doit() == \
        TransferFunction(a2*s + k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) + p + (-a2*p + s)*(s**2 + 2*s*wn*zeta + wn**2) \
            , (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert not Parallel(tf1, -tf2, evaluate=False) == Parallel(tf1, -tf2).doit()

    assert Parallel(Series(tf1, tf2), Series(tf2, tf3)).doit() == \
        TransferFunction(k*(a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) + k*(a2*s + p), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Parallel(-tf1, -tf2, -tf3).doit() == \
        TransferFunction(-a2*s - k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) - p + (-a2*p + s)*(s**2 + 2*s*wn*zeta + wn**2), \
            (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert -Parallel(tf1, tf2, tf3).doit() == \
        TransferFunction(-a2*s - k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) - p - (a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2), \
            (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Parallel(tf2, tf3, Series(tf2, -tf1), tf3).doit() == \
        TransferFunction(k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) - k*(a2*s + p) + (2*a2*p - 2*s)*(s**2 + 2*s*wn*zeta \
            + wn**2), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)

    assert Parallel(tf1, tf2).rewrite(TransferFunction) == \
        TransferFunction(k*(s**2 + 2*s*wn*zeta + wn**2) + 1, s**2 + 2*s*wn*zeta + wn**2, s)
    assert Parallel(tf2, tf1, -tf3).rewrite(TransferFunction) == \
        TransferFunction(a2*s + k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) + p + (-a2*p + s)*(s**2 + 2*s*wn*zeta + \
             wn**2), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)

    assert Parallel(tf1, Parallel(tf2, tf3)) == Parallel(tf1, tf2, tf3) == Parallel(Parallel(tf1, tf2), tf3)

    P1 = Parallel(Series(tf1, tf2), Series(tf2, tf3))
    assert P1.is_proper
    assert not P1.is_strictly_proper
    assert P1.is_biproper

    P2 = Parallel(tf1, -tf2, -tf3)
    assert P2.is_proper
    assert not P2.is_strictly_proper
    assert P2.is_biproper

    P3 = Parallel(tf1, -tf2, Series(tf1, tf3))
    assert P3.is_proper
    assert not P3.is_strictly_proper
    assert P3.is_biproper

    # discrete-time tests
    assert dtf1 + dtf2 + dtf3 == Parallel(dtf1, dtf2, dtf3)
    assert dtf1 + dtf2 + dtf3 + dtf5 == Parallel(dtf1, dtf2, dtf3, dtf5)
    assert dtf1 + dtf2 - dtf3 - dtf5 == Parallel(dtf1, dtf2, -dtf3, -dtf5)
    assert dtf1 + dtf2*dtf3 == Parallel(dtf1, Series(dtf2, dtf3))
    assert dtf1 - dtf2*dtf3 == Parallel(dtf1, -Series(dtf2,dtf3))
    assert -dtf1 - dtf2 == Parallel(-dtf1, -dtf2)
    assert -(dtf1 + dtf2) == Series(DiscreteTransferFunction(-1, 1, s, 0.01),
                                    Parallel(dtf1, dtf2))
    assert (dtf2 + dtf3)*dtf1 == Series(Parallel(dtf2, dtf3), dtf1)
    assert (dtf1 + dtf2)*(dtf3*dtf5) == Series(Parallel(dtf1, dtf2), dtf3, dtf5)
    assert -(dtf2 + dtf3)*-dtf5 == Series(DiscreteTransferFunction(-1, 1, s, 0.01),
                                       Parallel(dtf2, dtf3), -dtf5)
    assert dtf2 + dtf3 + dtf2*dtf1 + dtf5 == Parallel(dtf2, dtf3,
                                                      Series(dtf2, dtf1), dtf5)
    assert dtf2 + dtf3 + dtf2*dtf1 - dtf3 == \
        Parallel(dtf2, dtf3, Series(dtf2, dtf1), -dtf3)
    assert (dtf1 + dtf2 + dtf5)*(dtf3 + dtf5) == \
        Series(Parallel(dtf1, dtf2, dtf5), Parallel(dtf3, dtf5))

    raises(ValueError, lambda: dtf1 + dtf2 + dtf4)
    raises(ValueError, lambda: dtf1 - dtf2*dtf4)
    raises(ValueError, lambda: dtf3 + Matrix([1, 2, 3]))

    # evaluate=True -> doit()
    assert Parallel(dtf1, dtf2, evaluate=True) == \
        Parallel(dtf1, dtf2).doit() == \
        DiscreteTransferFunction(k*(s**2 + 2*s*wn*zeta + wn**2) + 1,
                           s**2 + 2*s*wn*zeta + wn**2, s, 0.01)
    assert Parallel(dtf1, dtf2, Series(-dtf1, dtf3), evaluate=True) == \
        Parallel(dtf1, dtf2, Series(-dtf1, dtf3)).doit() == \
        DiscreteTransferFunction(k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2)**2 + \
            (-a2*p + s)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*s + p)*\
                (s**2 + 2*s*wn*zeta + wn**2), (a2*s + p)*(s**2 + \
                2*s*wn*zeta + wn**2)**2, s, 0.01)
    assert Parallel(dtf2, dtf1, -dtf3, evaluate=True) == \
        Parallel(dtf2, dtf1, -dtf3).doit() == \
        DiscreteTransferFunction(a2*s + k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) +\
                            p + (-a2*p + s)*(s**2 + 2*s*wn*zeta + wn**2) \
                        , (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s, 0.01)
    assert not Parallel(dtf1, -dtf2, evaluate=False) == \
          Parallel(dtf1, -dtf2).doit()

    assert Parallel(Series(dtf1, dtf2), Series(dtf2, dtf3)).doit() == \
        DiscreteTransferFunction(k*(a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) +\
                            k*(a2*s + p),
                          (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s, 0.01)
    assert Parallel(-dtf1, -dtf2, -dtf3).doit() == \
        DiscreteTransferFunction(-a2*s - k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) - \
                            p + (-a2*p + s)*(s**2 + 2*s*wn*zeta + wn**2), \
            (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s, 0.01)
    assert -Parallel(dtf1, dtf2, dtf3).doit() == \
        DiscreteTransferFunction(-a2*s - k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) - \
                            p - (a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2), \
            (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s, 0.01)
    assert Parallel(dtf2, dtf3, Series(dtf2, -dtf1), dtf3).doit() == \
        DiscreteTransferFunction(k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) - \
                            k*(a2*s + p) + (2*a2*p - 2*s)*(s**2 + 2*s*wn*zeta \
            + wn**2), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s, 0.01)

    assert Parallel(dtf1, dtf2).rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction(k*(s**2 + 2*s*wn*zeta + wn**2) + 1,
                            s**2 + 2*s*wn*zeta + wn**2, s, 0.01)
    assert Parallel(dtf2, dtf1, -dtf3).rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction(a2*s + k*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) + \
                            p + (-a2*p + s)*(s**2 + 2*s*wn*zeta + \
             wn**2), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s, 0.01)

    assert Parallel(dtf1, Parallel(dtf2, dtf3)) == \
           Parallel(dtf1, dtf2, dtf3) == Parallel(Parallel(dtf1, dtf2), dtf3)

    P1 = Parallel(Series(dtf1, dtf2), Series(dtf2, dtf3))
    assert P1.is_proper
    assert not P1.is_strictly_proper
    assert P1.is_biproper

    P2 = Parallel(dtf1, -dtf2, -dtf3)
    assert P2.is_proper
    assert not P2.is_strictly_proper
    assert P2.is_biproper

    P3 = Parallel(dtf1, -dtf2, Series(dtf1, dtf3))
    assert P3.is_proper
    assert not P3.is_strictly_proper
    assert P3.is_biproper

    # test compatibility
    raises(TypeError, lambda:
           Parallel(dtf1, dtf2).rewrite(TransferFunction))
    raises(TypeError, lambda:
           Parallel(tf1, tf2).rewrite(DiscreteTransferFunction))


def test_MIMOParallel_functions():
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)

    dtf4 = DiscreteTransferFunction(a0*p + p**a1 - s, p, p, Rational(1,2))
    dtf5 = DiscreteTransferFunction(a1*s**2 + a2*s - a0, s + a0, s, Rational(1,2))

    tfm1 = TransferFunctionMatrix([[TF1], [TF2], [TF3]])
    tfm2 = TransferFunctionMatrix([[-TF2], [tf5], [-TF1]])
    tfm3 = TransferFunctionMatrix([[tf5], [-tf5], [TF2]])
    tfm4 = TransferFunctionMatrix([[TF2, -tf5], [TF1, tf5]])
    tfm5 = TransferFunctionMatrix([[TF1, TF2], [TF3, -tf5]])
    tfm6 = TransferFunctionMatrix([[-TF2]])
    tfm7 = TransferFunctionMatrix([[tf4], [-tf4], [tf4]])

    dTF1 = DiscreteTransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s, Rational(1,2))
    dTF2 = DiscreteTransferFunction(k, 1, s, Rational(1,2))
    dTF3 = DiscreteTransferFunction(a2*p - s, a2*s + p, s, Rational(1,2))

    dtfm1 = TransferFunctionMatrix([[dTF1], [dTF2], [dTF3]])
    dtfm2 = TransferFunctionMatrix([[-dTF2], [dtf5], [-dTF1]])
    dtfm3 = TransferFunctionMatrix([[dtf5], [-dtf5], [dTF2]])
    dtfm4 = TransferFunctionMatrix([[dTF2, -dtf5], [dTF1, dtf5]])
    dtfm5 = TransferFunctionMatrix([[dTF1, dTF2], [dTF3, -dtf5]])
    dtfm6 = TransferFunctionMatrix([[-dTF2]])
    dtfm7 = TransferFunctionMatrix([[dtf4], [-dtf4], [dtf4]])

    # continuous-time tests
    assert tfm1 + tfm2 + tfm3 == MIMOParallel(tfm1, tfm2, tfm3) == MIMOParallel(MIMOParallel(tfm1, tfm2), tfm3)
    assert tfm2 - tfm1 - tfm3 == MIMOParallel(tfm2, -tfm1, -tfm3)
    assert tfm2 - tfm3 + (-tfm1*tfm6*-tfm6) == MIMOParallel(tfm2, -tfm3, MIMOSeries(-tfm6, tfm6, -tfm1))
    assert tfm1 + tfm1 - (-tfm1*tfm6) == MIMOParallel(tfm1, tfm1, -MIMOSeries(tfm6, -tfm1))
    assert tfm2 - tfm3 - tfm1 + tfm2 == MIMOParallel(tfm2, -tfm3, -tfm1, tfm2)
    assert tfm1 + tfm2 - tfm3 - tfm1 == MIMOParallel(tfm1, tfm2, -tfm3, -tfm1)
    raises(ValueError, lambda: tfm1 + tfm2 + TF2)
    raises(TypeError, lambda: tfm1 - tfm2 - a1)
    raises(TypeError, lambda: tfm2 - tfm3 - (s - 1))
    raises(TypeError, lambda: -tfm3 - tfm2 - 9)
    raises(TypeError, lambda: (1 - p**3) - tfm3 - tfm2)
    # All TFMs must use the same complex var. tfm7 uses 'p'.
    raises(ValueError, lambda: tfm3 - tfm2 - tfm7)
    raises(ValueError, lambda: tfm2 - tfm1 + tfm7)
    # (tfm1 +/- tfm2) has (3, 1) shape while tfm4 has (2, 2) shape.
    raises(TypeError, lambda: tfm1 + tfm2 + tfm4)
    raises(TypeError, lambda: (tfm1 - tfm2) - tfm4)

    assert (tfm1 + tfm2)*tfm6 == MIMOSeries(tfm6, MIMOParallel(tfm1, tfm2))
    assert (tfm2 - tfm3)*tfm6*-tfm6 == MIMOSeries(-tfm6, tfm6, MIMOParallel(tfm2, -tfm3))
    assert (tfm2 - tfm1 - tfm3)*(tfm6 + tfm6) == MIMOSeries(MIMOParallel(tfm6, tfm6), MIMOParallel(tfm2, -tfm1, -tfm3))
    raises(ValueError, lambda: (tfm4 + tfm5)*TF1)
    raises(TypeError, lambda: (tfm2 - tfm3)*a2)
    raises(TypeError, lambda: (tfm3 + tfm2)*(s - 6))
    raises(TypeError, lambda: (tfm1 + tfm2 + tfm3)*0)
    raises(TypeError, lambda: (1 - p**3)*(tfm1 + tfm3))

    # (tfm3 - tfm2) has (3, 1) shape while tfm4*tfm5 has (2, 2) shape.
    raises(ValueError, lambda: (tfm3 - tfm2)*tfm4*tfm5)
    # (tfm1 - tfm2) has (3, 1) shape while tfm5 has (2, 2) shape.
    raises(ValueError, lambda: (tfm1 - tfm2)*tfm5)

    # TFM in the arguments.
    assert (MIMOParallel(tfm1, tfm2, evaluate=True) == MIMOParallel(tfm1, tfm2).doit()
    == MIMOParallel(tfm1, tfm2).rewrite(TransferFunctionMatrix)
    == TransferFunctionMatrix((
        (TransferFunction(-k*(s**2 + 2*s*wn*zeta + wn**2) + 1,
                          s**2 + 2*s*wn*zeta + wn**2, s),), \
        (TransferFunction(-a0 + a1*s**2 + a2*s + k*(a0 + s), a0 + s, s),),
        (TransferFunction(-a2*s - p + (a2*p - s)* \
        (s**2 + 2*s*wn*zeta + wn**2),
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s),))))

    # discrete-time tests
    assert dtfm1 + dtfm2 + dtfm3 == MIMOParallel(dtfm1, dtfm2, dtfm3) == \
        MIMOParallel(MIMOParallel(dtfm1, dtfm2), dtfm3)
    assert dtfm2 - dtfm1 - dtfm3 == MIMOParallel(dtfm2, -dtfm1, -dtfm3)
    assert dtfm2 - dtfm3 + (-dtfm1*dtfm6*-dtfm6) == \
        MIMOParallel(dtfm2, -dtfm3, MIMOSeries(-dtfm6, dtfm6, -dtfm1))
    assert dtfm1 + dtfm1 - (-dtfm1*dtfm6) == \
        MIMOParallel(dtfm1, dtfm1, -MIMOSeries(dtfm6, -dtfm1))
    assert dtfm2 - dtfm3 - dtfm1 + dtfm2 == MIMOParallel(dtfm2, -dtfm3,
                                                         -dtfm1, dtfm2)
    assert dtfm1 + dtfm2 - dtfm3 - dtfm1 == MIMOParallel(dtfm1, dtfm2,
                                                         -dtfm3, -dtfm1)
    raises(ValueError, lambda: dtfm1 + dtfm2 + dTF2)
    raises(TypeError, lambda: dtfm1 - dtfm2 - a1)
    raises(TypeError, lambda: dtfm2 - dtfm3 - (s - 1))
    raises(TypeError, lambda: -dtfm3 - dtfm2 - 9)
    raises(TypeError, lambda: (1 - p**3) - dtfm3 - dtfm2)
    # All TFMs must use the same complex var. dtfm7 uses 'p'.
    raises(ValueError, lambda: dtfm3 - dtfm2 - dtfm7)
    raises(ValueError, lambda: dtfm2 - dtfm1 + dtfm7)
    # (dtfm1 +/- dtfm2) has (3, 1) shape while dtfm4 has (2, 2) shape.
    raises(TypeError, lambda: dtfm1 + dtfm2 + dtfm4)
    raises(TypeError, lambda: (dtfm1 - dtfm2) - dtfm4)

    assert (dtfm1 + dtfm2)*dtfm6 == MIMOSeries(dtfm6, MIMOParallel(dtfm1,
                                                                   dtfm2))
    assert (dtfm2 - dtfm3)*dtfm6*-dtfm6 == \
        MIMOSeries(-dtfm6, dtfm6, MIMOParallel(dtfm2, -dtfm3))
    assert (dtfm2 - dtfm1 - dtfm3)*(dtfm6 + dtfm6) == \
        MIMOSeries(MIMOParallel(dtfm6, dtfm6), MIMOParallel(dtfm2, -dtfm1,
                                                            -dtfm3))
    raises(ValueError, lambda: (dtfm4 + dtfm5)*dTF1)
    raises(TypeError, lambda: (dtfm2 - dtfm3)*a2)
    raises(TypeError, lambda: (dtfm3 + dtfm2)*(s - 6))
    raises(TypeError, lambda: (dtfm1 + dtfm2 + dtfm3)*0)
    raises(TypeError, lambda: (1 - p**3)*(dtfm1 + dtfm3))

    # (dtfm3 - dtfm2) has (3, 1) shape while dtfm4*dtfm5 has (2, 2) shape.
    raises(ValueError, lambda: (dtfm3 - dtfm2)*dtfm4*dtfm5)
    # (dtfm1 - dtfm2) has (3, 1) shape while dtfm5 has (2, 2) shape.
    raises(ValueError, lambda: (dtfm1 - dtfm2)*dtfm5)

    # TFM in the arguments.
    assert (MIMOParallel(dtfm1, dtfm2, evaluate=True) == \
            MIMOParallel(dtfm1, dtfm2).doit()
    == MIMOParallel(dtfm1, dtfm2).rewrite(TransferFunctionMatrix)
    == TransferFunctionMatrix((
        (DiscreteTransferFunction(-k*(s**2 + 2*s*wn*zeta + wn**2) + 1,
                          s**2 + 2*s*wn*zeta + wn**2, s, Rational(1,2)),), \
        (DiscreteTransferFunction(-a0 + a1*s**2 + a2*s + k*(a0 + s), a0 + s, s,
         Rational(1,2)),),
        (DiscreteTransferFunction(-a2*s - p + (a2*p - s)* \
        (s**2 + 2*s*wn*zeta + wn**2),
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s, Rational(1,2)),))))


def test_Feedback_construction():
    tf1 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    tf2 = TransferFunction(k, 1, s)
    tf3 = TransferFunction(a2*p - s, a2*s + p, s)
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)
    tf6 = TransferFunction(s - p, p + s, p)

    dtf1 = DiscreteTransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s, 2)
    dtf2 = DiscreteTransferFunction(k, 1, s, 2)
    dtf3 = DiscreteTransferFunction(a2*p - s, a2*s + p, s, 2)
    dtf4 = DiscreteTransferFunction(a0*p + p**a1 - s, p, p, 2)
    dtf5 = DiscreteTransferFunction(a1*s**2 + a2*s - a0, s + a0, s, 2)
    dtf6 = DiscreteTransferFunction(s - p, p + s, p, 2)

    assert tf1.is_continuous is True
    assert dtf1.is_continuous is False

    # continuous-time tests
    f1 = Feedback(TransferFunction(1, 1, s), tf1*tf2*tf3)
    assert f1.args == (TransferFunction(1, 1, s), Series(tf1, tf2, tf3), -1)
    assert f1.sys1 == TransferFunction(1, 1, s)
    assert f1.sys2 == Series(tf1, tf2, tf3)
    assert f1.var == s

    f2 = Feedback(tf1, tf2*tf3)
    assert f2.args == (tf1, Series(tf2, tf3), -1)
    assert f2.sys1 == tf1
    assert f2.sys2 == Series(tf2, tf3)
    assert f2.var == s

    f3 = Feedback(tf1*tf2, tf5)
    assert f3.args == (Series(tf1, tf2), tf5, -1)
    assert f3.sys1 == Series(tf1, tf2)

    f4 = Feedback(tf4, tf6)
    assert f4.args == (tf4, tf6, -1)
    assert f4.sys1 == tf4
    assert f4.var == p

    f5 = Feedback(tf5, TransferFunction(1, 1, s))
    assert f5.args == (tf5, TransferFunction(1, 1, s), -1)
    assert f5.var == s
    assert f5 == Feedback(tf5)  # When sys2 is not passed explicitly, it is assumed to be unit tf.

    f6 = Feedback(TransferFunction(1, 1, p), tf4)
    assert f6.args == (TransferFunction(1, 1, p), tf4, -1)
    assert f6.var == p

    f7 = -Feedback(tf4*tf6, TransferFunction(1, 1, p))
    assert f7.args == (Series(TransferFunction(-1, 1, p), Series(tf4, tf6)), -TransferFunction(1, 1, p), -1)
    assert f7.sys1 == Series(TransferFunction(-1, 1, p), Series(tf4, tf6))

    # denominator can't be a Parallel instance
    raises(TypeError, lambda: Feedback(tf1, tf2 + tf3))
    raises(TypeError, lambda: Feedback(tf1, Matrix([1, 2, 3])))
    raises(TypeError, lambda: Feedback(TransferFunction(1, 1, s), s - 1))
    raises(TypeError, lambda: Feedback(1, 1))
    # raises(ValueError, lambda: Feedback(TransferFunction(1, 1, s), TransferFunction(1, 1, s)))
    raises(ValueError, lambda: Feedback(tf2, tf4*tf5))
    raises(ValueError, lambda: Feedback(tf2, tf1, 1.5))  # `sign` can only be -1 or 1
    raises(ValueError, lambda: Feedback(tf1, -tf1**-1))  # denominator can't be zero
    raises(ValueError, lambda: Feedback(tf4, tf5))  # Both systems should use the same `var`

    # discrete-time tests
    df1 = Feedback(DiscreteTransferFunction(1, 1, s, 2), dtf1*dtf2*dtf3)
    assert df1.args == (DiscreteTransferFunction(1, 1, s, 2), Series(dtf1, dtf2, dtf3),
                        -1)
    assert df1.sys1 == DiscreteTransferFunction(1, 1, s, 2)
    assert df1.sys2 == Series(dtf1, dtf2, dtf3)
    assert df1.var == s

    df2 = Feedback(dtf1, dtf2*dtf3)
    assert df2.args == (dtf1, Series(dtf2, dtf3), -1)
    assert df2.sys1 == dtf1
    assert df2.sys2 == Series(dtf2, dtf3)
    assert df2.var == s
    assert df2.sampling_time == 2

    df3 = Feedback(dtf1*dtf2, dtf5)
    assert df3.args == (Series(dtf1, dtf2), dtf5, -1)
    assert df3.sys1 == Series(dtf1, dtf2)
    assert df3.sampling_time == 2

    df4 = Feedback(dtf4, dtf6)
    assert df4.args == (dtf4, dtf6, -1)
    assert df4.sys1 == dtf4
    assert df4.var == p
    assert df4.sampling_time == 2

    df5 = Feedback(dtf5, DiscreteTransferFunction(1, 1, s, 2))
    assert df5.args == (dtf5, DiscreteTransferFunction(1, 1, s, 2), -1)
    assert df5.var == s
    assert df5 == Feedback(dtf5)  # When sys2 is not passed explicitly, it is assumed to be unit tf.
    assert df5.sampling_time == 2

    df6 = Feedback(DiscreteTransferFunction(1, 1, p, 2), dtf4)
    assert df6.args == (DiscreteTransferFunction(1, 1, p, 2), dtf4, -1)
    assert df6.var == p
    assert df6.sampling_time == 2

    df7 = -Feedback(dtf4*dtf6, DiscreteTransferFunction(1, 1, p, 2))
    assert df7.args == (Series(DiscreteTransferFunction(-1, 1, p, 2),
                               Series(dtf4, dtf6)),
                               -DiscreteTransferFunction(1, 1, p, 2), -1)
    assert df7.sys1 == Series(DiscreteTransferFunction(-1, 1, p, 2),
                              Series(dtf4, dtf6))
    assert df7.sampling_time == 2

    # denominator can't be a Parallel instance
    raises(TypeError, lambda: Feedback(dtf1, dtf2 + dtf3))
    raises(TypeError, lambda: Feedback(dtf1, Matrix([1, 2, 3])))
    raises(TypeError, lambda: Feedback(DiscreteTransferFunction(1, 1, s), s - 1))
    raises(TypeError, lambda: Feedback(1, 1))
    dtf_ = DiscreteTransferFunction(s, s + 1, s, 0.1)
    raises(TypeError, lambda: Feedback(dtf_, dtf1)) # can't do Feedback with different sampling time systems
    raises(TypeError, lambda: Feedback(f1, df1)) # can't do Feedback with continuous and discrete time systems
    raises(ValueError, lambda: Feedback(dtf2, dtf4*dtf5))
    raises(ValueError, lambda: Feedback(dtf2, dtf1, 1.5))  # `sign` can only be -1 or 1
    raises(ValueError, lambda: Feedback(dtf1, -dtf1**-1))  # denominator can't be zero
    raises(ValueError, lambda: Feedback(dtf4, dtf5))  # Both systems should use the same `var`


def test_Feedback_functions():
    tf = TransferFunction(1, 1, s)
    tf1 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    tf2 = TransferFunction(k, 1, s)
    tf3 = TransferFunction(a2*p - s, a2*s + p, s)
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)
    tf6 = TransferFunction(s - p, p + s, p)

    dtf = DiscreteTransferFunction(1, 1, s, 10)
    dtf1 = DiscreteTransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s, 10)
    dtf2 = DiscreteTransferFunction(k, 1, s, 10)
    dtf3 = DiscreteTransferFunction(a2*p - s, a2*s + p, s, 10)
    dtf4 = DiscreteTransferFunction(a0*p + p**a1 - s, p, p, 10)
    dtf5 = DiscreteTransferFunction(a1*s**2 + a2*s - a0, s + a0, s, 10)
    dtf6 = DiscreteTransferFunction(s - p, p + s, p, 10)

    # continuous-time tests
    assert (tf1*tf2*tf3 / tf3*tf5) == Series(tf1, tf2, tf3, pow(tf3, -1), tf5)
    assert (tf1*tf2*tf3) / (tf3*tf5) == Series((tf1*tf2*tf3).doit(), pow((tf3*tf5).doit(),-1))
    assert tf / (tf + tf1) == Feedback(tf, tf1)
    assert tf / (tf + tf1*tf2*tf3) == Feedback(tf, tf1*tf2*tf3)
    assert tf1 / (tf + tf1*tf2*tf3) == Feedback(tf1, tf2*tf3)
    assert (tf1*tf2) / (tf + tf1*tf2) == Feedback(tf1*tf2, tf)
    assert (tf1*tf2) / (tf + tf1*tf2*tf5) == Feedback(tf1*tf2, tf5)
    assert (tf1*tf2) / (tf + tf1*tf2*tf5*tf3) in (Feedback(tf1*tf2, tf5*tf3), Feedback(tf1*tf2, tf3*tf5))
    assert tf4 / (TransferFunction(1, 1, p) + tf4*tf6) == Feedback(tf4, tf6)
    assert tf5 / (tf + tf5) == Feedback(tf5, tf)

    raises(TypeError, lambda: tf1*tf2*tf3 / (1 + tf1*tf2*tf3))
    raises(ValueError, lambda: tf2*tf3 / (tf + tf2*tf3*tf4))

    assert Feedback(tf, tf1*tf2*tf3).doit() == \
        TransferFunction((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), k*(a2*p - s) + \
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(tf, tf1*tf2*tf3).sensitivity == \
        1/(k*(a2*p - s)/((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2)) + 1)
    assert Feedback(tf1, tf2*tf3).doit() == \
        TransferFunction((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), (k*(a2*p - s) + \
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2))*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(tf1, tf2*tf3).sensitivity == \
        1/(k*(a2*p - s)/((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2)) + 1)
    assert Feedback(tf1*tf2, tf5).doit() == \
        TransferFunction(k*(a0 + s)*(s**2 + 2*s*wn*zeta + wn**2), (k*(-a0 + a1*s**2 + a2*s) + \
        (a0 + s)*(s**2 + 2*s*wn*zeta + wn**2))*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(tf1*tf2, tf5, 1).sensitivity == \
        1/(-k*(-a0 + a1*s**2 + a2*s)/((a0 + s)*(s**2 + 2*s*wn*zeta + wn**2)) + 1)
    assert Feedback(tf4, tf6).doit() == \
        TransferFunction(p*(p + s)*(a0*p + p**a1 - s), p*(p*(p + s) + (-p + s)*(a0*p + p**a1 - s)), p)
    assert -Feedback(tf4*tf6, TransferFunction(1, 1, p)).doit() == \
        TransferFunction(-p*(-p + s)*(p + s)*(a0*p + p**a1 - s), p*(p + s)*(p*(p + s) + (-p + s)*(a0*p + p**a1 - s)), p)
    assert Feedback(tf, tf).doit() == TransferFunction(1, 2, s)

    assert Feedback(tf1, tf2*tf5).rewrite(TransferFunction) == \
        TransferFunction((a0 + s)*(s**2 + 2*s*wn*zeta + wn**2),
                         (k*(-a0 + a1*s**2 + a2*s) + (a0 + s)* \
                          (s**2 + 2*s*wn*zeta + wn**2))* \
                            (s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(TransferFunction(1, 1, p), tf4).\
        rewrite(TransferFunction) == TransferFunction(p, a0*p + p + p**a1 - s,
                                                      p)

    # discrete-time tests
    assert (dtf1*dtf2*dtf3 / dtf3*dtf5) == Series(dtf1, dtf2, dtf3,
                                                  pow(dtf3, -1), dtf5)
    assert (dtf1*dtf2*dtf3) / (dtf3*dtf5) == Series((dtf1*dtf2*dtf3).doit(),
                                                    pow((dtf3*dtf5).doit(),-1))
    assert dtf / (dtf + dtf1) == Feedback(dtf, dtf1)
    assert dtf / (dtf + dtf1*dtf2*dtf3) == Feedback(dtf, dtf1*dtf2*dtf3)
    assert dtf1 / (dtf + dtf1*dtf2*dtf3) == Feedback(dtf1, dtf2*dtf3)
    assert (dtf1*dtf2) / (dtf + dtf1*dtf2) == Feedback(dtf1*dtf2, dtf)
    assert (dtf1*dtf2) / (dtf + dtf1*dtf2*dtf5) == Feedback(dtf1*dtf2, dtf5)
    assert (dtf1*dtf2) / (dtf + dtf1*dtf2*dtf5*dtf3) in (
        Feedback(dtf1*dtf2, dtf5*dtf3), Feedback(dtf1*dtf2, dtf3*dtf5))
    assert dtf4 / (DiscreteTransferFunction(1, 1, p, 10) + dtf4*dtf6) == \
        Feedback(dtf4, dtf6)
    assert dtf5 / (dtf + dtf5) == Feedback(dtf5, dtf)

    raises(TypeError, lambda: dtf1*dtf2*dtf3 / (1 + dtf1*dtf2*dtf3))
    raises(ValueError, lambda: dtf2*dtf3 / (dtf + dtf2*dtf3*dtf4))

    assert Feedback(dtf, dtf1*dtf2*dtf3).doit() == \
        DiscreteTransferFunction((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2),
                           k*(a2*p - s) + (a2*s + p)* \
                            (s**2 + 2*s*wn*zeta + wn**2), s, 10)
    assert Feedback(dtf, dtf1*dtf2*dtf3).sensitivity == \
        1/(k*(a2*p - s)/((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2)) + 1)
    assert Feedback(dtf1, dtf2*dtf3).doit() == \
        DiscreteTransferFunction((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2),
                           (k*(a2*p - s) + (a2*s + p)* \
                            (s**2 + 2*s*wn*zeta + wn**2))*\
                            (s**2 + 2*s*wn*zeta + wn**2), s, 10)
    assert Feedback(dtf1, dtf2*dtf3).sensitivity == \
        1/(k*(a2*p - s)/((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2)) + 1)
    assert Feedback(dtf1*dtf2, dtf5).doit() == \
        DiscreteTransferFunction(k*(a0 + s)*(s**2 + 2*s*wn*zeta + wn**2),
                           (k*(-a0 + a1*s**2 + a2*s) + (a0 + s)* \
                            (s**2 + 2*s*wn*zeta + wn**2))* \
                                (s**2 + 2*s*wn*zeta + wn**2), s, 10)
    assert Feedback(dtf1*dtf2, dtf5, 1).sensitivity == \
        1/(-k*(-a0 + a1*s**2 + a2*s)/((a0 + s)* \
        (s**2 + 2*s*wn*zeta + wn**2)) + 1)
    assert Feedback(dtf4, dtf6).doit() == \
        DiscreteTransferFunction(p*(p + s)*(a0*p + p**a1 - s),
                           p*(p*(p + s) + (-p + s)*(a0*p + p**a1 - s)), p, 10)
    assert -Feedback(dtf4*dtf6, DiscreteTransferFunction(1, 1, p, 10)).doit() == \
        DiscreteTransferFunction(-p*(-p + s)*(p + s)*(a0*p + p**a1 - s),
                           p*(p + s)*(p*(p + s) + (-p + s)*(a0*p + p**a1 - s)),
                           p, 10)
    assert Feedback(dtf, dtf).doit() == DiscreteTransferFunction(1, 2, s, 10)

    assert Feedback(dtf1, dtf2*dtf5).rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction((a0 + s)*(s**2 + 2*s*wn*zeta + wn**2),
                           (k*(-a0 + a1*s**2 + a2*s) + (a0 + s)* \
                            (s**2 + 2*s*wn*zeta + wn**2))* \
                            (s**2 + 2*s*wn*zeta + wn**2), s, 10)
    assert Feedback(DiscreteTransferFunction(1, 1, p, 10),
                    dtf4).rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction(p, a0*p + p + p**a1 - s, p, 10)

    # test compatibility
    raises(TypeError, lambda:
           Feedback(dtf1, dtf2).rewrite(TransferFunction))
    raises(TypeError, lambda:
           Feedback(tf1, tf2).rewrite(DiscreteTransferFunction))

def test_Feedback_with_Series():
    # Solves issue https://github.com/sympy/sympy/issues/26161
    tf1 = TransferFunction(s+1, 1, s)
    tf2 = TransferFunction(s+2, 1, s)
    fd1 = Feedback(tf1, tf2, -1) # Negative Feedback system
    fd2 = Feedback(tf1, tf2, 1) # Positive Feedback system
    unit = TransferFunction(1, 1, s)

    dtf1 = DiscreteTransferFunction(s+1, 1, s)
    dtf2 = DiscreteTransferFunction(s+2, 1, s)
    dfd1 = Feedback(dtf1, dtf2, -1) # Negative Feedback system
    dfd2 = Feedback(dtf1, dtf2, 1) # Positive Feedback system
    dunit = DiscreteTransferFunction(1, 1, s)

    assert tf1.is_continuous is True
    assert dtf1.is_continuous is False

    # continuous-time tests
    # Checking the type
    assert isinstance(fd1, SISOLinearTimeInvariant)
    assert isinstance(fd1, Feedback)

    # Testing the numerator and denominator
    assert fd1.num == tf1
    assert fd2.num == tf1
    assert fd1.den == Parallel(unit, Series(tf2, tf1))
    assert fd2.den == Parallel(unit, -Series(tf2, tf1))

    # Testing the Series and Parallel Combination with Feedback and TransferFunction
    s1 = Series(tf1, fd1)
    p1 = Parallel(tf1, fd1)
    assert tf1 * fd1 == s1
    assert tf1 + fd1 == p1
    assert s1.doit() == TransferFunction((s + 1)**2, (s + 1)*(s + 2) + 1, s)
    assert p1.doit() == TransferFunction(s + (s + 1)*((s + 1)*(s + 2) + 1) + 1, (s + 1)*(s + 2) + 1, s)

    # Testing the use of Feedback and TransferFunction with Feedback
    fd3 = Feedback(tf1*fd1, tf2, -1)
    assert fd3 == Feedback(Series(tf1, fd1), tf2)
    assert fd3.num == tf1 * fd1
    assert fd3.den == Parallel(unit, Series(tf2, Series(tf1, fd1)))

    # Testing the use of Feedback and TransferFunction with TransferFunction
    tf3 = TransferFunction(tf1*fd1, tf2, s)
    assert tf3 == TransferFunction(Series(tf1, fd1), tf2, s)
    assert tf3.num == tf1*fd1

    # discrete-time tests
    # Checking the type
    assert isinstance(dfd1, SISOLinearTimeInvariant)
    assert isinstance(dfd1, Feedback)

    # Testing the numerator and denominator
    assert dfd1.num == dtf1
    assert dfd2.num == dtf1
    assert dfd1.den == Parallel(dunit, Series(dtf2, dtf1))
    assert dfd2.den == Parallel(dunit, -Series(dtf2, dtf1))

    # Testing the Series and Parallel Combination with Feedback and DiscreteTransferFunction
    ds1 = Series(dtf1, dfd1)
    dp1 = Parallel(dtf1, dfd1)
    assert dtf1 * dfd1 == ds1
    assert dtf1 + dfd1 == dp1
    assert ds1.doit() == DiscreteTransferFunction((s + 1)**2, (s + 1)*(s + 2) + 1, s)
    assert dp1.doit() == DiscreteTransferFunction(s + (s + 1)* \
                                            ((s + 1)*(s + 2) + 1) + 1,
                                            (s + 1)*(s + 2) + 1, s)

    # Testing the use of Feedback and DiscreteTransferFunction with Feedback
    dfd3 = Feedback(dtf1*dfd1, dtf2, -1)
    assert dfd3 == Feedback(Series(dtf1, dfd1), dtf2)
    assert dfd3.num == dtf1 * dfd1
    assert dfd3.den == Parallel(dunit, Series(dtf2, Series(dtf1, dfd1)))

    # Testing the use of Feedback and DiscreteTransferFunction with DiscreteTransferFunction
    dtf3 = DiscreteTransferFunction(dtf1*dfd1, dtf2, s)
    assert dtf3 == DiscreteTransferFunction(Series(dtf1, dfd1), dtf2, s)
    assert dtf3.num == dtf1*dfd1


def test_issue_26161():
    # Issue https://github.com/sympy/sympy/issues/26161
    Ib, Is, m, h, l2, l1 = symbols('I_b, I_s, m, h, l2, l1',
                                            real=True, nonnegative=True)
    KD, KP, v = symbols('K_D, K_P, v', real=True)

    tau1_sq = (Ib + m * h ** 2) / m / g / h
    tau2 = l2 / v
    tau3 = v / (l1 + l2)
    K = v ** 2 / g / (l1 + l2)

    Gtheta = TransferFunction(-K * (tau2 * s + 1), tau1_sq * s ** 2 - 1, s)
    Gdelta = TransferFunction(1, Is * s ** 2 + c * s, s)
    Gpsi = TransferFunction(1, tau3 * s, s)
    Dcont = TransferFunction(KD * s, 1, s)
    PIcont = TransferFunction(KP, s, s)
    Gunity = TransferFunction(1, 1, s)

    Ginner = Feedback(Dcont * Gdelta, Gtheta)
    Gouter = Feedback(PIcont * Ginner * Gpsi, Gunity)
    assert Gouter == Feedback(Series(PIcont, Series(Ginner, Gpsi)), Gunity)
    assert Gouter.num == Series(PIcont, Series(Ginner, Gpsi))
    assert Gouter.den == Parallel(Gunity, Series(Gunity, Series(PIcont, Series(Ginner, Gpsi))))
    expr = (KD*KP*g*s**3*v**2*(l1 + l2)*(Is*s**2 + c*s)**2*(-g*h*m + s**2*(Ib + h**2*m))*(-KD*g*h*m*s*v**2*(l2*s + v) + \
            g*v*(l1 + l2)*(Is*s**2 + c*s)*(-g*h*m + s**2*(Ib + h**2*m))))/((s**2*v*(Is*s**2 + c*s)*(-KD*g*h*m*s*v**2* \
            (l2*s + v) + g*v*(l1 + l2)*(Is*s**2 + c*s)*(-g*h*m + s**2*(Ib + h**2*m)))*(KD*KP*g*s*v*(l1 + l2)**2* \
            (Is*s**2 + c*s)*(-g*h*m + s**2*(Ib + h**2*m)) + s**2*v*(Is*s**2 + c*s)*(-KD*g*h*m*s*v**2*(l2*s + v) + \
            g*v*(l1 + l2)*(Is*s**2 + c*s)*(-g*h*m + s**2*(Ib + h**2*m))))/(l1 + l2)))

    assert (Gouter.to_expr() - expr).simplify() == 0


def test_MIMOFeedback_construction():
    tf1 = TransferFunction(1, s, s)
    tf2 = TransferFunction(s, s**3 - 1, s)
    tf3 = TransferFunction(s, s + 1, s)
    tf4 = TransferFunction(s, s**2 + 1, s)

    dtf1 = DiscreteTransferFunction(1, s, s, 0.5)
    dtf2 = DiscreteTransferFunction(s, s**3 - 1, s, 0.5)
    dtf3 = DiscreteTransferFunction(s, s + 1, s, 0.5)
    dtf4 = DiscreteTransferFunction(s, s**2 + 1, s, 0.5)

    tfm_1 = TransferFunctionMatrix([[tf1, tf2], [tf3, tf4]])
    tfm_2 = TransferFunctionMatrix([[tf2, tf3], [tf4, tf1]])
    tfm_3 = TransferFunctionMatrix([[tf3, tf4], [tf1, tf2]])

    dtfm_1 = TransferFunctionMatrix([[dtf1, dtf2], [dtf3, dtf4]])
    dtfm_2 = TransferFunctionMatrix([[dtf2, dtf3], [dtf4, dtf1]])
    dtfm_3 = TransferFunctionMatrix([[dtf3, dtf4], [dtf1, dtf2]])

    # continuous-time tests
    f1 = MIMOFeedback(tfm_1, tfm_2)
    assert f1.args == (tfm_1, tfm_2, -1)
    assert f1.sys1 == tfm_1
    assert f1.sys2 == tfm_2
    assert f1.var == s
    assert f1.sign == -1
    assert -(-f1) == f1
    assert f1.is_continuous is True

    f2 = MIMOFeedback(tfm_2, tfm_1, 1)
    assert f2.args == (tfm_2, tfm_1, 1)
    assert f2.sys1 == tfm_2
    assert f2.sys2 == tfm_1
    assert f2.var == s
    assert f2.sign == 1
    assert f2.is_continuous is True

    f3 = MIMOFeedback(tfm_1, MIMOSeries(tfm_3, tfm_2))
    assert f3.args == (tfm_1, MIMOSeries(tfm_3, tfm_2), -1)
    assert f3.sys1 == tfm_1
    assert f3.sys2 == MIMOSeries(tfm_3, tfm_2)
    assert f3.var == s
    assert f3.sign == -1
    assert f3.is_continuous is True

    mat = Matrix([[1, 1/s], [0, 1]])
    sys1 = controller = TransferFunctionMatrix.from_Matrix(mat, s)
    f4 = MIMOFeedback(sys1, controller)
    assert f4.args == (sys1, controller, -1)
    assert f4.sys1 == f4.sys2 == sys1
    assert f4.is_continuous is True

    # discrete-time tests
    df1 = MIMOFeedback(dtfm_1, dtfm_2)
    assert df1.args == (dtfm_1, dtfm_2, -1)
    assert df1.sys1 == dtfm_1
    assert df1.sys2 == dtfm_2
    assert df1.var == s
    assert df1.sign == -1
    assert df1.sampling_time == 0.5
    assert -(-df1) == df1
    assert df1.is_continuous is False

    df2 = MIMOFeedback(dtfm_2, dtfm_1, 1)
    assert df2.args == (dtfm_2, dtfm_1, 1)
    assert df2.sys1 == dtfm_2
    assert df2.sys2 == dtfm_1
    assert df2.var == s
    assert df2.sign == 1
    assert df2.sampling_time == 0.5
    assert df2.is_continuous is False

    df3 = MIMOFeedback(dtfm_1, MIMOSeries(dtfm_3, dtfm_2))
    assert df3.args == (dtfm_1, MIMOSeries(dtfm_3, dtfm_2), -1)
    assert df3.sys1 == dtfm_1
    assert df3.sys2 == MIMOSeries(dtfm_3, dtfm_2)
    assert df3.var == s
    assert df3.sign == -1
    assert df3.sampling_time == 0.5
    assert df3.is_continuous is False

    dmat = Matrix([[1, 1/s], [0, 1]])
    dsys1 = dcontroller = TransferFunctionMatrix.from_Matrix(dmat, s, 0.5)

    df4 = MIMOFeedback(dsys1, dcontroller)
    assert df4.args == (dsys1, dcontroller, -1)
    assert df4.sys1 == df4.sys2 == dsys1
    assert df4.sampling_time == 0.5
    assert df4.is_continuous is False

    # test compatibility
    raises(TypeError, lambda: MIMOFeedback(dtfm_1, tfm_2))
    dtfm_4 = TransferFunctionMatrix.from_Matrix(eye(2) * s, var=s,
                                                sampling_time=0.1)
    raises(TypeError, lambda: MIMOFeedback(dtfm_1, dtfm_4))


def test_MIMOFeedback_errors():
    tf1 = TransferFunction(1, s, s)
    tf2 = TransferFunction(s, s**3 - 1, s)
    tf3 = TransferFunction(s, s - 1, s)
    tf4 = TransferFunction(s, s**2 + 1, s)
    tf5 = TransferFunction(1, 1, s)
    tf6 = TransferFunction(-1, s - 1, s)

    dtf1 = DiscreteTransferFunction(1, s, s, 0.5)
    dtf2 = DiscreteTransferFunction(s, s**3 - 1, s, 0.5)
    dtf3 = DiscreteTransferFunction(s, s - 1, s, 0.5)
    dtf4 = DiscreteTransferFunction(s, s**2 + 1, s, 0.5)
    dtf5 = DiscreteTransferFunction(1, 1, s, 0.5)
    dtf6 = DiscreteTransferFunction(-1, s - 1, s, 0.5)

    tfm_1 = TransferFunctionMatrix([[tf1, tf2], [tf3, tf4]])
    tfm_2 = TransferFunctionMatrix([[tf2, tf3], [tf4, tf1]])
    tfm_3 = TransferFunctionMatrix.from_Matrix(eye(2), var=s)
    tfm_4 = TransferFunctionMatrix([[tf1, tf5], [tf5, tf5]])
    tfm_5 = TransferFunctionMatrix([[-tf3, tf3], [tf3, tf6]])
    # tfm_4 is inverse of tfm_5. Therefore tfm_5*tfm_4 = I
    tfm_6 = TransferFunctionMatrix([[-tf3]])
    tfm_7 = TransferFunctionMatrix([[tf3, tf4]])

    dtfm_1 = TransferFunctionMatrix([[dtf1, dtf2], [dtf3, dtf4]])
    dtfm_2 = TransferFunctionMatrix([[dtf2, dtf3], [dtf4, dtf1]])
    dtfm_3 = TransferFunctionMatrix.from_Matrix(eye(2), var=s)
    dtfm_4 = TransferFunctionMatrix([[dtf1, dtf5], [dtf5, dtf5]])
    dtfm_5 = TransferFunctionMatrix([[-dtf3, dtf3], [dtf3, dtf6]])
    # dtfm_4 is inverse of dtfm_5. Therefore dtfm_5*dtfm_4 = I
    dtfm_6 = TransferFunctionMatrix([[-dtf3]])
    dtfm_7 = TransferFunctionMatrix([[dtf3, dtf4]])

    # continuous-time tests
    # Unsupported Types
    raises(TypeError, lambda: MIMOFeedback(tf1, tf2))
    raises(TypeError, lambda: MIMOFeedback(MIMOParallel(tfm_1, tfm_2), tfm_3))
    # Shape Errors
    raises(ValueError, lambda: MIMOFeedback(tfm_1, tfm_6, 1))
    raises(ValueError, lambda: MIMOFeedback(tfm_7, tfm_7))
    # sign not 1/-1
    raises(ValueError, lambda: MIMOFeedback(tfm_1, tfm_2, -2))
    # Non-Invertible Systems
    raises(ValueError, lambda: MIMOFeedback(tfm_5, tfm_4, 1))
    raises(ValueError, lambda: MIMOFeedback(tfm_4, -tfm_5))
    raises(ValueError, lambda: MIMOFeedback(tfm_3, tfm_3, 1))
    # Variable not same in both the systems
    tfm_8 = TransferFunctionMatrix.from_Matrix(eye(2), var=p)
    raises(ValueError, lambda: MIMOFeedback(tfm_1, tfm_8, 1))

    # ciscrete-time tests
    # Unsupported Types
    raises(TypeError, lambda: MIMOFeedback(dtf1, dtf2))
    raises(TypeError, lambda: MIMOFeedback(MIMOParallel(dtfm_1, dtfm_2), dtfm_3))
    # Shape Errors
    raises(ValueError, lambda: MIMOFeedback(dtfm_1, dtfm_6, 1))
    raises(ValueError, lambda: MIMOFeedback(dtfm_7, dtfm_7))
    # sign not 1/-1
    raises(ValueError, lambda: MIMOFeedback(dtfm_1, dtfm_2, -2))
    # Non-Invertible Systems
    raises(ValueError, lambda: MIMOFeedback(dtfm_5, dtfm_4, 1))
    raises(ValueError, lambda: MIMOFeedback(dtfm_4, -dtfm_5))
    raises(ValueError, lambda: MIMOFeedback(dtfm_3, dtfm_3, 1))
    # Variable not same in both the systems
    dtfm_8 = TransferFunctionMatrix.from_Matrix(eye(2), var=p)
    raises(ValueError, lambda: MIMOFeedback(dtfm_1, dtfm_8, 1))
    # Mixed continuous/discrete-time systems
    raises(TypeError, lambda: MIMOFeedback(tfm_1, dtfm_2))
    raises(TypeError, lambda: MIMOFeedback(dtfm_1, tfm_2))

def test_MIMOFeedback_functions():
    tf1 = TransferFunction(1, s, s)
    tf2 = TransferFunction(s, s - 1, s)
    tf3 = TransferFunction(1, 1, s)
    tf4 = TransferFunction(-1, s - 1, s)

    dtf1 = DiscreteTransferFunction(1, s, s, 0.4)
    dtf2 = DiscreteTransferFunction(s, s - 1, s, 0.4)
    dtf3 = DiscreteTransferFunction(1, 1, s, 0.4)
    dtf4 = DiscreteTransferFunction(-1, s - 1, s, 0.4)

    tfm_1 = TransferFunctionMatrix.from_Matrix(eye(2), var=s)
    tfm_2 = TransferFunctionMatrix([[tf1, tf3], [tf3, tf3]])
    tfm_3 = TransferFunctionMatrix([[-tf2, tf2], [tf2, tf4]])
    tfm_4 = TransferFunctionMatrix([[tf1, tf2], [-tf2, tf1]])

    dtfm_1 = TransferFunctionMatrix.from_Matrix(eye(2), var=s,
                                                sampling_time=0.4)
    dtfm_2 = TransferFunctionMatrix([[dtf1, dtf3], [dtf3, dtf3]])
    dtfm_3 = TransferFunctionMatrix([[-dtf2, dtf2], [dtf2, dtf4]])
    dtfm_4 = TransferFunctionMatrix([[dtf1, dtf2], [-dtf2, dtf1]])

    # continuous-time tests
    # sensitivity, doit(), rewrite()
    F_1 = MIMOFeedback(tfm_2, tfm_3)
    F_2 = MIMOFeedback(tfm_2, MIMOSeries(tfm_4, -tfm_1), 1)

    assert F_1.sensitivity == Matrix([[S.Half, 0], [0, S.Half]])
    assert F_2.sensitivity == Matrix([[(-2*s**4 + s**2)/(s**2 - s + 1),
        (2*s**3 - s**2)/(s**2 - s + 1)], [-s**2, s]])

    assert F_1.doit() == \
        TransferFunctionMatrix(((TransferFunction(1, 2*s, s),
        TransferFunction(1, 2, s)), (TransferFunction(1, 2, s),
        TransferFunction(1, 2, s)))) == F_1.rewrite(TransferFunctionMatrix)
    assert F_2.doit(cancel=False, expand=True) == \
        TransferFunctionMatrix(((TransferFunction(-s**5 + 2*s**4 - 2*s**3 + s**2, s**5 - 2*s**4 + 3*s**3 - 2*s**2 + s, s),
        TransferFunction(-2*s**4 + 2*s**3, s**2 - s + 1, s)), (TransferFunction(0, 1, s), TransferFunction(-s**2 + s, 1, s))))
    assert F_2.doit(cancel=False) == \
        TransferFunctionMatrix(((TransferFunction(s*(2*s**3 - s**2)*(s**2 - s + 1) + \
        (-2*s**4 + s**2)*(s**2 - s + 1), s*(s**2 - s + 1)**2, s), TransferFunction(-2*s**4 + 2*s**3, s**2 - s + 1, s)),
        (TransferFunction(0, 1, s), TransferFunction(-s**2 + s, 1, s))))
    assert F_2.doit() == \
        TransferFunctionMatrix(((TransferFunction(s*(-2*s**2 + s*(2*s - 1) + 1), s**2 - s + 1, s),
        TransferFunction(-2*s**3*(s - 1), s**2 - s + 1, s)), (TransferFunction(0, 1, s), TransferFunction(s*(1 - s), 1, s))))
    assert F_2.doit(expand=True) == \
        TransferFunctionMatrix(((TransferFunction(-s**2 + s, s**2 - s + 1, s), TransferFunction(-2*s**4 + 2*s**3, s**2 - s + 1, s)),
        (TransferFunction(0, 1, s), TransferFunction(-s**2 + s, 1, s))))

    assert -(F_1.doit()) == (-F_1).doit()  # First negating then calculating vs calculating then negating.

    # continuous-time tests
    dF_1 = MIMOFeedback(dtfm_2, dtfm_3)
    dF_2 = MIMOFeedback(dtfm_2, MIMOSeries(dtfm_4, -dtfm_1), 1)

    assert dF_1.sampling_time == 0.4
    assert dF_2.sampling_time == 0.4

    assert dF_1.sensitivity == Matrix([[S.Half, 0], [0, S.Half]])
    assert dF_2.sensitivity == Matrix([[(-2*s**4 + s**2)/(s**2 - s + 1),
        (2*s**3 - s**2)/(s**2 - s + 1)], [-s**2, s]])

    assert dF_1.doit() == \
        TransferFunctionMatrix(((
            DiscreteTransferFunction(1, 2*s, s, 0.4),
            DiscreteTransferFunction(1, 2, s, 0.4)),
            (DiscreteTransferFunction(1, 2, s, 0.4),
             DiscreteTransferFunction(1, 2, s, 0.4)))) == \
        dF_1.rewrite(TransferFunctionMatrix)
    assert dF_2.doit(cancel=False, expand=True) == \
        TransferFunctionMatrix(((
            DiscreteTransferFunction(-s**5 + 2*s**4 - 2*s**3 + s**2, s**5 - 2*s**4 + \
                             3*s**3 - 2*s**2 + s, s, 0.4),
            DiscreteTransferFunction(-2*s**4 + 2*s**3, s**2 - s + 1, s, 0.4)),
            (DiscreteTransferFunction(0, 1, s, 0.4),
             DiscreteTransferFunction(-s**2 + s, 1, s, 0.4))))
    assert dF_2.doit(cancel=False) == \
        TransferFunctionMatrix(((
            DiscreteTransferFunction(s*(2*s**3 - s**2)*(s**2 - s + 1) + \
                             (-2*s**4 + s**2)*(s**2 - s + 1),
                             s*(s**2 - s + 1)**2, s, 0.4),
            DiscreteTransferFunction(-2*s**4 + 2*s**3, s**2 - s + 1, s, 0.4)),
            (DiscreteTransferFunction(0, 1, s, 0.4),
             DiscreteTransferFunction(-s**2 + s, 1, s, 0.4))))
    assert dF_2.doit() == \
        TransferFunctionMatrix(((
            DiscreteTransferFunction(s*(-2*s**2 + s*(2*s - 1) + 1),
                               s**2 - s + 1, s, 0.4),
            DiscreteTransferFunction(-2*s**3*(s - 1), s**2 - s + 1, s, 0.4)),
            (DiscreteTransferFunction(0, 1, s, 0.4), DiscreteTransferFunction(s*(1 - s),
                                                                  1, s, 0.4))))
    assert dF_2.doit(expand=True) == \
        TransferFunctionMatrix(((
            DiscreteTransferFunction(-s**2 + s, s**2 - s + 1, s, 0.4),
            DiscreteTransferFunction(-2*s**4 + 2*s**3, s**2 - s + 1, s, 0.4)),
            (DiscreteTransferFunction(0, 1, s, 0.4),
             DiscreteTransferFunction(-s**2 + s, 1, s, 0.4))))

    assert -(dF_1.doit()) == (-dF_1).doit()  # First negating then calculating vs calculating then negating.


def test_TransferFunctionMatrix_construction():
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)

    dtf1 = DiscreteTransferFunction(a1*z**2 + a2*z - a0, z + a0, z, 2)
    dtf2 = DiscreteTransferFunction(z, z - 1, z, 2)
    dtf3 = DiscreteTransferFunction(a0*p + p**a1 - s, p, p, 3)

    tfm3_ = TransferFunctionMatrix([[-TF3]])
    assert tfm3_.shape == (tfm3_.num_outputs, tfm3_.num_inputs) == (1, 1)
    assert tfm3_.args == Tuple(Tuple(Tuple(-TF3)))
    assert tfm3_.var == s
    assert tfm3_.sampling_time == 0

    tfm5 = TransferFunctionMatrix([[TF1, -TF2], [TF3, tf5]])
    assert tfm5.shape == (tfm5.num_outputs, tfm5.num_inputs) == (2, 2)
    assert tfm5.args == Tuple(Tuple(Tuple(TF1, -TF2), Tuple(TF3, tf5)))
    assert tfm5.var == s

    tfm7 = TransferFunctionMatrix([[TF1, TF2], [TF3, -tf5], [-tf5, TF2]])
    assert tfm7.shape == (tfm7.num_outputs, tfm7.num_inputs) == (3, 2)
    assert tfm7.args == Tuple(Tuple(Tuple(TF1, TF2), Tuple(TF3, -tf5), Tuple(-tf5, TF2)))
    assert tfm7.var == s

    dtfm1 = TransferFunctionMatrix([[dtf1, dtf2], [-dtf1 + dtf2, dtf2 * dtf1]])
    assert dtfm1.shape == (dtfm1.num_outputs, dtfm1.num_inputs) == (2, 2)
    assert dtfm1.args == Tuple(Tuple(Tuple(dtf1, dtf2),
                                     Tuple(-dtf1 + dtf2, dtf2 * dtf1)
                                    ))
    assert dtfm1.var == z
    assert dtfm1.sampling_time == 2

    assert tfm3_.is_continuous is True
    assert dtfm1.is_continuous is False

    # all transfer functions will use the same complex variable. tf4 uses 'p'.
    raises(ValueError, lambda: TransferFunctionMatrix([[TF1], [TF2], [tf4]]))
    raises(ValueError, lambda: TransferFunctionMatrix([[TF1, tf4], [TF3, tf5]]))

    # length of all the lists in the TFM should be equal.
    raises(ValueError, lambda: TransferFunctionMatrix([[TF1], [TF3, tf5]]))
    raises(ValueError, lambda: TransferFunctionMatrix([[TF1, TF3], [tf5]]))

    # lists should only support transfer functions in them.
    raises(TypeError, lambda: TransferFunctionMatrix([[TF1, TF2], [TF3, Matrix([1, 2])]]))
    raises(TypeError, lambda: TransferFunctionMatrix([[TF1, Matrix([1, 2])], [TF3, TF2]]))

    # `arg` should strictly be nested list of TransferFunction
    raises(ValueError, lambda: TransferFunctionMatrix([TF1, TF2, tf5]))
    raises(ValueError, lambda: TransferFunctionMatrix([TF1]))

    # all transfer function should be either continuous or discrete.
    raises(ValueError, lambda: TransferFunctionMatrix([[TF1, dtf2],
                                                       [TF3, tf5]]))
    # all transfer function should have the same sampling time.
    raises(ValueError, lambda: TransferFunctionMatrix([[dtf1, dtf2],
                                                       [dtf3, dtf1]]))

def test_TransferFunctionMatrix_functions():
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)

    #  Classmethod (from_matrix)

    mat_1 = ImmutableMatrix([
        [s*(s + 1)*(s - 3)/(s**4 + 1), 2],
        [p, p*(s + 1)/(s*(s**1 + 1))]
        ])
    mat_2 = ImmutableMatrix([[(2*s + 1)/(s**2 - 9)]])
    mat_3 = ImmutableMatrix([[1, 2], [3, 4]])
    assert TransferFunctionMatrix.from_Matrix(mat_1, s) == \
        TransferFunctionMatrix([[TransferFunction(s*(s - 3)*(s + 1), s**4 + 1, s), TransferFunction(2, 1, s)],
        [TransferFunction(p, 1, s), TransferFunction(p, s, s)]])
    assert TransferFunctionMatrix.from_Matrix(mat_2, s) == \
        TransferFunctionMatrix([[TransferFunction(2*s + 1, s**2 - 9, s)]])
    assert TransferFunctionMatrix.from_Matrix(mat_3, p) == \
        TransferFunctionMatrix([[TransferFunction(1, 1, p), TransferFunction(2, 1, p)],
        [TransferFunction(3, 1, p), TransferFunction(4, 1, p)]])

    assert TransferFunctionMatrix.from_Matrix(mat_1, s, 2) == \
        TransferFunctionMatrix([
            [DiscreteTransferFunction(s*(s - 3)*(s + 1), s**4 + 1, s, 2),
             DiscreteTransferFunction(2, 1, s, 2)],
            [DiscreteTransferFunction(p, 1, s, 2),
             DiscreteTransferFunction(p, s, s, 2)]])
    assert TransferFunctionMatrix.from_Matrix(mat_2, s, 0.1) == \
        TransferFunctionMatrix([[DiscreteTransferFunction(2*s + 1, s**2 - 9, s, 0.1)]])
    assert TransferFunctionMatrix.from_Matrix(mat_3, p, T) == \
        TransferFunctionMatrix([
            [DiscreteTransferFunction(1, 1, p, T), DiscreteTransferFunction(2, 1, p, T)],
            [DiscreteTransferFunction(3, 1, p, T), DiscreteTransferFunction(4, 1, p, T)]])

    #  Negating a TFM

    tfm1 = TransferFunctionMatrix([[TF1], [TF2]])
    assert -tfm1 == TransferFunctionMatrix([[-TF1], [-TF2]])

    tfm2 = TransferFunctionMatrix([[TF1, TF2, TF3], [tf5, -TF1, -TF3]])
    assert -tfm2 == TransferFunctionMatrix([[-TF1, -TF2, -TF3], [-tf5, TF1, TF3]])

    # subs()

    H_1 = TransferFunctionMatrix.from_Matrix(mat_1, s)
    H_2 = TransferFunctionMatrix([[TransferFunction(a*p*s, k*s**2, s), TransferFunction(p*s, k*(s**2 - a), s)]])
    assert H_1.subs(p, 1) == TransferFunctionMatrix([[TransferFunction(s*(s - 3)*(s + 1), s**4 + 1, s), TransferFunction(2, 1, s)], [TransferFunction(1, 1, s), TransferFunction(1, s, s)]])
    assert H_1.subs({p: 1}) == TransferFunctionMatrix([[TransferFunction(s*(s - 3)*(s + 1), s**4 + 1, s), TransferFunction(2, 1, s)], [TransferFunction(1, 1, s), TransferFunction(1, s, s)]])
    assert H_1.subs({p: 1, s: 1}) == TransferFunctionMatrix([[TransferFunction(s*(s - 3)*(s + 1), s**4 + 1, s), TransferFunction(2, 1, s)], [TransferFunction(1, 1, s), TransferFunction(1, s, s)]]) # This should ignore `s` as it is `var`
    assert H_2.subs(p, 2) == TransferFunctionMatrix([[TransferFunction(2*a*s, k*s**2, s), TransferFunction(2*s, k*(-a + s**2), s)]])
    assert H_2.subs(k, 1) == TransferFunctionMatrix([[TransferFunction(a*p*s, s**2, s), TransferFunction(p*s, -a + s**2, s)]])
    assert H_2.subs(a, 0) == TransferFunctionMatrix([[TransferFunction(0, k*s**2, s), TransferFunction(p*s, k*s**2, s)]])
    assert H_2.subs({p: 1, k: 1, a: a0}) == TransferFunctionMatrix([[TransferFunction(a0*s, s**2, s), TransferFunction(s, -a0 + s**2, s)]])

    # eval_frequency()
    assert H_2.eval_frequency(S(1)/2 + I) == Matrix([[2*a*p/(5*k) - 4*I*a*p/(5*k), I*p/(-a*k - 3*k/4 + I*k) + p/(-2*a*k - 3*k/2 + 2*I*k)]])

    # transpose()

    assert H_1.transpose() == TransferFunctionMatrix([[TransferFunction(s*(s - 3)*(s + 1), s**4 + 1, s), TransferFunction(p, 1, s)], [TransferFunction(2, 1, s), TransferFunction(p, s, s)]])
    assert H_2.transpose() == TransferFunctionMatrix([[TransferFunction(a*p*s, k*s**2, s)], [TransferFunction(p*s, k*(-a + s**2), s)]])
    assert H_1.transpose().transpose() == H_1
    assert H_2.transpose().transpose() == H_2

    # elem_poles()

    assert H_1.elem_poles() == [[[-sqrt(2)/2 - sqrt(2)*I/2, -sqrt(2)/2 + sqrt(2)*I/2, sqrt(2)/2 - sqrt(2)*I/2, sqrt(2)/2 + sqrt(2)*I/2], []],
        [[], [0]]]
    assert H_2.elem_poles() == [[[0, 0], [sqrt(a), -sqrt(a)]]]
    assert tfm2.elem_poles() == [[[wn*(-zeta + sqrt((zeta - 1)*(zeta + 1))), wn*(-zeta - sqrt((zeta - 1)*(zeta + 1)))], [], [-p/a2]],
        [[-a0], [wn*(-zeta + sqrt((zeta - 1)*(zeta + 1))), wn*(-zeta - sqrt((zeta - 1)*(zeta + 1)))], [-p/a2]]]

    # elem_zeros()

    assert H_1.elem_zeros() == [[[-1, 0, 3], []], [[], []]]
    assert H_2.elem_zeros() == [[[0], [0]]]
    assert tfm2.elem_zeros() == [[[], [], [a2*p]],
        [[-a2/(2*a1) - sqrt(4*a0*a1 + a2**2)/(2*a1), -a2/(2*a1) + sqrt(4*a0*a1 + a2**2)/(2*a1)], [], [a2*p]]]

    # doit()

    H_3 = TransferFunctionMatrix([[Series(TransferFunction(1, s**3 - 3, s), TransferFunction(s**2 - 2*s + 5, 1, s), TransferFunction(1, s, s))]])
    H_4 = TransferFunctionMatrix([[Parallel(TransferFunction(s**3 - 3, 4*s**4 - s**2 - 2*s + 5, s), TransferFunction(4 - s**3, 4*s**4 - s**2 - 2*s + 5, s))]])

    assert H_3.doit() == TransferFunctionMatrix([[TransferFunction(s**2 - 2*s + 5, s*(s**3 - 3), s)]])
    assert H_4.doit() == TransferFunctionMatrix([[TransferFunction(1, 4*s**4 - s**2 - 2*s + 5, s)]])

    dH_3 = TransferFunctionMatrix([
        [Series(DiscreteTransferFunction(1, s**3 - 3, s, 0.1),
            DiscreteTransferFunction(s**2 - 2*s + 5, 1, s, 0.1),
            DiscreteTransferFunction(1, s, s, 0.1))]])
    dH_4 = TransferFunctionMatrix([
        [Parallel(DiscreteTransferFunction(s**3 - 3, 4*s**4 - s**2 - 2*s + 5, s, 0.1),
                  DiscreteTransferFunction(4 - s**3, 4*s**4 - s**2 - 2*s + 5, s, 0.1))
                  ]])

    assert dH_3.doit() == TransferFunctionMatrix([
        [DiscreteTransferFunction(s**2 - 2*s + 5, s*(s**3 - 3), s, 0.1)]])
    assert dH_4.doit() == TransferFunctionMatrix([
        [DiscreteTransferFunction(1, 4*s**4 - s**2 - 2*s + 5, s, 0.1)]])

    # _flat()

    assert H_1._flat() == [TransferFunction(s*(s - 3)*(s + 1), s**4 + 1, s), TransferFunction(2, 1, s), TransferFunction(p, 1, s), TransferFunction(p, s, s)]
    assert H_2._flat() == [TransferFunction(a*p*s, k*s**2, s), TransferFunction(p*s, k*(-a + s**2), s)]
    assert H_3._flat() == [Series(TransferFunction(1, s**3 - 3, s), TransferFunction(s**2 - 2*s + 5, 1, s), TransferFunction(1, s, s))]
    assert H_4._flat() == [Parallel(TransferFunction(s**3 - 3, 4*s**4 - s**2 - 2*s + 5, s), TransferFunction(4 - s**3, 4*s**4 - s**2 - 2*s + 5, s))]

    # evalf()

    assert H_1.evalf() == \
        TransferFunctionMatrix(((TransferFunction(s*(s - 3.0)*(s + 1.0), s**4 + 1.0, s), TransferFunction(2.0, 1, s)), (TransferFunction(1.0*p, 1, s), TransferFunction(p, s, s))))
    assert H_2.subs({a:3.141, p:2.88, k:2}).evalf() == \
        TransferFunctionMatrix(((TransferFunction(4.5230399999999999494093572138808667659759521484375, s, s),
        TransferFunction(2.87999999999999989341858963598497211933135986328125*s, 2.0*s**2 - 6.282000000000000028421709430404007434844970703125, s)),))

    dH_1 = TransferFunctionMatrix.from_Matrix(mat_1, s, T)
    dH_2 = TransferFunctionMatrix([[DiscreteTransferFunction(a*p*s, k*s**2, s, T),
                                   DiscreteTransferFunction(p*s, k*(s**2 - a), s, T)
                                   ]])
    assert dH_1.evalf() == \
        TransferFunctionMatrix((
            (DiscreteTransferFunction(s*(s - 3.0)*(s + 1.0), s**4 + 1.0, s, T),
             DiscreteTransferFunction(2.0, 1, s, T)),
            (DiscreteTransferFunction(1.0*p, 1, s, T),
             DiscreteTransferFunction(p, s, s, T))))
    assert dH_2.subs({a:3.141, p:2.88, k:2}).evalf() == \
        TransferFunctionMatrix(((DiscreteTransferFunction(4.5230399999999999494093572138808667659759521484375, s, s, T),
        DiscreteTransferFunction(2.87999999999999989341858963598497211933135986328125*s, 2.0*s**2 - 6.282000000000000028421709430404007434844970703125, s, T)),))
    # simplify()

    H_5 = TransferFunctionMatrix([[TransferFunction(s**5 + s**3 + s, s - s**2, s),
        TransferFunction((s + 3)*(s - 1), (s - 1)*(s + 5), s)]])

    assert H_5.simplify() == simplify(H_5) == \
        TransferFunctionMatrix(((
            TransferFunction(-s**4 - s**2 - 1, s - 1, s),
            TransferFunction(s + 3, s + 5, s)),))

    dH_5 = TransferFunctionMatrix([[
        DiscreteTransferFunction(s**5 + s**3 + s, s - s**2, s, 2),
        DiscreteTransferFunction((s + 3)*(s - 1), (s - 1)*(s + 5), s, 2)]])

    assert dH_5.simplify() == simplify(dH_5) == \
        TransferFunctionMatrix(((
            DiscreteTransferFunction(-s**4 - s**2 - 1, s - 1, s, 2),
            DiscreteTransferFunction(s + 3, s + 5, s, 2)),))

    # expand()

    assert (H_1.expand()
            == TransferFunctionMatrix(((TransferFunction(s**3 - 2*s**2 - 3*s, s**4 + 1, s), TransferFunction(2, 1, s)),
            (TransferFunction(p, 1, s), TransferFunction(p, s, s)))))
    assert H_5.expand() == \
        TransferFunctionMatrix((
            (TransferFunction(s**5 + s**3 + s, -s**2 + s, s),
             TransferFunction(s**2 + 2*s - 3, s**2 + 4*s - 5, s)),))

    assert (dH_1.expand()
            == TransferFunctionMatrix((
                (DiscreteTransferFunction(s**3 - 2*s**2 - 3*s, s**4 + 1, s, T),
                 DiscreteTransferFunction(2, 1, s, T)),
                (DiscreteTransferFunction(p, 1, s, T), DiscreteTransferFunction(p, s, s, T))
                )))
    assert dH_5.expand() == \
        TransferFunctionMatrix((
            (DiscreteTransferFunction(s**5 + s**3 + s, -s**2 + s, s, 2),
             DiscreteTransferFunction(s**2 + 2*s - 3, s**2 + 4*s - 5, s, 2)),))

def test_TransferFunction_gbt():
    # simple transfer function, e.g. ohms law
    tf = TransferFunction(1, a*s+b, s)
    numZ, denZ = gbt(tf, T, 0.5)
    # discretized transfer function with coefs from tf.gbt()
    dtf_test_bilinear1 = DiscreteTransferFunction.from_coeff_lists(numZ,
                                                             denZ,
                                                             z, T).expand()
    dtf_test_bilinear2 = DiscreteTransferFunction.from_gbt(tf, T,
                                                           0.5, z).expand()
    # corresponding tf with manually calculated coefs
    num_manual = [T/(2*(a + b*T/2)), T/(2*(a + b*T/2))]
    den_manual = [1, (-a + b*T/2)/(a + b*T/2)]
    dtf_manual = DiscreteTransferFunction.\
        from_coeff_lists(num_manual, den_manual, z, T).expand()
    num_manual = [T/(2*(a + b*T/2)), T/(2*(a + b*T/2))]
    den_manual = [1, (-a + b*T/2)/(a + b*T/2)]
    dtf_manual = DiscreteTransferFunction.\
        from_coeff_lists(num_manual, den_manual, z, T).expand()

    assert dtf_test_bilinear1 == dtf_test_bilinear2 == dtf_manual
    assert dtf_test_bilinear1 == dtf_test_bilinear2 == dtf_manual

    tf = TransferFunction(1, a*s+b, s)
    numZ, denZ = gbt(tf, T, 0)
    # discretized transfer function with coefs from tf.from_gbt()
    dtf_test_forward1 = DiscreteTransferFunction.from_coeff_lists(numZ,
                                                           denZ,
                                                           z, T).expand()
    dtf_test_forward2 = DiscreteTransferFunction.from_gbt(tf, T, 0, z).expand()

    # corresponding tf with manually calculated coefs
    num_manual = [T/a]
    den_manual = [1, (-a + b*T)/a]
    dtf_manual = DiscreteTransferFunction.\
        from_coeff_lists(num_manual, den_manual, z, T).expand()

    assert dtf_test_forward1 == dtf_test_forward2 == dtf_manual

    tf = TransferFunction(1, a*s+b, s)
    numZ, denZ = gbt(tf, T, 1)
    # discretized transfer function with coefs from tf.from_gbt()
    dtf_test_backward1 = DiscreteTransferFunction.from_coeff_lists(numZ,
                                                             denZ,
                                                             z, T).expand()
    dtf_test_backward2 = DiscreteTransferFunction.from_gbt(tf, T, 1, z).expand()
    # corresponding tf with manually calculated coefs
    num_manual = [T/(a + b*T), 0]
    den_manual = [1, -a/(a + b*T)]
    dtf_manual = DiscreteTransferFunction.\
        from_coeff_lists(num_manual, den_manual, z, T).expand()

    assert dtf_test_backward1 == dtf_test_backward2 == dtf_manual

    tf = TransferFunction(1, a*s+b, s)
    numZ, denZ = gbt(tf, T, 0.3)
    # discretized transfer function with coefs from tf.from_gbt()
    dtf_test_gbt1 = DiscreteTransferFunction.from_coeff_lists(numZ,
                                                        denZ,
                                                        z, T).expand()
    dtf_test_gbt2 = DiscreteTransferFunction.from_gbt(tf, T, 0.3, z).expand()
    # corresponding tf with manually calculated coefs
    num_manual = [3*T/(10*(a + 3*b*T/10)), 7*T/(10*(a + 3*b*T/10))]
    den_manual = [1, (-a + 7*b*T/10)/(a + 3*b*T/10)]
    dtf_manual = DiscreteTransferFunction.\
        from_coeff_lists(num_manual, den_manual, z, T).expand()

    assert dtf_test_gbt1 == dtf_test_gbt2 == dtf_manual

def test_TransferFunction_bilinear():
    # simple transfer function, e.g. ohms law
    tf = TransferFunction(1, a*s+b, s)
    numZ, denZ = bilinear(tf, T)
    # discretized transfer function with coefs from tf.bilinear()
    dtf_test_bilinear1 = DiscreteTransferFunction.from_coeff_lists(numZ,
                                                             denZ,
                                                             z, T).expand()
    dtf_test_bilinear2 = \
        DiscreteTransferFunction.from_bilinear(tf, T, z).expand()
    # corresponding tf with manually calculated coefs
    num_manual = [T/(2*(a + b*T/2)), T/(2*(a + b*T/2))]
    den_manual = [1, (-a + b*T/2)/(a + b*T/2)]
    dtf_manual = DiscreteTransferFunction.\
        from_coeff_lists(num_manual, den_manual, z, T).expand()

    assert dtf_test_bilinear1 == dtf_test_bilinear2 == dtf_manual

def test_TransferFunction_forward_diff():
    # simple transfer function, e.g. ohms law
    tf = TransferFunction(1, a*s+b, s)
    numZ, denZ = forward_diff(tf, T)
    # discretized transfer function with coefs from tf.forward_diff()
    dtf_test_forward1 = DiscreteTransferFunction.from_coeff_lists(numZ,
                                                           denZ,
                                                           z, T).expand()
    dtf_test_forward2 = \
        DiscreteTransferFunction.from_forward_diff(tf, T, z).expand()
    # corresponding tf with manually calculated coefs
    num_manual = [T/a]
    den_manual = [1, (-a + b*T)/a]
    dtf_manual = DiscreteTransferFunction.\
        from_coeff_lists(num_manual, den_manual, z, T).expand()

    assert dtf_test_forward1 == dtf_test_forward2 == dtf_manual

def test_TransferFunction_backward_diff():
    # simple transfer function, e.g. ohms law
    tf = TransferFunction(1, a*s+b, s)
    numZ, denZ = backward_diff(tf, T)
    # discretized transfer function with coefs from tf.from_backward_diff()
    dtf_test_backward1 = DiscreteTransferFunction.from_coeff_lists(numZ,
                                                             denZ,
                                                             z, T).expand()
    dtf_test_backward2 = \
        DiscreteTransferFunction.from_backward_diff(tf, T, z).expand()
    # corresponding tf with manually calculated coefs
    num_manual = [T/(a + b*T), 0]
    den_manual = [1, -a/(a + b*T)]
    dtf_manual = DiscreteTransferFunction.\
        from_coeff_lists(num_manual, den_manual, z, T).expand()

    assert dtf_test_backward1 == dtf_test_backward2 == dtf_manual

def test_TransferFunction_phase_margin():
    # Test for phase margin
    tf1 = TransferFunction(10, p**3 + 1, p)
    tf2 = TransferFunction(s**2, 10, s)
    tf3 = TransferFunction(1, a*s+b, s)
    tf4 = TransferFunction((s + 1)*exp(s/tau), s**2 + 2, s)
    tf_m = TransferFunctionMatrix([[tf2],[tf3]])

    assert phase_margin(tf1) == -180 + 180*atan(3*sqrt(11))/pi
    assert phase_margin(tf2) == 0

    raises(NotImplementedError, lambda: phase_margin(tf4))
    raises(ValueError, lambda: phase_margin(tf3))
    raises(ValueError, lambda: phase_margin(MIMOSeries(tf_m)))

def test_TransferFunction_gain_margin():
    # Test for gain margin
    tf1 = TransferFunction(s**2, 5*(s+1)*(s-5)*(s-10), s)
    tf2 = TransferFunction(s**2 + 2*s + 1, 1, s)
    tf3 = TransferFunction(1, a*s+b, s)
    tf4 = TransferFunction((s + 1)*exp(s/tau), s**2 + 2, s)
    tf_m = TransferFunctionMatrix([[tf2],[tf3]])

    assert gain_margin(tf1) == -20*log(S(7)/540)/log(10)
    assert gain_margin(tf2) == oo

    raises(NotImplementedError, lambda: gain_margin(tf4))
    raises(ValueError, lambda: gain_margin(tf3))
    raises(ValueError, lambda: gain_margin(MIMOSeries(tf_m)))

def test_create_state_space():
    A = Matrix([[0, 1], [1, 0]])
    B = Matrix([1, 0])
    C = Matrix([[0, 1]])
    D = Matrix([0])

    cont_ss1 = create_state_space(A, B, C, D)
    assert isinstance(cont_ss1, StateSpace)

    cont_ss2 = create_state_space(A, B, C, D, sampling_time = 0)
    assert isinstance(cont_ss2, StateSpace)

    disc_ss1 = create_state_space(A, B, C, D, 0.1)
    assert isinstance(disc_ss1, DiscreteStateSpace)

    disc_ss2 = create_state_space(A, B, C, D, T)
    assert isinstance(disc_ss2, DiscreteStateSpace)

def test_StateSpace_construction():
    # using different numbers for a SISO system.
    A1 = Matrix([[0, 1], [1, 0]])
    B1 = Matrix([1, 0])
    C1 = Matrix([[0, 1]])
    D1 = Matrix([0])
    ss1 = StateSpace(A1, B1, C1, D1)

    assert ss1.state_matrix == Matrix([[0, 1], [1, 0]])
    assert ss1.input_matrix == Matrix([1, 0])
    assert ss1.output_matrix == Matrix([[0, 1]])
    assert ss1.feedforward_matrix == Matrix([0])
    assert ss1.sampling_time == 0
    assert ss1.args == (Matrix([[0, 1], [1, 0]]), Matrix([[1], [0]]), Matrix([[0, 1]]), Matrix([[0]]))

    # using different symbols for a SISO system.
    ss2 = StateSpace(Matrix([a0]), Matrix([a1]),
                    Matrix([a2]), Matrix([a3]))

    assert ss2.state_matrix == Matrix([[a0]])
    assert ss2.input_matrix == Matrix([[a1]])
    assert ss2.output_matrix == Matrix([[a2]])
    assert ss2.feedforward_matrix == Matrix([[a3]])
    assert ss2.args == (Matrix([[a0]]), Matrix([[a1]]), Matrix([[a2]]), Matrix([[a3]]))

    # using different numbers for a MIMO system.
    ss3 = StateSpace(Matrix([[-1.5, -2], [1, 0]]),
                    Matrix([[0.5, 0], [0, 1]]),
                    Matrix([[0, 1], [0, 2]]),
                    Matrix([[2, 2], [1, 1]]))

    assert ss3.state_matrix == Matrix([[-1.5, -2], [1,  0]])
    assert ss3.input_matrix == Matrix([[0.5, 0], [0, 1]])
    assert ss3.output_matrix == Matrix([[0, 1], [0, 2]])
    assert ss3.feedforward_matrix == Matrix([[2, 2], [1, 1]])
    assert ss3.args == (Matrix([[-1.5, -2],
                                [1,  0]]),
                        Matrix([[0.5, 0],
                                [0, 1]]),
                        Matrix([[0, 1],
                                [0, 2]]),
                        Matrix([[2, 2],
                                [1, 1]]))

    # using different symbols for a MIMO system.
    A4 = Matrix([[a0, a1], [a2, a3]])
    B4 = Matrix([[b0, b1], [b2, b3]])
    C4 = Matrix([[c0, c1], [c2, c3]])
    D4 = Matrix([[d0, d1], [d2, d3]])
    ss4 = StateSpace(A4, B4, C4, D4)

    assert ss4.state_matrix == Matrix([[a0, a1], [a2, a3]])
    assert ss4.input_matrix == Matrix([[b0, b1], [b2, b3]])
    assert ss4.output_matrix == Matrix([[c0, c1], [c2, c3]])
    assert ss4.feedforward_matrix == Matrix([[d0, d1], [d2, d3]])
    assert ss4.args == (Matrix([[a0, a1],
                                [a2, a3]]),
                        Matrix([[b0, b1],
                                [b2, b3]]),
                        Matrix([[c0, c1],
                                [c2, c3]]),
                        Matrix([[d0, d1],
                                [d2, d3]]))

    # using less matrices. Rest will be filled with a minimum of zeros.
    ss5 = StateSpace()
    assert ss5.args == (Matrix([[0]]), Matrix([[0]]), Matrix([[0]]), Matrix([[0]]))

    A6 = Matrix([[0, 1], [1, 0]])
    B6 = Matrix([1, 1])
    ss6 = StateSpace(A6, B6)

    assert ss6.state_matrix == Matrix([[0, 1], [1, 0]])
    assert ss6.input_matrix ==  Matrix([1, 1])
    assert ss6.output_matrix == Matrix([[0, 0]])
    assert ss6.feedforward_matrix == Matrix([[0]])
    assert ss6.args == (Matrix([[0, 1],
                                [1, 0]]),
                        Matrix([[1],
                                [1]]),
                        Matrix([[0, 0]]),
                        Matrix([[0]]))

    # Check if the system is SISO or MIMO.
    # If system is not SISO, then it is definitely MIMO.

    assert ss1.is_SISO == True
    assert ss2.is_SISO == True
    assert ss3.is_SISO == False
    assert ss4.is_SISO == False
    assert ss5.is_SISO == True
    assert ss6.is_SISO == True

    # ShapeError if matrices do not fit.
    raises(ShapeError, lambda: StateSpace(Matrix([s, (s+1)**2]), Matrix([s+1]),
                                          Matrix([s**2 - 1]), Matrix([2*s])))
    raises(ShapeError, lambda: StateSpace(Matrix([s]), Matrix([s+1, s**3 + 1]),
                                          Matrix([s**2 - 1]), Matrix([2*s])))
    raises(ShapeError, lambda: StateSpace(Matrix([s]), Matrix([s+1]),
                                          Matrix([[s**2 - 1], [s**2 + 2*s + 1]]), Matrix([2*s])))
    raises(ShapeError, lambda: StateSpace(Matrix([[-s, -s], [s, 0]]),
                                                Matrix([[s/2, 0], [0, s]]),
                                                Matrix([[0, s]]),
                                                Matrix([[2*s, 2*s], [s, s]])))

    # TypeError if arguments are not sympy matrices.
    raises(TypeError, lambda: StateSpace(s**2, s+1, 2*s, 1))
    raises(TypeError, lambda: StateSpace(Matrix([2, 0.5]), Matrix([-1]),
                                         Matrix([1]), 0))

def test_StateSpace_add_mul():
    A1 = Matrix([[4, 1],[2, -3]])
    B1 = Matrix([[5, 2],[-3, -3]])
    C1 = Matrix([[2, -4],[0, 1]])
    D1 = Matrix([[3, 2],[1, -1]])
    ss1 = StateSpace(A1, B1, C1, D1)

    A2 = Matrix([[-3, 4, 2],[-1, -3, 0],[2, 5, 3]])
    B2 = Matrix([[1, 4],[-3, -3],[-2, 1]])
    C2 = Matrix([[4, 2, -3],[1, 4, 3]])
    D2 = Matrix([[-2, 4],[0, 1]])
    ss2 = StateSpace(A2, B2, C2, D2)
    ss3 = StateSpace()
    ss4 = StateSpace(Matrix([1]), Matrix([2]), Matrix([3]), Matrix([4]))

    expected_add = \
        StateSpace(
        Matrix([
        [4,  1,  0,  0, 0],
        [2, -3,  0,  0, 0],
        [0,  0, -3,  4, 2],
        [0,  0, -1, -3, 0],
        [0,  0,  2,  5, 3]]),
        Matrix([
        [ 5,  2],
        [-3, -3],
        [ 1,  4],
        [-3, -3],
        [-2,  1]]),
        Matrix([
        [2, -4, 4, 2, -3],
        [0,  1, 1, 4,  3]]),
        Matrix([
        [1, 6],
        [1, 0]]))

    expected_mul = \
        StateSpace(
        Matrix([
        [ -3,   4,  2, 0,  0],
        [ -1,  -3,  0, 0,  0],
        [  2,   5,  3, 0,  0],
        [ 22,  18, -9, 4,  1],
        [-15, -18,  0, 2, -3]]),
        Matrix([
        [  1,   4],
        [ -3,  -3],
        [ -2,   1],
        [-10,  22],
        [  6, -15]]),
        Matrix([
        [14, 14, -3, 2, -4],
        [ 3, -2, -6, 0,  1]]),
        Matrix([
        [-6, 14],
        [-2,  3]]))

    assert ss1 + ss2 == expected_add
    assert ss1*ss2 == expected_mul
    assert ss3 + 1/2 == StateSpace(Matrix([[0]]), Matrix([[0]]), Matrix([[0]]), Matrix([[0.5]]))
    assert ss4*1.5 == StateSpace(Matrix([[1]]), Matrix([[2]]), Matrix([[4.5]]), Matrix([[6.0]]))
    assert 1.5*ss4 == StateSpace(Matrix([[1]]), Matrix([[3.0]]), Matrix([[3]]), Matrix([[6.0]]))
    raises(ShapeError, lambda: ss1 + ss3)
    raises(ShapeError, lambda: ss2*ss4)

def test_StateSpace_negation():
    A = Matrix([[a0, a1], [a2, a3]])
    B = Matrix([[b0, b1], [b2, b3]])
    C = Matrix([[c0, c1], [c1, c2], [c2, c3]])
    D = Matrix([[d0, d1], [d1, d2], [d2, d3]])
    SS = StateSpace(A, B, C, D)
    SS_neg = -SS

    state_mat = Matrix([[-1, 1], [1, -1]])
    input_mat = Matrix([1, -1])
    output_mat = Matrix([[-1, 1]])
    feedforward_mat = Matrix([1])
    system = StateSpace(state_mat, input_mat, output_mat, feedforward_mat)

    assert SS_neg == \
        StateSpace(Matrix([[a0, a1],
                           [a2, a3]]),
                   Matrix([[b0, b1],
                           [b2, b3]]),
                   Matrix([[-c0, -c1],
                           [-c1, -c2],
                           [-c2, -c3]]),
                   Matrix([[-d0, -d1],
                           [-d1, -d2],
                           [-d2, -d3]]))
    assert -system == \
        StateSpace(Matrix([[-1,  1],
                           [ 1, -1]]),
                   Matrix([[ 1],[-1]]),
                   Matrix([[1, -1]]),
                   Matrix([[-1]]))
    assert -SS_neg == SS
    assert -(-(-(-system))) == system

def test_SymPy_substitution_functions():
    # subs
    ss1 = StateSpace(Matrix([s]), Matrix([(s + 1)**2]), Matrix([s**2 - 1]), Matrix([2*s]))
    ss2 = StateSpace(Matrix([s + p]), Matrix([(s + 1)*(p - 1)]), Matrix([p**3 - s**3]), Matrix([s - p]))

    dss1 = DiscreteStateSpace(Matrix([s]), Matrix([(s + 1)**2]),
                              Matrix([s**2 - 1]), Matrix([2*s]), 0.1)
    dss2 = DiscreteStateSpace(Matrix([s + p]), Matrix([(s + 1)*(p - 1)]),
                              Matrix([p**3 - s**3]), Matrix([s - p]), 0.1)

    assert ss1.subs({s:5}) == StateSpace(Matrix([[5]]), Matrix([[36]]), Matrix([[24]]), Matrix([[10]]))
    assert ss2.subs({p:1}) == StateSpace(Matrix([[s + 1]]), Matrix([[0]]), Matrix([[1 - s**3]]), Matrix([[s - 1]]))

    assert dss1.subs({s:5}) == DiscreteStateSpace(Matrix([[5]]), Matrix([[36]]),
                                          Matrix([[24]]), Matrix([[10]]), 0.1)
    assert dss2.subs({p:1}) == DiscreteStateSpace(
        Matrix([[s + 1]]), Matrix([[0]]),
        Matrix([[1 - s**3]]), Matrix([[s - 1]]), 0.1)

    # xreplace
    assert ss1.xreplace({s:p}) == \
        StateSpace(Matrix([[p]]), Matrix([[(p + 1)**2]]), Matrix([[p**2 - 1]]), Matrix([[2*p]]))
    assert ss2.xreplace({s:a, p:b}) == \
        StateSpace(Matrix([[a + b]]), Matrix([[(a + 1)*(b - 1)]]), Matrix([[-a**3 + b**3]]), Matrix([[a - b]]))

    assert dss1.xreplace({s:p}) == \
        DiscreteStateSpace(Matrix([[p]]), Matrix([[(p + 1)**2]]),
                           Matrix([[p**2 - 1]]), Matrix([[2*p]]), 0.1)
    assert dss2.xreplace({s:a, p:b}) == \
        DiscreteStateSpace(Matrix([[a + b]]), Matrix([[(a + 1)*(b - 1)]]),
                           Matrix([[-a**3 + b**3]]), Matrix([[a - b]]), 0.1)

    # evalf
    p1 = a1*s + a0
    p2 = b2*s**2 + b1*s + b0
    G = StateSpace(Matrix([p1]), Matrix([p2]))

    dG = DiscreteStateSpace(Matrix([p1]), Matrix([p2]), sampling_time = 0.3)

    expect = StateSpace(Matrix([[2*s + 1]]), Matrix([[5*s**2 + 4*s + 3]]), Matrix([[0]]), Matrix([[0]]))
    expect_ = StateSpace(Matrix([[2.0*s + 1.0]]), Matrix([[5.0*s**2 + 4.0*s + 3.0]]), Matrix([[0]]), Matrix([[0]]))

    dexpect = DiscreteStateSpace(Matrix([[2*s + 1]]),
                                 Matrix([[5*s**2 + 4*s + 3]]), Matrix([[0]]),
                                 Matrix([[0]]), 0.3)
    dexpect_ = DiscreteStateSpace(Matrix([[2.0*s + 1.0]]),
                                  Matrix([[5.0*s**2 + 4.0*s + 3.0]]),
                                  Matrix([[0]]), Matrix([[0]]), 0.3)

    assert G.subs({a0: 1, a1: 2, b0: 3, b1: 4, b2: 5}) == expect
    assert G.subs({a0: 1, a1: 2, b0: 3, b1: 4, b2: 5}).evalf() == expect_
    assert expect.evalf() == expect_

    assert dG.subs({a0: 1, a1: 2, b0: 3, b1: 4, b2: 5}) == dexpect
    assert dG.subs({a0: 1, a1: 2, b0: 3, b1: 4, b2: 5}).evalf() == dexpect_
    assert dexpect.evalf() == dexpect_

def test_conversion():
    # StateSpace to TransferFunction for SISO
    A1 = Matrix([[-5, -1], [3, -1]])
    B1 = Matrix([2, 5])
    C1 = Matrix([[1, 2]])
    D1 = Matrix([0])
    H1 = StateSpace(A1, B1, C1, D1)
    H3 = StateSpace(Matrix([[a0, a1], [a2, a3]]), B = Matrix([[b1], [b2]]), C = Matrix([[c1, c2]]))
    tm1 = H1.rewrite(TransferFunction)
    raises(TypeError, lambda: H1.rewrite(DiscreteTransferFunction)) # H1 is not a discrete system
    tm2 = (-H1).rewrite(TransferFunction)

    tf1 = tm1[0][0]
    tf2 = tm2[0][0]

    assert tf1 == TransferFunction(12*s + 59, s**2 + 6*s + 8, s)
    assert tf2.num == -tf1.num
    assert tf2.den == tf1.den

    # StateSpace to TransferFunction for MIMO
    A2 = Matrix([[-1.5, -2, 3], [1, 0, 1], [2, 1, 1]])
    B2 = Matrix([[0.5, 0, 1], [0, 1, 2], [2, 2, 3]])
    C2 = Matrix([[0, 1, 0], [0, 2, 1], [1, 0, 2]])
    D2 = Matrix([[2, 2, 0], [1, 1, 1], [3, 2, 1]])
    H2 = StateSpace(A2, B2, C2, D2)
    tm3 = H2.rewrite(TransferFunction)

    # outputs for input i obtained at Index i-1. Consider input 1
    assert tm3[0][0] == TransferFunction(2.0*s**3 + 1.0*s**2 - 10.5*s + 4.5, 1.0*s**3 + 0.5*s**2 - 6.5*s - 2.5, s)
    assert tm3[0][1] == TransferFunction(2.0*s**3 + 2.0*s**2 - 10.5*s - 3.5, 1.0*s**3 + 0.5*s**2 - 6.5*s - 2.5, s)
    assert tm3[0][2] == TransferFunction(2.0*s**2 + 5.0*s - 0.5, 1.0*s**3 + 0.5*s**2 - 6.5*s - 2.5, s)
    assert H3.rewrite(TransferFunction) == [[TransferFunction(-c1*(a1*b2 - a3*b1 + b1*s) - c2*(-a0*b2 + a2*b1 + b2*s),
                                                              -a0*a3 + a0*s + a1*a2 + a3*s - s**2, s)]]
    # TransferFunction to StateSpace
    SS = TF1.rewrite(StateSpace)
    assert SS == \
        StateSpace(Matrix([[     0,          1],
                           [-wn**2, -2*wn*zeta]]),
                   Matrix([[0],
                           [1]]),
                   Matrix([[1, 0]]),
                   Matrix([[0]]))
    assert SS.rewrite(TransferFunction)[0][0] == TF1

    # TransferFunction cannot be converted to DiscreteStateSpace
    raises(TypeError, lambda: TF1.rewrite(DiscreteStateSpace))

    # Transfer function has to be proper
    raises(ValueError, lambda: TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s).rewrite(StateSpace))

    # DiscreteStateSpace to DiscreteTransferFunction for SISO
    DH1 = DiscreteStateSpace(A1, B1, C1, D1, 0.1)
    DH3 = DiscreteStateSpace(Matrix([[a0, a1], [a2, a3]]),
                             B = Matrix([[b1], [b2]]), C = Matrix([[c1, c2]]),
                             sampling_time = T)
    dtm1 = DH1.rewrite(DiscreteTransferFunction)
    raises(TypeError, lambda: DH1.rewrite(TransferFunction)) # H3 is not a continuous system
    dtm2 = (-DH1).rewrite(DiscreteTransferFunction)

    dtf1 = dtm1[0][0]
    dtf2 = dtm2[0][0]

    assert dtf1 == DiscreteTransferFunction(12*z + 59, z**2 + 6*z + 8, z, 0.1)
    assert dtf2.num == -dtf1.num
    assert dtf2.den == dtf1.den
    assert dtf1.sampling_time == 0.1
    assert dtf2.sampling_time == 0.1

    # DiscreteStateSpace to DiscreteTransferFunction for MIMO
    DH2 = DiscreteStateSpace(A2, B2, C2, D2, T)
    dtm3 = DH2.rewrite(DiscreteTransferFunction)

    # outputs for input i obtained at Index i-1. Consider input 1
    assert dtm3[0][0] == DiscreteTransferFunction(2.0*z**3 + 1.0*z**2 - 10.5*z + 4.5, 1.0*z**3 + 0.5*z**2 - 6.5*z - 2.5, z, T)
    assert dtm3[0][1] == DiscreteTransferFunction(2.0*z**3 + 2.0*z**2 - 10.5*z - 3.5, 1.0*z**3 + 0.5*z**2 - 6.5*z - 2.5, z, T)
    assert dtm3[0][2] == DiscreteTransferFunction(2.0*z**2 + 5.0*z - 0.5, 1.0*z**3 + 0.5*z**2 - 6.5*z - 2.5, z, T)
    assert DH3.rewrite(DiscreteTransferFunction) == [[DiscreteTransferFunction(-c1*(a1*b2 - a3*b1 + b1*z) - c2*(-a0*b2 + a2*b1 + b2*z),
                                                              -a0*a3 + a0*z + a1*a2 + a3*z - z**2, z, T)]]
    # DiscreteTransferFunction to DiscreteStateSpace
    dtf1 = DiscreteTransferFunction(z+5, z**3+2*z**2+4*z+3, z, 0.1)
    DSS = dtf1.rewrite(DiscreteStateSpace)
    A = Matrix([[0, 1, 0], [0, 0, 1], [-3, -4, -2]])
    B = Matrix([0, 0, 1])
    C = Matrix([[5, 1, 0]])
    D = Matrix([[0]])

    assert DSS == DiscreteStateSpace(A, B, C, D, 0.1)

    assert DSS.rewrite(DiscreteTransferFunction)[0][0] == dtf1

    # DiscreteTransferFunction cannot be converted to DiscreteStateSpace
    raises(TypeError, lambda: dtf1.rewrite(StateSpace))

    # discrete-time Transfer function has to be proper
    raises(ValueError, lambda: DiscreteTransferFunction(b2*z**2 + b1*z + b0, z + a0, z).rewrite(DiscreteStateSpace))

def test_StateSpace_dsolve():
    # https://web.mit.edu/2.14/www/Handouts/StateSpaceResponse.pdf
    # https://lpsa.swarthmore.edu/Transient/TransMethSS.html
    A1 = Matrix([[0, 1], [-2, -3]])
    B1 = Matrix([[0], [1]])
    C1 = Matrix([[1, -1]])
    D1 = Matrix([0])
    I1 = Matrix([[1], [2]])
    t = symbols('t')
    ss1 = StateSpace(A1, B1, C1, D1)

    # Zero input and Zero initial conditions
    assert ss1.dsolve() == Matrix([[0]])
    assert ss1.dsolve(initial_conditions=I1) == Matrix([[8*exp(-t) - 9*exp(-2*t)]])

    A2 = Matrix([[-2, 0], [1, -1]])
    C2 = eye(2,2)
    I2 = Matrix([2, 3])
    ss2 = StateSpace(A=A2, C=C2)
    assert ss2.dsolve(initial_conditions=I2) == Matrix([[2*exp(-2*t)], [5*exp(-t) - 2*exp(-2*t)]])

    A3 = Matrix([[-1, 1], [-4, -4]])
    B3 = Matrix([[0], [4]])
    C3 = Matrix([[0, 1]])
    D3 = Matrix([0])
    U3 = Matrix([10])
    ss3 = StateSpace(A3, B3, C3, D3)
    op = ss3.dsolve(input_vector=U3, var=t)
    assert str(op.simplify().expand().evalf()[0]) == str(5.0 + 20.7880460155075*exp(-5*t/2)*sin(sqrt(7)*t/2)
                                            - 5.0*exp(-5*t/2)*cos(sqrt(7)*t/2))

    # Test with Heaviside as input
    A4 = Matrix([[-1, 1], [-4, -4]])
    B4 = Matrix([[0], [4]])
    C4 = Matrix([[0, 1]])
    U4 = Matrix([[10*Heaviside(t)]])
    ss4 = StateSpace(A4, B4, C4)
    op4 = str(ss4.dsolve(var=t, input_vector=U4)[0].simplify().expand().evalf())
    assert op4 == str(5.0*Heaviside(t) + 20.7880460155075*exp(-5*t/2)*sin(sqrt(7)*t/2)*Heaviside(t)
                                            - 5.0*exp(-5*t/2)*cos(sqrt(7)*t/2)*Heaviside(t))

    # Test with Symbolic Matrices
    m, a, x0 = symbols('m a x_0')
    A5 = Matrix([[0, 1], [0, 0]])
    B5 = Matrix([[0], [1 / m]])
    C5 = Matrix([[1, 0]])
    I5 = Matrix([[x0], [0]])
    U5 = Matrix([[exp(-a * t)]])
    ss5 = StateSpace(A5, B5, C5)
    op5 = ss5.dsolve(initial_conditions=I5, input_vector=U5, var=t).simplify()
    assert op5[0].args[0][0] == x0 + t/(a*m) - 1/(a**2*m) + exp(-a*t)/(a**2*m)
    a11, a12, a21, a22, b1, b2, c1, c2, i1, i2 = symbols('a_11 a_12 a_21 a_22 b_1 b_2 c_1 c_2 i_1 i_2')
    A6 = Matrix([[a11, a12], [a21, a22]])
    B6 = Matrix([b1, b2])
    C6 = Matrix([[c1, c2]])
    I6 = Matrix([i1, i2])
    ss6 = StateSpace(A6, B6, C6)
    expr6 = ss6.dsolve(initial_conditions=I6)[0]
    expr6 = expr6.subs([(a11, 0), (a12, 1), (a21, -2), (a22, -3), (b1, 0), (b2, 1), (c1, 1), (c2, -1), (i1, 1), (i2, 2)])
    assert expr6 == 8*exp(-t) - 9*exp(-2*t)

def test_StateSpace_functions():
    # https://in.mathworks.com/help/control/ref/statespacemodel.obsv.html

    A_mat = Matrix([[-1.5, -2], [1, 0]])
    B_mat = Matrix([0.5, 0])
    C_mat = Matrix([[0, 1]])
    D_mat = Matrix([1])
    SS1 = StateSpace(A_mat, B_mat, C_mat, D_mat)
    SS2 = StateSpace(Matrix([[1, 1], [4, -2]]),Matrix([[0, 1], [0, 2]]),
                     Matrix([[-1, 1], [1, -1]]))
    SS3 = StateSpace(Matrix([[1, 1], [4, -2]]),Matrix([[1, -1], [1, -1]]))
    SS4 = StateSpace(Matrix([[a0, a1], [a2, a3]]), Matrix([[b1], [b2]]),
                     Matrix([[c1, c2]]))

    # Observability
    assert SS1.is_observable() == True
    assert SS2.is_observable() == False
    assert SS1.observability_matrix() == Matrix([[0, 1], [1, 0]])
    assert SS2.observability_matrix() == Matrix([[-1,  1], [ 1, -1], [ 3, -3], [-3,  3]])
    assert SS1.observable_subspace() == [Matrix([[1], [0]]), Matrix([[0], [1]])]
    assert SS2.observable_subspace() == [Matrix([[-1], [ 1]])]
    raises(NotImplementedError, lambda: SS4.observable_subspace())
    assert SS1.unobservable_subspace() == []
    assert SS2.unobservable_subspace() == [Matrix([[1],[1]])]
    raises(NotImplementedError, lambda: SS4.unobservable_subspace())

    Qo = SS4.observability_matrix().subs([(a0, 0), (a1, -6), (a2, 1), (a3, -5), (c1, 0), (c2, 1)])
    assert Qo == Matrix([[0, 1], [1, -5]])

    ss_obs = StateSpace(Matrix([[1, 0, 1], [0,0,0],[0,0,-2]]), Matrix([1,1,0]),
                        Matrix([1,1, Rational(1,3)]).T).to_observable_form()
    A_obs = ss_obs.A[:2, :2]
    A_nobs = ss_obs.A[2:, 2:]
    assert A_obs.eigenvals() == {0: 1, 1: 1}
    assert A_nobs.eigenvals() == {-2: 1}

    # Controllability
    assert SS1.is_controllable() == True
    assert SS3.is_controllable() == False
    assert SS1.controllability_matrix() ==  Matrix([[0.5, -0.75], [  0,   0.5]])
    assert SS3.controllability_matrix() == Matrix([[1, -1, 2, -2], [1, -1, 2, -2]])
    assert SS4.controllability_matrix() == \
        Matrix([[b1, a0*b1 + a1*b2], [b2, a2*b1 + a3*b2]])
    assert SS1.controllable_subspace() == [Matrix([[0.5], [  0]]), Matrix([[-0.75], [  0.5]])]
    assert SS3.controllable_subspace() == [Matrix([[1], [1]])]
    assert SS1.uncontrollable_subspace() == []
    assert SS3.uncontrollable_subspace() == [Matrix([[-1], [1]])]
    raises(NotImplementedError, lambda: SS4.uncontrollable_subspace()) # uncontrollable subspace fo symbols not implemented

    Qc = SS4.controllability_matrix().subs([(a0, 0), (a1, 1), (a2, -6), (a3, -5), (b1, 0), (b2, 1)])
    assert Qc == Matrix([[0, 1], [1, -5]])
    ss_contr = StateSpace(Matrix([[1, 0, 1], [0,0,0],[0,0,-2]]), Matrix([1,1,0]), Matrix([1,1,0]).T).to_controllable_form()
    A_contr = ss_contr.A[:2, :2]
    A_ncontr = ss_contr.A[2:, 2:]
    assert A_contr.eigenvals() == {0: 1, 1: 1}
    assert A_ncontr.eigenvals() == {-2: 1}

    # Append
    A1 = Matrix([[0, 1], [1, 0]])
    B1 = Matrix([[0], [1]])
    C1 = Matrix([[0, 1]])
    D1 = Matrix([[0]])
    ss1 = StateSpace(A1, B1, C1, D1)
    ss2 = StateSpace(Matrix([[1, 0], [0, 1]]), Matrix([[1], [0]]), Matrix([[1, 0]]), Matrix([[1]]))
    ss3 = ss1.append(ss2)
    ss4 = SS4.append(ss1)

    assert ss3.num_states == ss1.num_states + ss2.num_states
    assert ss3.num_inputs == ss1.num_inputs + ss2.num_inputs
    assert ss3.num_outputs == ss1.num_outputs + ss2.num_outputs
    assert ss3.state_matrix == Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    assert ss3.input_matrix == Matrix([[0, 0], [1, 0], [0, 1], [0, 0]])
    assert ss3.output_matrix == Matrix([[0, 1, 0, 0], [0, 0, 1, 0]])
    assert ss3.feedforward_matrix == Matrix([[0, 0], [0, 1]])

    # Using symbolic matrices
    assert ss4.num_states == SS4.num_states + ss1.num_states
    assert ss4.num_inputs == SS4.num_inputs + ss1.num_inputs
    assert ss4.num_outputs == SS4.num_outputs + ss1.num_outputs
    assert ss4.state_matrix == Matrix([[a0, a1, 0, 0], [a2, a3, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    assert ss4.input_matrix == Matrix([[b1, 0], [b2, 0], [0, 0], [0, 1]])
    assert ss4.output_matrix == Matrix([[c1, c2, 0, 0], [0, 0, 0, 1]])
    assert ss4.feedforward_matrix == Matrix([[0, 0], [0, 0]])

def test_StateSpace_series():
    # For SISO Systems
    a1 = Matrix([[0, 1], [1, 0]])
    b1 = Matrix([[0], [1]])
    c1 = Matrix([[0, 1]])
    d1 = Matrix([[0]])
    a2 = Matrix([[1, 0], [0, 1]])
    b2 = Matrix([[1], [0]])
    c2 = Matrix([[1, 0]])
    d2 = Matrix([[1]])

    ss1 = StateSpace(a1, b1, c1, d1)
    ss2 = StateSpace(a2, b2, c2, d2)
    tf1 = TransferFunction(s, s+1, s)
    ser1 = Series(ss1, ss2)
    assert ser1 == Series(StateSpace(Matrix([
                            [0, 1],
                            [1, 0]]), Matrix([
                            [0],
                            [1]]), Matrix([[0, 1]]), Matrix([[0]])), StateSpace(Matrix([
                            [1, 0],
                            [0, 1]]), Matrix([
                            [1],
                            [0]]), Matrix([[1, 0]]), Matrix([[1]])))
    assert ser1.doit() == StateSpace(
                            Matrix([
                            [0, 1, 0, 0],
                            [1, 0, 0, 0],
                            [0, 1, 1, 0],
                            [0, 0, 0, 1]]),
                            Matrix([
                            [0],
                            [1],
                            [0],
                            [0]]),
                            Matrix([[0, 1, 1, 0]]),
                            Matrix([[0]]))

    assert ser1.num_inputs == 1
    assert ser1.num_outputs == 1
    assert ser1.rewrite(TransferFunction) == TransferFunction(s**2, s**3 - s**2 - s + 1, s)
    raises(TypeError, lambda: ser1.rewrite(DiscreteTransferFunction))  # cannot convert to discrete transfer function

    ser2 = Series(ss1)
    ser3 = Series(ser2, ss2)
    assert ser3.doit() == ser1.doit()

    # test issue #28326
    ss_issue = StateSpace(a2, b1+b2, c1+c2, d1)
    assert Series(ss_issue, ss_issue, ss_issue).doit() == StateSpace(
        Matrix([
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [1, 1, 1, 0, 0, 0],
        [1, 1, 0, 1, 0, 0],
        [0, 0, 1, 1, 1, 0],
        [0, 0, 1, 1, 0, 1]]),
        Matrix([
        [1],
        [1],
        [0],
        [0],
        [0],
        [0]]),
        Matrix([[0, 0, 0, 0, 1, 1]]),
        Matrix([[0]]))

    # TransferFunction interconnection with StateSpace
    ser_tf = Series(tf1, ss1)
    assert ser_tf == Series(TransferFunction(s, s + 1, s), StateSpace(Matrix([
                            [0, 1],
                            [1, 0]]), Matrix([
                            [0],
                            [1]]), Matrix([[0, 1]]), Matrix([[0]])))
    assert ser_tf.doit() == StateSpace(
                            Matrix([
                            [-1, 0,  0],
                            [0, 0,  1],
                            [-1, 1, 0]]),
                            Matrix([
                            [1],
                            [0],
                            [1]]),
                            Matrix([[0, 0, 1]]),
                            Matrix([[0]]))
    assert ser_tf.rewrite(TransferFunction) == TransferFunction(s**2, s**3 + s**2 - s - 1, s)

    # For MIMO Systems
    a3 = Matrix([[4, 1], [2, -3]])
    b3 = Matrix([[5, 2], [-3, -3]])
    c3 = Matrix([[2, -4], [0, 1]])
    d3 = Matrix([[3, 2], [1, -1]])
    a4 = Matrix([[-3, 4, 2], [-1, -3, 0], [2, 5, 3]])
    b4 = Matrix([[1, 4], [-3, -3], [-2, 1]])
    c4 = Matrix([[4, 2, -3], [1, 4, 3]])
    d4 = Matrix([[-2, 4], [0, 1]])
    ss3 = StateSpace(a3, b3, c3, d3)
    ss4 = StateSpace(a4, b4, c4, d4)
    ser4 = MIMOSeries(ss3, ss4)
    assert ser4 == MIMOSeries(StateSpace(Matrix([
                    [4,  1],
                    [2, -3]]), Matrix([
                    [ 5,  2],
                    [-3, -3]]), Matrix([
                    [2, -4],
                    [0,  1]]), Matrix([
                    [3,  2],
                    [1, -1]])), StateSpace(Matrix([
                    [-3,  4, 2],
                    [-1, -3, 0],
                    [ 2,  5, 3]]), Matrix([
                    [ 1,  4],
                    [-3, -3],
                    [-2,  1]]), Matrix([
                    [4, 2, -3],
                    [1, 4,  3]]), Matrix([
                    [-2, 4],
                    [ 0, 1]])))
    assert ser4.doit() == StateSpace(
                        Matrix([
                        [4,   1,  0, 0,  0],
                        [2,  -3,  0, 0,  0],
                        [2,   0,  -3, 4,  2],
                        [-6,  9, -1, -3,  0],
                        [-4, 9,  2, 5, 3]]),
                        Matrix([
                        [5,   2],
                        [-3,  -3],
                        [7,   -2],
                        [-12,  -3],
                        [-5, -5]]),
                        Matrix([
                        [-4, 12, 4, 2, -3],
                        [0, 1, 1, 4, 3]]),
                        Matrix([
                        [-2, -8],
                        [1, -1]]))
    assert ser4.num_inputs == ss3.num_inputs
    assert ser4.num_outputs == ss4.num_outputs
    ser5 = MIMOSeries(ss3)
    ser6 = MIMOSeries(ser5, ss4)
    assert ser6.doit() == ser4.doit()
    assert ser6.rewrite(TransferFunctionMatrix) == ser4.rewrite(TransferFunctionMatrix)
    tf2 = TransferFunction(1, s, s)
    tf3 = TransferFunction(1, s+1, s)
    tf4 = TransferFunction(s, s+2, s)
    tfm = TransferFunctionMatrix([[tf1, tf2], [tf3, tf4]])
    ser6 = MIMOSeries(ss3, tfm)
    assert ser6 == MIMOSeries(StateSpace(Matrix([
                        [4,  1],
                        [2, -3]]), Matrix([
                        [ 5,  2],
                        [-3, -3]]), Matrix([
                        [2, -4],
                        [0,  1]]), Matrix([
                        [3,  2],
                        [1, -1]])), TransferFunctionMatrix((
                        (TransferFunction(s, s + 1, s), TransferFunction(1, s, s)),
                        (TransferFunction(1, s + 1, s), TransferFunction(s, s + 2, s)))))

def test_StateSpace_parallel():
    # For SISO system
    a1 = Matrix([[0, 1], [1, 0]])
    b1 = Matrix([[0], [1]])
    c1 = Matrix([[0, 1]])
    d1 = Matrix([[0]])
    a2 = Matrix([[1, 0], [0, 1]])
    b2 = Matrix([[1], [0]])
    c2 = Matrix([[1, 0]])
    d2 = Matrix([[1]])
    ss1 = StateSpace(a1, b1, c1, d1)
    ss2 = StateSpace(a2, b2, c2, d2)
    p1 = Parallel(ss1, ss2)
    assert p1 == Parallel(StateSpace(Matrix([[0, 1], [1, 0]]), Matrix([[0], [1]]), Matrix([[0, 1]]), Matrix([[0]])),
                          StateSpace(Matrix([[1, 0],[0, 1]]), Matrix([[1],[0]]), Matrix([[1, 0]]), Matrix([[1]])))
    assert p1.doit() == StateSpace(Matrix([
                        [0, 1, 0, 0],
                        [1, 0, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]]),
                        Matrix([
                        [0],
                        [1],
                        [1],
                        [0]]),
                        Matrix([[0, 1, 1, 0]]),
                        Matrix([[1]]))
    assert p1.rewrite(TransferFunction) == TransferFunction(s*(s + 2), s**2 - 1, s)
    raises(TypeError, lambda: p1.rewrite(DiscreteTransferFunction))  # cannot convert to discrete transfer function

    # Connecting StateSpace with TransferFunction
    tf1 = TransferFunction(s, s+1, s)
    p2 = Parallel(ss1, tf1)
    assert p2 == Parallel(StateSpace(Matrix([
                        [0, 1],
                        [1, 0]]), Matrix([
                        [0],
                        [1]]), Matrix([[0, 1]]), Matrix([[0]])), TransferFunction(s, s + 1, s))
    assert p2.doit() == StateSpace(
                        Matrix([
                        [0, 1,  0],
                        [1, 0,  0],
                        [0, 0, -1]]),
                        Matrix([
                        [0],
                        [1],
                        [1]]),
                        Matrix([[0, 1, -1]]),
                        Matrix([[1]]))
    assert p2.rewrite(TransferFunction) == TransferFunction(s**2, s**2 - 1, s)

    # For MIMO
    a3 = Matrix([[4, 1], [2, -3]])
    b3 = Matrix([[5, 2], [-3, -3]])
    c3 = Matrix([[2, -4], [0, 1]])
    d3 = Matrix([[3, 2], [1, -1]])
    a4 = Matrix([[-3, 4, 2], [-1, -3, 0], [2, 5, 3]])
    b4 = Matrix([[1, 4], [-3, -3], [-2, 1]])
    c4 = Matrix([[4, 2, -3], [1, 4, 3]])
    d4 = Matrix([[-2, 4], [0, 1]])
    ss3 = StateSpace(a3, b3, c3, d3)
    ss4 = StateSpace(a4, b4, c4, d4)
    p3 = MIMOParallel(ss3, ss4)
    assert p3 == MIMOParallel(StateSpace(Matrix([
                        [4,  1],
                        [2, -3]]), Matrix([
                        [ 5,  2],
                        [-3, -3]]), Matrix([
                        [2, -4],
                        [0,  1]]), Matrix([
                        [3,  2],
                        [1, -1]])), StateSpace(Matrix([
                        [-3,  4, 2],
                        [-1, -3, 0],
                        [ 2,  5, 3]]), Matrix([
                        [ 1,  4],
                        [-3, -3],
                        [-2,  1]]), Matrix([
                        [4, 2, -3],
                        [1, 4,  3]]), Matrix([
                        [-2, 4],
                        [ 0, 1]])))
    assert p3.doit() == StateSpace(Matrix([
                        [4, 1, 0, 0, 0],
                        [2, -3, 0, 0, 0],
                        [0, 0, -3, 4, 2],
                        [0, 0, -1, -3, 0],
                        [0, 0, 2, 5, 3]]),
                        Matrix([
                        [5, 2],
                        [-3, -3],
                        [1, 4],
                        [-3, -3],
                        [-2, 1]]),
                        Matrix([
                        [2, -4, 4, 2, -3],
                        [0, 1, 1, 4, 3]]),
                        Matrix([
                        [1, 6],
                        [1, 0]]))

    # Using StateSpace with MIMOParallel.
    tf2 = TransferFunction(1, s, s)
    tf3 = TransferFunction(1, s + 1, s)
    tf4 = TransferFunction(s, s + 2, s)
    tfm = TransferFunctionMatrix([[tf1, tf2], [tf3, tf4]])
    p4 = MIMOParallel(tfm, ss3)
    assert p4 == MIMOParallel(TransferFunctionMatrix((
                        (TransferFunction(s, s + 1, s), TransferFunction(1, s, s)),
                        (TransferFunction(1, s + 1, s), TransferFunction(s, s + 2, s)))),
                        StateSpace(Matrix([
                        [4, 1],
                        [2, -3]]), Matrix([
                        [5, 2],
                        [-3, -3]]), Matrix([
                        [2, -4],
                        [0, 1]]), Matrix([
                        [3, 2],
                        [1, -1]])))

def test_StateSpace_feedback():
    # For SISO
    a1 = Matrix([[0, 1], [1, 0]])
    b1 = Matrix([[0], [1]])
    c1 = Matrix([[0, 1]])
    d1 = Matrix([[0]])
    a2 = Matrix([[1, 0], [0, 1]])
    b2 = Matrix([[1], [0]])
    c2 = Matrix([[1, 0]])
    d2 = Matrix([[1]])
    ss1 = StateSpace(a1, b1, c1, d1)
    ss2 = StateSpace(a2, b2, c2, d2)
    fd1 = Feedback(ss1, ss2)

    # Negative feedback
    assert fd1 == Feedback(StateSpace(Matrix([[0, 1], [1, 0]]), Matrix([[0], [1]]), Matrix([[0, 1]]), Matrix([[0]])),
                          StateSpace(Matrix([[1, 0],[0, 1]]), Matrix([[1],[0]]), Matrix([[1, 0]]), Matrix([[1]])), -1)
    assert fd1.doit() == StateSpace(Matrix([
                            [0,  1,  0, 0],
                            [1, -1, -1, 0],
                            [0,  1,  1, 0],
                            [0,  0,  0, 1]]), Matrix([
                            [0],
                            [1],
                            [0],
                            [0]]), Matrix(
                            [[0, 1, 0, 0]]), Matrix(
                            [[0]]))
    assert fd1.rewrite(TransferFunction) == TransferFunction(s*(s - 1), s**3 - s + 1, s)
    raises(TypeError, lambda: fd1.rewrite(DiscreteTransferFunction))  # cannot convert to discrete transfer function

    # Positive Feedback
    fd2 = Feedback(ss1, ss2, 1)
    assert fd2.doit() == StateSpace(Matrix([
                            [0, 1, 0, 0],
                            [1, 1, 1, 0],
                            [0, 1, 1, 0],
                            [0, 0, 0, 1]]), Matrix([
                            [0],
                            [1],
                            [0],
                            [0]]), Matrix(
                            [[0, 1, 0, 0]]), Matrix(
                            [[0]]))
    assert fd2.rewrite(TransferFunction) == TransferFunction(s*(s - 1), s**3 - 2*s**2 - s + 1, s)

    # Connection with TransferFunction
    tf1 = TransferFunction(s, s+1, s)
    fd3 = Feedback(ss1, tf1)
    assert fd3 == Feedback(StateSpace(Matrix([
                            [0, 1],
                            [1, 0]]), Matrix([
                            [0],
                            [1]]), Matrix([[0, 1]]), Matrix([[0]])),
                            TransferFunction(s, s + 1, s), -1)
    assert fd3.doit() == StateSpace (Matrix([
                            [0,  1,  0],
                            [1, -1,  1],
                            [0,  1, -1]]), Matrix([
                            [0],
                            [1],
                            [0]]), Matrix(
                            [[0, 1, 0]]), Matrix(
                            [[0]]))

    # For MIMO
    a3 = Matrix([[4, 1], [2, -3]])
    b3 = Matrix([[5, 2], [-3, -3]])
    c3 = Matrix([[2, -4], [0, 1]])
    d3 = Matrix([[3, 2], [1, -1]])
    a4 = Matrix([[-3, 4, 2], [-1, -3, 0], [2, 5, 3]])
    b4 = Matrix([[1, 4], [-3, -3], [-2, 1]])
    c4 = Matrix([[4, 2, -3], [1, 4, 3]])
    d4 = Matrix([[-2, 4], [0, 1]])
    ss3 = StateSpace(a3, b3, c3, d3)
    ss4 = StateSpace(a4, b4, c4, d4)

    # Negative Feedback
    fd4 = MIMOFeedback(ss3, ss4)
    assert fd4 == MIMOFeedback(StateSpace(Matrix([
                            [4,  1],
                            [2, -3]]), Matrix([
                            [ 5,  2],
                            [-3, -3]]), Matrix([
                            [2, -4],
                            [0,  1]]), Matrix([
                            [3,  2],
                            [1, -1]])), StateSpace(Matrix([
                            [-3,  4, 2],
                            [-1, -3, 0],
                            [ 2,  5, 3]]), Matrix([
                            [ 1,  4],
                            [-3, -3],
                            [-2,  1]]), Matrix([
                            [4, 2, -3],
                            [1, 4,  3]]), Matrix([
                            [-2, 4],
                            [ 0, 1]])), -1)
    assert fd4.doit() == StateSpace(Matrix([
                            [Rational(3), Rational(-3, 4), Rational(-15, 4), Rational(-37, 2), Rational(-15)],
                            [Rational(7, 2), Rational(-39, 8), Rational(9, 8), Rational(39, 4), Rational(9)],
                            [Rational(3), Rational(-41, 4), Rational(-45, 4), Rational(-51, 2), Rational(-19)],
                            [Rational(-9, 2), Rational(129, 8), Rational(73, 8), Rational(171, 4), Rational(36)],
                            [Rational(-3, 2), Rational(47, 8), Rational(31, 8), Rational(85, 4), Rational(18)]]), Matrix([
                            [Rational(-1, 4), Rational(19, 4)],
                            [Rational(3, 8), Rational(-21, 8)],
                            [Rational(1, 4), Rational(29, 4)],
                            [Rational(3, 8), Rational(-93, 8)],
                            [Rational(5, 8), Rational(-35, 8)]]), Matrix([
                            [Rational(1), Rational(-15, 4), Rational(-7, 4), Rational(-21, 2), Rational(-9)],
                            [Rational(1, 2), Rational(-13, 8), Rational(-13, 8), Rational(-19, 4), Rational(-3)]]), Matrix([
                            [Rational(-1, 4), Rational(11, 4)],
                            [Rational(1, 8), Rational(9, 8)]]))

    # Positive Feedback
    fd5 = MIMOFeedback(ss3, ss4, 1)
    assert fd5.doit() == StateSpace(
        Matrix([
            [Rational(4, 7), Rational(62, 7), Rational(1), Rational(-8),
             Rational(-69, 7)],
            [Rational(32, 7), Rational(-135, 14),Rational(-3, 2), Rational(3),
             Rational(36, 7)],
            [Rational(-10, 7), Rational(41, 7), Rational(-4), Rational(-12),
             Rational(-97, 7)],
            [Rational(12, 7), Rational(-111, 14), Rational(-5, 2), Rational(18),
             Rational(171, 7)],
            [Rational(2, 7), Rational(-29, 14), Rational(-1, 2), Rational(10),
             Rational(81, 7)]]),
        Matrix([[Rational(6, 7), Rational(-17, 7)],
                [Rational(-9, 14), Rational(15, 14)],
                [Rational(6, 7), Rational(-31, 7)],
                [Rational(-27, 14), Rational(87, 14)],
                [Rational(-15, 14), Rational(25, 14)]]),
        Matrix([[Rational(-2, 7), Rational(11, 7),Rational(1), Rational(-4),
                 Rational(-39, 7)],
                [Rational(-2, 7), Rational(15, 14), Rational(-1, 2),
                 Rational(-3), Rational(-18, 7)]]),
        Matrix([[Rational(4, 7), Rational(-9, 7)],
                [Rational(1, 14), Rational(-11, 14)]]))

def test_StateSpace_stability():
    k = symbols('k')
    B = Matrix([1, 0, 0])
    C = Matrix([[0, 1, 0]])
    D = Matrix([0])

    A1 = Matrix([[0,1,0],[0,0,1], [k-1, -2*k, -1]])
    ss1 = StateSpace(A1, B, C, D)
    ineq = ss1.get_asymptotic_stability_conditions()
    assert ineq == [True, 3 * k - 1 > 0, 1 - k > 0]

    A2 = Matrix([[1,0,0], [0,-1,k], [0,0,-1]])
    ss2 = StateSpace(A2, B, C, D)
    ineq = ss2.get_asymptotic_stability_conditions()
    assert ineq == [False]

def test_DiscreteStateSpace_construction():
    # using different numbers for a SISO system.
    A1 = Matrix([[0, 1], [1, 0]])
    B1 = Matrix([1, 0])
    C1 = Matrix([[0, 1]])
    D1 = Matrix([0])
    ss1 = DiscreteStateSpace(A1, B1, C1, D1)

    assert ss1.state_matrix == Matrix([[0, 1], [1, 0]])
    assert ss1.input_matrix == Matrix([1, 0])
    assert ss1.output_matrix == Matrix([[0, 1]])
    assert ss1.feedforward_matrix == Matrix([0])
    assert ss1.sampling_time == 1
    assert ss1.args == (Matrix([[0, 1], [1, 0]]), Matrix([[1], [0]]),
                        Matrix([[0, 1]]), Matrix([[0]]), 1)

    # using different symbols for a SISO system.
    ss2 = DiscreteStateSpace(Matrix([a0]), Matrix([a1]),
                    Matrix([a2]), Matrix([a3]), 0.1)

    assert ss2.state_matrix == Matrix([[a0]])
    assert ss2.input_matrix == Matrix([[a1]])
    assert ss2.output_matrix == Matrix([[a2]])
    assert ss2.feedforward_matrix == Matrix([[a3]])
    assert ss2.sampling_time == 0.1
    assert isinstance(ss2.sampling_time, Number) is True # ensure it is a sympy object, not just a python float
    assert ss2.args == (Matrix([[a0]]), Matrix([[a1]]), Matrix([[a2]]),
                        Matrix([[a3]]), 0.1)

    # using different numbers for a MIMO system.
    ss3 = DiscreteStateSpace(Matrix([[-1.5, -2], [1, 0]]),
                    Matrix([[0.5, 0], [0, 1]]),
                    Matrix([[0, 1], [0, 2]]),
                    Matrix([[2, 2], [1, 1]]), T)

    assert ss3.state_matrix == Matrix([[-1.5, -2], [1,  0]])
    assert ss3.input_matrix == Matrix([[0.5, 0], [0, 1]])
    assert ss3.output_matrix == Matrix([[0, 1], [0, 2]])
    assert ss3.feedforward_matrix == Matrix([[2, 2], [1, 1]])
    assert ss3.sampling_time == T
    assert ss3.args == (Matrix([[-1.5, -2],
                                [1,  0]]),
                        Matrix([[0.5, 0],
                                [0, 1]]),
                        Matrix([[0, 1],
                                [0, 2]]),
                        Matrix([[2, 2],
                                [1, 1]]), T)

    # using different symbols for a MIMO system.
    A4 = Matrix([[a0, a1], [a2, a3]])
    B4 = Matrix([[b0, b1], [b2, b3]])
    C4 = Matrix([[c0, c1], [c2, c3]])
    D4 = Matrix([[d0, d1], [d2, d3]])
    ss4 = DiscreteStateSpace(A4, B4, C4, D4, 12)

    assert ss4.state_matrix == Matrix([[a0, a1], [a2, a3]])
    assert ss4.input_matrix == Matrix([[b0, b1], [b2, b3]])
    assert ss4.output_matrix == Matrix([[c0, c1], [c2, c3]])
    assert ss4.feedforward_matrix == Matrix([[d0, d1], [d2, d3]])
    assert ss4.sampling_time == 12
    assert ss4.args == (Matrix([[a0, a1],
                                [a2, a3]]),
                        Matrix([[b0, b1],
                                [b2, b3]]),
                        Matrix([[c0, c1],
                                [c2, c3]]),
                        Matrix([[d0, d1],
                                [d2, d3]]), 12)

    # using less matrices. Rest will be filled with a minimum of zeros.
    ss5 = DiscreteStateSpace()
    assert ss5.args == (Matrix([[0]]), Matrix([[0]]), Matrix([[0]]),
                        Matrix([[0]]), 1)

    A6 = Matrix([[0, 1], [1, 0]])
    B6 = Matrix([1, 1])
    ss6 = DiscreteStateSpace(A6, B6, sampling_time = 2)

    assert ss6.state_matrix == Matrix([[0, 1], [1, 0]])
    assert ss6.input_matrix ==  Matrix([1, 1])
    assert ss6.output_matrix == Matrix([[0, 0]])
    assert ss6.feedforward_matrix == Matrix([[0]])
    assert ss6.sampling_time == 2
    assert ss6.args == (Matrix([[0, 1],
                                [1, 0]]),
                        Matrix([[1],
                                [1]]),
                        Matrix([[0, 0]]),
                        Matrix([[0]]), 2)

    # Check if the system is SISO or MIMO.
    # If system is not SISO, then it is definitely MIMO.

    assert ss1.is_SISO == True
    assert ss2.is_SISO == True
    assert ss3.is_SISO == False
    assert ss4.is_SISO == False
    assert ss5.is_SISO == True
    assert ss6.is_SISO == True

    # ShapeError if matrices do not fit.
    raises(ShapeError, lambda: StateSpace(Matrix([s, (s+1)**2]), Matrix([s+1]),
                                          Matrix([s**2 - 1]), Matrix([2*s])))
    raises(ShapeError, lambda: StateSpace(Matrix([s]), Matrix([s+1, s**3 + 1]),
                                          Matrix([s**2 - 1]), Matrix([2*s])))
    raises(ShapeError, lambda: StateSpace(Matrix([s]), Matrix([s+1]),
                                          Matrix([[s**2 - 1], [s**2 + 2*s + 1]]), Matrix([2*s])))
    raises(ShapeError, lambda: StateSpace(Matrix([[-s, -s], [s, 0]]),
                                                Matrix([[s/2, 0], [0, s]]),
                                                Matrix([[0, s]]),
                                                Matrix([[2*s, 2*s], [s, s]])))

    # TypeError if arguments are not sympy matrices.
    raises(TypeError, lambda: StateSpace(s**2, s+1, 2*s, 1))
    raises(TypeError, lambda: StateSpace(Matrix([2, 0.5]), Matrix([-1]),
                                         Matrix([1]), 0))

    raises(ValueError, lambda: DiscreteStateSpace(A1, B1, C1, D1, 0)) # sampling time cannot be zero

def test_DiscreteStateSpace_add_mul():
    A1 = Matrix([[4, 1],[2, -3]])
    B1 = Matrix([[5, 2],[-3, -3]])
    C1 = Matrix([[2, -4],[0, 1]])
    D1 = Matrix([[3, 2],[1, -1]])
    ss1 = DiscreteStateSpace(A1, B1, C1, D1, T)

    A2 = Matrix([[-3, 4, 2],[-1, -3, 0],[2, 5, 3]])
    B2 = Matrix([[1, 4],[-3, -3],[-2, 1]])
    C2 = Matrix([[4, 2, -3],[1, 4, 3]])
    D2 = Matrix([[-2, 4],[0, 1]])
    ss2 = DiscreteStateSpace(A2, B2, C2, D2, T)
    ss3 = DiscreteStateSpace(sampling_time = T)
    ss4 = DiscreteStateSpace(Matrix([1]), Matrix([2]),
                             Matrix([3]), Matrix([4]), T)

    expected_add = \
        DiscreteStateSpace(
        Matrix([
        [4,  1,  0,  0, 0],
        [2, -3,  0,  0, 0],
        [0,  0, -3,  4, 2],
        [0,  0, -1, -3, 0],
        [0,  0,  2,  5, 3]]),
        Matrix([
        [ 5,  2],
        [-3, -3],
        [ 1,  4],
        [-3, -3],
        [-2,  1]]),
        Matrix([
        [2, -4, 4, 2, -3],
        [0,  1, 1, 4,  3]]),
        Matrix([
        [1, 6],
        [1, 0]]), T)

    expected_mul = \
        DiscreteStateSpace(
        Matrix([
        [ -3,   4,  2, 0,  0],
        [ -1,  -3,  0, 0,  0],
        [  2,   5,  3, 0,  0],
        [ 22,  18, -9, 4,  1],
        [-15, -18,  0, 2, -3]]),
        Matrix([
        [  1,   4],
        [ -3,  -3],
        [ -2,   1],
        [-10,  22],
        [  6, -15]]),
        Matrix([
        [14, 14, -3, 2, -4],
        [ 3, -2, -6, 0,  1]]),
        Matrix([
        [-6, 14],
        [-2,  3]]), T)

    assert ss1 + ss2 == expected_add
    assert ss1*ss2 == expected_mul
    assert ss3 + 1/2 == DiscreteStateSpace(Matrix([[0]]), Matrix([[0]]),
                                           Matrix([[0]]), Matrix([[0.5]]), T)
    assert ss4*1.5 == DiscreteStateSpace(Matrix([[1]]), Matrix([[2]]),
                                         Matrix([[4.5]]), Matrix([[6.0]]), T)
    assert 1.5*ss4 == DiscreteStateSpace(Matrix([[1]]), Matrix([[3.0]]),
                                         Matrix([[3]]), Matrix([[6.0]]), T)

    raises(ShapeError, lambda: ss1 + ss3)
    raises(ShapeError, lambda: ss2*ss4)

    # compatibility tests
    ss5 = DiscreteStateSpace(A1, B1, C1, D1, 0.1)
    cont_ss = StateSpace(A2, B2, C2, D2)

    raises(TypeError, lambda: ss5 + ss2) # sum of state spaces with different sampling time
    raises(TypeError, lambda: ss5 * ss2) # multiplication of state spaces with different sampling time
    raises(TypeError, lambda: ss1 + cont_ss) # sum of state space and continuous state space
    raises(TypeError, lambda: ss1 * cont_ss) # multiplication of state space and continuous state space

def test_DiscreteStateSpace_negation():
    A = Matrix([[a0, a1], [a2, a3]])
    B = Matrix([[b0, b1], [b2, b3]])
    C = Matrix([[c0, c1], [c1, c2], [c2, c3]])
    D = Matrix([[d0, d1], [d1, d2], [d2, d3]])
    SS = DiscreteStateSpace(A, B, C, D, T)
    SS_neg = -SS

    state_mat = Matrix([[-1, 1], [1, -1]])
    input_mat = Matrix([1, -1])
    output_mat = Matrix([[-1, 1]])
    feedforward_mat = Matrix([1])
    system = DiscreteStateSpace(state_mat, input_mat, output_mat,
                                feedforward_mat, T)

    assert SS_neg == \
        DiscreteStateSpace(Matrix([[a0, a1],
                           [a2, a3]]),
                   Matrix([[b0, b1],
                           [b2, b3]]),
                   Matrix([[-c0, -c1],
                           [-c1, -c2],
                           [-c2, -c3]]),
                   Matrix([[-d0, -d1],
                           [-d1, -d2],
                           [-d2, -d3]]), T)
    assert -system == \
        DiscreteStateSpace(Matrix([[-1,  1],
                           [ 1, -1]]),
                   Matrix([[ 1],[-1]]),
                   Matrix([[1, -1]]),
                   Matrix([[-1]]), T)
    assert -SS_neg == SS
    assert -(-(-(-system))) == system

def test_DiscreteStateSpace_functions():
    # https://in.mathworks.com/help/control/ref/statespacemodel.obsv.html
    A_mat = Matrix([[-1.5, -2], [1, 0]])
    B_mat = Matrix([0.5, 0])
    C_mat = Matrix([[0, 1]])
    D_mat = Matrix([1])
    SS1 = DiscreteStateSpace(A_mat, B_mat, C_mat, D_mat)
    SS2 = DiscreteStateSpace(Matrix([[1, 1], [4, -2]]),
                             Matrix([[0, 1], [0, 2]]),
                             Matrix([[-1, 1], [1, -1]]))
    SS3 = DiscreteStateSpace(Matrix([[1, 1], [4, -2]]),
                             Matrix([[1, -1], [1, -1]]))
    SS4 = DiscreteStateSpace(Matrix([[a0, a1], [a2, a3]]),
                             Matrix([[b1], [b2]]), Matrix([[c1, c2]]))

    # Observability
    assert SS1.is_observable() == True
    assert SS2.is_observable() == False
    assert SS1.observability_matrix() == Matrix([[0, 1], [1, 0]])
    assert SS2.observability_matrix() == Matrix([[-1,  1], [ 1, -1], [ 3, -3], [-3,  3]])
    assert SS1.observable_subspace() == [Matrix([[1], [0]]), Matrix([[0], [1]])]
    assert SS2.observable_subspace() == [Matrix([[-1], [ 1]])]
    raises(NotImplementedError, lambda: SS4.observable_subspace())
    assert SS1.unobservable_subspace() == []
    assert SS2.unobservable_subspace() == [Matrix([[1],[1]])]
    raises(NotImplementedError, lambda: SS4.unobservable_subspace())

    Qo = SS4.observability_matrix().subs([(a0, 0), (a1, -6), (a2, 1), (a3, -5), (c1, 0), (c2, 1)])
    assert Qo == Matrix([[0, 1], [1, -5]])

    ss_obs = StateSpace(Matrix([[1, 0, 1], [0,0,0],[0,0,-2]]), Matrix([1,1,0]),
                        Matrix([1,1, Rational(1,3)]).T).to_observable_form()
    A_obs = ss_obs.A[:2, :2]
    A_nobs = ss_obs.A[2:, 2:]
    assert A_obs.eigenvals() == {0: 1, 1: 1}
    assert A_nobs.eigenvals() == {-2: 1}

    # Controllability
    assert SS1.is_controllable() == True
    assert SS3.is_controllable() == False
    assert SS1.controllability_matrix() ==  Matrix([[0.5, -0.75], [  0,   0.5]])
    assert SS3.controllability_matrix() == Matrix([[1, -1, 2, -2], [1, -1, 2, -2]])
    assert SS4.controllability_matrix() == \
        Matrix([[b1, a0*b1 + a1*b2], [b2, a2*b1 + a3*b2]])
    assert SS1.controllable_subspace() == [Matrix([[0.5], [  0]]), Matrix([[-0.75], [  0.5]])]
    assert SS3.controllable_subspace() == [Matrix([[1], [1]])]
    assert SS1.uncontrollable_subspace() == []
    assert SS3.uncontrollable_subspace() == [Matrix([[-1], [1]])]
    raises(NotImplementedError, lambda: SS4.uncontrollable_subspace()) # uncontrollable subspace fo symbols not implemented

    Qc = SS4.controllability_matrix().subs([(a0, 0), (a1, 1), (a2, -6), (a3, -5), (b1, 0), (b2, 1)])
    assert Qc == Matrix([[0, 1], [1, -5]])
    ss_contr = StateSpace(Matrix([[1, 0, 1], [0,0,0],[0,0,-2]]), Matrix([1,1,0]), Matrix([1,1,0]).T).to_controllable_form()
    A_contr = ss_contr.A[:2, :2]
    A_ncontr = ss_contr.A[2:, 2:]
    assert A_contr.eigenvals() == {0: 1, 1: 1}
    assert A_ncontr.eigenvals() == {-2: 1}

    # Append
    A1 = Matrix([[0, 1], [1, 0]])
    B1 = Matrix([[0], [1]])
    C1 = Matrix([[0, 1]])
    D1 = Matrix([[0]])
    ss1 = DiscreteStateSpace(A1, B1, C1, D1)
    ss2 = DiscreteStateSpace(Matrix([[1, 0], [0, 1]]), Matrix([[1], [0]]),
                     Matrix([[1, 0]]), Matrix([[1]]))
    ss3 = ss1.append(ss2)
    ss4 = SS4.append(ss1)

    assert ss3.num_states == ss1.num_states + ss2.num_states
    assert ss3.num_inputs == ss1.num_inputs + ss2.num_inputs
    assert ss3.num_outputs == ss1.num_outputs + ss2.num_outputs
    assert ss3.state_matrix == Matrix([[0, 1, 0, 0], [1, 0, 0, 0],
                                       [0, 0, 1, 0], [0, 0, 0, 1]])
    assert ss3.input_matrix == Matrix([[0, 0], [1, 0], [0, 1], [0, 0]])
    assert ss3.output_matrix == Matrix([[0, 1, 0, 0], [0, 0, 1, 0]])
    assert ss3.feedforward_matrix == Matrix([[0, 0], [0, 1]])

    cont_ss1 = StateSpace(A1, B1, C1, D1)
    ss1_ = DiscreteStateSpace(A1, B1, C1, D1, sampling_time=2)
    raises(TypeError, lambda: ss1.append(cont_ss1))  # TypeError if appending a continuous system to a discrete system
    raises(TypeError, lambda: ss1.append(ss1_))  # TypeError if appending a discrete system with different sampling time

    # Using symbolic matrices
    assert ss4.num_states == SS4.num_states + ss1.num_states
    assert ss4.num_inputs == SS4.num_inputs + ss1.num_inputs
    assert ss4.num_outputs == SS4.num_outputs + ss1.num_outputs
    assert ss4.state_matrix == Matrix([[a0, a1, 0, 0], [a2, a3, 0, 0],
                                       [0, 0, 0, 1], [0, 0, 1, 0]])
    assert ss4.input_matrix == Matrix([[b1, 0], [b2, 0], [0, 0], [0, 1]])
    assert ss4.output_matrix == Matrix([[c1, c2, 0, 0], [0, 0, 0, 1]])
    assert ss4.feedforward_matrix == Matrix([[0, 0], [0, 0]])

def test_DiscreteStateSpace_series():
    # For SISO Systems
    a1 = Matrix([[0, 1], [1, 0]])
    b1 = Matrix([[0], [1]])
    c1 = Matrix([[0, 1]])
    d1 = Matrix([[0]])
    a2 = Matrix([[1, 0], [0, 1]])
    b2 = Matrix([[1], [0]])
    c2 = Matrix([[1, 0]])
    d2 = Matrix([[1]])

    ss1 = DiscreteStateSpace(a1, b1, c1, d1, T)
    ss2 = DiscreteStateSpace(a2, b2, c2, d2, T)
    tf1 = DiscreteTransferFunction(s, s+1, s, T)
    ser1 = Series(ss1, ss2)
    assert ser1 == Series(
        DiscreteStateSpace(
            Matrix([[0, 1], [1, 0]]), Matrix([[0], [1]]), Matrix([[0, 1]]),
            Matrix([[0]]), T),
        DiscreteStateSpace(
            Matrix([[1, 0], [0, 1]]), Matrix([[1], [0]]), Matrix([[1, 0]]),
            Matrix([[1]]), T),)
    assert ser1.doit() == DiscreteStateSpace(
        Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 1, 1, 0], [0, 0, 0, 1]]),
        Matrix([[0], [1], [0], [0]]), Matrix([[0, 1, 1, 0]]), Matrix([[0]]), T)

    assert ser1.num_inputs == 1
    assert ser1.num_outputs == 1
    assert ser1.rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction(z**2, z**3 - z**2 - z + 1, z, T)
    raises(TypeError, lambda: ser1.rewrite(TransferFunction))  # Cannot rewrite to TransferFunction
    ser2 = Series(ss1)
    ser3 = Series(ser2, ss2)
    assert ser3.doit() == ser1.doit()

    # DiscreteTransferFunction interconnection with DiscreteStateSpace
    ser_tf = Series(tf1, ss1)
    assert ser_tf == Series(DiscreteTransferFunction(s, s + 1, s, T),
                            DiscreteStateSpace(Matrix([
                            [0, 1],[1, 0]]), Matrix([[0], [1]]),
                            Matrix([[0, 1]]), Matrix([[0]]), T))
    assert ser_tf.doit() == DiscreteStateSpace(
                            Matrix([[-1, 0,  0], [0, 0,  1], [-1, 1, 0]]),
                            Matrix([[1], [0], [1]]), Matrix([[0, 0, 1]]),
                            Matrix([[0]]), T)
    assert ser_tf.rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction(z**2, z**3 + z**2 - z - 1, z, T)

    # For MIMO Systems
    a3 = Matrix([[4, 1], [2, -3]])
    b3 = Matrix([[5, 2], [-3, -3]])
    c3 = Matrix([[2, -4], [0, 1]])
    d3 = Matrix([[3, 2], [1, -1]])
    a4 = Matrix([[-3, 4, 2], [-1, -3, 0], [2, 5, 3]])
    b4 = Matrix([[1, 4], [-3, -3], [-2, 1]])
    c4 = Matrix([[4, 2, -3], [1, 4, 3]])
    d4 = Matrix([[-2, 4], [0, 1]])
    ss3 = DiscreteStateSpace(a3, b3, c3, d3, T)
    ss4 = DiscreteStateSpace(a4, b4, c4, d4, T)
    ser4 = MIMOSeries(ss3, ss4)
    assert ser4 == MIMOSeries(DiscreteStateSpace(Matrix([
                    [4,  1],
                    [2, -3]]), Matrix([
                    [ 5,  2],
                    [-3, -3]]), Matrix([
                    [2, -4],
                    [0,  1]]), Matrix([
                    [3,  2],
                    [1, -1]]), T), DiscreteStateSpace(Matrix([
                    [-3,  4, 2],
                    [-1, -3, 0],
                    [ 2,  5, 3]]), Matrix([
                    [ 1,  4],
                    [-3, -3],
                    [-2,  1]]), Matrix([
                    [4, 2, -3],
                    [1, 4,  3]]), Matrix([
                    [-2, 4],
                    [ 0, 1]]), T))
    assert ser4.doit() == DiscreteStateSpace(
                        Matrix([
                        [4,   1,  0, 0,  0],
                        [2,  -3,  0, 0,  0],
                        [2,   0,  -3, 4,  2],
                        [-6,  9, -1, -3,  0],
                        [-4, 9,  2, 5, 3]]),
                        Matrix([
                        [5,   2],
                        [-3,  -3],
                        [7,   -2],
                        [-12,  -3],
                        [-5, -5]]),
                        Matrix([
                        [-4, 12, 4, 2, -3],
                        [0, 1, 1, 4, 3]]),
                        Matrix([
                        [-2, -8],
                        [1, -1]]), T)
    assert ser4.num_inputs == ss3.num_inputs
    assert ser4.num_outputs == ss4.num_outputs
    ser5 = MIMOSeries(ss3)
    ser6 = MIMOSeries(ser5, ss4)
    assert ser6.doit() == ser4.doit()
    assert ser6.rewrite(TransferFunctionMatrix) == \
        ser4.rewrite(TransferFunctionMatrix)
    tf2 = DiscreteTransferFunction(1, s, s, T)
    tf3 = DiscreteTransferFunction(1, s+1, s, T)
    tf4 = DiscreteTransferFunction(s, s+2, s, T)
    tfm = TransferFunctionMatrix([[tf1, tf2], [tf3, tf4]])
    ser6 = MIMOSeries(ss3, tfm)
    assert ser6 == MIMOSeries(DiscreteStateSpace(Matrix([
                        [4,  1],
                        [2, -3]]), Matrix([
                        [ 5,  2],
                        [-3, -3]]), Matrix([
                        [2, -4],
                        [0,  1]]), Matrix([
                        [3,  2],
                        [1, -1]]), T), TransferFunctionMatrix((
                        (DiscreteTransferFunction(s, s + 1, s, T),
                         DiscreteTransferFunction(1, s, s, T)),
                        (DiscreteTransferFunction(1, s + 1, s, T),
                         DiscreteTransferFunction(s, s + 2, s, T)))))

    # compatibility tests
    ss5 = DiscreteStateSpace(a1, b1, c1, d1, 0.2)
    cont_ss1 = StateSpace(a2, b2, c2, d2)

    # SISO
    raises(TypeError, lambda: Series(ss5, ss2))  # series with different sampling time
    raises(TypeError, lambda: Series(ss1, cont_ss1))  # series with continuous state space

    # MIMO
    ss6 = DiscreteStateSpace(a3, b3, c3, d3)
    cont_ss2 = StateSpace(a4, b4, c4, d4)
    raises(TypeError, lambda: MIMOSeries(ss6, ss4))  # series with different sampling time
    raises(TypeError, lambda: MIMOSeries(ss3, cont_ss2))  # series with continuous state space

def test_DiscreteStateSpace_parallel():
    # For SISO system
    a1 = Matrix([[0, 1], [1, 0]])
    b1 = Matrix([[0], [1]])
    c1 = Matrix([[0, 1]])
    d1 = Matrix([[0]])
    a2 = Matrix([[1, 0], [0, 1]])
    b2 = Matrix([[1], [0]])
    c2 = Matrix([[1, 0]])
    d2 = Matrix([[1]])
    ss1 = DiscreteStateSpace(a1, b1, c1, d1, 0.12)
    ss2 = DiscreteStateSpace(a2, b2, c2, d2, 0.12)
    p1 = Parallel(ss1, ss2)
    assert p1 == Parallel(DiscreteStateSpace(
        Matrix([[0, 1], [1, 0]]), Matrix([[0], [1]]),
        Matrix([[0, 1]]), Matrix([[0]]), 0.12),
                          DiscreteStateSpace(
        Matrix([[1, 0],[0, 1]]), Matrix([[1],[0]]),
        Matrix([[1, 0]]), Matrix([[1]]), 0.12))
    assert p1.doit() == DiscreteStateSpace(Matrix([
                        [0, 1, 0, 0],
                        [1, 0, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]]),
                        Matrix([
                        [0],
                        [1],
                        [1],
                        [0]]),
                        Matrix([[0, 1, 1, 0]]),
                        Matrix([[1]]), 0.12)
    assert p1.rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction(z*(z + 2), z**2 - 1, z, 0.12)
    raises(TypeError, lambda: p1.rewrite(TransferFunction))  # cannot rewrite to TransferFunction

    # Connecting DiscreteStateSpace with DiscreteTransferFunction
    tf1 = DiscreteTransferFunction(s, s+1, s, 0.12)
    p2 = Parallel(ss1, tf1)
    assert p2 == Parallel(DiscreteStateSpace(Matrix([
                        [0, 1],
                        [1, 0]]), Matrix([
                        [0],
                        [1]]), Matrix([[0, 1]]), Matrix([[0]]), 0.12),
                        DiscreteTransferFunction(s, s + 1, s, 0.12))
    assert p2.doit() == DiscreteStateSpace(
                        Matrix([
                        [0, 1,  0],
                        [1, 0,  0],
                        [0, 0, -1]]),
                        Matrix([
                        [0],
                        [1],
                        [1]]),
                        Matrix([[0, 1, -1]]),
                        Matrix([[1]]), 0.12)
    assert p2.rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction(z**2, z**2 - 1, z, 0.12)

    # For MIMO
    a3 = Matrix([[4, 1], [2, -3]])
    b3 = Matrix([[5, 2], [-3, -3]])
    c3 = Matrix([[2, -4], [0, 1]])
    d3 = Matrix([[3, 2], [1, -1]])
    a4 = Matrix([[-3, 4, 2], [-1, -3, 0], [2, 5, 3]])
    b4 = Matrix([[1, 4], [-3, -3], [-2, 1]])
    c4 = Matrix([[4, 2, -3], [1, 4, 3]])
    d4 = Matrix([[-2, 4], [0, 1]])
    ss3 = DiscreteStateSpace(a3, b3, c3, d3, 0.12)
    ss4 = DiscreteStateSpace(a4, b4, c4, d4, 0.12)
    p3 = MIMOParallel(ss3, ss4)
    assert p3 == MIMOParallel(DiscreteStateSpace(Matrix([
                        [4,  1],
                        [2, -3]]), Matrix([
                        [ 5,  2],
                        [-3, -3]]), Matrix([
                        [2, -4],
                        [0,  1]]), Matrix([
                        [3,  2],
                        [1, -1]]), 0.12), DiscreteStateSpace(Matrix([
                        [-3,  4, 2],
                        [-1, -3, 0],
                        [ 2,  5, 3]]), Matrix([
                        [ 1,  4],
                        [-3, -3],
                        [-2,  1]]), Matrix([
                        [4, 2, -3],
                        [1, 4,  3]]), Matrix([
                        [-2, 4],
                        [ 0, 1]]), 0.12))
    assert p3.doit() == DiscreteStateSpace(Matrix([
                        [4, 1, 0, 0, 0],
                        [2, -3, 0, 0, 0],
                        [0, 0, -3, 4, 2],
                        [0, 0, -1, -3, 0],
                        [0, 0, 2, 5, 3]]),
                        Matrix([
                        [5, 2],
                        [-3, -3],
                        [1, 4],
                        [-3, -3],
                        [-2, 1]]),
                        Matrix([
                        [2, -4, 4, 2, -3],
                        [0, 1, 1, 4, 3]]),
                        Matrix([
                        [1, 6],
                        [1, 0]]), 0.12)

    # Using StateSpace with MIMOParallel.
    tf2 = DiscreteTransferFunction(1, s, s, 0.12)
    tf3 = DiscreteTransferFunction(1, s + 1, s, 0.12)
    tf4 = DiscreteTransferFunction(s, s + 2, s, 0.12)
    tfm = TransferFunctionMatrix([[tf1, tf2], [tf3, tf4]])
    p4 = MIMOParallel(tfm, ss3)
    assert p4 == MIMOParallel(TransferFunctionMatrix((
                        (DiscreteTransferFunction(s, s + 1, s, 00.12),
                         DiscreteTransferFunction(1, s, s, 0.12)),
                        (DiscreteTransferFunction(1, s + 1, s, 0.12),
                         DiscreteTransferFunction(s, s + 2, s, 0.12)))),
                        DiscreteStateSpace(Matrix([
                        [4, 1],
                        [2, -3]]), Matrix([
                        [5, 2],
                        [-3, -3]]), Matrix([
                        [2, -4],
                        [0, 1]]), Matrix([
                        [3, 2],
                        [1, -1]]), 0.12))

    # compatibility tests
    ss5 = DiscreteStateSpace(a1, b1, c1, d1, 0.33)
    cont_ss1 = DiscreteStateSpace(a2, b2, c2, d2)

    # SISO
    raises(TypeError, lambda: Parallel(ss5, ss2))  # parallel with different sampling time
    raises(TypeError, lambda: Parallel(ss1, cont_ss1))  # parallel with continuous state space

    # MIMO
    ss6 = DiscreteStateSpace(a3, b3, c3, d3, 0.33)
    cont_ss2 = DiscreteStateSpace(a4, b4, c4, d4)
    raises(TypeError, lambda: MIMOParallel(ss6, ss4))  # parallel with different sampling time
    raises(TypeError, lambda: MIMOParallel(ss3, cont_ss2))  # parallel with continuous state space

def test_DiscreteStateSpace_feedback():
    # For SISO
    a1 = Matrix([[0, 1], [1, 0]])
    b1 = Matrix([[0], [1]])
    c1 = Matrix([[0, 1]])
    d1 = Matrix([[0]])
    a2 = Matrix([[1, 0], [0, 1]])
    b2 = Matrix([[1], [0]])
    c2 = Matrix([[1, 0]])
    d2 = Matrix([[1]])
    ss1 = DiscreteStateSpace(a1, b1, c1, d1, T)
    ss2 = DiscreteStateSpace(a2, b2, c2, d2, T)
    fd1 = Feedback(ss1, ss2)

    # Negative feedback
    assert fd1 == Feedback(
        DiscreteStateSpace(
            Matrix([[0, 1], [1, 0]]), Matrix([[0], [1]]),
            Matrix([[0, 1]]), Matrix([[0]]), T),
        DiscreteStateSpace(
            Matrix([[1, 0],[0, 1]]), Matrix([[1],[0]]),
            Matrix([[1, 0]]), Matrix([[1]]), T), -1)
    assert fd1.doit() == DiscreteStateSpace(Matrix([
                            [0,  1,  0, 0],
                            [1, -1, -1, 0],
                            [0,  1,  1, 0],
                            [0,  0,  0, 1]]), Matrix([
                            [0],
                            [1],
                            [0],
                            [0]]), Matrix(
                            [[0, 1, 0, 0]]), Matrix(
                            [[0]]), T)
    assert fd1.rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction(z*(z - 1), z**3 - z + 1, z, T)
    raises(TypeError, lambda: fd1.rewrite(TransferFunction))  # Cannot rewrite to TransferFunction

    # Positive Feedback
    fd2 = Feedback(ss1, ss2, 1)
    assert fd2.doit() == DiscreteStateSpace(Matrix([
                            [0, 1, 0, 0],
                            [1, 1, 1, 0],
                            [0, 1, 1, 0],
                            [0, 0, 0, 1]]), Matrix([
                            [0],
                            [1],
                            [0],
                            [0]]), Matrix(
                            [[0, 1, 0, 0]]), Matrix(
                            [[0]]), T)
    assert fd2.rewrite(DiscreteTransferFunction) == \
        DiscreteTransferFunction(z*(z - 1), z**3 - 2*z**2 - z + 1, z, T)

    # Connection with DiscreteTransferFunction
    tf1 = DiscreteTransferFunction(s, s+1, s, T)
    fd3 = Feedback(ss1, tf1)
    assert fd3 == Feedback(DiscreteStateSpace(Matrix([
                            [0, 1],
                            [1, 0]]), Matrix([
                            [0],
                            [1]]), Matrix([[0, 1]]), Matrix([[0]]), T),
                            DiscreteTransferFunction(s, s + 1, s, T), -1)
    assert fd3.doit() == DiscreteStateSpace(Matrix([
                            [0,  1,  0],
                            [1, -1,  1],
                            [0,  1, -1]]), Matrix([
                            [0],
                            [1],
                            [0]]), Matrix(
                            [[0, 1, 0]]), Matrix(
                            [[0]]), T)

    # For MIMO
    a3 = Matrix([[4, 1], [2, -3]])
    b3 = Matrix([[5, 2], [-3, -3]])
    c3 = Matrix([[2, -4], [0, 1]])
    d3 = Matrix([[3, 2], [1, -1]])
    a4 = Matrix([[-3, 4, 2], [-1, -3, 0], [2, 5, 3]])
    b4 = Matrix([[1, 4], [-3, -3], [-2, 1]])
    c4 = Matrix([[4, 2, -3], [1, 4, 3]])
    d4 = Matrix([[-2, 4], [0, 1]])
    ss3 = DiscreteStateSpace(a3, b3, c3, d3, T)
    ss4 = DiscreteStateSpace(a4, b4, c4, d4, T)

    # Negative Feedback
    fd4 = MIMOFeedback(ss3, ss4)
    assert fd4 == MIMOFeedback(DiscreteStateSpace(Matrix([
                            [4,  1],
                            [2, -3]]), Matrix([
                            [ 5,  2],
                            [-3, -3]]), Matrix([
                            [2, -4],
                            [0,  1]]), Matrix([
                            [3,  2],
                            [1, -1]]), T), DiscreteStateSpace(Matrix([
                            [-3,  4, 2],
                            [-1, -3, 0],
                            [ 2,  5, 3]]), Matrix([
                            [ 1,  4],
                            [-3, -3],
                            [-2,  1]]), Matrix([
                            [4, 2, -3],
                            [1, 4,  3]]), Matrix([
                            [-2, 4],
                            [ 0, 1]]), T), -1)
    assert fd4.doit() == DiscreteStateSpace(Matrix([
                            [Rational(3), Rational(-3, 4), Rational(-15, 4),
                            Rational(-37, 2), Rational(-15)],
                            [Rational(7, 2), Rational(-39, 8), Rational(9, 8),
                             Rational(39, 4), Rational(9)],
                            [Rational(3), Rational(-41, 4), Rational(-45, 4),
                            Rational(-51, 2), Rational(-19)],
                            [Rational(-9, 2), Rational(129, 8), Rational(73, 8),
                            Rational(171, 4), Rational(36)],
                            [Rational(-3, 2), Rational(47, 8), Rational(31, 8),
                             Rational(85, 4), Rational(18)]]), Matrix([
                            [Rational(-1, 4), Rational(19, 4)],
                            [Rational(3, 8), Rational(-21, 8)],
                            [Rational(1, 4), Rational(29, 4)],
                            [Rational(3, 8), Rational(-93, 8)],
                            [Rational(5, 8), Rational(-35, 8)]]), Matrix([
                            [Rational(1), Rational(-15, 4), Rational(-7, 4),
                             Rational(-21, 2), Rational(-9)],
                            [Rational(1, 2), Rational(-13, 8), Rational(-13, 8),
                             Rational(-19, 4), Rational(-3)]]), Matrix([
                            [Rational(-1, 4), Rational(11, 4)],
                            [Rational(1, 8), Rational(9, 8)]]), T)

    # Positive Feedback
    fd5 = MIMOFeedback(ss3, ss4, 1)
    assert fd5.doit() == DiscreteStateSpace(Matrix([
                            [Rational(4, 7), Rational(62, 7), Rational(1),
                             Rational(-8), Rational(-69, 7)],
                            [Rational(32, 7), Rational(-135, 14),
                             Rational(-3, 2), Rational(3), Rational(36, 7)],
                            [Rational(-10, 7), Rational(41, 7), Rational(-4),
                             Rational(-12), Rational(-97, 7)],
                            [Rational(12, 7), Rational(-111, 14),
                             Rational(-5, 2), Rational(18), Rational(171, 7)],
                            [Rational(2, 7), Rational(-29, 14), Rational(-1, 2),
                             Rational(10), Rational(81, 7)]]), Matrix([
                            [Rational(6, 7), Rational(-17, 7)],
                            [Rational(-9, 14), Rational(15, 14)],
                            [Rational(6, 7), Rational(-31, 7)],
                            [Rational(-27, 14), Rational(87, 14)],
                            [Rational(-15, 14), Rational(25, 14)]]), Matrix([
                            [Rational(-2, 7), Rational(11, 7), Rational(1),
                             Rational(-4), Rational(-39, 7)],
                            [Rational(-2, 7), Rational(15, 14), Rational(-1, 2),
                             Rational(-3), Rational(-18, 7)]]), Matrix([
                            [Rational(4, 7), Rational(-9, 7)],
                            [Rational(1, 14), Rational(-11, 14)]]), T)

    # compatibility tests
    ss5 = DiscreteStateSpace(a1, b1, c1, d1, 0.12)
    cont_ss1 = DiscreteStateSpace(a2, b2, c2, d2)

    # SISO
    raises(TypeError, lambda: Feedback(ss5, ss2))  # feedback with different sampling time
    raises(TypeError, lambda: Feedback(ss1, cont_ss1))  # feedback with continuous state space

    # MIMO
    ss6 = DiscreteStateSpace(a3, b3, c3, d3, 0.12)
    cont_ss2 = DiscreteStateSpace(a4, b4, c4, d4)
    raises(TypeError, lambda: MIMOFeedback(ss6, ss4))  # feedback with different sampling time
    raises(TypeError, lambda: MIMOFeedback(ss3, cont_ss2))  # feedback with continuous state space

def test_DiscreteStateSpace_stability():
    k = symbols('k')
    B = Matrix([1, 0, 0])
    C = Matrix([[0, 1, 0]])
    D = Matrix([0])

    A1 = Matrix([[0,1,0],[0,0,1], [k-1, -2*k, -1]])
    ss1 = DiscreteStateSpace(A1, B, C, D)
    ineq = ss1.get_asymptotic_stability_conditions()
    assert ineq == [-15*k**2 + 20*k - 5 > 0, -8*k**2 - 8*k + 8 > 0,
                    3*k**2 + 8*k - 3 > 0]

    A2 = Matrix([[1,0,0], [0,-1,k], [0,0,-1]])
    ss2 = DiscreteStateSpace(A2, B, C, D)
    ineq = ss2.get_asymptotic_stability_conditions()
    assert ineq == [False]
