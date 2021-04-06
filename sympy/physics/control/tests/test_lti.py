from sympy import symbols, Matrix, factor, Function, simplify, exp, pi, oo, I, \
    Rational, sqrt, CRootOf
from sympy.physics.control.lti import TransferFunction, Series, Parallel, Feedback
from sympy.testing.pytest import raises

a, x, b, s, g, d, p, k, a0, a1, a2, b0, b1, b2 = symbols('a, x, b, s, g, d, p, k, a0:3, b0:3')


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
    tf18 = TransferFunction((k_p + k_o*s + k_i/s), s**2 + 2*omega_o*s + omega_o**2, s)
    assert tf18.num == k_i/s + k_o*s + k_p
    assert tf18.args == (k_i/s + k_o*s + k_p, omega_o**2 + 2*omega_o*s + s**2, s)

    # ValueError when denominator is zero.
    raises(ValueError, lambda: TransferFunction(4, 0, s))
    raises(ValueError, lambda: TransferFunction(s, 0, s))
    raises(ValueError, lambda: TransferFunction(0, 0, s))

    raises(TypeError, lambda: TransferFunction(Matrix([1, 2, 3]), s, s))
    raises(TypeError, lambda: TransferFunction(s**pi*exp(s), s, s))

    raises(TypeError, lambda: TransferFunction(s**2 + 2*s - 1, s + 3, 3))
    raises(TypeError, lambda: TransferFunction(p + 1, 5 - p, 4))
    raises(TypeError, lambda: TransferFunction(3, 4, 8))


def test_TransferFunction_functions():
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

    # purely symbolic polynomials.
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
    assert tf1.is_stable() is False   # as one pole is +ve, and the other is -ve.
    assert expect2.is_stable() is False
    assert expect1.is_stable() is True
    assert stable_tf.is_stable() is True
    assert stable_tf_.is_stable() is True
    assert TF_.is_stable() is False
    assert expect4_.is_stable() is None   # no assumption provided for the only pole 's'.
    assert SP4.is_stable() is None

    # Zeros of a transfer function.
    assert G1.zeros() == [1, 1]
    assert G2.zeros() == []
    assert tf1.zeros() == [-3, 1]
    assert expect4_.zeros() == [7**(Rational(2, 3))*(-s)**(Rational(1, 3))/7, -7**(Rational(2, 3))*(-s)**(Rational(1, 3))/14 -
        sqrt(3)*7**(Rational(2, 3))*I*(-s)**(Rational(1, 3))/14, -7**(Rational(2, 3))*(-s)**(Rational(1, 3))/14 + sqrt(3)*7**(Rational(2, 3))*I*(-s)**(Rational(1, 3))/14]
    assert SP4.zeros() == [(s/a1)**(Rational(1, 3)), -(s/a1)**(Rational(1, 3))/2 - sqrt(3)*I*(s/a1)**(Rational(1, 3))/2,
        -(s/a1)**(Rational(1, 3))/2 + sqrt(3)*I*(s/a1)**(Rational(1, 3))/2]
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

    # sympy's own functions.
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
    raises(ValueError, lambda: G5 / G6)
    raises(ValueError, lambda: -G3 /G4)
    raises(ValueError, lambda: G7 / (1 + G6))
    raises(ValueError, lambda: G7 / (G5 * G6))
    raises(ValueError, lambda: G7 / (G7 + (G5 + G6)))


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


def test_Series_construction():
    zeta, wn = symbols('zeta, wn')
    tf = TransferFunction(a0*s**3 + a1*s**2 - a2*s, b0*p**4 + b1*p**3 - b2*s*p, s)
    tf2 = TransferFunction(a2*p - s, a2*s + p, s)
    tf3 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf4 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    inp = Function('X_d')(s)
    out = Function('X')(s)

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


def test_Series_functions():
    zeta, wn = symbols('zeta, wn')
    tf1 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    tf2 = TransferFunction(k, 1, s)
    tf3 = TransferFunction(a2*p - s, a2*s + p, s)
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)

    assert tf1*tf2*tf3 == Series(tf1, tf2, tf3)
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


def test_Parallel_construction():
    zeta, wn = symbols('zeta, wn')
    tf = TransferFunction(a0*s**3 + a1*s**2 - a2*s, b0*p**4 + b1*p**3 - b2*s*p, s)
    tf2 = TransferFunction(a2*p - s, a2*s + p, s)
    tf3 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf4 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    inp = Function('X_d')(s)
    out = Function('X')(s)

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


def test_Parallel_functions():
    zeta, wn = symbols('zeta, wn')
    tf1 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    tf2 = TransferFunction(k, 1, s)
    tf3 = TransferFunction(a2*p - s, a2*s + p, s)
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)

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
        Parallel(tf1, tf2, Series(-tf1, tf3)).doit()== TransferFunction((-a2*p + s)*(s**2 + 2*s*wn*zeta + wn**2) + \
        (a2*s + p)*(k*(s**2 + 2*s*wn*zeta + wn**2) + 1)*(s**2 + 2*s*wn*zeta + wn**2), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2)**2, s)
    assert Parallel(tf2, tf1, -tf3, evaluate=True) == Parallel(tf2, tf1, -tf3).doit() == \
        TransferFunction(-(a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*s + p)*(k*(s**2 + 2*s*wn*zeta + wn**2) + 1), \
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert not Parallel(tf1, -tf2, evaluate=False) == Parallel(tf1, -tf2).doit()

    assert Parallel(Series(tf1, tf2), Series(tf2, tf3)).doit() == \
        TransferFunction(k*(a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) + k*(a2*s + p), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Parallel(-tf1, -tf2, -tf3).doit() == \
        TransferFunction(-(a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) + \
        (a2*s + p)*(-k*(s**2 + 2*s*wn*zeta + wn**2) - 1), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert -Parallel(tf1, tf2, tf3).doit() == \
        TransferFunction(-((a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*s + p)*(k*(s**2 + 2*s*wn*zeta + wn**2) + 1)),
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Parallel(tf2, tf3, Series(tf2, -tf1), tf3).doit() == \
        TransferFunction((a2*p - s)*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*s + p)*(-k*(a2*s + p) + \
        (s**2 + 2*s*wn*zeta + wn**2)*(a2*p + k*(a2*s + p) - s)), (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2), s)

    assert Parallel(tf1, tf2).rewrite(TransferFunction) == \
        TransferFunction(k*(s**2 + 2*s*wn*zeta + wn**2) + 1, s**2 + 2*s*wn*zeta + wn**2, s)
    assert Parallel(tf2, tf1, -tf3).rewrite(TransferFunction) == \
        TransferFunction(-(a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*s + p)*(k*(s**2 + 2*s*wn*zeta + wn**2) + 1), \
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)

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


def test_Feedback_construction():
    zeta, wn = symbols('zeta, wn')
    tf1 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    tf2 = TransferFunction(k, 1, s)
    tf3 = TransferFunction(a2*p - s, a2*s + p, s)
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)
    tf6 = TransferFunction(s - p, p + s, p)

    f1 = Feedback(TransferFunction(1, 1, s), tf1*tf2*tf3)
    assert f1.args == (TransferFunction(1, 1, s), Series(tf1, tf2, tf3))
    assert f1.num == TransferFunction(1, 1, s)
    assert f1.den == Series(tf1, tf2, tf3)
    assert f1.var == s

    f2 = Feedback(tf1, tf2*tf3)
    assert f2.args == (tf1, Series(tf2, tf3))
    assert f2.num == tf1
    assert f2.den == Series(tf2, tf3)
    assert f2.var == s

    f3 = Feedback(tf1*tf2, tf5)
    assert f3.args == (Series(tf1, tf2), tf5)
    assert f3.num == Series(tf1, tf2)

    f4 = Feedback(tf4, tf6)
    assert f4.args == (tf4, tf6)
    assert f4.num == tf4
    assert f4.var == p

    f5 = Feedback(tf5, TransferFunction(1, 1, s))
    assert f5.args == (tf5, TransferFunction(1, 1, s))
    assert f5.var == s

    f6 = Feedback(TransferFunction(1, 1, p), tf4)
    assert f6.args == (TransferFunction(1, 1, p), tf4)
    assert f6.var == p

    f7 = -Feedback(tf4*tf6, TransferFunction(1, 1, p))
    assert f7.args == (Series(TransferFunction(-1, 1, p), Series(tf4, tf6)), TransferFunction(1, 1, p))
    assert f7.num == Series(TransferFunction(-1, 1, p), Series(tf4, tf6))

    # denominator can't be a Parallel instance
    raises(TypeError, lambda: Feedback(tf1, tf2 + tf3))
    raises(TypeError, lambda: Feedback(tf1, Matrix([1, 2, 3])))
    raises(TypeError, lambda: Feedback(TransferFunction(1, 1, s), s - 1))
    raises(TypeError, lambda: Feedback(1, 1))
    raises(ValueError, lambda: Feedback(TransferFunction(1, 1, s), TransferFunction(1, 1, s)))
    raises(ValueError, lambda: Feedback(tf2, tf4*tf5))


def test_Feedback_functions():
    zeta, wn = symbols('zeta, wn')
    tf = TransferFunction(1, 1, s)
    tf1 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    tf2 = TransferFunction(k, 1, s)
    tf3 = TransferFunction(a2*p - s, a2*s + p, s)
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)
    tf6 = TransferFunction(s - p, p + s, p)

    assert tf / (tf + tf1) == Feedback(tf, tf1)
    assert tf / (tf + tf1*tf2*tf3) == Feedback(tf, tf1*tf2*tf3)
    assert tf1 / (tf + tf1*tf2*tf3) == Feedback(tf1, tf2*tf3)
    assert (tf1*tf2) / (tf + tf1*tf2) == Feedback(tf1*tf2, tf)
    assert (tf1*tf2) / (tf + tf1*tf2*tf5) == Feedback(tf1*tf2, tf5)
    assert (tf1*tf2) / (tf + tf1*tf2*tf5*tf3) in (Feedback(tf1*tf2, tf5*tf3), Feedback(tf1*tf2, tf3*tf5))
    assert tf4 / (TransferFunction(1, 1, p) + tf4*tf6) == Feedback(tf4, tf6)
    assert tf5 / (tf + tf5) == Feedback(tf5, tf)

    raises(ValueError, lambda: tf1*tf2*tf3 / (1 + tf1*tf2*tf3))
    raises(ValueError, lambda: tf1*tf2*tf3 / tf3*tf5)
    raises(ValueError, lambda: tf2*tf3 / (tf + tf2*tf3*tf4))

    assert Feedback(tf, tf1*tf2*tf3).doit() == \
        TransferFunction((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), k*(a2*p - s) + \
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(tf1, tf2*tf3).doit() == \
        TransferFunction((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), (k*(a2*p - s) + \
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2))*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(tf1*tf2, tf5).doit() == \
        TransferFunction(k*(a0 + s)*(s**2 + 2*s*wn*zeta + wn**2), (k*(-a0 + a1*s**2 + a2*s) + \
        (a0 + s)*(s**2 + 2*s*wn*zeta + wn**2))*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(tf4, tf6).doit() == \
        TransferFunction(p*(p + s)*(a0*p + p**a1 - s), p*(p*(p + s) + (-p + s)*(a0*p + p**a1 - s)), p)
    assert -Feedback(tf4*tf6, TransferFunction(1, 1, p)).doit() == \
        TransferFunction(-p*(-p + s)*(p + s)*(a0*p + p**a1 - s), p*(p + s)*(p*(p + s) + (-p + s)*(a0*p + p**a1 - s)), p)

    assert Feedback(tf1, tf2*tf5).rewrite(TransferFunction) == \
        TransferFunction((a0 + s)*(s**2 + 2*s*wn*zeta + wn**2), (k*(-a0 + a1*s**2 + a2*s) + \
        (a0 + s)*(s**2 + 2*s*wn*zeta + wn**2))*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(TransferFunction(1, 1, p), tf4).rewrite(TransferFunction) == \
        TransferFunction(p, a0*p + p + p**a1 - s, p)
