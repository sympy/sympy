from sympy import symbols, Matrix, factor, simplify, exp, log, pi
from sympy.physics.control.lti import TransferFunction
from sympy.testing.pytest import raises

a, b, s, g, d, p, k, a0, a1, a2, b0, b1, b2 = symbols('a, b, s, g, d, p, k, a0:3, b0:3')


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

    # when just the numerator is 0, leave the denominator along.
    tf12 = TransferFunction(0, p**2 - p + 1, p)
    assert tf12.args == (0, p**2 - p + 1, p)

    tf13 = TransferFunction(0, 1, s)
    assert tf13.args == (0, 1, s)

    # float exponents
    tf14 = TransferFunction(a0*s**0.5 + a2*s**0.6 - a1, a1*p**(-8.7), s)
    assert tf14.args == (a0*s**0.5 - a1 + a2*s**0.6, a1*p**(-8.7), s)

    tf15 = TransferFunction(a2**2*p**(1/4) + a1*s**(-4/5), a0*s - p, p)
    assert tf15.args == (a1*s**(-0.8) + a2**2*p**0.25, a0*s - p, p)

    tf16 = TransferFunction(s**pi*exp(s)/log(s), s, s)
    assert tf16.num == s**pi*exp(s)/log(s)
    assert tf16.args == (s**pi*exp(s)/log(s), s, s)

    tau, omega_o, k_p, k_o, k_i = symbols('tau, omega_o, k_p, k_o, k_i')
    tf17 = TransferFunction(exp(-tau*s), log(p)**s*s, s)
    assert tf17.num == exp(-s*tau)
    assert tf17.den == s*log(p)**s
    assert tf17.args == (exp(-s*tau), s*log(p)**s, s)

    tf18 = TransferFunction((k_p + k_o*s + k_i/s), s**2 + 2*omega_o*s + omega_o**2, s)
    assert tf18.num == k_i/s + k_o*s + k_p
    assert tf18.args == (k_i/s + k_o*s + k_p, omega_o**2 + 2*omega_o*s + s**2, s)

    # ValueError when denominator is zero.
    raises(ValueError, lambda: TransferFunction(4, 0, s))
    raises(ValueError, lambda: TransferFunction(s, 0, s))
    raises(ValueError, lambda: TransferFunction(0, 0, s))

    raises(TypeError, lambda: TransferFunction(Matrix([1, 2, 3]), s, s))

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
    assert SP1.subs({a0: 1, a1: 2, b0: 3, b1: 4, b2: 5}).evalf() == expect1
    assert expect1_.evalf() == expect1

    tau, c1, d0, d1, d2 = symbols('tau, c1, d0:3')
    p3, p4 = c1*p, d2*p**3 + d1*p**2 - d0
    SP2 = TransferFunction(p3, p4, p)
    expect2 = TransferFunction(2.0*p, 5.0*p**3 + 2.0*p**2 - 3.0, p)
    expect2_ = TransferFunction(2*p, 5*p**3 + 2*p**2 - 3, p)
    assert SP2.subs({c1: 2, d0: 3, d1: 2, d2: 5}).evalf() == expect2
    assert expect2_.evalf() == expect2

    SP3 = TransferFunction(a0*p**3 + a1*s**2 - b0*s + b1, a1*s + p, s)
    expect3 = TransferFunction(2.0*p**3 + 4.0*s**2 - s + 5.0, p + 4.0*s, s)
    expect3_ = TransferFunction(2*p**3 + 4*s**2 - s + 5, p + 4*s, s)
    assert SP3.subs({a0: 2, a1: 4, b0: 1, b1: 5}).evalf() == expect3
    assert expect3_.evalf() == expect3

    SP4 = TransferFunction(s - a1*p**3, a0*s + p, p)
    expect4 = TransferFunction(7.0*p**3 + s, p - s, p)
    expect4_ = TransferFunction(7*p**3 + s, p - s, p)
    assert SP4.subs({a0: -1, a1: -7}).evalf() == expect4
    assert expect4_.evalf() == expect4

    SP5 = TransferFunction(tau*p**s, exp(-tau*p)/s, s)
    expect5_ = TransferFunction(-4*4**s, exp(16)/s, s)
    expect5 = TransferFunction(-4.0*4.0**s, 8886110.52050787/s, s)
    assert str(SP5.subs({tau: -4, p: 4}).evalf()) == str(expect5)
    assert str(expect5_.evalf()) == str(expect5)

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
    assert tf4*tf4 == tf4**2 == pow(tf4, 2) == expect1
    assert tf5*tf5*tf5 == tf5**3 == pow(tf5, 3) == expect2
    assert tf5**0 == pow(tf5, 0) == TransferFunction(1, 1, s)
    assert tf4**-1 == pow(tf4, -1) == TransferFunction(p - 3, p + 4, p)
    assert tf5**-2 == pow(tf5, -2) == TransferFunction((1 - s)**2, (s**2 + 1)**2, s)

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
    expect1 = TransferFunction((s - 5)*(s + 3) + (s + 1)*(s + 6), (s - 5)*(s + 1), s)
    assert tf1 + tf2 == expect1
    assert tf1 + (s - 1) == TransferFunction(s + (s - 5)*(s - 1) + 6, s - 5, s)
    assert tf1 + 8 == TransferFunction(9*s - 34, s - 5, s)
    assert (1 - p**3) + tf1 == TransferFunction(s + (1 - p**3)*(s - 5) + 6, s - 5, s)

    c = symbols("c", commutative=False)
    raises(ValueError, lambda: tf1 + Matrix([1, 2, 3]))
    raises(ValueError, lambda: tf2 + c)
    raises(ValueError, lambda: tf3 + tf4)

    # subtraction
    expect1 = TransferFunction((s + 1)*(s + 6) - (s - 5)*(s + 3), (s - 5)*(s + 1), s)
    assert tf1 - tf2 == expect1
    assert tf3 - tf2 == \
        TransferFunction((s + 1)**2 - (s + 3)*(s**2 + s + 1), (s + 1)*(s**2 + s + 1), s)
    assert tf1 - (s - 1) == TransferFunction(s - (s - 5)*(s - 1) + 6, s - 5, s)
    assert tf1 - 8 == TransferFunction(-7*s + 46, s - 5, s)
    assert (s + 5) - tf2 == TransferFunction(-s + (s + 1)*(s + 5) - 3, s + 1, s)
    assert (1 + p**4) - tf1 == TransferFunction(-s + (p**4 + 1)*(s - 5) - 6, s - 5, s)

    raises(ValueError, lambda: tf1 - Matrix([1, 2, 3]))
    raises(ValueError, lambda: tf3 - tf4)

    p1 = (s - 5)*(s + 1)**2 + ((s + 1)*(s + 6) - (s - 5)*(s + 3))*(s**2 + s + 1)
    p2 = (s - 5)*(s + 1)*(s**2 + s + 1)
    assert tf1 - tf2 + tf3 == TransferFunction(p1, p2, s)


def test_TransferFunction_multiplication_and_division():
    G1 = TransferFunction(s + 3, -s**3 + 9, s)
    G2 = TransferFunction(s + 1, s - 5, s)
    G3 = TransferFunction(p, p**4 - 6, p)
    G4 = TransferFunction(p + 4, p - 5, p)
    G5 = TransferFunction(s + 6, s - 5, s)
    G6 = TransferFunction(s + 3, s + 1, s)

    # multiplication
    expect1 = TransferFunction((s + 1)*(s + 3), (9 - s**3)*(s - 5), s)
    assert G1*G2 == expect1

    expect2 = TransferFunction(p*(p + 4), (p - 5)*(p**4 - 6), p)
    assert G3*G4 == expect2
    assert G5*(s - 1) == TransferFunction((s - 1)*(s + 6), s - 5, s)
    assert 9*G5 == TransferFunction(9*s + 54, s - 5, s)

    c = symbols("c", commutative=False)
    raises(ValueError, lambda: G3 * Matrix([1, 2, 3]))
    raises(ValueError, lambda: G1 * c)
    raises(ValueError, lambda: G3 * G5)

    # division
    assert G5/G6 == TransferFunction((s + 1)*(s + 6), (s - 5)*(s + 3), s)
    assert G5/2 == TransferFunction(s + 6, 2*s - 10, s)
    assert G5/(s**2) == TransferFunction(s + 6, s**2*(s - 5), s)
    assert (s - 4*s**2)/G2 == TransferFunction((s - 5)*(-4*s**2 + s), s + 1, s)
    assert 0/G4 == TransferFunction(0, p + 4, p)

    raises(ValueError, lambda: G3 / Matrix([1, 2, 3]))
    raises(ValueError, lambda: G6 / 0)
    raises(ValueError, lambda: G3 / G5)

    p1 = (9 - s**3)*(s - 5)*(s + 1)*(s + 6) + (s - 5)*(s + 1)*(s + 3)**2
    p2 = (9 - s**3)*(s - 5)**2*(s + 3)
    assert G1*G2 + G5/G6 == TransferFunction(p1, p2, s)

    G3_ = TransferFunction(s**4 - 2*s**2 + 1, s - 1, s)
    p3 = (s - 5)*(s - 1)*(s + 3)
    p4 = (9 - s**3)*(s + 1)*(s**4 - 2*s**2 + 1)
    assert G1/G2/G3_ == TransferFunction(p3, p4, s)

    p5 = (s - 5)*(s - 1)*(s + 3)**2 - (9 - s**3)*(s + 1)**2*(s**4 - 2*s**2 + 1)
    p6 = (9 - s**3)*(s - 1)*(s + 1)*(s + 3)
    assert G1/G2 - G3_/G6 == TransferFunction(p5, p6, s)

    omega_o, zeta, tau = symbols('omega_o, zeta, tau')
    G7 = TransferFunction(a0*s**p + a1*p**s, a2*p*zeta - s, zeta)
    G8 = TransferFunction(omega_o**2, s**2 + p*omega_o*zeta*s + omega_o**2, zeta)
    G9 = TransferFunction(zeta - s**3, zeta + p**4, zeta)
    G10 = TransferFunction(exp(-tau*s*zeta), log(p)**s*s, zeta)
    num = (p**4 + zeta)*(a0*s**p + a1*p**s)*(omega_o**2 + omega_o*p*s*zeta + s**2)*exp(-s*tau*zeta)
    den = s*(omega_o**2*(p**4 + zeta) - (-s**3 + zeta)*(omega_o**2 + omega_o*p*s*zeta + s**2))
    assert G7*G10 / (G8 - G9) == TransferFunction(num, den*(a2*p*zeta - s)*log(p)**s, zeta)


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
    tau = symbols('tau')
    tf1 = TransferFunction(a0*s**5 + log(p)*s + a1, b0*s**4 + b1*s**2 - b2*s, s)
    tf2 = TransferFunction(tau - s**3, tau + p**4, tau)
    tf3 = TransferFunction(a*b*s**3 + s**2 - a*p + s, b - s*p**2, p)
    tf4 = TransferFunction(b*s**2 + p**2 - a*p + s, b - p**2, s)
    assert not tf1.is_biproper
    assert tf2.is_biproper
    assert not tf3.is_biproper
    assert not tf4.is_biproper
