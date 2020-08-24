from sympy import symbols, factor, Function, simplify, exp, pi, ShapeError
from sympy.matrices import Matrix
from sympy.physics.control import TransferFunction, Series, Parallel, \
    Feedback, TransferFunctionMatrix
from sympy.testing.pytest import raises

a, b, s, g, d, p, k, a0, a1, a2, b0, b1, b2, zeta, wn = symbols('a, b, s, g, d, p, k,\
    a0:3, b0:3, zeta, wn')
TF1 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
TF2 = TransferFunction(k, 1, s)
TF3 = TransferFunction(a2*p - s, a2*s + p, s)

def test_TransferFunction_construction():
    tf = TransferFunction(s + 1, s**2 + s + 1, s)
    assert tf.num == (s + 1)
    assert tf.den == (s**2 + s + 1)
    assert tf.args == (s + 1, s**2 + s + 1, s)
    assert tf.shape == (tf.num_outputs, tf.num_inputs) == (1, 1)

    tf1 = TransferFunction(s + 4, s - 5, s)
    assert tf1.num == (s + 4)
    assert tf1.den == (s - 5)
    assert tf1.args == (s + 4, s - 5, s)
    assert tf1.shape == (tf1.num_outputs, tf1.num_inputs) == (1, 1)

    # using different polynomial variables.
    tf2 = TransferFunction(p + 3, p**2 - 9, p)
    assert tf2.num == (p + 3)
    assert tf2.den == (p**2 - 9)
    assert tf2.args == (p + 3, p**2 - 9, p)
    assert tf2.shape == (tf2.num_outputs, tf2.num_inputs) == (1, 1)

    tf3 = TransferFunction(p**3 + 5*p**2 + 4, p**4 + 3*p + 1, p)
    assert tf3.args == (p**3 + 5*p**2 + 4, p**4 + 3*p + 1, p)

    # no pole-zero cancellation on its own.
    tf4 = TransferFunction((s + 3)*(s - 1), (s - 1)*(s + 5), s)
    assert tf4.den == (s - 1)*(s + 5)
    assert tf4.args == ((s + 3)*(s - 1), (s - 1)*(s + 5), s)
    assert tf4.shape == (tf4.num_outputs, tf4.num_inputs) == (1, 1)

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
    assert tf6_.shape == (tf6_.num_outputs, tf6_.num_inputs) == (1, 1)

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
    assert tf10_.shape == tf10.shape == (1, 1)

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
    tfm1 = TransferFunctionMatrix([tf1, tf2, tf3])

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
    raises(ValueError, lambda: tf1 + tfm1)
    raises(ValueError, lambda: tf1 - tf2 + tfm1)

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
    tf = TransferFunction(a0*s**3 + a1*s**2 - a2*s, b0*p**4 + b1*p**3 - b2*s*p, s)
    tf2 = TransferFunction(a2*p - s, a2*s + p, s)
    tf3 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf4 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    tf5 = TransferFunction(a1*p, p + a0, p)
    inp = Function('X_d')(s)
    out = Function('X')(s)

    # SISO transfer function in the arguments.
    s0 = Series(tf, tf2)
    assert s0.args == (tf, tf2)
    assert s0.var == s
    assert s0.is_SISO
    assert s0.shape == (s0.num_outputs, s0.num_inputs) == (1, 1)

    s1 = Series(Parallel(tf, -tf2), tf2)
    assert s1.args == (Parallel(tf, -tf2), tf2)
    assert s1.var == s
    assert s1.is_SISO
    assert s1.shape == (s1.num_outputs, s1.num_inputs) == (1, 1)

    tf3_ = TransferFunction(inp, 1, s)
    tf4_ = TransferFunction(-out, 1, s)
    s2 = Series(tf, Parallel(tf3_, tf4_), tf2)
    assert s2.args == (tf, Parallel(tf3_, tf4_), tf2)
    assert s2.is_SISO
    assert s2.shape == (s2.num_outputs, s2.num_inputs) == (1, 1)

    s3 = Series(tf, tf2, tf4)
    assert s3.args == (tf, tf2, tf4)
    assert s3.is_SISO
    assert s3.shape == (s3.num_outputs, s3.num_inputs) == (1, 1)

    s4 = Series(tf3_, tf4_)
    assert s4.args == (tf3_, tf4_)
    assert s4.var == s
    assert s4.is_SISO
    assert s4.shape == (s4.num_outputs, s4.num_inputs) == (1, 1)

    s6 = Series(tf2, tf4, Parallel(tf2, -tf), tf4)
    assert s6.args == (tf2, tf4, Parallel(tf2, -tf), tf4)
    assert s6.is_SISO
    assert s6.shape == (s6.num_outputs, s6.num_inputs) == (1, 1)

    s7 = Series(tf, tf2)
    assert s0 == s7
    assert not s0 == s2

    raises(ValueError, lambda: Series(tf, tf3))
    raises(ValueError, lambda: Series(tf, tf2, tf3, tf4))
    raises(ValueError, lambda: Series(-tf3, tf2))
    raises(TypeError, lambda: Series(2, tf, tf4))
    raises(TypeError, lambda: Series(s**2 + p*s, tf3, tf2))
    raises(TypeError, lambda: Series(tf3, Matrix([1, 2, 3, 4])))

    # Transfer function matrix in the arguments.
    tfm1 = TransferFunctionMatrix([[tf, tf2, tf4], [-tf4, -tf2, tf]], (2, 3), s)
    tfm2 = TransferFunctionMatrix([-tf, -tf2, -tf4], (3, 1), s)
    tfm3 = TransferFunctionMatrix([-tf4], (1, 1), s)
    tfm4 = TransferFunctionMatrix([[TF3, TF2], [-TF1, tf]])
    tfm5 = TransferFunctionMatrix([TF2, TF1, tf4], (3, 1), s)
    tfm6 = TransferFunctionMatrix([[-TF3, -TF2], [TF1, -tf]])
    tfm7 = TransferFunctionMatrix([[tf3, tf5], [-tf5, -tf3]], (2, 2), p)

    s8 = Series(tfm1, tfm2)
    assert s8.args == (tfm1, tfm2)
    assert s8.var == s
    assert not s8.is_SISO # .is_SISO gives either True or False, not None.
    assert s8.shape == (s8.num_outputs, s8.num_inputs) == (2, 1)

    s9 = Series(tfm1, tfm2, tfm3)
    assert s9.args == (tfm1, tfm2, tfm3)
    assert s9.var == s
    assert not s9.is_SISO
    # (2, 3) x (3, 1) x (1, 1) will give (2, 1)
    assert s9.shape == (s9.num_outputs, s9.num_inputs) == (2, 1)

    s10 = Series(Parallel(tfm2, tfm5), tfm3)
    assert s10.args == (Parallel(tfm2, tfm5), tfm3)
    assert s10.var == s
    assert not s10.is_SISO
    # ((3, 1) + (3, 1)) x (1, 1) will give (3, 1)
    assert s10.shape == (s10.num_outputs, s10.num_inputs) == (3, 1)

    s11 = Series(tfm1, Parallel(-tfm2, -tfm5), tfm3)
    assert s11.args == (tfm1, Parallel(-tfm2, -tfm5), tfm3)
    assert not s11.is_SISO
    # (2, 3) x ((3, 1) + (3, 1)) x (1, 1) will give (2, 1)
    assert s11.shape == (s11.num_outputs, s11.num_inputs) == (2, 1)

    s12 = Series(tfm4, tfm6)
    assert s12.args == (tfm4, tfm6)
    assert not s12.is_SISO
    # (2, 2) x (2, 2) will give (2, 2)
    assert s12.shape == (s12.num_outputs, s12.num_inputs) == (2, 2)

    s13 = Series(tfm4, -tfm1, Parallel(-tfm5, tfm2), -tfm3)
    assert s13.args == (tfm4, -tfm1, Parallel(-tfm5, tfm2), -tfm3)
    # (2, 2) x (2, 3) x ((3, 1) + (3, 1)) x (1, 1) will give (2, 1)
    assert s13.shape == (s13.num_outputs, s13.num_inputs) == (2, 1)

    # for all the adjascent transfer function matrices:
    # no. of inputs of first TFM must be equal to the no. of outputs of the second TFM.
    raises(ValueError, lambda: Series(tfm1, tfm2, -tfm1))

    # all the TFMs must use the same complex variable.
    raises(ValueError, lambda: Series(tfm4, tfm7))

    # Number or expression not allowed in the arguments.
    raises(TypeError, lambda: Series(2, tfm5, tfm3))
    raises(TypeError, lambda: Series(s**2 + p*s, -tfm5, tfm3))


def test_Series_functions():
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)

    tfm1 = TransferFunctionMatrix([[TF1, TF2, TF3], [-TF3, -TF2, TF1]], (2, 3), s)
    tfm2 = TransferFunctionMatrix([-TF1, -TF2, -TF3], (3, 1), s)
    tfm3 = TransferFunctionMatrix([-tf5], (1, 1), s)
    tfm4 = TransferFunctionMatrix([[-TF2, -TF3], [-TF1, TF2]])
    tfm5 = TransferFunctionMatrix([[tf5, -tf5], [-TF3, -TF2]])

    assert TF1*TF2*TF3 == Series(TF1, TF2, TF3)
    assert TF1*(TF2 + TF3) == Series(TF1, Parallel(TF2, TF3))
    assert TF1*TF2 + tf5 == Parallel(Series(TF1, TF2), tf5)
    assert TF1*TF2 - tf5 == Parallel(Series(TF1, TF2), -tf5)
    assert TF1*TF2 + TF3 + tf5 == Parallel(Series(TF1, TF2), TF3, tf5)
    assert TF1*TF2 - TF3 - tf5 == Parallel(Series(TF1, TF2), -TF3, -tf5)
    assert TF1*TF2 - TF3 + tf5 == Parallel(Series(TF1, TF2), -TF3, tf5)
    assert TF1*TF2 + TF3*tf5 == Parallel(Series(TF1, TF2), Series(TF3, tf5))
    assert TF2*TF3*(TF2 - TF1)*TF3 == Series(TF2, TF3, Parallel(TF2, -TF1), TF3)
    assert TF1*TF2 - TF3*tf5 == Parallel(Series(TF1, TF2), Series(TransferFunction(-1, 1, s), Series(TF3, tf5)))
    assert -TF1*TF2 == Series(-TF1, TF2)
    assert -(TF1*TF2) == Series(TransferFunction(-1, 1, s), Series(TF1, TF2))
    raises(ValueError, lambda: TF1*TF2*tf4)
    raises(ValueError, lambda: TF1*(TF2 - tf4))
    raises(ValueError, lambda: TF3*Matrix([1, 2, 3]))

    # SISO transfer function in the arguments.
    # evaluate=True -> doit()
    assert Series(TF1, TF2, evaluate=True) == Series(TF1, TF2).doit() == \
        TransferFunction(k, s**2 + 2*s*wn*zeta + wn**2, s)
    assert Series(TF1, TF2, Parallel(TF1, -TF3), evaluate=True) == Series(TF1, TF2, Parallel(TF1, -TF3)).doit() == \
        TransferFunction(k*(a2*s + p + (-a2*p + s)*(s**2 + 2*s*wn*zeta + wn**2)), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2)**2, s)
    assert Series(TF2, TF1, -TF3, evaluate=True) == Series(TF2, TF1, -TF3).doit() == \
        TransferFunction(k*(-a2*p + s), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert not Series(TF1, -TF2, evaluate=False) == Series(TF1, -TF2).doit()

    assert Series(Parallel(TF1, TF2), Parallel(TF2, -TF3)).doit() == \
        TransferFunction((k*(s**2 + 2*s*wn*zeta + wn**2) + 1)*(-a2*p + k*(a2*s + p) + s), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Series(-TF1, -TF2, -TF3).doit() == \
        TransferFunction(k*(-a2*p + s), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert -Series(TF1, TF2, TF3).doit() == \
        TransferFunction(-k*(a2*p - s), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Series(TF2, TF3, Parallel(TF2, -TF1), TF3).doit() == \
        TransferFunction(k*(a2*p - s)**2*(k*(s**2 + 2*s*wn*zeta + wn**2) - 1), (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2), s)

    assert Series(TF1, TF2).rewrite(TransferFunction) == TransferFunction(k, s**2 + 2*s*wn*zeta + wn**2, s)
    assert Series(TF2, TF1, -TF3).rewrite(TransferFunction) == \
        TransferFunction(k*(-a2*p + s), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    raises(ValueError, lambda: Series(TF1, TF2).rewrite(TransferFunctionMatrix))
    raises(ValueError, lambda: Series(TF2, TF1, -TF3).rewrite(TransferFunctionMatrix))

    S1 = Series(Parallel(TF1, TF2), Parallel(TF2, -TF3))
    assert S1.is_proper
    assert not S1.is_strictly_proper
    assert S1.is_biproper

    S2 = Series(TF1, TF2, TF3)
    assert S2.is_proper
    assert S2.is_strictly_proper
    assert not S2.is_biproper

    S3 = Series(TF1, -TF2, Parallel(TF1, -TF3))
    assert S3.is_proper
    assert S3.is_strictly_proper
    assert not S3.is_biproper

    # Transfer function matrix in the arguments.
    assert Series(tfm1, tfm2, evaluate=True) == Series(tfm1, tfm2).doit() == \
        TransferFunctionMatrix([Parallel(Series(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s)), Series(TransferFunction(k, 1, s), TransferFunction(-k, 1, s)), Series(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a2*p + s, a2*s + p, s))), Parallel(Series(TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s)), Series(TransferFunction(-k, 1, s), TransferFunction(-k, 1, s)), Series(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(-a2*p + s, a2*s + p, s)))])

    assert Series(tfm1, -tfm2, -tfm3, evaluate=True) == Series(tfm1, -tfm2, -tfm3).doit() == \
        TransferFunctionMatrix([Series(Parallel(Series(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s)), Series(TransferFunction(k, 1, s), TransferFunction(k, 1, s)), Series(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(a2*p - s, a2*s + p, s))), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s)), Series(Parallel(Series(TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s)), Series(TransferFunction(-k, 1, s), TransferFunction(k, 1, s)), Series(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(a2*p - s, a2*s + p, s))), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s))])

    assert Series(Parallel(tfm4, -tfm5), tfm1, evaluate=True) == \
        Series(Parallel(tfm4, -tfm5), tfm1).doit() == \
        TransferFunctionMatrix([[Parallel(Series(Parallel(TransferFunction(-k, 1, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s)), Series(Parallel(TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s)), TransferFunction(-a2*p + s, a2*s + p, s))), Parallel(Series(Parallel(TransferFunction(-k, 1, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), TransferFunction(k, 1, s)), Series(Parallel(TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s)), TransferFunction(-k, 1, s))), Parallel(Series(Parallel(TransferFunction(-k, 1, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), TransferFunction(a2*p - s, a2*s + p, s)), Series(Parallel(TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s)), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s)))], [Parallel(Series(Parallel(TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(a2*p - s, a2*s + p, s)), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s)), Series(Parallel(TransferFunction(k, 1, s), TransferFunction(k, 1, s)), TransferFunction(-a2*p + s, a2*s + p, s))), Parallel(Series(Parallel(TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(a2*p - s, a2*s + p, s)), TransferFunction(k, 1, s)), Series(Parallel(TransferFunction(k, 1, s), TransferFunction(k, 1, s)), TransferFunction(-k, 1, s))), Parallel(Series(Parallel(TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(a2*p - s, a2*s + p, s)), TransferFunction(a2*p - s, a2*s + p, s)), Series(Parallel(TransferFunction(k, 1, s), TransferFunction(k, 1, s)), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s)))]])

    assert Series(-tfm4, -tfm5, evaluate=True) == Series(-tfm4, -tfm5).doit() == \
        TransferFunctionMatrix([[Parallel(Series(TransferFunction(k, 1, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), Series(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(a2*p - s, a2*s + p, s))), Parallel(Series(TransferFunction(k, 1, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s)), Series(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(k, 1, s)))], [Parallel(Series(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), Series(TransferFunction(-k, 1, s), TransferFunction(a2*p - s, a2*s + p, s))), Parallel(Series(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s)), Series(TransferFunction(-k, 1, s), TransferFunction(k, 1, s)))]])

    assert Series(Parallel(tfm5, tfm4), Parallel(-tfm4, -tfm5)).doit() == \
        TransferFunctionMatrix([[Parallel(Series(Parallel(TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s), TransferFunction(-k, 1, s)), Parallel(TransferFunction(k, 1, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s))), Series(Parallel(TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s), TransferFunction(-a2*p + s, a2*s + p, s)), Parallel(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(a2*p - s, a2*s + p, s)))), Parallel(Series(Parallel(TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s), TransferFunction(-k, 1, s)), Parallel(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s))), Series(Parallel(TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s), TransferFunction(-a2*p + s, a2*s + p, s)), Parallel(TransferFunction(-k, 1, s), TransferFunction(k, 1, s))))], [Parallel(Series(Parallel(TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s)), Parallel(TransferFunction(k, 1, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s))), Series(Parallel(TransferFunction(-k, 1, s), TransferFunction(k, 1, s)), Parallel(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(a2*p - s, a2*s + p, s)))), Parallel(Series(Parallel(TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s)), Parallel(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s))), Series(Parallel(TransferFunction(-k, 1, s), TransferFunction(k, 1, s)), Parallel(TransferFunction(-k, 1, s), TransferFunction(k, 1, s))))]])

    assert Series(-tfm2, tfm3).rewrite(TransferFunctionMatrix) == \
        TransferFunctionMatrix([Series(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), Series(TransferFunction(k, 1, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), Series(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s))])
    assert Series(Parallel(tfm4, -tfm5), tfm1).doit() == Series(Parallel(tfm4, -tfm5), tfm1).rewrite(TransferFunctionMatrix)
    raises(ValueError, lambda: Series(Parallel(tfm5, tfm4), Parallel(-tfm4, -tfm5)).rewrite(TransferFunction))
    raises(ValueError, lambda: Series(-tfm2, tfm3).rewrite(TransferFunction))

    S4 = Series(Parallel(tfm5, tfm4), Parallel(-tfm4, -tfm5))
    assert not S4.is_proper
    assert not S4.is_strictly_proper
    assert not S4.is_biproper

    S5 = Series(tfm1, tfm2)
    assert S5.is_proper
    assert not S5.is_strictly_proper
    assert S5.is_biproper

    S6 = Series(tfm1, -tfm2, -tfm3)
    assert not S6.is_proper
    assert not S6.is_strictly_proper
    assert not S6.is_biproper


def test_Parallel_construction():
    tf = TransferFunction(a0*s**3 + a1*s**2 - a2*s, b0*p**4 + b1*p**3 - b2*s*p, s)
    tf2 = TransferFunction(a2*p - s, a2*s + p, s)
    tf3 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf4 = TransferFunction(1, s**2 + 2*zeta*wn*s + wn**2, s)
    tf5 = TransferFunction(p, a0 - p, p)
    inp = Function('X_d')(s)
    out = Function('X')(s)

    # SISO transfer function in the arguments.
    p0 = Parallel(tf, tf2)
    assert p0.args == (tf, tf2)
    assert p0.var == s
    assert p0.is_SISO
    assert p0.shape == (p0.num_outputs, p0.num_inputs) == (1, 1)

    p1 = Parallel(Series(tf, -tf2), tf2)
    assert p1.args == (Series(tf, -tf2), tf2)
    assert p1.var == s
    assert p1.is_SISO
    assert p1.shape == (p1.num_outputs, p1.num_inputs) == (1, 1)

    tf3_ = TransferFunction(inp, 1, s)
    tf4_ = TransferFunction(-out, 1, s)
    p2 = Parallel(tf, Series(tf3_, -tf4_), tf2)
    assert p2.args == (tf, Series(tf3_, -tf4_), tf2)
    assert p2.is_SISO
    assert p2.shape == (p2.num_outputs, p2.num_inputs) == (1, 1)

    p3 = Parallel(tf, tf2, tf4)
    assert p3.args == (tf, tf2, tf4)
    assert p3.is_SISO
    assert p3.shape == (p3.num_outputs, p3.num_inputs) == (1, 1)

    p4 = Parallel(tf3_, tf4_)
    assert p4.args == (tf3_, tf4_)
    assert p4.var == s
    assert p4.shape == (p4.num_outputs, p4.num_inputs) == (1, 1)

    p5 = Parallel(tf, tf2)
    assert p0 == p5
    assert not p0 == p1

    p6 = Parallel(tf2, tf4, Series(tf2, -tf4))
    assert p6.args == (tf2, tf4, Series(tf2, -tf4))
    assert p6.shape == (p6.num_outputs, p6.num_inputs) == (1, 1)

    p7 = Parallel(tf2, tf4, Series(tf2, -tf), tf4)
    assert p7.args == (tf2, tf4, Series(tf2, -tf), tf4)
    assert p7.is_SISO
    assert p7.shape == (p7.num_outputs, p7.num_inputs) == (1, 1)

    raises(ValueError, lambda: Parallel(tf, tf3))
    raises(ValueError, lambda: Parallel(tf, tf2, tf3, tf4))
    raises(ValueError, lambda: Parallel(-tf3, tf4))
    raises(TypeError, lambda: Parallel(2, tf, tf4))
    raises(TypeError, lambda: Parallel(s**2 + p*s, tf3, tf2))
    raises(TypeError, lambda: Parallel(tf3, Matrix([1, 2, 3, 4])))

    # Transfer function matrix in the arguments.
    tfm1 = TransferFunctionMatrix([tf, tf2, tf4], (3, 1), s)
    tfm2 = TransferFunctionMatrix((-tf2, tf4, tf), (3, 1), s)
    tfm3 = TransferFunctionMatrix([TF1])
    tfm4 = TransferFunctionMatrix((tf3_, tf4_, TF1), (3, 1), s)
    tfm5 = TransferFunctionMatrix([[tf, tf4], [tf2, tf]])
    tfm6 = TransferFunctionMatrix(((TF1, TF2), (tf, tf2)))
    tfm7 = TransferFunctionMatrix([tf3, tf5, -tf5], (3, 1), p)

    p8 = Parallel(tfm1, tfm2)
    assert p8.args == (tfm1, tfm2)
    assert p8.var == s
    assert not p8.is_SISO # .is_SISO gives either True or False, not None.
    assert p8.shape == (p8.num_outputs, p8.num_inputs) == (3, 1)

    p9 = Parallel(Series(tfm1, tfm3), tfm2)
    assert p9.args == (Series(tfm1, tfm3), tfm2)
    assert p9.var == s
    assert not p9.is_SISO
    assert p9.shape == (p9.num_outputs, p9.num_inputs) == (3, 1)

    p10 = Parallel(tfm1, Series(tfm4, tfm3), tfm2)
    assert p10.args == (tfm1, Series(tfm4, tfm3), tfm2)
    assert p10.var == s
    assert not p10.is_SISO
    assert p10.shape == (p10.num_outputs, p10.num_inputs) == (3, 1)

    p11 = Parallel(tfm2, tfm1, tfm4)
    assert p11.args == (tfm2, tfm1, tfm4)
    assert not p11.is_SISO
    assert p11.shape == (p11.num_outputs, p11.num_inputs) == (3, 1)

    p12 = Parallel(tfm6, tfm5)
    assert p12.args == (tfm6, tfm5)
    assert not p12.is_SISO
    assert p12.shape == (p12.num_outputs, p12.num_inputs) == (2, 2)

    p13 = Parallel(tfm2, tfm4, Series(tfm4, -tfm3), -tfm4)
    assert p13.args == (tfm2, tfm4, Series(tfm4, -tfm3), -tfm4)
    assert p13.shape == (p13.num_outputs, p13.num_inputs) == (3, 1)

    # all TFMs must have same shapes.
    raises(ShapeError, lambda: Parallel(tfm1, tfm3, tfm4))

    # all TFMs must be using the same complex variable.
    raises(ValueError, lambda: Parallel(tfm2, -tfm1, tfm7))

    # Number or expression not allowed in the arguments.
    raises(TypeError, lambda: Parallel(2, tfm1, tfm4))
    raises(TypeError, lambda: Parallel(s**2 + p*s, -tfm4, tfm2))


def test_Parallel_functions():
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)

    tfm1 = TransferFunctionMatrix([TF1, TF2, TF3], (3, 1), s)
    tfm2 = TransferFunctionMatrix((-TF2, tf5, -TF1), (3, 1), s)
    tfm3 = TransferFunctionMatrix((tf5, -tf5, TF2))
    tfm4 = TransferFunctionMatrix([[TF2, -tf5], [TF1, tf5]])
    tfm5 = TransferFunctionMatrix([[TF1, TF2], [TF3, -tf5]])

    assert TF1 + TF2 + TF3 == Parallel(TF1, TF2, TF3)
    assert TF1 + TF2 + TF3 + tf5 == Parallel(TF1, TF2, TF3, tf5)
    assert TF1 + TF2 - TF3 - tf5 == Parallel(TF1, TF2, -TF3, -tf5)
    assert TF1 + TF2*TF3 == Parallel(TF1, Series(TF2, TF3))
    assert TF1 - TF2*TF3 == Parallel(TF1, -Series(TF2,TF3))
    assert -TF1 - TF2 == Parallel(-TF1, -TF2)
    assert -(TF1 + TF2) == Series(TransferFunction(-1, 1, s), Parallel(TF1, TF2))
    assert (TF2 + TF3)*TF1 == Series(Parallel(TF2, TF3), TF1)
    assert (TF1 + TF2)*(TF3*tf5) == Series(Parallel(TF1, TF2), TF3, tf5)
    assert -(TF2 + TF3)*-tf5 == Series(TransferFunction(-1, 1, s), Parallel(TF2, TF3), -tf5)
    assert TF2 + TF3 + TF2*TF1 + tf5 == Parallel(TF2, TF3, Series(TF2, TF1), tf5)
    assert TF2 + TF3 + TF2*TF1 - TF3 == Parallel(TF2, TF3, Series(TF2, TF1), -TF3)
    assert (TF1 + TF2 + tf5)*(TF3 + tf5) == Series(Parallel(TF1, TF2, tf5), Parallel(TF3, tf5))
    raises(ValueError, lambda: TF1 + TF2 + tf4)
    raises(ValueError, lambda: TF1 - TF2*tf4)
    raises(ValueError, lambda: TF3 + Matrix([1, 2, 3]))

    # evaluate=True -> doit()
    assert Parallel(TF1, TF2, evaluate=True) == Parallel(TF1, TF2).doit() == \
        TransferFunction(k*(s**2 + 2*s*wn*zeta + wn**2) + 1, s**2 + 2*s*wn*zeta + wn**2, s)
    assert Parallel(TF1, TF2, Series(-TF1, TF3), evaluate=True) == \
        Parallel(TF1, TF2, Series(-TF1, TF3)).doit()== TransferFunction((-a2*p + s)*(s**2 + 2*s*wn*zeta + wn**2) + \
        (a2*s + p)*(k*(s**2 + 2*s*wn*zeta + wn**2) + 1)*(s**2 + 2*s*wn*zeta + wn**2), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2)**2, s)
    assert Parallel(TF2, TF1, -TF3, evaluate=True) == Parallel(TF2, TF1, -TF3).doit() == \
        TransferFunction(-(a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*s + p)*(k*(s**2 + 2*s*wn*zeta + wn**2) + 1), \
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert not Parallel(TF1, -TF2, evaluate=False) == Parallel(TF1, -TF2).doit()

    assert Parallel(Series(TF1, TF2), Series(TF2, TF3)).doit() == \
        TransferFunction(k*(a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) + k*(a2*s + p), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Parallel(-TF1, -TF2, -TF3).doit() == \
        TransferFunction(-(a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) + \
        (a2*s + p)*(-k*(s**2 + 2*s*wn*zeta + wn**2) - 1), (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert -Parallel(TF1, TF2, TF3).doit() == \
        TransferFunction(-((a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*s + p)*(k*(s**2 + 2*s*wn*zeta + wn**2) + 1)),
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Parallel(TF2, TF3, Series(TF2, -TF1), TF3).doit() == \
        TransferFunction((a2*p - s)*(a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*s + p)*(-k*(a2*s + p) + \
        (s**2 + 2*s*wn*zeta + wn**2)*(a2*p + k*(a2*s + p) - s)), (a2*s + p)**2*(s**2 + 2*s*wn*zeta + wn**2), s)

    assert Parallel(TF1, TF2).rewrite(TransferFunction) == \
        TransferFunction(k*(s**2 + 2*s*wn*zeta + wn**2) + 1, s**2 + 2*s*wn*zeta + wn**2, s)
    assert Parallel(TF2, TF1, -TF3).rewrite(TransferFunction) == \
        TransferFunction(-(a2*p - s)*(s**2 + 2*s*wn*zeta + wn**2) + (a2*s + p)*(k*(s**2 + 2*s*wn*zeta + wn**2) + 1), \
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    raises(ValueError, lambda:  Parallel(TF1, TF2).rewrite(TransferFunctionMatrix))

    P1 = Parallel(Series(TF1, TF2), Series(TF2, TF3))
    assert P1.is_proper
    assert not P1.is_strictly_proper
    assert P1.is_biproper

    P2 = Parallel(TF1, -TF2, -TF3)
    assert P2.is_proper
    assert not P2.is_strictly_proper
    assert P2.is_biproper

    P3 = Parallel(TF1, -TF2, Series(TF1, TF3))
    assert P3.is_proper
    assert not P3.is_strictly_proper
    assert P3.is_biproper

    # transfer function matrix in the arguments.
    assert Parallel(tfm1, tfm2, evaluate=True) == Parallel(tfm1, tfm2).doit() == \
        TransferFunctionMatrix([Parallel(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(-k, 1, s)), Parallel(TransferFunction(k, 1, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s)), Parallel(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s))])

    assert Parallel(tfm1, -tfm2, -tfm3, evaluate=True) == Parallel(tfm1, -tfm2, -tfm3).doit() == \
        TransferFunctionMatrix([Parallel(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(k, 1, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), Parallel(TransferFunction(k, 1, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s)), Parallel(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(-k, 1, s))])

    # both tfm4 and tfm5 have shape (2, 2).
    assert Parallel(tfm4, -tfm5, evaluate=True) == Parallel(tfm4, -tfm5).doit() == \
        TransferFunctionMatrix([[Parallel(TransferFunction(k, 1, s), TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s)), Parallel(TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s), TransferFunction(-k, 1, s))], [Parallel(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(-a2*p + s, a2*s + p, s)), Parallel(TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s))]])

    assert Parallel(tfm4, Series(-tfm4, tfm5), evaluate=True) == Parallel(tfm4, Series(-tfm4, tfm5)).doit() == \
        TransferFunctionMatrix([[Parallel(TransferFunction(k, 1, s), Series(TransferFunction(-k, 1, s), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s)), Series(TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s), TransferFunction(a2*p - s, a2*s + p, s))), Parallel(TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s), Series(TransferFunction(-k, 1, s), TransferFunction(k, 1, s)), Series(TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)))], [Parallel(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), Series(TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s)), Series(TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s), TransferFunction(a2*p - s, a2*s + p, s))), Parallel(TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s), Series(TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(k, 1, s)), Series(TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)))]])

    assert Parallel(Series(tfm4, tfm5), Series(-tfm5, tfm4)).doit() == \
        TransferFunctionMatrix([[Parallel(Series(TransferFunction(k, 1, s), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s)), Series(TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s), TransferFunction(a2*p - s, a2*s + p, s)), Series(TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(k, 1, s)), Series(TransferFunction(-k, 1, s), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s))), Parallel(Series(TransferFunction(k, 1, s), TransferFunction(k, 1, s)), Series(TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), Series(TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), Series(TransferFunction(-k, 1, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s)))], [Parallel(Series(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s)), Series(TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s), TransferFunction(a2*p - s, a2*s + p, s)), Series(TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(k, 1, s)), Series(TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s), TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s))), Parallel(Series(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(k, 1, s)), Series(TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), Series(TransferFunction(-a2*p + s, a2*s + p, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s)), Series(TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s)))]])

    assert Parallel(tfm1, -tfm3, tfm2).rewrite(TransferFunctionMatrix) == \
        TransferFunctionMatrix([Parallel(TransferFunction(1, s**2 + 2*s*wn*zeta + wn**2, s), TransferFunction(a0 - a1*s**2 - a2*s, a0 + s, s), TransferFunction(-k, 1, s)), Parallel(TransferFunction(k, 1, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s), TransferFunction(-a0 + a1*s**2 + a2*s, a0 + s, s)), Parallel(TransferFunction(a2*p - s, a2*s + p, s), TransferFunction(-k, 1, s), TransferFunction(-1, s**2 + 2*s*wn*zeta + wn**2, s))])
    assert Parallel(tfm4, Series(-tfm4, tfm5)).doit() == Parallel(tfm4, Series(-tfm4, tfm5)).rewrite(TransferFunctionMatrix)
    raises(ValueError, lambda: Parallel(tfm4, Series(-tfm4, tfm5)).rewrite(TransferFunction))
    raises(ValueError, lambda: Parallel(Series(tfm4, tfm5), Series(-tfm5, tfm4)).rewrite(TransferFunction))

    P4 = Parallel(tfm1, tfm3, -tfm2)
    assert not P4.is_proper
    assert not P4.is_strictly_proper
    assert not P4.is_biproper

    P5 = Parallel(Series(tfm4, tfm5), tfm4)
    assert not P5.is_proper
    assert not P5.is_strictly_proper
    assert not P5.is_biproper


def test_Feedback_construction():
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)
    tf6 = TransferFunction(s - p, p + s, p)

    f1 = Feedback(TransferFunction(1, 1, s), TF1*TF2*TF3)
    assert f1.args == (TransferFunction(1, 1, s), Series(TF1, TF2, TF3))
    assert f1.num == TransferFunction(1, 1, s)
    assert f1.den == Series(TF1, TF2, TF3)
    assert f1.var == s

    f2 = Feedback(TF1, TF2*TF3)
    assert f2.args == (TF1, Series(TF2, TF3))
    assert f2.num == TF1
    assert f2.den == Series(TF2, TF3)
    assert f2.var == s

    f3 = Feedback(TF1*TF2, tf5)
    assert f3.args == (Series(TF1, TF2), tf5)
    assert f3.num == Series(TF1, TF2)

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
    raises(TypeError, lambda: Feedback(TF1, TF2 + TF3))
    raises(TypeError, lambda: Feedback(TF1, Matrix([1, 2, 3])))
    raises(TypeError, lambda: Feedback(TransferFunction(1, 1, s), s - 1))
    raises(TypeError, lambda: Feedback(1, 1))
    raises(ValueError, lambda: Feedback(TransferFunction(1, 1, s), TransferFunction(1, 1, s)))
    raises(ValueError, lambda: Feedback(TF2, tf4*tf5))


def test_Feedback_functions():
    tf = TransferFunction(1, 1, s)
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)
    tf6 = TransferFunction(s - p, p + s, p)

    assert tf / (tf + TF1) == Feedback(tf, TF1)
    assert tf / (tf + TF1*TF2*TF3) == Feedback(tf, TF1*TF2*TF3)
    assert TF1 / (tf + TF1*TF2*TF3) == Feedback(TF1, TF2*TF3)
    assert (TF1*TF2) / (tf + TF1*TF2) == Feedback(TF1*TF2, tf)
    assert (TF1*TF2) / (tf + TF1*TF2*tf5) == Feedback(TF1*TF2, tf5)
    assert (TF1*TF2) / (tf + TF1*TF2*tf5*TF3) in (Feedback(TF1*TF2, tf5*TF3), Feedback(TF1*TF2, TF3*tf5))
    assert tf4 / (TransferFunction(1, 1, p) + tf4*tf6) == Feedback(tf4, tf6)
    assert tf5 / (tf + tf5) == Feedback(tf5, tf)

    raises(ValueError, lambda: TF1*TF2*TF3 / (1 + TF1*TF2*TF3))
    raises(ValueError, lambda: TF1*TF2*TF3 / TF3*tf5)
    raises(ValueError, lambda: TF2*TF3 / (tf + TF2*TF3*tf4))

    assert Feedback(tf, TF1*TF2*TF3).doit() == \
        TransferFunction((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), k*(a2*p - s) + \
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(TF1, TF2*TF3).doit() == \
        TransferFunction((a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2), (k*(a2*p - s) + \
        (a2*s + p)*(s**2 + 2*s*wn*zeta + wn**2))*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(TF1*TF2, tf5).doit() == \
        TransferFunction(k*(a0 + s)*(s**2 + 2*s*wn*zeta + wn**2), (k*(-a0 + a1*s**2 + a2*s) + \
        (a0 + s)*(s**2 + 2*s*wn*zeta + wn**2))*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(tf4, tf6).doit() == \
        TransferFunction(p*(p + s)*(a0*p + p**a1 - s), p*(p*(p + s) + (-p + s)*(a0*p + p**a1 - s)), p)
    assert -Feedback(tf4*tf6, TransferFunction(1, 1, p)).doit() == \
        TransferFunction(-p*(-p + s)*(p + s)*(a0*p + p**a1 - s), p*(p + s)*(p*(p + s) + (-p + s)*(a0*p + p**a1 - s)), p)

    assert Feedback(TF1, TF2*tf5).rewrite(TransferFunction) == \
        TransferFunction((a0 + s)*(s**2 + 2*s*wn*zeta + wn**2), (k*(-a0 + a1*s**2 + a2*s) + \
        (a0 + s)*(s**2 + 2*s*wn*zeta + wn**2))*(s**2 + 2*s*wn*zeta + wn**2), s)
    assert Feedback(TransferFunction(1, 1, p), tf4).rewrite(TransferFunction) == \
        TransferFunction(p, a0*p + p + p**a1 - s, p)


def test_TransferFunctionMatrix_construction():
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)
    tf4 = TransferFunction(a0*p + p**a1 - s, p, p)
    tf6 = TransferFunction(a0*s**3 + a1*s**2 - a2*s, b0*p**4 + b1*p**3 - b2*s*p, s)
    tf7 = TransferFunction(p, a0, p)

    tfm1 = TransferFunctionMatrix([TF1, TF2], (2, 1), s)
    assert tfm1.shape == (tfm1.num_outputs, tfm1.num_inputs) == (2, 1)
    assert tfm1.args == ([TF1, TF2],)
    assert tfm1.var == s

    tfm1_ = TransferFunctionMatrix((TF2, TF1, TF3), (3, 1), s)
    assert tfm1_.shape == (tfm1_.num_outputs, tfm1_.num_inputs) == (3, 1)
    assert tfm1_.args == ((TF2, TF1, TF3),)
    assert tfm1_.var == s

    tfm2 = TransferFunctionMatrix([-TF1, TF2])
    assert tfm2.shape == (tfm2.num_outputs, tfm2.num_inputs) == (2, 1)
    assert tfm2.args == ([-TF1, TF2],)
    assert tfm2.var == s

    tfm2_ = TransferFunctionMatrix((-TF1, TF2), (2, 1), s)
    assert tfm2_.shape == (tfm2_.num_outputs, tfm2_.num_inputs) == (2, 1)
    assert tfm2_.args == ((-TF1, TF2),)
    assert tfm2_.var == s

    tfm3 = TransferFunctionMatrix([tf7])
    assert tfm3.shape == (tfm3.num_outputs, tfm3.num_inputs) == (1, 1)
    assert tfm3.var == p
    assert tfm3.args == ([tf7],)

    tfm3_ = TransferFunctionMatrix((-TF3,), (1, 1), s)
    assert tfm3_.shape == (tfm3_.num_outputs, tfm3_.num_inputs) == (1, 1)
    assert tfm3_.args == ((-TF3,),)

    tfm4 = TransferFunctionMatrix([TF3, tf5, tf6])
    assert tfm4.shape == (tfm4.num_outputs, tfm4.num_inputs) == (3, 1)
    assert tfm4.args == ([TF3, tf5, tf6],)
    assert tfm4.var == s

    tfm4_ = TransferFunctionMatrix((TF3, tf5, tf6))
    assert tfm4_.shape == (tfm4_.num_outputs, tfm4_.num_inputs) == (3, 1)
    assert tfm4_.args == ((TF3, tf5, tf6),)

    tfm5 = TransferFunctionMatrix([[TF1, -TF2], [TF3, tf5]])
    assert tfm5.shape == (tfm5.num_outputs, tfm5.num_inputs) == (2, 2)
    assert tfm5.args == ([[TF1, -TF2], [TF3, tf5]],)

    tfm5_ = TransferFunctionMatrix(((TF1, -TF2), (TF3, tf5)), (2, 2), s)
    assert tfm5_.shape == (tfm5_.num_outputs, tfm5_.num_inputs) == (2, 2)
    assert tfm5_.args == (((TF1, -TF2), (TF3, tf5)),)
    assert tfm5_.var == s

    tfm6 = TransferFunctionMatrix([[TF1, TF2, TF3], [tf5, tf6, -tf6]], (2, 3), s)
    assert tfm6.shape == (tfm6.num_outputs, tfm6.num_inputs) == (2, 3)
    assert tfm6.args == ([[TF1, TF2, TF3], [tf5, tf6, -tf6]],)

    tfm6_ = TransferFunctionMatrix(((TF1, TF2, TF3), (tf5, tf6, -tf6)))
    assert tfm6_.shape == (tfm6_.num_outputs, tfm6_.num_inputs) == (2, 3)
    assert tfm6_.args == (((TF1, TF2, TF3), (tf5, tf6, -tf6)),)

    tfm7 = TransferFunctionMatrix([[TF1, TF2], [TF3, -tf5], [-tf5, TF2]])
    assert tfm7.shape == (tfm7.num_outputs, tfm7.num_inputs) == (3, 2)
    assert tfm7.args == ([[TF1, TF2], [TF3, -tf5], [-tf5, TF2]],)

    tfm7_ = TransferFunctionMatrix(((TF1, TF2), (TF3, -tf5), (-tf5, TF2)))
    assert tfm7_.shape == (tfm7_.num_outputs, tfm7_.num_inputs) == (3, 2)
    assert tfm7_.args == (((TF1, TF2), (TF3, -tf5), (-tf5, TF2)),)

    # all transfer functions will use the same complex variable. tf4 uses 'p'.
    raises(ValueError, lambda: TransferFunctionMatrix([TF1, TF2, tf4]))
    raises(ValueError, lambda: TransferFunctionMatrix((TF1, TF2, tf4), (3, 1), s))
    raises(ValueError, lambda: TransferFunctionMatrix([[TF1, tf4], [TF3, tf5]]))
    raises(ValueError, lambda: TransferFunctionMatrix(((TF1, tf4), (TF3, tf5))))

    # var provided should be the same complex variable used by all transfer functions.
    raises(ValueError, lambda: TransferFunctionMatrix(((TF1, TF2), (TF3, tf5)), (2, 2), p))
    raises(ValueError, lambda: TransferFunctionMatrix([[TF1, TF2], [TF3, tf5]], (2, 2), p))

    # length of all the lists/tuples in the first arg of TFM should be equal.
    raises(ValueError, lambda: TransferFunctionMatrix([[TF1], [TF3, tf5]]))
    raises(ValueError, lambda: TransferFunctionMatrix(((TF1,), (TF3, tf5))))

    # lists/tuples only support transfer functions in them.
    raises(TypeError, lambda: TransferFunctionMatrix([[TF1, TF2], [TF3, Matrix([1, 2])]]))
    raises(TypeError, lambda: TransferFunctionMatrix(((TF1, TF2), (TF3, Matrix([1, 2]))), (2, 2), s))

    # Shape provided should be equal to (len(args[0]), 1)
    raises(ValueError, lambda: TransferFunctionMatrix((-TF2, -TF1), (4, 1), s))
    # Shape provided should be equal to (len(args[0]), len(args[0][0]))
    raises(ValueError, lambda: TransferFunctionMatrix(((TF1, TF2), (TF3, tf5)), (3, 4), s))


def test_TransferFunctionMatrix_functions():
    tf5 = TransferFunction(a1*s**2 + a2*s - a0, s + a0, s)
    tf6 = TransferFunction(a0*s**3 + a1*s**2 - a2*s, b0*p**4 + b1*p**3 - b2*s*p, s)

    tfm1 = TransferFunctionMatrix([TF1, TF2])
    assert -tfm1 == TransferFunctionMatrix([-TF1, -TF2])

    tfm1_ = TransferFunctionMatrix((-TF1, TF2))
    assert -tfm1_ == TransferFunctionMatrix([TF1, -TF2])

    tfm2 = TransferFunctionMatrix([[TF1, -TF2], [TF3, tf5]])
    assert -tfm2 == TransferFunctionMatrix([[-TF1, TF2], [-TF3, -tf5]])

    tfm2_ = TransferFunctionMatrix(((TF1, -TF2), (TF3, tf5)), (2, 2), s)
    assert -tfm2_ == TransferFunctionMatrix([[-TF1, TF2], [-TF3, -tf5]])

    tfm3 = TransferFunctionMatrix([[TF1, TF2, TF3], [tf5, -TF1, -TF3]])
    assert -tfm3 == TransferFunctionMatrix([[-TF1, -TF2, -TF3], [-tf5, TF1, TF3]])

    tfm3_ = TransferFunctionMatrix(((TF1, TF2, TF3), (tf5, -TF1, -TF3)), (2, 3), s)
    assert -tfm3_ == TransferFunctionMatrix([[-TF1, -TF2, -TF3], [-tf5, TF1, TF3]])

    tfm4 = TransferFunctionMatrix([TF1, TF2, TF3])
    assert tfm4.is_proper
    assert not tfm4.is_strictly_proper
    assert not tfm4.is_biproper

    tfm4_ = TransferFunctionMatrix((TF1, TF2, TF3))
    assert tfm4_.is_proper
    assert not tfm4_.is_strictly_proper
    assert not tfm4_.is_biproper

    tfm5 = TransferFunctionMatrix([[TF1, TF2], [TF3, -tf5], [-tf5, TF2]])
    assert not tfm5.is_proper
    assert not tfm5.is_strictly_proper
    assert not tfm5.is_biproper

    tfm5_ = TransferFunctionMatrix(((TF1, TF2), (TF3, -tf5), (-tf5, TF2)), (3, 2), s)
    assert not tfm5_.is_proper
    assert not tfm5_.is_strictly_proper
    assert not tfm5_.is_biproper

    tfm6 = TransferFunctionMatrix([[TF1, TF2], [TF3, -TF3]])
    assert tfm6.is_proper
    assert not tfm6.is_strictly_proper
    assert not tfm6.is_biproper

    tfm6_ = TransferFunctionMatrix(((TF1, TF2), (TF3, -TF3)))
    assert tfm6_.is_proper
    assert not tfm6_.is_strictly_proper
    assert not tfm6_.is_biproper

    tfm7 = TransferFunctionMatrix([-TF1, -tf6])
    assert not tfm7.is_proper
    assert not tfm7.is_strictly_proper
    assert not tfm7.is_biproper


def test_TransferFunctionMatrix_addition_and_subtraction():
    tf = TransferFunction(a0*p, a1*p**2 + a2*p - a0, p)
    tf_ = TransferFunction(p**3 + a1, p - a2, p)

    tfm1 = TransferFunctionMatrix([TF1, TF2])
    tfm2 = TransferFunctionMatrix([-TF2, -TF1])
    tfm3 = TransferFunctionMatrix([TF1, TF1])
    tfm4 = TransferFunctionMatrix([tf])
    tfm5 = TransferFunctionMatrix([TF2, TF1, TF3])
    tfm6 = TransferFunctionMatrix([-TF1, TF2, -TF3])
    tfm7 = TransferFunctionMatrix([TF2])
    tfm8 = TransferFunctionMatrix([tf, tf_])

    # addition & subtraction.
    assert tfm1 + tfm2 == Parallel(tfm1, tfm2)
    assert tfm3 + tfm1 == Parallel(tfm3, tfm1)
    assert tfm1 - tfm2 == Parallel(tfm1, -tfm2)
    assert tfm3 - tfm1 == Parallel(tfm3, -tfm1)
    assert tfm2 + (tfm3 + tfm1) == Parallel(tfm2, tfm3, tfm1)
    assert tfm2 + (tfm1 - tfm3) == Parallel(tfm2, tfm1, -tfm3)
    assert tfm1 + tfm3 + tfm2 == Parallel(tfm1, tfm3, tfm2)
    assert tfm1 - tfm2 - tfm3 == Parallel(tfm1, -tfm2, -tfm3)
    assert tfm1 + tfm3*tfm7 == Parallel(tfm1, Series(tfm3, tfm7))
    assert tfm7 - tfm3*tfm7 == Parallel(tfm7, -Series(tfm3, tfm7))

    c = symbols("c", commutative=False)
    # Operation with a matrix not supported (for now).
    raises(ValueError, lambda: tfm1 + Matrix([1, 2, 3]))
    raises(ValueError, lambda: tfm1 - Matrix([1, 2, 3]))

    # shape should be equal.
    raises(ValueError, lambda: tfm3 + tfm4)
    raises(ValueError, lambda: tfm3 - tfm4)

    # Operation with a constant, expression or a symbol not allowed.
    raises(ValueError, lambda: tfm2 + c)
    raises(ValueError, lambda: tfm2 - c)
    raises(ValueError, lambda: tfm1 + (s - 1))
    raises(ValueError, lambda: tfm1 - (s - 1))
    raises(ValueError, lambda: (s + 5) - tfm2)
    raises(ValueError, lambda: tfm1 + 8)
    raises(ValueError, lambda: tfm1 - 8)
    raises(ValueError, lambda: (1 - p**3) + tfm1)
    raises(ValueError, lambda: (1 + p**4) - tfm1)

    # Parallel object should have all TFM as arguments.
    raises(ValueError, lambda: tfm1 + (TF1 + TF2))
    raises(ValueError, lambda: tfm2 - (TF2 + TF1))

    # Series object should have all TFM as arguments.
    raises(ValueError, lambda: tfm1 + TF1*TF2)
    raises(ValueError, lambda: tfm1 - TF1*TF2)

    # tfm2 has (2, 1) shape while (tfm5 +/- tfm6) has (3, 1) shape.
    raises(ShapeError, lambda: tfm1 + (tfm5 + tfm6))
    raises(ShapeError, lambda: tfm1 - (tfm6 - tfm5))
    # Both TFM should use the same complex variable.
    raises(ValueError, lambda: tfm1 + tfm8)
    # Both Parallel object and TFM should have the same shape for addition.
    raises(ShapeError, lambda: (tfm5 + tfm6) + tfm1)


def test_TransferFunctionMatrix_multiplication():
    tf = TransferFunction(a0*s**2, a1*s + a2, s)
    tf_ = TransferFunction(a0*p**2, a1 - p, p)
    tfm = TransferFunctionMatrix([tf_, -tf_], (2, 1), p)
    tfm1 = TransferFunctionMatrix([[TF1, TF2], [TF3, tf]])
    tfm2 = TransferFunctionMatrix([[-TF3, -tf], [TF2, TF1]])
    tfm3 = TransferFunctionMatrix([TF1, TF2, TF3])
    tfm5 = TransferFunctionMatrix((-TF2, -TF3))

    assert tfm1*tfm2 == Series(tfm1, tfm2)
    assert -tfm2*tfm1 == Series(-tfm2, tfm1)
    assert tfm1*(tfm2*tfm5) == Series(tfm1, tfm2, tfm5)
    assert -tfm2*(-tfm1*-tfm5) == Series(-tfm2, -tfm1, -tfm5)
    assert tfm1*(tfm2 + (-tfm1)) == Series(tfm1, Parallel(tfm2, -tfm1))
    assert tfm2*(-tfm1 - tfm2) == Series(tfm2, Parallel(-tfm1, -tfm2))

    raises(ValueError, lambda: tfm1 * TF2)

    # Only TFMs are allowed in the parallel/series configurations.
    raises(ValueError, lambda: tfm2 * Parallel(TF1, -TF2))
    raises(ValueError, lambda: tfm2 * Series(tf, TF2))

    # Multiplication of TFM with a Matrix not allowed (for now).
    raises(ValueError, lambda: tfm2 * Matrix([1, 2, 3]))

    # Operation with a constant, Symbol or expression also not supported.
    raises(ValueError, lambda: tfm2 * a0)
    raises(ValueError, lambda: 9 * tfm3)
    raises(ValueError, lambda: tfm2 * (s - 1))

    # tfm2 has (2, 2) shape while tfm3 has (3, 1) shape.
    # No. of inputs of the first TFM must be equal to
    # the no. of outputs of the second TFM.
    raises(ValueError, lambda: tfm2 * tfm3)

    # Both TFM should use the same complex variable.
    raises(ValueError, lambda: tfm1 * tfm)
