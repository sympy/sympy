from sympy import symbols, Matrix, simplify, factor, expand
from sympy.algebras.quaternion import Quaternion
from sympy.physics.control.lti import TransferFunction
from sympy.testing.pytest import raises

a, b, s, p, a0, a1, b0, b1, b2 = symbols('a, b, s, p, a0:2, b0:3')


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

    tf7 = TransferFunction(3*s**2 + 2*p + 4*s, 8*p**2 + 7*s, s)
    tf8 = TransferFunction(3*s**2 + 2*p + 4*s, 8*p**2 + 7*s, p)
    assert not tf7 == tf8

    # when just the numerator is 0, leave the denominator along.
    tf9 = TransferFunction(0, p**2 - p + 1, p)
    assert tf9.args == (0, p**2 - p + 1, p)

    tf10 = TransferFunction(0, 1, s)
    assert tf10.args == (0, 1, s)

    # ValueError when denominator is zero.
    raises(ValueError, lambda: TransferFunction(4, 0, s))
    raises(ValueError, lambda: TransferFunction(s, 0, s))
    raises(ValueError, lambda: TransferFunction(0, 0, s))

    raises(ValueError, lambda: TransferFunction(Matrix([1, 2, 3]), s, s))

    x, y = symbols('x, y', real=True)
    raises(TypeError, lambda: TransferFunction(x**2 + 2*x - 1, x + 3, x))
    raises(TypeError, lambda: TransferFunction(y + 1, 5 - y, y))
    raises(TypeError, lambda: TransferFunction(x - 1, 4 - y, x))


def test_TransferFunction_functions():
    # explicitly cancel poles and zeros.
    tf0 = TransferFunction(s**5 + s**3 + s, s - s**2, s)
    a = TransferFunction(-(s**4 + s**2 + 1), s - 1, s)
    assert tf0.simplify() == simplify(tf0) == a

    tf1 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5), p)
    b = TransferFunction(p + 3, p + 5, p)
    assert tf1.simplify() == simplify(tf1) == b

    # purely symbolic polynomials.
    p1 = a1*s + a0
    p2 = b2*s**2 + b1*s + b0
    SP1 = TransferFunction(p1, p2, s)
    expect = TransferFunction(2*s + 1, 5*s**2 + 4*s + 3, s)
    assert SP1.evalf(subs={a0: 1, a1: 2, b0: 3, b1: 4, b2: 5}) == expect

    c1, d0, d1, d2 = symbols('c1, d0:3')
    p3, p4 = c1*p, d2*p**3 + d1*p**2 - d0
    SP2 = TransferFunction(p3, p4, p)
    expect1 = TransferFunction(2*p, 5*p**3 + 2*p**2 - 3, p)
    assert SP2.evalf(subs={c1: 2, d0: 3, d1: 2, d2: 5}) == expect1

    p5 = a*s + b
    SP3 = TransferFunction(p5, p, s)
    expect2 = TransferFunction(3*s - 4, p, s)
    assert SP3.evalf(subs={a: 3, b: -4}) == expect2

    # negation of TF.
    tf2 = TransferFunction(s + 3, s**2 - s**3 + 9, s)
    tf3 = TransferFunction(-3*p + 3, 1 - p, p)
    assert -tf2 == TransferFunction(-s - 3, s**2 - s**3 + 9, s)
    assert -tf3 == TransferFunction(3*p - 3, 1 - p, p)

    # taking power of a TF.
    tf4 = TransferFunction(p + 4, p - 3, p)
    tf5 = TransferFunction(s**2 + 1, 1 - s, s)
    expect2 = TransferFunction(s**6 + 3*s**4 + 3*s**2 + 1, -s**3 + 3*s**2 - 3*s + 1, s)
    expect1 = TransferFunction(p**2 + 8*p + 16, p**2 - 6*p + 9, p)
    assert tf4*tf4 == tf4**2 == pow(tf4, 2) == expect1
    assert tf5*tf5*tf5 == tf5**3 == pow(tf5, 3) == expect2
    assert tf5**0 == pow(tf5, 0) == TransferFunction(1, 1, s)
    assert tf4**-1 == pow(tf4, -1) == TransferFunction(p - 3, p + 4, p)
    assert tf5**-2 == pow(tf5, -2) == TransferFunction(s**2 - 2*s + 1, s**4 + 2*s**2 + 1, s)

    raises(ValueError, lambda: tf4**(s**2 + s - 1))
    raises(ValueError, lambda: tf5**s)
    raises(ValueError, lambda: tf4**tf5)

    # sympy's own functions.
    assert tf3.xreplace({p: s}) == TransferFunction(-3*s + 3, 1 - s, s)
    tf = TransferFunction(s - 1, s**2 - 2*s + 1, s)
    tf_ = TransferFunction((s - 1)*(s + 3), s + 2, s)
    assert factor(tf) == TransferFunction(s - 1, (s - 1)**2, s)
    assert tf.num.subs(s, 2) == tf.den.subs(s, 2) == 1
    assert tf.subs(s, 2) == TransferFunction(1, 1, s)
    assert expand(tf_) == TransferFunction(s**2 + 2*s - 3, s + 2, s)


def test_TransferFunction_addition_and_subtraction():
    tf1 = TransferFunction(s + 6, s - 5, s)
    tf2 = TransferFunction(s + 3, s + 1, s)
    tf3 = TransferFunction(s + 1, s**2 + s + 1, s)
    tf4 = TransferFunction(p, 2 - p, p)

    # addition
    assert tf1 + tf2 == TransferFunction(2*s**2 + 5*s - 9, s**2 - 4*s - 5, s)
    assert tf1 + (s - 1) == TransferFunction(s**2 - 5*s + 11, s - 5, s)
    assert tf1 + 8 == TransferFunction(9*s - 34, s - 5, s)

    c = symbols("c", commutative=False)
    raises(ValueError, lambda: tf1 + Matrix([1, 2, 3]))
    raises(ValueError, lambda: tf4 + Quaternion(1, 2, 3, 4))
    raises(ValueError, lambda: tf2 + c)
    raises(ValueError, lambda: tf3 + tf4)

    # subtraction
    assert tf1 - tf2 == TransferFunction(9*s + 21, s**2 - 4*s - 5, s)
    assert tf3 - tf2 == \
        TransferFunction(-s**3 - 3*s**2 - 2*s - 2, s**3 + 2*s**2 + 2*s + 1, s)
    assert tf1 - (s - 1) == TransferFunction(-s**2 + 7*s + 1, s - 5, s)
    assert tf1 - 8 == TransferFunction(-7*s + 46, s - 5, s)

    raises(ValueError, lambda: tf1 - Matrix([1, 2, 3]))
    raises(ValueError, lambda: tf4 - Quaternion(3, 4, 5, 6))
    raises(ValueError, lambda: tf2 - c)
    raises(ValueError, lambda: tf3 - tf4)

    assert tf1 - tf2 + tf3 == \
        TransferFunction(10*s**3 + 27*s**2 + 21*s + 16, s**4 - 3*s**3 - 8*s**2 - 9*s - 5, s)


def test_TransferFunction_multiplication_and_division():
    G1 = TransferFunction(s + 3, -s**3 + 9, s)
    G2 = TransferFunction(s + 1, s - 5, s)
    G3 = TransferFunction(p, p**4 - 6, p)
    G4 = TransferFunction(p + 4, p - 5, p)
    G5 = TransferFunction(s + 6, s - 5, s)
    G6 = TransferFunction(s + 3, s + 1, s)

    # multiplication
    expect1 = TransferFunction(s**2 + 4*s + 3, -s**4 + 5*s**3 + 9*s - 45, s)
    assert G1*G2 == expect1

    expect2 = TransferFunction(p**2 + 4*p, p**5 - 5*p**4 - 6*p + 30, p)
    assert G3*G4 == expect2
    assert G5*(s - 1) == TransferFunction(s**2 + 5*s - 6, s - 5, s)
    assert 9*G5 == TransferFunction(9*s + 54, s - 5, s)

    c = symbols("c", commutative=False)
    raises(ValueError, lambda: G3 * Matrix([1, 2, 3]))
    raises(ValueError, lambda: G4 * Quaternion(1, 2, 3, 4))
    raises(ValueError, lambda: G1 * c)
    raises(ValueError, lambda: G3 * G5)

    # division
    assert G5/G6 == TransferFunction(s**2 + 7*s + 6, s**2 - 2*s - 15, s)
    assert G5/2 == TransferFunction(s + 6, 2*s - 10, s)
    assert G5/(s**2) == TransferFunction(s + 6, s**3 - 5*s**2, s)

    raises(ValueError, lambda: G3 / Matrix([1, 2, 3]))
    raises(ValueError, lambda: G4 / Quaternion(1, 2, 3, 4))
    raises(ValueError, lambda: G6 / 0)
    raises(ValueError, lambda: G3 / G5)

    p1 = -s**6 - 2*s**5 + 30*s**4 + 41*s**3 - 2*s**2 - 327*s - 315
    p2 = -s**6 + 7*s**5 + 5*s**4 - 66*s**3 - 63*s**2 - 45*s + 675
    assert G1*G2 + G5/G6 == TransferFunction(p1, p2, s)

    G3_ = TransferFunction(s**4 - 2*s**2 + 1, s - 1, s)
    p3 = s**3 - 3*s**2 - 13*s + 15
    p4 = -s**8 - s**7 + 2*s**6 + 11*s**5 + 8*s**4 - 19*s**3 - 18*s**2 + 9*s + 9
    assert G1/G2/G3_ == TransferFunction(p3, p4, s)

    p5 = s**9 + 2*s**8 - s**7 - 13*s**6 - 19*s**5 + 12*s**4 + 37*s**3 - 13*s**2 - 42*s + 36
    p6 = -s**6 - 3*s**5 + s**4 + 12*s**3 + 27*s**2 - 9*s - 27
    assert G1/G2 - G3_/G6 == TransferFunction(p5, p6, s)


def test_TransferFunction_is_proper():
    G1 = TransferFunction(s**4 + s**3 - 2*s, s**4 - 1, s)
    assert G1.is_proper

    G2 = TransferFunction(p**7, p - 4, p)
    assert not G2.is_proper


def test_TransferFunction_is_strictly_proper():
    tf1 = TransferFunction(s**3 - 2, s**4 + 5*s + 6, s)
    assert tf1.is_strictly_proper

    tf2 = TransferFunction(p + 3, p - 4, p)
    assert not tf2.is_strictly_proper


def test_TransferFunction_is_biproper():
    tf1 = TransferFunction(s - 1, s - 2, s)
    assert tf1.is_biproper

    tf2 = TransferFunction(p**3 - 2, p**4 + 5*p + 6, p)
    assert not tf2.is_biproper
