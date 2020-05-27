from sympy import symbols, Matrix, together, simplify, factor
from sympy.physics.control.lti import TransferFunctionMatrix, TransferFunction
from sympy.testing.pytest import raises

a, b, g, s, p = symbols('a, b, g, s, p')


def test_TransferFunctionMatrix_create():
    G = Matrix([1 / s, 1 / (s + 1)])
    tf1 = TransferFunctionMatrix(G)
    assert tf1.G == G

    G = Matrix([[s / (s**2 + 3), 1 / s],
            [1 / (s + g), s / (s**2 - g*s - 3)]])
    tf2 = TransferFunctionMatrix(G)
    assert tf2.G == G

    G = Matrix([[(4*p - 10) / (2*p + 1), 3 / (p + 2)],
            [1 / (2*p**2 + 5*p + 2), (p + 1) / (p**2 + 4*p + 4)]])
    tf3 = TransferFunctionMatrix(G)
    assert tf3.G == G


def test_TransferFunctionMatrix_solve():
    G1 = Matrix([[1 / s + 2 / s**2],
            [2*s / (s**2 + 3*s + 5)]])
    u1 = Matrix([s + 1])
    expect = G1*u1
    diff = Matrix([[0], [0]])
    assert simplify(TransferFunctionMatrix(G1).solve(u1) - expect) == diff

    G2 = Matrix([s / (1 + s**2), 1 / s])
    u2 = Matrix([1 / p])
    expect = G2*u2
    assert simplify(TransferFunctionMatrix(G2).solve(u2) - expect) == diff


def test_TransferFunctionMatrix_series():
    tfm1 = TransferFunctionMatrix(
        Matrix([1 / (a * s)])
    )
    tfm2 = TransferFunctionMatrix(
        Matrix([(b + s) / (a + s)])
    )
    expect = TransferFunctionMatrix(
        Matrix([(b + s) / (a * s * (a + s))])
    )
    assert tfm1.series(tfm2).G == expect.G

    H1 = TransferFunctionMatrix(
        Matrix([7 / (p**9 + 9*p)])
    )
    H2 = TransferFunctionMatrix(
        Matrix([5 / (s + 2)])
    )
    H = TransferFunctionMatrix(
        Matrix([35 / ((s + 2) * (p**9 + 9*p))])
    )
    assert H1.series(H2).G == H.G


def test_TransferFunctionMatrix_parallel():
    tfm1 = TransferFunctionMatrix(
        Matrix([1 / (a * s)])
    )
    tfm2 = TransferFunctionMatrix(
        Matrix([(b + s) / (a + s)])
    )
    expect = TransferFunctionMatrix(
        Matrix([(a + s + a * s * (b + s)) / (a * s * (a + s))])
    )
    assert together(tfm1.parallel(tfm2).G) == together(expect.G)

    H1 = TransferFunctionMatrix(
        Matrix([7 / (p**9 + 9*p)])
    )
    H2 = TransferFunctionMatrix(
        Matrix([5 / (s + 2)])
    )
    H = TransferFunctionMatrix(
        Matrix([7 / (p**9 + 9*p) + 5 / (s + 2)])
    )
    assert together(H1.parallel(H2).G) == together(H.G)


def test_TransferFunctionMatrix_negate():
    G1 = Matrix([s / (s + 1)])
    tf1 = TransferFunctionMatrix(G1)
    assert tf1.neg() == TransferFunctionMatrix(
        Matrix([-s / (s + 1)])
    )

    G2 = Matrix([[(4*p - 10) / (2*p + 1), 3 / (p + 2)],
            [1 / (2*p**2 + 5*p + 2), (p + 1) / (p**2 + 4*p + 4)]])
    tf2 = TransferFunctionMatrix(G2)
    assert tf2.neg() == TransferFunctionMatrix(
        Matrix([[(10 - 4*p) / (2*p + 1),  -3 / (p + 2)],
            [-1 / (2*p**2 + 5*p + 2), (-p - 1) / (p**2 + 4*p + 4)]])
    )


def test_TransferFunctionMatrix_is_proper():
    G1 = Matrix([2*p / (2 + p**2 - p), p / (1 + p**3)])
    assert TransferFunctionMatrix(G1).is_proper

    G2 = Matrix([2*s**6/(2 + s**2 - s), (s**4)/(1 + s**2)])
    assert not TransferFunctionMatrix(G2).is_proper


def test_TransferFunctionMatrix_is_strictly_proper():
    G1 = Matrix([2*s**6 / (2 + s**2 - s), (s**4) / (1 + s**2)])
    assert not TransferFunctionMatrix(G1).is_strictly_proper

    G2 = Matrix([[(4*p - 10) / (2*p**3 + 1), 3 / (p + 2)],
            [1 / (2*p**2 + 5*p + 2), (p + 1) / (p**2 + 4*p + 4)]])
    assert TransferFunctionMatrix(G2).is_strictly_proper


def test_TransferFunction_construction():
    tf = TransferFunction(s + 1, s**2 + s + 1)
    assert tf.num == (s + 1)
    assert tf.den == (s**2 + s + 1)
    assert tf.args == (s + 1, s**2 + s + 1)

    tf1 = TransferFunction(s + 4, s - 5)
    assert tf1.num == (s + 4)
    assert tf1.den == (s - 5)
    assert tf1.args == (s + 4, s - 5)

    # using different polynomial variables.
    tf2 = TransferFunction(p + 3, p**2 - 9)
    assert tf2.num == (p + 3)
    assert tf2.den == (p**2 - 9)
    assert tf2.args == (p + 3, p**2 - 9)

    tf3 = TransferFunction(p**3 + 5*p**2 + 4, p**4 + 3*p + 1)
    assert tf3.args == (p**3 + 5*p**2 + 4, p**4 + 3*p + 1)

    # no pole-zero cancellation on its own.
    tf4 = TransferFunction((s + 3)*(s - 1), (s - 1)*(s + 5))
    assert tf4.den == (s**2 + 4*s - 5)
    assert tf4.args == (s**2 + 2*s - 3, s**2 + 4*s - 5)

    # purely symbolic polynomials.
    tf5 = TransferFunction(p, p)
    assert tf5.num == tf5.den == p
    assert tf5.args == (p, p)

    tf6 = TransferFunction(5, 6)
    assert tf6.num == 5
    assert tf6.den == 6
    assert tf6.args == (5, 6)

    # ValueError when denominator is zero.
    raises(ValueError, lambda: TransferFunction(4, 0))
    raises(ValueError, lambda: TransferFunction(s, 0))
    raises(ValueError, lambda: TransferFunction(0, 0))


def test_TransferFunction_functions():
    # explicitly cancel poles and zeros.
    tf0 = TransferFunction(s**5 + s**3 + s, s - s**2)
    a = TransferFunction(-s**4 - s**2 - 1, s - 1)
    assert tf0.cancel_poles_and_zeros() == a

    tf1 = TransferFunction((p + 3)*(p - 1), (p - 1)*(p + 5))
    b = TransferFunction(p + 3, p + 5)
    assert tf1.cancel_poles_and_zeros() == b

    # negation of TF.
    tf2 = TransferFunction(s + 3, s**2 - s**3 + 9)
    tf3 = TransferFunction(-3*p + 3, 1 - p)
    assert tf2.__neg__() == TransferFunction(-s - 3, s**2 - s**3 + 9)
    assert tf3.__neg__() == TransferFunction(3*p - 3, 1 - p)

    # taking power of a TF.
    tf4 = TransferFunction(p + 4, p - 3)
    tf5 = TransferFunction(s**2 + 1, 1 - s)
    result = TransferFunction(s**6 + 3*s**4 + 3*s**2 + 1, -s**3 + 3*s**2 - 3*s + 1)
    assert tf4.pow(2) == TransferFunction(p**2 + 8*p + 16, p**2 - 6*p + 9)
    assert tf5.pow(3) == result

    # sympy's own functions.
    assert tf3.xreplace({p: s}) == TransferFunction(-3*s + 3, 1 - s)
    tf = TransferFunction(s - 1, s**2 - 2*s + 1)
    assert factor(tf) == TransferFunction(s - 1, (s - 1)**2)


def test_TransferFunction_addition():
    tf1 = TransferFunction(s + 6, s - 5)
    tf2 = TransferFunction(s + 3, s + 1)
    tf3 = TransferFunction(s + 1, s**2 + s + 1)
    tf4 = TransferFunction(p, 2 - p)
    assert tf1.add(tf2) == TransferFunction(2*s**2 + 5*s - 9, s**2 - 4*s - 5)
    assert tf3.add(tf4) == \
        TransferFunction(p*s**2 + 2*s + 2, -p*s**2 - p*s - p**2 + 4*s + 2)
    assert tf1.add(s - 1) == TransferFunction(s**2 - 5*s + 11, s - 5)
    assert tf1.add(8) == TransferFunction(9*s - 34, s - 5)


def test_TransferFunction_multiplication_and_division():
    G1 = TransferFunction(s + 3, -s**3 + 9)
    G2 = TransferFunction(s + 1, s - 5)
    expect = TransferFunction(-s**2 - 4*s - 3, s**4 - 5*s**3 -9*s + 45)
    assert G1.mul(G2) == expect

    G3 = TransferFunction(p, p**4 - 6)
    G4 = TransferFunction(p + 4, p - 5)
    expect = TransferFunction(p**2 + 4*p, p**5 - 5*p**4 - 6*p + 30)
    assert G3.mul(G4) == expect
    assert G3.mul(s) == TransferFunction(p*s, p**4 - 6)
    assert G4.mul(p**-1) == TransferFunction(p + 4, p**2 - 5*p)

    G5 = TransferFunction(s + 6, s - 5)
    G6 = TransferFunction(s + 3, s + 1)
    assert G5.div(G6) == TransferFunction(s**2 + 7*s + 6, s**2 - 2*s - 15)
    assert G5.div(2) == TransferFunction(s + 6, 2*s - 10)


def test_TransferFunction_is_proper():
    G1 = TransferFunction(s**4 + s**3 - 2*s, s**4 - 1)
    assert G1.is_proper

    G2 = TransferFunction(p**7, p - 4)
    assert not G2.is_proper


def test_TransferFunction_is_strictly_proper():
    tf1 = TransferFunction(s**3 - 2, s**4 + 5*s + 6)
    assert tf1.is_strictly_proper

    tf2 = TransferFunction(p + 3, p - 4)
    assert not tf2.is_strictly_proper
