from sympy import symbols, Matrix, together, simplify
from sympy.physics.control.lti import TransferFunction

a, b, g, s, p = symbols('a, b, g, s, p')


def test_TransferFunction_create():
    G = Matrix([1 / s, 1 / (s + 1)])
    tf1 = TransferFunction(G)
    assert tf1.G == G

    G = Matrix([[s / (s**2 + 3), 1 / s],
            [1 / (s + g), s / (s**2 - g*s - 3)]])
    tf2 = TransferFunction(G)
    assert tf2.G == G

    G = Matrix([[(4*p - 10) / (2*p + 1), 3 / (p + 2)],
            [1 / (2*p**2 + 5*p + 2), (p + 1) / (p**2 + 4*p + 4)]])
    tf3 = TransferFunction(G)
    assert tf3.G == G


def test_TransferFunction_solve():
    G1 = Matrix([[1 / s + 2 / s**2],
            [2*s / (s**2 + 3*s + 5)]])
    u1 = Matrix([s + 1])
    expect = G1*u1
    diff = Matrix([[0], [0]])
    assert simplify(TransferFunction(G1).solve(u1) - expect) == diff

    G2 = Matrix([s / (1 + s**2), 1 / s])
    u2 = Matrix([1 / p])
    expect = G2*u2
    assert simplify(TransferFunction(G2).solve(u2) - expect) == diff


def test_TransferFunction_series():
    tfm1 = TransferFunction(
        Matrix([1 / (a * s)])
    )
    tfm2 = TransferFunction(
        Matrix([(b + s) / (a + s)])
    )
    expect = TransferFunction(
        Matrix([(b + s) / (a * s * (a + s))])
    )
    assert tfm1.series(tfm2).G == expect.G

    H1 = TransferFunction(
        Matrix([7 / (p**9 + 9*p)])
    )
    H2 = TransferFunction(
        Matrix([5 / (s + 2)])
    )
    H = TransferFunction(
        Matrix([35 / ((s + 2) * (p**9 + 9*p))])
    )
    assert H1.series(H2).G == H.G


def test_TransferFunction_parallel():
    tfm1 = TransferFunction(
        Matrix([1 / (a * s)])
    )
    tfm2 = TransferFunction(
        Matrix([(b + s) / (a + s)])
    )
    expect = TransferFunction(
        Matrix([(a + s + a * s * (b + s)) / (a * s * (a + s))])
    )
    assert together(tfm1.parallel(tfm2).G) == together(expect.G)

    H1 = TransferFunction(
        Matrix([7 / (p**9 + 9*p)])
    )
    H2 = TransferFunction(
        Matrix([5 / (s + 2)])
    )
    H = TransferFunction(
        Matrix([7 / (p**9 + 9*p) + 5 / (s + 2)])
    )
    assert together(H1.parallel(H2).G) == together(H.G)


def test_TransferFunction_negate():
    G1 = Matrix([s / (s + 1)])
    tf1 = TransferFunction(G1)
    assert tf1.neg() == TransferFunction(
        Matrix([-s / (s + 1)])
    )

    G2 = Matrix([[(4*p - 10) / (2*p + 1), 3 / (p + 2)],
            [1 / (2*p**2 + 5*p + 2), (p + 1) / (p**2 + 4*p + 4)]])
    tf2 = TransferFunction(G2)
    assert tf2.neg() == TransferFunction(
        Matrix([[(10 - 4*p) / (2*p + 1),  -3 / (p + 2)],
            [-1 / (2*p**2 + 5*p + 2), (-p - 1) / (p**2 + 4*p + 4)]])
    )


def test_TransferFunction_is_proper():
    G1 = Matrix([2*p / (2 + p**2 - p), p / (1 + p**3)])
    assert TransferFunction(G1).is_proper

    G2 = Matrix([2*s**6/(2 + s**2 - s), (s**4)/(1 + s**2)])
    assert not TransferFunction(G2).is_proper


def test_TransferFunction_is_strictly_proper():
    G1 = Matrix([2*s**6/(2 + s**2 - s), (s**4)/(1 + s**2)])
    assert not TransferFunction(G1).is_strictly_proper

    G2 = Matrix([[(4*p - 10) / (2*p**3 + 1), 3 / (p + 2)],
            [1 / (2*p**2 + 5*p + 2), (p + 1) / (p**2 + 4*p + 4)]])
    assert TransferFunction(G2).is_strictly_proper
