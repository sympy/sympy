from __future__ import division
from sympy import (symbols, Matrix, zeros, simplify, eye, exp,
                   ones, cosh, sinh, Integral, expand, together
                   )
from sympy import Rational as R
from sympy.controlsys.lti import StateSpaceModel, TransferFunctionModel
from mpmath import e


def test_StateSpaceModel_create():
    a1, a2, a3, a4 = symbols('a1:5')
    ssm = StateSpaceModel(Matrix([a1]), Matrix([a2]),
                          Matrix([a3]), Matrix([a4]))
    assert ssm.represent[0] == Matrix([a1])
    assert ssm.represent[1] == Matrix([a2])
    assert ssm.represent[2] == Matrix([a3])
    assert ssm.represent[3] == Matrix([a4])

    b1, b2 = symbols('b1:3')
    c1, c2 = symbols('c1:3')
    d = symbols('d')

    A = Matrix([[a1, a2],
                [a3, a4]])
    B = Matrix([b1,
                b2])
    C = Matrix([[c1, c2]])
    D = Matrix([d])
    ssm = StateSpaceModel(A, B, C, D)
    assert ssm.represent[0] == A
    assert ssm.represent[1] == B
    assert ssm.represent[2] == C
    assert ssm.represent[3] == D

    s = symbols('s')
    G = Matrix([[(4 * s - 10) / (2 * s + 1), 3 / (s + 2)],
                [1 / ((2 * s + 1) * (s + 2)), (s + 1) / (s + 2)**2]])
    ssm = StateSpaceModel(TransferFunctionModel(G))
    assert ssm.represent[0] == Matrix([[-9 / 2, 0, -6, 0, -2, 0],
                                       [0, -9 / 2, 0, -6, 0, -2],
                                       [1, 0, 0, 0, 0, 0],
                                       [0, 1, 0, 0, 0, 0],
                                       [0, 0, 1, 0, 0, 0],
                                       [0, 0, 0, 1, 0, 0]])
    assert ssm.represent[1] == Matrix([[1, 0],
                                       [0, 1],
                                       [0, 0],
                                       [0, 0],
                                       [0, 0],
                                       [0, 0]])
    assert ssm.represent[2] == Matrix([[-6, 3, -24, 15 / 2, -24, 3],
                                       [0, 1, 1 / 2, 3 / 2, 1, 1 / 2]])
    assert ssm.represent[3] == Matrix([[2, 0],
                                       [0, 0]])


def test_TransferFunctionModel_create():
    s = symbols('s')
    G = Matrix([1 / s, 1 / (s + 1)])
    tfm = TransferFunctionModel(G)
    assert tfm.G == G

    g = symbols('g')
    G = Matrix([[s / (s**2 + 3), 1 / s],
                [1 / (s + g), s / (s**2 - g * s - 3)]])
    tfm = TransferFunctionModel(G)
    assert tfm.G == G

    A = Matrix([[-R(9, 2), 0, -6, 0, -2, 0],
                [0, -R(9, 2), 0, -6, 0, -2],
                [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0]])
    B = Matrix([[1, 0],
                [0, 1],
                [0, 0],
                [0, 0],
                [0, 0],
                [0, 0]])
    C = Matrix([[-6, 3, -24, R(15, 2), -24, 3],
                [0, 1, R(1, 2), R(3, 2), 1, R(1, 2)]])
    D = Matrix([[2, 0],
                [0, 0]])
    tfm = TransferFunctionModel(StateSpaceModel(A, B, C, D), s)
    G = Matrix(
        [[(4 * s - 10) / (2 * s + 1), 3 / (s + 2)],
         [1 / (2 * s**2 + 5 * s + 2), (s + 1) / (s**2 + 4 * s + 4)]])
    assert simplify(tfm.G - G) == zeros(2, 2)


def test_StateSpaceModel_eval():
    t, y1, y2, omega = symbols('t, Y1, Y2, omega')
    ssm = StateSpaceModel(eye(2) * omega, zeros(2, 1),
                          eye(2), zeros(2, 1))
    u = eye(1) * omega
    x0 = Matrix([y1, y2])
    assert ssm.evaluate(u, x0, t) == Matrix([y1 * exp(omega * t),
                                             y2 * exp(omega * t)])

    ssm = ssm.subs(omega, 1)
    u = u.subs(omega, 1)
    x0 = x0.subs([(y1, 1.0), (y2, 0.0)])
    sol_num = ssm.evaluate(u, x0, (t, [0, 1.0]))
    assert sol_num == [Matrix([1.0, 0.0]), Matrix([e, 0.0])]

    A = Matrix([[-1, 1], [1, -1]])
    ssm = StateSpaceModel(A, ones(2, 1),
                          eye(2), zeros(2, 1))
    u = Matrix([exp(2 * t)])
    x0 = Matrix([1, 0])
    sol = ssm.evaluate(u, x0, t, do_integrals=True, simplify=True)
    assert sol == Matrix([
                         [cosh(2 * t)],
                         [sinh(2 * t)]])

    sol_num = ssm.evaluate(u, x0, (t, [0]))
    assert sol_num == [Matrix([1.0, 0])]

    sol = ssm.evaluate(u, x0, t, do_integrals=False, simplify=True)
    tau = symbols('tau', positive=True)
    assert sol == Matrix([
        [Integral(exp(2 * tau), (tau, 0, t)) + 1 / 2 + exp(-2 * t) / 2],
        [Integral(exp(2 * tau), (tau, 0, t)) + 1 / 2 - exp(-2 * t) / 2]])

    sol = ssm.evaluate(u, x0, t, 1, do_integrals=False, simplify=True)
    assert sol == Matrix([
        [Integral(exp(2 * tau), (tau, 1, t)) + R(1, 2) + exp(-2 * t + 2) / 2],
        [Integral(exp(2 * tau), (tau, 1, t)) + R(1, 2) - exp(-2 * t + 2) / 2]])


def test_transferFunctionModel_eval():
    s = symbols('s')
    G = Matrix([[1 / s + 2 / s**2],
                [2 * s / (s**2 + 3 * s + 5)]])
    u = Matrix([[s + 1]])
    assert expand(TransferFunctionModel(G).evaluate(u, s) -
                  G * u).is_zero


def test_controllability_matrix():
    a1, a2, a3 = symbols('a1:4')
    A = Matrix([[a1, a2, a3],
                [1, 0, 0],
                [0, 1, 0]])
    B = Matrix([1, 0, 0])
    assert StateSpaceModel(A, B).controllability_matrix() == \
        Matrix([[1, a1, a1**2 + a2],
                [0, 1, a1],
                [0, 0, 1]])

    A = Matrix([[-1, 1, 0],
                [0, -1, 1],
                [0, 0, -1]])
    B = Matrix([[1, 0],
                [1, 0],
                [0, 1]])
    ctrb = StateSpaceModel(A, B).controllability_matrix()
    assert ctrb == Matrix([[1, 0, 0, 0, -1, 1],
                           [1, 0, -1, 1, 1, -2],
                           [0, 1, 0, -1, 0, 1]])


def test_controllable():
    a1, a2, b1, b2 = symbols('a1:3, b1:3')
    A = Matrix([[a1, 0],
                [0, a2]])
    B = Matrix([b1, b2])
    assert StateSpaceModel(A, B).is_controllable() == True

    B = B.subs(b2, 0)
    assert StateSpaceModel(A, B).is_controllable() == False

    A = Matrix([[-5., -2.],
                [6., 2.]])
    B = Matrix([1., -1.])
    assert StateSpaceModel(A, B).is_controllable() == True

    A = Matrix([[-1, 0],
                [0, -3]])
    B = eye(2)
    assert StateSpaceModel(A, B).is_controllable() == True


def test_StateSpaceModel_cascade():
    a0, a1, a2, b0, b1, b2 = symbols('a:3, b:3')
    ssm1 = StateSpaceModel(
        Matrix([[0, 1],
               [-a0, -a1]]),
        Matrix([0, 1]),
        Matrix([[b0, b1]]),
        zeros(1)
    )
    ssm2 = StateSpaceModel(
        Matrix([[0, 1, 0],
                [0, 0, 1],
                [-a0, -a1, -a2]]),
        Matrix([0, 0, 1]),
        Matrix([[b0, b1, b2]]),
        zeros(1)
    )
    expect = StateSpaceModel(
        Matrix([[0, 1, 0, 0, 0],
                [-a0, -a1, 0, 0, 0],
                [0, 0, 0, 1, 0],
                [0, 0, 0, 0, 1],
                [b0, b1, -a0, -a1, -a2]]),
        Matrix([0, 1, 0, 0, 0]),
        Matrix([[0, 0, b0, b1, b2]]),
        zeros(1)
    )
    assert ssm1.cascade(ssm2).represent == expect.represent


def test_TransferFunctionModel_cascade():
    a, b, s = symbols('a, b, s')
    tfm1 = TransferFunctionModel(
        Matrix([1 / (a * s)])
    )
    tfm2 = TransferFunctionModel(
        Matrix([(b + s) / (a + s)])
    )
    expect = TransferFunctionModel(
        Matrix([(b + s) / (a * s * (a + s))])
    )
    assert tfm1.cascade(tfm2).G == expect.G


def test_StateSpaceModel_parallel():
    a0, a1, a2, b0, b1, b2 = symbols('a:3, b:3')
    ssm1 = StateSpaceModel(
        Matrix([[0, 1],
               [-a0, -a1]]),
        Matrix([0, 1]),
        Matrix([[b0, b1]]),
        zeros(1)
    )
    ssm2 = StateSpaceModel(
        Matrix([[0, 1, 0],
                [0, 0, 1],
                [-a0, -a1, -a2]]),
        Matrix([0, 0, 1]),
        Matrix([[b0, b1, b2]]),
        zeros(1)
    )
    expect = StateSpaceModel(
        Matrix([[0, 1, 0, 0, 0],
                [-a0, -a1, 0, 0, 0],
                [0, 0, 0, 1, 0],
                [0, 0, 0, 0, 1],
                [0, 0, -a0, -a1, -a2]]),
        Matrix([0, 1, 0, 0, 1]),
        Matrix([[b0, b1, b0, b1, b2]]),
        zeros(1)
    )
    assert ssm1.parallel(ssm2).represent == expect.represent


def test_TransferFunctionModel_parallel():
    a, b, s = symbols('a, b, s')
    tfm1 = TransferFunctionModel(
        Matrix([1 / (a * s)])
    )
    tfm2 = TransferFunctionModel(
        Matrix([(b + s) / (a + s)])
    )
    expect = TransferFunctionModel(
        Matrix([(a + s + a * s * (b + s)) / (a * s * (a + s))])
    )
    assert together(tfm1.parallel(tfm2).G) == together(expect.G)


def test_equality():
    a0, a1, a2, b0, b1, b2 = symbols('a:3, b:3')
    ssm1 = StateSpaceModel(
        Matrix([[0, 1],
               [-a0, -a1]]),
        Matrix([0, 1]),
        Matrix([[b0, b1]]),
        zeros(1)
    )
    tfm = TransferFunctionModel(ssm1)
    ssm2 = StateSpaceModel(tfm)
    assert tfm == ssm1
    assert tfm == ssm2
    assert ssm1 == ssm2
