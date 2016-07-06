from sympy import cos, Matrix, simplify, sin, symbols, zeros
from sympy.physics.mechanics import dynamicsymbols
import sympy.physics.mechanics.joints as joints


def test_jcalc():
    # Obtain the results for jcalc for the three default joint types (revolute,
    # prismatic and helical)
    the_x, the_y, the_z, x, y, z = dynamicsymbols('the_x the_y the_z x y z')
    p = symbols('p')
    q = [the_x, the_y, the_z, x, y, z]
    [XJ_R, S_R, vJ_R, cJ_R] = joints.jcalc('R', q)
    [XJ_P, S_P, vJ_P, cJ_P] = joints.jcalc('P', q)
    [XJ_H, S_H, vJ_H, cJ_H] = joints.jcalc('H', q, jparam={"pitch": p})

    # Set up the expected results for a revolute joint
    XJ_exp_R = Matrix([[cos(the_z),  sin(the_z), 0, 0,           0,          0],
                       [-sin(the_z), cos(the_z), 0, 0,           0,          0],
                       [0,           0,          1, 0,           0,          0],
                       [0,           0,          0, cos(the_z),  sin(the_z), 0],
                       [0,           0,          0, -sin(the_z), cos(the_z), 0],
                       [0,           0,          0, 0,           0,         1]])

    S_exp_R = Matrix([0, 0, 1, 0, 0, 0])

    vJ_exp_R = Matrix([[the_z.diff()]])

    # Set up the expected results for a prismatic joint
    XJ_exp_P = Matrix([[1,  0, 0, 0, 0, 0],
                       [0,  1, 0, 0, 0, 0],
                       [0,  0, 1, 0, 0, 0],
                       [0,  z, 0, 1, 0, 0],
                       [-z, 0, 0, 0, 1, 0],
                       [0,  0, 0, 0, 0, 1]])

    S_exp_P = Matrix([0, 0, 0, 0, 0, 1])

    vJ_exp_P = Matrix([[z.diff()]])

    # Set up the expected results for a helical joint
    XJ_exp_H = Matrix([[cos(the_z), sin(the_z), 0, 0, 0, 0],
                       [-sin(the_z), cos(the_z), 0, 0, 0, 0],
                       [0, 0, 1, 0, 0, 0],
                       [-p*the_z*sin(the_z), p*the_z*cos(the_z), 0, cos(the_z),
                        sin(the_z), 0],
                       [-p*the_z*cos(the_z), -p*the_z*sin(the_z), 0,
                        -sin(the_z), cos(the_z), 0],
                       [0, 0, 0, 0, 0, 1]])

    S_exp_H = Matrix([0, 0, 1, 0, 0, p])

    vJ_exp_H = Matrix([[p*z.diff() + the_z.diff()]])

    # Test obtained results against expected results
    assert simplify(XJ_R - XJ_exp_R) == zeros(6)
    assert simplify(S_R - S_exp_R) == zeros(6, 1)
    assert simplify(vJ_R - vJ_exp_R) == zeros(1)

    assert simplify(XJ_P - XJ_exp_P) == zeros(6)
    assert simplify(S_P - S_exp_P) == zeros(6, 1)
    assert simplify(vJ_P - vJ_exp_P) == zeros(1)

    assert simplify(XJ_H - XJ_exp_H) == zeros(6)
    assert simplify(S_H - S_exp_H) == zeros(6, 1)
    assert simplify(vJ_H - vJ_exp_H) == zeros(1)


def test_revolute():
    # Obtain the return results for revolute
    theta = symbols('theta')
    [XJ, S, T] = joints.revolute(theta)

    # Set up the expected return results
    XJ_exp = Matrix([[cos(theta),  sin(theta), 0, 0,           0,          0],
                     [-sin(theta), cos(theta), 0, 0,           0,          0],
                     [0,           0,          1, 0,           0,          0],
                     [0,           0,          0, cos(theta),  sin(theta), 0],
                     [0,           0,          0, -sin(theta), cos(theta), 0],
                     [0,           0,          0, 0,           0,          1]])

    S_exp = Matrix([0, 0, 1, 0, 0, 0])

    T_exp = Matrix([[1, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0],
                    [0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 1]])

    # Test the obtained results against the expected results
    assert simplify(XJ - XJ_exp) == zeros(6)
    assert simplify(S - S_exp) == zeros(6, 1)
    assert simplify(T - T_exp) == zeros(6, 5)


def test_prismatic():
    # Obtain the return results for prismatic
    L = symbols('L')
    [XJ, S, T] = joints.prismatic(L)

    # Set up the expected return results
    XJ_exp = Matrix([[1,  0, 0, 0, 0, 0],
                     [0,  1, 0, 0, 0, 0],
                     [0,  0, 1, 0, 0, 0],
                     [0,  L, 0, 1, 0, 0],
                     [-L, 0, 0, 0, 1, 0],
                     [0,  0, 0, 0, 0, 1]])

    S_exp = Matrix([0, 0, 0, 0, 0, 1])

    T_exp = Matrix([[1, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0],
                    [0, 0, 1, 0, 0],
                    [0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0]])

    # Test the obtained results against the expected results
    assert simplify(XJ - XJ_exp) == zeros(6)
    assert simplify(S - S_exp) == zeros(6, 1)
    assert simplify(T - T_exp) == zeros(6, 5)


def test_helical():
    # Obtain the return results for a helical joint
    # Testing three different input options
    #    - Default w/ specifying the pitch
    #    - Defining the pitch with a symbol
    #    - Defining the pitch with a number
    theta, h, p = symbols('theta h p')
    [XJ1, S1, T1] = joints.helical(theta)
    [XJ2, S2, T2] = joints.helical(theta, pitch=p)
    [XJ3, S3, T3] = joints.helical(theta, pitch=5)

    # Set up the expected results for case 1
    XJ_exp1 = Matrix([[cos(theta), sin(theta), 0, 0, 0, 0],
                      [-sin(theta), cos(theta), 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [-h*theta*sin(theta), h*theta*cos(theta), 0, cos(theta),
                       sin(theta), 0],
                      [-h*theta*cos(theta), -h*theta*sin(theta), 0, -sin(theta),
                       cos(theta), 0],
                      [0, 0, 0, 0, 0, 1]])

    S_exp1 = Matrix([0, 0, 1, 0, 0, h])

    T_exp1 = Matrix([[1, 0, 0, 0,  0],
                     [0, 1, 0, 0,  0],
                     [0, 0, 0, 0, -h],
                     [0, 0, 1, 0,  0],
                     [0, 0, 0, 1,  0],
                     [0, 0, 0, 0,  1]])

    # Set up the expected results for case 2
    XJ_exp2 = Matrix([[cos(theta), sin(theta), 0, 0, 0, 0],
                      [-sin(theta), cos(theta), 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [-p*theta*sin(theta), p*theta*cos(theta), 0, cos(theta),
                       sin(theta), 0],
                      [-p*theta*cos(theta), -p*theta*sin(theta), 0, -sin(theta),
                       cos(theta), 0],
                      [0, 0, 0, 0, 0, 1]])

    S_exp2 = Matrix([0, 0, 1, 0, 0, p])

    T_exp2 = Matrix([[1, 0, 0, 0,  0],
                     [0, 1, 0, 0,  0],
                     [0, 0, 0, 0, -p],
                     [0, 0, 1, 0,  0],
                     [0, 0, 0, 1,  0],
                     [0, 0, 0, 0,  1]])

    # Set up the expected results for case 3
    XJ_exp3 = Matrix([[cos(theta), sin(theta), 0, 0, 0, 0],
                      [-sin(theta), cos(theta), 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [-5*theta*sin(theta), 5*theta*cos(theta), 0, cos(theta),
                       sin(theta), 0],
                      [-5*theta*cos(theta), -5*theta*sin(theta), 0, -sin(theta),
                       cos(theta), 0],
                      [0, 0, 0, 0, 0, 1]])

    S_exp3 = Matrix([0, 0, 1, 0, 0, 5])

    T_exp3 = Matrix([[1, 0, 0, 0,  0],
                     [0, 1, 0, 0,  0],
                     [0, 0, 0, 0, -5],
                     [0, 0, 1, 0,  0],
                     [0, 0, 0, 1,  0],
                     [0, 0, 0, 0,  1]])

    # Test the obtained results against the expected results
    assert simplify(XJ1 - XJ_exp1) == zeros(6)
    assert simplify(S1 - S_exp1) == zeros(6, 1)
    assert simplify(T1 - T_exp1) == zeros(6, 5)

    assert simplify(XJ2 - XJ_exp2) == zeros(6)
    assert simplify(S2 - S_exp2) == zeros(6, 1)
    assert simplify(T2 - T_exp2) == zeros(6, 5)

    assert simplify(XJ3 - XJ_exp3) == zeros(6)
    assert simplify(S3 - S_exp3) == zeros(6, 1)
    assert simplify(T3 - T_exp3) == zeros(6, 5)
