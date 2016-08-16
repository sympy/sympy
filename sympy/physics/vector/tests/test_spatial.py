from sympy import cos, Matrix, simplify, sin, symbols, zeros
import sympy.physics.vector.spatial as spatial


def test_rx():
    # Obtain the return results for rx
    theta = symbols('theta')
    E = spatial.rx(theta)

    # Set up the expected return results
    E_exp = Matrix([[1,           0,           0],
                    [0,  cos(theta),  sin(theta)],
                    [0, -sin(theta), cos(theta)]])

    # Test the obtained results against the expected return results
    assert simplify(E - E_exp) == zeros(3)


def test_ry():
    # Obtain the return results for ry
    theta = symbols('theta')
    E = spatial.ry(theta)

    # Set up the expected return results
    E_exp = Matrix([[cos(theta), 0, -sin(theta)],
                    [0,          1,           0],
                    [sin(theta), 0,  cos(theta)]])

    # Test the obtained results against the expected return results
    assert simplify(E - E_exp) == zeros(3)


def test_rz():
    # Obtain the return results for rz
    theta = symbols('theta')
    E = spatial.rz(theta)

    # Set up the expected return results
    E_exp = Matrix([[cos(theta),  sin(theta), 0],
                    [-sin(theta), cos(theta), 0],
                    [0,           0,          1]])

    # Test the obtained rsults against the expected return results
    assert simplify(E - E_exp) == zeros(3)


def test_rot():
    # Obtain the return results for rot for a rotation theta about the z axis
    theta = symbols('theta')
    E = spatial.rz(theta)
    X = spatial.rot(E)

    # Set up the expected return results
    X_exp = Matrix([[cos(theta),  sin(theta), 0, 0,           0,          0],
                    [-sin(theta), cos(theta), 0, 0,           0,          0],
                    [0,           0,          1, 0,           0,          0],
                    [0,           0,          0, cos(theta),  sin(theta), 0],
                    [0,           0,          0, -sin(theta), cos(theta), 0],
                    [0,           0,          0, 0,           0,          1]])

    # Test the obtained results against the expected return results
    assert simplify(X - X_exp) == zeros(6)


def test_xlt():
    # Obtain the return results for xlt
    L = symbols('L')
    r = Matrix([0, 0, L])
    X = spatial.xlt(r)

    # Set up the expected return results
    X_exp = Matrix([[1,  0, 0, 0, 0, 0],
                    [0,  1, 0, 0, 0, 0],
                    [0,  0, 1, 0, 0, 0],
                    [0,  L, 0, 1, 0, 0],
                    [-L, 0, 0, 0, 1, 0],
                    [0,  0, 0, 0, 0, 1]])

    # Test the obtained results against the expected results
    assert simplify(X - X_exp) == zeros(6)


def test_cross_f():
    # Obtain the return results for cross_f
    v1, v2, v3, v4, v5, v6 = symbols('v1 v2 v3 v4 v5 v6')
    v = Matrix([v1, v2, v3, v4, v5, v6])
    crossed_f = spatial.cross_f(v)

    # Set up the expected result
    crossed_f_exp = Matrix([[0,  -v3,  v2,   0, -v6,  v5],
                            [v3,   0, -v1,  v6,   0, -v4],
                            [-v2, v1,   0, -v5,  v4,   0],
                            [0,    0,   0,   0, -v3,  v2],
                            [0,    0,   0,  v3,   0, -v1],
                            [0,    0,   0, -v2,  v1,   0]])

    # Test the obtained results against the expected results
    assert simplify(crossed_f - crossed_f_exp) == zeros(6)


def test_cross_m():
    v1, v2, v3, v4, v5, v6 = symbols('v1 v2 v3 v4 v5 v6')
    # Test the 3x3 return results
    v = Matrix([v1, v2, v3])
    crossed_m = spatial.cross_m(v)

    # Set up the 3x3 expected result
    crossed_m_exp = Matrix([[0,  -v3,  v2],
                            [v3,   0, -v1],
                            [-v2, v1,   0]])

    # Test the obtained 3x3 results against the expected results
    assert simplify(crossed_m - crossed_m_exp) == zeros(3)

    # Set up the expected 3x3 results
    # Obtain the 6x6 return results for cross_m
    v = Matrix([v1, v2, v3, v4, v5, v6])
    crossed_m = spatial.cross_m(v)

    # Set up the 6x6 expected result
    crossed_m_exp = Matrix([[0,  -v3,  v2,  0,   0,   0],
                            [v3,   0, -v1,  0,   0,   0],
                            [-v2, v1,   0,  0,   0,   0],
                            [0,  -v6,  v5,  0, -v3,  v2],
                            [v6,   0, -v4, v3,   0, -v1],
                            [-v5, v4,   0, -v2, v1,   0]])

    # Test the obtained 6x6 results against the expected results
    assert simplify(crossed_m - crossed_m_exp) == zeros(6)


def test_X_to_Xstar():
    # Obtain the return results for X_to_Xstar
    theta = symbols('theta')
    E = spatial.rz(theta)
    X = spatial.rot(E)
    Xstar = spatial.X_to_Xstar(X)

    # Set up the expected return value
    Xstar_exp = Matrix([[-sin(theta)**2/cos(theta) + 1/cos(theta), sin(theta),
                         0, 0, 0, 0],
                        [-sin(theta), cos(theta), 0, 0, 0, 0],
                        [0, 0, 1, 0, 0, 0],
                        [0, 0, 0, -sin(theta)**2/cos(theta) + 1/cos(theta),
                         sin(theta), 0],
                        [0, 0, 0, -sin(theta), cos(theta), 0],
                        [0, 0, 0, 0, 0, 1]])

    # Test the obtained results against the expected results
    assert simplify(Xstar - Xstar_exp) == zeros(6)
