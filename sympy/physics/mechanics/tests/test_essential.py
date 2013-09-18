from sympy import cos, Matrix, sin, symbols, pi, S, Function, zeros
from sympy.abc import x, y, z
from sympy.physics.mechanics import Vector, ReferenceFrame, dot, dynamicsymbols
from sympy.physics.mechanics import Dyadic, CoordinateSym, express
from sympy.physics.mechanics.essential import MechanicsLatexPrinter
from sympy.utilities.pytest import raises

Vector.simp = True
A = ReferenceFrame('A')


def test_dyadic():
    d1 = A.x | A.x
    d2 = A.y | A.y
    d3 = A.x | A.y
    assert d1 * 0 == 0
    assert d1 != 0
    assert d1 * 2 == 2 * A.x | A.x
    assert d1 / 2. == 0.5 * d1
    assert d1 & (0 * d1) == 0
    assert d1 & d2 == 0
    assert d1 & A.x == A.x
    assert d1 ^ A.x == 0
    assert d1 ^ A.y == A.x | A.z
    assert d1 ^ A.z == - A.x | A.y
    assert d2 ^ A.x == - A.y | A.z
    assert A.x ^ d1 == 0
    assert A.y ^ d1 == - A.z | A.x
    assert A.z ^ d1 == A.y | A.x
    assert A.x & d1 == A.x
    assert A.y & d1 == 0
    assert A.y & d2 == A.y
    assert d1 & d3 == A.x | A.y
    assert d3 & d1 == 0
    assert d1.dt(A) == 0
    q = dynamicsymbols('q')
    qd = dynamicsymbols('q', 1)
    B = A.orientnew('B', 'Axis', [q, A.z])
    assert d1.express(B) == d1.express(B, B)
    assert d1.express(B) == ((cos(q)**2) * (B.x | B.x) + (-sin(q) * cos(q)) *
            (B.x | B.y) + (-sin(q) * cos(q)) * (B.y | B.x) + (sin(q)**2) *
            (B.y | B.y))
    assert d1.express(B, A) == (cos(q)) * (B.x | A.x) + (-sin(q)) * (B.y | A.x)
    assert d1.express(A, B) == (cos(q)) * (A.x | B.x) + (-sin(q)) * (A.x | B.y)
    assert d1.dt(B) == (-qd) * (A.y | A.x) + (-qd) * (A.x | A.y)


def test_coordinate_vars():
    """Tests the coordinate variables functionality"""
    assert CoordinateSym('Ax', A, 0) == A[0]
    assert CoordinateSym('Ax', A, 1) == A[1]
    assert CoordinateSym('Ax', A, 2) == A[2]
    q = dynamicsymbols('q')
    qd = dynamicsymbols('q', 1)
    assert isinstance(A[0], CoordinateSym) and \
           isinstance(A[0], CoordinateSym) and \
           isinstance(A[0], CoordinateSym)
    assert A.variable_map(A) == {A[0]:A[0], A[1]:A[1], A[2]:A[2]}
    assert A[0].frame == A
    B = A.orientnew('B', 'Axis', [q, A.z])
    assert B.variable_map(A) == {B[2]: A[2], B[1]: -A[0]*sin(q) + A[1]*cos(q),
                                 B[0]: A[0]*cos(q) + A[1]*sin(q)}
    assert A.variable_map(B) == {A[0]: B[0]*cos(q) - B[1]*sin(q),
                                 A[1]: B[0]*sin(q) + B[1]*cos(q), A[2]: B[2]}
    assert A.dt(B[0]) == -A[0]*sin(q)*qd + A[1]*cos(q)*qd
    assert A.dt(B[1]) == -A[0]*cos(q)*qd - A[1]*sin(q)*qd
    assert A.dt(B[2]) == 0
    assert express(B[0], A) == A[0]*cos(q) + A[1]*sin(q)
    assert express(B[1], A) == -A[0]*sin(q) + A[1]*cos(q)
    assert express(B[2], A) == A[2]
    assert B.dt(A[0]*A.x + A[1]*A.y + A[2]*A.z) == A[1]*qd*A.x - A[0]*qd*A.y
    assert A.dt(B[0]*B.x + B[1]*B.y + B[2]*B.z) == - B[1]*qd*B.x + B[0]*qd*B.y
    assert A.express(B[0]*B[1]*B[2]) == \
           A[2]*(-A[0]*sin(q) + A[1]*cos(q))*(A[0]*cos(q) + A[1]*sin(q))
    assert (A.dt(B[0]*B[1]*B[2]) -
            (A[2]*(-A[0]**2*cos(2*q) -
             2*A[0]*A[1]*sin(2*q) +
             A[1]**2*cos(2*q))*qd)).trigsimp() == 0
    assert A.express(B[0]*B.x + B[1]*B.y + B[2]*B.z) == \
           (B[0]*cos(q) - B[1]*sin(q))*A.x + (B[0]*sin(q) + \
           B[1]*cos(q))*A.y + B[2]*A.z
    assert A.express(B[0]*B.x + B[1]*B.y + B[2]*B.z, variables=True) == \
           A[0]*A.x + A[1]*A.y + A[2]*A.z
    assert B.express(A[0]*A.x + A[1]*A.y + A[2]*A.z) == \
           (A[0]*cos(q) + A[1]*sin(q))*B.x + \
           (-A[0]*sin(q) + A[1]*cos(q))*B.y + A[2]*B.z
    assert B.express(A[0]*A.x + A[1]*A.y + A[2]*A.z, variables=True) == \
           B[0]*B.x + B[1]*B.y + B[2]*B.z
    N = B.orientnew('N', 'Axis', [-q, B.z])
    assert N.variable_map(A) == {N[0]: A[0], N[2]: A[2], N[1]: A[1]}
    C = A.orientnew('C', 'Axis', [q, A.x + A.y + A.z])
    mapping = A.variable_map(C)
    assert mapping[A[0]] == 2*C[0]*cos(q)/3 + C[0]/3 - 2*C[1]*sin(q + pi/6)/3 +\
           C[1]/3 - 2*C[2]*cos(q + pi/3)/3 + C[2]/3
    assert mapping[A[1]] == -2*C[0]*cos(q + pi/3)/3 + \
           C[0]/3 + 2*C[1]*cos(q)/3 + C[1]/3 - 2*C[2]*sin(q + pi/6)/3 + C[2]/3
    assert mapping[A[2]] == -2*C[0]*sin(q + pi/6)/3 + C[0]/3 - \
           2*C[1]*cos(q + pi/3)/3 + C[1]/3 + 2*C[2]*cos(q)/3 + C[2]/3


def test_ang_vel():
    q1, q2, q3, q4 = dynamicsymbols('q1 q2 q3 q4')
    q1d, q2d, q3d, q4d = dynamicsymbols('q1 q2 q3 q4', 1)
    N = ReferenceFrame('N')
    A = N.orientnew('A', 'Axis', [q1, N.z])
    B = A.orientnew('B', 'Axis', [q2, A.x])
    C = B.orientnew('C', 'Axis', [q3, B.y])
    D = N.orientnew('D', 'Axis', [q4, N.y])
    u1, u2, u3 = dynamicsymbols('u1 u2 u3')
    assert A.ang_vel_in(N) == (q1d)*A.z
    assert B.ang_vel_in(N) == (q2d)*B.x + (q1d)*A.z
    assert C.ang_vel_in(N) == (q3d)*C.y + (q2d)*B.x + (q1d)*A.z

    A2 = N.orientnew('A2', 'Axis', [q4, N.y])
    assert N.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == -q1d*N.z
    assert N.ang_vel_in(B) == -q1d*A.z - q2d*B.x
    assert N.ang_vel_in(C) == -q1d*A.z - q2d*B.x - q3d*B.y
    assert N.ang_vel_in(A2) == -q4d*N.y

    assert A.ang_vel_in(N) == q1d*N.z
    assert A.ang_vel_in(A) == 0
    assert A.ang_vel_in(B) == - q2d*B.x
    assert A.ang_vel_in(C) == - q2d*B.x - q3d*B.y
    assert A.ang_vel_in(A2) == q1d*N.z - q4d*N.y

    assert B.ang_vel_in(N) == q1d*A.z + q2d*A.x
    assert B.ang_vel_in(A) == q2d*A.x
    assert B.ang_vel_in(B) == 0
    assert B.ang_vel_in(C) == -q3d*B.y
    assert B.ang_vel_in(A2) == q1d*A.z + q2d*A.x - q4d*N.y

    assert C.ang_vel_in(N) == q1d*A.z + q2d*A.x + q3d*B.y
    assert C.ang_vel_in(A) == q2d*A.x + q3d*C.y
    assert C.ang_vel_in(B) == q3d*B.y
    assert C.ang_vel_in(C) == 0
    assert C.ang_vel_in(A2) == q1d*A.z + q2d*A.x + q3d*B.y - q4d*N.y

    assert A2.ang_vel_in(N) == q4d*A2.y
    assert A2.ang_vel_in(A) == q4d*A2.y - q1d*N.z
    assert A2.ang_vel_in(B) == q4d*N.y - q1d*A.z - q2d*A.x
    assert A2.ang_vel_in(C) == q4d*N.y - q1d*A.z - q2d*A.x - q3d*B.y
    assert A2.ang_vel_in(A2) == 0

    C.set_ang_vel(N, u1*C.x + u2*C.y + u3*C.z)
    assert C.ang_vel_in(N) == (u1)*C.x + (u2)*C.y + (u3)*C.z
    assert N.ang_vel_in(C) == (-u1)*C.x + (-u2)*C.y + (-u3)*C.z
    assert C.ang_vel_in(D) == (u1)*C.x + (u2)*C.y + (u3)*C.z + (-q4d)*D.y
    assert D.ang_vel_in(C) == (-u1)*C.x + (-u2)*C.y + (-u3)*C.z + (q4d)*D.y

    q0 = dynamicsymbols('q0')
    q0d = dynamicsymbols('q0', 1)
    E = N.orientnew('E', 'Quaternion', (q0, q1, q2, q3))
    assert E.ang_vel_in(N) == (
        2 * (q1d * q0 + q2d * q3 - q3d * q2 - q0d * q1) * E.x +
        2 * (q2d * q0 + q3d * q1 - q1d * q3 - q0d * q2) * E.y +
        2 * (q3d * q0 + q1d * q2 - q2d * q1 - q0d * q3) * E.z)

    F = N.orientnew('F', 'Body', (q1, q2, q3), '313')
    assert F.ang_vel_in(N) == ((sin(q2)*sin(q3)*q1d + cos(q3)*q2d)*F.x +
        (sin(q2)*cos(q3)*q1d - sin(q3)*q2d)*F.y + (cos(q2)*q1d + q3d)*F.z)
    G = N.orientnew('G', 'Axis', (q1, N.x + N.y))
    assert G.ang_vel_in(N) == q1d * (N.x + N.y).normalize()
    assert N.ang_vel_in(G) == -q1d * (N.x + N.y).normalize()


def test_dcm():
    q1, q2, q3, q4 = dynamicsymbols('q1 q2 q3 q4')
    N = ReferenceFrame('N')
    A = N.orientnew('A', 'Axis', [q1, N.z])
    B = A.orientnew('B', 'Axis', [q2, A.x])
    C = B.orientnew('C', 'Axis', [q3, B.y])
    D = N.orientnew('D', 'Axis', [q4, N.y])
    E = N.orientnew('E', 'Space', [q1, q2, q3], '123')
    assert N.dcm(C) == Matrix([
        [- sin(q1) * sin(q2) * sin(q3) + cos(q1) * cos(q3), - sin(q1) *
        cos(q2), sin(q1) * sin(q2) * cos(q3) + sin(q3) * cos(q1)], [sin(q1) *
        cos(q3) + sin(q2) * sin(q3) * cos(q1), cos(q1) * cos(q2), sin(q1) *
            sin(q3) - sin(q2) * cos(q1) * cos(q3)], [- sin(q3) * cos(q2), sin(q2),
        cos(q2) * cos(q3)]])
    # This is a little touchy.  Is it ok to use simplify in assert?
    test_mat = D.dcm(C) - Matrix(
        [[cos(q1) * cos(q3) * cos(q4) - sin(q3) * (- sin(q4) * cos(q2) +
        sin(q1) * sin(q2) * cos(q4)), - sin(q2) * sin(q4) - sin(q1) *
            cos(q2) * cos(q4), sin(q3) * cos(q1) * cos(q4) + cos(q3) * (- sin(q4) *
        cos(q2) + sin(q1) * sin(q2) * cos(q4))], [sin(q1) * cos(q3) +
        sin(q2) * sin(q3) * cos(q1), cos(q1) * cos(q2), sin(q1) * sin(q3) -
            sin(q2) * cos(q1) * cos(q3)], [sin(q4) * cos(q1) * cos(q3) -
        sin(q3) * (cos(q2) * cos(q4) + sin(q1) * sin(q2) * sin(q4)), sin(q2) *
                cos(q4) - sin(q1) * sin(q4) * cos(q2), sin(q3) * sin(q4) * cos(q1) +
                cos(q3) * (cos(q2) * cos(q4) + sin(q1) * sin(q2) * sin(q4))]])
    assert test_mat.expand() == zeros(3, 3)
    assert E.dcm(N) == Matrix(
        [[cos(q2)*cos(q3), sin(q3)*cos(q2), -sin(q2)],
        [sin(q1)*sin(q2)*cos(q3) - sin(q3)*cos(q1), sin(q1)*sin(q2)*sin(q3) +
        cos(q1)*cos(q3), sin(q1)*cos(q2)], [sin(q1)*sin(q3) +
        sin(q2)*cos(q1)*cos(q3), - sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1),
         cos(q1)*cos(q2)]])


def test_Vector():
    assert A.x != A.y
    assert A.y != A.z
    assert A.z != A.x

    v1 = x*A.x + y*A.y + z*A.z
    v2 = x**2*A.x + y**2*A.y + z**2*A.z
    v3 = v1 + v2
    v4 = v1 - v2

    assert isinstance(v1, Vector)
    assert dot(v1, A.x) == x
    assert dot(v1, A.y) == y
    assert dot(v1, A.z) == z

    assert isinstance(v2, Vector)
    assert dot(v2, A.x) == x**2
    assert dot(v2, A.y) == y**2
    assert dot(v2, A.z) == z**2

    assert isinstance(v3, Vector)
    # We probably shouldn't be using simplify in dot...
    assert dot(v3, A.x) == x**2 + x
    assert dot(v3, A.y) == y**2 + y
    assert dot(v3, A.z) == z**2 + z

    assert isinstance(v4, Vector)
    # We probably shouldn't be using simplify in dot...
    assert dot(v4, A.x) == x - x**2
    assert dot(v4, A.y) == y - y**2
    assert dot(v4, A.z) == z - z**2


def test_Vector_diffs():
    q1, q2, q3, q4 = dynamicsymbols('q1 q2 q3 q4')
    q1d, q2d, q3d, q4d = dynamicsymbols('q1 q2 q3 q4', 1)
    q1dd, q2dd, q3dd, q4dd = dynamicsymbols('q1 q2 q3 q4', 2)
    N = ReferenceFrame('N')
    A = N.orientnew('A', 'Axis', [q3, N.z])
    B = A.orientnew('B', 'Axis', [q2, A.x])
    v1 = q2 * A.x + q3 * N.y
    v2 = q3 * B.x + v1
    v3 = v1.dt(B)
    v4 = v2.dt(B)
    v5 = q1*A.x + q2*A.y + q3*A.z

    assert v1.dt(N) == N.dt(v1) == q2d * A.x + q2 * q3d * A.y + q3d * N.y
    assert v1.dt(A) == A.dt(v1) == q2d * A.x + q3 * q3d * N.x + q3d * N.y
    assert v1.dt(B) == B.dt(v1) == (q2d * A.x + q3 * q3d * N.x + q3d *
                                    N.y - q3 * cos(q3) * q2d * N.z)
    assert v2.dt(N) == (q2d * A.x + (q2 + q3) * q3d * A.y + q3d * B.x + q3d *
                        N.y)
    assert v2.dt(A) == q2d * A.x + q3d * B.x + q3 * q3d * N.x + q3d * N.y
    assert v2.dt(B) == (q2d * A.x + q3d * B.x + q3 * q3d * N.x + q3d * N.y -
                        q3 * cos(q3) * q2d * N.z)
    assert v3.dt(N) == (q2dd * A.x + q2d * q3d * A.y + (q3d**2 + q3 * q3dd) *
                        N.x + q3dd * N.y + (q3 * sin(q3) * q2d * q3d -
                        cos(q3) * q2d * q3d - q3 * cos(q3) * q2dd) * N.z)
    assert v3.dt(A) == (q2dd * A.x + (2 * q3d**2 + q3 * q3dd) * N.x + (q3dd -
                        q3 * q3d**2) * N.y + (q3 * sin(q3) * q2d * q3d -
                        cos(q3) * q2d * q3d - q3 * cos(q3) * q2dd) * N.z)
    assert v3.dt(B) == (q2dd * A.x - q3 * cos(q3) * q2d**2 * A.y + (2 *
                        q3d**2 + q3 * q3dd) * N.x + (q3dd - q3 * q3d**2) *
                        N.y + (2 * q3 * sin(q3) * q2d * q3d - 2 * cos(q3) *
                        q2d * q3d - q3 * cos(q3) * q2dd) * N.z)
    assert v4.dt(N) == (q2dd * A.x + q3d * (q2d + q3d) * A.y + q3dd * B.x +
                        (q3d**2 + q3 * q3dd) * N.x + q3dd * N.y + (q3 *
                        sin(q3) * q2d * q3d - cos(q3) * q2d * q3d - q3 *
                        cos(q3) * q2dd) * N.z)
    assert v4.dt(A) == (q2dd * A.x + q3dd * B.x + (2 * q3d**2 + q3 * q3dd) *
                        N.x + (q3dd - q3 * q3d**2) * N.y + (q3 * sin(q3) *
                        q2d * q3d - cos(q3) * q2d * q3d - q3 * cos(q3) *
                        q2dd) * N.z)
    assert v4.dt(B) == (q2dd * A.x - q3 * cos(q3) * q2d**2 * A.y + q3dd * B.x +
                        (2 * q3d**2 + q3 * q3dd) * N.x + (q3dd - q3 * q3d**2) *
                        N.y + (2 * q3 * sin(q3) * q2d * q3d - 2 * cos(q3) *
                        q2d * q3d - q3 * cos(q3) * q2dd) * N.z)
    assert B.dt(v5) == v5.dt(B) == q1d*A.x + (q3*q2d + q2d)*A.y + (-q2*q2d + q3d)*A.z
    assert A.dt(v5) == v5.dt(A) == q1d*A.x + q2d*A.y + q3d*A.z
    assert N.dt(v5) == v5.dt(N) == (-q2*q3d + q1d)*A.x + (q1*q3d + q2d)*A.y + q3d*A.z
    assert v3.diff(q1d, N) == 0
    assert v3.diff(q2d, N) == A.x - q3 * cos(q3) * N.z
    assert v3.diff(q3d, N) == q3 * N.x + N.y
    assert v3.diff(q1d, A) == 0
    assert v3.diff(q2d, A) == A.x - q3 * cos(q3) * N.z
    assert v3.diff(q3d, A) == q3 * N.x + N.y
    assert v3.diff(q1d, B) == 0
    assert v3.diff(q2d, B) == A.x - q3 * cos(q3) * N.z
    assert v3.diff(q3d, B) == q3 * N.x + N.y
    assert v4.diff(q1d, N) == 0
    assert v4.diff(q2d, N) == A.x - q3 * cos(q3) * N.z
    assert v4.diff(q3d, N) == B.x + q3 * N.x + N.y
    assert v4.diff(q1d, A) == 0
    assert v4.diff(q2d, A) == A.x - q3 * cos(q3) * N.z
    assert v4.diff(q3d, A) == B.x + q3 * N.x + N.y
    assert v4.diff(q1d, B) == 0
    assert v4.diff(q2d, B) == A.x - q3 * cos(q3) * N.z
    assert v4.diff(q3d, B) == B.x + q3 * N.x + N.y


def test_vector_simplify():
    x, y, z, k, n, m, w, f, s, A = symbols('x, y, z, k, n, m, w, f, s, A')
    N = ReferenceFrame('N')

    test1 = (1 / x + 1 / y) * N.x
    assert (test1 & N.x) != (x + y) / (x * y)
    test1 = test1.simplify()
    assert (test1 & N.x) == (x + y) / (x * y)

    test2 = (A**2 * s**4 / (4 * pi * k * m**3)) * N.x
    test2 = test2.simplify()
    assert (test2 & N.x) == (A**2 * s**4 / (4 * pi * k * m**3))

    test3 = ((4 + 4 * x - 2 * (2 + 2 * x)) / (2 + 2 * x)) * N.x
    test3 = test3.simplify()
    assert (test3 & N.x) == 0

    test4 = ((-4 * x * y**2 - 2 * y**3 - 2 * x**2 * y) / (x + y)**2) * N.x
    test4 = test4.simplify()
    assert (test4 & N.x) == -2 * y


def test_dyadic_simplify():
    x, y, z, k, n, m, w, f, s, A = symbols('x, y, z, k, n, m, w, f, s, A')
    N = ReferenceFrame('N')

    dy = N.x | N.x
    test1 = (1 / x + 1 / y) * dy
    assert (N.x & test1 & N.x) != (x + y) / (x * y)
    test1 = test1.simplify()
    assert (N.x & test1 & N.x) == (x + y) / (x * y)

    test2 = (A**2 * s**4 / (4 * pi * k * m**3)) * dy
    test2 = test2.simplify()
    assert (N.x & test2 & N.x) == (A**2 * s**4 / (4 * pi * k * m**3))

    test3 = ((4 + 4 * x - 2 * (2 + 2 * x)) / (2 + 2 * x)) * dy
    test3 = test3.simplify()
    assert (N.x & test3 & N.x) == 0

    test4 = ((-4 * x * y**2 - 2 * y**3 - 2 * x**2 * y) / (x + y)**2) * dy
    test4 = test4.simplify()
    assert (N.x & test4 & N.x) == -2 * y


def test_latex_printer():
    r = Function('r')('t')
    assert MechanicsLatexPrinter().doprint(r**2) == "r^{2}"


def test_output_type():
    A = ReferenceFrame('A')
    v = A.x + A.y
    d = v | v
    zerov = Vector(0)
    zerod = Dyadic(0)

    # dot products
    assert isinstance(d & d, Dyadic)
    assert isinstance(d & zerod, Dyadic)
    assert isinstance(zerod & d, Dyadic)
    assert isinstance(d & v, Vector)
    assert isinstance(v & d, Vector)
    assert isinstance(d & zerov, Vector)
    assert isinstance(zerov & d, Vector)
    raises(TypeError, lambda: d & S(0))
    raises(TypeError, lambda: S(0) & d)
    raises(TypeError, lambda: d & 0)
    raises(TypeError, lambda: 0 & d)
    assert not isinstance(v & v, (Vector, Dyadic))
    assert not isinstance(v & zerov, (Vector, Dyadic))
    assert not isinstance(zerov & v, (Vector, Dyadic))
    raises(TypeError, lambda: v & S(0))
    raises(TypeError, lambda: S(0) & v)
    raises(TypeError, lambda: v & 0)
    raises(TypeError, lambda: 0 & v)

    # cross products
    raises(TypeError, lambda: d ^ d)
    raises(TypeError, lambda: d ^ zerod)
    raises(TypeError, lambda: zerod ^ d)
    assert isinstance(d ^ v, Dyadic)
    assert isinstance(v ^ d, Dyadic)
    assert isinstance(d ^ zerov, Dyadic)
    assert isinstance(zerov ^ d, Dyadic)
    assert isinstance(zerov ^ d, Dyadic)
    raises(TypeError, lambda: d ^ S(0))
    raises(TypeError, lambda: S(0) ^ d)
    raises(TypeError, lambda: d ^ 0)
    raises(TypeError, lambda: 0 ^ d)
    assert isinstance(v ^ v, Vector)
    assert isinstance(v ^ zerov, Vector)
    assert isinstance(zerov ^ v, Vector)
    raises(TypeError, lambda: v ^ S(0))
    raises(TypeError, lambda: S(0) ^ v)
    raises(TypeError, lambda: v ^ 0)
    raises(TypeError, lambda: 0 ^ v)

    # outer products
    raises(TypeError, lambda: d | d)
    raises(TypeError, lambda: d | zerod)
    raises(TypeError, lambda: zerod | d)
    raises(TypeError, lambda: d | v)
    raises(TypeError, lambda: v | d)
    raises(TypeError, lambda: d | zerov)
    raises(TypeError, lambda: zerov | d)
    raises(TypeError, lambda: zerov | d)
    raises(TypeError, lambda: d | S(0))
    raises(TypeError, lambda: S(0) | d)
    raises(TypeError, lambda: d | 0)
    raises(TypeError, lambda: 0 | d)
    assert isinstance(v | v, Dyadic)
    assert isinstance(v | zerov, Dyadic)
    assert isinstance(zerov | v, Dyadic)
    raises(TypeError, lambda: v | S(0))
    raises(TypeError, lambda: S(0) | v)
    raises(TypeError, lambda: v | 0)
    raises(TypeError, lambda: 0 | v)
