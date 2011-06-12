from sympy import symbols, Matrix, sin, cos
from sympy.physics.classical import (Vector, ReferenceFrame, dot, cross,
                                     dynamicsymbols)

phi, x, y, z = symbols('phi x y z')
A = ReferenceFrame('A')

def test_ang_vel():
    q1, q2, q3, q4, q1d, q2d, q3d, q4d = dynamicsymbols('q1 q2 q3 q4 q1d q2d \
                                                        q3d q4d')
    N = ReferenceFrame('N')
    A = N.orientnew('A', 'Simple', q1, 3)
    B = A.orientnew('B', 'Simple', q2, 1)
    C = B.orientnew('C', 'Simple', q3, 2)
    D = N.orientnew('D', 'Simple', q4, 2)
    u0, u1, u2 = dynamicsymbols('u1 u2 u3')
    assert A.ang_vel_in(N) == (q1d)*A.z
    assert B.ang_vel_in(N) == (q2d)*B.x + (q1d)*A.z
    assert C.ang_vel_in(N) == (q3d)*C.y + (q2d)*B.x + (q1d)*A.z

    A2 = N.orientnew('A2', 'Simple', q4, 2)
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

    C.set_ang_vel(N, u0*C.x + u1*C.y + u2*C.z)
    assert C.ang_vel_in(N) == (u0)*C.x + (u1)*C.y + (u2)*C.z
    assert N.ang_vel_in(C) == (-u0)*C.x + (-u1)*C.y + (-u2)*C.z
    assert C.ang_vel_in(D) == (u0)*C.x + (u1)*C.y + (u2)*C.z + (-q4d)*D.y
    assert D.ang_vel_in(C) == (-u0)*C.x + (-u1)*C.y + (-u2)*C.z + (q4d)*D.y

    q0, q0d = dynamicsymbols('q0 q0d')
    E = N.orientnew('E', 'Quaternion', (q0, q1, q2, q3))
    assert E.ang_vel_in(N) == (
        2 * (q1d * q0 + q2d * q3 - q3d * q2 - q0d * q1) * E.x +
        2 * (q2d * q0 + q3d * q1 - q1d * q3 - q0d * q2) * E.y +
        2 * (q3d * q0 + q1d * q2 - q2d * q1 - q0d * q3) * E.z)

    F = N.orientnew('F', 'Body', (q1, q2, q3), '123')
    assert F.ang_vel_in(N) == (
        (sin(q3) * q2d + cos(q2) * cos(q3) * q1d) * F.x + (cos(q3) * q2d -
        sin(q3) * cos(q2) * q1d) * F.y + (q3d + sin(q2) * q1d) * F.z)
    G = N.orientnew('G', 'Axis', (q1, N.x + N.y))
    assert G.ang_vel_in(N) == q1d * (N.x + N.y).unit
    assert N.ang_vel_in(G) == -q1d * (N.x + N.y).unit

def test_dcm():
    q1, q2, q3, q4 = dynamicsymbols('q1 q2 q3 q4')
    N = ReferenceFrame('N')
    A = N.orientnew('A', 'Simple', q1, 3)
    B = A.orientnew('B', 'Simple', q2, 1)
    C = B.orientnew('C', 'Simple', q3, 2)
    D = N.orientnew('D', 'Simple', q4, 2)
    assert N.dcm(C) == Matrix([
        [- sin(q1) * sin(q2) * sin(q3) + cos(q1) * cos(q3), - sin(q1) *
        cos(q2), sin(q1) * sin(q2) * cos(q3) + sin(q3) * cos(q1)], [sin(q1) *
        cos(q3) + sin(q2) * sin(q3) * cos(q1), cos(q1) * cos(q2), sin(q1) *
        sin(q3) - sin(q2) * cos(q1) * cos(q3)], [- sin(q3) * cos(q2), sin(q2),
        cos(q2) * cos(q3)]])
    # This is a little touchy.  Is it ok to use simplify in assert?
    assert D.dcm(C) == Matrix(
        [[cos(q1) * cos(q3) * cos(q4) - sin(q3) * (- sin(q4) * cos(q2) +
        sin(q1) * sin(q2) * cos(q4)), - sin(q2) * sin(q4) - sin(q1) *
        cos(q2) * cos(q4), sin(q3) * cos(q1) * cos(q4) + cos(q3) * (- sin(q4) *
        cos(q2) + sin(q1) * sin(q2) * cos(q4))], [sin(q1) * cos(q3) +
        sin(q2) * sin(q3) * cos(q1), cos(q1) * cos(q2), sin(q1) * sin(q3) -
        sin(q2) * cos(q1) * cos(q3)], [sin(q4) * cos(q1) * cos(q3) -
        sin(q3) * (cos(q2) * cos(q4) + sin(q1) * sin(q2) * sin(q4)), sin(q2) *
        cos(q4) - sin(q1) * sin(q4) * cos(q2), sin(q3) * sin(q4) * cos(q1) +
        cos(q3) * (cos(q2) * cos(q4) + sin(q1) * sin(q2) * sin(q4))]])

def test_Vector():
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
    assert dot(v3, A.x) == x + x**2
    assert dot(v3, A.y) == y + y**2
    assert dot(v3, A.z) == z + z**2

    assert isinstance(v4, Vector)
    # We probably shouldn't be using simplify in dot...
    assert dot(v4, A.x) == x - x**2
    assert dot(v4, A.y) == y - y**2
    assert dot(v4, A.z) == z - z**2

def test_Vector_diffs():
    q1, q2, q3, q4, q1d, q2d, q3d, q4d, q1dd, q2dd, q3dd, q4dd = (
        dynamicsymbols('q1 q2 q3 q4 q1d q2d q3d q4d q1dd q2dd q3dd q4dd'))
    N = ReferenceFrame('N')
    A = N.orientnew('A', 'Simple', q3, 3)
    B = A.orientnew('B', 'Simple', q2, 1)
    v1 = q2 * A.x + q3 * N.y
    v2 = q3 * B.x + v1
    v3 = v1.dt(B)
    v4 = v2.dt(B)

    assert v1.dt(N) == q2d * A.x + q2 * q3d * A.y + q3d * N.y
    assert v1.dt(A) == q2d * A.x + q3 * q3d * N.x + q3d * N.y
    assert v1.dt(B) == (q2d * A.x + q3 * q3d * N.x + q3d * N.y - q3 * cos(q3) *
                        q2d * N.z)
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




    #TODO: Put some tests here, mew.
