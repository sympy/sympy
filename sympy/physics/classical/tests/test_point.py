from sympy import symbols, Matrix, sin, cos
from sympy.physics.classical import (Vector, ReferenceFrame, dot, cross,
                                     dynamicsymbols, Point)

def test_point_v1pts():
    q, q2, qd, q2d, qdd, q2dd = dynamicsymbols('q q2 qd q2d qdd q2dd')
    N = ReferenceFrame('N')
    B = ReferenceFrame('B')
    B.set_ang_vel(N, qd * B.z)
    O = Point('O')
    P = O.newpoint('P', B.x)
    P.set_vel(0, B)
    O.set_vel(0, N)
    assert P.v1pt(O, N, B) == qd * B.y
    O.set_vel(N.x, N)
    assert P.v1pt(O, N, B) == N.x + qd * B.y
    P.set_vel(B.z, B)
    assert P.v1pt(O, N, B) == B.z + N.x + qd * B.y

def test_point_a1pts():
    q, q2, qd, q2d, qdd, q2dd = dynamicsymbols('q q2 qd q2d qdd q2dd')
    N = ReferenceFrame('N')
    B = ReferenceFrame('B')
    B.set_ang_vel(N, qd * B.z)
    O = Point('O')
    P = O.newpoint('P', B.x)
    P.set_vel(0, B)
    O.set_vel(0, N)
    assert P.a1pt(O, N, B) ==  -(qd**2) * B.x + qdd * B.y
    P.set_vel(q2d * B.z, B)
    assert P.a1pt(O, N, B) == -(qd**2) * B.x + qdd * B.y + q2dd * B.z
    O.set_vel(q2d * B.x, N)
    assert P.a1pt(O, N, B) == ((q2dd - qd**2) * B.x + (q2d * qd + qdd) * B.y +
                               q2dd * B.z)

def test_point_v2pts():
    pass

def test_point_a2pts():
    pass

def test_point_funcs():
    q, q2, qd, q2d, qdd, q2dd = dynamicsymbols('q q2 qd q2d qdd q2dd')
    N = ReferenceFrame('N')
    B = ReferenceFrame('B')
    B.set_ang_vel(N, 5 * B.y)
    O = Point('O')
    P = O.newpoint('P', q * B.x)
    assert P.pos_from(O) == q * B.x
    P.set_vel(qd * B.x + q2d * B.y, B)
    assert P.vel(B) == qd * B.x + q2d * B.y
    O.set_vel(0, N)
    assert O.vel(N) == 0
    assert P.a1pt(O, N, B) == ((-25 * q + qdd) * B.x + (q2dd) * B.y +
                               (-10 * qd) * B.z)

    B = N.orientnew('B', 'Simple', q, 3)
    O = Point('O')
    P = O.newpoint('P', 10 * B.x)
    O.set_vel(5 * N.x, N)
    assert O.vel(N) == 5 * N.x
    assert P.a2pt(O, N, B) == (-10 * qd**2) * B.x + (10 * qdd) * B.y

    B.set_ang_vel(N, 5 * B.y)
    O = Point('O')
    P = O.newpoint('P', q * B.x)
    P.set_vel(qd * B.x + q2d * B.y, B)
    O.set_vel(0, N)
    assert P.v1pt(O, N, B) == qd * B.x + q2d * B.y - 5 * q * B.z

def test_point_pos():
    q = dynamicsymbols('q')
    N = ReferenceFrame('N')
    B = N.orientnew('B', 'Simple', q, 3)
    O = Point('O')
    P = O.newpoint('P', 10 * N.x + 5 * B.x)
    assert P.pos_from(O) == 10 * N.x + 5 * B.x
    Q = P.newpoint('Q', 10 * N.y + 5 * B.y)
    assert Q.pos_from(P) == 10 * N.y + 5 * B.y
    assert Q.pos_from(O) == 10 * N.x + 10 * N.y + 5 * B.x + 5 * B.y
    assert O.pos_from(Q) == -10 * N.x - 10 * N.y - 5 * B.x - 5 * B.y


