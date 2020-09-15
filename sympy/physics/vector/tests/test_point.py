from sympy.physics.vector import dynamicsymbols, Point, ReferenceFrame
from sympy.testing.pytest import raises
from sympy import Derivative


def test_point_v1pt_theorys():
    q, q2 = dynamicsymbols('q q2')
    qd, q2d = dynamicsymbols('q q2', 1)
    qdd, q2dd = dynamicsymbols('q q2', 2)
    N = ReferenceFrame('N')
    B = ReferenceFrame('B')
    B.set_ang_vel(N, qd * B.z)
    O = Point('O')
    P = O.locatenew('P', B.x)
    P.set_vel(B, 0)
    O.set_vel(N, 0)
    assert P.v1pt_theory(O, N, B) == qd * B.y
    O.set_vel(N, N.x)
    assert P.v1pt_theory(O, N, B) == N.x + qd * B.y
    P.set_vel(B, B.z)
    assert P.v1pt_theory(O, N, B) == B.z + N.x + qd * B.y


def test_point_a1pt_theorys():
    q, q2 = dynamicsymbols('q q2')
    qd, q2d = dynamicsymbols('q q2', 1)
    qdd, q2dd = dynamicsymbols('q q2', 2)
    N = ReferenceFrame('N')
    B = ReferenceFrame('B')
    B.set_ang_vel(N, qd * B.z)
    O = Point('O')
    P = O.locatenew('P', B.x)
    P.set_vel(B, 0)
    O.set_vel(N, 0)
    assert P.a1pt_theory(O, N, B) == -(qd**2) * B.x + qdd * B.y
    P.set_vel(B, q2d * B.z)
    assert P.a1pt_theory(O, N, B) == -(qd**2) * B.x + qdd * B.y + q2dd * B.z
    O.set_vel(N, q2d * B.x)
    assert P.a1pt_theory(O, N, B) == ((q2dd - qd**2) * B.x + (q2d * qd + qdd) * B.y +
                               q2dd * B.z)


def test_point_v2pt_theorys():
    q = dynamicsymbols('q')
    qd = dynamicsymbols('q', 1)
    N = ReferenceFrame('N')
    B = N.orientnew('B', 'Axis', [q, N.z])
    O = Point('O')
    P = O.locatenew('P', 0)
    O.set_vel(N, 0)
    assert P.v2pt_theory(O, N, B) == 0
    P = O.locatenew('P', B.x)
    assert P.v2pt_theory(O, N, B) == (qd * B.z ^ B.x)
    O.set_vel(N, N.x)
    assert P.v2pt_theory(O, N, B) == N.x + qd * B.y


def test_point_a2pt_theorys():
    q = dynamicsymbols('q')
    qd = dynamicsymbols('q', 1)
    qdd = dynamicsymbols('q', 2)
    N = ReferenceFrame('N')
    B = N.orientnew('B', 'Axis', [q, N.z])
    O = Point('O')
    P = O.locatenew('P', 0)
    O.set_vel(N, 0)
    assert P.a2pt_theory(O, N, B) == 0
    P.set_pos(O, B.x)
    assert P.a2pt_theory(O, N, B) == (-qd**2) * B.x + (qdd) * B.y


def test_point_funcs():
    q, q2 = dynamicsymbols('q q2')
    qd, q2d = dynamicsymbols('q q2', 1)
    qdd, q2dd = dynamicsymbols('q q2', 2)
    N = ReferenceFrame('N')
    B = ReferenceFrame('B')
    B.set_ang_vel(N, 5 * B.y)
    O = Point('O')
    P = O.locatenew('P', q * B.x)
    assert P.pos_from(O) == q * B.x
    P.set_vel(B, qd * B.x + q2d * B.y)
    assert P.vel(B) == qd * B.x + q2d * B.y
    O.set_vel(N, 0)
    assert O.vel(N) == 0

    assert P.a1pt_theory(O, N, B) == ((-25 * q + qdd) * B.x + (q2dd) * B.y +
                               (-10 * qd) * B.z)

    B = N.orientnew('B', 'Axis', [q, N.z])
    O = Point('O')
    P = O.locatenew('P', 10 * B.x)
    O.set_vel(N, 5 * N.x)
    assert O.vel(N) == 5 * N.x
    assert P.a2pt_theory(O, N, B) == (-10 * qd**2) * B.x + (10 * qdd) * B.y

    B.set_ang_vel(N, 5 * B.y)
    O = Point('O')
    P = O.locatenew('P', q * B.x)
    P.set_vel(B, qd * B.x + q2d * B.y)
    O.set_vel(N, 0)
    assert P.v1pt_theory(O, N, B) == qd * B.x + q2d * B.y - 5 * q * B.z


def test_point_pos():
    q = dynamicsymbols('q')
    N = ReferenceFrame('N')
    B = N.orientnew('B', 'Axis', [q, N.z])
    O = Point('O')
    P = O.locatenew('P', 10 * N.x + 5 * B.x)
    assert P.pos_from(O) == 10 * N.x + 5 * B.x
    Q = P.locatenew('Q', 10 * N.y + 5 * B.y)
    assert Q.pos_from(P) == 10 * N.y + 5 * B.y
    assert Q.pos_from(O) == 10 * N.x + 10 * N.y + 5 * B.x + 5 * B.y
    assert O.pos_from(Q) == -10 * N.x - 10 * N.y - 5 * B.x - 5 * B.y

def test_point_partial_velocity():

    N = ReferenceFrame('N')
    A = ReferenceFrame('A')

    p = Point('p')

    u1, u2 = dynamicsymbols('u1, u2')

    p.set_vel(N, u1 * A.x + u2 * N.y)

    assert p.partial_velocity(N, u1) == A.x
    assert p.partial_velocity(N, u1, u2) == (A.x, N.y)
    raises(ValueError, lambda: p.partial_velocity(A, u1))

def test_point_vel():
    q1, q2 = dynamicsymbols('q1 q2')
    N = ReferenceFrame('N')
    B = ReferenceFrame('B')
    O = Point('O')
    Q = Point('Q')
    Q.set_pos(O, q1 * N.x)
    raises(ValueError , lambda : Q.vel(N))
    O.set_vel(N, q2 * N.y)
    assert Q.vel(N) == Derivative(q1) * N.x + q2 * N.y
    P1 = Point('P1')
    P1.set_pos(O, q1 * B.x)
    raises(ValueError, lambda : P1.vel(B))
    P2 = Point('P2')
    P2.set_pos(P1, q2 * B.z)
    raises(ValueError, lambda : P2.vel(B))
    P3 = Point('P3')
    P3.set_pos(P2, q1 * B.x)
    raises(ValueError, lambda : P3.vel(B))
    P1.set_vel(B, 10 * q1 * B.x)
    assert P3.vel(B) == (10 * q1 + Derivative(q1)) * B.x + Derivative(q2) * B.z
    Q1 = Point('Q1')
    Q1.set_vel(N, q1 * N.x)
    Q2 = Point('Q2')
    Q3 = Point('Q3')
    Q4 = Point('Q4')
    Q2.set_pos(Q1, q2 * N.x)
    Q3.set_pos(Q2, q1 * N.x)
    Q4.set_pos(Q3, q2 * N.x)
    assert Q4.vel(N) == (q1 + Derivative(q1) + 2 * Derivative(q2)) * N.x
    Q3.set_vel(N, q1 * N.x)
    assert Q3.vel(N) == q1 * N.x
    assert Q4.vel(N) == (q1 + Derivative(q2)) * N.x
