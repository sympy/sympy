from sympy.physics.vector import dynamicsymbols, Point, ReferenceFrame
from sympy.testing.pytest import raises


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

def test_point_vel(): #Basic functionality
    q1, q2 = dynamicsymbols('q1 q2')
    N = ReferenceFrame('N')
    B = ReferenceFrame('B')
    Q = Point('Q')
    O = Point('O')
    Q.set_pos(O, q1 * N.x)
    raises(ValueError , lambda: Q.vel(N)) #Velocity of Q is not defined
    O.set_vel(N, q2 * N.y)
    assert O.vel(N) == q2 * N.y
    raises(ValueError , lambda : O.vel(B)) #Velocity of O is not defined in B

def test_auto_point_vel():
    t = dynamicsymbols._t
    q1, q2 = dynamicsymbols('q1 q2')
    N = ReferenceFrame('N')
    B = ReferenceFrame('B')
    O = Point('O')
    Q = Point('Q')
    Q.set_pos(O, q1 * N.x)
    O.set_vel(N, q2 * N.y)
    assert Q.vel(N) == q1.diff(t) * N.x + q2 * N.y  # Velocity of Q using O
    P1 = Point('P1')
    P1.set_pos(O, q1 * B.x)
    P2 = Point('P2')
    P2.set_pos(P1, q2 * B.z)
    raises(ValueError, lambda : P2.vel(B)) # O's velocity is defined in different frame, and no
    #point in between has its velocity defined
    raises(ValueError, lambda: P2.vel(N)) # Pos vector defined in different frame

def test_auto_point_vel_shortest_path():
    t = dynamicsymbols._t
    q1, q2 = dynamicsymbols('q1 q2')
    B = ReferenceFrame('B')
    P = Point('P')
    P.set_vel(B, q1 * B.x)
    P1 = Point('P1')
    P1.set_pos(P, q2 * B.y)
    P1.set_vel(B, q1 * B.z)
    P2 = Point('P2')
    P2.set_pos(P1, q1 * B.z)
    P3 = Point('P3')
    P3.set_pos(P2, 10 * q1 * B.y)
    assert P3.vel(B) == 10 * q1.diff(t) * B.y + (q1 + q1.diff(t)) * B.z

def test_auto_vel_dont_overwrite():
    t = dynamicsymbols._t
    q1, q2 = dynamicsymbols('q1 q2')
    N = ReferenceFrame('N')
    P = Point('P1')
    P.set_vel(N, q1 * N.x)
    P1 = Point('P1')
    P1.set_pos(P, q2 * N.y)
    assert P1.vel(N) == q2.diff(t) * N.y + q1 * N.x
    assert P.vel(N) == q1 * N.x
    P1.set_vel(N, q1 * N.z)
    assert P1.vel(N) == q1 * N.z

def test_auto_point_vel_multiple_frames():
    t = dynamicsymbols._t
    q1, q2 = dynamicsymbols('q1 q2')
    B = ReferenceFrame('B')
    S = ReferenceFrame('S')
    N = ReferenceFrame('N')
    T1 = Point('T1')
    T1.set_vel(S, q1 * S.x)
    T2 = Point('T2')
    T2.set_pos(T1, q2 * S.x)
    T2.set_vel(N, q1 * N.y)
    T3 = Point('T3')
    T3.set_vel(B, q1 * B.z)
    T3.set_pos(T2, q1 * B.x + q2 * B.y)
    T3.set_pos(T1, q2 * S.x)
    T4 = Point('T4')
    T4.set_pos(T3, q1 * S.z)
    T5 = Point('T5')
    T5.set_pos(T4, q2 * S.y)
    T6 = Point('T6')
    T6.set_pos(T2, q1 * S.z)
    assert T6.vel(S) == (q1 + q2.diff(t)) * S.x + q1.diff(t) * S.z
    T6.set_pos(T3, q1 * B.z)
    assert T6.vel(B) == (q1 + q1.diff(t)) * B.z

def test_auto_point_vel_if_tree_has_vel_but_inappropriate_pos_vector():
    q1, q2 = dynamicsymbols('q1 q2')
    B = ReferenceFrame('B')
    S = ReferenceFrame('S')
    P = Point('P')
    P.set_vel(B, q1 * B.x)
    P1 = Point('P1')
    P1.set_pos(P, S.y)
    raises(ValueError, lambda : P1.vel(B)) # pos vector not according to frame
    raises(ValueError, lambda : P1.vel(S)) #velocity not defined wrt frame in tree
