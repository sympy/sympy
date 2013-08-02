from sympy import symbols, sqrt
from sympy.physics.mechanics import MovingRefFrame, dynamicsymbols


O = MovingRefFrame('O')
a, b, c = dynamicsymbols('a b c')
ad, bd, cd = dynamicsymbols('a b c', 1)
a2d, b2d, c2d = dynamicsymbols('a b c', 2)
v = a*O.x + b*O.y + c*O.z
vd = ad*O.x + bd*O.y + cd*O.z
v2d = a2d*O.x + b2d*O.y + c2d*O.z
v1 = 3*O.x + 4*O.y + 5*O.z
C, t1, t2 = symbols('C t1 t2')
#TODO- Add tests for tuple args

def test_pos_vector():
    """ Tests for relative position vectors """
    
    P = MovingRefFrame('P', pos_vector = v1, parentframe = O)
    assert P.pos_vector_in(O) == v1
    assert O.pos_vector_in(P) == -v1
    v2 = 6*P.x + 7*P.y + 8*P.z
    Q = MovingRefFrame('Q', pos_vector = v2, parentframe = P)
    assert Q.pos_vector_in(O) == v1 + v2
    assert O.pos_vector_in(Q) == - v1 - v2
    R = MovingRefFrame('R', pos_vector = v1 + v2, parentframe = O)
    assert Q.pos_vector_in(R) == 0
    S = MovingRefFrame('S', trans_vel = v, pos_vector_b = v1, t = t1,
                       parentframe = O)
    assert S.pos_vector_in(O) == (Integral(a, t) - Integral(a(t1), t) + 3)*O.x + \
           (Integral(b, t) - Integral(b(t1), t) + 4)*O.y + \
           (Integral(c, t) - Integral(c(t1), t) + 5)*O.z
    T = MovingRefFrame('T', trans_acc = v, trans_vel_b = v1, pos_vector_b = v2,
                       t1 = t1, t2 = t2, parentframe = O)
    #TODO- assert T.pos_vector_in(O) == 


def test_trans_vel():
    """ Tests for relative translational velocity """
    
    P = MovingRefFrame('P', trans_vel = v1, parentframe = O)
    assert P.trans_vel_in(O) == v1
    assert O.trans_vel_in(P) == -v1
    Q = MovingRefFrame('Q', ang_vel = C * O.z, parentframe = O)
    assert Q.trans_vel_in(O) == 0
    R = MovingRefFrame('R', pos_vector = Q.x, parentframe = Q)
    assert R.trans_vel_in(O) == C*Q.y
    assert R.trans_vel_in(P) == C*Q.y - v1
    S = MovingRefFrame('S', pos_vector = v, parentframe = O)
    assert S.trans_vel_in(O) == vd
    assert S.trans_vel_in(R) == C*b*O.x - C*a*O.y + vd
    T = MovingRefFrame('T', trans_acc = v, trans_vel_b = v1, t = t1, parentframe = O)
    assert T.trans_vel_in(O) == (Integral(a, t) - Integral(a(t1), t) + 3)*O.x + \
           (Integral(b, t) - Integral(b(t1), t) + 4)*O.y + \
           (Integral(c, t) - Integral(c(t1), t) + 5)*O.z


def test_trans_acc():
    """ Tests for relative translational acceleration """
    
    P = MovingRefFrame('P', trans_acc = v1, parentframe = O)
    assert P.trans_acc_in(O) == v1
    assert O.trans_acc_in(P) == -v1
    Q = MovingRefFrame('Q', ang_vel = C * O.z, parentframe = O)
    assert Q.trans_acc_in(O) == 0
    R = MovingRefFrame('R', pos_vector = Q.x, parentframe = Q)
    assert R.trans_acc_in(O) == - C**2*Q.x
    assert R.trans_acc_in(P) == - C**2*Q.x - v1
    S = MovingRefFrame('S', pos_vector = v, parentframe = O)
    assert S.trans_acc_in(O) == v2d
    assert S.trans_acc_in(R) == (C*(-C*a + bd) + C*bd)*O.x + \
           (C*(-C*b + ad) - C*ad)*O.y + v2d
    T = MovingRefFrame('T', trans_vel = v, parentframe = O)
    assert T.trans_acc_in(O) == vd


def test_ang_vel():
    """ Tests for relative angular velocity """
    
    P = MovingRefFrame('P', ang_vel = v1, parentframe = O)
    assert P.ang_vel_in(O) == v1
    assert O.ang_vel_in(P) == - v1
    Q = MovingRefFrame('Q', orient_type = 'Axis', orient_amount = [a, O.z],
                       parentframe = O)
    assert Q.ang_vel_in(O) == ad * Q.z
    assert Q.ang_vel_in(P) == ad*Q.z - v1
    v2 = 6*P.x + 7*P.y + 8*P.z
    R = MovingRefFrame('R', ang_acc = v1, ang_vel_b = v2, rt1 = t1, parentframe = Q)
    #TODO - assert R.ang_vel_in(O) == ...
    #TODO - test dcm with ang vel + rotation boundary condition


def test_ang_acc():
    """ Tests for relative angular acceleration """
    
    P = MovingRefFrame('P', ang_acc = v1, parentframe = O)
    assert P.ang_acc_in(O) == v1
    assert O.ang_acc_in(P) == -v1
    Q = MovingRefFrame('Q', orient_type = 'Axis', orient_amount = [a, O.z],
                       parentframe = O)
    assert Q.ang_acc_in(O) == a2d * O.z
    assert Q.ang_acc_in(P) == (-20*t**2 - 4*t*(-5*t + ad) - 3)*O.x + \
           (15*t**2 + 3*t*(-5*t + ad) - 4)*O.y + (a2d - 5)*O.z
    assert P.ang_acc_in(Q) == (4*t*ad + 3)*O.x + (-3*t*ad + 4)*O.y + (-a2d + 5)*O.z
    R = MovingRefFrame('R', ang_vel = v, parentframe = O)
    assert R.ang_acc_in(O) == vd
    assert P.ang_acc_in(R) == ((4*t - b)*c - (5*t - c)*b - ad + 3)*O.x + \
           (-(3*t - a)*c + (5*t - c)*a - bd + 4)*O.y + \
           ((3*t - a)*b - (4*t - b)*a - cd + 5)*O.z


def test_dt():
    """ Test time-differentiation """
    
    assert O.dt(v) == vd
    O1 = MovingRefFrame('O1', trans_acc = v1, parentframe = O)
    assert O1.dt(v) == vd
    P = MovingRefFrame('P', orient_type = 'Axis', orient_amount = [a, O.x + O.y + O.z],
                       parentframe = O)
    assert P.dt(v1) == - sqrt(3)*ad/3*O.x + 2*sqrt(3)*ad/3*O.y - sqrt(3)*ad/3*O.z
    Q = MovingRefFrame('Q', orient_type = 'Axis', orient_amount = [C, O.z],
                       parentframe = O)
    assert Q.dt(v) == vd
    R = MovingRefFrame('R', ang_vel = v1, parentframe = O)
    assert R.dt(v) == (5*b - 4*c + ad)*O.x + (-5*a + 3*c + bd)*O.y + \
           (4*a - 3*b + cd)*O.z
    assert R.dt(Q.x) == - 5*Q.y + (-3*sin(C) + 4*cos(C))*Q.z
    assert R.dt(Q.y) == 5*Q.x + (-4*sin(C) - 3*cos(C))*Q.z
    assert R.dt(Q.z) == (3*sin(C) - 4*cos(C))*Q.x + (4*sin(C) + 3*cos(C))*Q.y
    #TODO- Insert tests with coordinate variables
