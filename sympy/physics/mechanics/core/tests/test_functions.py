from sympy import symbols
from sympy.utilities.pytest import raises
from sympy.physics.mechanics import get_motion_pos, get_motion_vel, get_motion_acc,\
     dynamicsymbols

F = MovingRefFrame('F')
t = dynamicsymbols._t
s1, s2, s3 = symbols('s1 s2 s3')
S1, S2, S3 = symbols('S1 S2 S3')
S4, S5, S6 = symbols('S4 S5 S6')
t1, t2 = symbols('t1 t2')
a, b, c = dynamicsymbols('a b c')
ad, bd, cd = dynamicsymbols('a b c', 1)
a2d, b2d, c2d = dynamicsymbols('a b c', 2)
v0 = S1*F.x + S2*F.y + S3*F.z
v01 = S4*F.x + S5*F.y + S6*F.z
v1 = s1*F.x + s2*F.y + s3*F.z
v2 = a*F.x + b*F.y + c*F.z
v2d = ad*F.x + bd*F.y + cd*F.z
v2dd = a2d*F.x + b2d*F.y + c2d*F.z


def test_get_motion_pos():
    assert get_motion_pos(0, F) == [0, 0, 0]
    assert get_motion_pos(v1, F) == [0, 0, v1]
    assert get_motion_pos(v2, F) == [v2dd, v2d, v2]
    raises (TypeError, get_motion_pos())
    raises (TypeError, get_motion_pos(s1, F))


def test_get_motion_vel():
    assert get_motion_vel(frame = F) == [0, 0, 0]
    assert get_motion_vel(v1, frame = F) == [0, v1, v1 * t]
    assert get_motion_vel(v1, v0, t1, F) == [0, v1, v0 + v1*(t - t1)]
    assert get_motion_vel(v1, v2, t1, F) == [0, v1, v1*t - v1*t1 + v2.subs(t, t1)]
    assert get_motion_vel(v2, v0, t1, F) == \
           [v2d, v2,
            (S1 + Integral(a, t) - Integral(a(t1), t))*F.x + \
            (S2 + Integral(b, t) - Integral(b(t1), t))*F.y + \
            (S3 + Integral(c, t) - Integral(c(t1), t))*F.z]
    raises (TypeError, get_motion_vel())
    raises (TypeError, get_motion_vel(s1, v0, t1, F))
    raises (TypeError, get_motion_vel(v1, s1, t1, F))
    raises (ValueError, get_motion_vel(v1, v0, a, F))


def test_get_motion_acc():
    assert get_motion_acc(frame = F) == [0, 0, 0]
    assert get_motion_acc(v1, frame = F) == [v1, v1 * t, v1 * t**2/2]
    assert get_motion_acc(v1, v0, v2, t1, t2, F) == [v1, (v0 + v1*t - v1*t1),
                                                     -v0*t2 + v1*t**2/2 + v1*t1*t2 - \
                                                     v1*t2**2/2 + t*(v0 - v1*t1) + \
                                                     v2.subs(t, t2)]
    assert get_motion_acc(v1, v0, v01, t1, t2, F) == [v1, v0 + v1*t - v1*t1,
                                                      -v0*t2 + v01 + v1*t**2/2 + \
                                                      v1*t1*t2 - v1*t2**2/2 + \
                                                      t*(v0 - v1*t1)]
    assert get_motion_acc(a*F.x, S1*F.x,
                          S2*F.x, t1, t2, F) == [a*F.x, (S1 + Integral(a, t) - \
                                                         Integral(a(t1), t))*F.x,
                                                 (S2 + Integral(S1 + Integral(a, t) -\
                                                                Integral(a(t1), t),
                                                                t) - \
                                                  Integral(S1 - Integral(a(t1), t) + \
                                                           Integral(a(t2), t), t))*F.x]
     raises (TypeError, get_motion_acc())
     raises (TypeError, get_motion_acc(S1, v0, v01, t1, t2, F))
     raises (TypeError, get_motion_acc(v1, S1, v01, t1, t2, F))
     raises (TypeError, get_motion_acc(v1, v0, S1, t1, t2, F))
     raises (ValueError, get_motion_acc(v1, v0, v01, a, b, F))
