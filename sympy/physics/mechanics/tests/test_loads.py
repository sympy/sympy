from sympy.physics.mechanics.loads import Force, Torque, gravity
from sympy.physics.mechanics import (RigidBody, Particle, ReferenceFrame, Point,
                                     outer, dynamicsymbols)
from sympy.core.backend import symbols


def test_force_default():
    N = ReferenceFrame('N')
    Po = Point('Po')
    f1 = Force(Po, N.x)
    assert f1.point == Po
    assert f1.force == N.x
    assert f1.__repr__() == 'Force(point=Po, force=N.x)'
    # Test tuple behaviour
    assert isinstance(f1, tuple)
    assert f1[0] == Po
    assert f1[1] == N.x
    assert f1 == (Po, N.x)
    assert f1 != (N.x, Po)
    assert f1 != (Po, N.x + N.y)
    assert f1 != (Point('Co'), N.x)
    # Test body as input
    P = Particle('P', Po)
    f2 = Force(P, N.x)
    assert f1 == f2


def test_torque_default():
    N = ReferenceFrame('N')
    f1 = Torque(N, N.x)
    assert f1.frame == N
    assert f1.torque == N.x
    assert f1.__repr__() == 'Torque(frame=N, torque=N.x)'
    # Test tuple behaviour
    assert isinstance(f1, tuple)
    assert f1[0] == N
    assert f1[1] == N.x
    assert f1 == (N, N.x)
    assert f1 != (N.x, N)
    assert f1 != (N, N.x + N.y)
    assert f1 != (ReferenceFrame('A'), N.x)
    # Test body as input
    rb = RigidBody('P', frame=N)
    f2 = Torque(rb, N.x)
    assert f1 == f2


def test_gravity():
    N = ReferenceFrame('N')
    m, M, g = symbols('m M g')
    F1, F2 = dynamicsymbols('F1 F2')
    po = Point('po')
    pa = Particle('pa', po, m)
    A = ReferenceFrame('A')
    P = Point('P')
    I = outer(A.x, A.x)
    B = RigidBody('B', P, A, M, (I, P))
    forceList = [(po, F1), (P, F2)]
    forceList.extend(gravity(g * N.y, pa, B))
    l = [(po, F1), (P, F2), (po, g * m * N.y), (P, g * M * N.y)]

    for i in range(len(l)):
        for j in range(len(l[i])):
            assert forceList[i][j] == l[i][j]
