from sympy.vector.coordsysrect import CoordSysCartesian
from sympy.vector.coordsystem import (CurvilinearCoordSys, CoordSysSpherical,
                                      CoordSysCylindrical)
from sympy.vector.scalar import BaseScalar
from sympy import sin, cos, pi, ImmutableMatrix as Matrix, oo


def test_coord_system_superclass():
    assert issubclass(CoordSysCartesian, CurvilinearCoordSys)
    assert issubclass(CoordSysSpherical, CurvilinearCoordSys)
    assert issubclass(CoordSysCylindrical, CurvilinearCoordSys)

    assert issubclass(CurvilinearCoordSys, Basic)


def test_spherical_coord_system():
    S = CoordSysSpherical('S')
    assert S.base_scalars() == (S.r, S.theta, S.phi)
    assert S.base_vectors() == (S.rhat, S.thetahat, S.phihat)
    assert S.limits() == ((0, oo), (0, pi), (0, 2*pi))
    
    assert S.coord_relations() == (S.r*sin(S.theta)*cos(S.phi),
                                   S.r*sin(S.theta)*sin(S.phi),
                                   S.r*cos(S.theta))
    assert S.metric() == Matrix([[1, 0, 0],
                                 [0, S.r**2, 0],
                                 [0, 0, S.r**2 * sin(S.theta)**2]])

    assert BaseScalar('r', 0, S, ' ', ' ') == S.r
    assert BaseScalar('theta', 1, S, ' ', ' ') == S.theta
    assert BaseScalar('phi', 2, S, ' ', ' ') == S.phi
    assert isinstance(S.r, BaseScalar) and \
        isinstance(S.theta, BaseScalar) and \
        isinstance(S.phi, BaseScalar)


def test_cylindrical_coord_system():
    C = CoordSysCylindrical('C')
    assert C.base_scalars() == (S.rho, S.theta, S.z)
    assert C.base_vectors() == (S.rhohat, S.phihat, S.zhat)
    assert C.limits() == ((0, oo), (0, 2*pi), (-oo, oo))

    assert C.coord_relations() == (C.rho*cos(C.theta),
                                   C.rho*sin(C.theta),
                                   C.z)
    assert C.metric() == Matrix([[1, 0, 0],
                                 [0, C.rho**2, 0],
                                 [0, 0, 1]])

    assert BaseScalar('rho', 0, C, ' ', ' ') == C.rho
    assert BaseScalar('theta', 1, C, ' ', ' ') == C.theta
    assert BaseScalar('z', 2, C, ' ', ' ') == C.z
    assert isinstance(C.rho, BaseScalar) and \
        isinstance(C.theta, BaseScalar) and \
        isinstance(C.z, BaseScalar)


def test_general_coordinate_setup():
    S = CoordSysSpherical('S')
    gen = CurvilinearCoordSys('gen', coord_type='spherical')
    assert gen.base_scalars() == (gen.q1, gen.q2, gen.q3)
    assert gen.base_vectors() == (gen.e1, gen.e2, gen.e3)
    assert gen.limits() == ((-oo, oo), (-oo, oo), (-oo, oo))

    repl = [(S.r, gen.q1), (S.theta, gen.q2), (S.phi, gen.q3)]
    assert Matrix(gen.coord_relations()) == \
           Matrix(S.coord_relations()).subs(repl)
    assert gen.metric() == S.metric().subs(repl)

    assert not hasattr(gen, 'rotation_matrix')
    assert not hasattr(gen, 'locate_new')
    assert not hasattr(gen, 'orient_new')
