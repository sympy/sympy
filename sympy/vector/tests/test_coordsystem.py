from sympy.vector.coordsys import CoordinateSystem, CartesianCoordinateSystem, SphericalCoordinateSystem
from sympy.vector.scalar import BaseScalar
from sympy import sin, cos, pi, ImmutableMatrix as Matrix, oo

def test_coord_system_superclass():
    assert issubclass(CartesianCoordinateSystem, CoordinateSystem)
    assert issubclass(SphericalCoordinateSystem, CoordinateSystem)


def test_spherical_coord_system():
    S = SphericalCoordinateSystem('S')
    assert S.base_scalars() == (S.r, S.theta, S.phi)
    assert S.base_vectors() == (S.e_r, S.e_theta, S.e_phi)
    assert S.coordinate_relations() == (S.r*sin(S.theta)*cos(S.phi),
                                   S.r*sin(S.theta)*sin(S.phi),
                                   S.r*cos(S.theta))

    assert S.coordinate_metric() == Matrix([[1, 0, 0],
                                 [0, S.r**2, 0],
                                 [0, 0, S.r**2 * sin(S.theta)**2]])

    assert BaseScalar('S.r', 0, S, 'S_r', r'\mathbf{{r}_{S}}') == S.r
    assert BaseScalar('S.theta', 1, S, 'S_theta', r'\mathbf{{theta}_{S}}') == S.theta
    assert BaseScalar('S.phi', 2, S, 'S_phi', r'\mathbf{{phi}_{S}}') == S.phi

    assert isinstance(S.r, BaseScalar) and \
        isinstance(S.theta, BaseScalar) and \
        isinstance(S.phi, BaseScalar)
