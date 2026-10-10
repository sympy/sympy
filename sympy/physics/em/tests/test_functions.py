from sympy.physics.mechanics import ReferenceFrame
from sympy.physics.em import gradient, curl, is_conservative
from sympy.utilities.pytest import raises

R = ReferenceFrame('R')
scalar_field = 2*R[0]**2*R[1]*R[2]
grad_field = gradient(scalar_field, R)
vector_field = R[1]**2*R.x + 3*R[0]*R.y + 5*R[1]*R[2]*R.z
curl_field = curl(vector_field, R)


def test_conservative():
    assert is_conservative(0) is True
    assert is_conservative(R.x) is True
    assert is_conservative(2 * R.x + 3 * R.y + 4 * R.z) is True
    assert is_conservative(R[1]*R[2]*R.x + R[0]*R[2]*R.y + R[0]*R[1]*R.z) is \
           True
    assert is_conservative(R[0] * R.y) is False
    assert is_conservative(grad_field) is True
    assert is_conservative(curl_field) is False
    assert is_conservative(4*R[0]*R[1]*R[2]*R.x + 2*R[0]**2*R[2]*R.y) is \
                           False


def test_solenoidal():
    assert is_solenoidal(0) is True
    assert is_solenidal(R.x) is True
    assert is_solenoidal(2 * R.x + 3 * R.y + 4 * R.z) is True
    assert is_solenoidal(R[1]*R[2]*R.x + R[0]*R[2]*R.y + R[0]*R[1]*R.z) is \
           True
    assert is_solenoidal(R[1] * R.y) is False
    assert is_solenoidal(grad_field) is False
    assert is_solenoidal(curl_field) is True
    assert is_solenoidal((-2*R[1] + 3)*R.z) is False


def test_scalar_potential():
    assert scalar_potential(0, R) == 0
    assert scalar_potential(R.x, R) == R[0]
    assert scalar_potential(R.y, R) == R[1]
    assert scalar_potential(R.z, R) == R[2]
    assert scalar_potential(R[1]*R[2]*R.x + R[0]*R[2]*R.y + \
                            R[0]*R[1]*R.z, R) == R[0]*R[1]*R[2]
    assert scalar_potential(grad_field, R) == scalar_field
    assert raises(ValueError, scalar_potential(curl_field))


def test_scalarpotential_difference():
    origin = 0
    point1 = 1*R.x + 2*R.y + 3*R.z
    point2 = 4*R.x + 5*R.y + 6*R.z
    pos_vect = R[0]*R.x + R[1]*R.y + R[2]*R.z
    assert scalar_potential_difference(0, R, point1, point2) == 0
    assert scalar_potential_difference(scalar_field, R, 0, pos_vect) == \
                                              scalar_field
    assert scalar_potential_difference(grad_field, R, 0, pos_vect) == \
                                              scalar_field
    assert scalar_potential_difference(grad_field, R, point1, point2) == \
                                              948
    assert scalar_potential_difference(R[1]*R[2]*R.x + R[0]*R[2]*R.y + \
                                       R[0]*R[1]*R.z, R, point1, pos_vect) == \
                                       R[0]*R[1]*R[2] - 6
    assert raises(ValueError, scalar_potential_difference(curl_field, \
                                                          R, point1, point2))
