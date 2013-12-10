from sympy import S, Symbol, sin, cos
from sympy.physics.mechanics import ReferenceFrame, Vector
from sympy.physics.em import divergence, gradient, curl, is_conservative, \
     is_solenoidal, scalar_potential, scalar_potential_difference
from sympy.utilities.pytest import raises

R = ReferenceFrame('R')

def test_curl():
    assert curl(Vector(0), R) == Vector(0)
    assert curl(R.x, R) == Vector(0)
    assert curl(2*R[1]**2*R.y, R) == Vector(0)
    assert curl(R[0]*R[1]*R.z, R) == R[0]*R.x - R[1]*R.y
    assert curl(R[0]*R[1]*R[2] * (R.x+R.y+R.z), R) == \
           (-R[0]*R[1] + R[0]*R[2])*R.x + (R[0]*R[1] - R[1]*R[2])*R.y + \
           (-R[0]*R[2] + R[1]*R[2])*R.z
    assert curl(2*R[0]**2*R.y, R) == 4*R[0]*R.z


def test_divergence():
    assert divergence(Vector(0), R) == S(0)
    assert divergence(R.x, R) == S(0)
    assert divergence(R[0]**2*R.x, R) == 2*R[0]
    assert divergence(R[0]*R[1]*R[2] * (R.x+R.y+R.z), R) == \
           R[0]*R[1] + R[0]*R[2] + R[1]*R[2]
    assert divergence((1/(R[0]*R[1]*R[2])) * (R.x+R.y+R.z), R) == \
           -1/(R[0]*R[1]*R[2]**2) - 1/(R[0]*R[1]**2*R[2]) - 1/(R[0]**2*R[1]*R[2])


def test_gradient():
    a = Symbol('a')
    assert gradient(0, R) == Vector(0)
    assert gradient(R[0], R) == R.x
    assert gradient(R[0]*R[1]*R[2], R) == \
           R[1]*R[2]*R.x + R[0]*R[2]*R.y + R[0]*R[1]*R.z
    assert gradient(2*R[0]**2, R) == 4*R[0]*R.x
    assert gradient(a*sin(R[1])/R[0], R) == \
           - a*sin(R[1])/R[0]**2*R.x + a*cos(R[1])/R[0]*R.y


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
    assert is_solenoidal(R.x) is True
    assert is_solenoidal(2 * R.x + 3 * R.y + 4 * R.z) is True
    assert is_solenoidal(R[1]*R[2]*R.x + R[0]*R[2]*R.y + R[0]*R[1]*R.z) is \
           True
    assert is_solenoidal(R[1] * R.y) is False
    assert is_solenoidal(grad_field) is False
    assert is_solenoidal(curl_field) is True
    assert is_solenoidal((-2*R[1] + 3)*R.z) is True


def test_scalar_potential():
    assert scalar_potential(0, R) == 0
    assert scalar_potential(R.x, R) == R[0]
    assert scalar_potential(R.y, R) == R[1]
    assert scalar_potential(R.z, R) == R[2]
    assert scalar_potential(R[1]*R[2]*R.x + R[0]*R[2]*R.y + \
                            R[0]*R[1]*R.z, R) == R[0]*R[1]*R[2]
    assert scalar_potential(grad_field, R) == scalar_field


def test_scalar_potential_difference():
    origin = 0
    point1 = 1*R.x + 2*R.y + 3*R.z
    point2 = 4*R.x + 5*R.y + 6*R.z
    pos_vect = R[0]*R.x + R[1]*R.y + R[2]*R.z
    assert scalar_potential_difference(S(0), R, point1, point2) == 0
    assert scalar_potential_difference(scalar_field, R, Vector(0), pos_vect) == \
                                              scalar_field
    assert scalar_potential_difference(grad_field, R, Vector(0), pos_vect) == \
                                              scalar_field
    assert scalar_potential_difference(grad_field, R, point1, point2) == \
                                              948
    assert scalar_potential_difference(R[1]*R[2]*R.x + R[0]*R[2]*R.y + \
                                       R[0]*R[1]*R.z, R, point1, pos_vect) == \
                                       R[0]*R[1]*R[2] - 6
