from sympy.core.function import Derivative
from sympy.vector.vector import Vector
from sympy.vector.coordsysrect import CoordSysCartesian
from sympy.simplify import simplify
from sympy.core.symbol import symbols
from sympy.core import S
from sympy import sin, cos
from sympy.vector.functions import (curl, divergence, gradient,
                                    is_conservative, is_solenoidal,
                                    scalar_potential,
                                    scalar_potential_difference)
from sympy.utilities.pytest import raises

C = CoordSysCartesian('C')
i, j, k = C.base_vectors()
x, y, z = C.base_scalars()
delop = C.delop
a, b, c, q = symbols('a b c q')

def test_del_operator():

    #Tests for curl
    assert (delop ^ Vector.zero ==
            (Derivative(0, C.y) - Derivative(0, C.z))*C.i +
            (-Derivative(0, C.x) + Derivative(0, C.z))*C.j +
            (Derivative(0, C.x) - Derivative(0, C.y))*C.k)
    assert ((delop ^ Vector.zero).doit() == Vector.zero ==
            curl(Vector.zero, C))
    assert delop.cross(Vector.zero) == delop ^ Vector.zero
    assert (delop ^ i).doit() == Vector.zero
    assert delop.cross(2*y**2*j, doit = True) == Vector.zero
    assert delop.cross(2*y**2*j) == delop ^ 2*y**2*j
    v = x*y*z * (i + j + k)
    assert ((delop ^ v).doit() ==
            (-x*y + x*z)*i + (x*y - y*z)*j + (-x*z + y*z)*k ==
            curl(v, C))
    assert delop ^ v == delop.cross(v)
    assert (delop.cross(2*x**2*j) ==
            (Derivative(0, C.y) - Derivative(2*C.x**2, C.z))*C.i +
            (-Derivative(0, C.x) + Derivative(0, C.z))*C.j +
            (-Derivative(0, C.y) + Derivative(2*C.x**2, C.x))*C.k)
    assert (delop.cross(2*x**2*j, doit = True) == 4*x*k ==
            curl(2*x**2*j, C))

    #Tests for divergence
    assert delop & Vector.zero == S(0) == divergence(Vector.zero, C)
    assert (delop & Vector.zero).doit() == S(0)
    assert delop.dot(Vector.zero) == delop & Vector.zero
    assert (delop & i).doit() == S(0)
    assert (delop & x**2*i).doit() == 2*x == divergence(x**2*i, C)
    assert (delop.dot(v, doit = True) == x*y + y*z + z*x ==
            divergence(v, C))
    assert delop & v == delop.dot(v)
    assert delop.dot(1/(x*y*z) * (i + j + k), doit = True) == \
           - 1 / (x*y*z**2) - 1 / (x*y**2*z) - 1 / (x**2*y*z)
    v = x*i + y*j + z*k
    assert (delop & v == Derivative(C.x, C.x) +
            Derivative(C.y, C.y) + Derivative(C.z, C.z))
    assert delop.dot(v, doit = True) == 3 == divergence(v, C)
    assert delop & v == delop.dot(v)
    assert simplify((delop & v).doit()) == 3

    #Tests for gradient
    assert (delop.gradient(0, doit = True) == Vector.zero ==
            gradient(0, C))
    assert delop.gradient(0) == delop(0)
    assert (delop(S(0))).doit() == Vector.zero
    assert (delop(x) == (Derivative(C.x, C.x))*C.i +
            (Derivative(C.x, C.y))*C.j + (Derivative(C.x, C.z))*C.k)
    assert (delop(x)).doit() == i == gradient(x, C)
    assert (delop(x*y*z) ==
            (Derivative(C.x*C.y*C.z, C.x))*C.i +
            (Derivative(C.x*C.y*C.z, C.y))*C.j +
            (Derivative(C.x*C.y*C.z, C.z))*C.k)
    assert (delop.gradient(x*y*z, doit = True) ==
            y*z*i + z*x*j + x*y*k ==
            gradient(x*y*z, C))
    assert delop(x*y*z) == delop.gradient(x*y*z)
    assert (delop(2*x**2)).doit() == 4*x*i
    assert ((delop(a*sin(y) / x)).doit() ==
            -a*sin(y)/x**2 * i + a*cos(y)/x * j)

    #Tests for directional derivative
    assert (Vector.zero & delop)(a) == S(0)
    assert ((Vector.zero & delop)(a)).doit() == S(0)
    assert ((v & delop)(Vector.zero)).doit() == Vector.zero
    assert ((v & delop)(S(0))).doit() == S(0)
    assert ((i & delop)(x)).doit() == 1
    assert ((j & delop)(y)).doit() == 1
    assert ((k & delop)(z)).doit() == 1
    assert ((i & delop)(x*y*z)).doit() == y*z
    assert ((v & delop)(x)).doit() == x
    assert ((v & delop)(x*y*z)).doit() == 3*x*y*z
    assert (v & delop)(x + y + z) == C.x + C.y + C.z
    assert ((v & delop)(x + y + z)).doit() == x + y + z
    assert ((v & delop)(v)).doit() == v
    assert ((i & delop)(v)).doit() == i
    assert ((j & delop)(v)).doit() == j
    assert ((k & delop)(v)).doit() == k
    assert ((v & delop)(Vector.zero)).doit() == Vector.zero


def test_product_rules():
    """
    Tests the six product rules defined with respect to the Del
    operator

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Del

    """

    #Define the scalar and vector functions
    f = 2*x*y*z
    g = x*y + y*z + z*x
    u = x**2*i + 4*j - y**2*z*k
    v = 4*i + x*y*z*k

    #First product rule
    lhs = delop(f * g, doit = True)
    rhs = (f * delop(g) + g * delop(f)).doit()
    assert simplify(lhs) == simplify(rhs)

    #Second product rule
    lhs = delop(u & v).doit()
    rhs = ((u ^ (delop ^ v)) + (v ^ (delop ^ u)) + \
          ((u & delop)(v)) + ((v & delop)(u))).doit()
    assert simplify(lhs) == simplify(rhs)

    #Third product rule
    lhs = (delop & (f*v)).doit()
    rhs = ((f * (delop & v)) + (v & (delop(f)))).doit()
    assert simplify(lhs) == simplify(rhs)

    #Fourth product rule
    lhs = (delop & (u ^ v)).doit()
    rhs = ((v & (delop ^ u)) - (u & (delop ^ v))).doit()
    assert simplify(lhs) == simplify(rhs)

    #Fifth product rule
    lhs = (delop ^ (f * v)).doit()
    rhs = (((delop(f)) ^ v) + (f * (delop ^ v))).doit()
    assert simplify(lhs) == simplify(rhs)

    #Sixth product rule
    lhs = (delop ^ (u ^ v)).doit()
    rhs = ((u * (delop & v) - v * (delop & u) +
           (v & delop)(u) - (u & delop)(v))).doit()
    assert simplify(lhs) == simplify(rhs)


P = C.orient_new_axis('P', q, C.k)
scalar_field = 2*x**2*y*z
grad_field = gradient(scalar_field, C)
vector_field = y**2*i + 3*x*j + 5*y*z*k
curl_field = curl(vector_field, C)


def test_conservative():
    assert is_conservative(Vector.zero) is True
    assert is_conservative(i) is True
    assert is_conservative(2 * i + 3 * j + 4 * k) is True
    assert (is_conservative(y*z*i + x*z*j + x*y*k) is
            True)
    assert is_conservative(x * j) is False
    assert is_conservative(grad_field) is True
    assert is_conservative(curl_field) is False
    assert (is_conservative(4*x*y*z*i + 2*x**2*z*j) is
            False)
    assert is_conservative(z*P.i + P.x*k) is True


def test_solenoidal():
    assert is_solenoidal(Vector.zero) is True
    assert is_solenoidal(i) is True
    assert is_solenoidal(2 * i + 3 * j + 4 * k) is True
    assert (is_solenoidal(y*z*i + x*z*j + x*y*k) is
            True)
    assert is_solenoidal(y * j) is False
    assert is_solenoidal(grad_field) is False
    assert is_solenoidal(curl_field) is True
    assert is_solenoidal((-2*y + 3)*k) is True
    assert is_solenoidal(cos(q)*i + sin(q)*j + cos(q)*P.k) is True
    assert is_solenoidal(z*P.i + P.x*k) is True


def test_scalar_potential():
    assert scalar_potential(Vector.zero, C) == 0
    assert scalar_potential(i, C) == x
    assert scalar_potential(j, C) == y
    assert scalar_potential(k, C) == z
    assert scalar_potential(y*z*i + x*z*j + x*y*k, C) == x*y*z
    assert scalar_potential(grad_field, C) == scalar_field
    assert scalar_potential(z*P.i + P.x*k, C) == x*z*cos(q) + y*z*sin(q)
    assert scalar_potential(z*P.i + P.x*k, P) == P.x*P.z
    raises(ValueError, lambda: scalar_potential(x*j, C))


def test_scalar_potential_difference():
    point1 = C.origin.locate_new('P1', 1*i + 2*j + 3*k)
    point2 = C.origin.locate_new('P2', 4*i + 5*j + 6*k)
    genericpointC = C.origin.locate_new('RP', x*i + y*j + z*k)
    genericpointP = P.origin.locate_new('PP', P.x*P.i + P.y*P.j + P.z*P.k)
    assert scalar_potential_difference(S(0), C, point1, point2) == 0
    assert (scalar_potential_difference(scalar_field, C, C.origin,
                                        genericpointC) ==
            scalar_field)
    assert (scalar_potential_difference(grad_field, C, C.origin,
                                        genericpointC) ==
            scalar_field)
    assert scalar_potential_difference(grad_field, C, point1, point2) == 948
    assert (scalar_potential_difference(y*z*i + x*z*j +
                                        x*y*k, C, point1,
                                        genericpointC) ==
            x*y*z - 6)
    potential_diff_P = (2*P.z*(P.x*sin(q) + P.y*cos(q))*
                        (P.x*cos(q) - P.y*sin(q))**2)
    assert (scalar_potential_difference(grad_field, P, P.origin,
                                        genericpointP).simplify() ==
            potential_diff_P)
