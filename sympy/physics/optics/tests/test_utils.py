from sympy.physics.optics.utils import refraction_angle
from sympy.physics.optics.medium import Medium
from sympy.physics.optics.waves import TWave

from sympy import symbols, sqrt, Matrix
from sympy.geometry.point3d import Point3D
from sympy.geometry.line3d import Ray3D
from sympy.geometry.plane import Plane


def test_refraction_angle():
    n1, n2 = symbols('n1, n2')
    m1 = Medium('m1')
    m2 = Medium('m2')

    r1 = Ray3D(Point3D(-1, -1, 1), Point3D(0, 0, 0))
    n = Matrix([0, 0, 1])
    P = Plane(Point3D(0, 0, 0), normal_vector=[0, 0, 1])
    assert refraction_angle(r1, 1, 1, n) == Matrix([
                                            [ 1],
                                            [ 1],
                                            [-1]])
    assert refraction_angle(r1, 1, 1, plane=P) ==\
        Ray3D(Point3D(0, 0, 0), Point3D(1, 1, -1))
