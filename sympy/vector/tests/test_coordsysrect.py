from sympy.vector.coordsysrect import CoordSysCartesian
from sympy.vector.scalar import BaseScalar
from sympy import Symbol, sin, cos, pi, ImmutableMatrix as Matrix, \
     symbols, simplify, sqrt, zeros
from sympy.vector.functions import express
from sympy.vector.point import Point
from sympy.vector.vector import Vector


A = CoordSysCartesian('A')
a, b, c, q = symbols('a b c q')
q1, q2, q3, q4 = symbols('q1 q2 q3 q4')


def test_coordinate_vars():
    """
    Tests the coordinate variables functionality with respect to
    reorientation of coordinate systems.
    """
    assert BaseScalar('Ax', 0, A) == A.x
    assert BaseScalar('Ay', 1, A) == A.y
    assert BaseScalar('Az', 2, A) == A.z
    assert BaseScalar('Ax', 0, A).__hash__() == A.x.__hash__()
    assert isinstance(A.x, BaseScalar) and \
           isinstance(A.y, BaseScalar) and \
           isinstance(A.z, BaseScalar)
    assert A.scalar_map(A) == {A.x: A.x, A.y: A.y, A.z: A.z}
    assert A.x.system == A
    B = A.orient_new('B', 'Axis', [q, A.k])
    assert B.scalar_map(A) == {B.z: A.z, B.y: -A.x*sin(q) + A.y*cos(q),
                                 B.x: A.x*cos(q) + A.y*sin(q)}
    assert A.scalar_map(B) == {A.x: B.x*cos(q) - B.y*sin(q),
                                 A.y: B.x*sin(q) + B.y*cos(q), A.z: B.z}
    assert express(B.x, A, variables=True) == A.x*cos(q) + A.y*sin(q)
    assert express(B.y, A, variables=True) == -A.x*sin(q) + A.y*cos(q)
    assert express(B.z, A, variables=True) == A.z
    assert express(B.x*B.y*B.z, A, variables=True) == \
           A.z*(-A.x*sin(q) + A.y*cos(q))*(A.x*cos(q) + A.y*sin(q))
    assert express(B.x*B.i + B.y*B.j + B.z*B.k, A) == \
           (B.x*cos(q) - B.y*sin(q))*A.i + (B.x*sin(q) + \
           B.y*cos(q))*A.j + B.z*A.k
    assert simplify(express(B.x*B.i + B.y*B.j + B.z*B.k, A, \
                            variables=True)) == \
           A.x*A.i + A.y*A.j + A.z*A.k
    assert express(A.x*A.i + A.y*A.j + A.z*A.k, B) == \
           (A.x*cos(q) + A.y*sin(q))*B.i + \
           (-A.x*sin(q) + A.y*cos(q))*B.j + A.z*B.k
    assert simplify(express(A.x*A.i + A.y*A.j + A.z*A.k, B, \
                            variables=True)) == \
           B.x*B.i + B.y*B.j + B.z*B.k
    N = B.orient_new('N', 'Axis', [-q, B.k])
    assert N.scalar_map(A) == \
           {N.x: A.x, N.z: A.z, N.y: A.y}
    C = A.orient_new('C', 'Axis', [q, A.i + A.j + A.k])
    mapping = A.scalar_map(C)
    assert mapping[A.x] == 2*C.x*cos(q)/3 + C.x/3 - \
           2*C.y*sin(q + pi/6)/3 + C.y/3 - 2*C.z*cos(q + pi/3)/3 + C.z/3
    assert mapping[A.y] == -2*C.x*cos(q + pi/3)/3 + \
           C.x/3 + 2*C.y*cos(q)/3 + C.y/3 - 2*C.z*sin(q + pi/6)/3 + C.z/3
    assert mapping[A.z] == -2*C.x*sin(q + pi/6)/3 + C.x/3 - \
           2*C.y*cos(q + pi/3)/3 + C.y/3 + 2*C.z*cos(q)/3 + C.z/3


def test_rotation_matrix():
    N = CoordSysCartesian('N')
    A = N.orient_new('A', 'Axis', [q1, N.k])
    B = A.orient_new('B', 'Axis', [q2, A.i])
    C = B.orient_new('C', 'Axis', [q3, B.j])
    D = N.orient_new('D', 'Axis', [q4, N.j])
    E = N.orient_new('E', 'Space', [q1, q2, q3], '123')
    assert N.rotation_matrix(C) == Matrix([
        [- sin(q1) * sin(q2) * sin(q3) + cos(q1) * cos(q3), - sin(q1) *
        cos(q2), sin(q1) * sin(q2) * cos(q3) + sin(q3) * cos(q1)], \
        [sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1), \
         cos(q1) * cos(q2), sin(q1) * sin(q3) - sin(q2) * cos(q1) * \
         cos(q3)], [- sin(q3) * cos(q2), sin(q2), cos(q2) * cos(q3)]])
    test_mat = D.rotation_matrix(C) - Matrix(
        [[cos(q1) * cos(q3) * cos(q4) - sin(q3) * (- sin(q4) * cos(q2) +
        sin(q1) * sin(q2) * cos(q4)), - sin(q2) * sin(q4) - sin(q1) *
            cos(q2) * cos(q4), sin(q3) * cos(q1) * cos(q4) + cos(q3) * \
          (- sin(q4) * cos(q2) + sin(q1) * sin(q2) * cos(q4))], \
         [sin(q1) * cos(q3) + sin(q2) * sin(q3) * cos(q1), cos(q1) * \
          cos(q2), sin(q1) * sin(q3) - sin(q2) * cos(q1) * cos(q3)], \
         [sin(q4) * cos(q1) * cos(q3) - sin(q3) * (cos(q2) * cos(q4) + \
                                                   sin(q1) * sin(q2) * \
                                                   sin(q4)), sin(q2) *
                cos(q4) - sin(q1) * sin(q4) * cos(q2), sin(q3) * \
          sin(q4) * cos(q1) + cos(q3) * (cos(q2) * cos(q4) + \
                                         sin(q1) * sin(q2) * sin(q4))]])
    assert test_mat.expand() == zeros(3, 3)
    assert E.rotation_matrix(N) == Matrix(
        [[cos(q2)*cos(q3), sin(q3)*cos(q2), -sin(q2)],
        [sin(q1)*sin(q2)*cos(q3) - sin(q3)*cos(q1), \
         sin(q1)*sin(q2)*sin(q3) + cos(q1)*cos(q3), sin(q1)*cos(q2)], \
         [sin(q1)*sin(q3) + sin(q2)*cos(q1)*cos(q3), - \
          sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1), cos(q1)*cos(q2)]])


def test_vector():
    """
    Tests the effects of orientation of coordinate systems on
    basic vector operations.
    """
    N = CoordSysCartesian('N')
    A = N.orient_new('A', 'Axis', [q1, N.k])
    B = A.orient_new('B', 'Axis', [q2, A.i])
    C = B.orient_new('C', 'Axis', [q3, B.j])

    #Test to_matrix
    v1 = a*N.i + b*N.j + c*N.k
    assert v1.to_matrix(A) == Matrix([[ a*cos(q1) + b*sin(q1)],
                                      [-a*sin(q1) + b*cos(q1)],
                                      [                     c]])

    #Test dot
    assert N.i.dot(A.i) == cos(q1)
    assert N.i.dot(A.j) == -sin(q1)
    assert N.i.dot(A.k) == 0
    assert N.j.dot(A.i) == sin(q1)
    assert N.j.dot(A.j) == cos(q1)
    assert N.j.dot(A.k) == 0
    assert N.k.dot(A.i) == 0
    assert N.k.dot(A.j) == 0
    assert N.k.dot(A.k) == 1

    assert N.i.dot(A.i + A.j) == -sin(q1) + cos(q1) == \
           (A.i + A.j).dot(N.i)

    assert A.i.dot(C.i) == cos(q3)
    assert A.i.dot(C.j) == 0
    assert A.i.dot(C.k) == sin(q3)
    assert A.j.dot(C.i) == sin(q2)*sin(q3)
    assert A.j.dot(C.j) == cos(q2)
    assert A.j.dot(C.k) == -sin(q2)*cos(q3)
    assert A.k.dot(C.i) == -cos(q2)*sin(q3)
    assert A.k.dot(C.j) == sin(q2)
    assert A.k.dot(C.k) == cos(q2)*cos(q3)

    #Test cross
    assert N.i.cross(A.i) == sin(q1)*A.k
    assert N.i.cross(A.j) == cos(q1)*A.k
    assert N.i.cross(A.k) == -sin(q1)*A.i - cos(q1)*A.j
    assert N.j.cross(A.i) == -cos(q1)*A.k
    assert N.j.cross(A.j) == sin(q1)*A.k
    assert N.j.cross(A.k) == cos(q1)*A.i - sin(q1)*A.j
    assert N.k.cross(A.i) == A.j
    assert N.k.cross(A.j) == -A.i
    assert N.k.cross(A.k) == Vector.zero

    assert N.i.cross(A.i) == sin(q1)*A.k
    assert N.i.cross(A.j) == cos(q1)*A.k
    assert N.i.cross(A.i + A.j) == sin(q1)*A.k + cos(q1)*A.k
    assert (A.i + A.j).cross(N.i) == (-sin(q1) - cos(q1))*N.k

    assert A.i.cross(C.i) == sin(q3)*C.j
    assert A.i.cross(C.j) == -sin(q3)*C.i + cos(q3)*C.k
    assert A.i.cross(C.k) == -cos(q3)*C.j
    assert C.i.cross(A.i) == (-sin(q3)*cos(q2))*A.j + \
           (-sin(q2)*sin(q3))*A.k
    assert C.j.cross(A.i) == (sin(q2))*A.j + (-cos(q2))*A.k
    assert express(C.k.cross(A.i), C).trigsimp() == cos(q3)*C.j


def test_locatenew_point():
    """
    Tests Point class, and locate_new method in CoordSysCartesian.
    """
    assert isinstance(A.origin, Point)
    v = a*A.i + b*A.j + c*A.k
    C = A.locate_new('C', v)
    assert C.origin.position_wrt(A) == \
           C.origin.position_wrt(A.origin) == v
    assert A.origin.position_wrt(C) == \
           A.origin.position_wrt(C.origin) == -v
    assert A.origin.express_coordinates(C) == (-a, -b, -c)
    p = A.origin.locate_new('p', -v)
    assert p.express_coordinates(A) == (-a, -b, -c)
    assert p.position_wrt(C.origin) == p.position_wrt(C) == \
           -2 * v
    p1 = p.locate_new('p1', 2*v)
    assert p1.position_wrt(C.origin) == Vector.zero
    assert p1.express_coordinates(C) == (0, 0, 0)
    p2 = p.locate_new('p2', A.i)
    assert p1.position_wrt(p2) == 2*v - A.i
    assert p2.express_coordinates(C) == (-2*a + 1, -2*b, -2*c)
