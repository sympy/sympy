"""Tests for kinematics.py"""

from sympy.utilities.pytest import raises
#from sympy.physics.classical.kinematics import *
from kinematics import *

q1, q2, q3, q4, q5 = symbols('q1 q2 q3 q4 q5')
N = ReferenceFrame('N')
A = N.orientnew('A', 'Simple', q1, 3)
B = A.orientnew('B', 'Simple', q2, 1)
C = B.orientnew('C', 'Simple', q3, 2)


def test_dot_cross():
    assert dot(A.x, A.x) == 1
    assert dot(A.x, A.y) == 0
    assert dot(A.x, A.z) == 0

    assert dot(A.y, A.x) == 0
    assert dot(A.y, A.y) == 1
    assert dot(A.y, A.z) == 0

    assert dot(A.z, A.x) == 0
    assert dot(A.z, A.y) == 0
    assert dot(A.z, A.z) == 1

    assert cross(A.x, A.x) == 0
    assert cross(A.x, A.y) == A.z
    assert cross(A.x, A.z) == -A.y

    assert cross(A.y, A.x) == -A.z
    assert cross(A.y, A.y) == 0
    assert cross(A.y, A.z) == A.x

    assert cross(A.z, A.x) == A.y
    assert cross(A.z, A.y) == -A.x
    assert cross(A.z, A.z) == 0


def test_ReferenceFrame():
    A = ReferenceFrame('A')
    phi = Symbol("phi")
    B = A.orientnew('B', 'Simple', phi, 1)
    assert B.parent is not None
    assert B.parent_orient is not None
    B = A.orientnew('B', 'Simple', phi, 2)
    assert B.parent is not None
    assert B.parent_orient is not None
    B = A.orientnew('B', 'Simple', phi, 3)
    assert B.parent is not None
    assert B.parent_orient is not None

def test_cross_different_frames1():
    assert cross(N.x, A.x) == sin(q1)*A.z
    assert cross(N.x, A.y) == cos(q1)*A.z
    assert cross(N.x, A.z) == -sin(q1)*A.x-cos(q1)*A.y
    assert cross(N.y, A.x) == -cos(q1)*A.z
    assert cross(N.y, A.y) == sin(q1)*A.z
    assert cross(N.y, A.z) == cos(q1)*A.x-sin(q1)*A.y
    assert cross(N.z, A.x) == A.y
    assert cross(N.z, A.y) == -A.x
    assert cross(N.z, A.z) == 0

def test_cross_method():
    q1, q2, q3 = symbols('q1 q2 q3')
    N = ReferenceFrame('N')
    A = N.orientnew('A', 'Simple', q1, 1)
    B = N.orientnew('B', 'Simple', q2, 2)
    C = N.orientnew('C', 'Simple', q3, 3)
    assert cross(N.x, N.x) == 0
    assert cross(N.x, N.y) == N.z
    assert N.x.cross(N.z) == -N.y

    assert N.y.cross(N.x) == -N.z
    assert N.y.cross(N.y) == 0
    assert N.y.cross(N.z) == N.x

    assert N.z.cross(N.x) == N.y
    assert N.z.cross(N.y) == -N.x
    assert N.z.cross(N.z) == 0

    assert N.x.cross(A.x) == 0
    assert N.x.cross(A.y) == A.z
    assert N.x.cross(A.z) == -A.y

    assert N.y.cross(A.x) == -N.z
    assert N.y.cross(A.y) == sin(q1)*N.x
    assert N.y.cross(A.z) == cos(q1)*N.x

    assert N.x.cross(B.x) == sin(q2)*N.y
    assert N.x.cross(B.y) == N.z
    assert N.x.cross(B.z) == -cos(q2)*N.y

def test_cross_different_frames2():
    assert cross(N.x, A.x) == sin(q1)*A.z
    assert cross(N.x, A.y) == cos(q1)*A.z
    assert cross(N.x, A.x + A.y) == sin(q1)*A.z + cos(q1)*A.z
    assert cross(A.x + A.y, N.x) == -sin(q1)*A.z - cos(q1)*A.z


def test_cross_different_frames3():
    assert cross(A.x, C.x) == sin(q3)*C.y
    assert cross(A.x, C.y) == -sin(q3)*C.x + cos(q3)*C.z
    assert cross(A.x, C.z) == -cos(q3)*C.y
    assert cross(C.x, A.x) == -sin(q3)*C.y
    assert cross(C.y, A.x) == sin(q3)*C.x - cos(q3)*C.z
    assert cross(C.z, A.x) == cos(q3)*C.y

def test_express1():
    assert express(A.x, C) == cos(q3)*C.x + sin(q3)*C.z
    assert express(A.y, C) == sin(q2)*sin(q3)*C.x + cos(q2)*C.y - \
            sin(q2)*cos(q3)*C.z
    assert express(A.z, C) == -sin(q3)*cos(q2)*C.x + sin(q2)*C.y + \
            cos(q2)*cos(q3)*C.z

def test_express2():
    assert A.x.express(N) == cos(q1)*N.x + sin(q1)*N.y
    assert A.y.express(N) == -sin(q1)*N.x + cos(q1)*N.y
    assert A.z.express(N) == N.z
    assert A.x.express(A) == A.x
    assert A.y.express(A) == A.y
    assert A.z.express(A) == A.z
    assert A.x.express(B) == B.x
    assert A.y.express(B) == cos(q2)*B.y - sin(q2)*B.z
    assert A.z.express(B) == sin(q2)*B.y + cos(q2)*B.z
    assert A.x.express(C) == cos(q3)*C.x + sin(q3)*C.z
    assert A.y.express(C) == sin(q2)*sin(q3)*C.x + cos(q2)*C.y - \
            sin(q2)*cos(q3)*C.z
    assert A.z.express(C) == -sin(q3)*cos(q2)*C.x + sin(q2)*C.y + \
            cos(q2)*cos(q3)*C.z

def test_express3():
    # Check to make sure UnitVectors get converted properly
    assert express(N.x, N) == N.x
    assert express(N.y, N) == N.y
    assert express(N.z, N) == N.z
    assert express(N.x, A) == (cos(q1)*A.x - sin(q1)*A.y)
    assert express(N.y, A) == (sin(q1)*A.x + cos(q1)*A.y)
    assert express(N.z, A) == A.z
    assert express(N.x, B) == (cos(q1)*B.x - sin(q1)*cos(q2)*B.y + \
            sin(q1)*sin(q2)*B.z)
    assert express(N.y, B) == (sin(q1)*B.x + cos(q1)*cos(q2)*B.y - \
            sin(q2)*cos(q1)*B.z)
    assert express(N.z, B) == (sin(q2)*B.y + cos(q2)*B.z)
    assert express(N.x, C) == (
            (cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*C.x -
            sin(q1)*cos(q2)*C.y +
            (sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))*C.z)
    assert express(N.y, C) == (
            (sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1))*C.x +
            cos(q1)*cos(q2)*C.y +
            (sin(q1)*sin(q3) - sin(q2)*cos(q1)*cos(q3))*C.z)
    assert express(N.z, C) == (-sin(q3)*cos(q2)*C.x + sin(q2)*C.y +
            cos(q2)*cos(q3)*C.z)

    assert express(A.x, N) == (cos(q1)*N.x + sin(q1)*N.y)
    assert express(A.y, N) == (-sin(q1)*N.x + cos(q1)*N.y)
    assert express(A.z, N) == N.z
    assert express(A.x, A) == A.x
    assert express(A.y, A) == A.y
    assert express(A.z, A) == A.z
    assert express(A.x, B) == B.x
    assert express(A.y, B) == (cos(q2)*B.y - sin(q2)*B.z)
    assert express(A.z, B) == (sin(q2)*B.y + cos(q2)*B.z)
    assert express(A.x, C) == (cos(q3)*C.x + sin(q3)*C.z)
    assert express(A.y, C) == (sin(q2)*sin(q3)*C.x + cos(q2)*C.y -
            sin(q2)*cos(q3)*C.z)
    assert express(A.z, C) == (-sin(q3)*cos(q2)*C.x + sin(q2)*C.y +
            cos(q2)*cos(q3)*C.z)

    assert express(B.x, N) == (cos(q1)*N.x + sin(q1)*N.y)
    assert express(B.y, N) == (-sin(q1)*cos(q2)*N.x +
            cos(q1)*cos(q2)*N.y + sin(q2)*N.z)
    assert express(B.z, N) == (sin(q1)*sin(q2)*N.x -
            sin(q2)*cos(q1)*N.y + cos(q2)*N.z)
    assert express(B.x, A) == A.x
    assert express(B.y, A) == (cos(q2)*A.y + sin(q2)*A.z)
    assert express(B.z, A) == (-sin(q2)*A.y + cos(q2)*A.z)
    assert express(B.x, B) == B.x
    assert express(B.y, B) == B.y
    assert express(B.z, B) == B.z
    assert express(B.x, C) == (cos(q3)*C.x + sin(q3)*C.z)
    assert express(B.y, C) == C.y
    assert express(B.z, C) == (-sin(q3)*C.x + cos(q3)*C.z)

    assert express(C.x, N) == (
            (cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*N.x +
            (sin(q1)*cos(q3)+sin(q2)*sin(q3)*cos(q1))*N.y -
                sin(q3)*cos(q2)*N.z)
    assert express(C.y, N) == (
            -sin(q1)*cos(q2)*N.x + cos(q1)*cos(q2)*N.y + sin(q2)*N.z)
    assert express(C.z, N) == (
            (sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))*N.x +
            (sin(q1)*sin(q3)-sin(q2)*cos(q1)*cos(q3))*N.y +
            cos(q2)*cos(q3)*N.z)
    assert express(C.x, A) == (cos(q3)*A.x + sin(q2)*sin(q3)*A.y -
            sin(q3)*cos(q2)*A.z)
    assert express(C.y, A) == (cos(q2)*A.y + sin(q2)*A.z)
    assert express(C.z, A) == (sin(q3)*A.x - sin(q2)*cos(q3)*A.y +
            cos(q2)*cos(q3)*A.z)
    assert express(C.x, B) == (cos(q3)*B.x - sin(q3)*B.z)
    assert express(C.y, B) == B.y
    assert express(C.z, B) == (sin(q3)*B.x + cos(q3)*B.z)
    assert express(C.x, C) == C.x
    assert express(C.y, C) == C.y
    assert express(C.z, C) == C.z == (C.z)

    #  Check to make sure Vectors get converted back to UnitVectors
    assert N.x == express((cos(q1)*A.x - sin(q1)*A.y), N)
    assert N.y == express((sin(q1)*A.x + cos(q1)*A.y), N)
    assert N.x == express((cos(q1)*B.x - sin(q1)*cos(q2)*B.y +
            sin(q1)*sin(q2)*B.z), N)
    assert N.y == express((sin(q1)*B.x + cos(q1)*cos(q2)*B.y -
        sin(q2)*cos(q1)*B.z), N)
    assert N.z == express((sin(q2)*B.y + cos(q2)*B.z), N)


    assert N.x == express((
            (cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*C.x -
            sin(q1)*cos(q2)*C.y +
            (sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))*C.z), N)
    assert N.y == express((
            (sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1))*C.x +
            cos(q1)*cos(q2)*C.y +
            (sin(q1)*sin(q3) - sin(q2)*cos(q1)*cos(q3))*C.z), N)
    assert N.z == express((-sin(q3)*cos(q2)*C.x + sin(q2)*C.y +
            cos(q2)*cos(q3)*C.z), N)

    assert A.x == express((cos(q1)*N.x + sin(q1)*N.y), A)
    assert A.y == express((-sin(q1)*N.x + cos(q1)*N.y), A)

    assert A.y == express((cos(q2)*B.y - sin(q2)*B.z), A)
    assert A.z == express((sin(q2)*B.y + cos(q2)*B.z), A)

    assert A.x == express((cos(q3)*C.x + sin(q3)*C.z), A)

    # Tripsimp messes up here too.
    #print express((sin(q2)*sin(q3)*C.x + cos(q2)*C.y -
    #        sin(q2)*cos(q3)*C.z), A)
    assert A.y == express((sin(q2)*sin(q3)*C.x + cos(q2)*C.y -
            sin(q2)*cos(q3)*C.z), A)

    assert A.z == express((-sin(q3)*cos(q2)*C.x + sin(q2)*C.y +
            cos(q2)*cos(q3)*C.z), A)
    assert B.x == express((cos(q1)*N.x + sin(q1)*N.y), B)
    assert B.y == express((-sin(q1)*cos(q2)*N.x +
            cos(q1)*cos(q2)*N.y + sin(q2)*N.z), B)

    assert B.z == express((sin(q1)*sin(q2)*N.x -
            sin(q2)*cos(q1)*N.y + cos(q2)*N.z), B)

    assert B.y == express((cos(q2)*A.y + sin(q2)*A.z), B)
    assert B.z == express((-sin(q2)*A.y + cos(q2)*A.z), B)
    assert B.x == express((cos(q3)*C.x + sin(q3)*C.z), B)
    assert B.z == express((-sin(q3)*C.x + cos(q3)*C.z), B)

    assert C.x == express((
            (cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*N.x +
            (sin(q1)*cos(q3)+sin(q2)*sin(q3)*cos(q1))*N.y -
                sin(q3)*cos(q2)*N.z), C)
    assert C.y == express((
            -sin(q1)*cos(q2)*N.x + cos(q1)*cos(q2)*N.y + sin(q2)*N.z), C)
    assert C.z == express((
            (sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))*N.x +
            (sin(q1)*sin(q3)-sin(q2)*cos(q1)*cos(q3))*N.y +
            cos(q2)*cos(q3)*N.z), C)
    assert C.x == express((cos(q3)*A.x + sin(q2)*sin(q3)*A.y -
            sin(q3)*cos(q2)*A.z), C)
    assert C.y == express((cos(q2)*A.y + sin(q2)*A.z), C)
    assert C.z == express((sin(q3)*A.x - sin(q2)*cos(q3)*A.y +
            cos(q2)*cos(q3)*A.z), C)
    assert C.x == express((cos(q3)*B.x - sin(q3)*B.z), C)
    assert C.z == express((sin(q3)*B.x + cos(q3)*B.z), C)

def test_ang_vel():
    A2 = N.rotate('A2', 2, q4)
    assert N.ang_vel(N) == 0
    assert N.ang_vel(A) == -q1p*N.z
    assert N.ang_vel(B) == -q1p*A.z - q2p*B.x
    assert N.ang_vel(C) == -q1p*A.z - q2p*B.x - q3p*B.y
    assert N.ang_vel(A2) == -q4p*N.y

    assert A.ang_vel(N) == q1p*N.z
    assert A.ang_vel(A) == 0
    assert A.ang_vel(B) == - q2p*B.x
    assert A.ang_vel(C) == - q2p*B.x - q3p*B.y
    assert A.ang_vel(A2) == q1p*N.z - q4p*N.y

    assert B.ang_vel(N) == q1p*A.z + q2p*A.x
    assert B.ang_vel(A) == q2p*A.x
    assert B.ang_vel(B) == 0
    assert B.ang_vel(C) == -q3p*B.y
    assert B.ang_vel(A2) == q1p*A.z + q2p*A.x - q4p*N.y

    assert C.ang_vel(N) == q1p*A.z + q2p*A.x + q3p*B.y
    assert C.ang_vel(A) == q2p*A.x + q3p*C.y
    assert C.ang_vel(B) == q3p*B.y
    assert C.ang_vel(C) == 0
    assert C.ang_vel(A2) == q1p*A.z + q2p*A.x + q3p*B.y - q4p*N.y

    assert A2.ang_vel(N) == q4p*A2.y
    assert A2.ang_vel(A) == q4p*A2.y - q1p*N.z
    assert A2.ang_vel(B) == q4p*N.y - q1p*A.z - q2p*A.x
    assert A2.ang_vel(C) == q4p*N.y - q1p*A.z - q2p*A.x - q3p*B.y
    assert A2.ang_vel(A2) == 0


"""
class Vector(object):
    def __init__(self, inlist):
    def __str__(self):
    def __repr__(self):
    def __add__(self, other):
    def __and__(self, other):
    def __div__(self, other):
    def __eq__(self, other):
    def __mul__(self, other):
    def __neg__(self):
    def __rmul__(self, other):
    def __sub__(self, other):
    def __xor__(self, other):
    def dot(self, other):
    def cross(self, other):
    def diff(self, wrt, otherframe):
    def dt(self, otherframe):
    def express(self, otherframe):
    def mag(self):
    def unit(self):

class ReferenceFrame(object):
    def __init__(self, name=''):
    def __iter__(self):
    def __str__(self):
    def __repr__(self):
    def next(self):
    def _common_frame(self,other):
    def ang_vel(self, other):
    def dcm(self, other):
    def orientnew(self, newname, rot_type, amounts, rot_order):
    def orient(self, parent, rot_type, amounts, rot_order):
    def set_ang_vel(self, value, other):
    def x(self):
    def y(self):
    def z(self):
"""
