import pytest

from sympy.core.function import expand_mul
from sympy.core.numbers import pi
from sympy.core.singleton import S
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import (cos, sin)
from sympy.core.backend import Matrix, _simplify_matrix, eye, zeros
from sympy.core.symbol import symbols
from sympy.physics.mechanics import (dynamicsymbols, Body, JointsMethod,
                                     PinJoint, PrismaticJoint, CylindricalJoint,
                                     PlanarJoint, SphericalJoint, WeldJoint)
from sympy.physics.mechanics.joint import Joint
from sympy.physics.vector import Vector, ReferenceFrame, Point
from sympy.testing.pytest import raises, warns_deprecated_sympy


Vector.simp = True
t = dynamicsymbols._t  # type: ignore


def _generate_body(interframe=False):
    N = ReferenceFrame('N')
    A = ReferenceFrame('A')
    P = Body('P', frame=N)
    C = Body('C', frame=A)
    if interframe:
        Pint, Cint = ReferenceFrame('P_int'), ReferenceFrame('C_int')
        Pint.orient_axis(N, N.x, pi)
        Cint.orient_axis(A, A.y, -pi / 2)
        return N, A, P, C, Pint, Cint
    return N, A, P, C


def test_Joint():
    parent = Body('parent')
    child = Body('child')
    with pytest.raises(TypeError):
        Joint('J', parent, child)


class TestCoordinateGeneration:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.q, self.u, self.qj, self.uj = dynamicsymbols('q u q_J u_J')
        self.q0j, self.q1j, self.q2j, self.q3j = dynamicsymbols('q0:4_J')
        self.u0j, self.u1j, self.u2j, self.u3j = dynamicsymbols('u0:4_J')
        self.q0, self.q1, self.q2, self.q3 = dynamicsymbols('q0:4')
        self.u0, self.u1, self.u2, self.u3 = dynamicsymbols('u0:4')
        _, _, self.P, self.C = _generate_body()
        # Using PinJoint to access Joint's coordinate generation method
        self.J = PinJoint('J', self.P, self.C)

    def test_single_given(self):
        assert self.J._fill_coordinate_list(self.q, 1) == Matrix([self.q])
        assert self.J._fill_coordinate_list([self.u], 1) == Matrix([self.u])
        assert self.J._fill_coordinate_list([self.u], 1, offset=2) == Matrix(
            [self.u])

    def test_none(self):
        assert self.J._fill_coordinate_list(None, 1) == Matrix([self.qj])
        assert self.J._fill_coordinate_list([None], 1) == Matrix([self.qj])
        assert self.J._fill_coordinate_list([self.q0, None], 3) == Matrix(
            [self.q0, self.q1j, self.q2j])

    def test_autofill(self):
        assert self.J._fill_coordinate_list(None, 3) == Matrix(
            [self.q0j, self.q1j, self.q2j])
        assert self.J._fill_coordinate_list([], 3) == Matrix(
            [self.q0j, self.q1j, self.q2j])

    def test_offset(self):
        assert self.J._fill_coordinate_list([], 3, offset=1) == Matrix(
            [self.q1j, self.q2j, self.q3j])
        assert self.J._fill_coordinate_list(self.q1, 3, offset=1) == Matrix(
            [self.q1, self.q2j, self.q3j])
        assert self.J._fill_coordinate_list(
            [self.q1, None, self.q3], 3, offset=1
        ) == Matrix([self.q1, self.q2j, self.q3])
        assert self.J._fill_coordinate_list(None, 2, offset=2) == Matrix(
            [self.q2j, self.q3j])

    def test_label(self):
        assert self.J._fill_coordinate_list(None, 1, 'u') == Matrix([self.uj])
        assert self.J._fill_coordinate_list([], 3, 'u') == Matrix(
            [self.u0j, self.u1j, self.u2j])
        assert self.J._fill_coordinate_list([self.u0], 3, 'u', 1) == Matrix(
            [self.u0, self.u2j, self.u3j])

    def test_single_numbering(self):
        assert self.J._fill_coordinate_list(
            None, 1, number_single=True) == Matrix([self.q0j])
        assert self.J._fill_coordinate_list([], 1, 'u', 2, True) == Matrix(
            [self.u2j])
        assert self.J._fill_coordinate_list([], 3, 'q') == Matrix(
            [self.q0j, self.q1j, self.q2j])

    def test_too_many_coordinates_supplied(self):
        with pytest.raises(ValueError):
            self.J._fill_coordinate_list([self.q0, self.q1], 1)
        with pytest.raises(ValueError):
            self.J._fill_coordinate_list([self.u0, self.u1, None], 2, 'u')

    def test_incorrect_coordinate_type(self):
        with pytest.raises(TypeError):
            self.J._fill_coordinate_list([self.q0, symbols('q1')], 2)
        with pytest.raises(TypeError):
            self.J._fill_coordinate_list([self.q0 + self.q1, self.q1], 2)


class TestPinJoint1:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.P = Body('P')
        self.C = Body('C')
        self.l, self.m = symbols('l m')
        self.q, self.u = dynamicsymbols('q_J, u_J')
        self.Pj = PinJoint('J', self.P, self.C)

    def test_attributes(self):
        assert self.Pj.name == 'J'
        assert self.Pj.parent == self.P
        assert self.Pj.child == self.C

    def test_coordinates_and_kdes(self):
        assert self.Pj.coordinates == Matrix([self.q])
        assert self.Pj.speeds == Matrix([self.u])
        assert self.Pj.kdes == Matrix([self.u - self.q.diff(t)])

    def test_joint_axis(self):
        assert self.Pj.joint_axis == self.P.frame.x

    def test_points(self):
        assert self.Pj.child_point.pos_from(self.C.masscenter) == Vector(0)
        assert self.Pj.parent_point.pos_from(self.P.masscenter) == Vector(0)
        assert self.Pj.parent_point.pos_from(self.Pj._child_point) == Vector(0)
        assert self.C.masscenter.pos_from(self.P.masscenter) == Vector(0)

    def test_interframes(self):
        assert self.Pj.parent_interframe == self.P.frame
        assert self.Pj.child_interframe == self.C.frame

    def test_str(self):
        assert self.Pj.__str__() == 'PinJoint: J  parent: P  child: C'


class TestPinJoint2:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.P1 = Body('P1')
        self.C1 = Body('C1')
        self.l, self.m = symbols('l m')
        self.Pint = ReferenceFrame('P_int')
        self.Pint.orient_axis(self.P1.frame, self.P1.y, pi / 2)
        self.J1 = PinJoint('J1', self.P1, self.C1,
                           parent_point=self.l*self.P1.frame.x,
                           child_point=self.m*self.C1.frame.y,
                           joint_axis=self.P1.frame.z,
                           parent_interframe=self.Pint)

    def test_joint_axis(self):
        assert self.J1._joint_axis == self.P1.frame.z

    def test_points(self):
        assert self.J1._child_point.pos_from(
            self.C1.masscenter) == self.m * self.C1.frame.y
        assert self.J1._parent_point.pos_from(
            self.P1.masscenter) == self.l * self.P1.frame.x
        assert self.J1._parent_point.pos_from(self.J1._child_point) == Vector(0)
        assert self.P1.masscenter.pos_from(
            self.C1.masscenter
        ) == -self.l*self.P1.frame.x + self.m*self.C1.frame.y

    def test_interframes(self):
        assert self.J1.parent_interframe == self.Pint
        assert self.J1.child_interframe == self.C1.frame


class TestPinJoint3:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.q, self.u = dynamicsymbols('q, u')
        self.N, self.A, self.P, self.C, self.Pint, self.Cint = _generate_body(
            True)
        self.parent_point = self.P.masscenter.locatenew('parent_point',
                                                        self.N.x + self.N.y)
        self.child_point = self.C.masscenter.locatenew('child_point',
                                                       self.C.y + self.C.z)
        self.J = PinJoint('J', self.P, self.C, self.q, self.u,
                          parent_point=self.parent_point,
                          child_point=self.child_point,
                          parent_interframe=self.Pint,
                          child_interframe=self.Cint, joint_axis=self.N.z)
    def test_joint_axis(self):
        assert self.J.joint_axis == self.N.z

    def test_points(self):
        assert self.J.parent_point.vel(self.N) == 0
        assert self.J.parent_point == self.parent_point
        assert self.J.child_point == self.child_point
        assert self.J.child_point.pos_from(
            self.P.masscenter) == self.N.x + self.N.y
        assert self.J.parent_point.pos_from(
            self.C.masscenter) == self.C.y + self.C.z
        assert self.C.masscenter.pos_from(
            self.P.masscenter) == self.N.x + self.N.y - self.C.y - self.C.z
        assert self.C.masscenter.vel(self.N).express(self.N) == (
            self.u * sin(self.q) - self.u * cos(self.q)) * self.N.x + (
            -self.u * sin(self.q) - self.u * cos(self.q)) * self.N.y

    def test_interframes(self):
        assert self.J.parent_interframe == self.Pint
        assert self.J.child_interframe == self.Cint


class TestPinJointDoublePendulum:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.q1, self.q2 = dynamicsymbols('q1 q2')
        self.u1, self.u2 = dynamicsymbols('u1 u2')
        self.m, self.l = symbols('m l')
        self.N = ReferenceFrame('N')
        self.A = ReferenceFrame('A')
        self.B = ReferenceFrame('B')
        self.C = Body('C', frame=self.N)  # ceiling
        self.PartP = Body('P', frame=self.A, mass=self.m)
        self.PartR = Body('R', frame=self.B, mass=self.m)
        self.J1 = PinJoint('J1', self.C, self.PartP, speeds=self.u1,
                           coordinates=self.q1, child_point=-self.l*self.A.x,
                           joint_axis=self.C.frame.z)
        self.J2 = PinJoint('J2', self.PartP, self.PartR, speeds=self.u2,
                           coordinates=self.q2, child_point=-self.l*self.B.x,
                           joint_axis=self.PartP.frame.z)

    def test_orientation(self):
        assert self.N.dcm(self.A) == Matrix(
            [[cos(self.q1), -sin(self.q1), 0],
             [sin(self.q1), cos(self.q1), 0], [0, 0, 1]])
        assert self.A.dcm(self.B) == Matrix(
            [[cos(self.q2), -sin(self.q2), 0],
             [sin(self.q2), cos(self.q2), 0], [0, 0, 1]])
        assert _simplify_matrix(self.N.dcm(self.B)) == Matrix(
            [[cos(self.q1 + self.q2), -sin(self.q1 + self.q2), 0],
             [sin(self.q1 + self.q2), cos(self.q1 + self.q2), 0],
             [0, 0, 1]])

    def test_angular_velocity(self):
        assert self.A.ang_vel_in(self.N) == self.u1 * self.N.z
        assert self.B.ang_vel_in(self.A) == self.u2 * self.A.z
        assert self.B.ang_vel_in(self.N) == (
            self.u1 * self.N.z + self.u2 * self.A.z)

    def test_kdes(self):
        assert self.J1.kdes == Matrix([self.u1 - self.q1.diff(t)])
        assert self.J2.kdes == Matrix([self.u2 - self.q2.diff(t)])

    def test_linear_velocity(self):
        assert self.PartP.masscenter.vel(self.N) == self.l*self.u1*self.A.y
        assert self.PartR.masscenter.vel(self.A) == self.l*self.u2*self.B.y
        assert self.PartR.masscenter.vel(self.N) == (
            self.l*self.u1*self.A.y + self.l*(self.u1 + self.u2)*self.B.y)


class TestPinJointChaosPendulum():

    @pytest.fixture(autouse=True)
    def setup(self):
        self.mA, self.mB = symbols('mA, mB')
        self.lA, self.lB, self.h = symbols('lA, lB, h')
        self.theta, self.phi = dynamicsymbols('theta phi')
        self.omega, self.alpha = dynamicsymbols('omega alpha')
        self.N = ReferenceFrame('N')
        self.A = ReferenceFrame('A')
        self.B = ReferenceFrame('B')
        self.lA = (self.lB - self.h / 2) / 2
        self.lC = (self.lB/2 + self.h/4)
        self.rod = Body('rod', frame=self.A, mass=self.mA)
        self.plate = Body('plate', mass=self.mB, frame=self.B)
        self.C = Body('C', frame=self.N)
        self.J1 = PinJoint('J1', self.C, self.rod, coordinates=self.theta,
                           speeds=self.omega, child_point=self.lA*self.A.z,
                           joint_axis=self.N.y)
        self.J2 = PinJoint('J2', self.rod, self.plate, coordinates=self.phi,
                           speeds=self.alpha, parent_point=self.lC*self.A.z,
                           joint_axis=self.A.z)

    def test_orientation(self):
        assert self.A.dcm(self.N) == Matrix([
            [cos(self.theta), 0, -sin(self.theta)],
            [0, 1, 0],
            [sin(self.theta), 0, cos(self.theta)]])
        assert self.A.dcm(self.B) == Matrix([
            [cos(self.phi), -sin(self.phi), 0],
            [sin(self.phi), cos(self.phi), 0],
            [0, 0, 1]])
        assert self.B.dcm(self.N) == Matrix([
            [cos(self.phi)*cos(self.theta), sin(self.phi),
             -sin(self.theta)*cos(self.phi)],
            [-sin(self.phi)*cos(self.theta), cos(self.phi),
             sin(self.phi)*sin(self.theta)],
            [sin(self.theta), 0, cos(self.theta)]])

    def test_angular_velocity(self):
        assert self.A.ang_vel_in(self.N) == self.omega*self.N.y
        assert self.A.ang_vel_in(self.B) == -self.alpha*self.A.z
        assert self.N.ang_vel_in(self.B) == (
            -self.omega*self.N.y - self.alpha*self.A.z)

    def test_kdes(self):
        assert self.J1.kdes == Matrix([self.omega - self.theta.diff(t)])
        assert self.J2.kdes == Matrix([self.alpha - self.phi.diff(t)])

    def test_masscenters_positions(self):
        assert self.C.masscenter.pos_from(
            self.rod.masscenter) == self.lA*self.A.z
        assert self.rod.masscenter.pos_from(
            self.plate.masscenter) == - self.lC * self.A.z

    def test_linear_velocities(self):
        assert self.rod.masscenter.vel(self.N) == (
            (self.h/4 - self.lB/2)*self.omega*self.A.x)
        assert self.plate.masscenter.vel(self.N) == (
            ((self.h/4 - self.lB/2)*self.omega +
             (self.h/4 + self.lB/2)*self.omega)*self.A.x)


class TestPinJointInterframe:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.q, self.u = dynamicsymbols('q, u')

    def test_not_connected(self):
        N, A, P, C = _generate_body()
        Pint, Cint = ReferenceFrame('Pint'), ReferenceFrame('Cint')
        with pytest.raises(ValueError):
            PinJoint('J', P, C, parent_interframe=Pint)
        with pytest.raises(ValueError):
            PinJoint('J', P, C, child_interframe=Cint)

    def test_not_fixed_interframe(self):
        N, A, P, C = _generate_body()
        Pint, Cint = ReferenceFrame('Pint'), ReferenceFrame('Cint')
        Pint.orient_axis(N, N.z, self.q)
        Cint.orient_axis(A, A.z, self.q)
        with pytest.raises(ValueError):
            PinJoint('J', P, C, parent_interframe=Pint)
        with pytest.raises(ValueError):
            PinJoint('J', P, C, child_interframe=Cint)

    def test_parent_interframe_only(self):
        N, A, P, C = _generate_body()
        Pint = ReferenceFrame('Pint')
        Pint.orient_body_fixed(N, (pi / 4, pi, pi / 3), 'xyz')
        PinJoint('J', P, C, self.q, self.u, parent_point=N.x, child_point=-C.y,
                 parent_interframe=Pint, joint_axis=Pint.x)
        assert _simplify_matrix(N.dcm(A)) - Matrix([
            [-1 / 2, sqrt(3) * cos(self.q) / 2, -sqrt(3) * sin(self.q) / 2],
            [sqrt(6) / 4, sqrt(2) * (2 * sin(self.q) + cos(self.q)) / 4,
             sqrt(2) * (-sin(self.q) + 2 * cos(self.q)) / 4],
            [sqrt(6) / 4, sqrt(2) * (-2 * sin(self.q) + cos(self.q)) / 4,
             -sqrt(2) * (sin(self.q) + 2 * cos(self.q)) / 4]]) == zeros(3)
        assert A.ang_vel_in(N) == self.u * Pint.x
        assert C.masscenter.pos_from(P.masscenter) == N.x + A.y
        assert C.masscenter.vel(N) == self.u * A.z
        assert P.masscenter.vel(Pint) == Vector(0)
        assert C.masscenter.vel(Pint) == self.u * A.z

    def test_child_interframe_only(self):
        N, A, P, C = _generate_body()
        Cint = ReferenceFrame('Cint')
        Cint.orient_body_fixed(A, (2 * pi / 3, -pi, pi / 2), 'xyz')
        PinJoint('J', P, C, self.q, self.u, parent_point=-N.z, child_point=C.x,
                 child_interframe=Cint, joint_axis=P.x + P.z)
        assert _simplify_matrix(N.dcm(A)) == Matrix([
            [-sqrt(2) * sin(self.q) / 2,
             -sqrt(3) * (cos(self.q) - 1) / 4 - cos(self.q) / 4 - S(1) / 4,
             sqrt(3) * (cos(self.q) + 1) / 4 - cos(self.q) / 4 + S(1) / 4],
            [cos(self.q), (sqrt(2) + sqrt(6)) * -sin(self.q) / 4,
             (-sqrt(2) + sqrt(6)) * sin(self.q) / 4],
            [sqrt(2) * sin(self.q) / 2,
             sqrt(3) * (cos(self.q) + 1) / 4 + cos(self.q) / 4 - S(1) / 4,
             sqrt(3) * (1 - cos(self.q)) / 4 + cos(self.q) / 4 + S(1) / 4]])
        assert A.ang_vel_in(N) == (
            sqrt(2) * self.u / 2 * N.x + sqrt(2) * self.u / 2 * N.z)
        assert C.masscenter.pos_from(P.masscenter) == - N.z - A.x
        assert C.masscenter.vel(N).simplify() == (
            -sqrt(6) - sqrt(2)) * self.u / 4 * A.y + (
                   -sqrt(2) + sqrt(6)) * self.u / 4 * A.z
        assert C.masscenter.vel(Cint) == Vector(0)

    def test_combination(self):
        N, A, P, C = _generate_body()
        Pint, Cint = ReferenceFrame('Pint'), ReferenceFrame('Cint')
        Pint.orient_body_fixed(N, (-pi / 2, pi, pi / 2), 'xyz')
        Cint.orient_body_fixed(A, (2 * pi / 3, -pi, pi / 2), 'xyz')
        PinJoint('J', P, C, self.q, self.u, parent_point=N.x - N.y,
                 child_point=-C.z, parent_interframe=Pint,
                 child_interframe=Cint, joint_axis=Pint.x + Pint.z)
        assert _simplify_matrix(N.dcm(A)) == Matrix([
            [cos(self.q), (sqrt(2) + sqrt(6)) * -sin(self.q) / 4,
             (-sqrt(2) + sqrt(6)) * sin(self.q) / 4],
            [-sqrt(2) * sin(self.q) / 2,
             -sqrt(3) * (cos(self.q) + 1) / 4 - cos(self.q) / 4 + S(1) / 4,
             sqrt(3) * (cos(self.q) - 1) / 4 - cos(self.q) / 4 - S(1) / 4],
            [sqrt(2) * sin(self.q) / 2,
             sqrt(3) * (cos(self.q) - 1) / 4 + cos(self.q) / 4 + S(1) / 4,
             -sqrt(3) * (cos(self.q) + 1) / 4 + cos(self.q) / 4 - S(1) / 4]])
        assert A.ang_vel_in(N) == sqrt(2) * self.u / 2 * Pint.x + sqrt(
            2) * self.u / 2 * Pint.z
        assert C.masscenter.pos_from(P.masscenter) == N.x - N.y + A.z
        N_v_C = (-sqrt(2) + sqrt(6)) * self.u / 4 * A.x
        assert C.masscenter.vel(N).simplify() == N_v_C
        assert C.masscenter.vel(Pint).simplify() == N_v_C
        assert C.masscenter.vel(Cint) == Vector(0)


class TestPinJointJointAxis():

    @pytest.fixture(autouse=True)
    def setup(self):
        self.q, self.u = dynamicsymbols('q, u')
        self.N, self.A, self.P, self.C, self.Pint, self.Cint = _generate_body(
            True)

    def test_parent_as_reference(self):
        pin = PinJoint('J', self.P, self.C, self.q, self.u,
                       parent_interframe=self.Pint, child_interframe=self.Cint,
                       joint_axis=self.P.y)
        assert pin.joint_axis == self.P.y
        assert self.N.dcm(self.A) == Matrix([[sin(self.q), 0, cos(self.q)],
                                             [0, -1, 0],
                                             [cos(self.q), 0, -sin(self.q)]])

    def test_parent_interframe_as_reference(self):
        pin = PinJoint('J', self.P, self.C, self.q, self.u,
                       parent_interframe=self.Pint, child_interframe=self.Cint,
                       joint_axis=self.Pint.y)
        assert pin.joint_axis == self.Pint.y
        assert self.N.dcm(self.A) == Matrix([[-sin(self.q), 0, cos(self.q)],
                                             [0, -1, 0],
                                             [cos(self.q), 0, sin(self.q)]])

    def test_joint_axis_with_interframes_as_vectors_1(self):
        N, A, P, C = _generate_body()
        pin = PinJoint('J', P, C, self.q, self.u, parent_interframe=N.z,
                       child_interframe=-C.z, joint_axis=N.z)
        assert pin.joint_axis == N.z
        assert N.dcm(A) == Matrix([[-cos(self.q), -sin(self.q), 0],
                                   [-sin(self.q), cos(self.q), 0],
                                   [0, 0, -1]])

    def test_joint_axis_with_interframes_as_vectors_2(self):
        N, A, P, C = _generate_body()
        pin = PinJoint('J', P, C, self.q, self.u, parent_interframe=N.z,
                       child_interframe=-C.z, joint_axis=N.x)
        assert pin.joint_axis == N.x
        assert N.dcm(A) == Matrix([[-1, 0, 0], [0, cos(self.q), sin(self.q)],
                                   [0, sin(self.q), -cos(self.q)]])

    def test_time_varying_axis(self):
        with pytest.raises(ValueError):
            PinJoint('J', self.P, self.C,
                     joint_axis=cos(self.q)*self.N.x + sin(self.q)*self.N.y)

    def test_joint_axis_provided_in_child_frame(self):
        with pytest.raises(ValueError):
            PinJoint('J', self.P, self.C, joint_axis=self.C.x)

    def test_invalid_combinations(self):
        with pytest.raises(ValueError):
            PinJoint('J', self.P, self.C, joint_axis=self.P.x + self.C.y)
        with pytest.raises(ValueError):
            PinJoint('J', self.P, self.C, parent_interframe=self.Pint,
                     child_interframe=self.Cint,
                     joint_axis=self.Pint.x + self.C.y)
        with pytest.raises(ValueError):
            PinJoint('J', self.P, self.C, parent_interframe=self.Pint,
                     child_interframe=self.Cint,
                     joint_axis=self.P.x + self.Cint.y)

    def test_valid_special_combination(self):
        PinJoint('J', self.P, self.C, parent_interframe=self.Pint,
                 child_interframe=self.Cint, joint_axis=self.Pint.x + self.P.y)

    def test_zero_vector_invalid(self):
        with pytest.raises(Exception):
            PinJoint('J', self.P, self.C, parent_interframe=self.Pint,
                     child_interframe=self.Cint, joint_axis=Vector(0))
        with pytest.raises(Exception):
            PinJoint('J', self.P, self.C, parent_interframe=self.Pint,
                     child_interframe=self.Cint,
                     joint_axis=self.P.y + self.Pint.y)


class TestPinJointArbitraryAxis:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.q, self.u = dynamicsymbols('q_J, u_J')
        self.N, self.A, self.P, self.C = _generate_body()

    def test_attach_through_mass_centers_with_opposite_axes(self):
        PinJoint('J', self.P, self.C, child_interframe=-self.A.x)
        assert (-self.A.x).angle_between(self.N.x) == 0
        assert -self.A.x.express(self.N) == self.N.x
        assert self.A.dcm(self.N) == Matrix([[-1, 0, 0],
                                             [0, -cos(self.q), -sin(self.q)],
                                             [0, -sin(self.q), cos(self.q)]])
        assert self.A.ang_vel_in(self.N) == self.u*self.N.x
        assert self.A.ang_vel_in(self.N).magnitude() == sqrt(self.u**2)
        assert self.C.masscenter.pos_from(self.P.masscenter) == 0
        assert self.C.masscenter.pos_from(self.P.masscenter).express(
            self.N).simplify() == 0
        assert self.C.masscenter.vel(self.N) == 0

    def test_child_joint_unit_vector_from_masscenter(self):
        """Axes are different and parent joint is at masscenter but child joint
        is at a unit vector from child masscenter."""
        PinJoint('J', self.P, self.C, child_interframe=self.A.y,
                 child_point=self.A.x)
        assert self.A.y.angle_between(self.N.x) == 0  # Axis are aligned
        assert self.A.y.express(self.N) == self.N.x
        assert self.A.dcm(self.N) == Matrix([[0, -cos(self.q), -sin(self.q)],
                                            [1, 0, 0],
                                            [0, -sin(self.q), cos(self.q)]])
        assert self.A.ang_vel_in(self.N) == self.u * self.N.x
        assert self.A.ang_vel_in(self.N).express(self.A) == self.u * self.A.y
        assert self.A.ang_vel_in(self.N).magnitude() == sqrt(self.u**2)
        assert self.A.ang_vel_in(self.N).cross(self.A.y) == 0
        assert self.C.masscenter.vel(self.N) == self.u * self.A.z
        assert self.C.masscenter.pos_from(self.P.masscenter) == -self.A.x
        assert self.C.masscenter.pos_from(self.P.masscenter).express(
            self.N).simplify() == cos(self.q)*self.N.y + sin(self.q)*self.N.z
        assert self.C.masscenter.vel(self.N).angle_between(self.A.x) == pi/2

    def test_parent_joint_unit_vector_from_masscenter(self):
        """Similar to previous case but wrt parent body."""
        PinJoint('J', self.P, self.C, parent_interframe=self.N.y,
                 parent_point=self.N.x)
        assert self.N.y.angle_between(self.A.x) == 0  # Axis are aligned
        assert self.N.y.express(self.A) == self.A.x
        assert self.A.dcm(self.N) == Matrix([[0, 1, 0],
                                             [-cos(self.q), 0, sin(self.q)],
                                             [sin(self.q), 0, cos(self.q)]])
        assert self.A.ang_vel_in(self.N) == self.u*self.N.y
        assert self.A.ang_vel_in(self.N).express(self.A) == self.u*self.A.x
        assert self.A.ang_vel_in(self.N).magnitude() == sqrt(self.u**2)
        angle = self.A.ang_vel_in(self.N).angle_between(self.A.x)
        assert angle.xreplace({self.u: 1}) == 0
        assert self.C.masscenter.vel(self.N) == 0
        assert self.C.masscenter.pos_from(self.P.masscenter) == self.N.x

    def test_joint_pos_with_different_axes_1(self):
        PinJoint('J', self.P, self.C, parent_point=self.N.x,
                 child_point=self.A.x, child_interframe=self.A.x + self.A.y)
        assert expand_mul(self.N.x.angle_between(self.A.x + self.A.y)) == 0
        assert (self.A.x + self.A.y).express(
            self.N).simplify() == sqrt(2)*self.N.x
        assert _simplify_matrix(self.A.dcm(self.N)) == Matrix([
            [sqrt(2)/2, -sqrt(2)*cos(self.q)/2, -sqrt(2)*sin(self.q)/2],
            [sqrt(2)/2, sqrt(2)*cos(self.q)/2, sqrt(2)*sin(self.q)/2],
            [0, -sin(self.q), cos(self.q)]])
        assert self.A.ang_vel_in(self.N) == self.u*self.N.x
        assert (self.A.ang_vel_in(self.N).express(self.A).simplify() ==
                (self.u*self.A.x + self.u*self.A.y)/sqrt(2))
        assert self.A.ang_vel_in(self.N).magnitude() == sqrt(self.u**2)
        angle = self.A.ang_vel_in(self.N).angle_between(self.A.x + self.A.y)
        assert angle.xreplace({self.u: 1}) == 0
        assert self.C.masscenter.vel(
            self.N).simplify() == (self.u * self.A.z)/sqrt(2)
        assert self.C.masscenter.pos_from(
            self.P.masscenter) == self.N.x - self.A.x
        assert self.C.masscenter.pos_from(self.P.masscenter).express(
            self.N).simplify() == ((1 - sqrt(2)/2)*self.N.x +
                                    sqrt(2)*cos(self.q)/2*self.N.y +
                                    sqrt(2)*sin(self.q)/2*self.N.z)
        assert self.C.masscenter.vel(self.N).express(
            self.N).simplify() == (-sqrt(2)*self.u*sin(self.q)/2*self.N.y +
                                   sqrt(2)*self.u*cos(self.q)/2*self.N.z)
        assert self.C.masscenter.vel(self.N).angle_between(self.A.x) == pi/2

    def test_joint_pos_with_different_axes_2(self):
        PinJoint('J', self.P, self.C, parent_point=self.N.x,
                 child_point=self.A.x,
                 child_interframe=self.A.x + self.A.y - self.A.z)
        assert expand_mul(self.N.x.angle_between(
            self.A.x + self.A.y - self.A.z)) == 0
        assert (self.A.x + self.A.y - self.A.z).express(
            self.N).simplify() == sqrt(3)*self.N.x
        assert _simplify_matrix(self.A.dcm(self.N)) == Matrix([
            [sqrt(3)/3, -sqrt(6)*sin(self.q + pi/4)/3,
             sqrt(6)*cos(self.q + pi/4)/3],
            [sqrt(3)/3, sqrt(6)*cos(self.q + pi/12)/3,
             sqrt(6)*sin(self.q + pi/12)/3],
            [-sqrt(3)/3, sqrt(6)*cos(self.q + 5*pi/12)/3,
             sqrt(6)*sin(self.q + 5*pi/12)/3]])
        assert self.A.ang_vel_in(self.N) == self.u*self.N.x
        assert self.A.ang_vel_in(self.N).express(self.A).simplify() == (
            self.u*self.A.x + self.u*self.A.y - self.u*self.A.z)/sqrt(3)
        assert self.A.ang_vel_in(self.N).magnitude() == sqrt(self.u**2)
        angle = self.A.ang_vel_in(self.N).angle_between(
            self.A.x + self.A.y - self.A.z)
        assert angle.xreplace({self.u: 1}) == 0
        assert self.C.masscenter.vel(self.N).simplify() == (
            self.u*self.A.y + self.u*self.A.z) / sqrt(3)
        assert self.C.masscenter.pos_from(self.P.masscenter) == (
            self.N.x - self.A.x)
        assert self.C.masscenter.pos_from(self.P.masscenter).express(
            self.N).simplify() == ((1 - sqrt(3)/3)*self.N.x +
                                    sqrt(6)*sin(self.q + pi/4)/3*self.N.y -
                                    sqrt(6)*cos(self.q + pi/4)/3*self.N.z)
        assert self.C.masscenter.vel(self.N).express(self.N).simplify() == (
                sqrt(6)*self.u*cos(self.q + pi/4)/3*self.N.y +
                sqrt(6)*self.u*sin(self.q + pi/4)/3*self.N.z)
        assert self.C.masscenter.vel(self.N).angle_between(self.A.x) == pi/2

    def test_joint_pos_with_different_axes_3(self):
        m, n = symbols('m n')
        PinJoint('J', self.P, self.C, parent_point=m * self.N.x,
                 child_point=n * self.A.x,
                 child_interframe=self.A.x + self.A.y - self.A.z,
                 parent_interframe=self.N.x - self.N.y + self.N.z)
        angle = (self.N.x - self.N.y + self.N.z).angle_between(
            self.A.x + self.A.y - self.A.z)
        assert expand_mul(angle) == 0
        assert ((self.A.x - self.A.y + self.A.z).express(self.N).simplify() ==
                (-4*cos(self.q)/3 - S(1)/3)*self.N.x + (S(1)/3 -
                 4*sin(self.q + pi/6)/3)*self.N.y +
                (4*cos(self.q + pi/3)/3 - S(1)/3)*self.N.z)
        assert _simplify_matrix(self.A.dcm(self.N)) == Matrix([
            [S(1)/3 - 2*cos(self.q)/3, -2*sin(self.q + pi/6)/3 - S(1)/3,
             2*cos(self.q + pi/3)/3 + S(1)/3],
            [2*cos(self.q + pi/3)/3 + S(1)/3, 2*cos(self.q)/3 - S(1)/3,
             2*sin(self.q + pi/6)/3 + S(1)/3],
            [-2*sin(self.q + pi/6)/3 - S(1)/3, 2*cos(self.q + pi/3)/3 + S(1)/3,
             2*cos(self.q)/3 - S(1)/3]])
        assert self.A.ang_vel_in(self.N) == (
            self.u*self.N.x - self.u*self.N.y + self.u*self.N.z)/sqrt(3)
        assert self.A.ang_vel_in(self.N).express(self.A).simplify() == (
            self.u*self.A.x + self.u*self.A.y - self.u*self.A.z) / sqrt(3)
        assert self.A.ang_vel_in(self.N).magnitude() == sqrt(self.u**2)
        angle = self.A.ang_vel_in(self.N).angle_between(
            self.A.x + self.A.y - self.A.z)
        assert angle.xreplace({self.u: 1}) == 0
        assert self.C.masscenter.vel(self.N).simplify() == (
            sqrt(3)*n*self.u/3*self.A.y + sqrt(3)*n*self.u/3*self.A.z)
        assert self.C.masscenter.pos_from(self.P.masscenter) == (
            m*self.N.x - n*self.A.x)
        assert (self.C.masscenter.pos_from(self.P.masscenter).express(
            self.N).simplify() == ((m + n*(2*cos(self.q) - 1)/3)*self.N.x +
                                    n*(2*sin(self.q + pi/6) +
                                    1)/3*self.N.y - n*(2*cos(self.q + pi/3) +
                                    1)/3*self.N.z))
        assert self.C.masscenter.vel(self.N).express(self.N).simplify() == (
            - 2*n*self.u*sin(self.q)/3*self.N.x +
            2*n*self.u*cos(self.q + pi/6)/3*self.N.y +
            2*n*self.u*sin(self.q + pi/3)/3*self.N.z)
        assert self.C.masscenter.vel(self.N).dot(
            self.N.x - self.N.y + self.N.z).simplify() == 0


class TestCreateAlignedFramePi():

    @pytest.fixture(autouse=True)
    def setup(self):
        self.N, self.A, self.P, self.C = _generate_body()

    def test_case_1(self):
        f = Joint._create_aligned_interframe(self.P, -self.P.x, self.P.x)
        assert f.z == self.P.z

    def test_case_2(self):
        f = Joint._create_aligned_interframe(self.P, -self.P.y, self.P.y)
        assert f.x == self.P.x

    def test_case_3(self):
        f = Joint._create_aligned_interframe(self.P, -self.P.z, self.P.z)
        assert f.y == self.P.y

    def test_case_4(self):
        f = Joint._create_aligned_interframe(self.P,
                                             -self.P.x - self.P.y,
                                             self.P.x + self.P.y)
        assert f.z == self.P.z

    def test_case_5(self):
        f = Joint._create_aligned_interframe(self.P,
                                             -self.P.y - self.P.z,
                                             self.P.y + self.P.z)
        assert f.x == self.P.x

    def test_case_6(self):
        f = Joint._create_aligned_interframe(self.P,
                                             -self.P.x - self.P.z,
                                             self.P.x + self.P.z)
        assert f.y == self.P.y

    def test_case_7(self):
        f = Joint._create_aligned_interframe(self.P,
                                             -self.P.x - self.P.y - self.P.z,
                                             self.P.x + self.P.y + self.P.z)
        assert f.y - f.z == self.P.y - self.P.z


class TestPinJointAxis:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.q, self.u = dynamicsymbols('q u')
        self.N, self.A, self.P, self.C, self.Pint, self.Cint = _generate_body(
            True)
        self.N_R_A = Matrix([[0, sin(self.q), cos(self.q)],
                             [0, -cos(self.q), sin(self.q)],
                             [1, 0, 0]])

    def test_default_joint_axis(self):
        J = PinJoint('J', self.P, self.C, self.q, self.u,
                     parent_interframe=self.Pint, child_interframe=self.Cint)
        assert J.joint_axis == self.Pint.x

    def test_same_joint_axis_expressed_in_different_frames_1(self):
        PinJoint('J', self.P, self.C, self.q, self.u,
                 parent_interframe=self.Pint, child_interframe=self.Cint,
                 joint_axis=self.N.z)
        assert self.N.dcm(self.A) == self.N_R_A

    def test_same_joint_axis_expressed_in_different_frames_2(self):
        PinJoint('J', self.P, self.C, self.q, self.u,
                 parent_interframe=self.Pint, child_interframe=self.Cint,
                 joint_axis=-self.Pint.z)
        assert self.N.dcm(self.A) == self.N_R_A

    def test_time_varying_joint_axis(self):
        with pytest.raises(ValueError):
            PinJoint('J', self.P, self.C, joint_axis=self.q * self.N.z)


class TestLocateJointPos:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.N, self.A, self.P, self.C = _generate_body()
        self.q = dynamicsymbols('q')

    def test_vector_and_default(self):
        joint = PinJoint('J', self.P, self.C, parent_point=self.N.y + self.N.z)
        assert joint.parent_point.name == 'J_P_joint'
        assert joint.parent_point.pos_from(
            self.P.masscenter) == self.N.y + self.N.z
        assert joint.child_point == self.C.masscenter

    def test_point_objects(self):
        parent_point = self.P.masscenter.locatenew('p', self.N.y + self.N.z)
        joint = PinJoint('J', self.P, self.C, parent_point=parent_point,
                         child_point=self.C.masscenter)
        assert joint.parent_point == parent_point
        assert joint.child_point == self.C.masscenter

    def test_invalid_type(self):
        with pytest.raises(TypeError):
            PinJoint('J', self.P, self.C,
                     parent_point=self.N.x.to_matrix(self.N))

    def test_time_varying_positions_1(self):
        with pytest.raises(ValueError):
            PinJoint('J', self.P, self.C, parent_point=self.q * self.N.x)

    def test_time_varying_positions_2(self):
        child_point = self.C.masscenter.locatenew('p', self.q * self.A.y)
        with pytest.raises(ValueError):
            PinJoint('J', self.P, self.C, child_point=child_point)

    def test_undefined_position(self):
        child_point = Point('p')
        with pytest.raises(ValueError):
            PinJoint('J', self.P, self.C, child_point=child_point)


class TestLocateJointFrame:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.q = dynamicsymbols('q')
        self.N, self.A, self.P, self.C = _generate_body()
        self.parent_interframe = ReferenceFrame('parent_int_frame')
        self.child_interframe = ReferenceFrame('child_int_frame')

    def test_rotated_frame_and_default(self):
        self.parent_interframe.orient_axis(self.N, self.N.z, 1)
        joint = PinJoint('J', self.P, self.C,
                         parent_interframe=self.parent_interframe)
        assert joint.parent_interframe == self.parent_interframe
        assert joint.parent_interframe.ang_vel_in(self.N) == 0
        assert joint.child_interframe == self.A

    def test_time_varying_orientations(self):
        self.parent_interframe.orient_axis(self.N, self.N.z, self.q)
        with pytest.raises(ValueError):
            PinJoint('J', self.P, self.C,
                     parent_interframe=self.parent_interframe)

    def test_undefined_frame(self):
        """Child interframe is defined with respect to parent."""
        self.child_interframe.orient_axis(self.N, self.N.z, 1)
        with pytest.raises(ValueError):
            PinJoint('J', self.P, self.C,
                     child_interframe=self.child_interframe)


def test_sliding_joint():
    _, _, P, C = _generate_body()
    q, u = dynamicsymbols('q_S, u_S')
    S = PrismaticJoint('S', P, C)
    assert S.name == 'S'
    assert S.parent == P
    assert S.child == C
    assert S.coordinates == Matrix([q])
    assert S.speeds == Matrix([u])
    assert S.kdes == Matrix([u - q.diff(t)])
    assert S.joint_axis == P.frame.x
    assert S.child_point.pos_from(C.masscenter) == Vector(0)
    assert S.parent_point.pos_from(P.masscenter) == Vector(0)
    assert S.parent_point.pos_from(S.child_point) == - q * P.frame.x
    assert P.masscenter.pos_from(C.masscenter) == - q * P.frame.x
    assert C.masscenter.vel(P.frame) == u * P.frame.x
    assert P.ang_vel_in(C) == 0
    assert C.ang_vel_in(P) == 0
    assert S.__str__() == 'PrismaticJoint: S  parent: P  child: C'

    N, A, P, C = _generate_body()
    l, m = symbols('l m')
    Pint = ReferenceFrame('P_int')
    Pint.orient_axis(P.frame, P.y, pi / 2)
    S = PrismaticJoint('S', P, C, parent_point=l * P.frame.x,
                       child_point=m * C.frame.y, joint_axis=P.frame.z,
                       parent_interframe=Pint)

    assert S.joint_axis == P.frame.z
    assert S.child_point.pos_from(C.masscenter) == m * C.frame.y
    assert S.parent_point.pos_from(P.masscenter) == l * P.frame.x
    assert S.parent_point.pos_from(S.child_point) == - q * P.frame.z
    assert P.masscenter.pos_from(C.masscenter) == - l*N.x - q*N.z + m*A.y
    assert C.masscenter.vel(P.frame) == u * P.frame.z
    assert P.masscenter.vel(Pint) == Vector(0)
    assert C.ang_vel_in(P) == 0
    assert P.ang_vel_in(C) == 0

    _, _, P, C = _generate_body()
    Pint = ReferenceFrame('P_int')
    Pint.orient_axis(P.frame, P.y, pi / 2)
    S = PrismaticJoint('S', P, C, parent_point=l * P.frame.z,
                       child_point=m * C.frame.x, joint_axis=P.frame.z,
                       parent_interframe=Pint)
    assert S.joint_axis == P.frame.z
    assert S.child_point.pos_from(C.masscenter) == m * C.frame.x
    assert S.parent_point.pos_from(P.masscenter) == l * P.frame.z
    assert S.parent_point.pos_from(S.child_point) == - q * P.frame.z
    assert P.masscenter.pos_from(C.masscenter) == (-l - q)*P.frame.z + m*C.frame.x
    assert C.masscenter.vel(P.frame) == u * P.frame.z
    assert C.ang_vel_in(P) == 0
    assert P.ang_vel_in(C) == 0


def test_sliding_joint_arbitrary_axis():
    q, u = dynamicsymbols('q_S, u_S')

    N, A, P, C = _generate_body()
    PrismaticJoint('S', P, C, child_interframe=-A.x)

    assert (-A.x).angle_between(N.x) == 0
    assert -A.x.express(N) == N.x
    assert A.dcm(N) == Matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
    assert C.masscenter.pos_from(P.masscenter) == q * N.x
    assert C.masscenter.pos_from(P.masscenter).express(A).simplify() == -q * A.x
    assert C.masscenter.vel(N) == u * N.x
    assert C.masscenter.vel(N).express(A) == -u * A.x
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0

    #When axes are different and parent joint is at masscenter but child joint is at a unit vector from
    #child masscenter.
    N, A, P, C = _generate_body()
    PrismaticJoint('S', P, C, child_interframe=A.y, child_point=A.x)

    assert A.y.angle_between(N.x) == 0 #Axis are aligned
    assert A.y.express(N) == N.x
    assert A.dcm(N) == Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    assert C.masscenter.vel(N) == u * N.x
    assert C.masscenter.vel(N).express(A) == u * A.y
    assert C.masscenter.pos_from(P.masscenter) == q*N.x - A.x
    assert C.masscenter.pos_from(P.masscenter).express(N).simplify() == q*N.x + N.y
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0

    #Similar to previous case but wrt parent body
    N, A, P, C = _generate_body()
    PrismaticJoint('S', P, C, parent_interframe=N.y, parent_point=N.x)

    assert N.y.angle_between(A.x) == 0 #Axis are aligned
    assert N.y.express(A) ==  A.x
    assert A.dcm(N) == Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    assert C.masscenter.vel(N) == u * N.y
    assert C.masscenter.vel(N).express(A) == u * A.x
    assert C.masscenter.pos_from(P.masscenter) == N.x + q*N.y
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0

    #Both joint pos is defined but different axes
    N, A, P, C = _generate_body()
    PrismaticJoint('S', P, C, parent_point=N.x, child_point=A.x,
                   child_interframe=A.x + A.y)
    assert N.x.angle_between(A.x + A.y) == 0 #Axis are aligned
    assert (A.x + A.y).express(N) == sqrt(2)*N.x
    assert A.dcm(N) == Matrix([[sqrt(2)/2, -sqrt(2)/2, 0], [sqrt(2)/2, sqrt(2)/2, 0], [0, 0, 1]])
    assert C.masscenter.pos_from(P.masscenter) == (q + 1)*N.x - A.x
    assert C.masscenter.pos_from(P.masscenter).express(N) == (q - sqrt(2)/2 + 1)*N.x + sqrt(2)/2*N.y
    assert C.masscenter.vel(N).express(A) == u * (A.x + A.y)/sqrt(2)
    assert C.masscenter.vel(N) == u*N.x
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0

    N, A, P, C = _generate_body()
    PrismaticJoint('S', P, C, parent_point=N.x, child_point=A.x,
                   child_interframe=A.x + A.y - A.z)
    assert N.x.angle_between(A.x + A.y - A.z) == 0 #Axis are aligned
    assert (A.x + A.y - A.z).express(N) == sqrt(3)*N.x
    assert _simplify_matrix(A.dcm(N)) == Matrix([[sqrt(3)/3, -sqrt(3)/3, sqrt(3)/3],
                                                 [sqrt(3)/3, sqrt(3)/6 + S(1)/2, S(1)/2 - sqrt(3)/6],
                                                 [-sqrt(3)/3, S(1)/2 - sqrt(3)/6, sqrt(3)/6 + S(1)/2]])
    assert C.masscenter.pos_from(P.masscenter) == (q + 1)*N.x - A.x
    assert C.masscenter.pos_from(P.masscenter).express(N) == \
        (q - sqrt(3)/3 + 1)*N.x + sqrt(3)/3*N.y - sqrt(3)/3*N.z
    assert C.masscenter.vel(N) == u*N.x
    assert C.masscenter.vel(N).express(A) == sqrt(3)*u/3*A.x + sqrt(3)*u/3*A.y - sqrt(3)*u/3*A.z
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0

    N, A, P, C = _generate_body()
    m, n = symbols('m n')
    PrismaticJoint('S', P, C, parent_point=m*N.x, child_point=n*A.x,
                   child_interframe=A.x + A.y - A.z,
                   parent_interframe=N.x - N.y + N.z)
    # 0 angle means that the axis are aligned
    assert (N.x-N.y+N.z).angle_between(A.x+A.y-A.z).simplify() == 0
    assert (A.x+A.y-A.z).express(N) == N.x - N.y + N.z
    assert _simplify_matrix(A.dcm(N)) == Matrix([[-S(1)/3, -S(2)/3, S(2)/3],
                                                 [S(2)/3, S(1)/3, S(2)/3],
                                                 [-S(2)/3, S(2)/3, S(1)/3]])
    assert C.masscenter.pos_from(P.masscenter) == \
        (m + sqrt(3)*q/3)*N.x - sqrt(3)*q/3*N.y + sqrt(3)*q/3*N.z - n*A.x
    assert C.masscenter.pos_from(P.masscenter).express(N) == \
        (m + n/3 + sqrt(3)*q/3)*N.x + (2*n/3 - sqrt(3)*q/3)*N.y + (-2*n/3 + sqrt(3)*q/3)*N.z
    assert C.masscenter.vel(N) == sqrt(3)*u/3*N.x - sqrt(3)*u/3*N.y + sqrt(3)*u/3*N.z
    assert C.masscenter.vel(N).express(A) == sqrt(3)*u/3*A.x + sqrt(3)*u/3*A.y - sqrt(3)*u/3*A.z
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0


def test_cylindrical_joint():
    N, A, P, C = _generate_body()
    q0_def, q1_def, u0_def, u1_def = dynamicsymbols('q0:2_J, u0:2_J')
    Cj = CylindricalJoint('J', P, C)
    assert Cj.name == 'J'
    assert Cj.parent == P
    assert Cj.child == C
    assert Cj.coordinates == Matrix([q0_def, q1_def])
    assert Cj.speeds == Matrix([u0_def, u1_def])
    assert Cj.rotation_coordinate == q0_def
    assert Cj.translation_coordinate == q1_def
    assert Cj.rotation_speed == u0_def
    assert Cj.translation_speed == u1_def
    assert Cj.kdes == Matrix([u0_def - q0_def.diff(t), u1_def - q1_def.diff(t)])
    assert Cj.joint_axis == N.x
    assert Cj.child_point.pos_from(C.masscenter) == Vector(0)
    assert Cj.parent_point.pos_from(P.masscenter) == Vector(0)
    assert Cj.parent_point.pos_from(Cj._child_point) == -q1_def * N.x
    assert C.masscenter.pos_from(P.masscenter) == q1_def * N.x
    assert Cj.child_point.vel(N) == u1_def * N.x
    assert A.ang_vel_in(N) == u0_def * N.x
    assert Cj.parent_interframe == N
    assert Cj.child_interframe == A
    assert Cj.__str__() == 'CylindricalJoint: J  parent: P  child: C'

    q0, q1, u0, u1 = dynamicsymbols('q0:2, u0:2')
    l, m = symbols('l, m')
    N, A, P, C, Pint, Cint = _generate_body(True)
    Cj = CylindricalJoint('J', P, C, rotation_coordinate=q0, rotation_speed=u0,
                          translation_speed=u1, parent_point=m * N.x,
                          child_point=l * A.y, parent_interframe=Pint,
                          child_interframe=Cint, joint_axis=2 * N.z)
    assert Cj.coordinates == Matrix([q0, q1_def])
    assert Cj.speeds == Matrix([u0, u1])
    assert Cj.rotation_coordinate == q0
    assert Cj.translation_coordinate == q1_def
    assert Cj.rotation_speed == u0
    assert Cj.translation_speed == u1
    assert Cj.kdes == Matrix([u0 - q0.diff(t), u1 - q1_def.diff(t)])
    assert Cj.joint_axis == 2 * N.z
    assert Cj.child_point.pos_from(C.masscenter) == l * A.y
    assert Cj.parent_point.pos_from(P.masscenter) == m * N.x
    assert Cj.parent_point.pos_from(Cj._child_point) == -q1_def * N.z
    assert C.masscenter.pos_from(
        P.masscenter) == m * N.x + q1_def * N.z - l * A.y
    assert C.masscenter.vel(N) == u1 * N.z - u0 * l * A.z
    assert A.ang_vel_in(N) == u0 * N.z


def test_planar_joint():
    N, A, P, C = _generate_body()
    q0_def, q1_def, q2_def = dynamicsymbols('q0:3_J')
    u0_def, u1_def, u2_def = dynamicsymbols('u0:3_J')
    Cj = PlanarJoint('J', P, C)
    assert Cj.name == 'J'
    assert Cj.parent == P
    assert Cj.child == C
    assert Cj.coordinates == Matrix([q0_def, q1_def, q2_def])
    assert Cj.speeds == Matrix([u0_def, u1_def, u2_def])
    assert Cj.rotation_coordinate == q0_def
    assert Cj.planar_coordinates == Matrix([q1_def, q2_def])
    assert Cj.rotation_speed == u0_def
    assert Cj.planar_speeds == Matrix([u1_def, u2_def])
    assert Cj.kdes == Matrix([u0_def - q0_def.diff(t), u1_def - q1_def.diff(t),
                              u2_def - q2_def.diff(t)])
    assert Cj.rotation_axis == N.x
    assert Cj.planar_vectors == [N.y, N.z]
    assert Cj.child_point.pos_from(C.masscenter) == Vector(0)
    assert Cj.parent_point.pos_from(P.masscenter) == Vector(0)
    r_P_C = q1_def * N.y + q2_def * N.z
    assert Cj.parent_point.pos_from(Cj.child_point) == -r_P_C
    assert C.masscenter.pos_from(P.masscenter) == r_P_C
    assert Cj.child_point.vel(N) == u1_def * N.y + u2_def * N.z
    assert A.ang_vel_in(N) == u0_def * N.x
    assert Cj.parent_interframe == N
    assert Cj.child_interframe == A
    assert Cj.__str__() == 'PlanarJoint: J  parent: P  child: C'

    q0, q1, q2, u0, u1, u2 = dynamicsymbols('q0:3, u0:3')
    l, m = symbols('l, m')
    N, A, P, C, Pint, Cint = _generate_body(True)
    Cj = PlanarJoint('J', P, C, rotation_coordinate=q0, planar_coordinates=q1,
                     planar_speeds=[u1, u2], parent_point=m * N.x,
                     child_point=l * A.y, parent_interframe=Pint,
                     child_interframe=Cint)
    assert Cj.coordinates == Matrix([q0, q1, q2_def])
    assert Cj.speeds == Matrix([u0_def, u1, u2])
    assert Cj.rotation_coordinate == q0
    assert Cj.planar_coordinates == Matrix([q1, q2_def])
    assert Cj.rotation_speed == u0_def
    assert Cj.planar_speeds == Matrix([u1, u2])
    assert Cj.kdes == Matrix([u0_def - q0.diff(t), u1 - q1.diff(t),
                              u2 - q2_def.diff(t)])
    assert Cj.rotation_axis == Pint.x
    assert Cj.planar_vectors == [Pint.y, Pint.z]
    assert Cj.child_point.pos_from(C.masscenter) == l * A.y
    assert Cj.parent_point.pos_from(P.masscenter) == m * N.x
    assert Cj.parent_point.pos_from(Cj.child_point) == q1 * N.y + q2_def * N.z
    assert C.masscenter.pos_from(
        P.masscenter) == m * N.x - q1 * N.y - q2_def * N.z - l * A.y
    assert C.masscenter.vel(N) == -u1 * N.y - u2 * N.z + u0_def * l * A.x
    assert A.ang_vel_in(N) == u0_def * N.x


def test_planar_joint_advanced():
    # Tests whether someone is able to just specify two normals, which will form
    # the rotation axis seen from the parent and child body.
    # This specific example is a block on a slope, which has that same slope of
    # 30 degrees, so in the zero configuration the frames of the parent and
    # child are actually aligned.
    q0, q1, q2, u0, u1, u2 = dynamicsymbols('q0:3, u0:3')
    l1, l2 = symbols('l1:3')
    N, A, P, C = _generate_body()
    J = PlanarJoint('J', P, C, q0, [q1, q2], u0, [u1, u2],
                    parent_point=l1 * N.z,
                    child_point=-l2 * C.z,
                    parent_interframe=N.z + N.y / sqrt(3),
                    child_interframe=A.z + A.y / sqrt(3))
    assert J.rotation_axis.express(N) == (N.z + N.y / sqrt(3)).normalize()
    assert J.rotation_axis.express(A) == (A.z + A.y / sqrt(3)).normalize()
    assert J.rotation_axis.angle_between(N.z) == pi / 6
    assert N.dcm(A).xreplace({q0: 0, q1: 0, q2: 0}) == eye(3)
    N_R_A = Matrix([
        [cos(q0), -sqrt(3) * sin(q0) / 2, sin(q0) / 2],
        [sqrt(3) * sin(q0) / 2, 3 * cos(q0) / 4 + 1 / 4,
         sqrt(3) * (1 - cos(q0)) / 4],
        [-sin(q0) / 2, sqrt(3) * (1 - cos(q0)) / 4, cos(q0) / 4 + 3 / 4]])
    # N.dcm(A) == N_R_A did not work
    assert _simplify_matrix(N.dcm(A) - N_R_A) == zeros(3)


def test_spherical_joint():
    N, A, P, C = _generate_body()
    q0, q1, q2, u0, u1, u2 = dynamicsymbols('q0:3_S, u0:3_S')
    S = SphericalJoint('S', P, C)
    assert S.name == 'S'
    assert S.parent == P
    assert S.child == C
    assert S.coordinates == Matrix([q0, q1, q2])
    assert S.speeds == Matrix([u0, u1, u2])
    assert S.kdes == Matrix([u0 - q0.diff(t), u1 - q1.diff(t), u2 - q2.diff(t)])
    assert S.child_point.pos_from(C.masscenter) == Vector(0)
    assert S.parent_point.pos_from(P.masscenter) == Vector(0)
    assert S.parent_point.pos_from(S.child_point) == Vector(0)
    assert P.masscenter.pos_from(C.masscenter) == Vector(0)
    assert C.masscenter.vel(N) == Vector(0)
    assert P.ang_vel_in(C) == (-u0 * cos(q1) * cos(q2) - u1 * sin(q2)) * A.x + (
            u0 * sin(q2) * cos(q1) - u1 * cos(q2)) * A.y + (
                   -u0 * sin(q1) - u2) * A.z
    assert C.ang_vel_in(P) == (u0 * cos(q1) * cos(q2) + u1 * sin(q2)) * A.x + (
            -u0 * sin(q2) * cos(q1) + u1 * cos(q2)) * A.y + (
                   u0 * sin(q1) + u2) * A.z
    assert S.__str__() == 'SphericalJoint: S  parent: P  child: C'
    assert S._rot_type == 'BODY'
    assert S._rot_order == 123
    assert S._amounts is None


def test_spherical_joint_speeds_as_derivative_terms():
    # This tests checks whether the system remains valid if the user chooses to
    # pass the derivative of the generalized coordinates as generalized speeds
    q0, q1, q2 = dynamicsymbols('q0:3')
    u0, u1, u2 = dynamicsymbols('q0:3', 1)
    N, A, P, C = _generate_body()
    S = SphericalJoint('S', P, C, coordinates=[q0, q1, q2], speeds=[u0, u1, u2])
    assert S.coordinates == Matrix([q0, q1, q2])
    assert S.speeds == Matrix([u0, u1, u2])
    assert S.kdes == Matrix([0, 0, 0])
    assert P.ang_vel_in(C) == (-u0 * cos(q1) * cos(q2) - u1 * sin(q2)) * A.x + (
        u0 * sin(q2) * cos(q1) - u1 * cos(q2)) * A.y + (
               -u0 * sin(q1) - u2) * A.z


def test_spherical_joint_coords():
    q0s, q1s, q2s, u0s, u1s, u2s = dynamicsymbols('q0:3_S, u0:3_S')
    q0, q1, q2, q3, u0, u1, u2, u4 = dynamicsymbols('q0:4, u0:4')
    # Test coordinates as list
    N, A, P, C = _generate_body()
    S = SphericalJoint('S', P, C, [q0, q1, q2], [u0, u1, u2])
    assert S.coordinates == Matrix([q0, q1, q2])
    assert S.speeds == Matrix([u0, u1, u2])
    # Test coordinates as Matrix
    N, A, P, C = _generate_body()
    S = SphericalJoint('S', P, C, Matrix([q0, q1, q2]),
                       Matrix([u0, u1, u2]))
    assert S.coordinates == Matrix([q0, q1, q2])
    assert S.speeds == Matrix([u0, u1, u2])
    # Test too few generalized coordinates
    N, A, P, C = _generate_body()
    S = SphericalJoint('S', P, C, Matrix([q0, q1]), Matrix([u0]))
    assert S.coordinates == Matrix([q0, q1, q2s])
    assert S.speeds == Matrix([u0, u1s, u2s])
    # Test too many generalized coordinates
    N, A, P, C = _generate_body()
    raises(ValueError, lambda: SphericalJoint(
        'S', P, C, Matrix([q0, q1, q2, q3]), Matrix([u0, u1, u2])))
    raises(ValueError, lambda: SphericalJoint(
        'S', P, C, Matrix([q0, q1, q2]), Matrix([u0, u1, u2, u4])))


def test_spherical_joint_orient_body():
    q0, q1, q2, u0, u1, u2 = dynamicsymbols('q0:3, u0:3')
    N_R_A = Matrix([
        [-sin(q1), -sin(q2) * cos(q1), cos(q1) * cos(q2)],
        [-sin(q0) * cos(q1), sin(q0) * sin(q1) * sin(q2) - cos(q0) * cos(q2),
         -sin(q0) * sin(q1) * cos(q2) - sin(q2) * cos(q0)],
        [cos(q0) * cos(q1), -sin(q0) * cos(q2) - sin(q1) * sin(q2) * cos(q0),
         -sin(q0) * sin(q2) + sin(q1) * cos(q0) * cos(q2)]])
    N_w_A = Matrix([[-u0 * sin(q1) - u2],
                    [-u0 * sin(q2) * cos(q1) + u1 * cos(q2)],
                    [u0 * cos(q1) * cos(q2) + u1 * sin(q2)]])
    N_v_Co = Matrix([
        [-sqrt(2) * (u0 * cos(q2 + pi / 4) * cos(q1) + u1 * sin(q2 + pi / 4))],
        [-u0 * sin(q1) - u2], [-u0 * sin(q1) - u2]])
    # Test default rot_type='BODY', rot_order=123
    N, A, P, C, Pint, Cint = _generate_body(True)
    S = SphericalJoint('S', P, C, coordinates=[q0, q1, q2], speeds=[u0, u1, u2],
                       parent_point=N.x + N.y, child_point=-A.y + A.z,
                       parent_interframe=Pint, child_interframe=Cint,
                       rot_type='body', rot_order=123)
    assert S._rot_type.upper() == 'BODY'
    assert S._rot_order == 123
    assert _simplify_matrix(N.dcm(A) - N_R_A) == zeros(3)
    assert A.ang_vel_in(N).to_matrix(A) == N_w_A
    assert C.masscenter.vel(N).to_matrix(A) == N_v_Co
    # Test change of amounts
    N, A, P, C, Pint, Cint = _generate_body(True)
    S = SphericalJoint('S', P, C, coordinates=[q0, q1, q2], speeds=[u0, u1, u2],
                       parent_point=N.x + N.y, child_point=-A.y + A.z,
                       parent_interframe=Pint, child_interframe=Cint,
                       rot_type='BODY', amounts=(q1, q0, q2), rot_order=123)
    switch_order = lambda expr: expr.xreplace(
        {q0: q1, q1: q0, q2: q2, u0: u1, u1: u0, u2: u2})
    assert S._rot_type.upper() == 'BODY'
    assert S._rot_order == 123
    assert _simplify_matrix(N.dcm(A) - switch_order(N_R_A)) == zeros(3)
    assert A.ang_vel_in(N).to_matrix(A) == switch_order(N_w_A)
    assert C.masscenter.vel(N).to_matrix(A) == switch_order(N_v_Co)
    # Test different rot_order
    N, A, P, C, Pint, Cint = _generate_body(True)
    S = SphericalJoint('S', P, C, coordinates=[q0, q1, q2], speeds=[u0, u1, u2],
                       parent_point=N.x + N.y, child_point=-A.y + A.z,
                       parent_interframe=Pint, child_interframe=Cint,
                       rot_type='BodY', rot_order='yxz')
    assert S._rot_type.upper() == 'BODY'
    assert S._rot_order == 'yxz'
    assert _simplify_matrix(N.dcm(A) - Matrix([
        [-sin(q0) * cos(q1), sin(q0) * sin(q1) * cos(q2) - sin(q2) * cos(q0),
         sin(q0) * sin(q1) * sin(q2) + cos(q0) * cos(q2)],
        [-sin(q1), -cos(q1) * cos(q2), -sin(q2) * cos(q1)],
        [cos(q0) * cos(q1), -sin(q0) * sin(q2) - sin(q1) * cos(q0) * cos(q2),
         sin(q0) * cos(q2) - sin(q1) * sin(q2) * cos(q0)]])) == zeros(3)
    assert A.ang_vel_in(N).to_matrix(A) == Matrix([
        [u0 * sin(q1) - u2], [u0 * cos(q1) * cos(q2) - u1 * sin(q2)],
        [u0 * sin(q2) * cos(q1) + u1 * cos(q2)]])
    assert C.masscenter.vel(N).to_matrix(A) == Matrix([
        [-sqrt(2) * (u0 * sin(q2 + pi / 4) * cos(q1) + u1 * cos(q2 + pi / 4))],
        [u0 * sin(q1) - u2], [u0 * sin(q1) - u2]])


def test_spherical_joint_orient_space():
    q0, q1, q2, u0, u1, u2 = dynamicsymbols('q0:3, u0:3')
    N_R_A = Matrix([
        [-sin(q0) * sin(q2) - sin(q1) * cos(q0) * cos(q2),
         sin(q0) * sin(q1) * cos(q2) - sin(q2) * cos(q0), cos(q1) * cos(q2)],
        [-sin(q0) * cos(q2) + sin(q1) * sin(q2) * cos(q0),
         -sin(q0) * sin(q1) * sin(q2) - cos(q0) * cos(q2), -sin(q2) * cos(q1)],
        [cos(q0) * cos(q1), -sin(q0) * cos(q1), sin(q1)]])
    N_w_A = Matrix([
        [u1 * sin(q0) - u2 * cos(q0) * cos(q1)],
        [u1 * cos(q0) + u2 * sin(q0) * cos(q1)], [u0 - u2 * sin(q1)]])
    N_v_Co = Matrix([
        [u0 - u2 * sin(q1)], [u0 - u2 * sin(q1)],
        [sqrt(2) * (-u1 * sin(q0 + pi / 4) + u2 * cos(q0 + pi / 4) * cos(q1))]])
    # Test default rot_type='BODY', rot_order=123
    N, A, P, C, Pint, Cint = _generate_body(True)
    S = SphericalJoint('S', P, C, coordinates=[q0, q1, q2], speeds=[u0, u1, u2],
                       parent_point=N.x + N.z, child_point=-A.x + A.y,
                       parent_interframe=Pint, child_interframe=Cint,
                       rot_type='space', rot_order=123)
    assert S._rot_type.upper() == 'SPACE'
    assert S._rot_order == 123
    assert _simplify_matrix(N.dcm(A) - N_R_A) == zeros(3)
    assert _simplify_matrix(A.ang_vel_in(N).to_matrix(A)) == N_w_A
    assert _simplify_matrix(C.masscenter.vel(N).to_matrix(A)) == N_v_Co
    # Test change of amounts
    switch_order = lambda expr: expr.xreplace(
        {q0: q1, q1: q0, q2: q2, u0: u1, u1: u0, u2: u2})
    N, A, P, C, Pint, Cint = _generate_body(True)
    S = SphericalJoint('S', P, C, coordinates=[q0, q1, q2], speeds=[u0, u1, u2],
                       parent_point=N.x + N.z, child_point=-A.x + A.y,
                       parent_interframe=Pint, child_interframe=Cint,
                       rot_type='SPACE', amounts=(q1, q0, q2), rot_order=123)
    assert S._rot_type.upper() == 'SPACE'
    assert S._rot_order == 123
    assert _simplify_matrix(N.dcm(A) - switch_order(N_R_A)) == zeros(3)
    assert _simplify_matrix(A.ang_vel_in(N).to_matrix(A)) == switch_order(N_w_A)
    assert _simplify_matrix(C.masscenter.vel(N).to_matrix(A)) == switch_order(N_v_Co)
    # Test different rot_order
    N, A, P, C, Pint, Cint = _generate_body(True)
    S = SphericalJoint('S', P, C, coordinates=[q0, q1, q2], speeds=[u0, u1, u2],
                       parent_point=N.x + N.z, child_point=-A.x + A.y,
                       parent_interframe=Pint, child_interframe=Cint,
                       rot_type='SPaCe', rot_order='zxy')
    assert S._rot_type.upper() == 'SPACE'
    assert S._rot_order == 'zxy'
    assert _simplify_matrix(N.dcm(A) - Matrix([
        [-sin(q2) * cos(q1), -sin(q0) * cos(q2) + sin(q1) * sin(q2) * cos(q0),
         sin(q0) * sin(q1) * sin(q2) + cos(q0) * cos(q2)],
        [-sin(q1), -cos(q0) * cos(q1), -sin(q0) * cos(q1)],
        [cos(q1) * cos(q2), -sin(q0) * sin(q2) - sin(q1) * cos(q0) * cos(q2),
         -sin(q0) * sin(q1) * cos(q2) + sin(q2) * cos(q0)]]))
    assert _simplify_matrix(A.ang_vel_in(N).to_matrix(A) - Matrix([
        [-u0 + u2 * sin(q1)], [-u1 * sin(q0) + u2 * cos(q0) * cos(q1)],
        [u1 * cos(q0) + u2 * sin(q0) * cos(q1)]])) == zeros(3, 1)
    assert _simplify_matrix(C.masscenter.vel(N).to_matrix(A) - Matrix([
        [u1 * cos(q0) + u2 * sin(q0) * cos(q1)],
        [u1 * cos(q0) + u2 * sin(q0) * cos(q1)],
        [u0 + u1 * sin(q0) - u2 * sin(q1) -
         u2 * cos(q0) * cos(q1)]])) == zeros(3, 1)


def test_weld_joint():
    _, _, P, C = _generate_body()
    W = WeldJoint('W', P, C)
    assert W.name == 'W'
    assert W.parent == P
    assert W.child == C
    assert W.coordinates == Matrix()
    assert W.speeds == Matrix()
    assert W.kdes == Matrix(1, 0, []).T
    assert P.dcm(C) == eye(3)
    assert W.child_point.pos_from(C.masscenter) == Vector(0)
    assert W.parent_point.pos_from(P.masscenter) == Vector(0)
    assert W.parent_point.pos_from(W.child_point) == Vector(0)
    assert P.masscenter.pos_from(C.masscenter) == Vector(0)
    assert C.masscenter.vel(P.frame) == Vector(0)
    assert P.ang_vel_in(C) == 0
    assert C.ang_vel_in(P) == 0
    assert W.__str__() == 'WeldJoint: W  parent: P  child: C'

    N, A, P, C = _generate_body()
    l, m = symbols('l m')
    Pint = ReferenceFrame('P_int')
    Pint.orient_axis(P.frame, P.y, pi / 2)
    W = WeldJoint('W', P, C, parent_point=l * P.frame.x,
                  child_point=m * C.frame.y, parent_interframe=Pint)

    assert W.child_point.pos_from(C.masscenter) == m * C.frame.y
    assert W.parent_point.pos_from(P.masscenter) == l * P.frame.x
    assert W.parent_point.pos_from(W.child_point) == Vector(0)
    assert P.masscenter.pos_from(C.masscenter) == - l * N.x + m * A.y
    assert C.masscenter.vel(P.frame) == Vector(0)
    assert P.masscenter.vel(Pint) == Vector(0)
    assert C.ang_vel_in(P) == 0
    assert P.ang_vel_in(C) == 0
    assert P.x == A.z

    JointsMethod(P, W)  # Tests #10770


def test_deprecated_parent_child_axis():
    q, u = dynamicsymbols('q_J, u_J')
    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PinJoint('J', P, C, child_axis=-A.x)
    assert (-A.x).angle_between(N.x) == 0
    assert -A.x.express(N) == N.x
    assert A.dcm(N) == Matrix([[-1, 0, 0],
                               [0, -cos(q), -sin(q)],
                               [0, -sin(q), cos(q)]])
    assert A.ang_vel_in(N) == u * N.x
    assert A.ang_vel_in(N).magnitude() == sqrt(u ** 2)

    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PrismaticJoint('J', P, C, parent_axis=P.x + P.y)
    assert (A.x).angle_between(N.x + N.y) == 0
    assert A.x.express(N) == (N.x + N.y) / sqrt(2)
    assert A.dcm(N) == Matrix([[sqrt(2) / 2, sqrt(2) / 2, 0],
                               [-sqrt(2) / 2, sqrt(2) / 2, 0], [0, 0, 1]])
    assert A.ang_vel_in(N) == Vector(0)


def test_deprecated_joint_pos():
    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        pin = PinJoint('J', P, C, parent_joint_pos=N.x + N.y,
                       child_joint_pos=C.y - C.z)
    assert pin.parent_point.pos_from(P.masscenter) == N.x + N.y
    assert pin.child_point.pos_from(C.masscenter) == C.y - C.z

    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        slider = PrismaticJoint('J', P, C, parent_joint_pos=N.z + N.y,
                                child_joint_pos=C.y - C.x)
    assert slider.parent_point.pos_from(P.masscenter) == N.z + N.y
    assert slider.child_point.pos_from(C.masscenter) == C.y - C.x
