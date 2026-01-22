"""Tests for the ``sympy.physics.mechanics.wrapping_geometry.py`` module."""

import pytest

from sympy import (
    Integer,
    Rational,
    S,
    Symbol,
    acos,
    cos,
    pi,
    sin,
    sqrt,
    diff,
    Piecewise,
    symbols,
)
from sympy.core.relational import Eq
from sympy.physics.mechanics import (
    Point,
    ReferenceFrame,
    WrappingCylinder,
    WrappingSphere,
    WrappingCone,
    dynamicsymbols,
    LagrangesMethod,
    Lagrangian,
    Particle,
)
from sympy.simplify.simplify import simplify
from sympy.polys.polytools import cancel


r = Symbol('r', positive=True)
x = Symbol('x')
q = dynamicsymbols('q')
N = ReferenceFrame('N')


class TestWrappingSphere:

    @staticmethod
    def test_valid_constructor():
        r = Symbol('r', positive=True)
        pO = Point('pO')
        sphere = WrappingSphere(r, pO)
        assert isinstance(sphere, WrappingSphere)
        assert hasattr(sphere, 'radius')
        assert sphere.radius == r
        assert hasattr(sphere, 'point')
        assert sphere.point == pO

    @staticmethod
    @pytest.mark.parametrize(
        'position_1, position_2, expected',
        [
            (r*N.x, r*N.x, S.Zero),
            (r*N.x, r*N.y, S.Half*pi*r),
            (r*N.x, r*-N.x, pi*r),
            (r*-N.x, r*N.x, pi*r),
            (r*N.x, r*sqrt(2)*S.Half*(N.x + N.y), Rational(1, 4)*pi*r),
            (
                r*sqrt(2)*S.Half*(N.x + N.y),
                r*sqrt(3)*Rational(1, 3)*(N.x + N.y + N.z),
                r*acos(sqrt(6)*Rational(1, 3)),
            ),
        ]
    )
    def test_geodesic_length(position_1, position_2, expected):
        r = Symbol('r', positive=True)
        pO = Point('pO')
        sphere = WrappingSphere(r, pO)

        p1 = Point('p1')
        p1.set_pos(pO, position_1)
        p2 = Point('p2')
        p2.set_pos(pO, position_2)

        assert simplify(Eq(sphere.geodesic_length(p1, p2), expected))

    @staticmethod
    @pytest.mark.parametrize(
        'position_1, position_2, vector_1, vector_2',
        [
            (r * N.x, r * N.y, N.y, N.x),
            (r * N.x, -r * N.y, -N.y, N.x),
            (
                r * N.y,
                sqrt(2)/2 * r * N.x - sqrt(2)/2 * r * N.y,
                N.x,
                sqrt(2)/2 * N.x + sqrt(2)/2 * N.y,
            ),
            (
                r * N.x,
                r / 2 * N.x + sqrt(3)/2 * r * N.y,
                N.y,
                sqrt(3)/2 * N.x - 1/2 * N.y,
            ),
            (
                r * N.x,
                sqrt(2)/2 * r * N.x + sqrt(2)/2 * r * N.y,
                N.y,
                sqrt(2)/2 * N.x - sqrt(2)/2 * N.y,
            ),
        ]
    )
    def test_geodesic_end_vectors(position_1, position_2, vector_1, vector_2):
        r = Symbol('r', positive=True)
        pO = Point('pO')
        sphere = WrappingSphere(r, pO)

        p1 = Point('p1')
        p1.set_pos(pO, position_1)
        p2 = Point('p2')
        p2.set_pos(pO, position_2)

        expected = (vector_1, vector_2)

        assert sphere.geodesic_end_vectors(p1, p2) == expected

    @staticmethod
    @pytest.mark.parametrize(
        'position',
        [r * N.x, r * cos(q) * N.x + r * sin(q) * N.y]
    )
    def test_geodesic_end_vectors_invalid_coincident(position):
        r = Symbol('r', positive=True)
        pO = Point('pO')
        sphere = WrappingSphere(r, pO)

        p1 = Point('p1')
        p1.set_pos(pO, position)
        p2 = Point('p2')
        p2.set_pos(pO, position)

        with pytest.raises(ValueError):
            _ = sphere.geodesic_end_vectors(p1, p2)

    @staticmethod
    @pytest.mark.parametrize(
        'position_1, position_2',
        [
            (r * N.x, -r * N.x),
            (-r * N.y, r * N.y),
            (
                r * cos(q) * N.x + r * sin(q) * N.y,
                -r * cos(q) * N.x - r * sin(q) * N.y,
            )
        ]
    )
    def test_geodesic_end_vectors_invalid_diametrically_opposite(
        position_1,
        position_2,
    ):
        r = Symbol('r', positive=True)
        pO = Point('pO')
        sphere = WrappingSphere(r, pO)

        p1 = Point('p1')
        p1.set_pos(pO, position_1)
        p2 = Point('p2')
        p2.set_pos(pO, position_2)

        with pytest.raises(ValueError):
            _ = sphere.geodesic_end_vectors(p1, p2)


class TestWrappingCylinder:

    @staticmethod
    def test_valid_constructor():
        N = ReferenceFrame('N')
        r = Symbol('r', positive=True)
        pO = Point('pO')
        cylinder = WrappingCylinder(r, pO, N.x)
        assert isinstance(cylinder, WrappingCylinder)
        assert hasattr(cylinder, 'radius')
        assert cylinder.radius == r
        assert hasattr(cylinder, 'point')
        assert cylinder.point == pO
        assert hasattr(cylinder, 'axis')
        assert cylinder.axis == N.x

    @staticmethod
    @pytest.mark.parametrize(
        'position, expected',
        [
            (S.Zero, Eq(0, r**2, evaluate=False)),
            (r*N.y, Eq(r**2, r**2, evaluate=False)),
            (r*N.z, Eq(r**2, r**2, evaluate=False)),
            (r*(N.y + N.z).normalize(), Eq(r**2, r**2, evaluate=False)),
            (Integer(2)*r*N.y, Eq(4*r**2, r**2, evaluate=False)),
            (r*(N.x + N.y), Eq(r**2, r**2, evaluate=False)),
            (r*(Integer(2)*N.x + N.y), Eq(r**2, r**2, evaluate=False)),
            (Integer(2)*N.x + r*(Integer(2)*N.y + N.z).normalize(), Eq(r**2, r**2, evaluate=False)),
            (r*(cos(q)*N.y + sin(q)*N.z), Eq(r**2, r**2, evaluate=False))
        ]
    )
    def test_point_is_on_surface(position, expected):
        r = Symbol('r', positive=True)
        pO = Point('pO')
        cylinder = WrappingCylinder(r, pO, N.x)

        p1 = Point('p1')
        p1.set_pos(pO, position)

        assert cylinder.point_on_surface(p1) == expected

    @staticmethod
    @pytest.mark.parametrize(
        'axis, position_1, position_2, expected',
        [
            (N.x, r*N.y, r*N.y, S.Zero),
            (N.x, r*N.y, N.x + r*N.y, S.One),
            (N.x, r*N.y, -x*N.x + r*N.y, sqrt(x**2)),
            (-N.x, r*N.y, x*N.x + r*N.y, sqrt(x**2)),
            (N.x, r*N.y, r*N.z, S.Half*pi*sqrt(r**2)),
            (-N.x, r*N.y, r*N.z, Integer(3)*S.Half*pi*sqrt(r**2)),
            (N.x, r*N.z, r*N.y, Integer(3)*S.Half*pi*sqrt(r**2)),
            (-N.x, r*N.z, r*N.y, S.Half*pi*sqrt(r**2)),
            (N.x, r*N.y, r*(cos(q)*N.y + sin(q)*N.z), sqrt(r**2*q**2)),
            (
                -N.x, r*N.y,
                r*(cos(q)*N.y + sin(q)*N.z),
                sqrt(r**2*(Integer(2)*pi - q)**2),
            ),
        ]
    )
    def test_geodesic_length(axis, position_1, position_2, expected):
        r = Symbol('r', positive=True)
        pO = Point('pO')
        cylinder = WrappingCylinder(r, pO, axis)

        p1 = Point('p1')
        p1.set_pos(pO, position_1)
        p2 = Point('p2')
        p2.set_pos(pO, position_2)

        assert simplify(Eq(cylinder.geodesic_length(p1, p2), expected))

    @staticmethod
    @pytest.mark.parametrize(
        'axis, position_1, position_2, vector_1, vector_2',
        [
            (N.z, r * N.x, r * N.y, N.y, N.x),
            (N.z, r * N.x, -r * N.x, N.y, N.y),
            (N.z, -r * N.x, r * N.x, -N.y, -N.y),
            (-N.z, r * N.x, -r * N.x, -N.y, -N.y),
            (-N.z, -r * N.x, r * N.x, N.y, N.y),
            (N.z, r * N.x, -r * N.y, N.y, -N.x),
            (
                N.z,
                r * N.y,
                sqrt(2)/2 * r * N.x - sqrt(2)/2 * r * N.y,
                - N.x,
                - sqrt(2)/2 * N.x - sqrt(2)/2 * N.y,
            ),
            (
                N.z,
                r * N.x,
                r / 2 * N.x + sqrt(3)/2 * r * N.y,
                N.y,
                sqrt(3)/2 * N.x - 1/2 * N.y,
            ),
            (
                N.z,
                r * N.x,
                sqrt(2)/2 * r * N.x + sqrt(2)/2 * r * N.y,
                N.y,
                sqrt(2)/2 * N.x - sqrt(2)/2 * N.y,
            ),
            (
                N.z,
                r * N.x,
                r * N.x + N.z,
                N.z,
                -N.z,
            ),
            (
                N.z,
                r * N.x,
                r * N.y + pi/2 * r * N.z,
                sqrt(2)/2 * N.y + sqrt(2)/2 * N.z,
                sqrt(2)/2 * N.x - sqrt(2)/2 * N.z,
            ),
            (
                N.z,
                r * N.x,
                r * cos(q) * N.x + r * sin(q) * N.y,
                N.y,
                sin(q) * N.x - cos(q) * N.y,
            ),
        ]
    )
    def test_geodesic_end_vectors(
        axis,
        position_1,
        position_2,
        vector_1,
        vector_2,
    ):
        r = Symbol('r', positive=True)
        pO = Point('pO')
        cylinder = WrappingCylinder(r, pO, axis)

        p1 = Point('p1')
        p1.set_pos(pO, position_1)
        p2 = Point('p2')
        p2.set_pos(pO, position_2)

        expected = (vector_1, vector_2)
        end_vectors = tuple(
            end_vector.simplify()
            for end_vector in cylinder.geodesic_end_vectors(p1, p2)
        )

        assert end_vectors == expected

    @staticmethod
    @pytest.mark.parametrize(
        'axis, position',
        [
            (N.z, r * N.x),
            (N.z, r * cos(q) * N.x + r * sin(q) * N.y + N.z),
        ]
    )
    def test_geodesic_end_vectors_invalid_coincident(axis, position):
        r = Symbol('r', positive=True)
        pO = Point('pO')
        cylinder = WrappingCylinder(r, pO, axis)

        p1 = Point('p1')
        p1.set_pos(pO, position)
        p2 = Point('p2')
        p2.set_pos(pO, position)

        with pytest.raises(ValueError):
            _ = cylinder.geodesic_end_vectors(p1, p2)


class TestWrappingCone:

    @staticmethod
    def test_valid_constructor():
        N = ReferenceFrame('N')
        alpha, apex, axis = Symbol('alpha'), Point('p0'), N.z
        cone = WrappingCone(alpha, apex, axis)
        assert isinstance(cone, WrappingCone)
        assert hasattr(cone, 'alpha')
        assert cone.alpha == alpha
        assert hasattr(cone, 'apex')
        assert cone.apex == apex
        assert hasattr(cone, 'axis')
        assert cone.axis == axis

    @staticmethod
    @pytest.mark.parametrize(
        'position, expected',
        [
            (S.Zero, Eq(0, 0, evaluate=False)),
            (N.x + N.y + N.z, Eq(2, Rational(1, 3), evaluate=False)),
            (N.x + sqrt(3) * N.z, Eq(1, 1, evaluate=False)),
            (N.y + sqrt(3) * N.z, Eq(1, 1, evaluate=False)),
            ((N.x + N.y) / sqrt(2) + sqrt(3) * N.z, Eq(1, 1, evaluate=False)),
            (N.x / sqrt(3) + sqrt(2) * N.y / sqrt(3) + sqrt(3) * N.z, Eq(1, 1, evaluate=False)),
            (2 * N.x + sqrt(12) * N.z, Eq(4, 4, evaluate=False)),
            (2 * N.y + sqrt(12) * N.z, Eq(4, 4, evaluate=False)),
            (5 * N.x + sqrt(12) * N.z, Eq(25, 4, evaluate=False)),
            (sqrt(2) * (N.x + N.y) + sqrt(12) * N.z, Eq(4, 4, evaluate=False))
        ]
    )
    def test_point_on_surface(position, expected):
        axis = N.z
        apex = Point('p0')
        alpha = pi/6
        cone = WrappingCone(alpha, apex, axis)

        p1 = Point('p1')
        p1.set_pos(apex, position)

        assert cone.point_on_surface(p1) == expected

    @staticmethod
    @pytest.mark.parametrize(
        'axis, alpha, position_1, position_2, expected',
        [
            (N.z, pi/4, (N.x + N.z)/sqrt(2), (N.y + N.z)/sqrt(2), sqrt(2 - 2*cos(pi/(2*sqrt(2))))),
            (N.z, pi/6, N.x/sqrt(3) + N.z, N.y/sqrt(3) + N.z, sqrt(Rational(8, 3) - 4*sqrt(2)/3)),
            (N.z, pi/4, (N.x + N.z)/sqrt(2), (2*N.x + 2*N.z)/sqrt(2), 1),
            (N.z, pi/4, (N.x + N.z)/sqrt(2), (-N.x + N.z)/sqrt(2), sqrt(2 - 2*cos(pi/sqrt(2)))),
            (N.x, pi/3, (N.y + N.x)/2, (N.z + N.x)/2, sqrt(2 - 2*cos(pi*sqrt(3)/4))),
            (N.z, pi/6, (N.x + N.z*sqrt(3))/2, (N.y + N.z*sqrt(3))/2, sqrt(2 - sqrt(2))),
            (N.z, pi/3, (N.x*sqrt(3) + N.z)/2, (N.y*sqrt(3) + N.z)/2, sqrt(2 - 2*cos(pi*sqrt(3)/4))),
            (N.z, pi/4, (N.x + N.z)/sqrt(11), (N.x + N.z)/sqrt(11), 0),
            (N.z, pi/6, N.x/sqrt(3) + N.z, 2*N.y/sqrt(3) + 2*N.z, sqrt(Rational(20, 3) - 8*sqrt(2)/3)),
            (N.z, pi/6, (N.x + N.y)/(sqrt(2)*sqrt(3)) + N.z, (3*N.x - 3*N.y)/(sqrt(2)*sqrt(3)) + 3*N.z, sqrt(Rational(40, 3) - 4*sqrt(2))),
            (N.z, pi/4, (cos(pi/6)*N.x + sin(pi/6)*N.y + N.z)/sqrt(2), (3*cos(2*pi/3)*N.x + 3*sin(2*pi/3)*N.y + 3*N.z)/sqrt(2), sqrt(10 - 6*cos(pi/(2*sqrt(2))))),
            (N.z, pi/4, (N.x + N.z)/sqrt(2), (2*N.y + 2*N.z)/sqrt(2), sqrt(5 - 4*cos(pi/(2*sqrt(2))))),
        ]
    )
    def test_geodesic_length(axis, alpha, position_1, position_2, expected):
        apex = Point('p0')
        cone = WrappingCone(alpha, apex, axis)

        p1 = Point('p1')
        p1.set_pos(apex, position_1)

        p2 = Point('p2')
        p2.set_pos(apex, position_2)

        assert cone.geodesic_length(p1, p2) == expected


    @staticmethod
    @pytest.mark.parametrize(
        'axis, alpha, position_1, position_2',
        [
            (N.z, pi/4, N.x + N.z, N.y + N.z),
            (N.z, pi/6, N.x/sqrt(3) + N.z, 2*N.x/sqrt(3) + 2*N.z),
            (N.x, pi/4, N.y + N.x, N.z + 2*N.x),
            (N.z, pi/6, N.x/sqrt(3) + N.z, 2*N.y/sqrt(3) + 2*N.z),
        ]
    )
    def test_geodesic_end_vectors(axis, alpha, position_1, position_2):
        apex = Point('p0')
        cone = WrappingCone(alpha, apex, axis)
        p1 = Point('p1'); p1.set_pos(apex, position_1)
        p2 = Point('p2'); p2.set_pos(apex, position_2)

        v1_calc, v2_calc = cone.geodesic_end_vectors(p1, p2)

        pos1 = p1.pos_from(apex)
        pos2 = p2.pos_from(apex)
        z1 = pos1.dot(axis)
        z2 = pos2.dot(axis)
        s1 = z1 / cos(alpha)
        s2 = z2 / cos(alpha)
        L = cone.geodesic_length(p1, p2)
        n1 = (pos1 - z1*axis).normalize()
        n2 = (pos2 - z2*axis).normalize()

        central = acos(cancel(n1.dot(n2)))
        delta_u = central * sin(alpha)

        g1 = pos1.normalize()
        c1 = axis.cross(n1)
        g2 = pos2.normalize()
        c2 = axis.cross(n2)

        v1_radial_comp = (s2 * cos(delta_u) - s1) / L
        v1_circ_comp = (s2 * sin(delta_u)) / L
        expected_v1 = v1_radial_comp * g1 + v1_circ_comp * c1

        v2_radial_comp = (s1 * cos(delta_u) - s2) / L
        v2_circ_comp = (-s1 * sin(delta_u)) / L
        expected_v2 = v2_radial_comp * g2 + v2_circ_comp * c2

        assert v1_calc == expected_v1
        assert v2_calc == expected_v2

    @staticmethod
    @pytest.mark.parametrize(
        'position',
        [
            (N.x + N.z),
            (N.y/sqrt(3) + 2*N.z),
        ]
    )
    def test_geodesic_end_vectors_invalid_coincident(position):
        apex = Point('p0')
        cone = WrappingCone(pi/6, apex, N.z)

        p1 = Point('p1')
        p1.set_pos(apex, position)

        with pytest.raises(ValueError,
                           match='No unique geodesic exists for coincident points'):
            cone.geodesic_end_vectors(p1, p1)


class TestElasticConeModel:
    """
    A test class to verify the Lagrangian mechanics model of a particle
    tethered by an elastic cable that is constrained to move on the surface
    of a fixed cone.

    The physical system consists of a particle of mass `m` which can slide
    freely on a cone with a half-angle `alpha`. This particle is attached to a
    fixed anchor point `A` (also on the cone) by a massless, elastic cable
    with a spring constant `k` and a natural rest length `L0`.

    A key feature of this model is that the cable always follows the shortest
    path, the geodesic, between the anchor point and the particle. The system has
    two degrees of freedom, described by the generalized coordinates `r(t)`
    (the slant distance from the cone's apex) and `phi(t)` (the azimuthal angle).
    """

    @classmethod
    def setup_class(cls):
        # Setup basic symbols
        cls.t = symbols('t')
        cls.m, cls.g, cls.alpha, cls.k, cls.L0 = symbols('m g alpha k L0', positive=True)
        cls.rA, cls.phiA = symbols('rA phiA', real=True)
        cls.r, cls.phi = dynamicsymbols('r phi')
        cls.dr, cls.dphi = dynamicsymbols('r phi', 1)

        cls.N = ReferenceFrame('N')
        cls.O = Point('O')
        cls.O.set_vel(cls.N, 0)
        cls.cone = WrappingCone(cls.alpha, cls.O, cls.N.z)

        # Generalized point on the cone
        P = Point('P')
        x = cls.r * sin(cls.alpha) * cos(cls.phi)
        y = cls.r * sin(cls.alpha) * sin(cls.phi)
        z = cls.r * cos(cls.alpha)

        P.set_pos(cls.O, x * cls.N.x + y * cls.N.y + z * cls.N.z)
        P.set_vel(cls.N, P.pos_from(cls.O).diff(cls.t, cls.N))
        P_part = Particle('P_part', P, cls.m)

        A = Point('A')
        xA = cls.rA * sin(cls.alpha) * cos(cls.phiA)
        yA = cls.rA * sin(cls.alpha) * sin(cls.phiA)
        zA = cls.rA * cos(cls.alpha)
        A.set_pos(cls.O, xA * cls.N.x + yA * cls.N.y + zA * cls.N.z)

        cls.T = P_part.kinetic_energy(cls.N)
        cls.L_geo = cls.cone.geodesic_length(A, P)

        # Gravitational potential energy
        h = P.pos_from(cls.O).dot(cls.N.z)
        Vg = cls.m * cls.g * h

        # Elastic potential energy
        Ve = Rational(1, 2) * cls.k * (cls.L_geo - cls.L0)**2

        # Total potential energy
        cls.V = Vg + Ve

        P_part.potential_energy = cls.V

        # Get EOMs
        cls.L = Lagrangian(cls.N, P_part)
        LM = LagrangesMethod(cls.L, [cls.r, cls.phi])
        cls.eom = LM.form_lagranges_equations()
        for i in range(len(cls.eom)):
            cls.eom[i] = simplify(cls.eom[i])

    def test_kinetic_energy(self):
        expected_T = Rational(1, 2) * self.m * (self.dr**2 + self.r**2 * sin(self.alpha)**2 * self.dphi**2)
        assert simplify(self.T - expected_T) == 0

    def test_geodesic_length(self):
        expected_L_geo = sqrt(self.rA**2 - 2*self.rA*self.r*cos(Piecewise((-self.phiA + self.phi + 2*pi, -self.phiA + self.phi + 2*pi <= pi), (self.phiA - self.phi, True))*sin(self.alpha)) + self.r**2)
        assert simplify(self.L_geo - expected_L_geo) == 0

    def test_potential_energy(self):
        expected_V = (self.m * self.g * self.r * cos(self.alpha) + Rational(1, 2) * self.k * (self.L_geo - self.L0)**2)
        assert simplify(self.V - expected_V) == 0

    def test_lagrangian(self):
        assert simplify(self.L - (self.T - self.V)) == 0

    def test_equation_of_motion_r(self):
        # Derived by hand: EOM for r: m*r_ddot - m*r*sin(alpha)**2*phi_dot**2 + dV/dr = 0
        dV_dr = diff(self.V, self.r)
        expected_eom_r = (self.m * diff(self.r, self.t, 2) -
                          self.m * self.r * sin(self.alpha)**2 * self.dphi**2 +
                          dV_dr)

        assert simplify(self.eom[0] - expected_eom_r) == 0

    def test_equation_of_motion_phi(self):
        # Derived by hand: EOM for phi: d/dt(m*r**2*sin(alpha)**2*phi_dot) + dV/dphi = 0
        dV_dphi = diff(self.V, self.phi)
        expected_eom_phi = (diff(self.m * self.r**2 * sin(self.alpha)**2 * self.dphi, self.t) + dV_dphi)

        assert simplify(self.eom[1] - expected_eom_phi) == 0
