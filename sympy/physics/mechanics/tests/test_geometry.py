"""Tests for the ``sympy.physics.mechanics._geometry.py`` module."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from sympy.core.backend import (
    Integer,
    Rational,
    S,
    Symbol,
    acos,
    cos,
    pi,
    sin,
    sqrt,
)
from sympy.core.relational import Eq
from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.mechanics._geometry import Cylinder, Sphere
from sympy.simplify.simplify import simplify

if TYPE_CHECKING:
    from sympy.core.backend import USE_SYMENGINE
    from sympy.physics.mechanics import Vector

    if USE_SYMENGINE:
        from sympy.core.backend import Basic as ExprType
    else:
        from sympy.core.expr import Expr as ExprType


r = Symbol('r')
x = Symbol('x')
q = dynamicsymbols('q')
N = ReferenceFrame('N')


class TestSphere:

    @staticmethod
    def test_valid_constructor() -> None:
        r = Symbol('r')
        pO = Point('pO')
        sphere = Sphere(r, pO)
        assert isinstance(sphere, Sphere)
        assert hasattr(sphere, 'radius')
        assert sphere.radius == r
        assert hasattr(sphere, 'point')
        assert sphere.point == pO

    @staticmethod
    @pytest.mark.parametrize('position', [S.Zero, Integer(2) * r * N.x])
    def test_geodesic_length_point_not_on_surface_invalid(position: Vector) -> None:
        r = Symbol('r')
        pO = Point('pO')
        sphere = Sphere(r, pO)

        p1 = Point('p1')
        p1.set_pos(pO, position)
        p2 = Point('p2')
        p2.set_pos(pO, position)

        error_msg = r'point .* does not lie on the surface of'
        with pytest.raises(ValueError, match=error_msg):
            sphere.geodesic_length(p1, p2)

    @staticmethod
    @pytest.mark.parametrize(
        'position_1, position_2, expected',
        [
            (r * N.x, r * N.x, S.Zero),
            (r * N.x, r * N.y, S.Half * pi * r),
            (r * N.x, r * -N.x, pi * r),
            (r * -N.x, r * N.x, pi * r),
            (r * N.x, r * sqrt(2) * S.Half * (N.x + N.y), Rational(1, 4) * pi * r),
            (
                r * sqrt(2) * S.Half * (N.x + N.y),
                r * sqrt(3) * Rational(1, 3) * (N.x + N.y + N.z),
                r * acos(sqrt(6) * Rational(1, 3)),
            ),
        ]
    )
    def test_geodesic_length(position_1: Vector, position_2: Vector, expected: ExprType) -> None:
        r = Symbol('r')
        pO = Point('pO')
        sphere = Sphere(r, pO)

        p1 = Point('p1')
        p1.set_pos(pO, position_1)
        p2 = Point('p2')
        p2.set_pos(pO, position_2)

        assert simplify(Eq(sphere.geodesic_length(p1, p2), expected))


class TestCylinder:

    @staticmethod
    def test_valid_constructor() -> None:
        N = ReferenceFrame('N')
        r = Symbol('r')
        pO = Point('pO')
        cylinder = Cylinder(r, pO, N.x)
        assert isinstance(cylinder, Cylinder)
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
            (S.Zero, False),
            (r * N.y, True),
            (r * N.z, True),
            (r * (N.y + N.z).normalize(), True),
            (Integer(2) * r * N.y, False),
            (r * (N.x + N.y), True),
            (r * (Integer(2) * N.x + N.y), True),
            (Integer(2) * N.x + r * (Integer(2) * N.y + N.z).normalize(), True),
            (r * (cos(q) * N.y + sin(q) * N.z), True)
        ]
    )
    def test_point_is_on_surface(position: Vector, expected: bool) -> None:
        r = Symbol('r')
        pO = Point('pO')
        cylinder = Cylinder(r, pO, N.x)

        p1 = Point('p1')
        p1.set_pos(pO, position)

        assert cylinder._point_is_on_surface(p1) is expected

    @staticmethod
    @pytest.mark.parametrize('position', [S.Zero, Integer(2) * r * N.y])
    def test_geodesic_length_point_not_on_surface_invalid(position: Vector) -> None:
        r = Symbol('r')
        pO = Point('pO')
        cylinder = Cylinder(r, pO, N.x)

        p1 = Point('p1')
        p1.set_pos(pO, position)
        p2 = Point('p2')
        p2.set_pos(pO, position)

        error_msg = r'point .* does not lie on the surface of'
        with pytest.raises(ValueError, match=error_msg):
            cylinder.geodesic_length(p1, p2)

    @staticmethod
    @pytest.mark.parametrize(
        'axis, position_1, position_2, expected',
        [
            (N.x, r * N.y, r * N.y, S.Zero),
            (N.x, r * N.y, N.x + r * N.y, S.One),
            (N.x, r * N.y, -x * N.x + r * N.y, sqrt(x**2)),
            (-N.x, r * N.y, x * N.x + r * N.y, sqrt(x**2)),
            (N.x, r * N.y, r * N.z, S.Half * pi * sqrt(r**2)),
            (-N.x, r * N.y, r * N.z, Integer(3) * S.Half * pi * sqrt(r**2)),
            (N.x, r * N.z, r * N.y, Integer(3) * S.Half * pi * sqrt(r**2)),
            (-N.x, r * N.z, r * N.y, S.Half * pi * sqrt(r**2)),
            (N.x, r * N.y, r * (cos(q) * N.y + sin(q) * N.z), sqrt(r**2 * q**2)),
            (
                -N.x, r * N.y,
                r * (cos(q) * N.y + sin(q) * N.z),
                sqrt(r**2 * (Integer(2) * pi - q)**2),
            ),
        ]
    )
    def test_geodesic_length(
        axis: Vector,
        position_1: Vector,
        position_2: Vector,
        expected: ExprType,
    ) -> None:
        r = Symbol('r')
        pO = Point('pO')
        cylinder = Cylinder(r, pO, axis)

        p1 = Point('p1')
        p1.set_pos(pO, position_1)
        p2 = Point('p2')
        p2.set_pos(pO, position_2)

        assert simplify(Eq(cylinder.geodesic_length(p1, p2), expected))
