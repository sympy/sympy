"""Tests for the ``sympy.physics.mechanics._geometry.py`` module."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from sympy.core.backend import Integer, S, Symbol
from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.mechanics._geometry import Sphere

if TYPE_CHECKING:
    from sympy.physics.mechanics import Vector


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
