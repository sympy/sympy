"""Tests for the ``sympy.physics.mechanics._geometry.py`` module."""

from __future__ import annotations

from sympy.core.backend import Symbol
from sympy.physics.mechanics import Point
from sympy.physics.mechanics._geometry import Sphere


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
