"""Tests for the ``sympy.physics.mechanics._pathway.py`` module."""

from __future__ import annotations

from typing import Any, Sequence

import pytest

from sympy.core.backend import Symbol, sqrt
from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.mechanics._pathway import LinearPathway
from sympy.simplify.simplify import simplify


class TestLinearPathway:

    def test_is_pathway_base_subclass(self) -> None:
        assert issubclass(LinearPathway, PathwayBase)

    @staticmethod
    @pytest.mark.parametrize(
        'args, kwargs',
        [
            ((Point('pA'), Point('pB')), {}),
        ]
    )
    def test_valid_constructor(args: tuple, kwargs: dict) -> None:
        instance = LinearPathway(*args, **kwargs)
        assert isinstance(instance, LinearPathway)
        assert hasattr(instance, 'attachments')
        assert len(instance.attachments) == 2
        assert isinstance(instance.attachments[0], Point)
        assert instance.attachments[0].name == 'pA'
        assert isinstance(instance.attachments[1], Point)
        assert instance.attachments[1].name == 'pB'

    @staticmethod
    @pytest.mark.parametrize(
        'attachments',
        [
            (Point('pA'), ),
            (Point('pA'), Point('pB'), Point('pZ')),
        ]
    )
    def test_invalid_attachments_incorrect_number(
        attachments: tuple[Point, ...],
    ) -> None:
        with pytest.raises(ValueError):
            _ = LinearPathway(*attachments)  # type: ignore

    @staticmethod
    @pytest.mark.parametrize(
        'attachments',
        [
            (None, Point('pB')),
            (Point('pA'), None),
        ]
    )
    def test_invalid_attachments_not_point(attachments: Sequence[Any]) -> None:
        with pytest.raises(TypeError):
            _ = LinearPathway(*attachments)  # type: ignore

    @pytest.fixture(autouse=True)
    def _linear_pathway_fixture(self) -> None:
        self.N = ReferenceFrame('N')
        self.pA = Point('pA')
        self.pB = Point('pB')
        self.pathway = LinearPathway(self.pA, self.pB)
        self.q1 = dynamicsymbols('q1')
        self.q2 = dynamicsymbols('q2')
        self.q3 = dynamicsymbols('q3')
        self.q1d = dynamicsymbols('q1', 1)
        self.q2d = dynamicsymbols('q2', 1)
        self.q3d = dynamicsymbols('q3', 1)
        self.F = Symbol('F')

    def test_properties_are_immutable(self) -> None:
        instance = LinearPathway(self.pA, self.pB)
        with pytest.raises(AttributeError):
            instance.attachments = None  # type: ignore
        with pytest.raises(TypeError):
            instance.attachments[0] = None  # type: ignore
        with pytest.raises(TypeError):
            instance.attachments[1] = None  # type: ignore

    def test_repr(self) -> None:
        pathway = LinearPathway(self.pA, self.pB)
        expected = 'LinearPathway(pA, pB)'
        assert repr(pathway) == expected

    def test_static_pathway_length(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.N.x)
        assert self.pathway.length == 2

    def test_static_pathway_extension_velocity(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.N.x)
        assert self.pathway.extension_velocity == 0

    def test_static_pathway_compute_loads(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.N.x)
        expected = [
            (self.pA, - self.F * self.N.x),
            (self.pB, self.F * self.N.x),
        ]
        assert self.pathway.compute_loads(self.F) == expected

    def test_2D_pathway_length(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.q1 * self.N.x)
        expected = 2 * sqrt(self.q1**2)
        assert self.pathway.length == expected

    def test_2D_pathway_extension_velocity(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.q1 * self.N.x)
        expected = 2 * self.q1 * self.q1d / sqrt(self.q1**2)
        assert self.pathway.extension_velocity == expected

    def test_2D_pathway_compute_loads(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.q1 * self.N.x)
        expected = [
            (self.pA, - self.F * (self.q1 / sqrt(self.q1**2)) * self.N.x),
            (self.pB, self.F * (self.q1 / sqrt(self.q1**2)) * self.N.x),
        ]
        assert self.pathway.compute_loads(self.F) == expected

    def test_3D_pathway_length(self) -> None:
        self.pB.set_pos(
            self.pA,
            self.q1*self.N.x - self.q2*self.N.y + 2*self.q3*self.N.z,
        )
        expected = sqrt(self.q1**2 + self.q2**2 + 4*self.q3**2)
        assert simplify(self.pathway.length - expected) == 0

    def test_3D_pathway_extension_velocity(self) -> None:
        self.pB.set_pos(
            self.pA,
            self.q1*self.N.x - self.q2*self.N.y + 2*self.q3*self.N.z,
        )
        length = sqrt(self.q1**2 + self.q2**2 + 4*self.q3**2)
        expected = (
            self.q1 * self.q1d / length
            + self.q2 * self.q2d / length
            + 4 * self.q3 * self.q3d / length
        )
        assert simplify(self.pathway.extension_velocity - expected) == 0

    def test_3D_pathway_compute_loads(self) -> None:
        self.pB.set_pos(
            self.pA,
            self.q1*self.N.x - self.q2*self.N.y + 2*self.q3*self.N.z,
        )
        length = sqrt(self.q1**2 + self.q2**2 + 4*self.q3**2)
        pO_force = (
            - self.F * self.q1 * self.N.x / length
            + self.F * self.q2 * self.N.y / length
            - 2 * self.F * self.q3 * self.N.z / length
        )
        pI_force = (
            self.F * self.q1 * self.N.x / length
            - self.F * self.q2 * self.N.y / length
            + 2 * self.F * self.q3 * self.N.z / length
        )
        expected = [
            (self.pA, pO_force),
            (self.pB, pI_force),
        ]
        assert self.pathway.compute_loads(self.F) == expected


class TestWrappingPathway:

    def test_is_pathway_base_subclass(self) -> None:
        assert issubclass(WrappingPathway, PathwayBase)

    @pytest.fixture(autouse=True)
    def _wrapping_pathway_fixture(self) -> None:
        self.pA = Point('pA')
        self.pB = Point('pB')
        self.r = Symbol('r', positive=True)
        self.pO = Point('pO')
        self.N = ReferenceFrame('N')
        self.ax = self.N.z
        self.sphere = Sphere(self.r, self.pO)
        self.cylinder = Cylinder(self.r, self.pO, self.ax)
        self.pathway = WrappingPathway(self.pA, self.pB, self.cylinder)
        self.F = Symbol('F')

    def test_valid_constructor(self) -> None:
        instance = WrappingPathway(self.pA, self.pB, self.cylinder)
        assert isinstance(instance, WrappingPathway)
        assert hasattr(instance, 'attachments')
        assert len(instance.attachments) == 2
        assert isinstance(instance.attachments[0], Point)
        assert instance.attachments[0] == self.pA
        assert isinstance(instance.attachments[1], Point)
        assert instance.attachments[1] == self.pB
        assert hasattr(instance, 'geometry')
        assert isinstance(instance.geometry, GeometryBase)
        assert instance.geometry == self.cylinder
