"""Tests for the ``sympy.physics.mechanics._pathway.py`` module."""

from __future__ import annotations

from typing import Any, Sequence

import pytest

from sympy.core.backend import Symbol, sqrt
from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.mechanics._pathway import LinearPathway


class TestLinearPathway:

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

    def test_static_pathway_length(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.N.x)
        assert self.pathway.length == 2

    def test_static_pathway_shortening_velocity(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.N.x)
        assert self.pathway.shortening_velocity == 0

    def test_static_pathway_compute_loads(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.N.x)
        expected = [
            (self.pA, self.F * self.N.x),
            (self.pB, - self.F * self.N.x),
        ]
        assert self.pathway.compute_loads(self.F) == expected

    def test_2D_pathway_length(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.q1 * self.N.x)
        expected = 2 * sqrt(self.q1**2)
        assert self.pathway.length == expected

    def test_2D_pathway_shortening_velocity(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.q1 * self.N.x)
        expected = -2 * self.q1 * self.q1d / sqrt(self.q1**2)
        assert self.pathway.shortening_velocity == expected

    def test_2D_pathway_compute_loads(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.q1 * self.N.x)
        expected = [
            (self.pA, self.F * (self.q1 / sqrt(self.q1**2)) * self.N.x),
            (self.pB, - self.F * (self.q1 / sqrt(self.q1**2)) * self.N.x),
        ]
        assert self.pathway.compute_loads(self.F) == expected

    def test_3D_pathway_length(self) -> None:
        self.pB.set_pos(
            self.pA,
            self.q1*self.N.x - self.q2*self.N.y + 2*self.q3*self.N.z,
        )
        expected = sqrt(self.q1**2 + self.q2**2 + 4*self.q3**2)
        assert self.pathway.length == expected

    def test_3D_pathway_shortening_velocity(self) -> None:
        self.pB.set_pos(
            self.pA,
            self.q1*self.N.x - self.q2*self.N.y + 2*self.q3*self.N.z,
        )
        length = sqrt(self.q1**2 + self.q2**2 + 4*self.q3**2)
        expected = (
            - self.q1 * self.q1d / length
            - self.q2 * self.q2d / length
            - 4 * self.q3 * self.q3d / length
        )
        assert self.pathway.shortening_velocity == expected

    def test_3D_pathway_compute_loads(self) -> None:
        self.pB.set_pos(
            self.pA,
            self.q1*self.N.x - self.q2*self.N.y + 2*self.q3*self.N.z,
        )
        length = sqrt(self.q1**2 + self.q2**2 + 4*self.q3**2)
        pO_force = (
            self.F * self.q1 * self.N.x / length
            - self.F * self.q2 * self.N.y / length
            + 2 * self.F * self.q3 * self.N.z / length
        )
        pI_force = (
            - self.F * self.q1 * self.N.x / length
            + self.F * self.q2 * self.N.y / length
            - 2 * self.F * self.q3 * self.N.z / length
        )
        expected = [
            (self.pA, pO_force),
            (self.pB, pI_force),
        ]
        assert self.pathway.compute_loads(self.F) == expected
