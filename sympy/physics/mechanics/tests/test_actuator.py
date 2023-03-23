"""Tests for the ``sympy.physics.mechanics._actuator.py`` module."""

from __future__ import annotations

import pytest

from sympy.core.backend import Symbol
from sympy.core.expr import Expr
from sympy.physics.mechanics import (
    Point,
    ReferenceFrame,
    dynamicsymbols,
)
from sympy.physics.mechanics._actuator import ForceActuator
from sympy.physics.mechanics._pathway import LinearPathway, PathwayBase


class TestForceActuator:

    @pytest.fixture(autouse=True)
    def _linear_pathway_fixture(self) -> None:
        self.force = Symbol('F')
        self.pA = Point('pA')
        self.pB = Point('pB')
        self.pathway = LinearPathway(self.pA, self.pB)
        self.q1 = dynamicsymbols('q1')
        self.q2 = dynamicsymbols('q2')
        self.q3 = dynamicsymbols('q3')
        self.q1d = dynamicsymbols('q1', 1)
        self.q2d = dynamicsymbols('q2', 1)
        self.q3d = dynamicsymbols('q3', 1)
        self.N = ReferenceFrame('N')

    @pytest.mark.parametrize(
        'force',
        [
            Symbol('F'),
            dynamicsymbols('F'),
            Symbol('F')**2 + Symbol('F'),
        ]
    )
    def test_valid_constructor_force(self, force: Expr) -> None:
        instance = ForceActuator(force, self.pathway)
        assert isinstance(instance, ForceActuator)
        assert hasattr(instance, 'force')
        assert isinstance(instance.force, Expr)

    def test_invalid_constructor_force_not_expr(self) -> None:
        with pytest.raises(TypeError):
            _ = ForceActuator(None, self.pathway)  # type: ignore

    @pytest.mark.parametrize(
        'pathway',
        [
            LinearPathway(Point('pA'), Point('pB')),
        ]
    )
    def test_valid_constructor_pathway(self, pathway: PathwayBase) -> None:
        instance = ForceActuator(self.force, pathway)
        assert isinstance(instance, ForceActuator)
        assert hasattr(instance, 'pathway')
        assert isinstance(instance.pathway, LinearPathway)

    def test_invalid_constructor_pathway_not_pathway_base(self) -> None:
        with pytest.raises(TypeError):
            _ = ForceActuator(self.force, None)  # type: ignore

    def test_compute_loads_static_pathway(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.N.x)
        actuator = ForceActuator(self.force, self.pathway)
        expected = [
            (self.pA, self.force * self.N.x),
            (self.pB, - self.force * self.N.x),
        ]
        assert actuator.compute_loads() == expected
