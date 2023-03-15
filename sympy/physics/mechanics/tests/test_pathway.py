"""Tests for the ``sympy.physics.mechanics._pathway.py`` module."""

import pytest

from sympy.core.symbol import Symbol
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.mechanics._pathway import LinearPathway


class TestLinearPathway:

    @pytest.fixture(autouse=True)
    def _linear_pathway_fixture(self) -> None:
        self.N = ReferenceFrame('N')
        self.pO = Point('pO')
        self.pI = Point('pI')
        self.pathway = LinearPathway(origin=self.pO, insertion=self.pI)
        self.F = Symbol('F')

    def test_static_pathway_length(self) -> None:
        self.pI.set_pos(self.pO, 2 * self.N.x)
        assert self.pathway.length == 2

    def test_static_pathway_shortening_velocity(self) -> None:
        self.pI.set_pos(self.pO, 2 * self.N.x)
        assert self.pathway.shortening_velocity == 0

    def test_static_pathway_forces(self) -> None:
        self.pI.set_pos(self.pO, 2 * self.N.x)
        expected = [
            (self.pO, self.F * self.N.x),
            (self.pI, - self.F * self.N.x),
        ]
        assert self.pathway.forces(self.F) == expected
