"""Tests for the ``sympy.physics.mechanics._pathway.py`` module."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Sequence

import pytest

from sympy.core.backend import (
    Rational,
    Symbol,
    USE_SYMENGINE,
    cos,
    pi,
    sin,
    sqrt,
)
from sympy.physics.mechanics import (
    Force,
    Point,
    ReferenceFrame,
    WrappingCylinder,
    WrappingGeometryBase,
    WrappingSphere,
    dynamicsymbols,
)
from sympy.physics.mechanics._pathway import (
    LinearPathway,
    PathwayBase,
    WrappingPathway,
)
from sympy.simplify.simplify import simplify

if TYPE_CHECKING:
    from sympy.core.expr import Expr
    from sympy.core.numbers import Number
    from sympy.physics.mechanics.loads import LoadBase


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
        self.sphere = WrappingSphere(self.r, self.pO)
        self.cylinder = WrappingCylinder(self.r, self.pO, self.ax)
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
        assert isinstance(instance.geometry, WrappingGeometryBase)
        assert instance.geometry == self.cylinder

    @pytest.mark.parametrize(
        'attachments',
        [
            (Point('pA'), ),
            (Point('pA'), Point('pB'), Point('pZ')),
        ]
    )
    def test_invalid_constructor_attachments_incorrect_number(
        self,
        attachments: Sequence[Point],
    ) -> None:
        with pytest.raises(TypeError):
            _ = WrappingPathway(*attachments, self.cylinder)  # type: ignore

    @staticmethod
    @pytest.mark.parametrize(
        'attachments',
        [
            (None, Point('pB')),
            (Point('pA'), None),
        ]
    )
    def test_invalid_constructor_attachments_not_point(
        attachments: Sequence[Any],
    ) -> None:
        with pytest.raises(TypeError):
            _ = WrappingPathway(*attachments)  # type: ignore

    def test_invalid_constructor_geometry_is_not_supplied(self) -> None:
        with pytest.raises(TypeError):
            _ = WrappingPathway(self.pA, self.pB)  # type: ignore

    @pytest.mark.parametrize(
        'geometry',
        [
            Symbol('r'),
            dynamicsymbols('q'),
            ReferenceFrame('N'),
            ReferenceFrame('N').x,
        ]
    )
    def test_invalid_geometry_not_geometry(
        self,
        geometry: Any,
    ) -> None:
        with pytest.raises(TypeError):
            _ = WrappingPathway(self.pA, self.pB, geometry)

    def test_attachments_property_is_immutable(self) -> None:
        with pytest.raises(TypeError):
            self.pathway.attachments[0] = self.pB  # type: ignore
        with pytest.raises(TypeError):
            self.pathway.attachments[1] = self.pA  # type: ignore

    def test_geometry_property_is_immutable(self) -> None:
        with pytest.raises(AttributeError):
            self.pathway.geometry = None  # type: ignore

    def test_repr(self) -> None:
        expected = (
            f'WrappingPathway(pA, pB, '
            f'geometry={self.cylinder!r})'
        )
        assert repr(self.pathway) == expected

    @staticmethod
    def _expand_pos_to_vec(pos, frame):
        return sum(mag * unit for (mag, unit) in zip(pos, frame))

    @pytest.mark.parametrize(
        'pA_vec, pB_vec, expected_factor',
        [
            ((1, 0, 0), (0, 1, 0), pi / 2),
            ((0, 1, 0), (sqrt(2) / 2, -sqrt(2) / 2, 0), 3 * pi / 4),
            ((1, 0, 0), (Rational(1, 2), sqrt(3) / 2, 0), pi / 3),
        ]
    )
    def test_static_pathway_on_sphere_length(
        self,
        pA_vec: tuple[int | Number, int | Number, int | Number],
        pB_vec: tuple[int | Number, int | Number, int | Number],
        expected_factor: Expr,
    ) -> None:
        pA_vec = self._expand_pos_to_vec(pA_vec, self.N)
        pB_vec = self._expand_pos_to_vec(pB_vec, self.N)
        self.pA.set_pos(self.pO, self.r * pA_vec)
        self.pB.set_pos(self.pO, self.r * pB_vec)
        pathway = WrappingPathway(self.pA, self.pB, self.sphere)
        expected = expected_factor * self.r
        assert simplify(pathway.length - expected) == 0

    @pytest.mark.parametrize(
        'pA_vec, pB_vec, expected_factor',
        [
            ((1, 0, 0), (0, 1, 0), Rational(1, 2) * pi),
            ((1, 0, 0), (-1, 0, 0), pi),
            ((-1, 0, 0), (1, 0, 0), pi),
            ((0, 1, 0), (sqrt(2) / 2, -sqrt(2) / 2, 0), 5 * pi / 4),
            ((1, 0, 0), (Rational(1, 2), sqrt(3) / 2, 0), pi / 3),
            (
                (0, 1, 0),
                (sqrt(2) * Rational(1, 2), -sqrt(2)* Rational(1, 2), 1),
                sqrt(1 + (Rational(5, 4) * pi)**2),
            ),
            (
                (1, 0, 0),
                (Rational(1, 2), sqrt(3) * Rational(1, 2), 1),
                sqrt(1 + (Rational(1, 3) * pi)**2),
            ),
        ]
    )
    def test_static_pathway_on_cylinder_length(
        self,
        pA_vec: tuple[int | Number, int | Number, int | Number],
        pB_vec: tuple[int | Number, int | Number, int | Number],
        expected_factor: Expr,
    ) -> None:
        pA_vec = self._expand_pos_to_vec(pA_vec, self.N)
        pB_vec = self._expand_pos_to_vec(pB_vec, self.N)
        self.pA.set_pos(self.pO, self.r * pA_vec)
        self.pB.set_pos(self.pO, self.r * pB_vec)
        pathway = WrappingPathway(self.pA, self.pB, self.cylinder)
        expected = expected_factor * sqrt(self.r**2)
        assert simplify(pathway.length - expected) == 0

    @pytest.mark.parametrize(
        'pA_vec, pB_vec',
        [
            ((1, 0, 0), (0, 1, 0)),
            ((0, 1, 0), (sqrt(2) * Rational(1, 2), -sqrt(2)* Rational(1, 2), 0)),
            ((1, 0, 0), (Rational(1, 2), sqrt(3) * Rational(1, 2), 0)),
        ]
    )
    def test_static_pathway_on_sphere_extension_velocity(
        self,
        pA_vec: tuple[int | Number, int | Number, int | Number],
        pB_vec: tuple[int | Number, int | Number, int | Number],
    ) -> None:
        pA_vec = self._expand_pos_to_vec(pA_vec, self.N)
        pB_vec = self._expand_pos_to_vec(pB_vec, self.N)
        self.pA.set_pos(self.pO, self.r * pA_vec)
        self.pB.set_pos(self.pO, self.r * pB_vec)
        pathway = WrappingPathway(self.pA, self.pB, self.sphere)
        assert pathway.extension_velocity == 0

    @pytest.mark.parametrize(
        'pA_vec, pB_vec',
        [
            ((1, 0, 0), (0, 1, 0)),
            ((1, 0, 0), (-1, 0, 0)),
            ((-1, 0, 0), (1, 0, 0)),
            ((0, 1, 0), (sqrt(2) / 2, -sqrt(2) / 2, 0)),
            ((1, 0, 0), (Rational(1, 2), sqrt(3) / 2, 0)),
            ((0, 1, 0), (sqrt(2) * Rational(1, 2), -sqrt(2) / 2, 1)),
            ((1, 0, 0), (Rational(1, 2), sqrt(3) / 2, 1)),
        ]
    )
    def test_static_pathway_on_cylinder_extension_velocity(
        self,
        pA_vec: tuple[int | Number, int | Number, int | Number],
        pB_vec: tuple[int | Number, int | Number, int | Number],
    ) -> None:
        pA_vec = self._expand_pos_to_vec(pA_vec, self.N)
        pB_vec = self._expand_pos_to_vec(pB_vec, self.N)
        self.pA.set_pos(self.pO, self.r * pA_vec)
        self.pB.set_pos(self.pO, self.r * pB_vec)
        pathway = WrappingPathway(self.pA, self.pB, self.cylinder)
        assert pathway.extension_velocity == 0

    @pytest.mark.parametrize(
        'pA_vec, pB_vec, pA_vec_expected, pB_vec_expected, pO_vec_expected',
        (
            ((1, 0, 0), (0, 1, 0), (0, 1, 0), (1, 0, 0), (-1, -1, 0)),
            (
                (0, 1, 0),
                (sqrt(2) / 2, -sqrt(2) / 2, 0),
                (1, 0, 0),
                (sqrt(2) / 2, sqrt(2) / 2, 0),
                (-1 - sqrt(2) / 2, -sqrt(2) / 2, 0)
            ),
            (
                (1, 0, 0),
                (Rational(1, 2), sqrt(3) / 2, 0),
                (0, 1, 0),
                (sqrt(3) / 2, -Rational(1, 2), 0),
                (-sqrt(3) / 2, Rational(1, 2) - 1, 0),
            ),
        )
    )
    def test_static_pathway_on_sphere_compute_loads(
        self,
        pA_vec: tuple[int | Number, int | Number, int | Number],
        pB_vec: tuple[int | Number, int | Number, int | Number],
        pA_vec_expected: tuple[int | Number, int | Number, int | Number],
        pB_vec_expected: tuple[int | Number, int | Number, int | Number],
        pO_vec_expected: tuple[int | Number, int | Number, int | Number],
    ) -> None:
        pA_vec = self._expand_pos_to_vec(pA_vec, self.N)
        pB_vec = self._expand_pos_to_vec(pB_vec, self.N)
        self.pA.set_pos(self.pO, self.r * pA_vec)
        self.pB.set_pos(self.pO, self.r * pB_vec)
        pathway = WrappingPathway(self.pA, self.pB, self.sphere)

        pA_vec_expected = sum(
            mag * unit for (mag, unit) in zip(pA_vec_expected, self.N)
        )
        pB_vec_expected = sum(
            mag * unit for (mag, unit) in zip(pB_vec_expected, self.N)
        )
        pO_vec_expected = sum(
            mag * unit for (mag, unit) in zip(pO_vec_expected, self.N)
        )
        expected = [
            Force(self.pA, self.F * (self.r**3 / sqrt(self.r**6)) * pA_vec_expected),
            Force(self.pB, self.F * (self.r**3 / sqrt(self.r**6)) * pB_vec_expected),
            Force(self.pO, self.F * (self.r**3 / sqrt(self.r**6)) * pO_vec_expected),
        ]
        assert pathway.compute_loads(self.F) == expected

    @pytest.mark.skipif(USE_SYMENGINE, reason='SymEngine does not simplify')
    @pytest.mark.parametrize(
        'pA_vec, pB_vec, pA_vec_expected, pB_vec_expected, pO_vec_expected',
        (
            ((1, 0, 0), (0, 1, 0), (0, 1, 0), (1, 0, 0), (-1, -1, 0)),
            ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, 1, 0), (0, -2, 0)),
            ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, -1, 0), (0, 2, 0)),
            (
                (0, 1, 0),
                (sqrt(2) / 2, -sqrt(2) / 2, 0),
                (-1, 0, 0),
                (-sqrt(2) / 2, -sqrt(2) / 2, 0),
                (1 + sqrt(2) / 2, sqrt(2) / 2, 0)
            ),
            (
                (1, 0, 0),
                (Rational(1, 2), sqrt(3) / 2, 0),
                (0, 1, 0),
                (sqrt(3) / 2, -Rational(1, 2), 0),
                (-sqrt(3) / 2, Rational(1, 2) - 1, 0),
            ),
            (
                (1, 0, 0),
                (sqrt(2) / 2, sqrt(2) / 2, 0),
                (0, 1, 0),
                (sqrt(2) / 2, -sqrt(2) / 2, 0),
                (-sqrt(2) / 2, sqrt(2) / 2 - 1, 0),
            ),
            ((0, 1, 0), (0, 1, 1), (0, 0, 1), (0, 0, -1), (0, 0, 0)),
            (
                (0, 1, 0),
                (sqrt(2) / 2, -sqrt(2) / 2, 1),
                (-5 * pi / sqrt(16 + 25*pi**2), 0, 4 / sqrt(16 + 25*pi**2)),
                (
                    -5 * sqrt(2) * pi / (2 * sqrt(16 + 25*pi**2)),
                    -5 * sqrt(2) * pi / (2 * sqrt(16 + 25*pi**2)),
                    -4 / sqrt(16 + 25*pi**2),
                ),
                (
                    5 * (sqrt(2) + 2) * pi / (2 * sqrt(16 + 25*pi**2)),
                    5 * sqrt(2) * pi / (2 * sqrt(16 + 25*pi**2)),
                    0,
                ),
            ),
        )
    )
    def test_static_pathway_on_cylinder_compute_loads(
        self,
        pA_vec: tuple[int | Number, int | Number, int | Number],
        pB_vec: tuple[int | Number, int | Number, int | Number],
        pA_vec_expected: tuple[int | Number, int | Number, int | Number],
        pB_vec_expected: tuple[int | Number, int | Number, int | Number],
        pO_vec_expected: tuple[int | Number, int | Number, int | Number],
    ) -> None:
        pA_vec = self._expand_pos_to_vec(pA_vec, self.N)
        pB_vec = self._expand_pos_to_vec(pB_vec, self.N)
        self.pA.set_pos(self.pO, self.r * pA_vec)
        self.pB.set_pos(self.pO, self.r * pB_vec)
        pathway = WrappingPathway(self.pA, self.pB, self.cylinder)

        pA_force_expected = self.F * self._expand_pos_to_vec(pA_vec_expected, self.N)
        pB_force_expected = self.F * self._expand_pos_to_vec(pB_vec_expected, self.N)
        pO_force_expected = self.F * self._expand_pos_to_vec(pO_vec_expected, self.N)
        expected = [
            Force(self.pA, pA_force_expected),
            Force(self.pB, pB_force_expected),
            Force(self.pO, pO_force_expected),
        ]
        assert self._simplify_loads(pathway.compute_loads(self.F)) == expected

    @pytest.mark.skipif(USE_SYMENGINE, reason='SymEngine does not simplify')
    def test_2D_pathway_on_cylinder_length(self) -> None:
        q = dynamicsymbols('q')
        pA_pos = self.r * self.N.x
        pB_pos = self.r * (cos(q) * self.N.x + sin(q) * self.N.y)
        self.pA.set_pos(self.pO, pA_pos)
        self.pB.set_pos(self.pO, pB_pos)
        expected = self.r * sqrt(q**2)
        assert simplify(self.pathway.length - expected) == 0

    @pytest.mark.skipif(USE_SYMENGINE, reason='SymEngine does not simplify')
    def test_2D_pathway_on_cylinder_extension_velocity(self) -> None:
        q = dynamicsymbols('q')
        qd = dynamicsymbols('q', 1)
        pA_pos = self.r * self.N.x
        pB_pos = self.r * (cos(q) * self.N.x + sin(q) * self.N.y)
        self.pA.set_pos(self.pO, pA_pos)
        self.pB.set_pos(self.pO, pB_pos)
        expected = self.r * (sqrt(q**2) / q) * qd
        assert simplify(self.pathway.extension_velocity - expected) == 0

    @pytest.mark.skipif(USE_SYMENGINE, reason='SymEngine does not simplify')
    def test_2D_pathway_on_cylinder_compute_loads(self) -> None:
        q = dynamicsymbols('q')
        pA_pos = self.r * self.N.x
        pB_pos = self.r * (cos(q) * self.N.x + sin(q) * self.N.y)
        self.pA.set_pos(self.pO, pA_pos)
        self.pB.set_pos(self.pO, pB_pos)

        pA_force = self.F * self.N.y
        pB_force = self.F * (sin(q) * self.N.x - cos(q) * self.N.y)
        pO_force = self.F * (-sin(q) * self.N.x + (cos(q) - 1) * self.N.y)
        expected = [
            Force(self.pA, pA_force),
            Force(self.pB, pB_force),
            Force(self.pO, pO_force),
        ]

        loads = self._simplify_loads(self.pathway.compute_loads(self.F))
        assert loads == expected

    @staticmethod
    def _simplify_loads(loads: list[LoadBase]) -> list[LoadBase]:
        return [
            load.__class__(load.location, load.vector.simplify())
            for load in loads
        ]
