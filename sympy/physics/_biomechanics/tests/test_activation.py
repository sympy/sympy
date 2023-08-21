"""Tests for the ``sympy.physics._biomechanics.activation.py`` module."""

import pytest

from sympy.matrices import Matrix
from sympy.matrices.dense import zeros
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics._biomechanics import (
    ActivationBase,
    ZerothOrderActivation,
)
from sympy.physics._biomechanics._mixin import _NamedMixin


class TestZerothOrderActivation:

    @staticmethod
    def test_class():
        assert issubclass(ZerothOrderActivation, ActivationBase)
        assert issubclass(ZerothOrderActivation, _NamedMixin)
        assert ZerothOrderActivation.__name__ == 'ZerothOrderActivation'

    @pytest.fixture(autouse=True)
    def _zeroth_order_activation_fixture(self):
        self.name = 'name'
        self.e = dynamicsymbols('e_name')
        self.instance = ZerothOrderActivation(self.name)

    def test_instance(self):
        instance = ZerothOrderActivation(self.name)
        assert isinstance(instance, ZerothOrderActivation)

    def test_with_default_constants(self):
        instance = ZerothOrderActivation.with_default_constants(self.name)
        assert isinstance(instance, ZerothOrderActivation)
        assert instance == ZerothOrderActivation(self.name)

    def test_name(self):
        assert hasattr(self.instance, 'name')
        assert self.instance.name == self.name

    def test_order(self):
        assert hasattr(self.instance, 'order')
        assert self.instance.order == 0

    def test_excitation(self):
        assert hasattr(self.instance, 'e')
        assert hasattr(self.instance, 'excitation')
        e_expected = dynamicsymbols('e_name')
        assert self.instance.e == e_expected
        assert self.instance.excitation == e_expected
        assert self.instance.e is self.instance.excitation

    def test_activation(self):
        assert hasattr(self.instance, 'a')
        assert hasattr(self.instance, 'activation')
        a_expected = dynamicsymbols('e_name')
        assert self.instance.a == a_expected
        assert self.instance.activation == a_expected
        assert self.instance.a is self.instance.activation is self.instance.e

    def test_state_vars(self):
        assert hasattr(self.instance, 'x')
        assert hasattr(self.instance, 'state_vars')
        assert self.instance.x == self.instance.state_vars
        x_expected = zeros(0, 1)
        assert self.instance.x == x_expected
        assert self.instance.state_vars == x_expected
        assert isinstance(self.instance.x, Matrix)
        assert isinstance(self.instance.state_vars, Matrix)
        assert self.instance.x.shape == (0, 1)
        assert self.instance.state_vars.shape == (0, 1)

    def test_input_vars(self):
        assert hasattr(self.instance, 'r')
        assert hasattr(self.instance, 'input_vars')
        assert self.instance.r == self.instance.input_vars
        r_expected = Matrix([self.e])
        assert self.instance.r == r_expected
        assert self.instance.input_vars == r_expected
        assert isinstance(self.instance.r, Matrix)
        assert isinstance(self.instance.input_vars, Matrix)
        assert self.instance.r.shape == (1, 1)
        assert self.instance.input_vars.shape == (1, 1)

    def test_constants(self):
        assert hasattr(self.instance, 'p')
        assert hasattr(self.instance, 'constants')
        assert self.instance.p == self.instance.constants
        p_expected = zeros(0, 1)
        assert self.instance.p == p_expected
        assert self.instance.constants == p_expected
        assert isinstance(self.instance.p, Matrix)
        assert isinstance(self.instance.constants, Matrix)
        assert self.instance.p.shape == (0, 1)
        assert self.instance.constants.shape == (0, 1)

    def test_M(self):
        assert hasattr(self.instance, 'M')
        M_expected = Matrix([])
        assert self.instance.M == M_expected
        assert isinstance(self.instance.M, Matrix)
        assert self.instance.M.shape == (0, 0)

    def test_F(self):
        assert hasattr(self.instance, 'F')
        F_expected = zeros(0, 1)
        assert self.instance.F == F_expected
        assert isinstance(self.instance.F, Matrix)
        assert self.instance.F.shape == (0, 1)

    def test_rhs(self):
        assert hasattr(self.instance, 'rhs')
        rhs_expected = zeros(0, 1)
        rhs = self.instance.rhs()
        assert rhs == rhs_expected
        assert isinstance(rhs, Matrix)
        assert rhs.shape == (0, 1)

    def test_repr(self):
        expected = 'ZerothOrderActivation(\'name\')'
        assert repr(self.instance) == expected
