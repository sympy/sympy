"""Tests for objects from `sympy/physics/_biomechanics/mixin.py`."""


import pytest

import sympy as sm
from sympy.physics._biomechanics.mixin import _NamedMixin


class TestNamedMixin:

    class NamedMixinClass(_NamedMixin):
        """Utility class inheriting `_NamedMixin` that can be instantiated."""

        def __init__(self, name):
            self.name = name

    @pytest.mark.parametrize('name', ['a', 'valid_name', '_valid_name'])
    def test_subclasses_get_and_set_correctly(self, name):
        """Subclasses of `_NamedMixin` set and get valid names correctly.

        Valid names are `str` with nonzero length.

        """
        instance = self.NamedMixinClass(name)
        assert instance.name == name
        assert isinstance(instance.name, str)

    def test_zero_length_name_raises_value_error(self):
        """A `ValueError` is raised is the supplied name is `''`."""
        with pytest.raises(ValueError):
            _ = self.NamedMixinClass('')

    @pytest.mark.parametrize('name', [0, sm.Symbol('name')])
    def test_non_string_name_raises_type_error(self, name):
        """A `TypeError` is raised if `name` is not a `str`."""
        with pytest.raises(TypeError):
            _ = self.NamedMixinClass(name)

    def test_name_immutable_after_initialization(self):
        """An `AttributeError` is raised if mutating `name` is attempted."""
        instance = self.NamedMixinClass('name')
        with pytest.raises(AttributeError):
            instance.name = 'new_name'
