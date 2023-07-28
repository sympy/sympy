"""Tests for the ``sympy.physics._biomechanics.characteristic.py`` module."""

import pytest

from sympy.core.backend import Symbol
from sympy.physics._biomechanics.characteristic import (
    CharacteristicCurveFunction,
    fl_T_de_groote_2016,
)


class TestTendonForceLengthDeGroote2016:

    @pytest.fixture(autouse=True)
    def _fl_T_de_groote_2016_fixture(self) -> None:
        self.l_T_tilde = Symbol(r'l_T_tilde')
        self.c0 = Symbol(r'c_0')
        self.c1 = Symbol(r'c_1')
        self.c2 = Symbol(r'c_2')
        self.c3 = Symbol(r'c_3')
        self.constants = (self.c0, self.c1, self.c2, self.c3)

    @staticmethod
    def test_class() -> None:
        assert issubclass(fl_T_de_groote_2016, Function)
        assert issubclass(fl_T_de_groote_2016, CharacteristicCurveFunction)
        assert fl_T_de_groote_2016.__name__ == 'fl_T_de_groote_2016'
