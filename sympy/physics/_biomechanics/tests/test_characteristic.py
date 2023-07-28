"""Tests for the ``sympy.physics._biomechanics.characteristic.py`` module."""

from sympy.core.function import Function
from sympy.physics._biomechanics.characteristic import (
    CharacteristicCurveFunction,
    fl_T_de_groote_2016,
)


class TestTendonForceLengthDeGroote2016:

    @staticmethod
    def test_class() -> None:
        assert issubclass(fl_T_de_groote_2016, Function)
        assert issubclass(fl_T_de_groote_2016, CharacteristicCurveFunction)
        assert fl_T_de_groote_2016.__name__ == 'fl_T_de_groote_2016'
