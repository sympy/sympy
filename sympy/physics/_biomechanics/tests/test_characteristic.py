"""Tests for the ``sympy.physics._biomechanics.characteristic.py`` module."""

import pytest

from sympy.core.backend import Function, Symbol, exp
from sympy.physics._biomechanics.characteristic import (
    CharacteristicCurveFunction,
    fl_T_de_groote_2016,
)
from sympy.printing.latex import LatexPrinter


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

    def test_instance(self) -> None:
        fl_T = fl_T_de_groote_2016(self.l_T_tilde, *self.constants)
        assert isinstance(fl_T, fl_T_de_groote_2016)
        assert str(fl_T) == 'fl_T_de_groote_2016(l_T_tilde, c_0, c_1, c_2, c_3)'

    def test_doit(self) -> None:
        fl_T = fl_T_de_groote_2016(self.l_T_tilde, *self.constants).doit()
        assert fl_T == self.c0*exp(self.c3*(self.l_T_tilde - self.c1)) - self.c2

    def test_function_print_latex(self) -> None:
        fl_T = fl_T_de_groote_2016(self.l_T_tilde, *self.constants)
        printer = LatexPrinter()
        output = printer.doprint(fl_T)
        expected = r'\operatorname{fl}^T \left( l_{T tilde} \right)'
        assert output == expected
