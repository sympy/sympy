"""Tests for the ``sympy.physics._biomechanics.characteristic.py`` module."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from sympy.core.backend import Function, Symbol, exp
from sympy.external import import_module
from sympy.physics._biomechanics.characteristic import (
    CharacteristicCurveFunction,
    fl_T_de_groote_2016,
)
from sympy.printing.c import C89CodePrinter, C99CodePrinter, C11CodePrinter
from sympy.printing.cxx import (
    CXX98CodePrinter,
    CXX11CodePrinter,
    CXX17CodePrinter,
)
from sympy.printing.fortran import FCodePrinter
from sympy.printing.latex import LatexPrinter
from sympy.printing.octave import OctaveCodePrinter
from sympy.printing.numpy import (
    CuPyPrinter,
    JaxPrinter,
    NumPyPrinter,
    SciPyPrinter,
)
from sympy.printing.pycode import PythonCodePrinter
from sympy.utilities.lambdify import lambdify

if TYPE_CHECKING:
    from sympy.printing.codeprinter import CodePrinter

jax = import_module('jax')
numpy = import_module('numpy')

if jax:
    jax.config.update('jax_enable_x64', True)


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

    def test_expression_print_latex(self) -> None:
        fl_T = fl_T_de_groote_2016(self.l_T_tilde, *self.constants)
        printer = LatexPrinter()
        output = printer.doprint(fl_T.doit())
        expected = r'c_{0} e^{c_{3} \left(- c_{1} + l_{T tilde}\right)} - c_{2}'
        assert output == expected

    @pytest.mark.parametrize(
        'code_printer, expected',
        [
            (
                C89CodePrinter,
                '0.20000000000000001*exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
            (
                C99CodePrinter,
                '0.20000000000000001*exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
            (
                C11CodePrinter,
                '0.20000000000000001*exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
            (
                CXX98CodePrinter,
                '0.20000000000000001*exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
            (
                CXX11CodePrinter,
                '0.20000000000000001*std::exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
            (
                CXX17CodePrinter,
                '0.20000000000000001*std::exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
            (
                FCodePrinter,
                '      0.2d0*2.1635839670136945d-15*exp(33.93669377311689d0*l_T_tilde) -\n      @ 0.25d0',
            ),
            (
                OctaveCodePrinter,
                '0.2*exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
            (
                PythonCodePrinter,
                '0.2*math.exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
            (
                NumPyPrinter,
                '0.2*numpy.exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
            (
                SciPyPrinter,
                '0.2*numpy.exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
            (
                CuPyPrinter,
                '0.2*cupy.exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
            (
                JaxPrinter,
                '0.2*jax.numpy.exp(33.93669377311689*(l_T_tilde - 0.995)) - 0.25',
            ),
        ]
    )
    def test_print_code(self, code_printer: type[CodePrinter], expected: str) -> None:
        fl_T = fl_T_de_groote_2016.with_default_constants(self.l_T_tilde)
        assert code_printer().doprint(fl_T) == expected

    def test_lambdify(self) -> None:
        fl_T = fl_T_de_groote_2016.with_default_constants(self.l_T_tilde)
        fl_T_callable = lambdify(self.l_T_tilde, fl_T)
        assert fl_T_callable(1.0) == pytest.approx(-0.013014055039221595)

    @pytest.mark.skipif(numpy is None, reason='NumPy not installed')
    def test_lambdify_numpy(self) -> None:
        fl_T = fl_T_de_groote_2016.with_default_constants(self.l_T_tilde)
        fl_T_callable = lambdify(self.l_T_tilde, fl_T, 'numpy')
        l_T_tilde = numpy.array([0.95, 1.0, 1.01, 1.05])
        expected = numpy.array([
            -0.2065693181344816,
            -0.0130140550392216,
            0.0827421191989246,
            1.04314889144172,
        ])
        numpy.testing.assert_allclose(fl_T_callable(l_T_tilde), expected)
