"""Tests for the ``sympy.physics._biomechanics.characteristic.py`` module."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from sympy.core.expr import UnevaluatedExpr
from sympy.core.function import Function
from sympy.core.numbers import Float, Integer
from sympy.core.symbol import Symbol
from sympy.external.importtools import import_module
from sympy.functions.elementary.exponential import exp
from sympy.physics._biomechanics.characteristic import (
    CharacteristicCurveFunction,
    TendonForceLengthDeGroote2016,
    TendonForceLengthInverseDeGroote2016,
)
from sympy.printing.c import C89CodePrinter, C99CodePrinter, C11CodePrinter
from sympy.printing.cxx import (
    CXX98CodePrinter,
    CXX11CodePrinter,
    CXX17CodePrinter,
)
from sympy.printing.fortran import FCodePrinter
from sympy.printing.lambdarepr import LambdaPrinter
from sympy.printing.latex import LatexPrinter
from sympy.printing.octave import OctaveCodePrinter
from sympy.printing.numpy import (
    CuPyPrinter,
    JaxPrinter,
    NumPyPrinter,
    SciPyPrinter,
)
from sympy.printing.pycode import MpmathPrinter, PythonCodePrinter
from sympy.utilities.lambdify import lambdify

if TYPE_CHECKING:
    from sympy.printing.codeprinter import CodePrinter

jax = import_module('jax')
numpy = import_module('numpy')

if jax:
    jax.config.update('jax_enable_x64', True)


class TestTendonForceLengthDeGroote2016:

    @pytest.fixture(autouse=True)
    def _tendon_force_length_arguments_fixture(self) -> None:
        self.l_T_tilde = Symbol('l_T_tilde')
        self.c0 = Symbol('c_0')
        self.c1 = Symbol('c_1')
        self.c2 = Symbol('c_2')
        self.c3 = Symbol('c_3')
        self.constants = (self.c0, self.c1, self.c2, self.c3)

    @staticmethod
    def test_class() -> None:
        assert issubclass(TendonForceLengthDeGroote2016, Function)
        assert issubclass(TendonForceLengthDeGroote2016, CharacteristicCurveFunction)
        assert TendonForceLengthDeGroote2016.__name__ == 'TendonForceLengthDeGroote2016'

    def test_instance(self) -> None:
        fl_T = TendonForceLengthDeGroote2016(self.l_T_tilde, *self.constants)
        assert isinstance(fl_T, TendonForceLengthDeGroote2016)
        assert str(fl_T) == 'TendonForceLengthDeGroote2016(l_T_tilde, c_0, c_1, c_2, c_3)'

    def test_doit(self) -> None:
        fl_T = TendonForceLengthDeGroote2016(self.l_T_tilde, *self.constants).doit()
        assert fl_T == self.c0*exp(self.c3*(self.l_T_tilde - self.c1)) - self.c2

    def test_doit_evaluate_false(self) -> None:
        fl_T = TendonForceLengthDeGroote2016(self.l_T_tilde, *self.constants).doit(evaluate=False)
        assert fl_T == self.c0*exp(self.c3*UnevaluatedExpr(self.l_T_tilde - self.c1)) - self.c2

    def test_with_default_constants(self) -> None:
        constants = (
            Float('0.2'),
            Float('0.995'),
            Float('0.25'),
            Float('33.93669377311689'),
        )
        fl_T_manual = TendonForceLengthDeGroote2016(self.l_T_tilde, *constants)
        fl_T_constants = TendonForceLengthDeGroote2016.with_default_constants(self.l_T_tilde)
        assert fl_T_manual == fl_T_constants

    def test_differentiate_wrt_l_T_tilde(self) -> None:
        fl_T = TendonForceLengthDeGroote2016(self.l_T_tilde, *self.constants)
        expected = self.c0 * self.c3 * exp(self.c3 * UnevaluatedExpr(-self.c1 + self.l_T_tilde))
        assert fl_T.diff(self.l_T_tilde) == expected

    def test_differentiate_wrt_c0(self) -> None:
        fl_T = TendonForceLengthDeGroote2016(self.l_T_tilde, *self.constants)
        expected = exp(self.c3 * UnevaluatedExpr(-self.c1 + self.l_T_tilde))
        assert fl_T.diff(self.c0) == expected

    def test_differentiate_wrt_c1(self) -> None:
        fl_T = TendonForceLengthDeGroote2016(self.l_T_tilde, *self.constants)
        expected = -self.c0 * self.c3 * exp(self.c3 * UnevaluatedExpr(self.l_T_tilde - self.c1))
        assert fl_T.diff(self.c1) == expected

    def test_differentiate_wrt_c2(self) -> None:
        fl_T = TendonForceLengthDeGroote2016(self.l_T_tilde, *self.constants)
        expected = Integer(-1)
        assert fl_T.diff(self.c2) == expected

    def test_differentiate_wrt_c3(self) -> None:
        fl_T = TendonForceLengthDeGroote2016(self.l_T_tilde, *self.constants)
        expected = self.c0 * (self.l_T_tilde - self.c1) * exp(self.c3 * UnevaluatedExpr(self.l_T_tilde - self.c1))
        assert fl_T.diff(self.c3) == expected

    def test_inverse(self) -> None:
        fl_T = TendonForceLengthDeGroote2016(self.l_T_tilde, *self.constants)
        assert fl_T.inverse() is TendonForceLengthInverseDeGroote2016

    def test_function_print_latex(self) -> None:
        fl_T = TendonForceLengthDeGroote2016(self.l_T_tilde, *self.constants)
        expected = r'\operatorname{fl}^T \left( l_{T tilde} \right)'
        assert LatexPrinter().doprint(fl_T) == expected

    def test_expression_print_latex(self) -> None:
        fl_T = TendonForceLengthDeGroote2016(self.l_T_tilde, *self.constants)
        expected = r'c_{0} e^{c_{3} \left(- c_{1} + l_{T tilde}\right)} - c_{2}'
        assert LatexPrinter().doprint(fl_T.doit()) == expected

    @pytest.mark.parametrize(
        'code_printer, expected',
        [
            (
                C89CodePrinter,
                '-0.25 + 0.20000000000000001*exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                C99CodePrinter,
                '-0.25 + 0.20000000000000001*exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                C11CodePrinter,
                '-0.25 + 0.20000000000000001*exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                CXX98CodePrinter,
                '-0.25 + 0.20000000000000001*exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                CXX11CodePrinter,
                '-0.25 + 0.20000000000000001*std::exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                CXX17CodePrinter,
                '-0.25 + 0.20000000000000001*std::exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                FCodePrinter,
                '      -0.25d0 + 0.2d0*exp(33.93669377311689d0*(l_T_tilde - 0.995d0))',
            ),
            (
                OctaveCodePrinter,
                '-0.25 + 0.2*exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                PythonCodePrinter,
                '-0.25 + 0.2*math.exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                NumPyPrinter,
                '-0.25 + 0.2*numpy.exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                SciPyPrinter,
                '-0.25 + 0.2*numpy.exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                CuPyPrinter,
                '-0.25 + 0.2*cupy.exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                JaxPrinter,
                '-0.25 + 0.2*jax.numpy.exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
            (
                MpmathPrinter,
                'mpmath.mpf((1, 1, -2, 1)) + mpmath.mpf((0, 3602879701896397, -54, 52))'
                '*mpmath.exp(mpmath.mpf((0, 9552330089424741, -48, 54))*(l_T_tilde + '
                'mpmath.mpf((1, 8962163258467287, -53, 53))))',
            ),
            (
                LambdaPrinter,
                '-0.25 + 0.2*math.exp(33.93669377311689*(l_T_tilde - 0.995))',
            ),
        ]
    )
    def test_print_code(self, code_printer: type[CodePrinter], expected: str) -> None:
        fl_T = TendonForceLengthDeGroote2016.with_default_constants(self.l_T_tilde)
        assert code_printer().doprint(fl_T) == expected

    def test_derivative_print_code(self) -> None:
        fl_T = TendonForceLengthDeGroote2016.with_default_constants(self.l_T_tilde)
        dfl_T_dl_T_tilde = fl_T.diff(self.l_T_tilde)
        expected = '6.787338754623378*math.exp(33.93669377311689*(l_T_tilde - 0.995))'
        assert PythonCodePrinter().doprint(dfl_T_dl_T_tilde) == expected

    def test_lambdify(self) -> None:
        fl_T = TendonForceLengthDeGroote2016.with_default_constants(self.l_T_tilde)
        fl_T_callable = lambdify(self.l_T_tilde, fl_T)
        assert fl_T_callable(1.0) == pytest.approx(-0.013014055039221595)

    @pytest.mark.skipif(numpy is None, reason='NumPy not installed')
    def test_lambdify_numpy(self) -> None:
        fl_T = TendonForceLengthDeGroote2016.with_default_constants(self.l_T_tilde)
        fl_T_callable = lambdify(self.l_T_tilde, fl_T, 'numpy')
        l_T_tilde = numpy.array([0.95, 1.0, 1.01, 1.05])
        expected = numpy.array([
            -0.2065693181344816,
            -0.0130140550392216,
            0.0827421191989246,
            1.04314889144172,
        ])
        numpy.testing.assert_allclose(fl_T_callable(l_T_tilde), expected)

    @pytest.mark.skipif(jax is None, reason='JAX not installed')
    def test_lambdify_jax(self) -> None:
        fl_T = TendonForceLengthDeGroote2016.with_default_constants(self.l_T_tilde)
        fl_T_callable = jax.jit(lambdify(self.l_T_tilde, fl_T, 'jax'))
        l_T_tilde = jax.numpy.array([0.95, 1.0, 1.01, 1.05])
        expected = jax.numpy.array([
            -0.2065693181344816,
            -0.0130140550392216,
            0.0827421191989246,
            1.04314889144172,
        ])
        numpy.testing.assert_allclose(fl_T_callable(l_T_tilde), expected)


class TestTendonForceLengthInverseDeGroote2016:

    @pytest.fixture(autouse=True)
    def _tendon_force_length_inverse_arguments_fixture(self) -> None:
        self.fl_T = Symbol('fl_T')
        self.c0 = Symbol('c_0')
        self.c1 = Symbol('c_1')
        self.c2 = Symbol('c_2')
        self.c3 = Symbol('c_3')
        self.constants = (self.c0, self.c1, self.c2, self.c3)

    @staticmethod
    def test_class() -> None:
        assert issubclass(TendonForceLengthInverseDeGroote2016, Function)
        assert issubclass(TendonForceLengthInverseDeGroote2016, CharacteristicCurveFunction)
        assert TendonForceLengthInverseDeGroote2016.__name__ == 'TendonForceLengthInverseDeGroote2016'

    def test_instance(self) -> None:
        fl_T_inv = TendonForceLengthInverseDeGroote2016(self.fl_T, *self.constants)
        assert isinstance(fl_T_inv, TendonForceLengthInverseDeGroote2016)
        assert str(fl_T_inv) == 'TendonForceLengthInverseDeGroote2016(fl_T, c_0, c_1, c_2, c_3)'
