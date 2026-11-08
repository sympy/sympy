from __future__ import annotations

import os
from doctest import ELLIPSIS, NORMALIZE_WHITESPACE

from sympy.testing.runtests import SymPyOutputChecker, convert_to_native_paths


def test_convert_to_native_paths():
    expected = os.path.normcase(os.path.join("sympy", "solvers"))
    paths = [
        "sympy/solvers",
        r"sympy\solvers",
        "./sympy/solvers",
        r".\sympy\solvers",
    ]
    assert set(convert_to_native_paths(paths)) == {expected}


def test_sympy_output_checker_ellipsis():
    checker = SymPyOutputChecker()

    # An ellipsis can hide several floating-point coefficients.
    got = 'Force(0.814, 1.06, 0.162, 0.0633, 0.433, -8.149, -0.374, 0.886)'
    want = 'Force(0.814, ..., -8.149, -0.374, 0.886)'
    assert checker.check_output(want, got, ELLIPSIS)

    # Without ELLIPSIS the output must not match, and an incorrect
    # visible number must not be hidden by the wildcard.
    assert not checker.check_output(want, got, 0)
    assert not checker.check_output(
        'Force(0.814, ..., -8.149, -0.374, 0.999)', got, ELLIPSIS)

    # The standard doctest matcher also normalizes whitespace.
    assert checker.check_output(
        'Force(0.814, ...,\n -8.149, -0.374, 0.886)',
        got, ELLIPSIS | NORMALIZE_WHITESPACE)

    # Approximate float matching without an ellipsis is unchanged.
    assert checker.check_output('0.60000', '0.6000001', 0)
