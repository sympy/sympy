from __future__ import annotations

import os

from sympy.testing.runtests import convert_to_native_paths


def test_convert_to_native_paths():
    expected = os.path.normcase(os.path.join("sympy", "solvers"))
    paths = [
        "sympy/solvers",
        r"sympy\solvers",
        "./sympy/solvers",
        r".\sympy\solvers",
    ]
    assert set(convert_to_native_paths(paths)) == {expected}
