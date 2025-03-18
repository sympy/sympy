#!/usr/bin/env python
"""
Run tests involving pytorch

These are separate from the other optional dependency tests because pytorch
may have specific version requirements.
"""

TEST_LIST = DOCTEST_LIST = [
    'sympy/printing/pytorch.py',
    'sympy/printing/tests/test_torch.py',
    'sympy/utilities/lambdify.py',
    'sympy/utilities/tests/test_lambdify.py',
]


if __name__ == "__main__":

    import sys

    # Add the local SymPy to sys.path (needed for CI)
    from get_sympy import path_hack
    path_hack()
    import sympy

    # Note: The doctests are not tested here but there are many failures when
    # running them with symengine.
    args = TEST_LIST
    test_exit_code = sympy.test(*args, verbose=True)
    if test_exit_code != 0:
        sys.exit(test_exit_code)
    doctest_exit_code = sympy.doctest(*DOCTEST_LIST)
    exit_code = 0 if doctest_exit_code is True else 1
    sys.exit(exit_code)
