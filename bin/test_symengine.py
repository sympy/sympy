#!/usr/bin/env python
"""
Run tests involving SymEngine

These are separate from the other optional dependency tests because they need
to be run with the `USE_SYMENGINE=1` environment variable set.

Run this as:

    $ USE_SYMENGINE=1 bin/test_symengine.py

"""

TEST_LIST = [
    'sympy/physics/mechanics',
    'sympy/liealgebras',
]


if __name__ == "__main__":

    import os
    import sys
    os.environ["USE_SYMENGINE"] = "1"

    # Add the local SymPy to sys.path (needed for CI)
    from get_sympy import path_hack
    path_hack()
    import sympy

    # Note: The doctests are not tested here but there are many failures when
    # running them with symengine.
    args = TEST_LIST
    exit_code = sympy.test(*args, verbose=True)
    sys.exit(exit_code)
