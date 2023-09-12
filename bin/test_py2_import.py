#!/usr/bin/env python
#
# Tests that a useful message is give in the ImportError when trying to import
# sympy from Python 2. This ensures that we don't get a Py2 SyntaxError from
# sympy/__init__.py

import sys
assert sys.version_info[:2] == (2, 7), "This test is for Python 2.7 only"

import os
thisdir = os.path.dirname(__file__)
parentdir = os.path.normpath(os.path.join(thisdir, '..'))

# Append the SymPy root directory to path
sys.path.append(parentdir)

try:
    import sympy
except ImportError as exc:
    message = str(exc)
    # e.g. "Python version 3.5 or above is required for SymPy."
    assert message.startswith("Python version")
    assert message.endswith(" or above is required for SymPy.")
else:
    raise AssertionError("import sympy should give ImportError on Python 2.7")
