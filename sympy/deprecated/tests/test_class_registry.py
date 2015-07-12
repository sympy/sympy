from __future__ import division, print_function

from sympy.deprecated.class_registry import C
from sympy.utilities.exceptions import SymPyDeprecationWarning
from sympy.utilities.pytest import raises

def test_C():
    with raises(SymPyDeprecationWarning):
        C.Add
