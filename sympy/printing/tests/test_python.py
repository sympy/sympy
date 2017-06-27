# -*- coding: utf-8 -*-
import warnings
from sympy import limit, oo, Symbol
from sympy.printing.python import python
from sympy.utilities.pytest import raises
from sympy.utilities.exceptions import SymPyDeprecationWarning

x = Symbol

def test_python_limits():
    with warnings.catch_warnings():
        warnings.filterwarnings("error", category=SymPyDeprecationWarning)
        with raises(SymPyDeprecationWarning):
            python(limit(x, x, oo))
