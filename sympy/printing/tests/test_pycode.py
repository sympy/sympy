# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)
from sympy.core import symbols
from sympy.printing.pycode import PythonCodePrinter

x, y, z = symbols('x y z')

def test_PythonCodePrinter():
    prntr = PythonCodePrinter()
    assert prntr.doprint(x**y) == 'x**y'
