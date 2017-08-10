from __future__ import (absolute_import, division, print_function)

import pytest
import sympy as sp
from sympy.codegen.ast import Assignment
from sympy.codegen.algorithms import newtons_method, newtons_method_function
from sympy.printing.ccode import ccode
from sympy.utilities._compilation import compile_link_import_strings


def test_newtons_method():
    x, dx, atol = sp.symbols('x dx atol')
    expr = sp.cos(x) - x**3
    algo = newtons_method(expr, x, atol, dx)
    assert algo.has(Assignment(dx, -expr/expr.diff(x)))


def test_newtons_method_function():
    x = sp.Symbol('x')
    expr = sp.cos(x) - x**3
    func = newtons_method_function(expr, x)
    compile_kw = dict(std='c99')
    mod = compile_link_import_strings([
        ('newton.c', ('#include <math.h>\n'
                      '#include <stdio.h>\n') + ccode(func)),
        ('_newton.pyx', ("cdef extern double newton(double)\n"
                         "def py_newton(x):\n"
                         "    return newton(x)\n"))
    ], compile_kwargs=compile_kw)
    assert abs(mod.py_newton(0.5) - 0.865474033102) < 1e-12


def test_newtons_method_function_parameters():
    args = x, A, k, p = sp.symbols('x A k p')
    expr = A*sp.cos(k*x) - p*x**3
    with pytest.raises(ValueError):
        newtons_method_function(expr, x)
    func = newtons_method_function(expr, x, args, debug=True)
    compile_kw = dict(std='c99')
    mod = compile_link_import_strings([
        ('newton.c', ('#include <math.h>\n'
                      '#include <stdio.h>\n') + ccode(func)),
        ('_newton.pyx', ("cdef extern double newton(double, double, double, double)\n"
                         "def py_newton(x, A=1, k=1, p=1):\n"
                         "    return newton(x, A, k, p)\n"))
    ], compile_kwargs=compile_kw)
    assert abs(mod.py_newton(0.5) - 0.865474033102) < 1e-12
