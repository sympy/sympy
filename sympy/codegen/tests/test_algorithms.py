from __future__ import (absolute_import, division, print_function)

import sys
import pytest
import sympy as sp

from sympy.core.compatibility import exec_
from sympy.codegen.ast import Assignment
from sympy.codegen.algorithms import newtons_method, newtons_method_function
from sympy.codegen.fnodes import bind_C
from sympy.codegen.futils import render_as_module as f_module
from sympy.codegen.pyutils import render_as_module as py_module
from sympy.external import import_module
from sympy.printing.ccode import ccode
from sympy.utilities._compilation import compile_link_import_strings
from sympy.utilities.pytest import skip, USE_PYTEST

cython = import_module('cython')
wurlitzer = import_module('wurlitzer')

def test_newtons_method():
    x, dx, atol = sp.symbols('x dx atol')
    expr = sp.cos(x) - x**3
    algo = newtons_method(expr, x, atol, dx)
    assert algo.has(Assignment(dx, -expr/expr.diff(x)))


def test_newtons_method_function__ccode():
    x = sp.Symbol('x')
    expr = sp.cos(x) - x**3
    func = newtons_method_function(expr, x)

    if not cython:
        skip("cython not installed.")

    compile_kw = dict(std='c99')
    mod = compile_link_import_strings([
        ('newton.c', ('#include <math.h>\n'
                      '#include <stdio.h>\n') + ccode(func)),
        ('_newton.pyx', ("cdef extern double newton(double)\n"
                         "def py_newton(x):\n"
                         "    return newton(x)\n"))
    ], compile_kwargs=compile_kw)
    assert abs(mod.py_newton(0.5) - 0.865474033102) < 1e-12


def test_newtons_method_function__fcode():
    x = sp.Symbol('x')
    expr = sp.cos(x) - x**3
    func = newtons_method_function(expr, x, attrs=[bind_C(name='newton')])

    if not cython:
        skip("cython not installed.")

    f_mod = f_module([func], 'mod_newton')
    mod = compile_link_import_strings([
        ('newton.f90', f_mod),
        ('_newton.pyx', ("cdef extern double newton(double*)\n"
                         "def py_newton(double x):\n"
                         "    return newton(&x)\n"))
    ])
    assert abs(mod.py_newton(0.5) - 0.865474033102) < 1e-12


def test_newtons_method_function__pycode():
    x = sp.Symbol('x')
    expr = sp.cos(x) - x**3
    func = newtons_method_function(expr, x)
    py_mod = py_module(func)
    namespace = {}
    exec_(py_mod, namespace, namespace)
    res = eval('newton(0.5)', namespace)
    assert abs(res - 0.865474033102) < 1e-12


def test_newtons_method_function__ccode_parameters():
    args = x, A, k, p = sp.symbols('x A k p')
    expr = A*sp.cos(k*x) - p*x**3
    with pytest.raises(ValueError):
        newtons_method_function(expr, x)
    use_wurlitzer = wurlitzer and sys.version_info[0] == 2  # weird threading issues with Python 3

    func = newtons_method_function(expr, x, args, debug=use_wurlitzer)

    if not cython:
        skip("cython not installed.")

    compile_kw = dict(std='c99')
    mod = compile_link_import_strings([
        ('newton.c', ('#include <math.h>\n'
                      '#include <stdio.h>\n') + ccode(func)),
        ('_newton.pyx', ("cdef extern double newton(double, double, double, double)\n"
                         "def py_newton(x, A=1, k=1, p=1):\n"
                         "    return newton(x, A, k, p)\n"))
    ], compile_kwargs=compile_kw)

    if use_wurlitzer:
        with wurlitzer.pipes() as (out, err):
            result = mod.py_newton(0.5)
    else:
        result = mod.py_newton(0.5)

    assert abs(result - 0.865474033102) < 1e-12

    if not use_wurlitzer:
        skip("C-level output only tested on Python 2 when package wurlitzer is available.")

    assert err.read() == ''
    assert out.read() == """\
x=         0.5 d_x=     0.61214
x=      1.1121 d_x=    -0.20247
x=     0.90967 d_x=   -0.042409
x=     0.86726 d_x=  -0.0017867
x=     0.86548 d_x= -3.1022e-06
x=     0.86547 d_x= -9.3421e-12
x=     0.86547 d_x=  3.6902e-17
"""
