from __future__ import annotations
from sympy.abc import x
from sympy.external import import_module
from sympy.functions.elementary.trigonometric import sin
from sympy.printing.gtk import print_gtk
from sympy.testing.pytest import raises, skip, warns_deprecated_sympy

lxml = import_module('lxml')


def test_1():
    if not lxml:
        skip("lxml not installed.")

    with warns_deprecated_sympy():
        print_gtk(x**2, start_viewer=False)
    with warns_deprecated_sympy():
        print_gtk(x**2 + sin(x)/4, start_viewer=False)


def test_settings():
    with warns_deprecated_sympy():
        raises(TypeError, lambda: print_gtk(x, method="garbage"))
