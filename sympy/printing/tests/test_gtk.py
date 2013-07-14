from sympy import print_gtk, sin
from sympy.external import import_module
from sympy.utilities.pytest import XFAIL, raises

if not import_module('lxml'):
    disabled = True

def test_1():
    from sympy.abc import x
    print_gtk(x**2, start_viewer=False)
    print_gtk(x**2 + sin(x)/4, start_viewer=False)


def test_settings():
    from sympy.abc import x
    raises(TypeError, lambda: print_gtk(x, method="garbage"))
