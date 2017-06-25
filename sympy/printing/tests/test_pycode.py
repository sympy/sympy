from sympy.core import symbols
from ..pycode import PyCodePrinter

x, y, z = symbols('x y z')

def test_PyCodePrinter():
    prntr = PyCodePrinter()
    assert prntr.doprint(x**y) == 'x**y'
