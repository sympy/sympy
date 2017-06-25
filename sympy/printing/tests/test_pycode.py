from sympy.core import symbols
from ..pycode import PythonCodePrinter

x, y, z = symbols('x y z')

def test_PythonCodePrinter():
    prntr = PythonCodePrinter()
    assert prntr.doprint(x**y) == 'x**y'
