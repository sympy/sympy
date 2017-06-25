from sympy.core import symbols
from ..scipycode import SciPyCodePrinter

x, y, z = symbols('x y z')

def test_SciPyCodePrinter():
    prntr = SciPyCodePrinter()
    assert prntr.doprint(x**y) == 'x**y'
