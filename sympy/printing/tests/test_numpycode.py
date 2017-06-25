from sympy.core import symbols
from ..numpycode import NumPyCodePrinter

x, y, z = symbols('x y z')

def test_NumPyCodePrinter():
    prntr = NumPyCodePrinter()
    assert prntr.doprint(x**y) == 'x**y'
