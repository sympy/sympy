# This test file tests the SymPy function interface, that people use to create
# their own new functions. It should be as easy as possible.

from sympy import Function, sympify, sin, cos
from sympy.abc import x

def test_function_series1():
    """Create our new "sin" function."""

    class my_function(Function):
        nofargs = 1

        def fdiff(self, argindex = 1):
            return cos(self[0])

        @classmethod
        def canonize(cls, arg):
            arg = sympify(arg)
            if arg == 0:
                return sympify(0)

    #Test that the taylor series is correct
    assert my_function(x).series(x, 10) == sin(x).series(x, 10)

def test_function_series2():
    """Create our new "cos" function."""

    class my_function2(Function):
        nofargs = 1

        def fdiff(self, argindex = 1):
            return -sin(self[0])

        @classmethod
        def canonize(cls, arg):
            arg = sympify(arg)
            if arg == 0:
                return sympify(1)

    #Test that the taylor series is correct
    assert my_function2(x).series(x, 10) == cos(x).series(x, 10)
