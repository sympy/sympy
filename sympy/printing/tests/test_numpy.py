from sympy import Piecewise
from sympy.abc import x
from sympy.printing.lambdarepr import NumPyPrinter

def test_numpy_piecewise_regression():
    """
    NumPyPrinter needs to print Piecewise()'s choicelist as a list to avoid
    breaking compatibility with numpy 1.8. This is not necessary in numpy 1.9+.
    See gh-9747 and gh-9749 for details.
    """
    p = Piecewise((1, x < 0), (0, True))
    assert NumPyPrinter().doprint(p) == 'select([x < 0,True], [1,0], default=nan)'
