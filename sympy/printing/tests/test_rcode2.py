from sympy.core import pi, oo, symbols, Function, Rational, Integer, GoldenRatio, EulerGamma, Catalan, Lambda, Dummy, Eq
from sympy.functions import Piecewise, sin, cos, Abs, exp, ceiling, sqrt, gamma
from sympy.utilities.pytest import raises
from sympy.printing.rcode import RCodePrinter
from sympy.utilities.lambdify import implemented_function
from sympy.tensor import IndexedBase, Idx

# import test
from sympy import rcode

x, y, z = symbols('x,y,z')
g = Function('g')



#mm
def test_rcode_constants_other():
    assert rcode(GoldenRatio) == "GoldenRatio = 1.61803398874989;\nGoldenRatio"
    #assert rcode(2*GoldenRatio) == "double const GoldenRatio = 1.61803398874989;\n2*GoldenRatio"
    #assert rcode(
    #    2*Catalan) == "double const Catalan = 0.915965594177219;\n2*Catalan"
    #assert rcode(2*EulerGamma) == "double const EulerGamma = 0.577215664901533;\n2*EulerGamma"

