from sympy import Q, symbols
from sympy.printing import sstr, pretty, latex

x,y,z = symbols('x y z')

def test_printing():
    eq = Q.eq(x,y)
    assert sstr(eq) == "x = y"
    assert pretty(eq) == "x = y"
    assert latex(eq) == "x = y"
