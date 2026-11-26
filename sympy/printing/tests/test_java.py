from sympy.core import oo
from sympy.printing.java import javacode

def test_constants():
    assert javacode(oo) == "Double.POSITIVE_INFINITY"
