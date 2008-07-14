from sympy.core.symbol import symbols
from sympy.printing.precedence import precedence, PRECEDENCE
x, y ,z = symbols("xyz")

def test_add():
    assert precedence(x+y) == PRECEDENCE["Add"]
    assert precedence(x*z+1) == PRECEDENCE["Add"]

def test_mul():
    assert precedence(x*y) == PRECEDENCE["Mul"]
    assert precedence(-x*y) == PRECEDENCE["Add"]

def test_symbol():
    assert precedence(x) == PRECEDENCE["Atom"]

def test_pow():
    assert precedence(x**y) == PRECEDENCE["Pow"]