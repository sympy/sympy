from sympy import Symbol, normal
from sympy.abc import n

def test_normal():
    assert normal(4*n+5, 2*(4*n+1)*(2*n+3), n)

def test_gosper():
    pass
