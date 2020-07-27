from sympy import (
    Ring, S, scalar_add, scalar_mul, symbols
)

x, y = symbols('x y') 
R = Ring('R', (S.Complexes,), (scalar_add, scalar_mul))

def test_Ring():
    assert R.add(x, x, evaluate=True) == scalar_mul(2, x)
