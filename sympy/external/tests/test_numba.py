# This testfile tests SymPy's compatibility with Numba

from sympy.external import import_module
from sympy import sin, exp, symbols, lambdify, Matrix

numba = import_module('numba')
if not numba:
    #bin/test will not execute any tests now
    disabled = True

def test_jit():
    # Test that numba.jit works with lambdify
    x = symbols('x')
    f1 = lambdify(x, sin(x), 'numpy')
    assert numba.jit(f1)(1) == f1(1)

    y = exp(-x)
    tanh = (1.0 - y) / (1.0 + y)
    ddx_tanh = tanh.diff(x)
    f2 = lambdify(x, ddx_tanh, 'numpy')
    assert numba.jit(f2)(1) == f2(1)

    from itertools import chain
    X = Matrix(symbols('x:3')).T
    Y = Matrix(symbols('y:3')).T
    f3 = lambdify((*X, *Y), X+Y)
    xs = (1,2,3,4,5,6)
    assert all(chain.from_iterable(f3(*xs) == numba.njit(f3)(*xs)))

