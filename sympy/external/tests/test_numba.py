# This testfile tests SymPy <-> Numba compatibility

from sympy.external import import_module


numba = import_module('numba')
if not numba:
    #bin/test will not execute any tests now
    disabled = True

from sympy import lambdify, Symbol, sin, Matrix, pi
from sympy.testing.pytest import skip

numpy = import_module('numpy')

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

# test that numba.jit works with lambdify

def test_trig():
    # test with a literal
    f1 = lambdify(x, sin(x), 'numpy')
    assert abs(numba.jit(f1)(1) - f1(1)) < 1e-15

    # test with a symbol
    f2 = lambdify(x, sin(x), 'math')
    assert abs(numba.jit(f2, forceobj=True)(pi) - f2(pi)) < 1e-15

def test_matrix():
    f = lambdify(x, Matrix([[x, 2*x], [1, 2]]), 'numpy')
    assert abs(numba.jit(f)(3) - f(3)).all() < 1e-15

def test_numpy_matrix():
    if not numpy:
        skip("numpy not installed.")
    A = Matrix([[x, x*y], [sin(z) + 4, x**z]])
    f = lambdify((x, y, z), A, ['numpy'])
    assert abs(numba.jit(f, forceobj=True)(1, 2, 3) - f(1, 2, 3)).all() < 1e-15
