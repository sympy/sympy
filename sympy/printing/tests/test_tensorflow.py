from sympy.printing.lambdarepr import TensorflowPrinter
from sympy.printing.lambdarepr import _lambdify_tensorflow
from sympy import (eye, symbols, MatrixSymbol, Symbol)

from sympy.utilities.pytest import skip
from sympy.external import import_module

tf = import_module("tensorflow")

def tensorflow_code(expr):
    printer = TensorflowPrinter()
    return printer.doprint(expr)


def test_tensorflow_matrix():
    if not tf:
        skip("TensorFlow not installed")

    assert tensorflow_code(eye(3)) == "tensorflow.Variable([[1, 0, 0], [0, 1, 0], [0, 0, 1]])"

    M = MatrixSymbol("M", 3, 3)
    N = MatrixSymbol("N", 3, 3)

    expr = M
    assert tensorflow_code(expr) == "M"

    assert _lambdify_tensorflow((M,), expr) == """\
def _lambdified(M_):
    M = tensorflow.placeholder(tensorflow.float32, (3, 3))
    feed_dict = {M: M_}
    return tensorflow.InteractiveSession().run(M, feed_dict=feed_dict)
"""
    expr = M + N
    assert _lambdify_tensorflow((M, N), expr) == """\
def _lambdified(M_, N_):
    M = tensorflow.placeholder(tensorflow.float32, (3, 3))
    N = tensorflow.placeholder(tensorflow.float32, (3, 3))
    feed_dict = {M: M_, N: N_}
    return tensorflow.InteractiveSession().run(M + N, feed_dict=feed_dict)
"""
