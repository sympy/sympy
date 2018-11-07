from sympy.printing.lambdarepr import TensorflowPrinter
from sympy.printing.lambdarepr import _lambdify_tensorflow
from sympy import (eye, symbols, MatrixSymbol, Symbol)
from sympy.codegen.array_utils import (CodegenArrayContraction,
        CodegenArrayTensorProduct, CodegenArrayElementwiseAdd,
        CodegenArrayPermuteDims, CodegenArrayDiagonal)
from sympy.utilities.lambdify import lambdify

from sympy.utilities.pytest import skip
from sympy.external import import_module

tf = import_module("tensorflow")
np = import_module("numpy")

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

def test_codegen_einsum():
    if not np:
        skip("NumPy not installed")

    M = MatrixSymbol("M", 2, 2)
    N = MatrixSymbol("N", 2, 2)

    cg = CodegenArrayContraction.from_MatMul(M*N)
    f = lambdify((M, N), cg, 'numpy')

    ma = np.matrix([[1, 2], [3, 4]])
    mb = np.matrix([[1,-2], [-1, 3]])
    assert (f(ma, mb) == ma*mb).all()


def test_codegen_extra():
    if not np:
        skip("NumPy not installed")

    M = MatrixSymbol("M", 2, 2)
    N = MatrixSymbol("N", 2, 2)
    P = MatrixSymbol("P", 2, 2)
    Q = MatrixSymbol("Q", 2, 2)
    ma = np.matrix([[1, 2], [3, 4]])
    mb = np.matrix([[1,-2], [-1, 3]])
    mc = np.matrix([[2, 0], [1, 2]])
    md = np.matrix([[1,-1], [4, 7]])

    cg = CodegenArrayTensorProduct(M, N)
    assert tensorflow_code(cg) == 'tensorflow.einsum("ab,cd", M, N)'

    cg = CodegenArrayElementwiseAdd(M, N)
    assert tensorflow_code(cg) == 'tensorflow.add(M, N)'

    cg = CodegenArrayElementwiseAdd(M, N, P)
    assert tensorflow_code(cg) == 'tensorflow.add(tensorflow.add(M, N), P)'

    cg = CodegenArrayElementwiseAdd(M, N, P, Q)
    assert tensorflow_code(cg) == 'tensorflow.add(tensorflow.add(tensorflow.add(M, N), P), Q)'

    cg = CodegenArrayPermuteDims(M, [1, 0])
    assert tensorflow_code(cg) == 'tensorflow.transpose(M, [1, 0])'

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), [1, 2, 3, 0])
    assert tensorflow_code(cg) == 'tensorflow.transpose(tensorflow.einsum("ab,cd", M, N), [1, 2, 3, 0])'
    #f = lambdify((M, N), cg, 'tensorflow')
    #assert (f(ma, mb) == np.transpose(np.einsum(ma, [0, 1], mb, [2, 3]), (1, 2, 3, 0))).all()

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N), (1, 2))
    assert tensorflow_code(cg) == 'tensorflow.einsum("ab,bc->acb", M, N)'
    #f = lambdify((M, N), cg, 'tensorflow')
    #assert (f(ma, mb) == np.diagonal(np.einsum(ma, [0, 1], mb, [2, 3]), axis1=1, axis2=2)).all()
