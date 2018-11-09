from sympy.printing.tensorflow import TensorflowPrinter
from sympy.printing.tensorflow import tensorflow_code
from sympy import (eye, symbols, MatrixSymbol, Symbol, Matrix, symbols, sin,
        exp, Function, Derivative, Trace)
from sympy.codegen.array_utils import (CodegenArrayContraction,
        CodegenArrayTensorProduct, CodegenArrayElementwiseAdd,
        CodegenArrayPermuteDims, CodegenArrayDiagonal)
from sympy.utilities.lambdify import lambdify

from sympy.utilities.pytest import skip
from sympy.external import import_module

tf = import_module("tensorflow")

M = MatrixSymbol("M", 3, 3)
N = MatrixSymbol("N", 3, 3)
P = MatrixSymbol("P", 3, 3)
Q = MatrixSymbol("Q", 3, 3)

x, y, z, t = symbols("x y z t")


def test_tensorflow_matrix():
    if not tf:
        skip("TensorFlow not installed")

    assert tensorflow_code(eye(3)) == "tensorflow.constant([[1, 0, 0], [0, 1, 0], [0, 0, 1]])"

    expr = Matrix([[x, sin(y)], [exp(z), -t]])
    assert tensorflow_code(expr) == "tensorflow.Variable([[x, tensorflow.sin(y)], [tensorflow.exp(z), -t]])"

    expr = M
    assert tensorflow_code(expr) == "M"
    expr = M + N
    assert tensorflow_code(expr) == "M + N"
    expr = M*N
    assert tensorflow_code(expr) == "tensorflow.matmul(M, N)"
    expr = M*N*P*Q
    assert tensorflow_code(expr) == "tensorflow.matmul(tensorflow.matmul(tensorflow.matmul(M, N), P), Q)"
    expr = M**3
    assert tensorflow_code(expr) == "tensorflow.matmul(tensorflow.matmul(M, M), M)"
    expr = M.T
    assert tensorflow_code(expr) == "tensorflow.matrix_transpose(M)"
    expr = Trace(M)
    assert tensorflow_code(expr) == "tensorflow.trace(M)"


def test_codegen_einsum():
    if not tf:
        skip("TensorFlow not installed")

    session = tf.Session()

    M = MatrixSymbol("M", 2, 2)
    N = MatrixSymbol("N", 2, 2)

    cg = CodegenArrayContraction.from_MatMul(M*N)
    f = lambdify((M, N), cg, 'tensorflow')

    ma = tf.constant([[1, 2], [3, 4]])
    mb = tf.constant([[1,-2], [-1, 3]])
    y = session.run(f(ma, mb))
    c = session.run(tf.matmul(ma, mb))
    assert (y == c).all()


def test_codegen_extra():
    if not tf:
        skip("TensorFlow not installed")

    session = tf.Session()

    M = MatrixSymbol("M", 2, 2)
    N = MatrixSymbol("N", 2, 2)
    P = MatrixSymbol("P", 2, 2)
    Q = MatrixSymbol("Q", 2, 2)
    ma = tf.constant([[1, 2], [3, 4]])
    mb = tf.constant([[1,-2], [-1, 3]])
    mc = tf.constant([[2, 0], [1, 2]])
    md = tf.constant([[1,-1], [4, 7]])

    cg = CodegenArrayTensorProduct(M, N)
    assert tensorflow_code(cg) == 'tensorflow.einsum("ab,cd", M, N)'
    f = lambdify((M, N), cg, 'tensorflow')
    y = session.run(f(ma, mb))
    c = session.run(tf.einsum("ij,kl", ma, mb))
    assert (y == c).all()

    cg = CodegenArrayElementwiseAdd(M, N)
    assert tensorflow_code(cg) == 'tensorflow.add(M, N)'
    f = lambdify((M, N), cg, 'tensorflow')
    y = session.run(f(ma, mb))
    c = session.run(ma + mb)
    assert (y == c).all()

    cg = CodegenArrayElementwiseAdd(M, N, P)
    assert tensorflow_code(cg) == 'tensorflow.add(tensorflow.add(M, N), P)'
    f = lambdify((M, N, P), cg, 'tensorflow')
    y = session.run(f(ma, mb, mc))
    c = session.run(ma + mb + mc)
    assert (y == c).all()

    cg = CodegenArrayElementwiseAdd(M, N, P, Q)
    assert tensorflow_code(cg) == 'tensorflow.add(tensorflow.add(tensorflow.add(M, N), P), Q)'
    f = lambdify((M, N, P, Q), cg, 'tensorflow')
    y = session.run(f(ma, mb, mc, md))
    c = session.run(ma + mb + mc + md)
    assert (y == c).all()

    cg = CodegenArrayPermuteDims(M, [1, 0])
    assert tensorflow_code(cg) == 'tensorflow.transpose(M, [1, 0])'
    f = lambdify((M,), cg, 'tensorflow')
    y = session.run(f(ma))
    c = session.run(tf.transpose(ma))
    assert (y == c).all()

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), [1, 2, 3, 0])
    assert tensorflow_code(cg) == 'tensorflow.transpose(tensorflow.einsum("ab,cd", M, N), [1, 2, 3, 0])'
    f = lambdify((M, N), cg, 'tensorflow')
    y = session.run(f(ma, mb))
    c = session.run(tf.transpose(tf.einsum("ab,cd", ma, mb), [1, 2, 3, 0]))
    assert (y == c).all()

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N), (1, 2))
    assert tensorflow_code(cg) == 'tensorflow.einsum("ab,bc->acb", M, N)'
    f = lambdify((M, N), cg, 'tensorflow')
    y = session.run(f(ma, mb))
    c = session.run(tf.einsum("ab,bc->acb", ma, mb))
    assert (y == c).all()


def test_MatrixElement_printing():
    A = MatrixSymbol("A", 1, 3)
    B = MatrixSymbol("B", 1, 3)
    C = MatrixSymbol("C", 1, 3)

    assert tensorflow_code(A[0, 0]) == "A[0, 0]"
    assert tensorflow_code(3 * A[0, 0]) == "3*A[0, 0]"

    F = C[0, 0].subs(C, A - B)
    assert tensorflow_code(F) == "((-1)*B + A)[0, 0]"


def test_tensorflow_Derivative():
    f = Function("f")

    expr = Derivative(sin(x), x)
    assert tensorflow_code(expr) == "tensorflow.gradients(tensorflow.sin(x), x)[0]"
