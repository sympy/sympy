import random

from sympy.printing.torch import torch_code
from sympy import (eye, symbols, MatrixSymbol, Symbol, Matrix)
from sympy.codegen.array_utils import (CodegenArrayContraction,
        CodegenArrayTensorProduct, CodegenArrayElementwiseAdd,
        CodegenArrayPermuteDims, CodegenArrayDiagonal)
from sympy.utilities.lambdify import lambdify

from sympy.utilities.pytest import skip
from sympy.external import import_module

torch = import_module("torch")

M = MatrixSymbol("M", 3, 3)
N = MatrixSymbol("N", 3, 3)
P = MatrixSymbol("P", 3, 3)
Q = MatrixSymbol("Q", 3, 3)


def _compare_torch_matrix(variables, expr):
    f = lambdify(variables, expr, 'torch')
    random_matrices = [Matrix([[random.randint(0, 10) for k in
        range(i.shape[1])] for j in range(i.shape[0])]) for i in variables]
    random_variables = [eval(torch_code(i)) for i in
            random_matrices]
    r = f(*random_variables)
    e = expr.subs({k: v for k, v in zip(variables, random_matrices)}).doit()
    if e.is_Matrix:
        e = torch.FloatTensor(e.tolist())
    assert (r == e).all()


def test_torch_matrix():
    if torch is None:
        skip("Torch not installed")

    expr = M
    assert torch_code(expr) == "M"
    f = lambdify((M,), expr, "torch")
    assert f(eye(3)) == eye(3)

    expr = M*N
    assert torch_code(expr) == "torch.mm(M, N)"
    _compare_torch_matrix((M, N), expr)

    expr = M**3
    assert torch_code(expr) == "torch.mm(torch.mm(M, M), M)"
    _compare_torch_matrix((M,), expr)

    expr = M*N*P*Q
    assert torch_code(expr) == "torch.mm(torch.mm(torch.mm(M, N), P), Q)"
    _compare_torch_matrix((M, N, P, Q), expr)
