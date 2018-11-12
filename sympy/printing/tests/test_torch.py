import random

from sympy.printing.torch import torch_code
from sympy import (eye, symbols, MatrixSymbol, Symbol, Matrix, tensorproduct)
from sympy.tensor.array import NDimArray
from sympy.codegen.array_utils import (CodegenArrayContraction,
        CodegenArrayTensorProduct, CodegenArrayElementwiseAdd,
        CodegenArrayPermuteDims, CodegenArrayDiagonal, _CodegenArrayAbstract)
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
    if isinstance(e, _CodegenArrayAbstract):
        e = e.doit()
    if e.is_Matrix or isinstance(e, NDimArray):
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


def test_codegen_extra():
    if not torch:
        skip("Torch not installed")

    M = MatrixSymbol("M", 2, 2)
    N = MatrixSymbol("N", 2, 2)
    P = MatrixSymbol("P", 2, 2)
    Q = MatrixSymbol("Q", 2, 2)
    ma = torch.FloatTensor([[1, 2], [3, 4]])
    mb = torch.FloatTensor([[1,-2], [-1, 3]])
    mc = torch.FloatTensor([[2, 0], [1, 2]])
    md = torch.FloatTensor([[1,-1], [4, 7]])

    cg = CodegenArrayTensorProduct(M, N)
    assert torch_code(cg) == 'torch.einsum("ab,cd", [M, N])'
    _compare_torch_matrix((M, N), cg)

    cg = CodegenArrayElementwiseAdd(M, N)
    assert torch_code(cg) == 'torch.add(M, N)'
    _compare_torch_matrix((M, N), cg)

    cg = CodegenArrayElementwiseAdd(M, N, P)
    assert torch_code(cg) == 'torch.add(torch.add(M, N), P)'
    _compare_torch_matrix((M, N, P), cg)

    cg = CodegenArrayElementwiseAdd(M, N, P, Q)
    assert torch_code(cg) == 'torch.add(torch.add(torch.add(M, N), P), Q)'
    _compare_torch_matrix((M, N, P, Q), cg)

    cg = CodegenArrayPermuteDims(M, [1, 0])
    assert torch_code(cg) == 'M.permute(1, 0)'
    _compare_torch_matrix((M,), cg)

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), [1, 2, 3, 0])
    assert torch_code(cg) == 'torch.einsum("ab,cd", [M, N]).permute(1, 2, 3, 0)'
    _compare_torch_matrix((M, N), cg)

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N), (1, 2))
    assert torch_code(cg) == 'torch.einsum("ab,bc->acb", [M, N])'
    # TODO: no implemented for diagonal
    #_compare_torch_matrix((M, N), cg)
