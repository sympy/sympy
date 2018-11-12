from sympy.printing.torch import torch_code
from sympy import MatrixSymbol
from sympy import (eye, symbols, MatrixSymbol, Symbol)
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


def test_torch_matrix():
    if torch is None:
        skip("Torch not installed")

    expr = M
    assert torch_code(expr) == "M"
    expr = M*N
    assert torch_code(expr) == "torch.mm(M, N)"
    expr = M**3
    assert torch_code(expr) == "torch.mm(torch.mm(M, M), M)"
