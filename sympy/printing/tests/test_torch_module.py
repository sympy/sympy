"""
Unit tests for the torch_module function and lambdify with "torch_module" in sympy.printing.pytorch.
"""

import random
import math
import numpy as np

import sympy
from sympy.testing.pytest import skip
from sympy.external import import_module

from sympy import symbols, Derivative
from sympy.printing.pytorch import torch_module
from sympy.utilities.lambdify import lambdify
from sympy import (eye, MatrixSymbol, Matrix)
from sympy.tensor.array import NDimArray
from sympy.tensor.array.expressions.array_expressions import (
    ArrayTensorProduct, ArrayAdd, PermuteDims, ArrayDiagonal, _CodegenArrayAbstract)
from sympy.core.relational import Eq, Ne, Ge, Gt, Le, Lt
from sympy.functions import (
    Abs, ceiling, exp, floor, sign, sin, asin, cos, acos, tan, atan, atan2,
    cosh, acosh, sinh, asinh, tanh, atanh, re, im, arg, erf, loggamma, sqrt)
from sympy.matrices.expressions import Determinant, HadamardProduct, Inverse, Trace
from sympy.matrices import randMatrix
from sympy.matrices import Identity, ZeroMatrix, OneMatrix
from sympy import conjugate, I
from sympy import Heaviside, gamma, polygamma


torch = import_module("torch")

if torch is not None:
    m3x3sympy = Matrix([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
    m3x3 = torch.tensor(m3x3sympy.tolist(), dtype=torch.float64)


M = MatrixSymbol("M", 3, 3)
N = MatrixSymbol("N", 3, 3)
P = MatrixSymbol("P", 3, 3)
Q = MatrixSymbol("Q", 3, 3)

x, y, z, t = symbols("x y z t")


def _compare_torch_module_matrix(variables, expr, use_lambdify=False):
    if not torch:
        skip("PyTorch not installed")

    if use_lambdify:
        mod = lambdify(variables, expr, "torch_module")
    else:
        mod = torch_module(variables, expr)
    random_matrices = [randMatrix(i.shape[0], i.shape[1]) for i in variables]
    random_tensors = [torch.tensor(i.tolist(), dtype=torch.float64) for i in random_matrices]
    result = mod(*random_tensors)
    expected = expr.subs(dict(zip(variables, random_matrices))).doit()

    if isinstance(expected, _CodegenArrayAbstract):
        expected = expected.doit()

    if hasattr(expected, 'is_number') and expected.is_number:
        if isinstance(result, torch.Tensor) and result.dim() == 0:
            result = result.item()
            expected = float(expected)
            assert abs(result - expected) < 1e-6
            return

    if expected.is_Matrix or isinstance(expected, NDimArray):
        expected = torch.tensor(expected.tolist(), dtype=torch.float64)
        assert torch.allclose(result, expected, atol=1e-6)
    else:
        raise TypeError(f"Cannot compare {type(result)} with {type(expected)}")

def _compare_torch_module_scalar(variables, expr, rng=lambda: random.uniform(-5, 5), use_lambdify=False):
    if not torch:
        skip("PyTorch not installed")

    if use_lambdify:
        mod = lambdify(variables, expr, "torch_module")
    else:
        mod = torch_module(variables, expr)
    rvs = [rng() for v in variables]
    t_rvs = [torch.tensor(i, dtype=torch.float64) for i in rvs]
    result = mod(*t_rvs)
    if isinstance(result, torch.Tensor) and result.dim() == 0:
        result = result.item()
    expected = expr.subs(dict(zip(variables, rvs))).doit()
    assert abs(result - expected) < 1e-6

def _compare_torch_module_relational(variables, expr, rng=lambda: random.randint(0, 10), use_lambdify=False):
    if not torch:
        skip("PyTorch not installed")

    if use_lambdify:
        mod = lambdify(variables, expr, "torch_module")
    else:
        mod = torch_module(variables, expr)
    rvs = [rng() for v in variables]
    t_rvs = [torch.tensor(i, dtype=torch.float64) for i in rvs]
    result = mod(*t_rvs)
    expected = bool(expr.subs(dict(zip(variables, rvs))).doit())
    assert result.item() == expected

def test_torch_module_math():
    if not torch:
        skip("PyTorch not installed")

    expr = Abs(x)
    mod = torch_module(x, expr)
    ma = torch.tensor([[-1, 2, -3, -4]], dtype=torch.float64)
    result = mod(ma)
    expected = torch.abs(ma)
    assert torch.all(result == expected)

    expr = sign(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-10, 10))

    expr = ceiling(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random())

    expr = floor(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random())

    expr = exp(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-2, 2))

    expr = sqrt(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random())

    expr = x ** 4
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random())

    expr = cos(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random())

    expr = acos(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-0.99, 0.99))

    expr = sin(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random())

    expr = asin(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-0.99, 0.99))

    expr = tan(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-1.5, 1.5))

    expr = atan(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-5, 5))

    expr = atan2(y, x)
    _compare_torch_module_scalar((y, x), expr, rng=lambda: random.uniform(-5, 5))

    expr = cosh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-2, 2))

    expr = acosh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(1.1, 5))

    expr = sinh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-2, 2))

    expr = asinh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-5, 5))

    expr = tanh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-2, 2))

    expr = atanh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-0.9, 0.9))

    expr = erf(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-2, 2))

    expr = loggamma(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(0.5, 5))

def test_torch_module_lambdify_math():
    if not torch:
        skip("PyTorch not installed")

    expr = Abs(x)
    mod = lambdify(x, expr, "torch_module")
    ma = torch.tensor([[-1, 2, -3, -4]], dtype=torch.float64)
    result = mod(ma)
    expected = torch.abs(ma)
    assert torch.all(result == expected)

    expr = sign(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-10, 10), use_lambdify=True)

    expr = ceiling(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random(), use_lambdify=True)

    expr = floor(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random(), use_lambdify=True)

    expr = exp(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-2, 2), use_lambdify=True)

    expr = sqrt(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random(), use_lambdify=True)

    expr = x ** 4
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random(), use_lambdify=True)

    expr = cos(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random(), use_lambdify=True)

    expr = acos(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-0.99, 0.99), use_lambdify=True)

    expr = sin(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.random(), use_lambdify=True)

    expr = asin(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-0.99, 0.99), use_lambdify=True)

    expr = tan(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-1.5, 1.5), use_lambdify=True)

    expr = atan(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-5, 5), use_lambdify=True)

    expr = atan2(y, x)
    _compare_torch_module_scalar((y, x), expr, rng=lambda: random.uniform(-5, 5), use_lambdify=True)

    expr = cosh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-2, 2), use_lambdify=True)

    expr = acosh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(1.1, 5), use_lambdify=True)

    expr = sinh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-2, 2), use_lambdify=True)

    expr = asinh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-5, 5), use_lambdify=True)

    expr = tanh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-2, 2), use_lambdify=True)

    expr = atanh(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-0.9, 0.9), use_lambdify=True)

    expr = erf(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(-2, 2), use_lambdify=True)

    expr = loggamma(x)
    _compare_torch_module_scalar((x,), expr, rng=lambda: random.uniform(0.5, 5), use_lambdify=True)

def test_torch_module_complexes():
    if not torch:
        skip("PyTorch not installed")

    expr = re(x)
    mod = torch_module(x, expr)
    input_tensor = torch.tensor([1.0 + 2.0j], dtype=torch.complex128)
    result = mod(input_tensor)
    assert torch.allclose(result, torch.real(input_tensor))

    expr = im(x)
    mod = torch_module(x, expr)
    result = mod(input_tensor)
    assert torch.allclose(result, torch.imag(input_tensor))

    expr = arg(x)
    mod = torch_module(x, expr)
    result = mod(input_tensor)
    assert torch.allclose(result, torch.angle(input_tensor))

def test_torch_module_lambdify_complexes():
    if not torch:
        skip("PyTorch not installed")

    expr = re(x)
    mod = lambdify(x, expr, "torch_module")
    input_tensor = torch.tensor([1.0 + 2.0j], dtype=torch.complex128)
    result = mod(input_tensor)
    assert torch.allclose(result, torch.real(input_tensor))

    expr = im(x)
    mod = lambdify(x, expr, "torch_module")
    result = mod(input_tensor)
    assert torch.allclose(result, torch.imag(input_tensor))

    expr = arg(x)
    mod = lambdify(x, expr, "torch_module")
    result = mod(input_tensor)
    assert torch.allclose(result, torch.angle(input_tensor))

def test_torch_module_relational():
    if not torch:
        skip("PyTorch not installed")

    expr = Eq(x, y)
    _compare_torch_module_relational((x, y), expr)

    expr = Ne(x, y)
    _compare_torch_module_relational((x, y), expr)

    expr = Ge(x, y)
    _compare_torch_module_relational((x, y), expr)

    expr = Gt(x, y)
    _compare_torch_module_relational((x, y), expr)

    expr = Le(x, y)
    _compare_torch_module_relational((x, y), expr)

    expr = Lt(x, y)
    _compare_torch_module_relational((x, y), expr)

def test_torch_module_lambdify_relational():
    if not torch:
        skip("PyTorch not installed")

    expr = Eq(x, y)
    _compare_torch_module_relational((x, y), expr, use_lambdify=True)

    expr = Ne(x, y)
    _compare_torch_module_relational((x, y), expr, use_lambdify=True)

    expr = Ge(x, y)
    _compare_torch_module_relational((x, y), expr, use_lambdify=True)

    expr = Gt(x, y)
    _compare_torch_module_relational((x, y), expr, use_lambdify=True)

    expr = Le(x, y)
    _compare_torch_module_relational((x, y), expr, use_lambdify=True)

    expr = Lt(x, y)
    _compare_torch_module_relational((x, y), expr, use_lambdify=True)

def test_torch_module_matrix():
    if not torch:
        skip("PyTorch not installed")

    expr = M
    mod = torch_module((M,), expr)
    eye_tensor = torch.tensor(eye(3).tolist(), dtype=torch.float64)
    result = mod(eye_tensor)
    assert torch.allclose(result, eye_tensor)

    expr = M * N
    _compare_torch_module_matrix((M, N), expr)

    expr = M ** 3
    _compare_torch_module_matrix((M,), expr)

    expr = M * N * P * Q
    _compare_torch_module_matrix((M, N, P, Q), expr)

    expr = Trace(M)
    _compare_torch_module_matrix((M,), expr)

    expr = Determinant(M)
    _compare_torch_module_matrix((M,), expr)

    expr = HadamardProduct(M, N)
    _compare_torch_module_matrix((M, N), expr)

    expr = Inverse(M)
    mod = torch_module((M,), expr)
    result = mod(eye_tensor)
    expected = torch.linalg.inv(eye_tensor)
    assert torch.allclose(result, expected)

def test_torch_module_lambdify_matrix():
    if not torch:
        skip("PyTorch not installed")

    expr = M
    mod = lambdify((M,), expr, "torch_module")
    eye_tensor = torch.tensor(eye(3).tolist(), dtype=torch.float64)
    result = mod(eye_tensor)
    assert torch.allclose(result, eye_tensor)

    expr = M * N
    _compare_torch_module_matrix((M, N), expr, use_lambdify=True)

    expr = M ** 3
    _compare_torch_module_matrix((M,), expr, use_lambdify=True)

    expr = M * N * P * Q
    _compare_torch_module_matrix((M, N, P, Q), expr, use_lambdify=True)

    expr = Trace(M)
    _compare_torch_module_matrix((M,), expr, use_lambdify=True)

    expr = Determinant(M)
    _compare_torch_module_matrix((M,), expr, use_lambdify=True)

    expr = HadamardProduct(M, N)
    _compare_torch_module_matrix((M, N), expr, use_lambdify=True)

    expr = Inverse(M)
    mod = lambdify((M,), expr, "torch_module")
    result = mod(eye_tensor)
    expected = torch.linalg.inv(eye_tensor)
    assert torch.allclose(result, expected)

def test_torch_module_array_operations():
    if not torch:
        skip("PyTorch not installed")

    M = MatrixSymbol("M", 2, 2)
    N = MatrixSymbol("N", 2, 2)
    P = MatrixSymbol("P", 2, 2)
    Q = MatrixSymbol("Q", 2, 2)

    ma = torch.tensor([[1., 2.], [3., 4.]], dtype=torch.float64)
    mb = torch.tensor([[1., -2.], [-1., 3.]], dtype=torch.float64)
    mc = torch.tensor([[2., 0.], [1., 2.]], dtype=torch.float64)
    md = torch.tensor([[1., -1.], [4., 7.]], dtype=torch.float64)

    cg = ArrayTensorProduct(M, N)
    mod = torch_module((M, N), cg)
    result = mod(ma, mb)
    expected = torch.einsum("ij,kl", ma, mb)

    assert torch.allclose(result, expected)

    cg = ArrayAdd(M, N)
    mod = torch_module((M, N), cg)
    result = mod(ma, mb)
    expected = ma + mb

    assert torch.allclose(result, expected)

    cg = ArrayAdd(M, N, P)
    mod = torch_module((M, N, P), cg)
    result = mod(ma, mb, mc)
    expected = ma + mb + mc
    assert torch.allclose(result, expected)

    cg = ArrayAdd(M, N, P, Q)
    mod = torch_module((M, N, P, Q), cg)
    result = mod(ma, mb, mc, md)
    expected = ma + mb + mc + md
    assert torch.allclose(result, expected)

    cg = PermuteDims(M, [1, 0])
    mod = torch_module((M,), cg)
    result = mod(ma)
    expected = ma.T
    assert torch.allclose(result, expected)

    cg = PermuteDims(ArrayTensorProduct(M, N), [1, 2, 3, 0])
    mod = torch_module((M, N), cg)
    result = mod(ma, mb)
    expected = torch.einsum("ab,cd", ma, mb).permute(1, 2, 3, 0)
    assert torch.allclose(result, expected)

    cg = ArrayDiagonal(ArrayTensorProduct(M, N), (1, 2))
    mod = torch_module((M, N), cg)
    result = mod(ma, mb)
    expected = torch.tensor([[[1., -2.], [-2., 6.]], [[3., -6.], [-4., 12.]]], dtype=torch.float64)
    assert torch.allclose(result, expected)

    cg = ArrayDiagonal(ArrayTensorProduct(M, N), (1, 2))
    mod = lambdify((M, N), cg, "torch_module")
    result = mod(ma, mb)
    expected = torch.einsum("ab,bc->abc", ma, mb)
    assert torch.allclose(result, expected)

def test_torch_module_lambdify_array_operations():
    if not torch:
        skip("PyTorch not installed")

    M = MatrixSymbol("M", 2, 2)
    N = MatrixSymbol("N", 2, 2)
    P = MatrixSymbol("P", 2, 2)
    Q = MatrixSymbol("Q", 2, 2)

    ma = torch.tensor([[1., 2.], [3., 4.]], dtype=torch.float64)
    mb = torch.tensor([[1., -2.], [-1., 3.]], dtype=torch.float64)
    mc = torch.tensor([[2., 0.], [1., 2.]], dtype=torch.float64)
    md = torch.tensor([[1., -1.], [4., 7.]], dtype=torch.float64)

    cg = ArrayTensorProduct(M, N)
    mod = lambdify((M, N), cg, "torch_module")
    result = mod(ma, mb)
    expected = torch.einsum("ij,kl", ma, mb)
    assert torch.allclose(result, expected)

    cg = ArrayAdd(M, N)
    mod = lambdify((M, N), cg, "torch_module")
    result = mod(ma, mb)
    expected = ma + mb
    assert torch.allclose(result, expected)

    cg = ArrayAdd(M, N, P)
    mod = lambdify((M, N, P), cg, "torch_module")
    result = mod(ma, mb, mc)
    expected = ma + mb + mc
    assert torch.allclose(result, expected)

    cg = ArrayAdd(M, N, P, Q)
    mod = lambdify((M, N, P, Q), cg, "torch_module")
    result = mod(ma, mb, mc, md)
    expected = ma + mb + mc + md
    assert torch.allclose(result, expected)

    cg = PermuteDims(M, [1, 0])
    mod = lambdify((M,), cg, "torch_module")
    result = mod(ma)
    expected = ma.T
    assert torch.allclose(result, expected)

    cg = PermuteDims(ArrayTensorProduct(M, N), [1, 2, 3, 0])
    mod = lambdify((M, N), cg, "torch_module")
    result = mod(ma, mb)
    expected = torch.einsum("ab,cd", ma, mb).permute(1, 2, 3, 0)
    assert torch.allclose(result, expected)

    cg = ArrayDiagonal(ArrayTensorProduct(M, N), (1, 2))
    mod = lambdify((M, N), cg, "torch_module")
    result = mod(ma, mb)
    expected = torch.einsum("ab,bc->abc", ma, mb)
    assert torch.allclose(result, expected)

def test_torch_module_derivative():
    if not torch:
        skip("PyTorch not installed")

    expr = Derivative(sin(x), x)
    try:
        torch_module(x, expr)
        assert False, "Expected ValueError for unsupported Derivative"
    except ValueError as e:
        assert str(e) == "Unsupported SymPy function: Derivative"

def test_torch_module_lambdify_derivative():
    if not torch:
        skip("PyTorch not installed")

    expr = Derivative(sin(x), x)
    try:
        lambdify(x, expr, "torch_module")
        assert False, "Expected ValueError for unsupported Derivative"
    except ValueError as e:
        assert str(e) == "Unsupported SymPy function: Derivative"

def test_torch_module_requires_grad():
    if not torch:
        skip("PyTorch not installed")

    expr = sin(x) + cos(y)
    mod = torch_module([x, y], expr)
    x_val = torch.tensor(1.0, requires_grad=True, dtype=torch.float64)
    y_val = torch.tensor(2.0, requires_grad=True, dtype=torch.float64)
    result = mod(x_val, y_val)
    assert result.requires_grad
    result.backward()
    assert torch.isclose(x_val.grad, torch.cos(x_val), atol=1e-6)
    assert torch.isclose(y_val.grad, -torch.sin(y_val), atol=1e-6)

def test_torch_module_lambdify_requires_grad():
    if not torch:
        skip("PyTorch not installed")

    expr = sin(x) + cos(y)
    mod = lambdify([x, y], expr, "torch_module")
    x_val = torch.tensor(1.0, requires_grad=True, dtype=torch.float64)
    y_val = torch.tensor(2.0, requires_grad=True, dtype=torch.float64)
    result = mod(x_val, y_val)
    assert result.requires_grad
    result.backward()
    assert torch.isclose(x_val.grad, torch.cos(x_val), atol=1e-6)
    assert torch.isclose(y_val.grad, -torch.sin(y_val), atol=1e-6)

def test_torch_module_lambdify_requires_grad_low_precision():
    if not torch:
        skip("PyTorch not installed")

    expr = sin(x) + cos(y)
    mod = lambdify([x, y], expr, "torch_module")
    x_val = torch.tensor(1.0, requires_grad=True, dtype=torch.float32)
    y_val = torch.tensor(2.0, requires_grad=True, dtype=torch.float32)
    result = mod(x_val, y_val)
    assert result.requires_grad
    result.backward()
    assert torch.isclose(x_val.grad, torch.cos(x_val), atol=1e-6)
    assert torch.isclose(y_val.grad, -torch.sin(y_val), atol=1e-6)

def test_torch_module_special_matrices():
    if not torch:
        skip("PyTorch not installed")

    expr = Identity(3)
    mod = torch_module([], expr)
    result = mod()
    expected = torch.eye(3, dtype=torch.float64)
    assert torch.allclose(result, expected)

    expr = ZeroMatrix(2, 3)
    mod = torch_module([], expr)
    result = mod()
    expected = torch.zeros((2, 3), dtype=torch.float64)
    assert torch.allclose(result, expected)

    expr = OneMatrix(2, 3)
    mod = torch_module([], expr)
    result = mod()
    expected = torch.ones((2, 3), dtype=torch.float64)
    assert torch.allclose(result, expected)

def test_torch_module_lambdify_special_matrices():
    if not torch:
        skip("PyTorch not installed")

    expr = Identity(3)
    mod = lambdify([], expr, "torch_module")
    result = mod()
    expected = torch.eye(3, dtype=torch.float64)
    assert torch.allclose(result, expected)

    expr = ZeroMatrix(2, 3)
    mod = lambdify([], expr, "torch_module")
    result = mod()
    expected = torch.zeros((2, 3), dtype=torch.float64)
    assert torch.allclose(result, expected)

    expr = OneMatrix(2, 3)
    mod = lambdify([], expr, "torch_module")
    result = mod()
    expected = torch.ones((2, 3), dtype=torch.float64)
    assert torch.allclose(result, expected)

def test_torch_module_complex_operations():
    if not torch:
        skip("PyTorch not installed")

    expr = conjugate(x)
    mod = torch_module(x, expr)
    input_tensor = torch.tensor([1.0 + 2.0j], dtype=torch.complex64)
    result = mod(input_tensor)
    assert torch.allclose(result, torch.conj(input_tensor))

    expr = conjugate(sin(x) + I * cos(y))
    mod = torch_module([x, y], expr)
    x_val = torch.tensor(1.0, dtype=torch.float64)
    y_val = torch.tensor(2.0, dtype=torch.float64)
    result = mod(x_val, y_val)
    expected = torch.sin(x_val) - 1j * torch.cos(y_val)
    assert torch.allclose(result, expected, atol=1e-6)

    expr = I
    mod = torch_module([], expr)
    result = mod()
    expected = torch.tensor(1j, dtype=torch.complex64)
    assert torch.allclose(result, expected)

    expr = 2 * I + x
    mod = torch_module(x, expr)
    x_val = torch.tensor(3.0, dtype=torch.float64)
    result = mod(x_val)
    expected = x_val + 2 * 1j
    assert torch.allclose(result, expected)

    expr = exp(I * x)
    mod = torch_module(x, expr)
    x_val = torch.tensor(1.0, dtype=torch.float64)
    result = mod(x_val)
    expected = torch.exp(1j * x_val)
    assert torch.allclose(result, expected)

def test_torch_module_lambdify_complex_operations():
    if not torch:
        skip("PyTorch not installed")

    expr = conjugate(x)
    mod = lambdify(x, expr, "torch_module")
    input_tensor = torch.tensor([1.0 + 2.0j], dtype=torch.complex64)
    result = mod(input_tensor)
    assert torch.allclose(result, torch.conj(input_tensor))

    expr = conjugate(sin(x) + I * cos(y))
    mod = lambdify([x, y], expr, "torch_module")
    x_val = torch.tensor(1.0, dtype=torch.float64)
    y_val = torch.tensor(2.0, dtype=torch.float64)
    result = mod(x_val, y_val)
    expected = torch.sin(x_val) - 1j * torch.cos(y_val)
    assert torch.allclose(result, expected, atol=1e-6)

    expr = I
    mod = lambdify([], expr, "torch_module")
    result = mod()
    expected = torch.tensor(1j, dtype=torch.complex64)
    assert torch.allclose(result, expected)

    expr = 2 * I + x
    mod = lambdify(x, expr, "torch_module")
    x_val = torch.tensor(3.0, dtype=torch.float64)
    result = mod(x_val)
    expected = x_val + 2 * 1j
    assert torch.allclose(result, expected)

    expr = exp(I * x)
    mod = lambdify(x, expr, "torch_module")
    x_val = torch.tensor(1.0, dtype=torch.float64)
    result = mod(x_val)
    expected = torch.exp(1j * x_val)
    assert torch.allclose(result, expected)

def test_torch_module_special_functions():
    if not torch:
        skip("PyTorch not installed")

    expr = Heaviside(x)
    mod = torch_module(x, expr)
    x_val = torch.tensor([-1.0, 0.0, 1.0], dtype=torch.float64)
    result = mod(x_val)
    expected = torch.heaviside(x_val, torch.tensor(0.5, dtype=torch.float64))
    assert torch.allclose(result, expected)

    expr = Heaviside(x, 0)
    mod = torch_module(x, expr)
    result = mod(x_val)
    expected = torch.heaviside(x_val, torch.tensor(0.0, dtype=torch.float64))
    assert torch.allclose(result, expected)

    expr = gamma(x)
    try:
        torch_module(x, expr)
        assert False, "Expected ValueError for unsupported gamma"
    except ValueError as e:
        assert str(e) == "Unsupported SymPy function: gamma"

    expr = polygamma(0, x)
    try:
        torch_module(x, expr)
        assert False, "Expected ValueError for unsupported polygamma"
    except ValueError as e:
        assert str(e) == "Unsupported SymPy function: polygamma"

def test_torch_module_lambdify_special_functions():
    if not torch:
        skip("PyTorch not installed")

    expr = Heaviside(x)
    mod = lambdify(x, expr, "torch_module")
    x_val = torch.tensor([-1.0, 0.0, 1.0], dtype=torch.float64)
    result = mod(x_val)
    expected = torch.heaviside(x_val, torch.tensor(0.5, dtype=torch.float64))
    assert torch.allclose(result, expected)

    expr = Heaviside(x, 0)
    mod = lambdify(x, expr, "torch_module")
    result = mod(x_val)
    expected = torch.heaviside(x_val, torch.tensor(0.0, dtype=torch.float64))
    assert torch.allclose(result, expected)

    expr = gamma(x)
    try:
        lambdify(x, expr, "torch_module")
        assert False, "Expected ValueError for unsupported gamma"
    except ValueError as e:
        assert str(e) == "Unsupported SymPy function: gamma"

    expr = polygamma(0, x)
    try:
        lambdify(x, expr, "torch_module")
        assert False, "Expected ValueError for unsupported polygamma"
    except ValueError as e:
        assert str(e) == "Unsupported SymPy function: polygamma"

def test_example():
    x = sympy.symbols("x_name")
    cosx = 1.0 * sympy.cos(x)
    sinx = 2.0 * sympy.sin(x)

    mod = torch_module([x], [cosx, sinx])

    x_ = torch.rand(3)
    out = mod(x_)

    assert torch.equal(out[:, 0], x_.cos())
    assert torch.equal(out[:, 1], 2 * x_.sin())
    assert out.requires_grad

def test_grad():
    x = sympy.symbols("x_name")
    y = 1.0 * x
    mod = torch_module([x], [y])
    x_val = torch.ones((), requires_grad=True)
    out = mod(x_val)
    out.backward()
    # torch_module doesn't expose parameters like sympytorch, so we test gradient w.r.t. input
    assert torch.isclose(x_val.grad, torch.tensor(1.0), atol=1e-6)

def test_reduce():
    x, y = sympy.symbols("x y")
    z = 2 * x * y
    mod = torch_module([x, y], [z])
    result = mod(torch.rand(2), torch.rand(2))
    assert result.shape == (2,)

    z = 2 + x + y
    mod = torch_module([x, y], [z])
    result = mod(torch.rand(2), torch.rand(2))
    assert result.shape == (2,)

def test_special_subclasses():
    x, y = sympy.symbols("x y")
    z = x - 1
    w = y * 0
    u = sympy.Integer(1)

    mod = torch_module([x, y], [z, w, u])
    x_val = torch.tensor([2.0], dtype=torch.float64)
    y_val = torch.tensor([3.0], dtype=torch.float64)
    out = mod(x_val, y_val)
    # Adjust expected shape to (1, 3) since expressions are stacked horizontally
    assert torch.allclose(out, torch.tensor([[1.0, 0.0, 1.0]], dtype=torch.float64), atol=1e-6)

def test_constants():
    x = sympy.symbols("x")
    y = 2.0 * x + sympy.UnevaluatedExpr(1.0)
    mod = torch_module([x], [y])
    x_val = torch.tensor([3.0], dtype=torch.float64)
    out = mod(x_val)
    assert torch.allclose(out, torch.tensor([[7.0]], dtype=torch.float64), atol=1e-6)

def test_custom_function():
    x, y = sympy.symbols("x y")
    f = sympy.Function("f")
    z = x + f(y)
    extra_funcs = {f: lambda y_: y_**2}
    mod = torch_module([x, y], [z], extra_funcs=extra_funcs)
    result = mod(torch.tensor(1.0), torch.tensor(2.0))
    assert torch.isclose(result, torch.tensor([[5.0]]), atol=1e-6)  # 1 + 2^2 = 5

def test_rationals():
    xvals = np.random.randn(100)
    x = sympy.symbols("x")
    y = x * sympy.Rational(2, 7)
    mod = torch_module([x], [y])
    y_tilde = mod(torch.tensor(xvals, dtype=torch.float64))
    error = y_tilde.detach().numpy().flatten() - xvals * 2 / 7
    assert (error**2).mean() < 1e-10

def test_half1():
    xvals = np.random.randn(100)
    x = sympy.symbols("x")
    y = abs(x) ** sympy.S.Half
    mod = torch_module([x], [y])
    y_tilde = mod(torch.tensor(xvals, dtype=torch.float64))
    error = y_tilde.detach().numpy().flatten() - np.abs(xvals) ** 0.5
    assert (error**2).mean() < 1e-10

def test_half2():
    xvals = np.random.randn(100)
    y = sympy.parse_expr("sqrt(Abs(x))")
    mod = torch_module([x], [y])
    y_tilde = mod(torch.tensor(xvals, dtype=torch.float64))
    error = y_tilde.detach().numpy().flatten() - np.abs(xvals) ** 0.5
    assert (error**2).mean() < 1e-10

def test_constants2():
    constants = [
        sympy.pi,
        sympy.E,
        sympy.GoldenRatio,
        sympy.TribonacciConstant,
        sympy.EulerGamma,
        sympy.Catalan,
    ]
    # Explicitly set dtype to torch.float64
    mod = torch_module([], constants, dtype=torch.float32)
    out = mod()
    torch.testing.assert_close(out.flatten(), torch.tensor([float(c) for c in constants], dtype=torch.float32))

def test_complex():
    x = sympy.symbols("x")
    complex_func_torch = torch_module(
        [x],
        [
            x * sympy.I,
            sympy.conjugate(x),
            sympy.sqrt(sympy.conjugate(x * sympy.I) * x * sympy.I),
        ]
    )

    out = complex_func_torch(torch.tensor(2.0, dtype=torch.double)).detach().numpy()
    assert np.isclose(out[0].item(), 2.0j, atol=1e-6)
    assert np.isclose(out[1].item(), 2.0, atol=1e-6)
    assert np.isclose(out[2].item(), 2.0, atol

=1e-6)

    out = complex_func_torch(torch.tensor(2.0j, dtype=torch.complex128)).detach().numpy()
    assert np.isclose(out[0].item(), -2.0, atol=1e-6)
    assert np.isclose(out[1].item(), -2.0j, atol=1e-6)
    assert np.isclose(out[2].item(), 2.0, atol=1e-6)

    theta, phi = sympy.symbols("theta phi")
    max_l = 2
    m_list = range(-max_l, max_l + 1)
    func_list = [
        sympy.simplify(
            sympy.functions.special.spherical_harmonics.Znm(max_l, m, theta, phi).expand(func=True)
        ).evalf()
        for m in m_list
    ]
    func_list_np = [sympy.lambdify([theta, phi], func, "numpy") for func in func_list]
    func_list_torch = torch_module([theta, phi], func_list)

    np_eval = np.array([func_list_np[i](np.pi / 3, np.pi / 3) for i in range(5)])
    torch_eval = func_list_torch(torch.tensor(np.pi / 3, dtype=torch.double), torch.tensor(np.pi / 3, dtype=torch.double))
    error = np.sum(np.abs(torch_eval.detach().numpy() - np_eval))
    assert error < 1e-7

def test_integers():
    mod = torch_module([], [sympy.core.numbers.Zero()])
    assert torch.equal(mod(), torch.tensor(0.0))  # Adjust to scalar tensor

    mod = torch_module([], [sympy.core.numbers.One()])
    assert torch.equal(mod(), torch.tensor(1.0))  # Adjust to scalar tensor

    mod = torch_module([], [sympy.core.numbers.NegativeOne()])
    assert torch.equal(mod(), torch.tensor(-1.0))  # Adjust to scalar tensor

    for i in range(-10, 10):
        mod = torch_module([], [sympy.core.numbers.Integer(i)])
        assert torch.equal(mod(), torch.tensor(float(i)))


def test_torch_module_complex_system():
    x, y, z = symbols('x y z')
    A00, A01, A10, A11 = symbols('A00 A01 A10 A11')

    # Output will be a 3-vector: [f1, f2, f3]
    system = Matrix([
        sin(x) + cos(y) * A00 + sqrt(abs(z)),  # f1: Nonlinear mix of scalars and matrix element
        A01 * x ** 2 - A10 * y + A11 * z,  # f2: Polynomial with matrix elements
        (A00 + A11) * sin(x * y)  # f3: Diagonal trace-like term with nonlinearity
    ])

    mod = torch_module([x, y, z, A00, A01, A10, A11], system)

    x_val = torch.tensor(1.0, dtype=torch.float64)
    y_val = torch.tensor(2.0, dtype=torch.float64)
    z_val = torch.tensor(4.0, dtype=torch.float64)
    a_val = torch.tensor([1.0, -1.0, 0.5, 2.0], dtype=torch.float64)  # [A00, A01, A10, A11]

    result = mod(x_val, y_val, z_val, a_val[0], a_val[1], a_val[2], a_val[3])

    expected = torch.tensor([
        [0.8414709848078965 - 0.4161468365471424 + 2.0],
        [-1.0 - 1.0 + 8.0],
        [3.0 * 0.9092974268256817]
    ], dtype=torch.float64).squeeze()

    assert torch.allclose(result, expected, atol=1e-6), \
        f"Evaluation failed: Expected {expected}, got {result}"

    # Test case 2: Gradient computation
    x_val = torch.tensor(1.0, dtype=torch.float64, requires_grad=True)
    y_val = torch.tensor(2.0, dtype=torch.float64, requires_grad=True)
    z_val = torch.tensor(4.0, dtype=torch.float64, requires_grad=True)
    a_val = torch.tensor([1.0, -1.0, 0.5, 2.0], dtype=torch.float64)

    result = mod(x_val, y_val, z_val, a_val[0], a_val[1], a_val[2], a_val[3])
    assert result.requires_grad, "Output should support gradients"

    result.sum().backward()
    assert x_val.grad is not None, "Gradient w.r.t. x should be computed"
    assert y_val.grad is not None, "Gradient w.r.t. y should be computed"
    assert z_val.grad is not None, "Gradient w.r.t. z should be computed"

    expected_x_grad = (
        cos(1.0) +
        2 * (-1.0) * 1.0 +
        (1.0 + 2.0) * cos(1.0 * 2.0) * 2.0
    )
    expected_y_grad = (
        -sin(2.0) * 1.0 +
        -0.5 +
        (1.0 + 2.0) * cos(1.0 * 2.0) * 1.0
    )
    expected_z_grad = (
        1.0 / (2 * sqrt(4.0)) * 1.0 +
        2.0 +
        0.0
    )

    expected_grad = torch.tensor([expected_x_grad, expected_y_grad, expected_z_grad], dtype=torch.float64)
    computed_grad = torch.tensor([x_val.grad.item(), y_val.grad.item(), z_val.grad.item()], dtype=torch.float64)

    assert torch.allclose(computed_grad, expected_grad, atol=1e-6), \
        f"Gradient failed: Expected {expected_grad}, got {computed_grad}"

    x_val = torch.tensor(0.5, dtype=torch.float64)
    y_val = torch.tensor(1.0, dtype=torch.float64)
    z_val = torch.tensor(-9.0, dtype=torch.float64)  # Negative z to test abs
    a_val = torch.tensor([2.0, 0.5, -1.0, 3.0], dtype=torch.float64)

    result = mod(x_val, y_val, z_val, a_val[0], a_val[1], a_val[2], a_val[3])

    expected = torch.tensor([
        [0.479425538604203 + 0.5403023058681398 * 2 + 3.0],
        [0.5 * 0.25 + 1.0 - 27.0],
        [5.0 * 0.479425538604203]
    ], dtype=torch.float64).squeeze()

    assert torch.allclose(result, expected, atol=1e-6), \
        f"Different input failed: Expected {expected}, got {result}"
