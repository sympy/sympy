from __future__ import annotations

from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol
from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
from sympy.tensor.algebraic.pure_tensor import PureTensor


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
I3 = MatrixSymbol("I3", 3, 3)


def test_pure_tensor_str_with_coefficient():
    pt = PureTensor(A, C)
    scaled = 3 * pt
    s = str(scaled)
    assert "3" in s
    assert "A" in s
    assert "C" in s
    assert "\u2297" in s


def test_pure_tensor_str_with_negative_coefficient():
    pt = PureTensor(A, C)
    neg = -3 * pt
    s = str(neg)
    assert s.startswith("-")
    assert "3" in s
    assert "A" in s
    assert "C" in s


def test_pure_tensor_str_with_negative_one():
    pt = PureTensor(A, C)
    neg = -pt
    s = str(neg)
    assert s.startswith("-")
    assert "A" in s
    assert "C" in s
    assert "-1" not in s


def test_pure_tensor_repr_with_coefficient():
    pt = PureTensor(A, C)
    scaled = 3 * pt
    r = repr(scaled)
    assert "3" in r
    assert "PureTensor" in r


def test_pure_tensor_repr_without_coefficient():
    pt = PureTensor(A, C)
    r = repr(pt)
    assert "PureTensor(A, C)" in r or "PureTensor(A,C)" in r


def test_pure_tensor_str_coefficient_one_omitted():
    pt = PureTensor(A, C)
    s = str(pt)
    assert not s.startswith("1")
    assert "*A" not in s.split("\u2297")[0] or s.startswith("A")


def test_pure_tensor_double_scaled_str():
    pt = PureTensor(A, C)
    result = (2 * pt) * 3
    s = str(result)
    assert "6" in s
    assert isinstance(result, PureTensor)


def test_algebraic_tensor_str_with_coefficients():
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    at = AlgebraicTensor(2 * pt_a, 3 * pt_b)
    s = str(at)
    assert "2" in s
    assert "3" in s
    assert "+" in s


def test_algebraic_tensor_str_with_negative_term():
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    at = AlgebraicTensor(2 * pt_a, -3 * pt_b)
    s = str(at)
    assert "-" in s
    assert "+ -" not in s


def test_linear_combination_str_display():
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    combo = 2 * pt_a + 3 * pt_b
    s = str(combo)
    assert isinstance(combo, AlgebraicTensor)
    assert "2" in s
    assert "3" in s
    assert "+" in s


def test_linear_combination_subtraction_str():
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    combo = 2 * pt_a - 3 * pt_b
    s = str(combo)
    assert isinstance(combo, AlgebraicTensor)
    assert "2" in s
    assert "3" in s
    assert "-" in s
    assert "+ -" not in s


def test_linear_combination_three_terms_str():
    D2 = MatrixSymbol("D2", 3, 4)
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    pt_d = PureTensor(D2, C)
    combo = pt_a + 2 * pt_b - 3 * pt_d
    s = str(combo)
    assert isinstance(combo, AlgebraicTensor)
    assert "2" in s
    assert "3" in s
    assert "-" in s


def test_linear_combination_symbol_coefficient_str():
    from sympy import Symbol
    x = Symbol("x")
    y = Symbol("y")
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    scaled_a = pt_a * x
    scaled_b = pt_b * y
    assert isinstance(scaled_a, PureTensor)
    assert isinstance(scaled_b, PureTensor)
    combo = AlgebraicTensor(scaled_a, scaled_b)
    s = str(combo)
    assert isinstance(combo, AlgebraicTensor)
    assert "x" in s
    assert "y" in s
    assert "+" in s


def test_linear_combination_repr():
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    combo = 2 * pt_a + 3 * pt_b
    r = repr(combo)
    assert "AlgebraicTensor" in r
    assert "2" in r
    assert "3" in r


def test_algebraic_tensor_mul_scalar_str():
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    at = pt_a + pt_b
    scaled = at * 5
    s = str(scaled)
    assert isinstance(scaled, AlgebraicTensor)
    assert "5" in s


def test_algebraic_tensor_rmul_scalar_str():
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    at = pt_a + pt_b
    scaled = 5 * at
    s = str(scaled)
    assert isinstance(scaled, AlgebraicTensor)
    assert "5" in s


def test_pure_tensor_coeff_preserved_through_add():
    pt = PureTensor(A, C)
    scaled = 3 * pt
    at = AlgebraicTensor(scaled, PureTensor(B, C))
    for arg in at.args:
        if isinstance(arg, PureTensor) and arg.factors == (A, C):
            assert arg._get_coeff() == S(3)


def test_linear_combination_with_zerotensor_str():
    from sympy.tensor.algebraic.zero_tensor import ZeroTensor
    z = ZeroTensor(((3, 4), (4, 5)))
    pt = PureTensor(A, C)
    at = AlgebraicTensor(2 * pt, z)
    s = str(at)
    assert "2" in s
    assert "0_" in s
