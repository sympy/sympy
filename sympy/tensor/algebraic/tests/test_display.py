"""Tests for LaTeX display of algebraic tensors via _repr_latex_()."""

from sympy.core.add import Add
from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol
from sympy.printing.latex import latex
from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
)
from sympy import Symbol


class TestZeroTensorDisplay:
    """Test LaTeX display of AlgebraicZeroTensor."""

    def test_zero_tensor_two_factors(self):
        z = AlgebraicZeroTensor(((3, 4), (4, 5)))
        s = latex(z, mode="plain")
        assert "0_" in s
        assert "3 \\times 4" in s
        assert "4 \\times 5" in s
        assert "\\otimes" in s

    def test_zero_tensor_single_factor(self):
        z = AlgebraicZeroTensor(((3, 4),))
        s = latex(z, mode="plain")
        assert "0_" in s
        assert "3 \\times 4" in s
        # Single factor should not have \otimes
        assert "\\otimes" not in s

    def test_repr_latex(self):
        z = AlgebraicZeroTensor(((2, 3),))
        r = z._repr_latex_()
        assert r.startswith("$\\displaystyle ")
        assert "0_" in r
        assert "2 \\times 3" in r


class TestPureTensorDisplay:
    """Test LaTeX display of AlgebraicPureTensor."""

    def test_basic_tensor_product(self):
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        pt = AlgebraicPureTensor(A, B)
        s = latex(pt, mode="plain")
        assert "A \\otimes B" == s

    def test_three_factors(self):
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        C = MatrixSymbol("C", 4, 5)
        pt = AlgebraicPureTensor(A, B, C)
        s = latex(pt, mode="plain")
        assert "A \\otimes B \\otimes C" == s

    def test_positive_coefficient(self):
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        pt = AlgebraicPureTensor(2, A, B)
        s = latex(pt, mode="plain")
        assert "2 A \\otimes B" == s

    def test_negative_one_coefficient(self):
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        pt = AlgebraicPureTensor(S.NegativeOne, A, B)
        s = latex(pt, mode="plain")
        assert "-A \\otimes B" == s

    def test_negative_coefficient(self):
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        pt = AlgebraicPureTensor(-3, A, B)
        s = latex(pt, mode="plain")
        assert "-3 A \\otimes B" == s

    def test_repr_latex(self):
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        pt = AlgebraicPureTensor(A, B)
        r = pt._repr_latex_()
        assert r.startswith("$\\displaystyle ")
        assert "A \\otimes B" in r
        assert r.endswith("$")


class TestAlgebraicTensorDisplay:
    """Test LaTeX display of AlgebraicTensor."""

    def test_sum_two_terms(self):
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        C = MatrixSymbol("C", 3, 4)
        D = MatrixSymbol("D", 4, 5)
        at = AlgebraicTensor(
            AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D)
        )
        s = latex(at, mode="plain")
        # Both terms should be present, separated by +
        assert "A \\otimes B" in s
        assert "C \\otimes D" in s
        assert "+" in s

    def test_subtraction(self):
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        C = MatrixSymbol("C", 3, 4)
        D = MatrixSymbol("D", 4, 5)
        at = AlgebraicTensor(
            AlgebraicPureTensor(A, B),
            AlgebraicPureTensor(S.NegativeOne, C, D),
        )
        s = latex(at, mode="plain")
        assert "A \\otimes B" in s
        assert "C \\otimes D" in s
        assert "-" in s
        # Should not have double minus
        assert "--" not in s

    def test_mixed_coefficients(self):
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        C = MatrixSymbol("C", 3, 4)
        D = MatrixSymbol("D", 4, 5)
        at = AlgebraicTensor(
            AlgebraicPureTensor(2, A, B),
            AlgebraicPureTensor(S.NegativeOne, C, D),
        )
        s = latex(at, mode="plain")
        assert "2 A \\otimes B" in s
        assert "C \\otimes D" in s
        assert "-" in s

    def test_negated_term_first(self):
        """When the negative term is ordered first."""
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        C = MatrixSymbol("C", 3, 4)
        D = MatrixSymbol("D", 4, 5)
        at = AlgebraicTensor(
            AlgebraicPureTensor(S.NegativeOne, C, D),
            AlgebraicPureTensor(2, A, B),
        )
        s = latex(at, mode="plain")
        assert "--" not in s
        assert "2 A \\otimes B" in s

    def test_repr_latex(self):
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        C = MatrixSymbol("C", 3, 4)
        D = MatrixSymbol("D", 4, 5)
        at = AlgebraicTensor(
            AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D)
        )
        r = at._repr_latex_()
        assert r.startswith("$\\displaystyle ")
        assert "\\otimes" in r
        assert r.endswith("$")


class TestAddPrefactorBracketing:
    """Test that Add prefactors are wrapped in brackets."""

    def test_add_coeff_pure_tensor(self):
        r"""An Add coefficient gets \left...\right) brackets."""
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        a = Symbol("a")
        b = Symbol("b")
        pt = AlgebraicPureTensor(a + b, A, B)
        s = latex(pt, mode="plain")
        assert "\\left(" in s
        assert "\\right)" in s
        assert "a + b" in s

    def test_mul_coeff_no_brackets(self):
        """A simple Mul coefficient does not get brackets."""
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        a = Symbol("a")
        b = Symbol("b")
        pt = AlgebraicPureTensor(a * b, A, B)
        s = latex(pt, mode="plain")
        assert "\\left(" not in s

    def test_add_coeff_algebraic_tensor(self):
        """Add coefficients in AlgebraicTensor get brackets."""
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        a = Symbol("a")
        b = Symbol("b")
        at = AlgebraicTensor(
            AlgebraicPureTensor(a + b, A, B),
            AlgebraicPureTensor(A, B),
        )
        s = latex(at, mode="plain")
        assert "\\left(" in s
        assert "\\right)" in s

    def test_display_method_exists(self):
        """All tensor types have a display() method."""
        A = MatrixSymbol("A", 3, 4)
        B = MatrixSymbol("B", 4, 5)
        pt = AlgebraicPureTensor(A, B)
        z = AlgebraicZeroTensor(((3, 4), (4, 5)))
        at = AlgebraicTensor(AlgebraicPureTensor(A, B))

        assert hasattr(pt, "display")
        assert hasattr(z, "display")
        assert hasattr(at, "display")
        assert callable(pt.display)
        assert callable(z.display)
        assert callable(at.display)
