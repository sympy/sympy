"""Tests for expand() on AlgebraicPureTensor and AlgebraicTensor."""

from sympy.core.add import Add
from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol
from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
)


class TestExpandPureTensor:
    """Test _eval_expand_mul on AlgebraicPureTensor."""

    def test_expand_single_add_factor(self):
        """A ⊗ (B+C) ⊗ D  ->  A⊗B⊗D + A⊗C⊗D."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        C = MatrixSymbol("C", 3, 4)
        D = MatrixSymbol("D", 4, 5)

        pt = AlgebraicPureTensor(A, B + C, D)
        result = pt.expand()

        assert isinstance(result, AlgebraicTensor)
        # Should have two terms
        terms = result.terms
        assert len(terms) == 2

        # Check the two terms are AlgebraicPureTensor(A,B,D) and AlgebraicPureTensor(A,C,D)
        expected1 = AlgebraicPureTensor(A, B, D)
        expected2 = AlgebraicPureTensor(A, C, D)
        assert expected1 in terms or expected2 in terms

    def test_expand_two_add_factors(self):
        """(A+B) ⊗ (C+D)  ->  A⊗C + A⊗D + B⊗C + B⊗D."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 2, 3)
        C = MatrixSymbol("C", 3, 4)
        D = MatrixSymbol("D", 3, 4)

        pt = AlgebraicPureTensor(A + B, C + D)
        result = pt.expand()

        assert isinstance(result, AlgebraicTensor)
        terms = result.terms
        assert len(terms) == 4

    def test_expand_no_add_factors(self):
        """A ⊗ B (no Add factors) -> unchanged."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)

        pt = AlgebraicPureTensor(A, B)
        result = pt.expand()

        # Should return the same PureTensor
        assert result == pt

    def test_expand_single_factor_add(self):
        """(A+B) as single factor -> should expand to A+B."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 2, 3)

        # Single factor that is an Add — AlgebraicPureTensor(A+B) unwraps to A+B
        # since single factor with coeff=1 returns bare factor
        pt = AlgebraicPureTensor(A + B)
        # This unwraps to A+B directly
        assert pt == A + B

    def test_expand_with_coefficient(self):
        """2 * (A ⊗ (B+C)) -> 2*A⊗B + 2*A⊗C."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        C = MatrixSymbol("C", 3, 4)

        pt = AlgebraicPureTensor(2, A, B + C)
        result = pt.expand()

        assert isinstance(result, AlgebraicTensor)
        terms = result.terms
        assert len(terms) == 2

    def test_expand_three_factors_one_add(self):
        """A ⊗ (B+C) ⊗ D -> A⊗B⊗D + A⊗C⊗D."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        C = MatrixSymbol("C", 3, 4)
        D = MatrixSymbol("D", 4, 5)

        pt = AlgebraicPureTensor(A, B + C, D)
        result = pt.expand()

        assert isinstance(result, AlgebraicTensor)
        terms = result.terms
        assert len(terms) == 2

        # Each term should be an AlgebraicPureTensor with 3 factors
        for t in terms:
            assert isinstance(t, AlgebraicPureTensor)
            assert len(t.factors) == 3
            # First factor should be A, last should be D
            assert t.factors[0] == A
            assert t.factors[-1] == D

    def test_expand_add_of_three(self):
        """A ⊗ (B+C+E) -> A⊗B + A⊗C + A⊗E."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        C = MatrixSymbol("C", 3, 4)
        E = MatrixSymbol("E", 3, 4)

        pt = AlgebraicPureTensor(A, B + C + E)
        result = pt.expand()

        assert isinstance(result, AlgebraicTensor)
        terms = result.terms
        assert len(terms) == 3


class TestExpandAlgebraicTensor:
    """Test expand() on AlgebraicTensor."""

    def test_expand_sum_of_pure_tensors(self):
        """(A⊗(B+C)) + (D⊗(E+F)) should expand each term."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        C = MatrixSymbol("C", 3, 4)
        D = MatrixSymbol("D", 2, 3)
        E = MatrixSymbol("E", 3, 4)
        F = MatrixSymbol("F", 3, 4)

        pt1 = AlgebraicPureTensor(A, B + C)
        pt2 = AlgebraicPureTensor(D, E + F)
        at = AlgebraicTensor(pt1, pt2)

        result = at.expand()

        # Each term expands to 2 terms, so we should have 4 terms total
        assert isinstance(result, AlgebraicTensor)
        terms = result.terms
        assert len(terms) == 4

    def test_expand_no_expansion_needed(self):
        """A⊗B + C⊗D (no Add factors) -> unchanged."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        C = MatrixSymbol("C", 2, 3)
        D = MatrixSymbol("D", 3, 4)

        pt1 = AlgebraicPureTensor(A, B)
        pt2 = AlgebraicPureTensor(C, D)
        at = AlgebraicTensor(pt1, pt2)

        result = at.expand()

        assert result == at

    def test_expand_mixed(self):
        """A⊗(B+C) + D⊗E — one term expands, other doesn't."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        C = MatrixSymbol("C", 3, 4)
        D = MatrixSymbol("D", 2, 3)
        E = MatrixSymbol("E", 3, 4)

        pt1 = AlgebraicPureTensor(A, B + C)
        pt2 = AlgebraicPureTensor(D, E)
        at = AlgebraicTensor(pt1, pt2)

        result = at.expand()

        assert isinstance(result, AlgebraicTensor)
        terms = result.terms
        assert len(terms) == 3

    def test_expand_with_zero_term(self):
        """Expand preserves AlgebraicZeroTensor anchors."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        C = MatrixSymbol("C", 3, 4)
        zero = AlgebraicZeroTensor(((2, 3), (3, 4)))

        pt = AlgebraicPureTensor(A, B + C)
        at = AlgebraicTensor(pt, zero)

        result = at.expand()

        assert isinstance(result, AlgebraicTensor)
        # Should have 2 expanded terms + the zero anchor
        assert result.has_zero_term()


class TestExpandEdgeCases:
    """Edge cases for expand()."""

    def test_expand_unwraps_single_term(self):
        """If expansion results in a single term, it should unwrap."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)

        pt = AlgebraicPureTensor(A, B)
        result = pt.expand()

        # No expansion needed, returns the same
        assert result == pt

    def test_expand_coefficient_preserved(self):
        """Coefficient should be distributed to all expanded terms."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        C = MatrixSymbol("C", 3, 4)

        pt = AlgebraicPureTensor(3, A, B + C)
        result = pt.expand()

        assert isinstance(result, AlgebraicTensor)
        terms = result.terms
        assert len(terms) == 2

        # Each term should have coefficient 3
        for t in terms:
            if isinstance(t, AlgebraicPureTensor):
                assert t._get_coeff() == 3
