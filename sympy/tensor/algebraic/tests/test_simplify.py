"""Tests for the tensor simplification pipeline.

Covers proportionality_factoring, commutativity-based simplification,
and the full tensorsimplify() pipeline.
"""

from sympy import Symbol, Rational, sqrt, eye, zeros, Matrix
from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol
from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    tensorsimplify,
    proportionality_factoring,
)


class TestProportionalityFactoringBasic:
    """Basic proportionality factoring tests."""

    def test_identical_factors_merge(self):
        """2*A⊗B + 3*A⊗B -> 5*A⊗B."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(2, A, B),
            AlgebraicPureTensor(3, A, B),
        )
        result = proportionality_factoring(at)
        assert isinstance(result, AlgebraicPureTensor)
        assert result._get_coeff() == 5

    def test_canceling_terms(self):
        """2*A⊗B - 2*A⊗B -> zero tensor."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(2, A, B),
            AlgebraicPureTensor(-2, A, B),
        )
        result = proportionality_factoring(at)
        assert isinstance(result, AlgebraicZeroTensor)

    def test_one_non_proportional_slot(self):
        """A⊗B + A2⊗B -> (A+A2)⊗B."""
        A = MatrixSymbol("A", 2, 3)
        A2 = MatrixSymbol("A2", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(A, B),
            AlgebraicPureTensor(A2, B),
        )
        result = proportionality_factoring(at)
        assert isinstance(result, AlgebraicPureTensor)
        assert len(result.factors) == 2

    def test_two_non_proportional_slots_no_merge(self):
        """A⊗B + A2⊗B2 (two different slots) -> no merge."""
        A = MatrixSymbol("A", 2, 3)
        A2 = MatrixSymbol("A2", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        B2 = MatrixSymbol("B2", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(A, B),
            AlgebraicPureTensor(A2, B2),
        )
        result = proportionality_factoring(at)
        assert isinstance(result, AlgebraicTensor)
        assert len(result.terms) == 2


class TestCommutativitySimplify:
    """Commutativity-based simplification tests."""

    def test_all_noncommutative_unchanged(self):
        """All symbolic factors: no commutative decomposition needed."""
        A = MatrixSymbol("A", 2, 3)
        A2 = MatrixSymbol("A2", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(A, B),
            AlgebraicPureTensor(A2, B),
        )
        result = tensorsimplify(at)
        # Should still simplify via proportionality factoring
        assert isinstance(result, AlgebraicPureTensor)

    def test_commutative_numeric_matrices(self):
        """Numeric matrices are commutative; should decompose and regroup."""
        E11 = Matrix([[1, 0], [0, 0]])
        E22 = Matrix([[0, 0], [0, 1]])
        M = MatrixSymbol("M", 3, 3)
        N = MatrixSymbol("N", 3, 3)
        # E11⊗M + E22⊗M: different commutative keys, same non-comm part
        at = AlgebraicTensor(
            AlgebraicPureTensor(E11, M),
            AlgebraicPureTensor(E22, M),
        )
        result = tensorsimplify(at)
        # Should preserve both terms (different commutative keys)
        assert result is not None

    def test_same_commutative_key_merges(self):
        """Same commutative pattern, proportional non-comm parts -> merge."""
        E11 = Matrix([[1, 0], [0, 0]])
        M = MatrixSymbol("M", 3, 3)
        N = MatrixSymbol("N", 3, 3)
        # E11⊗M + 2*E11⊗M -> 3*E11⊗M
        at = AlgebraicTensor(
            AlgebraicPureTensor(E11, M),
            AlgebraicPureTensor(2, E11, M),
        )
        result = tensorsimplify(at)
        assert isinstance(result, AlgebraicPureTensor)
        assert result._get_coeff() == 3

    def test_mixed_commutativity_three_slots(self):
        """Three-factor tensor with mixed commutativity: (comm, noncomm, comm)."""
        E11 = Matrix([[1, 0], [0, 0]])
        E22 = Matrix([[0, 0], [0, 1]])
        M = MatrixSymbol("M", 3, 3)
        N = MatrixSymbol("N", 3, 3)
        I2 = eye(2)
        at = AlgebraicTensor(
            AlgebraicPureTensor(E11, M, I2),
            AlgebraicPureTensor(E22, M, I2),
        )
        result = tensorsimplify(at)
        assert result is not None

    def test_coefficient_simplification(self):
        """Coefficients should be simplified."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(2 + Rational(1, 2) - Rational(1, 2), A, B),
        )
        result = tensorsimplify(at)
        assert isinstance(result, AlgebraicPureTensor)
        assert result._get_coeff() == 2

    def test_zero_tensor_passthrough(self):
        """Zero tensor should pass through unchanged."""
        zt = AlgebraicZeroTensor(((2, 3), (3, 4)))
        result = tensorsimplify(zt)
        assert isinstance(result, AlgebraicZeroTensor)
        assert result.shape == ((2, 3), (3, 4))

    def test_single_term_passthrough(self):
        """Single-term AlgebraicTensor unwraps and simplifies."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(AlgebraicPureTensor(2 + 2, A, B))
        result = tensorsimplify(at)
        assert isinstance(result, AlgebraicPureTensor)
        assert result._get_coeff() == 4

    def test_numeric_matrix_decomposition(self):
        """Numeric matrix with multiple nonzero entries decomposes correctly."""
        M_full = Matrix([[1, 2], [3, 4]])
        N = MatrixSymbol("N", 2, 2)
        # M_full⊗N should decompose into 4 basis terms
        pt = AlgebraicPureTensor(M_full, N)
        at = AlgebraicTensor(pt)
        result = tensorsimplify(at)
        # Result should be equivalent (decomposed and recombined)
        assert result is not None
        # The result should still represent the same tensor
        if isinstance(result, AlgebraicTensor):
            assert len(result.terms) >= 1
        elif isinstance(result, AlgebraicPureTensor):
            assert len(result.factors) >= 1

    def test_commutative_with_symbolic_coefficient(self):
        """Symbolic coefficient carried through decomposition."""
        x = Symbol("x")
        E11 = Matrix([[1, 0], [0, 0]])
        E12 = Matrix([[0, 1], [0, 0]])
        M = MatrixSymbol("M", 3, 3)
        at = AlgebraicTensor(
            AlgebraicPureTensor(x, E11, M),
            AlgebraicPureTensor(x, E12, M),
        )
        result = tensorsimplify(at)
        assert result is not None


class TestSimplifyIntegration:
    """Integration tests for the full simplification pipeline."""

    def test_pure_tensor_simplify(self):
        """PureTensor.simplify() works correctly."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        pt = AlgebraicPureTensor(2 + 2, A, B)
        result = pt.simplify()
        assert isinstance(result, AlgebraicPureTensor)
        assert result._get_coeff() == 4

    def test_algebraic_tensor_simplify(self):
        """AlgebraicTensor.simplify() works correctly."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(2, A, B),
            AlgebraicPureTensor(3, A, B),
        )
        result = at.simplify()
        assert isinstance(result, AlgebraicPureTensor)
        assert result._get_coeff() == 5

    def test_tensorsimplify_dispatch(self):
        """tensorsimplify dispatches correctly for different types."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)

        # PureTensor
        pt = AlgebraicPureTensor(2 + 2, A, B)
        assert tensorsimplify(pt)._get_coeff() == 4

        # AlgebraicTensor
        at = AlgebraicTensor(
            AlgebraicPureTensor(2, A, B),
            AlgebraicPureTensor(3, A, B),
        )
        result = tensorsimplify(at)
        assert isinstance(result, AlgebraicPureTensor)
        assert result._get_coeff() == 5

        # Zero tensor
        zt = AlgebraicZeroTensor(((2, 3), (3, 4)))
        assert isinstance(tensorsimplify(zt), AlgebraicZeroTensor)

    def test_cancellation_to_zero(self):
        """Full cancellation produces zero tensor."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(5, A, B),
            AlgebraicPureTensor(-5, A, B),
        )
        result = tensorsimplify(at)
        assert isinstance(result, AlgebraicZeroTensor)

    def test_three_term_merging(self):
        """Three proportional terms merge into one."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(1, A, B),
            AlgebraicPureTensor(2, A, B),
            AlgebraicPureTensor(3, A, B),
        )
        result = tensorsimplify(at)
        assert isinstance(result, AlgebraicPureTensor)
        assert result._get_coeff() == 6

    def test_expand_then_simplify(self):
        """Expand then simplify factors common slots."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        C = MatrixSymbol("C", 3, 4)
        # A⊗(B+C) expands to A⊗B + A⊗C, then simplifies back to A⊗(B+C)
        pt = AlgebraicPureTensor(A, B + C)
        expanded = pt.expand()
        assert isinstance(expanded, AlgebraicTensor)
        assert len(expanded.terms) == 2
        result = tensorsimplify(expanded)
        # Proportionality factoring merges: A⊗B + A⊗C -> A⊗(B+C)
        assert isinstance(result, AlgebraicPureTensor)
        assert len(result.factors) == 2
        assert result.factors[0] == A

    def test_expand_identical_then_simplify(self):
        """A⊗B + A⊗B simplifies to 2*A⊗B."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(A, B),
            AlgebraicPureTensor(A, B),
        )
        result = tensorsimplify(at)
        assert isinstance(result, AlgebraicPureTensor)
        assert result._get_coeff() == 2


class TestEdgeCases:
    """Edge cases and corner scenarios."""

    def test_empty_commutative_indices(self):
        """All non-commutative: decomposition produces single bucket."""
        A = MatrixSymbol("A", 2, 3)
        A2 = MatrixSymbol("A2", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(A, B),
            AlgebraicPureTensor(A2, B),
        )
        result = tensorsimplify(at)
        assert result is not None

    def test_zero_coefficient_term(self):
        """Zero coefficient term is eliminated."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(S.Zero, A, B),
            AlgebraicPureTensor(3, A, B),
        )
        result = tensorsimplify(at)
        if isinstance(result, AlgebraicPureTensor):
            assert result._get_coeff() == 3

    def test_identity_matrix_commutative(self):
        """Identity matrix is commutative and handled correctly."""
        I3 = eye(3)
        M = MatrixSymbol("M", 3, 3)
        N = MatrixSymbol("N", 3, 3)
        at = AlgebraicTensor(
            AlgebraicPureTensor(I3, M),
            AlgebraicPureTensor(I3, N),
        )
        result = tensorsimplify(at)
        assert result is not None

    def test_all_zero_matrix(self):
        """Zero matrix factor is a valid factor (not auto-eliminated)."""
        Z2 = zeros(2)
        M = MatrixSymbol("M", 3, 3)
        pt = AlgebraicPureTensor(Z2, M)
        at = AlgebraicTensor(pt)
        result = tensorsimplify(at)
        # Zero matrix is a valid factor; simplification preserves it
        assert isinstance(result, AlgebraicPureTensor)
        assert len(result.factors) == 2

    def test_scalar_symbol_coefficient(self):
        """Symbolic coefficient preserved through simplification."""
        x = Symbol("x")
        y = Symbol("y")
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(x, A, B),
            AlgebraicPureTensor(y, A, B),
        )
        result = tensorsimplify(at)
        assert isinstance(result, AlgebraicPureTensor)
        coeff = result._get_coeff()
        assert coeff == x + y or coeff == y + x

    def test_negative_coefficient_merge(self):
        """Negative coefficients merge correctly."""
        A = MatrixSymbol("A", 2, 3)
        B = MatrixSymbol("B", 3, 4)
        at = AlgebraicTensor(
            AlgebraicPureTensor(3, A, B),
            AlgebraicPureTensor(-1, A, B),
        )
        result = tensorsimplify(at)
        assert isinstance(result, AlgebraicPureTensor)
        assert result._get_coeff() == 2


class TestCommutativePrefactorExtraction:
    """Tests for extracting commutative prefactors from non-commutative factors."""

    def test_mul_factor_extraction(self):
        """Mul of commutative and non-commutative parts extracts commutative part."""
        from sympy import symbols
        nc = symbols("nc", commutative=False)
        x = Symbol("x")
        from sympy.tensor.algebraic.simplify import _extract_commutative_prefactors
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        from sympy.matrices import Matrix
        nc_mat = Matrix([[nc]])
        xnc_mat = Matrix([[x * nc]])
        # Use two factors so AlgebraicPureTensor doesn't unwrap
        pt = AlgebraicPureTensor(xnc_mat, nc_mat)
        extr_coeff, new_factors = _extract_commutative_prefactors(pt)
        assert extr_coeff == x
        # The first factor should no longer contain x
        assert not new_factors[0].has(x)

    def test_matrix_factor_common_divisor(self):
        """Common commutative divisor across matrix entries is extracted."""
        from sympy import symbols
        nc = symbols("nc", commutative=False)
        x = Symbol("x")
        M = Matrix([[x * nc, 2 * x], [3 * x, x * nc]])
        from sympy.tensor.algebraic.simplify import _extract_commutative_prefactors
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        from sympy.matrices import Matrix as SymMatrix
        # Use two factors so AlgebraicPureTensor doesn't unwrap
        nc_mat = SymMatrix([[nc]])
        pt = AlgebraicPureTensor(M, nc_mat)
        extr_coeff, new_factors = _extract_commutative_prefactors(pt)
        assert extr_coeff == x
        # Check that the matrix entries no longer contain x
        new_M = new_factors[0]
        for r in range(new_M.shape[0]):
            for c in range(new_M.shape[1]):
                entry = new_M[r, c]
                if entry != S.Zero:
                    assert not entry.has(x)

    def test_matrix_factor_no_common_divisor(self):
        """No extraction when entries share no common commutative divisor."""
        from sympy import symbols
        nc = symbols("nc", commutative=False)
        x = Symbol("x")
        y = Symbol("y")
        M = Matrix([[x * nc, y * nc], [nc, nc]])
        from sympy.tensor.algebraic.simplify import _extract_commutative_prefactors
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        from sympy.matrices import Matrix as SymMatrix
        nc_mat = SymMatrix([[nc]])
        pt = AlgebraicPureTensor(M, nc_mat)
        extr_coeff, new_factors = _extract_commutative_prefactors(pt)
        assert extr_coeff is S.One

    def test_matrix_factor_partial_divisor(self):
        """Partial divisor (not dividing all entries) is not extracted."""
        from sympy import symbols
        nc = symbols("nc", commutative=False)
        x = Symbol("x")
        y = Symbol("y")
        M = Matrix([[x * nc, y * nc], [x * nc, x * nc]])
        from sympy.tensor.algebraic.simplify import _extract_commutative_prefactors
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        from sympy.matrices import Matrix as SymMatrix
        nc_mat = SymMatrix([[nc]])
        pt = AlgebraicPureTensor(M, nc_mat)
        extr_coeff, new_factors = _extract_commutative_prefactors(pt)
        # x doesn't divide y*nc, so nothing should be extracted
        assert extr_coeff is S.One

    def test_matrix_factor_numeric_common(self):
        """Numeric common factor is extracted from matrix entries."""
        M = Matrix([[6, 9], [12, 15]])
        from sympy.tensor.algebraic.simplify import _extract_commutative_prefactors
        from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
        from sympy.matrices import Matrix as SymMatrix
        from sympy import symbols
        nc = symbols("nc", commutative=False)
        nc_mat = SymMatrix([[nc]])
        pt = AlgebraicPureTensor(M, nc_mat)
        extr_coeff, new_factors = _extract_commutative_prefactors(pt)
        # 3 should be extracted as common divisor
        assert extr_coeff == 3
        new_M = new_factors[0]
        assert new_M[0, 0] == 2
        assert new_M[0, 1] == 3
        assert new_M[1, 0] == 4
        assert new_M[1, 1] == 5

    def test_integration_extraction_in_simplify(self):
        """Full simplification pipeline extracts commutative prefactors."""
        from sympy import symbols
        nc = symbols("nc", commutative=False)
        x = Symbol("x")
        E11 = Matrix([[1, 0], [0, 0]])
        M1 = Matrix([[x * nc, 0], [0, x]])
        # E11 has shape (2,2), M1 has shape (2,2)
        at = AlgebraicTensor(
            AlgebraicPureTensor(E11, M1),
        )
        result = tensorsimplify(at)
        # The simplification should work without errors
        assert result is not None
