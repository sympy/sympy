from __future__ import annotations

import pickle

from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.abc import x, y
from sympy.matrices.expressions import MatrixSymbol
from sympy.testing.pytest import raises

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicZeroTensor,
    algebraic_zero_tensor,
)


A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 4, 5)
C = MatrixSymbol("C", 3, 4)
D = MatrixSymbol("D", 4, 5)
E = MatrixSymbol("E", 5, 3)


# ---------------------------------------------------------------------------
# Constructor and shape normalization
# ---------------------------------------------------------------------------

def test_constructor_bare_pair():
    """Bare (m, n) pair is wrapped as ((m, n),)."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.shape == ((3, 4),)


def test_constructor_list_pair():
    """List [m, n] is wrapped as ((m, n),)."""
    Z = AlgebraicZeroTensor([3, 4])
    assert Z.shape == ((3, 4),)


def test_constructor_single_factor_tuple():
    """Single-factor tuple ((m, n),) is preserved."""
    Z = AlgebraicZeroTensor(((3, 4),))
    assert Z.shape == ((3, 4),)


def test_constructor_multi_factor():
    """Multi-factor shape is preserved."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z.shape == ((3, 4), (4, 5))


def test_constructor_list_of_lists():
    """List of lists is converted to tuple of tuples."""
    Z = AlgebraicZeroTensor([[3, 4], [4, 5]])
    assert Z.shape == ((3, 4), (4, 5))


def test_constructor_three_factors():
    """Three-factor shape is handled correctly."""
    Z = AlgebraicZeroTensor(((2, 3), (3, 4), (4, 5)))
    assert Z.shape == ((2, 3), (3, 4), (4, 5))


def test_constructor_list_of_lists_three_factors():
    """List of three lists is converted to tuple of tuples."""
    Z = AlgebraicZeroTensor([[1, 2], [2, 3], [4, 5]])
    assert Z.shape == ((1, 2), (2, 3), (4, 5))


def test_convenience_constructor():
    """algebraic_zero_tensor convenience function works."""
    Z1 = algebraic_zero_tensor((3, 4))
    assert Z1.shape == ((3, 4),)
    Z2 = algebraic_zero_tensor(((3, 4), (4, 5)))
    assert Z2.shape == ((3, 4), (4, 5))
    assert isinstance(Z1, AlgebraicZeroTensor)
    assert isinstance(Z2, AlgebraicZeroTensor)


# ---------------------------------------------------------------------------
# Type flags and attributes
# ---------------------------------------------------------------------------

def test_type_flags():
    """is_AlgebraicZeroTensor, is_commutative, is_zero flags."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.is_AlgebraicZeroTensor is True
    assert Z.is_commutative is True
    assert Z.is_zero is True


def test_op_priority():
    """_op_priority is set higher than Expr default."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z._op_priority == 11


def test_free_symbols():
    """Zero tensor has no free symbols."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z.free_symbols == set()


# ---------------------------------------------------------------------------
# Hashing and equality
# ---------------------------------------------------------------------------

def test_hashable_content():
    """_hashable_content returns the expected tuple."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z._hashable_content() == (((3, 4), (4, 5)),)


def test_equality_same_shape():
    """Two zero tensors with the same shape are equal."""
    Z1 = AlgebraicZeroTensor((3, 4))
    Z2 = AlgebraicZeroTensor((3, 4))
    assert Z1 == Z2


def test_equality_different_shape():
    """Two zero tensors with different shapes are not equal."""
    Z1 = AlgebraicZeroTensor((3, 4))
    Z2 = AlgebraicZeroTensor((2, 3))
    assert Z1 != Z2


def test_equality_multi_factor():
    """Multi-factor zero tensors compare correctly."""
    Z1 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    Z2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    Z3 = AlgebraicZeroTensor(((3, 4), (5, 6)))
    assert Z1 == Z2
    assert Z1 != Z3


def test_hash_consistency():
    """Hash is consistent with equality."""
    Z1 = AlgebraicZeroTensor((3, 4))
    Z2 = AlgebraicZeroTensor((3, 4))
    assert hash(Z1) == hash(Z2)
    assert Z1 == Z2


def test_set_membership():
    """Zero tensors work correctly in sets."""
    Z1 = AlgebraicZeroTensor((3, 4))
    Z2 = AlgebraicZeroTensor((3, 4))
    Z3 = AlgebraicZeroTensor((2, 3))
    s = {Z1, Z2, Z3}
    assert len(s) == 2


# ---------------------------------------------------------------------------
# Pickle
# ---------------------------------------------------------------------------

def test_pickle():
    """Zero tensor can be pickled and unpickled."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    data = pickle.dumps(Z)
    Z2 = pickle.loads(data)
    assert Z2 == Z
    assert Z2.shape == Z.shape


def test_getnewargs():
    """__getnewargs__ returns the shape for pickle reconstruction."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    args = Z.__getnewargs__()
    assert args == (((3, 4), (4, 5)),)


# ---------------------------------------------------------------------------
# commutativity_pattern
# ---------------------------------------------------------------------------

def test_commutativity_pattern_single():
    """Single-factor zero tensor has all-1s pattern."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.commutativity_pattern == (1,)


def test_commutativity_pattern_multi():
    """Multi-factor zero tensor has all-1s pattern."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z.commutativity_pattern == (1, 1)


def test_commutativity_pattern_three():
    """Three-factor zero tensor has all-1s pattern."""
    Z = AlgebraicZeroTensor(((2, 3), (3, 4), (4, 5)))
    assert Z.commutativity_pattern == (1, 1, 1)


# ---------------------------------------------------------------------------
# Negation
# ---------------------------------------------------------------------------

def test_negation():
    """Negation of zero is itself."""
    Z = AlgebraicZeroTensor((3, 4))
    assert -Z is Z


def test_negation_multi_factor():
    """Negation of multi-factor zero is itself."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert -Z is Z


# ---------------------------------------------------------------------------
# Addition (additive identity)
# ---------------------------------------------------------------------------

def test_add_pure_tensor():
    """Zero + PureTensor returns the PureTensor."""
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = Z + T
    assert result == T


def test_radd_pure_tensor():
    """PureTensor + Zero returns the PureTensor."""
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = T + Z
    assert result == T


def test_add_zero_with_zero():
    """Zero + Zero returns the other zero."""
    Z1 = AlgebraicZeroTensor((3, 4))
    Z2 = AlgebraicZeroTensor((3, 4))
    result = Z1 + Z2
    assert result == Z2


def test_add_bare_matrix():
    """Zero + bare matrix returns the bare matrix."""
    Z = AlgebraicZeroTensor(((3, 4),))
    result = Z + A
    assert result == A


def test_radd_bare_matrix():
    """Bare matrix + Zero goes through MatrixSymbol.__add__ which raises TypeError
    because AlgebraicZeroTensor is not a MatrixExpr. The AlgebraicZeroTensor.__radd__
    is not called because MatrixSymbol has its own __add__ that doesn't return NotImplemented."""
    Z = AlgebraicZeroTensor(((3, 4),))
    raises(TypeError, lambda: A + Z)


# ---------------------------------------------------------------------------
# Subtraction
# ---------------------------------------------------------------------------

def test_sub_pure_tensor():
    """Zero - PureTensor returns negated PureTensor (zero dropped)."""
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = Z - T
    assert isinstance(result, AlgebraicPureTensor)
    assert result == -T


def test_rsub_pure_tensor():
    """PureTensor - Zero returns PureTensor (zero dropped)."""
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = T - Z
    assert isinstance(result, AlgebraicPureTensor)
    assert result == T


def test_sub_zero_with_zero():
    """Zero - Zero returns AlgebraicZeroTensor (AlgebraicTensor simplifies)."""
    Z1 = AlgebraicZeroTensor((3, 4))
    Z2 = AlgebraicZeroTensor((3, 4))
    result = Z1 - Z2
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4),)


def test_rsub_zero_with_zero():
    """Zero - Zero (reversed) returns AlgebraicZeroTensor."""
    Z1 = AlgebraicZeroTensor((3, 4))
    Z2 = AlgebraicZeroTensor((3, 4))
    result = Z2 - Z1
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4),)


def test_rsub_zero_with_zero_dupe():
    """Zero - Zero (reversed) returns AlgebraicZeroTensor."""
    Z1 = AlgebraicZeroTensor((3, 4))
    Z2 = AlgebraicZeroTensor((3, 4))
    result = Z2 - Z1
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4),)


# ---------------------------------------------------------------------------
# Multiplication (scalar scaling)
# ---------------------------------------------------------------------------

def test_mul_number():
    """Zero * number returns self."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z * 5 is Z


def test_mul_zero():
    """Zero * S.Zero returns self."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z * S.Zero is Z


def test_mul_symbol():
    """Zero * commutative symbol returns self."""
    x = Symbol("x")
    Z = AlgebraicZeroTensor((3, 4))
    assert Z * x is Z


def test_rmul_number():
    """Number * Zero returns self."""
    Z = AlgebraicZeroTensor((3, 4))
    assert 5 * Z is Z


def test_rmul_zero():
    """S.Zero * Zero returns self."""
    Z = AlgebraicZeroTensor((3, 4))
    assert S.Zero * Z is Z


def test_rmul_symbol():
    """Commutative symbol * Zero returns self."""
    x = Symbol("x")
    Z = AlgebraicZeroTensor((3, 4))
    assert x * Z is Z


def test_mul_pure_tensor():
    """Zero * PureTensor delegates to composition."""
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = Z * T
    # Composition with zero returns a zero tensor
    assert isinstance(result, AlgebraicZeroTensor)


def test_rmul_pure_tensor():
    """PureTensor * Zero delegates to composition."""
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = T * Z
    # Composition with zero returns a zero tensor
    assert isinstance(result, AlgebraicZeroTensor)


# ---------------------------------------------------------------------------
# Boolean
# ---------------------------------------------------------------------------

def test_bool():
    """Zero tensor evaluates to False in boolean context."""
    Z = AlgebraicZeroTensor((3, 4))
    assert not Z
    assert bool(Z) is False


def test_bool_in_condition():
    """Zero tensor works correctly in if statement."""
    Z = AlgebraicZeroTensor((3, 4))
    if Z:
        assert False, "Zero tensor should be falsy"


# ---------------------------------------------------------------------------
# Transpose
# ---------------------------------------------------------------------------

def test_transpose_single():
    """Transpose of single-factor zero tensor reverses shape."""
    Z = AlgebraicZeroTensor((3, 4))
    ZT = Z.T
    assert ZT.shape == ((4, 3),)
    assert ZT.commutativity_pattern == (1,)
    assert isinstance(ZT, AlgebraicZeroTensor)


def test_transpose_multi():
    """Transpose of multi-factor zero tensor reverses each factor."""
    Z = AlgebraicZeroTensor(((1, 2), (3, 4)))
    ZT = Z.T
    assert ZT.shape == ((2, 1), (4, 3))
    assert ZT.commutativity_pattern == (1, 1)


def test_transpose_double():
    """Double transpose returns original shape."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z.T.T.shape == Z.shape


def test_transpose_square():
    """Transpose of square factors keeps same shape."""
    Z = AlgebraicZeroTensor(((3, 3), (4, 4)))
    ZT = Z.T
    assert ZT.shape == ((3, 3), (4, 4))
    assert ZT == Z


# ---------------------------------------------------------------------------
# Copy
# ---------------------------------------------------------------------------

def test_copy():
    """copy() returns self (immutable atom)."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.copy() is Z


# ---------------------------------------------------------------------------
# Expand
# ---------------------------------------------------------------------------

def test_expand():
    """expand() returns self (already expanded)."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z.expand() is Z


def test_expand_with_kwargs():
    """expand() with kwargs returns self."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.expand(deep=True) is Z


# ---------------------------------------------------------------------------
# SymPy integration
# ---------------------------------------------------------------------------

def test_doit():
    """doit() returns self (zero tensor has nothing to evaluate)."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.doit() is Z


def test_doit_multi_factor():
    """doit() on multi-factor zero tensor returns self."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z.doit() is Z


def test_doit_with_hints():
    """doit() with hints returns self."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.doit(deep=True) is Z
    assert Z.doit(deep=False) is Z


def test_subs():
    """subs() returns self (no symbols to substitute)."""
    x = Symbol("x")
    y = Symbol("y")
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.subs(x, y) is Z


def test_xreplace():
    """xreplace() with no matching rule returns self."""
    Z = AlgebraicZeroTensor((3, 4))
    result = Z.xreplace({})
    assert result is Z


def test_atoms():
    """atoms() returns self (it is an atom)."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.atoms() == {Z}


def test_as_coeff_Add():
    """as_coeff_Add() works (inherited from Expr)."""
    Z = AlgebraicZeroTensor((3, 4))
    coeff, addends = Z.as_coeff_Add()
    # Zero-like: coefficient is 0, remaining is the zero tensor itself
    assert coeff == 0
    assert addends == Z


def test_as_coeff_Mul():
    """as_coeff_Mul() works (inherited from Expr)."""
    Z = AlgebraicZeroTensor((3, 4))
    coeff, mul_args = Z.as_coeff_Mul()
    assert coeff == S.One


def test_as_base_exp():
    """as_base_exp() works (inherited from Expr)."""
    Z = AlgebraicZeroTensor((3, 4))
    base, exp = Z.as_base_exp()
    assert base == Z
    assert exp == S.One


# ---------------------------------------------------------------------------
# Identity preservation across operations
# ---------------------------------------------------------------------------

def test_identity_after_negation():
    """Type is preserved after negation."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = -Z
    assert isinstance(result, AlgebraicZeroTensor)


def test_identity_after_transpose():
    """Type is preserved after transpose."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = Z.T
    assert isinstance(result, AlgebraicZeroTensor)


def test_identity_after_expand():
    """Type is preserved after expand."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = Z.expand()
    assert isinstance(result, AlgebraicZeroTensor)


def test_identity_after_copy():
    """Type is preserved after copy."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = Z.copy()
    assert isinstance(result, AlgebraicZeroTensor)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_shape_with_symbolic_dimensions():
    """Zero tensor can carry symbolic dimensions."""
    n = Symbol("n", integer=True)
    m = Symbol("m", integer=True)
    Z = AlgebraicZeroTensor((n, m))
    assert Z.shape == ((n, m),)


def test_different_shapes_not_equal():
    """Different shapes produce distinct objects."""
    shapes = [
        ((2, 3),),
        ((2, 3), (3, 4)),
        ((2, 3), (3, 4), (4, 5)),
        ((5, 6),),
    ]
    zeros = [AlgebraicZeroTensor(s) for s in shapes]
    for i in range(len(zeros)):
        for j in range(i + 1, len(zeros)):
            assert zeros[i] != zeros[j], f"Shapes {shapes[i]} and {shapes[j]} should differ"


def test_hash_different_shapes():
    """Different shapes produce different hashes."""
    Z1 = AlgebraicZeroTensor((3, 4))
    Z2 = AlgebraicZeroTensor((2, 3))
    assert hash(Z1) != hash(Z2)


# ---------------------------------------------------------------------------
# Diff
# ---------------------------------------------------------------------------

def test_diff_returns_self():
    """diff() returns self for any symbol."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.diff(x) is Z


def test_diff_multi_symbol():
    """diff() with multiple symbols returns self."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z.diff(x, y) is Z


def test_diff_second_order():
    """Second-order diff() returns self."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.diff(x, x) is Z


def test_diff_with_assumptions():
    """diff() with assumptions kwargs returns self."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.diff(x, evaluate=True) is Z


# ---------------------------------------------------------------------------
# Conjugate
# ---------------------------------------------------------------------------

def test_conjugate_returns_self():
    """conjugate() returns self for a single-factor zero tensor."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.conjugate() is Z


def test_conjugate_multi_factor():
    """conjugate() returns self for a multi-factor zero tensor."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z.conjugate() is Z


def test_conjugate_preserves_shape():
    """conjugate() preserves the zero tensor shape."""
    Z = AlgebraicZeroTensor(((2, 3), (3, 4), (4, 5)))
    assert Z.conjugate().shape == Z.shape


def test_conjugate_preserves_type():
    """conjugate() preserves the AlgebraicZeroTensor type."""
    Z = AlgebraicZeroTensor((3, 4))
    assert isinstance(Z.conjugate(), AlgebraicZeroTensor)
