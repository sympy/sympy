"""Playground for testing the algebraic tensor module.

Run from the sympy root directory:
    PYTHONPATH=. python sympy/tensor/algebraic/tests/playground_functions.py
"""

import os, sys
# Allow running from any directory by adding the sympy package root to sys.path.
_HERE = os.path.dirname(os.path.abspath(__file__))
_SYMPY_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(_HERE))))
if _SYMPY_ROOT not in sys.path:
    sys.path.insert(0, _SYMPY_ROOT)

from sympy import Symbol, Rational, sqrt, pi, I
from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol, MatAdd, MatMul

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    algebraic_tensor_product,
    algebraic_zero_tensor,
    compose_algebraic_pure_tensors,
    compose_algebraic_tensors,
    tensorsimplify,
    proportionality_factoring,
    ShapeMismatchError,
)

# ---------------------------------------------------------------------------
# Commutative symbols
# ---------------------------------------------------------------------------
a, b, c, x, y, z = Symbol("a"), Symbol("b"), Symbol("c"), Symbol("x"), Symbol("y"), Symbol("z")

# ---------------------------------------------------------------------------
# Matrix symbols – various shapes for testing
# ---------------------------------------------------------------------------
# 2x3 matrices
A = MatrixSymbol("A", 2, 3)
A2 = MatrixSymbol("A2", 2, 3)

# 3x4 matrices
B = MatrixSymbol("B", 3, 4)
B2 = MatrixSymbol("B2", 3, 4)

# 3x3 matrices (square, for composition with self-shape)
M = MatrixSymbol("M", 3, 3)
N = MatrixSymbol("N", 3, 3)

# 4x5 matrices
C = MatrixSymbol("C", 4, 5)
C2 = MatrixSymbol("C2", 4, 5)

# 2x2 identity-like square
P = MatrixSymbol("P", 2, 2)
Q = MatrixSymbol("Q", 2, 2)

# 4x2 matrices (for composition chain: 2x3 * 3x4 * 4x2)
R = MatrixSymbol("R", 4, 2)
R2 = MatrixSymbol("R2", 4, 2)

# 2x2 for chaining after R
S_mat = MatrixSymbol("S", 2, 2)
S2 = MatrixSymbol("S2", 2, 2)

# ---------------------------------------------------------------------------
# 1. LINEAR COMBINATIONS
# ---------------------------------------------------------------------------
print("=" * 60)
print("1. LINEAR COMBINATIONS")
print("=" * 60)

# --- Pure tensors ---
print("\n--- Pure Tensors ---")

# Single tensor product
pt1 = AlgebraicPureTensor(A, B)
print(f"pt1 = A ⊗ B              : {pt1}")
print(f"  type   : {type(pt1).__name__}")
print(f"  shape  : {pt1.tensor_shape}")

# With integer coefficient
pt2 = AlgebraicPureTensor(2, A, B)
print(f"\npt2 = 2*(A ⊗ B)         : {pt2}")
print(f"  coeff  : {pt2._get_coeff()}")
print(f"  factors: {[str(f) for f in pt2.factors]}")

# With symbolic coefficient
pt3 = AlgebraicPureTensor(x, A, B)
print(f"\npt3 = x*(A ⊗ B)         : {pt3}")
print(f"  coeff  : {pt3._get_coeff()}")

# Three-factor tensor product
pt4 = AlgebraicPureTensor(A, B, C)
print(f"\npt4 = A ⊗ B ⊗ C         : {pt4}")
print(f"  shape  : {pt4.tensor_shape}")

# Scalar multiplication from left
pt5 = 3 * AlgebraicPureTensor(A, B)
print(f"\npt5 = 3 * (A ⊗ B)       : {pt5}")
print(f"  type   : {type(pt5).__name__}")
print(f"  coeff  : {pt5._get_coeff()}")

# Scalar multiplication from right
pt6 = AlgebraicPureTensor(A, B) * 3
print(f"\npt6 = (A ⊗ B) * 3       : {pt6}")
print(f"  coeff  : {pt6._get_coeff()}")

# Symbolic scalar from left
pt7 = x * AlgebraicPureTensor(A, B)
print(f"\npt7 = x * (A ⊗ B)       : {pt7}")
print(f"  coeff  : {pt7._get_coeff()}")

# Negation
pt8 = -AlgebraicPureTensor(A, B)
print(f"\npt8 = -(A ⊗ B)          : {pt8}")
print(f"  coeff  : {pt8._get_coeff()}")

# --- Addition (creates AlgebraicTensor) ---
print("\n--- Addition ---")

at1 = AlgebraicPureTensor(A, B) + AlgebraicPureTensor(A2, B2)
print(f"at1 = A⊗B + A2⊗B2        : {at1}")
print(f"  type   : {type(at1).__name__}")
print(f"  shape  : {at1.tensor_shape}")

at2 = AlgebraicPureTensor(2, A, B) + AlgebraicPureTensor(3, A2, B2)
print(f"\nat2 = 2*A⊗B + 3*A2⊗B2  : {at2}")
print(f"  type   : {type(at2).__name__}")

at3 = AlgebraicPureTensor(A, B) + AlgebraicPureTensor(A2, B2) + AlgebraicPureTensor(A, B2)
print(f"\nat3 = A⊗B + A2⊗B2 + A⊗B2 : {at3}")
print(f"  type   : {type(at3).__name__}")
print(f"  terms  : {len(at3.terms)}")

# --- Subtraction ---
print("\n--- Subtraction ---")

at4 = AlgebraicPureTensor(A, B) - AlgebraicPureTensor(A2, B2)
print(f"at4 = A⊗B - A2⊗B2        : {at4}")

at5 = AlgebraicPureTensor(2, A, B) - AlgebraicPureTensor(3, A, B)
print(f"\nat5 = 2*A⊗B - 3*A⊗B     : {at5}")

# --- Scalar times AlgebraicTensor ---
print("\n--- Scalar * AlgebraicTensor ---")

at6 = x * AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(A2, B2))
print(f"at6 = x*(A⊗B + A2⊗B2)    : {at6}")
print(f"  type   : {type(at6).__name__}")

at7 = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(A2, B2)) * y
print(f"\nat7 = (A⊗B + A2⊗B2)*y   : {at7}")

# --- Zero tensor ---
print("\n--- Zero Tensor ---")

z1 = AlgebraicZeroTensor(((2, 3), (3, 4)))
print(f"z1 = zero tensor          : {z1}")
print(f"  shape  : {z1.shape}")

at8 = AlgebraicPureTensor(A, B) + z1
print(f"\nat8 = A⊗B + zero         : {at8}")
print(f"  type   : {type(at8).__name__}")
print(f"  has_zero: {at8.has_zero_term()}")

# --- Unwrapping: single-term AlgebraicTensor ---
print("\n--- Unwrapping ---")

at9 = AlgebraicTensor(AlgebraicPureTensor(A, B))
print(f"at9 = AlgebraicTensor(A⊗B)  : {at9}")
print(f"  type   : {type(at9).__name__}  (should unwrap to AlgebraicPureTensor)")

# --- Complex linear combination ---
print("\n--- Complex Linear Combination ---")

at10 = (
    AlgebraicPureTensor(2, A, B)
    + AlgebraicPureTensor(x, A2, B2)
    - AlgebraicPureTensor(3, A, B2)
)
print(f"at10 = 2*A⊗B + x*A2⊗B2 - 3*A⊗B2")
print(f"  result: {at10}")
print(f"  type   : {type(at10).__name__}")
print(f"  terms  : {len(at10.terms)}")

# ---------------------------------------------------------------------------
# 2. COMPOSITION (* operation)
# ---------------------------------------------------------------------------
print("=" * 60)
print("2. COMPOSITION")
print("=" * 60)

# --- PureTensor * PureTensor ---
print("\n--- PureTensor * PureTensor ---")

# A(2,3) * A^T-like: need compatible shapes
# M(3,3) * N(3,3) -> MatMul(M,N)(3,3)
comp1 = AlgebraicPureTensor(M) * AlgebraicPureTensor(N)
print(f"comp1 = M * N              : {comp1}")
print(f"  type   : {type(comp1).__name__}")

# Two-factor composition
comp2 = AlgebraicPureTensor(M, N) * AlgebraicPureTensor(N, M)
print(f"\ncomp2 = (M⊗N) * (N⊗M)    : {comp2}")
print(f"  type   : {type(comp2).__name__}")

# --- PureTensor with coefficient * PureTensor ---
print("\n--- Coefficient handling in composition ---")

comp3 = AlgebraicPureTensor(2, M, N) * AlgebraicPureTensor(3, N, M)
print(f"comp3 = (2*M⊗N) * (3*N⊗M) : {comp3}")
print(f"  coeff  : {comp3._get_coeff()}")

comp4 = AlgebraicPureTensor(x, M) * AlgebraicPureTensor(y, N)
print(f"\ncomp4 = (x*M) * (y*N)     : {comp4}")
print(f"  coeff  : {comp4._get_coeff()}")

# --- AlgebraicTensor * AlgebraicTensor ---
print("\n--- AlgebraicTensor * AlgebraicTensor ---")

at_comp1 = AlgebraicTensor(AlgebraicPureTensor(M), AlgebraicPureTensor(N))
at_comp2 = AlgebraicTensor(AlgebraicPureTensor(N), AlgebraicPureTensor(M))
comp5 = at_comp1 * at_comp2
print(f"comp5 = (M+N)*(N+M)        : {comp5}")
print(f"  type   : {type(comp5).__name__}")

# --- AlgebraicTensor * PureTensor ---
print("\n--- AlgebraicTensor * PureTensor ---")

comp6 = AlgebraicTensor(AlgebraicPureTensor(M), AlgebraicPureTensor(N)) * AlgebraicPureTensor(M)
print(f"comp6 = (M+N)*M            : {comp6}")
print(f"  type   : {type(comp6).__name__}")

# --- PureTensor * AlgebraicTensor ---
print("\n--- PureTensor * AlgebraicTensor ---")

# Use 2-factor PureTensor so it doesn't unwrap to bare matrix
comp7 = AlgebraicPureTensor(M, N) * AlgebraicTensor(AlgebraicPureTensor(N, M), AlgebraicPureTensor(M, N))
print(f"comp7 = (M⊗N)*(N⊗M + M⊗N)  : {comp7}")
print(f"  type   : {type(comp7).__name__}")

# Also test with compose_algebraic_tensors directly for bare matrix * AlgebraicTensor
comp7b = compose_algebraic_tensors(M, AlgebraicTensor(AlgebraicPureTensor(M), AlgebraicPureTensor(N)))
print(f"\ncomp7b = compose(M, M+N)   : {comp7b}")
print(f"  type   : {type(comp7b).__name__}")

# --- Composition with zero tensor ---
print("\n--- Composition with Zero Tensor ---")

z2 = AlgebraicZeroTensor(((3, 3),))
comp8 = AlgebraicPureTensor(M) * z2
print(f"comp8 = M * zero           : {comp8}")
print(f"  type   : {type(comp8).__name__}")

comp9 = z2 * AlgebraicPureTensor(N)
print(f"\ncomp9 = zero * N          : {comp9}")
print(f"  type   : {type(comp9).__name__}")

# --- Multi-factor composition chain ---
print("\n--- Multi-factor composition ---")

# Need compatible shapes: M(3,3)*N(3,3), N(3,3)*M(3,3), M(3,3)*N(3,3)
comp10 = AlgebraicPureTensor(M, N, M) * AlgebraicPureTensor(N, M, N)
print(f"comp10 = (M⊗N⊗M)*(N⊗M⊗N)    : {comp10}")
print(f"  type   : {type(comp10).__name__}")

# Composition with coefficient
comp11 = AlgebraicPureTensor(2, M, N) * AlgebraicPureTensor(3, N, M)
print(f"\ncomp11 = (2*M⊗N)*(3*N⊗M)     : {comp11}")
print(f"  coeff  : {comp11._get_coeff()}")

# ---------------------------------------------------------------------------
# 3. SIMPLIFICATION AND EXPANSION
# ---------------------------------------------------------------------------
print("=" * 60)
print("3. SIMPLIFICATION AND EXPANSION")
print("=" * 60)

# --- PureTensor .simplify() ---
print("\n--- PureTensor .simplify() ---")

pt_simp1 = AlgebraicPureTensor(2 + 2, A, B)
print(f"pt_simp1 = (2+2)*(A⊗B)     : {pt_simp1}")
simp1 = pt_simp1.simplify()
print(f"  simplified              : {simp1}")
print(f"  coeff  : {simp1._get_coeff()}")

pt_simp2 = AlgebraicPureTensor(x + x, M, N)
print(f"\npt_simp2 = (x+x)*(M⊗N)    : {pt_simp2}")
simp2 = pt_simp2.simplify()
print(f"  simplified              : {simp2}")
print(f"  coeff  : {simp2._get_coeff()}")

# --- PureTensor .expand() ---
print("\n--- PureTensor .expand() ---")

pt_exp1 = AlgebraicPureTensor(A, B + B2)
print(f"pt_exp1 = A ⊗ (B+B2)       : {pt_exp1}")
exp1 = pt_exp1.expand()
print(f"  expanded                : {exp1}")
print(f"  type   : {type(exp1).__name__}")
print(f"  terms  : {len(exp1.terms) if hasattr(exp1, 'terms') else 'N/A'}")

pt_exp2 = AlgebraicPureTensor(A + A2, B + B2)
print(f"\npt_exp2 = (A+A2) ⊗ (B+B2) : {pt_exp2}")
exp2 = pt_exp2.expand()
print(f"  expanded                : {exp2}")
print(f"  type   : {type(exp2).__name__}")
print(f"  terms  : {len(exp2.terms) if hasattr(exp2, 'terms') else 'N/A'}")

pt_exp3 = AlgebraicPureTensor(3, A, B + B2)
print(f"\npt_exp3 = 3*(A ⊗ (B+B2))  : {pt_exp3}")
exp3 = pt_exp3.expand()
print(f"  expanded                : {exp3}")
print(f"  type   : {type(exp3).__name__}")

# --- AlgebraicTensor .simplify() ---
print("\n--- AlgebraicTensor .simplify() ---")

at_simp1 = AlgebraicTensor(
    AlgebraicPureTensor(2, A, B),
    AlgebraicPureTensor(3, A, B),
)
print(f"at_simp1 = 2*A⊗B + 3*A⊗B   : {at_simp1}")
simp_at1 = at_simp1.simplify()
print(f"  simplified              : {simp_at1}")
print(f"  type   : {type(simp_at1).__name__}")

at_simp2 = AlgebraicPureTensor(2, A, B) -2*AlgebraicPureTensor(A, B)

print(f"\nat_simp2 = 2*A⊗B - 2*A⊗B  : {at_simp2}")
simp_at2 = at_simp2.simplify()
print(f"  simplified              : {simp_at2}")
print(f"  type   : {type(simp_at2).__name__}")

at_simp3 = AlgebraicTensor(
    AlgebraicPureTensor(A, B),
    AlgebraicPureTensor(A2, B),
)
print(f"\nat_simp3 = A⊗B + A2⊗B     : {at_simp3}")
simp_at3 = at_simp3.simplify()
print(f"  simplified              : {simp_at3}")

# --- AlgebraicTensor .expand() ---
print("\n--- AlgebraicTensor .expand() ---")

at_exp1 = AlgebraicTensor(
    AlgebraicPureTensor(A, B + B2),
    AlgebraicPureTensor(A2, B),
)
print(f"at_exp1 = A⊗(B+B2) + A2⊗B  : {at_exp1}")
exp_at1 = at_exp1.expand()
print(f"  expanded                : {exp_at1}")
print(f"  type   : {type(exp_at1).__name__}")
print(f"  terms  : {len(exp_at1.terms) if hasattr(exp_at1, 'terms') else 'N/A'}")

# --- tensorsimplify() public API ---
print("\n--- tensorsimplify() public API ---")

ts1 = tensorsimplify(AlgebraicPureTensor(2 + Rational(1, 2) - Rational(1, 2), A, B))
print(f"ts1 = simplify((2+1/2-1/2)*A⊗B) : {ts1}")

ts2 = tensorsimplify(AlgebraicTensor(
    AlgebraicPureTensor(2, A, B),
    AlgebraicPureTensor(3, A, B),
))
print(f"\nts2 = simplify(2*A⊗B + 3*A⊗B) : {ts2}")

# --- proportionality_factoring ---
print("\n--- proportionality_factoring() ---")

pf1 = proportionality_factoring(AlgebraicTensor(
    AlgebraicPureTensor(2, A, B),
    AlgebraicPureTensor(3, A, B),
))
print(f"pf1 = factor(2*A⊗B + 3*A⊗B)  : {pf1}")
print(f"  type   : {type(pf1).__name__}")

pf2 = proportionality_factoring(AlgebraicTensor(
    AlgebraicPureTensor(A, B),
    AlgebraicPureTensor(A2, B),
))
print(f"\npf2 = factor(A⊗B + A2⊗B)    : {pf2}")
print(f"  type   : {type(pf2).__name__}")

# ---------------------------------------------------------------------------
# 4. EDGE CASES
# ---------------------------------------------------------------------------
print("=" * 60)
print("4. EDGE CASES")
print("=" * 60)

# Unwrap: single factor with coeff=1
print("\n--- Unwrapping ---")
uw1 = AlgebraicPureTensor(M)
print(f"uw1 = AlgebraicPureTensor(M) : {uw1}")
print(f"  type   : {type(uw1).__name__}  (should unwrap to MatrixSymbol)")

# Zero coefficient
uw2 = AlgebraicPureTensor(S.Zero, A, B)
print(f"\nuw2 = AlgebraicPureTensor(0, A, B) : {uw2}")
print(f"  type   : {type(uw2).__name__}")

# Zero tensor expand
z_exp = AlgebraicZeroTensor(((2, 3), (3, 4))).expand()
print(f"\nz_exp = zero.expand()        : {z_exp}")
print(f"  type   : {type(z_exp).__name__}")

# Identity expand
id_exp = AlgebraicPureTensor(A, B).expand()
print(f"\nid_exp = (A⊗B).expand()      : {id_exp}")
print(f"  is_self: {id_exp is AlgebraicPureTensor(A, B)}")

print("\n" + "=" * 60)
print("READY — modify this file to add your own tests.")
print("=" * 60)
