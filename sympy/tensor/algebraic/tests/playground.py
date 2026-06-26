"""Interactive playground for algebraic tensors.

Run with:
    python3 tensor/algebraic/tests/playground.py
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..")))

from sympy.matrices.expressions import MatrixSymbol, Transpose, MatAdd

from sympy.tensor.algebraic.zero_tensor import ZeroTensor, zero_tensor
from sympy.tensor.algebraic.pure_tensor import PureTensor, tensor_product
from sympy.tensor.algebraic.algebraic_tensor import (
    AlgebraicTensor,
    ShapeMismatchError,
)

# ---------------------------------------------------------------------------
# Matrix symbols of various shapes
# ---------------------------------------------------------------------------

# Regular matrices
A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)

# Vectors
x = MatrixSymbol("x", 5, 1)       # column vector
y = MatrixSymbol("y", 1, 5)       # row vector
v = MatrixSymbol("v", 4, 1)       # column vector (compatible with A's columns)

# Transposed vectors
xT = Transpose(x)                 # 1x5 row
vT = Transpose(v)                 # 1x4 row

# Identity matrices
I3 = MatrixSymbol("I3", 3, 3)
I4 = MatrixSymbol("I4", 4, 4)
I5 = MatrixSymbol("I5", 5, 5)

# 1x1 noncommutative matrix
u = MatrixSymbol("u", 1, 1)

# ---------------------------------------------------------------------------
# Nontrivial PureTensors
# ---------------------------------------------------------------------------

# --- Basic two-factor product ---
t1 = tensor_product(A, C)
print(f"t1 = tensor_product(A, C)")
print(f"  type   : {type(t1).__name__}")
print(f"  shape  : {t1.tensor_shape}")
print(f"  str    : {t1}")
print()

# --- Four-factor chain ---
t2 = PureTensor(A, I4, C, x)
print(f"t2 = PureTensor(A, I4, C, x)")
print(f"  type   : {type(t2).__name__}")
print(f"  shape  : {t2.tensor_shape}")
print(f"  str    : {t2}")
print()

# --- Vector ⊗ matrix ⊗ row-vector ---
t3 = PureTensor(x, C, vT)
print(f"t3 = PureTensor(x, C, vT)")
print(f"  type   : {type(t3).__name__}")
print(f"  shape  : {t3.tensor_shape}")
print(f"  str    : {t3}")
print()

# --- Row ⊗ column ⊗ matrix (mixed spaces) ---
t4 = PureTensor(y, x, A)
print(f"t4 = PureTensor(y, x, A)")
print(f"  type   : {type(t4).__name__}")
print(f"  shape  : {t4.tensor_shape}")
print(f"  str    : {t4}")
print()

# --- 1x1 noncommutative matrix in the middle ---
t5 = PureTensor(A, u, C)
print(f"t5 = PureTensor(A, u, C)")
print(f"  type   : {type(t5).__name__}")
print(f"  shape  : {t5.tensor_shape}")
print(f"  str    : {t5}")
print()

# --- Scalar-scaled PureTensor ---
t6 = 3 * t1
print(f"t6 = 3 * t1")
print(f"  type   : {type(t6).__name__}")
print(f"  str    : {t6}")
print()

# ---------------------------------------------------------------------------
# AlgebraicTensor: sums of same-shape PureTensors
# ---------------------------------------------------------------------------

s1 = t1 + PureTensor(B, C)
print(f"s1 = tensor_product(A, C) + tensor_product(B, C)")
print(f"  type   : {type(s1).__name__}")
print(f"  shape  : {s1.tensor_shape}")
print(f"  str    : {s1}")
print()

s2 = PureTensor(A, I3) + PureTensor(B, I3) + PureTensor(A, I3)
print(f"s2 = PureTensor(A, I3) + PureTensor(B, I3) + PureTensor(A, I3)")
print(f"  type   : {type(s2).__name__}")
print(f"  shape  : {s2.tensor_shape}")
print(f"  str    : {s2}")
print()

# --- With coefficients ---
s3 = AlgebraicTensor(2 * t1, -3 * PureTensor(B, C))
print(f"s3 = AlgebraicTensor(2*t1, -3*PureTensor(B, C))")
print(f"  type   : {type(s3).__name__}")
print(f"  shape  : {s3.tensor_shape}")
print(f"  str    : {s3}")
print()

# ---------------------------------------------------------------------------
# ZeroTensor
# ---------------------------------------------------------------------------

z1 = zero_tensor(((3, 4), (4, 5)))
print(f"z1 = zero_tensor(((3,4), (4,5)))")
print(f"  type   : {type(z1).__name__}")
print(f"  shape  : {z1.shape}")
print(f"  str    : {z1}")
print(f"  bool   : {bool(z1)}")
print()

z2 = z1 + t1
print(f"z2 = z1 + t1  (ZeroTensor + PureTensor)")
print(f"  type   : {type(z2).__name__}")
print(f"  str    : {z2}")
print()

# ---------------------------------------------------------------------------
# Shape mismatch (commented out — raises ShapeMismatchError)
# ---------------------------------------------------------------------------

# PureTensor(A, C) + PureTensor(C, D)   # ((3,4),(4,5)) vs ((4,5),(3,5))
# PureTensor(A, I4, C) + PureTensor(A, C)  # 3 factors vs 2 factors

# ---------------------------------------------------------------------------
# Factorization helpers
# ---------------------------------------------------------------------------

s4 = PureTensor(A, I4, C) + PureTensor(A, I4, C)
print(f"s4 = PureTensor(A, I4, C) + PureTensor(A, I4, C)")
left, mid, right = s4.as_common_factors()
print(f"  as_common_factors: left={left}, mid={mid}, right={right}")
print()

# --- Factorization with differing middle factor ---
J4 = MatrixSymbol("J4", 4, 4)
s5 = PureTensor(A, I4, C) + PureTensor(A, J4, C)
print(f"s5 = PureTensor(A, I4, C) + PureTensor(A, J4, C)")
left, mid, right = s5.as_common_factors()
print(f"  as_common_factors: left={left}, mid={mid}, right={right}")
print()


t8 = t6 + 2*t6
print(f"  type   : {type(t8).__name__}")
print(f"  shape  : {t8.tensor_shape}")
print(f"  str    : {t8}")
print()
