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

from sympy import Symbol, Rational, sqrt, pi, I, symbols
from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol, MatAdd, MatMul
from sympy.matrices import Matrix, zeros, eye, diag

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    compose_algebraic_tensors,
    tensorsimplify,
    ShapeMismatchError,
    algebraic_tensor_product
)

from sympy.printing.latex import latex
from sympy.printing.repr import srepr

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
print(f"  shape  : {pt1.shape}")

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
print(f"  shape  : {pt4.shape}")

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
print(f"  shape  : {at1.shape}")

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




def direct_sum0(matrix1, matrix2):
    result = Matrix(zeros(sqrt(len(matrix1))+sqrt(len(matrix2)),sqrt(len(matrix1))+sqrt(len(matrix2))))
    for i in range (sqrt(len(matrix1))):
        for j in range (sqrt(len(matrix1))):
            result[i,j] = matrix1[i,j]

    for i in range (sqrt(len(matrix2))):
        for j in range (sqrt(len(matrix2))):
            result[i + sqrt(len(matrix1)),j + sqrt(len(matrix1))] = matrix2[i,j]

    return result


### Direct sum of arbitraty number of matrices

def direct_sum(matrices):
    result = matrices[0]
    for i in range (1,len(matrices)):
        result = direct_sum0(result, matrices[i])
    return result


upsilon_u, upsilon_d, upsilon_e, upsilon_nu, upsilon_R = symbols(\
        r"""\Upsilon_u, \Upsilon_d, \Upsilon_e, \Upsilon_\nu, \Upsilon_R""", commutative = False)
upsilont_u, upsilont_d, upsilont_e, upsilont_nu, upsilont_R = symbols(\
        r"""\Upsilon^t_u, \Upsilon^t_d, \Upsilon^t_e, \Upsilon^t_\nu, \Upsilon^t_R""", commutative = False)
upsilonc_u, upsilonc_d, upsilonc_e, upsilonc_nu, upsilonc_R = symbols(\
        r"""\Upsilon^*_u, \Upsilon^*_d, \Upsilon^*_e, \Upsilon^*_\nu, \Upsilon^*_R""", commutative = False)
cupsilon_u, cupsilon_d, cupsilon_e, cupsilon_nu, cupsilon_R = symbols(\
        r"""\overline{\Upsilon}_u, \overline{\Upsilon}_d, \overline{\Upsilon}_e,"""\
       +r"""\overline{\Upsilon}_\nu, \overline{\Upsilon}_R""", commutative = False)


D1 = AlgebraicPureTensor(Matrix([[1,0],[0,0]]), direct_sum([Matrix([1]), Matrix(zeros(3,3))]),Matrix(\
     [[0,0,upsilon_R, upsilonc_nu],[0,0,cupsilon_nu,0],[upsilonc_R,upsilont_nu,0,0],[upsilon_nu,0,0,0]]))

D2 = AlgebraicPureTensor(Matrix([[1,0],[0,0]]),direct_sum([Matrix([0]), eye(3)]), Matrix(\
     [[0,0,0,upsilonc_u],[0,0,cupsilon_u,0],[0,upsilont_u,0,0],[upsilon_u,0,0,0]]))

D3 = AlgebraicPureTensor(Matrix([[0,0],[0,1]]), direct_sum([Matrix([1]),Matrix(zeros(3,3))]),Matrix(\
     [[0,0,0,upsilonc_e],[0,0,cupsilon_e,0],[0,upsilont_e,0,0],[upsilon_e,0,0,0]]))

D4 = AlgebraicPureTensor(Matrix([[0,0],[0,1]]), direct_sum([Matrix([0]),eye(3)]), Matrix(\
     [[0,0,0,upsilonc_d],[0,0,cupsilon_d,0],[0,upsilont_d,0,0],[upsilon_d,0,0,0]]))

Dirac = D1 + D2 + D3 + D4

### Defining representation pi, algebra elements a0,...,a5 and the lists of symbols, list of d(symbol)s, list of -simbols (they are relevant)

def rep(z, w, alpha, beta, gamma, delta, m):
    pi1 = AlgebraicPureTensor(Matrix([[z,0],[0,w]]), eye(4), direct_sum([Matrix([1]), Matrix(zeros(3,3))]))
    pi2 = AlgebraicPureTensor(Matrix([[alpha, beta],[gamma,delta]]), eye(4),\
           direct_sum([Matrix(zeros(3,3)),Matrix([1])]))
    pi3 = AlgebraicPureTensor(eye(2), direct_sum([Matrix([z]), m]), Matrix(diag(0,1,1,0)))

    return pi1 + pi2 + pi3



### Defining two algebra elements, a1 and a2:
z0,z1,z2 = symbols(r"""z,z_1,z_2""", complex = True)
w0,w1,w2 = symbols(r"""w,w_1,w_2""", complex = True)
z3,z4,z5 = symbols(r"""z_3,z_4,z_5""", complex = True)
w3,w4,w5 = symbols(r"""w_3,w_4,w_5""", complex = True)
z6,w6 = symbols(r"""z_6, w_6""", complex = True)
z7,z8, w7, w8 = symbols(r"""z_7, z_8, w_7, w_8""", complex = True)

alpha0,alpha1, alpha2,beta0, beta1, beta2 = symbols(r"""\alpha,\alpha_1, \alpha_2,\beta, \beta_1, \beta_2""", complex = True)
gamma0,gamma1, gamma2,delta0, delta1, delta2 = symbols(r"""\gamma,\gamma_1, \gamma_2,\delta, \delta_1, \delta_2""", complex = True)
alpha3,alpha4, alpha5,beta3, beta4, beta5 = symbols(r"""\alpha_3,\alpha_4, \alpha_5,\beta_3, \beta_4, \beta_5""", complex = True)
gamma3,gamma4, gamma5,delta3, delta4, delta5 = symbols(r"""\gamma_3,\gamma_4, \gamma_5,\delta_3, \delta_4, \delta_5""", complex = True)
alpha5, alpha6, beta5, beta6 = symbols(r"""\alpha_5, \alpha_6, \beta_5, \beta_6""", complex = True)
gamma5, gamma6, delta5, delta6 = symbols(r"""\gamma_5, \gamma_6, \delta_5, \delta_6""", complex = True)
alpha7, alpha8, beta7, beta8 = symbols(r"""\alpha_7, \alpha_8, \beta_7, \beta_8""", complex = True)
gamma7, gamma8, delta7, delta8 = symbols(r"""\gamma_7, \gamma_8, \delta_7, \delta_8""", complex = True)


m111,m112,m113, m121,m122,m123, m131,m132,m133 = symbols(\
    r"""m^1_{11},m^1_{12},m^1_{13},m^1_{21},m^1_{22},m^1_{23},m^1_{31},m^1_{32},m^1_{33}""", complex=True)
m011,m012,m013, m021,m022,m023, m031,m032,m033 = symbols(\
    r"""m_{11},m_{12},m_{13},m_{21},m_{22},m_{23},m_{31},m_{32},m_{33}""", complex=True)
m211,m212,m213, m221,m222,m223, m231,m232,m233 = symbols(\
    r"""m^2_{11},m^2_{12},m^2_{13},m^2_{21},m^2_{22},m^2_{23},m^2_{31},m^2_{32},m^2_{33}""", complex=True)
m311,m312,m313, m321,m322,m323, m331,m332,m333 = symbols(\
    r"""m^3_{11},m^3_{12},m^3_{13},m^3_{21},m^3_{22},m^3_{23},m^3_{31},m^3_{32},m^3_{33}""", complex=True)
m411,m412,m413, m421,m422,m423, m431,m432,m433 = symbols(\
    r"""m^4_{11},m^4_{12},m^4_{13},m^4_{21},m^4_{22},m^4_{23},m^4_{31},m^4_{32},m^4_{33}""", complex=True)
m511,m512,m513, m521,m522,m523, m531,m532,m533 = symbols(\
    r"""m^5_{11},m^5_{12},m^5_{13},m^5_{21},m^5_{22},m^5_{23},m^5_{31},m^5_{32},m^5_{33}""", complex=True)
m611,m612,m613, m621,m622,m623, m631,m632,m633 = symbols(\
    r"""m^6_{11},m^6_{12},m^6_{13},m^6_{21},m^6_{22},m^6_{23},m^6_{31},m^6_{32},m^6_{33}""", complex=True)
m711,m712,m713, m721,m722,m723, m731,m732,m733 = symbols(\
    r"""m^7_{11},m^7_{12},m^7_{13},m^7_{21},m^7_{22},m^7_{23},m^7_{31},m^7_{32},m^7_{33}""", complex=True)
m811,m812,m813, m821,m822,m823, m831,m832,m833 = symbols(\
    r"""m^8_{11},m^8_{12},m^8_{13},m^8_{21},m^8_{22},m^8_{23},m^8_{31},m^8_{32},m^8_{33}""", complex=True)

m0 = Matrix([[m011,m012,m013],[m021,m022,m023],[m031,m032,m033]])
m1 = Matrix([[m111,m112,m113],[m121,m122,m123],[m131,m132,m133]])
m2 = Matrix([[m211,m212,m213],[m221,m222,m223],[m231,m232,m233]])
m3 = Matrix([[m311,m312,m313],[m321,m322,m323],[m331,m332,m333]])
m4 = Matrix([[m411,m412,m413],[m421,m422,m423],[m431,m432,m433]])
m5 = Matrix([[m511,m512,m513],[m521,m522,m523],[m531,m532,m533]])
m6 = Matrix([[m611,m612,m613],[m621,m622,m623],[m631,m632,m633]])
m7 = Matrix([[m711,m712,m713],[m721,m722,m723],[m731,m732,m733]])
m8 = Matrix([[m811,m812,m813],[m821,m822,m823],[m831,m832,m833]])


a0 = rep(z0, w0, alpha0, beta0, gamma0, delta0, m0)
a1 = rep(z1, w1, alpha1, beta1, gamma1, delta1, m1)
a2 = rep(z2, w2, alpha2, beta2, gamma2, delta2, m2)
a3 = rep(z3, w3, alpha3, beta3, gamma3, delta3, m3)
a4 = rep(z4, w4, alpha4, beta4, gamma4, delta4, m4)
a5 = rep(z5, w5, alpha5, beta5, gamma5, delta5, m5)
a6 = rep(z6, w6, alpha6, beta6, gamma6, delta6, m6)
a7 = rep(z7, w7, alpha7, beta7, gamma7, delta7, m7)
a8 = rep(z8, w8, alpha8, beta8, gamma8, delta8, m8)

C1,C2,C3,C4,C5,C6,C7,C8,\
C9,C10,C11,C12,C13,C14,C15 = symbols(r"""C_1,C_2,C_3,C_4,C_5,C_6,C_7,C_8,C_9,C_{10},C_{11},C_{12},C_{13},C_{14},C_{15}""", complex = True)
C_m = Matrix([[C7,C8,C9],[C10,C11,C12],[C13,C14,C15]])

general_a = rep(C1,C2,C3,C4,C5,C6,C_m)

test_a1 = rep(0, 0, 1, 0, 0, 0, 0*m1)
test_b1 = rep(0, 0, 0, 1, 0, 0, 0*m1)

test_a2 = rep(0, 0, 0, 1, 0, 0, 0*m1)
test_b2 = rep(0, 1, 0, 0, 0, 0, 0*m1)

minus_simbolovi = [-1*z1,-1*w1,-1*z2,-1*w2,-1*alpha1,-1*alpha2,-1*beta1,-1*beta2,\
             -1*gamma1,-1*gamma2,\
             -1*delta1,-1*delta2,
             -1*m111,-1*m112,-1*m113,-1*m121,-1*m122,-1*m123,-1*m131,-1*m132,-1*m133,\
             -1*m211,-1*m212,-1*m213,-1*m221,-1*m222,-1*m223,-1*m231,-1*m232,-1*m233]


simbolovi0 = [z0, alpha0, beta0, w0, gamma0, delta0,\
              m011,m012,m013, m021,m022,m023, m031,m032,m033]
simbolovi1 = [z1,w1,alpha1,beta1,gamma1,delta1,\
              m111,m112,m113,m121,m122,m123,m131,m132,m133]
simbolovi2 = [z2,w2,alpha2,beta2,gamma2,delta2,\
              m211,m212,m213,m221,m222,m223,m231,m232,m233]
simbolovi3 = [z3,w3,alpha3,beta3,gamma3,delta3,\
              m311,m312,m313,m321,m322,m323,m331,m332,m333]
simbolovi4 = [z4,w4,alpha4,beta4,gamma4,delta4,\
              m411,m412,m413,m421,m422,m423,m431,m432,m433]
simbolovi5 = [z5,w5,alpha5,beta5,gamma5,delta5,\
              m511,m512,m513,m521,m522,m523,m531,m532,m533]
simbolovi6 = [z6,w6,alpha6,beta6,gamma6,delta6,\
              m611,m612,m613,m621,m622,m623,m631,m632,m633]
simbolovi7 = [z7,w7,alpha7,beta7,gamma7,delta7,\
              m711,m712,m713,m721,m722,m723,m731,m732,m733]
simbolovi8 = [z8,w8,alpha8,beta8,gamma8,delta8,\
              m811,m812,m813,m821,m822,m823,m831,m832,m833]

simboloviC = [C1,C2,C3,C4,C5,C6,C7,C8, C9,C10,C11,C12,C13,C14,C15]

simbolovi = simbolovi0 + simbolovi1 + simbolovi2 + simbolovi3 + simbolovi4 + simbolovi5 + simbolovi6 + simbolovi7 + simbolovi8
simbolovi08 = simbolovi0 + simbolovi1 + simbolovi2 + simbolovi3 + simbolovi4 + simbolovi5 + simbolovi6 + simbolovi7 + simbolovi8

da1 = Dirac*a1 - a1*Dirac
da1 = da1.expand()
da1 = da1.simplify()
print(da1.args[0])
print(len(da1.args))

da1 = x*da1
# ---------------------------------------------------------------------------
# LaTeX display
# ---------------------------------------------------------------------------
print("=" * 60)
print("LATEX DISPLAY")
print("=" * 60)

print("\n--- da1 str ---")
print(da1)

print("\n--- da1 latex ---")
print(latex(da1))

print("\n--- da1 _repr_latex_ ---")
print(da1._repr_latex_())

da1.display()
print("da1 args = ")
print(len(da1.args))

print("")
print(srepr(da1))


testoni_100 = AlgebraicPureTensor(2+x,P+Q,Q)
testoni_200 = x* AlgebraicPureTensor(Q,P)
testoni_300 = testoni_100*testoni_200
testoni_300.display()
print(srepr(testoni_100))
print(srepr(testoni_200))
print(testoni_100.expand().simplify())

algebraic_da1 = algebraic_tensor_product(da1, da1)
print(len(algebraic_da1.args))
print(algebraic_da1.args[0])