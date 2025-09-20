from sympy.assumptions.assume import recursive_ask
'\nThis module contains query handlers responsible for Matrices queries:\nSquare, Symmetric, Invertible etc.\n'
from sympy.logic.boolalg import conjuncts
from sympy.assumptions import Q
from sympy.assumptions.handlers import test_closed_group
from sympy.matrices import MatrixBase
from sympy.matrices.expressions import BlockMatrix, BlockDiagMatrix, Determinant, DiagMatrix, DiagonalMatrix, HadamardProduct, Identity, Inverse, MatAdd, MatMul, MatPow, MatrixExpr, MatrixSlice, MatrixSymbol, OneMatrix, Trace, Transpose, ZeroMatrix
from sympy.matrices.expressions.blockmatrix import reblock_2x2
from sympy.matrices.expressions.factorizations import Factorization
from sympy.matrices.expressions.fourier import DFT
from sympy.core.logic import fuzzy_and
from sympy.utilities.iterables import sift
from sympy.core import Basic
from ..predicates.matrices import SquarePredicate, SymmetricPredicate, InvertiblePredicate, OrthogonalPredicate, UnitaryPredicate, FullRankPredicate, PositiveDefinitePredicate, UpperTriangularPredicate, LowerTriangularPredicate, DiagonalPredicate, IntegerElementsPredicate, RealElementsPredicate, ComplexElementsPredicate

def _Factorization(predicate, expr, assumptions):
    if predicate in expr.predicates:
        return True

@SquarePredicate.register(MatrixExpr)
def _(expr, assumptions, rec):
    return expr.shape[0] == expr.shape[1]

@SymmetricPredicate.register(MatMul)
def _(expr, assumptions, rec):
    factor, mmul = expr.as_coeff_mmul()
    if all((recursive_ask(Q.symmetric(arg), assumptions=assumptions, rec=rec) for arg in mmul.args)):
        return True
    if recursive_ask(Q.diagonal(expr), assumptions=assumptions, rec=rec):
        return True
    if len(mmul.args) >= 2 and mmul.args[0] == mmul.args[-1].T:
        if len(mmul.args) == 2:
            return True
        return recursive_ask(Q.symmetric(MatMul(*mmul.args[1:-1])), assumptions=assumptions, rec=rec)

@SymmetricPredicate.register(MatPow)
def _(expr, assumptions, rec):
    base, exp = expr.args
    int_exp = recursive_ask(Q.integer(exp), assumptions=assumptions, rec=rec)
    if not int_exp:
        return None
    non_negative = recursive_ask(~Q.negative(exp), assumptions=assumptions, rec=rec)
    if non_negative or (non_negative == False and recursive_ask(Q.invertible(base), assumptions=assumptions, rec=rec)):
        return recursive_ask(Q.symmetric(base), assumptions=assumptions, rec=rec)
    return None

@SymmetricPredicate.register(MatAdd)
def _(expr, assumptions, rec):
    return all((recursive_ask(Q.symmetric(arg), assumptions=assumptions, rec=rec) for arg in expr.args))

@SymmetricPredicate.register(MatrixSymbol)
def _(expr, assumptions, rec):
    if not expr.is_square:
        return False
    if recursive_ask(Q.diagonal(expr), assumptions=assumptions, rec=rec):
        return True
    if Q.symmetric(expr) in conjuncts(assumptions):
        return True

@SymmetricPredicate.register_many(OneMatrix, ZeroMatrix)
def _(expr, assumptions, rec):
    return recursive_ask(Q.square(expr), assumptions=assumptions, rec=rec)

@SymmetricPredicate.register_many(Inverse, Transpose)
def _(expr, assumptions, rec):
    return recursive_ask(Q.symmetric(expr.arg), assumptions=assumptions, rec=rec)

@SymmetricPredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    if recursive_ask(Q.diagonal(expr), assumptions=assumptions, rec=rec):
        return True
    if not expr.on_diag:
        return None
    else:
        return recursive_ask(Q.symmetric(expr.parent), assumptions=assumptions, rec=rec)

@SymmetricPredicate.register(Identity)
def _(expr, assumptions, rec):
    return True

@InvertiblePredicate.register(MatMul)
def _(expr, assumptions, rec):
    factor, mmul = expr.as_coeff_mmul()
    if all((recursive_ask(Q.invertible(arg), assumptions=assumptions, rec=rec) for arg in mmul.args)):
        return True
    if any((recursive_ask(Q.invertible(arg), assumptions=assumptions, rec=rec) is False for arg in mmul.args)):
        return False

@InvertiblePredicate.register(MatPow)
def _(expr, assumptions, rec):
    base, exp = expr.args
    int_exp = recursive_ask(Q.integer(exp), assumptions=assumptions, rec=rec)
    if not int_exp:
        return None
    if exp.is_negative == False:
        return recursive_ask(Q.invertible(base), assumptions=assumptions, rec=rec)
    return None

@InvertiblePredicate.register(MatAdd)
def _(expr, assumptions, rec):
    return None

@InvertiblePredicate.register(MatrixSymbol)
def _(expr, assumptions, rec):
    if not expr.is_square:
        return False
    if Q.invertible(expr) in conjuncts(assumptions):
        return True

@InvertiblePredicate.register_many(Identity, Inverse)
def _(expr, assumptions, rec):
    return True

@InvertiblePredicate.register(ZeroMatrix)
def _(expr, assumptions, rec):
    return False

@InvertiblePredicate.register(OneMatrix)
def _(expr, assumptions, rec):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@InvertiblePredicate.register(Transpose)
def _(expr, assumptions, rec):
    return recursive_ask(Q.invertible(expr.arg), assumptions=assumptions, rec=rec)

@InvertiblePredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    if not expr.on_diag:
        return None
    else:
        return recursive_ask(Q.invertible(expr.parent), assumptions=assumptions, rec=rec)

@InvertiblePredicate.register(MatrixBase)
def _(expr, assumptions, rec):
    if not expr.is_square:
        return False
    return expr.rank() == expr.rows

@InvertiblePredicate.register(MatrixExpr)
def _(expr, assumptions, rec):
    if not expr.is_square:
        return False
    return None

@InvertiblePredicate.register(BlockMatrix)
def _(expr, assumptions, rec):
    if not expr.is_square:
        return False
    if expr.blockshape == (1, 1):
        return recursive_ask(Q.invertible(expr.blocks[0, 0]), assumptions=assumptions, rec=rec)
    expr = reblock_2x2(expr)
    if expr.blockshape == (2, 2):
        [[A, B], [C, D]] = expr.blocks.tolist()
        if recursive_ask(Q.invertible(A), assumptions=assumptions, rec=rec) == True:
            invertible = recursive_ask(Q.invertible(D - C * A.I * B), assumptions=assumptions, rec=rec)
            if invertible is not None:
                return invertible
        if recursive_ask(Q.invertible(B), assumptions=assumptions, rec=rec) == True:
            invertible = recursive_ask(Q.invertible(C - D * B.I * A), assumptions=assumptions, rec=rec)
            if invertible is not None:
                return invertible
        if recursive_ask(Q.invertible(C), assumptions=assumptions, rec=rec) == True:
            invertible = recursive_ask(Q.invertible(B - A * C.I * D), assumptions=assumptions, rec=rec)
            if invertible is not None:
                return invertible
        if recursive_ask(Q.invertible(D), assumptions=assumptions, rec=rec) == True:
            invertible = recursive_ask(Q.invertible(A - B * D.I * C), assumptions=assumptions, rec=rec)
            if invertible is not None:
                return invertible
    return None

@InvertiblePredicate.register(BlockDiagMatrix)
def _(expr, assumptions, rec):
    if expr.rowblocksizes != expr.colblocksizes:
        return None
    return fuzzy_and([recursive_ask(Q.invertible(a), assumptions=assumptions, rec=rec) for a in expr.diag])

@OrthogonalPredicate.register(MatMul)
def _(expr, assumptions, rec):
    factor, mmul = expr.as_coeff_mmul()
    if all((recursive_ask(Q.orthogonal(arg), assumptions=assumptions, rec=rec) for arg in mmul.args)) and factor == 1:
        return True
    if any((recursive_ask(Q.invertible(arg), assumptions=assumptions, rec=rec) is False for arg in mmul.args)):
        return False

@OrthogonalPredicate.register(MatPow)
def _(expr, assumptions, rec):
    base, exp = expr.args
    int_exp = recursive_ask(Q.integer(exp), assumptions=assumptions, rec=rec)
    if int_exp:
        return recursive_ask(Q.orthogonal(base), assumptions=assumptions, rec=rec)
    return None

@OrthogonalPredicate.register(MatAdd)
def _(expr, assumptions, rec):
    if len(expr.args) == 1 and recursive_ask(Q.orthogonal(expr.args[0]), assumptions=assumptions, rec=rec):
        return True

@OrthogonalPredicate.register(MatrixSymbol)
def _(expr, assumptions, rec):
    if not expr.is_square or recursive_ask(Q.invertible(expr), assumptions=assumptions, rec=rec) is False:
        return False
    if Q.orthogonal(expr) in conjuncts(assumptions):
        return True

@OrthogonalPredicate.register(Identity)
def _(expr, assumptions, rec):
    return True

@OrthogonalPredicate.register(ZeroMatrix)
def _(expr, assumptions, rec):
    return False

@OrthogonalPredicate.register_many(Inverse, Transpose)
def _(expr, assumptions, rec):
    return recursive_ask(Q.orthogonal(expr.arg), assumptions=assumptions, rec=rec)

@OrthogonalPredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    if not expr.on_diag:
        return None
    else:
        return recursive_ask(Q.orthogonal(expr.parent), assumptions=assumptions, rec=rec)

@OrthogonalPredicate.register(Factorization)
def _(expr, assumptions, rec):
    return _Factorization(Q.orthogonal, expr, assumptions)

@UnitaryPredicate.register(MatMul)
def _(expr, assumptions, rec):
    factor, mmul = expr.as_coeff_mmul()
    if all((recursive_ask(Q.unitary(arg), assumptions=assumptions, rec=rec) for arg in mmul.args)) and abs(factor) == 1:
        return True
    if any((recursive_ask(Q.invertible(arg), assumptions=assumptions, rec=rec) is False for arg in mmul.args)):
        return False

@UnitaryPredicate.register(MatPow)
def _(expr, assumptions, rec):
    base, exp = expr.args
    int_exp = recursive_ask(Q.integer(exp), assumptions=assumptions, rec=rec)
    if int_exp:
        return recursive_ask(Q.unitary(base), assumptions=assumptions, rec=rec)
    return None

@UnitaryPredicate.register(MatrixSymbol)
def _(expr, assumptions, rec):
    if not expr.is_square or recursive_ask(Q.invertible(expr), assumptions=assumptions, rec=rec) is False:
        return False
    if Q.unitary(expr) in conjuncts(assumptions):
        return True

@UnitaryPredicate.register_many(Inverse, Transpose)
def _(expr, assumptions, rec):
    return recursive_ask(Q.unitary(expr.arg), assumptions=assumptions, rec=rec)

@UnitaryPredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    if not expr.on_diag:
        return None
    else:
        return recursive_ask(Q.unitary(expr.parent), assumptions=assumptions, rec=rec)

@UnitaryPredicate.register_many(DFT, Identity)
def _(expr, assumptions, rec):
    return True

@UnitaryPredicate.register(ZeroMatrix)
def _(expr, assumptions, rec):
    return False

@UnitaryPredicate.register(Factorization)
def _(expr, assumptions, rec):
    return _Factorization(Q.unitary, expr, assumptions)

@FullRankPredicate.register(MatMul)
def _(expr, assumptions, rec):
    if all((recursive_ask(Q.fullrank(arg), assumptions=assumptions, rec=rec) for arg in expr.args)):
        return True

@FullRankPredicate.register(MatPow)
def _(expr, assumptions, rec):
    base, exp = expr.args
    int_exp = recursive_ask(Q.integer(exp), assumptions=assumptions, rec=rec)
    if int_exp and recursive_ask(~Q.negative(exp), assumptions=assumptions, rec=rec):
        return recursive_ask(Q.fullrank(base), assumptions=assumptions, rec=rec)
    return None

@FullRankPredicate.register(Identity)
def _(expr, assumptions, rec):
    return True

@FullRankPredicate.register(ZeroMatrix)
def _(expr, assumptions, rec):
    return False

@FullRankPredicate.register(OneMatrix)
def _(expr, assumptions, rec):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@FullRankPredicate.register_many(Inverse, Transpose)
def _(expr, assumptions, rec):
    return recursive_ask(Q.fullrank(expr.arg), assumptions=assumptions, rec=rec)

@FullRankPredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    if recursive_ask(Q.orthogonal(expr.parent), assumptions=assumptions, rec=rec):
        return True

@PositiveDefinitePredicate.register(MatMul)
def _(expr, assumptions, rec):
    factor, mmul = expr.as_coeff_mmul()
    if all((recursive_ask(Q.positive_definite(arg), assumptions=assumptions, rec=rec) for arg in mmul.args)) and factor > 0:
        return True
    if len(mmul.args) >= 2 and mmul.args[0] == mmul.args[-1].T and recursive_ask(Q.fullrank(mmul.args[0]), assumptions=assumptions, rec=rec):
        return recursive_ask(Q.positive_definite(MatMul(*mmul.args[1:-1])), assumptions=assumptions, rec=rec)

@PositiveDefinitePredicate.register(MatPow)
def _(expr, assumptions, rec):
    if recursive_ask(Q.positive_definite(expr.args[0]), assumptions=assumptions, rec=rec):
        return True

@PositiveDefinitePredicate.register(MatAdd)
def _(expr, assumptions, rec):
    if all((recursive_ask(Q.positive_definite(arg), assumptions=assumptions, rec=rec) for arg in expr.args)):
        return True

@PositiveDefinitePredicate.register(MatrixSymbol)
def _(expr, assumptions, rec):
    if not expr.is_square:
        return False
    if Q.positive_definite(expr) in conjuncts(assumptions):
        return True

@PositiveDefinitePredicate.register(Identity)
def _(expr, assumptions, rec):
    return True

@PositiveDefinitePredicate.register(ZeroMatrix)
def _(expr, assumptions, rec):
    return False

@PositiveDefinitePredicate.register(OneMatrix)
def _(expr, assumptions, rec):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@PositiveDefinitePredicate.register_many(Inverse, Transpose)
def _(expr, assumptions, rec):
    return recursive_ask(Q.positive_definite(expr.arg), assumptions=assumptions, rec=rec)

@PositiveDefinitePredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    if not expr.on_diag:
        return None
    else:
        return recursive_ask(Q.positive_definite(expr.parent), assumptions=assumptions, rec=rec)

@UpperTriangularPredicate.register(MatMul)
def _(expr, assumptions, rec):
    factor, matrices = expr.as_coeff_matrices()
    if all((recursive_ask(Q.upper_triangular(m), assumptions=assumptions, rec=rec) for m in matrices)):
        return True

@UpperTriangularPredicate.register(MatAdd)
def _(expr, assumptions, rec):
    if all((recursive_ask(Q.upper_triangular(arg), assumptions=assumptions, rec=rec) for arg in expr.args)):
        return True

@UpperTriangularPredicate.register(MatPow)
def _(expr, assumptions, rec):
    base, exp = expr.args
    int_exp = recursive_ask(Q.integer(exp), assumptions=assumptions, rec=rec)
    if not int_exp:
        return None
    non_negative = recursive_ask(~Q.negative(exp), assumptions=assumptions, rec=rec)
    if non_negative or (non_negative == False and recursive_ask(Q.invertible(base), assumptions=assumptions, rec=rec)):
        return recursive_ask(Q.upper_triangular(base), assumptions=assumptions, rec=rec)
    return None

@UpperTriangularPredicate.register(MatrixSymbol)
def _(expr, assumptions, rec):
    if Q.upper_triangular(expr) in conjuncts(assumptions):
        return True

@UpperTriangularPredicate.register_many(Identity, ZeroMatrix)
def _(expr, assumptions, rec):
    return True

@UpperTriangularPredicate.register(OneMatrix)
def _(expr, assumptions, rec):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@UpperTriangularPredicate.register(Transpose)
def _(expr, assumptions, rec):
    return recursive_ask(Q.lower_triangular(expr.arg), assumptions=assumptions, rec=rec)

@UpperTriangularPredicate.register(Inverse)
def _(expr, assumptions, rec):
    return recursive_ask(Q.upper_triangular(expr.arg), assumptions=assumptions, rec=rec)

@UpperTriangularPredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    if not expr.on_diag:
        return None
    else:
        return recursive_ask(Q.upper_triangular(expr.parent), assumptions=assumptions, rec=rec)

@UpperTriangularPredicate.register(Factorization)
def _(expr, assumptions, rec):
    return _Factorization(Q.upper_triangular, expr, assumptions)

@LowerTriangularPredicate.register(MatMul)
def _(expr, assumptions, rec):
    factor, matrices = expr.as_coeff_matrices()
    if all((recursive_ask(Q.lower_triangular(m), assumptions=assumptions, rec=rec) for m in matrices)):
        return True

@LowerTriangularPredicate.register(MatAdd)
def _(expr, assumptions, rec):
    if all((recursive_ask(Q.lower_triangular(arg), assumptions=assumptions, rec=rec) for arg in expr.args)):
        return True

@LowerTriangularPredicate.register(MatPow)
def _(expr, assumptions, rec):
    base, exp = expr.args
    int_exp = recursive_ask(Q.integer(exp), assumptions=assumptions, rec=rec)
    if not int_exp:
        return None
    non_negative = recursive_ask(~Q.negative(exp), assumptions=assumptions, rec=rec)
    if non_negative or (non_negative == False and recursive_ask(Q.invertible(base), assumptions=assumptions, rec=rec)):
        return recursive_ask(Q.lower_triangular(base), assumptions=assumptions, rec=rec)
    return None

@LowerTriangularPredicate.register(MatrixSymbol)
def _(expr, assumptions, rec):
    if Q.lower_triangular(expr) in conjuncts(assumptions):
        return True

@LowerTriangularPredicate.register_many(Identity, ZeroMatrix)
def _(expr, assumptions, rec):
    return True

@LowerTriangularPredicate.register(OneMatrix)
def _(expr, assumptions, rec):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@LowerTriangularPredicate.register(Transpose)
def _(expr, assumptions, rec):
    return recursive_ask(Q.upper_triangular(expr.arg), assumptions=assumptions, rec=rec)

@LowerTriangularPredicate.register(Inverse)
def _(expr, assumptions, rec):
    return recursive_ask(Q.lower_triangular(expr.arg), assumptions=assumptions, rec=rec)

@LowerTriangularPredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    if not expr.on_diag:
        return None
    else:
        return recursive_ask(Q.lower_triangular(expr.parent), assumptions=assumptions, rec=rec)

@LowerTriangularPredicate.register(Factorization)
def _(expr, assumptions, rec):
    return _Factorization(Q.lower_triangular, expr, assumptions)

def _is_empty_or_1x1(expr):
    return expr.shape in ((0, 0), (1, 1))

@DiagonalPredicate.register(MatMul)
def _(expr, assumptions, rec):
    if _is_empty_or_1x1(expr):
        return True
    factor, matrices = expr.as_coeff_matrices()
    if all((recursive_ask(Q.diagonal(m), assumptions=assumptions, rec=rec) for m in matrices)):
        return True

@DiagonalPredicate.register(MatPow)
def _(expr, assumptions, rec):
    base, exp = expr.args
    int_exp = recursive_ask(Q.integer(exp), assumptions=assumptions, rec=rec)
    if not int_exp:
        return None
    non_negative = recursive_ask(~Q.negative(exp), assumptions=assumptions, rec=rec)
    if non_negative or (non_negative == False and recursive_ask(Q.invertible(base), assumptions=assumptions, rec=rec)):
        return recursive_ask(Q.diagonal(base), assumptions=assumptions, rec=rec)
    return None

@DiagonalPredicate.register(MatAdd)
def _(expr, assumptions, rec):
    if all((recursive_ask(Q.diagonal(arg), assumptions=assumptions, rec=rec) for arg in expr.args)):
        return True

@DiagonalPredicate.register(MatrixSymbol)
def _(expr, assumptions, rec):
    if _is_empty_or_1x1(expr):
        return True
    if Q.diagonal(expr) in conjuncts(assumptions):
        return True

@DiagonalPredicate.register(OneMatrix)
def _(expr, assumptions, rec):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@DiagonalPredicate.register_many(Inverse, Transpose)
def _(expr, assumptions, rec):
    return recursive_ask(Q.diagonal(expr.arg), assumptions=assumptions, rec=rec)

@DiagonalPredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    if _is_empty_or_1x1(expr):
        return True
    if not expr.on_diag:
        return None
    else:
        return recursive_ask(Q.diagonal(expr.parent), assumptions=assumptions, rec=rec)

@DiagonalPredicate.register_many(DiagonalMatrix, DiagMatrix, Identity, ZeroMatrix)
def _(expr, assumptions, rec):
    return True

@DiagonalPredicate.register(Factorization)
def _(expr, assumptions, rec):
    return _Factorization(Q.diagonal, expr, assumptions)

def BM_elements(predicate, expr, assumptions, rec):
    """ Block Matrix elements. """
    return all((recursive_ask(predicate(b), assumptions=assumptions, rec=rec) for b in expr.blocks))

def MS_elements(predicate, expr, assumptions, rec):
    """ Matrix Slice elements. """
    return recursive_ask(predicate(expr.parent), assumptions=assumptions, rec=rec)

def MatMul_elements(matrix_predicate, scalar_predicate, expr, assumptions, rec):
    d = sift(expr.args, lambda x: isinstance(x, MatrixExpr))
    factors, matrices = (d[False], d[True])
    return fuzzy_and([test_closed_group(Basic(*factors), assumptions, scalar_predicate, rec=rec), test_closed_group(Basic(*matrices), assumptions, matrix_predicate, rec=rec)])

@IntegerElementsPredicate.register_many(Determinant, HadamardProduct, MatAdd, Trace, Transpose)
def _(expr, assumptions, rec):
    return test_closed_group(expr, assumptions, Q.integer_elements, rec=rec)

@IntegerElementsPredicate.register(MatPow)
def _(expr, assumptions, rec):
    base, exp = expr.args
    int_exp = recursive_ask(Q.integer(exp), assumptions=assumptions, rec=rec)
    if not int_exp:
        return None
    if exp.is_negative == False:
        return recursive_ask(Q.integer_elements(base), assumptions=assumptions, rec=rec)
    return None

@IntegerElementsPredicate.register_many(Identity, OneMatrix, ZeroMatrix)
def _(expr, assumptions, rec):
    return True

@IntegerElementsPredicate.register(MatMul)
def _(expr, assumptions, rec):
    return MatMul_elements(Q.integer_elements, Q.integer, expr, assumptions, rec=rec)

@IntegerElementsPredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    return MS_elements(Q.integer_elements, expr, assumptions, rec=rec)

@IntegerElementsPredicate.register(BlockMatrix)
def _(expr, assumptions, rec):
    return BM_elements(Q.integer_elements, expr, assumptions, rec=rec)

@RealElementsPredicate.register_many(Determinant, Factorization, HadamardProduct, MatAdd, Trace, Transpose)
def _(expr, assumptions, rec):
    return test_closed_group(expr, assumptions, Q.real_elements, rec=rec)

@RealElementsPredicate.register(MatPow)
def _(expr, assumptions, rec):
    base, exp = expr.args
    int_exp = recursive_ask(Q.integer(exp), assumptions=assumptions, rec=rec)
    if not int_exp:
        return None
    non_negative = recursive_ask(~Q.negative(exp), assumptions=assumptions, rec=rec)
    if non_negative or (non_negative == False and recursive_ask(Q.invertible(base), assumptions=assumptions, rec=rec)):
        return recursive_ask(Q.real_elements(base), assumptions=assumptions, rec=rec)
    return None

@RealElementsPredicate.register(MatMul)
def _(expr, assumptions, rec):
    return MatMul_elements(Q.real_elements, Q.real, expr, assumptions, rec=rec)

@RealElementsPredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    return MS_elements(Q.real_elements, expr, assumptions, rec=rec)

@RealElementsPredicate.register(BlockMatrix)
def _(expr, assumptions, rec):
    return BM_elements(Q.real_elements, expr, assumptions, rec=rec)

@ComplexElementsPredicate.register_many(Determinant, Factorization, HadamardProduct, Inverse, MatAdd, Trace, Transpose)
def _(expr, assumptions, rec):
    return test_closed_group(expr, assumptions, Q.complex_elements, rec=rec)

@ComplexElementsPredicate.register(MatPow)
def _(expr, assumptions, rec):
    base, exp = expr.args
    int_exp = recursive_ask(Q.integer(exp), assumptions=assumptions, rec=rec)
    if not int_exp:
        return None
    non_negative = recursive_ask(~Q.negative(exp), assumptions=assumptions, rec=rec)
    if non_negative or (non_negative == False and recursive_ask(Q.invertible(base), assumptions=assumptions, rec=rec)):
        return recursive_ask(Q.complex_elements(base), assumptions=assumptions, rec=rec)
    return None

@ComplexElementsPredicate.register(MatMul)
def _(expr, assumptions, rec):
    return MatMul_elements(Q.complex_elements, Q.complex, expr, assumptions, rec=rec)

@ComplexElementsPredicate.register(MatrixSlice)
def _(expr, assumptions, rec):
    return MS_elements(Q.complex_elements, expr, assumptions, rec=rec)

@ComplexElementsPredicate.register(BlockMatrix)
def _(expr, assumptions, rec):
    return BM_elements(Q.complex_elements, expr, assumptions, rec=rec)

@ComplexElementsPredicate.register(DFT)
def _(expr, assumptions, rec):
    return True