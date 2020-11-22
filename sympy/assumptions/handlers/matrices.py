"""
This module contains query handlers responsible for calculus queries:
infinitesimal, bounded, etc.
"""

from sympy.logic.boolalg import conjuncts
from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import (
    AskHandlerClass, CommonHandler, test_closed_group
)
from sympy.matrices import (
    MatrixExpr, MatMul, MatPow, MatAdd, MatrixSymbol, ZeroMatrix, OneMatrix,
    Transpose, Inverse, MatrixSlice, Identity,
    MatrixBase, BlockMatrix, BlockDiagMatrix,
    DiagonalMatrix, DiagMatrix,
    HadamardProduct, Determinant, Trace
)
from sympy.matrices.expressions.factorizations import Factorization
from sympy.matrices.expressions.fourier import DFT
from sympy.core.logic import fuzzy_and
from sympy.utilities.iterables import sift
from sympy.core import Basic
from functools import partial


### Helper functions ###

def _Factorization(predicate, expr, assumptions):
    if predicate in expr.predicates:
        return True

def BM_elements(predicate, expr, assumptions):
    """ Block Matrix elements """
    return all(ask(predicate(b), assumptions) for b in expr.blocks)

def MS_elements(predicate, expr, assumptions):
    """ Matrix Slice elements """
    return ask(predicate(expr.parent), assumptions)

def MatMul_elements(matrix_predicate, scalar_predicate, expr, assumptions):
    d = sift(expr.args, lambda x: isinstance(x, MatrixExpr))
    factors, matrices = d[False], d[True]
    return fuzzy_and([
        test_closed_group(Basic(*factors), assumptions, scalar_predicate),
        test_closed_group(Basic(*matrices), assumptions, matrix_predicate)])


### AskSquareHandler ###

AskSquareHandler = AskHandlerClass(
    'AskSquareHandler',
    doc="""
    Handler for key 'square'
    """,
    base=CommonHandler
)

@AskSquareHandler.register(MatrixExpr)
def _(expr, assumptions):
    return expr.shape[0] == expr.shape[1]


### AskSymmetricHandler ###

AskSymmetricHandler = AskHandlerClass(
    'AskSymmetricHandler',
    doc="""
    Handler for key 'symmetric'
    """,
    base=CommonHandler
)

AskSymmetricHandler.register(Identity)(AskSymmetricHandler.AlwaysTrue)

@AskSymmetricHandler.register(MatMul)
def _(expr, assumptions):
    factor, mmul = expr.as_coeff_mmul()
    if all(ask(Q.symmetric(arg), assumptions) for arg in mmul.args):
        return True
    # TODO: implement sathandlers system for the matrices.
    # Now it duplicates the general fact: Implies(Q.diagonal, Q.symmetric).
    if ask(Q.diagonal(expr), assumptions):
        return True
    if len(mmul.args) >= 2 and mmul.args[0] == mmul.args[-1].T:
        if len(mmul.args) == 2:
            return True
        return ask(Q.symmetric(MatMul(*mmul.args[1:-1])), assumptions)

@AskSymmetricHandler.register(MatPow)
def _(expr, assumptions):
    # only for integer powers
    base, exp = expr.args
    int_exp = ask(Q.integer(exp), assumptions)
    if not int_exp:
        return None
    non_negative = ask(~Q.negative(exp), assumptions)
    if (non_negative or non_negative == False
                        and ask(Q.invertible(base), assumptions)):
        return ask(Q.symmetric(base), assumptions)
    return None

@AskSymmetricHandler.register(MatAdd)
def _(expr, assumptions):
    return all(ask(Q.symmetric(arg), assumptions) for arg in expr.args)

@AskSymmetricHandler.register(MatrixSymbol)
def _(expr, assumptions):
    if not expr.is_square:
        return False
    # TODO: implement sathandlers system for the matrices.
    # Now it duplicates the general fact: Implies(Q.diagonal, Q.symmetric).
    if ask(Q.diagonal(expr), assumptions):
        return True
    if Q.symmetric(expr) in conjuncts(assumptions):
        return True

@AskSymmetricHandler.register(MatrixSlice)
def _(expr, assumptions):
    # TODO: implement sathandlers system for the matrices.
    # Now it duplicates the general fact: Implies(Q.diagonal, Q.symmetric).
    if ask(Q.diagonal(expr), assumptions):
        return True
    if not expr.on_diag:
        return None
    else:
        return ask(Q.symmetric(expr.parent), assumptions)

@AskSymmetricHandler.register(ZeroMatrix)
@AskSymmetricHandler.register(OneMatrix)
def _(expr, assumptions):
    return ask(Q.square(expr), assumptions)

@AskSymmetricHandler.register(Transpose)
@AskSymmetricHandler.register(Inverse)
def _(expr, assumptions):
    return ask(Q.symmetric(expr.arg), assumptions)


### AskInvertibleHandler ###

AskInvertibleHandler = AskHandlerClass(
    'AskInvertibleHandler',
    doc="""
    Handler for key 'invertible'
    """,
    base=CommonHandler
)

AskInvertibleHandler.register(ZeroMatrix)(AskInvertibleHandler.AlwaysFalse)

for sig in (Identity, Inverse):
    AskInvertibleHandler.register(sig)(AskInvertibleHandler.AlwaysTrue)

@AskInvertibleHandler.register(MatMul)
def _(expr, assumptions):
    factor, mmul = expr.as_coeff_mmul()
    if all(ask(Q.invertible(arg), assumptions) for arg in mmul.args):
        return True
    if any(ask(Q.invertible(arg), assumptions) is False
            for arg in mmul.args):
        return False

@AskInvertibleHandler.register(MatPow)
def _(expr, assumptions):
    # only for integer powers
    base, exp = expr.args
    int_exp = ask(Q.integer(exp), assumptions)
    if not int_exp:
        return None
    if exp.is_negative == False:
        return ask(Q.invertible(base), assumptions)
    return None

@AskInvertibleHandler.register(MatAdd)
def _(expr, assumptions):
    return None

@AskInvertibleHandler.register(MatrixSymbol)
def _(expr, assumptions):
    if not expr.is_square:
        return False
    if Q.invertible(expr) in conjuncts(assumptions):
        return True

@AskInvertibleHandler.register(OneMatrix)
def _(expr, assumptions):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@AskInvertibleHandler.register(Transpose)
def _(expr, assumptions):
    return ask(Q.invertible(expr.arg), assumptions)

@AskInvertibleHandler.register(MatrixSlice)
def _(expr, assumptions):
    if not expr.on_diag:
        return None
    else:
        return ask(Q.invertible(expr.parent), assumptions)

@AskInvertibleHandler.register(MatrixBase)
def _(expr, assumptions):
    if not expr.is_square:
        return False
    return expr.rank() == expr.rows

@AskInvertibleHandler.register(MatrixExpr)
def _(expr, assumptions):
    if not expr.is_square:
        return False
    return None

@AskInvertibleHandler.register(BlockMatrix)
def _(expr, assumptions):
    from sympy.matrices.expressions.blockmatrix import reblock_2x2
    if not expr.is_square:
        return False
    if expr.blockshape == (1, 1):
        return ask(Q.invertible(expr.blocks[0, 0]), assumptions)
    expr = reblock_2x2(expr)
    if expr.blockshape == (2, 2):
        [[A, B], [C, D]] = expr.blocks.tolist()
        if ask(Q.invertible(A), assumptions) == True:
            invertible = ask(Q.invertible(D - C * A.I * B), assumptions)
            if invertible is not None:
                return invertible
        if ask(Q.invertible(B), assumptions) == True:
            invertible = ask(Q.invertible(C - D * B.I * A), assumptions)
            if invertible is not None:
                return invertible
        if ask(Q.invertible(C), assumptions) == True:
            invertible = ask(Q.invertible(B - A * C.I * D), assumptions)
            if invertible is not None:
                return invertible
        if ask(Q.invertible(D), assumptions) == True:
            invertible = ask(Q.invertible(A - B * D.I * C), assumptions)
            if invertible is not None:
                return invertible
    return None

@AskInvertibleHandler.register(BlockDiagMatrix)
def _(expr, assumptions):
    if expr.rowblocksizes != expr.colblocksizes:
        return None
    return fuzzy_and([ask(Q.invertible(a), assumptions) for a in expr.diag])


### AskOrthogonalHandler ###

AskOrthogonalHandler = AskHandlerClass(
    'AskOrthogonalHandler',
    doc="""
    Handler for key 'orthogonal'
    """,
    base=CommonHandler
)

AskOrthogonalHandler.register(Identity)(AskOrthogonalHandler.AlwaysTrue)

AskOrthogonalHandler.register(ZeroMatrix)(AskOrthogonalHandler.AlwaysFalse)

AskOrthogonalHandler.register(Factorization)(partial(_Factorization, Q.orthogonal))

@AskOrthogonalHandler.register(MatMul)
def _(expr, assumptions):
    factor, mmul = expr.as_coeff_mmul()
    if (all(ask(Q.orthogonal(arg), assumptions) for arg in mmul.args) and
            factor == 1):
        return True
    if any(ask(Q.invertible(arg), assumptions) is False
            for arg in mmul.args):
        return False

@AskOrthogonalHandler.register(MatPow)
def _(expr, assumptions):
    # only for integer powers
    base, exp = expr.args
    int_exp = ask(Q.integer(exp), assumptions)
    if int_exp:
        return ask(Q.orthogonal(base), assumptions)
    return None

@AskOrthogonalHandler.register(MatAdd)
def _(expr, assumptions):
    if (len(expr.args) == 1 and
            ask(Q.orthogonal(expr.args[0]), assumptions)):
        return True

@AskOrthogonalHandler.register(MatrixSymbol)
def _(expr, assumptions):
    if (not expr.is_square or
                    ask(Q.invertible(expr), assumptions) is False):
        return False
    if Q.orthogonal(expr) in conjuncts(assumptions):
        return True

@AskOrthogonalHandler.register(MatrixSlice)
def _(expr, assumptions):
    if not expr.on_diag:
        return None
    else:
        return ask(Q.orthogonal(expr.parent), assumptions)

@AskOrthogonalHandler.register(Transpose)
@AskOrthogonalHandler.register(Inverse)
def _(expr, assumptions):
    return ask(Q.orthogonal(expr.arg), assumptions)


### AskUnitaryHandler ###

AskUnitaryHandler = AskHandlerClass(
    'AskUnitaryHandler',
    doc="""
    Handler for key 'unitary'
    """,
    base=CommonHandler
)

for sig in (DFT, Identity):
    AskUnitaryHandler.register(sig)(AskUnitaryHandler.AlwaysTrue)

AskUnitaryHandler.register(ZeroMatrix)(AskUnitaryHandler.AlwaysFalse)

AskUnitaryHandler.register(Factorization)(partial(_Factorization, Q.unitary))

@AskUnitaryHandler.register(MatMul)
def _(expr, assumptions):
    factor, mmul = expr.as_coeff_mmul()
    if (all(ask(Q.unitary(arg), assumptions) for arg in mmul.args) and
            abs(factor) == 1):
        return True
    if any(ask(Q.invertible(arg), assumptions) is False
            for arg in mmul.args):
        return False

@AskUnitaryHandler.register(MatPow)
def _(expr, assumptions):
    # only for integer powers
    base, exp = expr.args
    int_exp = ask(Q.integer(exp), assumptions)
    if int_exp:
        return ask(Q.unitary(base), assumptions)
    return None

@AskUnitaryHandler.register(MatrixSymbol)
def _(expr, assumptions):
    if (not expr.is_square or
                    ask(Q.invertible(expr), assumptions) is False):
        return False
    if Q.unitary(expr) in conjuncts(assumptions):
        return True

@AskUnitaryHandler.register(MatrixSlice)
def _(expr, assumptions):
    if not expr.on_diag:
        return None
    else:
        return ask(Q.unitary(expr.parent), assumptions)

@AskUnitaryHandler.register(Transpose)
@AskUnitaryHandler.register(Inverse)
def _(expr, assumptions):
    return ask(Q.unitary(expr.arg), assumptions)


### AskFullRankHandler ###

AskFullRankHandler = AskHandlerClass(
    'AskFullRankHandler',
    doc="""
    Handler for key 'fullrank'
    """,
    base=CommonHandler
)

AskFullRankHandler.register(Identity)(AskFullRankHandler.AlwaysTrue)

AskFullRankHandler.register(ZeroMatrix)(AskFullRankHandler.AlwaysFalse)

@AskFullRankHandler.register(MatMul)
def _(expr, assumptions):
    if all(ask(Q.fullrank(arg), assumptions) for arg in expr.args):
        return True

@AskFullRankHandler.register(MatPow)
def _(expr, assumptions):
    # only for integer powers
    base, exp = expr.args
    int_exp = ask(Q.integer(exp), assumptions)
    if int_exp and ask(~Q.negative(exp), assumptions):
        return ask(Q.fullrank(base), assumptions)
    return None

@AskFullRankHandler.register(OneMatrix)
def _(expr, assumptions):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@AskFullRankHandler.register(MatrixSlice)
def _(expr, assumptions):
    if ask(Q.orthogonal(expr.parent), assumptions):
        return True

@AskFullRankHandler.register(Transpose)
@AskFullRankHandler.register(Inverse)
def _(expr, assumptions):
    return ask(Q.fullrank(expr.arg), assumptions)


### AskPositiveDefiniteHandler ###

AskPositiveDefiniteHandler = AskHandlerClass(
    'AskPositiveDefiniteHandler',
    doc="""
    Handler for key 'positive_definite'
    """,
    base=CommonHandler
)

AskPositiveDefiniteHandler.register(Identity)(AskPositiveDefiniteHandler.AlwaysTrue)

AskPositiveDefiniteHandler.register(ZeroMatrix)(AskPositiveDefiniteHandler.AlwaysFalse)

@AskPositiveDefiniteHandler.register(MatMul)
def _(expr, assumptions):
    factor, mmul = expr.as_coeff_mmul()
    if (all(ask(Q.positive_definite(arg), assumptions)
            for arg in mmul.args) and factor > 0):
        return True
    if (len(mmul.args) >= 2
            and mmul.args[0] == mmul.args[-1].T
            and ask(Q.fullrank(mmul.args[0]), assumptions)):
        return ask(Q.positive_definite(
            MatMul(*mmul.args[1:-1])), assumptions)

@AskPositiveDefiniteHandler.register(MatPow)
def _(expr, assumptions):
    # a power of a positive definite matrix is positive definite
    if ask(Q.positive_definite(expr.args[0]), assumptions):
        return True

@AskPositiveDefiniteHandler.register(MatAdd)
def _(expr, assumptions):
    if all(ask(Q.positive_definite(arg), assumptions)
            for arg in expr.args):
        return True

@AskPositiveDefiniteHandler.register(MatrixSymbol)
def _(expr, assumptions):
    if not expr.is_square:
        return False
    if Q.positive_definite(expr) in conjuncts(assumptions):
        return True

@AskPositiveDefiniteHandler.register(OneMatrix)
def _(expr, assumptions):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@AskPositiveDefiniteHandler.register(MatrixSlice)
def _(expr, assumptions):
    if not expr.on_diag:
        return None
    else:
        return ask(Q.positive_definite(expr.parent), assumptions)

@AskPositiveDefiniteHandler.register(Transpose)
@AskPositiveDefiniteHandler.register(Inverse)
def _(expr, assumptions):
    return ask(Q.positive_definite(expr.arg), assumptions)


### AskUpperTriangularHandler ###

AskUpperTriangularHandler = AskHandlerClass(
    'AskUpperTriangularHandler',
    doc="""
    Handler for key 'upper_triangular'
    """,
    base=CommonHandler
)

for sig in (Identity, ZeroMatrix):
    AskUpperTriangularHandler.register(sig)(AskUpperTriangularHandler.AlwaysTrue)

AskUpperTriangularHandler.register(Factorization)(partial(_Factorization, Q.upper_triangular))

@AskUpperTriangularHandler.register(MatMul)
def _(expr, assumptions):
    factor, matrices = expr.as_coeff_matrices()
    if all(ask(Q.upper_triangular(m), assumptions) for m in matrices):
        return True

@AskUpperTriangularHandler.register(MatAdd)
def _(expr, assumptions):
    if all(ask(Q.upper_triangular(arg), assumptions) for arg in expr.args):
        return True

@AskUpperTriangularHandler.register(MatPow)
def _(expr, assumptions):
    # only for integer powers
    base, exp = expr.args
    int_exp = ask(Q.integer(exp), assumptions)
    if not int_exp:
        return None
    non_negative = ask(~Q.negative(exp), assumptions)
    if (non_negative or non_negative == False
                        and ask(Q.invertible(base), assumptions)):
        return ask(Q.upper_triangular(base), assumptions)
    return None

@AskUpperTriangularHandler.register(MatrixSymbol)
def _(expr, assumptions):
    if Q.upper_triangular(expr) in conjuncts(assumptions):
        return True

@AskUpperTriangularHandler.register(OneMatrix)
def _(expr, assumptions):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@AskUpperTriangularHandler.register(Transpose)
def _(expr, assumptions):
    return ask(Q.lower_triangular(expr.arg), assumptions)

@AskUpperTriangularHandler.register(Inverse)
def _(expr, assumptions):
    return ask(Q.upper_triangular(expr.arg), assumptions)

@AskUpperTriangularHandler.register(MatrixSlice)
def _(expr, assumptions):
    if not expr.on_diag:
        return None
    else:
        return ask(Q.upper_triangular(expr.parent), assumptions)


### AskLowerTriangularHandler ###

AskLowerTriangularHandler = AskHandlerClass(
    'AskLowerTriangularHandler',
    doc="""
    Handler for key 'lower_triangular'
    """,
    base=CommonHandler
)

for sig in (Identity, ZeroMatrix):
    AskLowerTriangularHandler.register(sig)(AskLowerTriangularHandler.AlwaysTrue)

AskLowerTriangularHandler.register(Factorization)(partial(_Factorization, Q.lower_triangular))

@AskLowerTriangularHandler.register(MatMul)
def _(expr, assumptions):
    factor, matrices = expr.as_coeff_matrices()
    if all(ask(Q.lower_triangular(m), assumptions) for m in matrices):
        return True

@AskLowerTriangularHandler.register(MatAdd)
def _(expr, assumptions):
    if all(ask(Q.lower_triangular(arg), assumptions) for arg in expr.args):
        return True

@AskLowerTriangularHandler.register(MatPow)
def _(expr, assumptions):
    # only for integer powers
    base, exp = expr.args
    int_exp = ask(Q.integer(exp), assumptions)
    if not int_exp:
        return None
    non_negative = ask(~Q.negative(exp), assumptions)
    if (non_negative or non_negative == False
                        and ask(Q.invertible(base), assumptions)):
        return ask(Q.lower_triangular(base), assumptions)
    return None

@AskLowerTriangularHandler.register(MatrixSymbol)
def _(expr, assumptions):
    if Q.lower_triangular(expr) in conjuncts(assumptions):
        return True

@AskLowerTriangularHandler.register(OneMatrix)
def _(expr, assumptions):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@AskLowerTriangularHandler.register(Transpose)
def _(expr, assumptions):
    return ask(Q.upper_triangular(expr.arg), assumptions)

@AskLowerTriangularHandler.register(Inverse)
def _(expr, assumptions):
    return ask(Q.lower_triangular(expr.arg), assumptions)

@AskLowerTriangularHandler.register(MatrixSlice)
def _(expr, assumptions):
    if not expr.on_diag:
        return None
    else:
        return ask(Q.lower_triangular(expr.parent), assumptions)


### AskDiagonalHandler ###

AskDiagonalHandler = AskHandlerClass(
    'AskDiagonalHandler',
    doc="""
    Handler for key 'diagonal'
    """,
    base=CommonHandler
)

for sig in (Identity, ZeroMatrix, DiagonalMatrix, DiagMatrix):
    AskDiagonalHandler.register(sig)(AskDiagonalHandler.AlwaysTrue)

AskDiagonalHandler.register(Factorization)(partial(_Factorization, Q.diagonal))

def _is_empty_or_1x1(expr):
    return expr.shape == (0, 0) or expr.shape == (1, 1)

@AskDiagonalHandler.register(MatMul)
def _(expr, assumptions):
    if _is_empty_or_1x1(expr):
        return True
    factor, matrices = expr.as_coeff_matrices()
    if all(ask(Q.diagonal(m), assumptions) for m in matrices):
        return True

@AskDiagonalHandler.register(MatPow)
def _(expr, assumptions):
    # only for integer powers
    base, exp = expr.args
    int_exp = ask(Q.integer(exp), assumptions)
    if not int_exp:
        return None
    non_negative = ask(~Q.negative(exp), assumptions)
    if (non_negative or non_negative == False
                        and ask(Q.invertible(base), assumptions)):
        return ask(Q.diagonal(base), assumptions)
    return None

@AskDiagonalHandler.register(MatAdd)
def _(expr, assumptions):
    if all(ask(Q.diagonal(arg), assumptions) for arg in expr.args):
        return True

@AskDiagonalHandler.register(MatrixSymbol)
def _(expr, assumptions):
    if _is_empty_or_1x1(expr):
        return True
    if Q.diagonal(expr) in conjuncts(assumptions):
        return True

@AskDiagonalHandler.register(OneMatrix)
def _(expr, assumptions):
    return expr.shape[0] == 1 and expr.shape[1] == 1

@AskDiagonalHandler.register(Transpose)
def _(expr, assumptions):
    return ask(Q.diagonal(expr.arg), assumptions)

@AskDiagonalHandler.register(Inverse)
def _(expr, assumptions):
    return ask(Q.diagonal(expr.arg), assumptions)

@AskDiagonalHandler.register(MatrixSlice)
def _(expr, assumptions):
    if _is_empty_or_1x1(expr):
        return True
    if not expr.on_diag:
        return None
    else:
        return ask(Q.diagonal(expr.parent), assumptions)


### AskIntegerElementsHandler ###

AskIntegerElementsHandler = AskHandlerClass(
    'AskIntegerElementsHandler',
    doc="""
    Handler for key 'integer_elements'
    """,
    base=CommonHandler
)

for sig in (ZeroMatrix, OneMatrix, Identity):
    AskIntegerElementsHandler.register(sig)(AskIntegerElementsHandler.AlwaysTrue)

AskIntegerElementsHandler.register(MatMul)(partial(MatMul_elements, Q.integer_elements, Q.integer))

AskIntegerElementsHandler.register(MatrixSlice)(partial(MS_elements, Q.integer_elements))

AskIntegerElementsHandler.register(BlockMatrix)(partial(BM_elements, Q.integer_elements))

@AskIntegerElementsHandler.register(MatAdd)
@AskIntegerElementsHandler.register(HadamardProduct)
@AskIntegerElementsHandler.register(Determinant)
@AskIntegerElementsHandler.register(Trace)
@AskIntegerElementsHandler.register(Transpose)
def _(expr, assumptions):
    return test_closed_group(expr, assumptions, Q.integer_elements)

@AskIntegerElementsHandler.register(MatPow)
def _(expr, assumptions):
    # only for integer powers
    base, exp = expr.args
    int_exp = ask(Q.integer(exp), assumptions)
    if not int_exp:
        return None
    if exp.is_negative == False:
        return ask(Q.integer_elements(base), assumptions)
    return None


### AskRealElementsHandler ###

AskRealElementsHandler = AskHandlerClass(
    'AskRealElementsHandler',
    doc="""
    Handler for key 'real_elements'
    """,
    base=CommonHandler
)

AskRealElementsHandler.register(MatMul)(partial(MatMul_elements, Q.real_elements, Q.real))

AskRealElementsHandler.register(MatrixSlice)(partial(MS_elements, Q.real_elements))

AskRealElementsHandler.register(BlockMatrix)(partial(BM_elements, Q.real_elements))

@AskRealElementsHandler.register(MatAdd)
@AskRealElementsHandler.register(HadamardProduct)
@AskRealElementsHandler.register(Determinant)
@AskRealElementsHandler.register(Trace)
@AskRealElementsHandler.register(Transpose)
@AskRealElementsHandler.register(Factorization)
def _(expr, assumptions):
    return test_closed_group(expr, assumptions, Q.real_elements)

@AskRealElementsHandler.register(MatPow)
def _(expr, assumptions):
    # only for integer powers
    base, exp = expr.args
    int_exp = ask(Q.integer(exp), assumptions)
    if not int_exp:
        return None
    non_negative = ask(~Q.negative(exp), assumptions)
    if (non_negative or non_negative == False
                        and ask(Q.invertible(base), assumptions)):
        return ask(Q.real_elements(base), assumptions)
    return None


### AskComplexElementsHandler ###

AskComplexElementsHandler = AskHandlerClass(
    'AskComplexElementsHandler',
    doc="""
    Handler for key 'complex_elements'
    """,
    base=CommonHandler
)

AskComplexElementsHandler.register(MatMul)(partial(MatMul_elements, Q.complex_elements, Q.complex))

AskComplexElementsHandler.register(MatrixSlice)(partial(MS_elements, Q.complex_elements))

AskComplexElementsHandler.register(BlockMatrix)(partial(BM_elements, Q.complex_elements))

AskComplexElementsHandler.register(DFT)(AskComplexElementsHandler.AlwaysTrue)

@AskComplexElementsHandler.register(MatAdd)
@AskComplexElementsHandler.register(HadamardProduct)
@AskComplexElementsHandler.register(Determinant)
@AskComplexElementsHandler.register(Trace)
@AskComplexElementsHandler.register(Transpose)
@AskComplexElementsHandler.register(Inverse)
@AskComplexElementsHandler.register(Factorization)
def _(expr, assumptions):
    return test_closed_group(expr, assumptions, Q.complex_elements)

@AskComplexElementsHandler.register(MatPow)
def _(expr, assumptions):
    # only for integer powers
    base, exp = expr.args
    int_exp = ask(Q.integer(exp), assumptions)
    if not int_exp:
        return None
    non_negative = ask(~Q.negative(exp), assumptions)
    if (non_negative or non_negative == False
                        and ask(Q.invertible(base), assumptions)):
        return ask(Q.complex_elements(base), assumptions)
    return None
