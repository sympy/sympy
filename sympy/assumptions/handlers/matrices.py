"""
This module contains query handlers responsible for calculus queries:
infinitesimal, bounded, etc.
"""
from sympy.logic.boolalg import conjuncts
from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler
from sympy.matrices.expressions import MatMul

class AskSquareHandler(CommonHandler):
    """
    Handler for key 'square'
    """

    @staticmethod
    def MatrixExpr(expr, assumptions):
        return expr.shape[0] == expr.shape[1]


class AskSymmetricHandler(CommonHandler):
    """
    Handler for key 'symmetric'
    """

    @staticmethod
    def MatMul(expr, assumptions):
        factor, mmul = expr.as_coeff_mmul()
        if all(ask(Q.symmetric(arg), assumptions) for arg in mmul.args):
            return True
        if len(mmul.args) >= 2 and mmul.args[0] == mmul.args[-1].T:
            return ask(Q.symmetric(MatMul(*mmul.args[1:-1])), assumptions)

    @staticmethod
    def MatAdd(expr, assumptions):
        return all(ask(Q.symmetric(arg), assumptions) for arg in expr.args)

    @staticmethod
    def MatrixSymbol(expr, assumptions):
        if not expr.is_square:
            return False
        if Q.symmetric(expr) in conjuncts(assumptions):
            return True

    ZeroMatrix, Identity = [staticmethod(CommonHandler.AlwaysTrue)]*2

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.symmetric(expr.arg), assumptions)
    Inverse = Transpose

    @staticmethod
    def MatrixSlice(expr, assumptions):
        if not expr.on_diag:
            return None
        else:
            return ask(Q.symmetric(expr.parent), assumptions)


class AskInvertibleHandler(CommonHandler):
    """
    Handler for key 'invertible'
    """

    @staticmethod
    def MatMul(expr, assumptions):
        factor, mmul = expr.as_coeff_mmul()
        if all(ask(Q.invertible(arg), assumptions) for arg in mmul.args):
            return True
        if any(ask(Q.invertible(arg), assumptions) is False
               for arg in mmul.args):
            return False

    @staticmethod
    def MatAdd(expr, assumptions):
        return None

    @staticmethod
    def MatrixSymbol(expr, assumptions):
        if not expr.is_square:
            return False
        if Q.invertible(expr) in conjuncts(assumptions):
            return True

    Identity, Inverse = [staticmethod(CommonHandler.AlwaysTrue)]*2

    ZeroMatrix = staticmethod(CommonHandler.AlwaysFalse)

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.invertible(expr.arg), assumptions)

    @staticmethod
    def MatrixSlice(expr, assumptions):
        if not expr.on_diag:
            return None
        else:
            return ask(Q.invertible(expr.parent), assumptions)

class AskOrthogonalHandler(CommonHandler):
    """
    Handler for key 'orthogonal'
    """
    @staticmethod
    def MatMul(expr, assumptions):
        factor, mmul = expr.as_coeff_mmul()
        if (all(ask(Q.orthogonal(arg), assumptions) for arg in mmul.args) and
                factor == 1):
            return True
        if any(ask(Q.invertible(arg), assumptions) is False
                for arg in mmul.args):
            return False

    @staticmethod
    def MatAdd(expr, assumptions):
        if (len(expr.args) == 1 and
                ask(Q.orthogonal(expr.args[0]), assumptions)):
            return True

    @staticmethod
    def MatrixSymbol(expr, assumptions):
        if not expr.is_square:
            return False
        if Q.orthogonal(expr) in conjuncts(assumptions):
            return True

    Identity = staticmethod(CommonHandler.AlwaysTrue)

    ZeroMatrix = staticmethod(CommonHandler.AlwaysFalse)

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.orthogonal(expr.arg), assumptions)
    Inverse = Transpose

    @staticmethod
    def MatrixSlice(expr, assumptions):
        if not expr.on_diag:
            return None
        else:
            return ask(Q.orthogonal(expr.parent), assumptions)


class AskFullRankHandler(CommonHandler):
    """
    Handler for key 'fullrank'
    """
    @staticmethod
    def MatMul(expr, assumptions):
        if all(ask(Q.fullrank(arg), assumptions) for arg in expr.args):
            return True

    Identity = staticmethod(CommonHandler.AlwaysTrue)

    ZeroMatrix = staticmethod(CommonHandler.AlwaysFalse)

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.fullrank(expr.arg), assumptions)
    Inverse = Transpose

    @staticmethod
    def MatrixSlice(expr, assumptions):
        if ask(Q.orthogonal(expr.parent), assumptions):
            return True

class AskPositiveDefiniteHandler(CommonHandler):
    """
    Handler for key 'positive_definite'
    """
    @staticmethod
    def MatMul(expr, assumptions):
        factor, mmul = expr.as_coeff_mmul()
        if (all(ask(Q.positive_definite(arg), assumptions)
                for arg in mmul.args) and factor > 0):
            return True
        if (len(mmul.args) >= 2
                and mmul.args[0] == mmul.args[-1].T
                and ask(Q.fullrank(mmul.args[0]), assumptions)):
            return ask(Q.positive_definite(
                MatMul(*mmul.args[1:-1])), assumptions)

    @staticmethod
    def MatAdd(expr, assumptions):
        if all(ask(Q.positive_definite(arg), assumptions)
                for arg in expr.args):
            return True

    @staticmethod
    def MatrixSymbol(expr, assumptions):
        if not expr.is_square:
            return False
        if Q.positive_definite(expr) in conjuncts(assumptions):
            return True

    Identity = staticmethod(CommonHandler.AlwaysTrue)

    ZeroMatrix = staticmethod(CommonHandler.AlwaysFalse)

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.positive_definite(expr.arg), assumptions)
    Inverse = Transpose

    @staticmethod
    def MatrixSlice(expr, assumptions):
        if not expr.on_diag:
            return None
        else:
            return ask(Q.positive_definite(expr.parent), assumptions)

class AskUpperTriangularHandler(CommonHandler):
    """
    Handler for key 'upper_triangular'
    """
    @staticmethod
    def MatMul(expr, assumptions):
        factor, matrices = expr.as_coeff_matrices()
        if all(ask(Q.upper_triangular(m), assumptions) for m in matrices):
            return True

    @staticmethod
    def MatAdd(expr, assumptions):
        if all(ask(Q.upper_triangular(arg), assumptions) for arg in expr.args):
            return True

    @staticmethod
    def MatrixSymbol(expr, assumptions):
        if Q.upper_triangular(expr) in conjuncts(assumptions):
            return True

    Identity, ZeroMatrix = [staticmethod(CommonHandler.AlwaysTrue)]*2

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.lower_triangular(expr.arg), assumptions)

    @staticmethod
    def Inverse(expr, assumptions):
        return ask(Q.upper_triangular(expr.arg), assumptions)

    @staticmethod
    def MatrixSlice(expr, assumptions):
        if not expr.on_diag:
            return None
        else:
            return ask(Q.upper_triangular(expr.parent), assumptions)

class AskLowerTriangularHandler(CommonHandler):
    """
    Handler for key 'lower_triangular'
    """
    @staticmethod
    def MatMul(expr, assumptions):
        factor, matrices = expr.as_coeff_matrices()
        if all(ask(Q.lower_triangular(m), assumptions) for m in matrices):
            return True

    @staticmethod
    def MatAdd(expr, assumptions):
        if all(ask(Q.lower_triangular(arg), assumptions) for arg in expr.args):
            return True

    @staticmethod
    def MatrixSymbol(expr, assumptions):
        if Q.lower_triangular(expr) in conjuncts(assumptions):
            return True

    Identity, ZeroMatrix = [staticmethod(CommonHandler.AlwaysTrue)]*2

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.upper_triangular(expr.arg), assumptions)

    @staticmethod
    def Inverse(expr, assumptions):
        return ask(Q.lower_triangular(expr.arg), assumptions)

    @staticmethod
    def MatrixSlice(expr, assumptions):
        if not expr.on_diag:
            return None
        else:
            return ask(Q.lower_triangular(expr.parent), assumptions)

class AskDiagonalHandler(CommonHandler):
    """
    Handler for key 'diagonal'
    """
    @staticmethod
    def MatMul(expr, assumptions):
        factor, matrices = expr.as_coeff_matrices()
        if all(ask(Q.diagonal(m), assumptions) for m in matrices):
            return True

    @staticmethod
    def MatAdd(expr, assumptions):
        if all(ask(Q.diagonal(arg), assumptions) for arg in expr.args):
            return True

    @staticmethod
    def MatrixSymbol(expr, assumptions):
        if Q.diagonal(expr) in conjuncts(assumptions):
            return True

    Identity, ZeroMatrix = [staticmethod(CommonHandler.AlwaysTrue)]*2

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.diagonal(expr.arg), assumptions)

    @staticmethod
    def Inverse(expr, assumptions):
        return ask(Q.diagonal(expr.arg), assumptions)

    @staticmethod
    def MatrixSlice(expr, assumptions):
        if not expr.on_diag:
            return None
        else:
            return ask(Q.diagonal(expr.parent), assumptions)
