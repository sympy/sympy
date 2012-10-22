"""
This module contains query handlers responsible for calculus queries:
infinitesimal, bounded, etc.
"""
from sympy.logic.boolalg import conjuncts
from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler
from sympy.matrices.expressions import MatMul

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
            return ask(Q.symmetric(MatMul(*mmul.args[1:-1], evaluate=False)),
                                   assumptions)

    @staticmethod
    def MatAdd(expr, assumptions):
        return all(ask(Q.symmetric(arg), assumptions) for arg in expr.args)

    @staticmethod
    def MatrixSymbol(expr, assumptions):
        if not expr.is_square:
            return False
        if Q.symmetric(expr) in conjuncts(assumptions):
            return True

    @staticmethod
    def Identity(expr, assumptions):
        return True
    ZeroMatrix = Identity

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.symmetric(expr.arg), assumptions)
    Inverse = Transpose

class AskInvertibleHandler(CommonHandler):
    """
    Handler for key 'invertible'
    """

    @staticmethod
    def MatMul(expr, assumptions):
        factor, mmul = expr.as_coeff_mmul()
        if all(ask(Q.invertible(arg), assumptions) for arg in mmul.args):
            return True
        if any(ask(Q.invertible(arg), assumptions)==False for arg in mmul.args):
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

    @staticmethod
    def Identity(expr, assumptions):
        return True

    @staticmethod
    def ZeroMatrix(expr, assumptions):
        return False

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.invertible(expr.arg), assumptions)

    @staticmethod
    def Inverse(expr, assumptions):
        return True

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
        if any(ask(Q.invertible(arg), assumptions) == False
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

    @staticmethod
    def Identity(expr, assumptions):
        return True

    @staticmethod
    def ZeroMatrix(expr, assumptions):
        return False

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.orthogonal(expr.arg), assumptions)
    Inverse = Transpose


class AskPositiveDefiniteHandler(CommonHandler):
    """
    Handler for key 'positive_definite'
    """
    @staticmethod
    def MatMul(expr, assumptions):
        factor, mmul = expr.as_coeff_mmul()
        if (all(ask(Q.positive_definite(arg), assumptions) for arg in mmul.args)
                and factor > 0):
            return True
        if len(mmul.args) >= 2 and mmul.args[0] == mmul.args[-1].T:
            return ask(Q.positive_definite(
                                MatMul(*mmul.args[1:-1], evaluate=False)),
                                assumptions)

    @staticmethod
    def MatAdd(expr, assumptions):
        if (all(ask(Q.positive_definite(arg), assumptions) for arg in mmul.args)
                and factor > 0):
            return True

    @staticmethod
    def MatrixSymbol(expr, assumptions):
        if not expr.is_square:
            return False
        if Q.positive_definite(expr) in conjuncts(assumptions):
            return True

    @staticmethod
    def Identity(expr, assumptions):
        return True

    @staticmethod
    def ZeroMatrix(expr, assumptions):
        return False

    @staticmethod
    def Transpose(expr, assumptions):
        return ask(Q.positive_definite(expr.arg), assumptions)
    Inverse = Transpose
