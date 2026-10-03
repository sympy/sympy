from __future__ import annotations
from sympy.matrices.expressions import MatrixExpr
from sympy.assumptions.ask import Q

class Factorization(MatrixExpr):
    arg = property(lambda self: self.args[0])
    shape = property(lambda self: self.arg.shape)  # type: ignore

class LofLU(Factorization):
    @property
    def predicates(self):
        return (Q.lower_triangular,)
class UofLU(Factorization):
    @property
    def predicates(self):
        return (Q.upper_triangular,)

class LofCholesky(LofLU): pass
class UofCholesky(UofLU): pass

class QofQR(Factorization):
    @property
    def predicates(self):
        return (Q.orthogonal,)
class RofQR(Factorization):
    @property
    def predicates(self):
        return (Q.upper_triangular,)

class EigenVectors(Factorization):
    @property
    def predicates(self):
        return (Q.orthogonal,)
class EigenValues(Factorization):
    @property
    def predicates(self):
        return (Q.diagonal,)

class UofSVD(Factorization):
    @property
    def predicates(self):
        return (Q.orthogonal,)
class SofSVD(Factorization):
    @property
    def predicates(self):
        return (Q.diagonal,)
class VofSVD(Factorization):
    @property
    def predicates(self):
        return (Q.orthogonal,)


def lu(expr):
    return LofLU(expr), UofLU(expr)

def qr(expr):
    return QofQR(expr), RofQR(expr)

def eig(expr):
    return EigenValues(expr), EigenVectors(expr)

def svd(expr):
    return UofSVD(expr), SofSVD(expr), VofSVD(expr)


class JordanVectors(Factorization):
    @property
    def predicates(self):
        return (Q.invertible,)


class JordanBlocks(Factorization):
    @property
    def predicates(self):
        return (Q.upper_triangular,)


def jordan(expr):
    """Jordan decomposition of a matrix expression.

    Returns JordanVectors(expr), JordanBlocks(expr).
    """
    return JordanVectors(expr), JordanBlocks(expr)


def JordanForm(expr):
    """Construct symbolic similarity transform representing the Jordan normal form.

    Returns P * J * P**(-1) where P = JordanVectors(expr) and J = JordanBlocks(expr).
    """
    from sympy.matrices.expressions.inverse import Inverse
    from sympy.matrices.expressions.matmul import MatMul
    P, J = jordan(expr)
    return MatMul(P, J, Inverse(P))
