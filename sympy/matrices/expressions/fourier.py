from sympy.core.sympify import _sympify
from sympy.matrices.expressions import MatrixExpr
from sympy.matrices.expressions.determinant import Determinant
from sympy import S, I, sqrt, exp


class DFT(MatrixExpr):
    """ Discrete Fourier Transform """
    def __new__(cls, n):
        n = _sympify(n)
        cls._check_dim(n)

        obj = super().__new__(cls, n)
        return obj

    n = property(lambda self: self.args[0])  # type: ignore
    shape = property(lambda self: (self.n, self.n))  # type: ignore

    def _entry(self, i, j, **kwargs):
        w = exp(-2*S.Pi*I/self.n)
        return w**(i*j) / sqrt(self.n)

    def _eval_inverse(self):
        return IDFT(self.n)

    def _eval_determinant(self):
        eigs = self._eval_eigenvals()
        if eigs:
            tmp = 1
            for e, cnt in eigs.items():
                tmp *= e**cnt
            return tmp
        return Determinant(self)

    def _eval_eigenvals(self):
        # Based on Dickinson and Steiglitz, "Eigenvectors and Functions of the
        # Discrete Fourier Transform", 1982
        if self.n.is_Integer:
            m = self.n // 4
            r = self.n % 4
            if r == 0:
                ret = {1: m+1, -1: m, I: m-1, -I: m}
            elif r == 1:
                ret = {1: m+1, -1: m, I: m, -I: m}
            elif r == 2:
                ret = {1: m+1, -1: m+1, I: m, -I: m}
            elif r == 3:
                ret = {1: m+1, -1: m+1, I: m, -I: m+1}
            return {eig: ret[eig] for eig in ret if ret[eig] > 0}


class IDFT(DFT):
    """ Inverse Discrete Fourier Transform """
    def _entry(self, i, j, **kwargs):
        w = exp(-2*S.Pi*I/self.n)
        return w**(-i*j) / sqrt(self.n)

    def _eval_inverse(self):
        return DFT(self.n)

    def _eval_eigenvals(self):
        # Based on Dickinson and Steiglitz, "Eigenvectors and Functions of the
        # Discrete Fourier Transform", 1982
        if self.n.is_Integer:
            m = self.n // 4
            r = self.n % 4
            if r == 0:
                ret = {S.One: m+1, -1: m, I: m, -I: m-1}
            elif r == 1:
                ret = {S.One: m+1, -1: m, I: m, -I: m}
            elif r == 2:
                ret = {S.One: m+1, -1: m+1, I: m, -I: m}
            elif r == 3:
                ret = {S.One: m+1, -1: m+1, I: m+1, -I: m}
            return {eig: ret[eig] for eig in ret if ret[eig] > 0}
