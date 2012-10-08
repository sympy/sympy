from matexpr import MatrixExpr, ShapeError, Identity, ZeroMatrix
from sympy.core import Mul, Add, Basic, sympify

class MatMul(MatrixExpr):
    """A Product of Matrix Expressions

    MatMul inherits from and operates like SymPy Mul

    >>> from sympy import MatMul, MatrixSymbol
    >>> A = MatrixSymbol('A', 5, 4)
    >>> B = MatrixSymbol('B', 4, 3)
    >>> C = MatrixSymbol('C', 3, 6)
    >>> MatMul(A, B, C)
    A*B*C
    """
    is_MatMul = True

    def __new__(cls, *args, **kwargs):
        simplify = kwargs.get('simplify', True)
        check    = kwargs.get('check'   , True)

        args = map(sympify, args)
        obj = Basic.__new__(cls, *args)
        factor, matrices = obj.as_coeff_matrices()
        if check:
            validate(*matrices)
        if simplify:
            return canonicalize(obj)
        return obj

    @property
    def shape(self):
        matrices = [arg for arg in self.args if arg.is_Matrix]
        return (matrices[0].rows, matrices[-1].cols)

    def _entry(self, i, j):
        coeff, matmul = self.as_coeff_mmul()

        if not matmul.is_MatMul: # situation like 2*X, matmul is just X
            return coeff * matmul[i, j]

        head, tail = matmul.args[0], matmul.args[1:]
        assert len(tail) != 0

        X = head
        Y = MatMul(*tail)

        if X.shape[1].is_Number:
            # Numeric shape like (3,5)
            return coeff*Add(*[X[i, k]*Y[k, j] for k in range(X.shape[1])])
        else:
            # Symbolic shape like (n, m)
            from sympy import Dummy, summation
            k = Dummy('k', integer=True)
            return summation(coeff*X[i, k]*Y[k, j], (k, 0, X.cols - 1))

    def as_coeff_matrices(self):
        scalars = [x for x in self.args if not x.is_Matrix]
        matrices = [x for x in self.args if x.is_Matrix]
        coeff = Mul(*scalars)

        return coeff, matrices

    def as_coeff_mmul(self):
        coeff, matrices = self.as_coeff_matrices()
        return coeff, MatMul(*matrices)

    def _eval_transpose(self):
        from transpose import Transpose
        return MatMul(*[Transpose(arg) for arg in self.args[::-1]])

    def _eval_trace(self):
        factor, mmul = self.as_coeff_mmul()
        if factor != 1:
            from trace import Trace
            return factor * Trace(mmul)
        else:
            raise NotImplementedError("Can't simplify any further")

    def _eval_inverse(self):
        from inverse import Inverse
        try:
            return MatMul(*[Inverse(arg) for arg in self.args[::-1]])
        except ShapeError:
            raise NotImplementedError("Can not decompose this Inverse")

    def canonicalize(self):
        return canonicalize(self)

def validate(*matrices):
    """ Checks for valid shapes for args of MatMul """
    for i in range(len(matrices)-1):
        A,B = matrices[i:i+2]
        if A.cols != B.rows:
            raise ShapeError("Matrices %s and %s are not aligned"%(A, B))

# Rules

from sympy.rr import (rmid, unpack, canon, condition, debug, flatten, chain,
        exhaust, do_one)

def newmul(*args):
    if args[0] == 1:
        args = args[1:]
    return Basic.__new__(MatMul, *args)

def any_zeros(mul):
    if any([arg.is_zero or (arg.is_Matrix and arg.is_ZeroMatrix)
                       for arg in mul.args]):
        matrices = [arg for arg in mul.args if arg.is_Matrix]
        return ZeroMatrix(matrices[0].rows, matrices[-1].cols)
    return mul

def xxinv(mul):
    from sympy.matrices.expressions import Inverse
    """ Y * X * X.I -> Y """
    factor, matrices = mul.as_coeff_matrices()
    for i, (X, Y) in enumerate(zip(matrices[:-1], matrices[1:])):
        if X.is_square and Y.is_square and X == Inverse(Y):
            I = Identity(X.rows)
            return newmul(factor, *(matrices[:i] + [I] + matrices[i+2:]))
    return mul

def remove_ids(mul):
    factor, matrices = mul.as_coeff_matrices()
    if not any(m.is_Identity for m in matrices) or len(matrices) == 1:
        return mul

    non_ids = [x for x in matrices if not x.is_Identity]
    return newmul(factor, *non_ids)

def factor_in_front(mul):
    factor, matrices = mul.as_coeff_matrices()
    if factor == 1:
        return mul
    else:
        return newmul(factor, *matrices)

def condition_matmul(rule):
    is_matmul = lambda x: x.is_Matrix and x.is_MatMul
    return condition(is_matmul, rule)

rules = (any_zeros, remove_ids, xxinv, unpack, rmid(lambda x: x == 1),
         factor_in_front, flatten)

canonicalize = exhaust(condition_matmul(do_one(*rules)))

from matadd import MatAdd
