from matexpr import MatrixExpr, ShapeError, Identity, ZeroMatrix
from sympy.core import Mul, Add, Basic
from sympy import sympify
from sympy.rules import (rm_id, unpack, condition, debug, flatten, exhaust,
        do_one)

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
        evaluate = kwargs.get('evaluate', True)
        check    = kwargs.get('check', True)

        args = map(sympify, args)
        obj = Basic.__new__(cls, *args)
        factor, matrices = obj.as_coeff_matrices()
        if check:
            validate(*matrices)
        if evaluate:
            return canonicalize(obj)
        return obj

    @property
    def shape(self):
        matrices = [arg for arg in self.args if arg.is_Matrix]
        return (matrices[0].rows, matrices[-1].cols)

    def _entry(self, i, j):
        coeff, matmul = self.as_coeff_mmul()
        if not matmul.is_MatMul: # situation like 2*X, matmul is just X
            return coeff * matmul[i,j]

        head, tail = matmul.args[0], matmul.args[1:]
        assert len(tail) != 0

        X = head
        Y = MatMul(*tail)

        if X.shape[1].is_Number:
            # Numeric shape like (3,5)
            return coeff*Add(*[X[i,k]*Y[k,j] for k in range(X.shape[1])])
        else:
            # Symbolic shape like (n, m)
            from sympy import Dummy, summation
            k = Dummy('k', integer=True)
            return summation(coeff*X[i,k]*Y[k,j], (k, 0, X.cols-1))

    def as_coeff_matrices(self):
        scalars = [x for x in self.args if not x.is_Matrix]
        matrices = [x for x in self.args if x.is_Matrix]
        coeff = Mul(*scalars)

        return coeff, matrices

    def as_coeff_mmul(self):
        coeff, matrices = self.as_coeff_matrices()
        return coeff, Basic.__new__(MatMul, *matrices)

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
    """ Y * X * X.I -> Y """
    from sympy.matrices.expressions import Inverse
    factor, matrices = mul.as_coeff_matrices()
    for i, (X, Y) in enumerate(zip(matrices[:-1], matrices[1:])):
        if X.is_square and Y.is_square and X == Inverse(Y):
            I = Identity(X.rows)
            return newmul(factor, *(matrices[:i] + [I] + matrices[i+2:]))
    return mul

def remove_ids(mul):
    """ Remove Identities from a MatMul

    This is a modified version of sympy.rules.rm_id.
    This is necesssary because MatMul may contain both MatrixExprs and Exprs
    as args.

    See Also
    --------
        sympy.rules.rm_id
    """
    # Separate Exprs from MatrixExprs in args
    factor, mmul = mul.as_coeff_mmul()
    # Apply standard rm_id for MatMuls
    result = rm_id(lambda x: x.is_Identity == True)(mmul)
    if result != mmul:
        return newmul(factor, *result.args) # Recombine and return
    else:
        return mul

def factor_in_front(mul):
    factor, matrices = mul.as_coeff_matrices()
    if factor != 1:
        return newmul(factor, *matrices)
    return mul

rules = (any_zeros, remove_ids, xxinv, unpack, rm_id(lambda x: x == 1),
         factor_in_front, flatten)

canonicalize = exhaust(condition(lambda x: isinstance(x, MatMul),
                                 do_one(*rules)))
