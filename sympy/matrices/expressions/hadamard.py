from __future__ import print_function, division

from sympy.core import Mul, Basic, sympify
from sympy.strategies import unpack, flatten, sort, condition, exhaust, do_one

from sympy.matrices.expressions.matexpr import MatrixExpr, ShapeError
from sympy.matrices import matrix_multiply_elementwise

def hadamard_product(*matrices):
    """
    Return the elementwise (aka Hadamard) product of matrices.

    Examples
    --------
    >>> from sympy.matrices import hadamard_product, MatrixSymbol
    >>> A = MatrixSymbol('A', 2, 3)
    >>> B = MatrixSymbol('B', 2, 3)
    >>> hadamard_product(A)
    A
    >>> hadamard_product(A, B)
    A.*B
    >>> hadamard_product(A, B)[0, 1]
    A[0, 1]*B[0, 1]
    """
    if not matrices:
        raise TypeError("Empty Hadamard product is undefined")
    validate(*matrices)
    return HadamardProduct(*matrices).doit()

class HadamardProduct(MatrixExpr):
    """
    Elementwise product of matrix expressions

    This is a symbolic object that simply stores its argument without
    evaluating it. To actually compute the product, use the function
    ``hadamard_product()``.

    >>> from sympy.matrices import hadamard_product, HadamardProduct, MatrixSymbol
    >>> A = MatrixSymbol('A', 5, 5)
    >>> B = MatrixSymbol('B', 5, 5)
    >>> isinstance(hadamard_product(A, B), HadamardProduct)
    True
    """
    is_HadamardProduct = True

    def __new__(cls, *args, **kwargs):
        args = list(map(sympify, args))
        check = kwargs.get('check'   , True)
        if check:
            validate(*args)
        return super(HadamardProduct, cls).__new__(cls, *args)

    @property
    def shape(self):
        return self.args[0].shape

    def _entry(self, i, j):
        return Mul(*[arg._entry(i, j) for arg in self.args])

    def _eval_transpose(self):
        from sympy.matrices.expressions.transpose import transpose
        return HadamardProduct(*list(map(transpose, self.args)))

    def doit(self, **ignored):
        from sympy.matrices import ImmutableMatrix
        matrices = sorted(self.args, key=lambda args: type(args))
        ans = matrices[0]
        index = 0
        if isinstance(ans, ImmutableMatrix):
            for index in range(1, len(matrices)):
                if not isinstance(matrices[index], ImmutableMatrix):
                    break
                else:
                    ans = matrix_multiply_elementwise(ans, matrices[index])
        if index == 0 and isinstance(matrices[0], ImmutableMatrix):
            return HadamardProduct(ans)
        elif index == 0:
            return canonicalize(self)
        elif index == len(matrices)-1 and isinstance(matrices[-1], ImmutableMatrix):
            return HadamardProduct(ans)
        else:
            matrices = tuple([ans] + matrices[index:])
            return HadamardProduct(*matrices)

def validate(*args):
    if not all(arg.is_Matrix for arg in args):
        raise TypeError("Mix of Matrix and Scalar symbols")
    A = args[0]
    for B in args[1:]:
        if A.shape != B.shape:
            raise ShapeError("Matrices %s and %s are not aligned" % (A, B))

rules = (unpack,
         flatten)

canonicalize = exhaust(condition(lambda x: isinstance(x, HadamardProduct),
                                 do_one(*rules)))
