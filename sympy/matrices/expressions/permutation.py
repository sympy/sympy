from sympy import Basic, S
from matexpr import MatrixExpr
from sympy.core.compatibility import iterable

class PermutationMatrix(MatrixExpr):
    """ Permutation Matrix

    >>> from sympy.matrices import Matrix, PermutationMatrix
    >>> P = PermutationMatrix((1, 2, 0))
    >>> Matrix(P)
    [0, 1, 0]
    [0, 0, 1]
    [1, 0, 0]

    >>> Matrix(3, 3, range(9))
    [0, 1, 2]
    [3, 4, 5]
    [6, 7, 8]

    >>> P*Matrix(3, 3, range(9))
    [3, 4, 5]
    [6, 7, 8]
    [0, 1, 2]
    """

    def __new__(cls, *args):
        from sympy.matrices import ImmutableMatrix, MutableMatrix

        if len(args) == 1 and iterable(args[0]):
            arg = args[0]
        else:
            arg = args

        if not isinstance(arg, MatrixExpr):
            arg = ImmutableMatrix(list(arg)).T

        if arg.shape[0] != 1:
            raise ValueError("Input must be row vector")

        return Basic.__new__(cls, arg)

    arg = property(lambda self: self.args[0])
    shape = property(lambda self: (S(self.arg.cols),)*2)

    @property
    def is_Identity(self):
        return (S(self.arg.cols).is_Number and
                all(self.arg[0, i] == i for i in range(self.arg.cols)))

    def _entry(self, i, j):
        return 1 if self.arg[i] == j else 0
