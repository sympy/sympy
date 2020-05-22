from sympy import (
    Basic, Matrix, latex, ShapeError, Mul,
    ImmutableMatrix, MutableMatrix, MutableDenseMatrix
)
from sympy.polys.polytools import degree

__all__ = ['TransferFunction']

_matrixTypes = (
    Matrix, ImmutableMatrix, MutableMatrix, MutableDenseMatrix)


class TransferFunction(Basic):
    def __new__(cls, args):
        # constructor is from a given transfer function.
        if  isinstance(args, (Matrix, ImmutableMatrix, MutableMatrix)):
            G = args
            obj = Basic.__new__(cls, G)
            obj.G = args
            return obj
        else:
            raise TypeError("Provided argument is of unsupported type.")

    def series(self, other):
        if not isinstance(other, TransferFunction):
            raise TypeError("Argument must be of type TransferFunction, not {}.".
                format(type(other)))
        # assert matching shapes.
        if not self.G.shape[0] == other.G.shape[1]:
            raise ShapeError("Dimensions of the input of the argument and the output of"
                "the system must match.")

        return TransferFunction(self.G * other.G)

    def parallel(self, other):
        if not isinstance(other, TransferFunction):
            raise TypeError("Argument must be of type TransferFunction, not {}.".
                format(type(other)))
        # assert matching shapes.
        if not ((self.G.shape[1] == other.G.shape[1]) and
                (self.G.shape[0] == other.G.shape[0])):
            raise ShapeError("Dimensions of the input and the output must match.")

        return TransferFunction(self.G + other.G)

    def neg(self):
        neg_G = self.G
        for i, j in enumerate(neg_G):
            num, denom = j.as_numer_denom()
            neg_G[i] = Mul(-num, 1/denom)

        return TransferFunction(neg_G)

    def solve(self, u):
        # assert the right shape of u
        if not u.shape[1] == 1:
            raise ShapeError("u must be a column vector, not a matrix.")
        if not self.G.shape[1] == u.shape[0]:
            raise ShapeError("u must have a length of {}.".format(self.G.shape[1]))

        return self.G * u

    @property
    def is_proper(self):
        return all(degree(tf.as_numer_denom()[0]) <=
                   degree(tf.as_numer_denom()[1]) for tf in self.G)

    @property
    def is_strictly_proper(self):
        return all(degree(tf.as_numer_denom()[0]) <
                   degree(tf.as_numer_denom()[1]) for tf in self.G)

    def _repr_latex_(self):
        return '$' + latex(self.G) + '$'
