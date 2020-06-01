from sympy import (
    Basic, Matrix, ShapeError, Mul,
    ImmutableMatrix, MutableMatrix, MutableDenseMatrix,
)
from sympy.polys.polytools import degree

__all__ = ['TransferFunctionMatrix', 'TransferFunction']

_matrixTypes = (
    Matrix, ImmutableMatrix, MutableMatrix, MutableDenseMatrix)


class TransferFunctionMatrix(Basic):
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
        if not isinstance(other, TransferFunctionMatrix):
            raise TypeError("Argument must be of type TransferFunctionMatrix, not {}.".
                format(type(other)))
        # assert matching shapes.
        if not self.G.shape[0] == other.G.shape[1]:
            raise ShapeError("Dimensions of the input of the argument and the output of"
                "the system must match.")

        return TransferFunctionMatrix(self.G * other.G)

    def parallel(self, other):
        if not isinstance(other, TransferFunctionMatrix):
            raise TypeError("Argument must be of type TransferFunctionMatrix, not {}.".
                format(type(other)))
        # assert matching shapes.
        if not ((self.G.shape[1] == other.G.shape[1]) and
                (self.G.shape[0] == other.G.shape[0])):
            raise ShapeError("Dimensions of the input and the output must match.")

        return TransferFunctionMatrix(self.G + other.G)

    def neg(self):
        neg_G = self.G
        for i, j in enumerate(neg_G):
            num, denom = j.as_numer_denom()
            neg_G[i] = Mul(-num, 1/denom)

        return TransferFunctionMatrix(neg_G)

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


class TransferFunction(Mul):
    # this will subclass from `Basic`.
    def __new__(cls, num, den):
        obj = Mul.__new__(cls, num, 1/den)
        obj.num = num
        obj.den = den
        return obj

    def cancel_poles_and_zeros(self):
        pass

    def __add__(self, other):
        return self.add(other)

    def __radd__(self, other):
        return self.add(other)

    def __mul__(self, other):
        return self.mul(other)

    def __div__(self, other):
        return self.div(other)

    def __pow__(self, other):
        return self.pow(other)

    def __neg__(self):
        return TransferFunction(-self.num, self.den)

    def add(self, other):
        pass

    def mul(self, other):
        pass

    def div(self, other):
        pass

    def pow(self, other):
        pass

    @property
    def is_proper(self):
        return degree(self.num) <= degree(self.den)

    @property
    def is_strictly_proper(self):
        return degree(self.num) < degree(self.den)
