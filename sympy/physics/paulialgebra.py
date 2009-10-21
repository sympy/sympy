from sympy import Symbol,I,Integer

"""
This module implements Pauli algebra by subclassing Symbol. Only aglebraic
properties of Pauli matrices are used (we don't use the Matrix class).

See the documentation to the class Pauli for examples.

See also:
    http://en.wikipedia.org/wiki/Pauli_matrices
"""

def delta(i,j):
    if i==j:
        return 1
    else:
        return 0

def epsilon(i,j,k):
    if (i,j,k) in [(1,2,3), (2,3,1), (3,1,2)]:
        return 1
    elif (i,j,k) in [(1,3,2), (3,2,1), (2,1,3)]:
        return -1
    else:
        return 0

class Pauli(Symbol):
    """
    >>> from sympy import *
    >>> Pauli(1)
    sigma1
    >>> Pauli(1)*Pauli(2)
    I*sigma3
    >>> Pauli(1)*Pauli(1)
    1
    >>> Pauli(3)**4
    1
    >>> Pauli(1)*Pauli(2)*Pauli(3)
    I

    """

    __slots__ = ["i"]

    def __new__(cls, i):
        if not i in [1,2,3]:
            raise IndexError("Invalid Pauli index")
        obj = Symbol.__new__(cls, "sigma%d"%i, commutative=False)
        obj.i=i
        return obj

    def __getnewargs__(self):
        return (self.i,)

    # FIXME don't work for -I*Pauli(2)*Pauli(3)
    def __mul__(self, other):
        if isinstance(other, Pauli):
            j=self.i
            k=other.i
            return delta(j,k) \
                +I*epsilon(j,k,1)*Pauli(1) \
                +I*epsilon(j,k,2)*Pauli(2) \
                +I*epsilon(j,k,3)*Pauli(3)
        return super(Pauli, self).__mul__(other)

    def _eval_power(b, e):
        if e.is_Integer and e.is_positive:
            return super(Pauli, b).__pow__(int(e) % 2)
