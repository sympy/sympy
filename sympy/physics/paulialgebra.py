"""
This module implements Pauli algebra by subclassing Symbol. Only algebraic
properties of Pauli matrices are used (we don't use the Matrix class).

See the documentation to the class Pauli for examples.

References
~~~~~~~~~~
.. [1] http://en.wikipedia.org/wiki/Pauli_matrices
"""

from __future__ import print_function, division

from sympy import Symbol, I

__all__ = ['evaluate_pauli_product']


def delta(i, j):
    """
    Returns 1 if i == j, else 0.

    This is used in the multiplication of Pauli matrices.

    Examples
    ========

    >>> from sympy.physics.paulialgebra import delta
    >>> delta(1, 1)
    1
    >>> delta(2, 3)
    0
    """
    if i == j:
        return 1
    else:
        return 0


def epsilon(i, j, k):
    """
    Return 1 if i,j,k is equal to (1,2,3), (2,3,1), or (3,1,2);
    -1 if i,j,k is equal to (1,3,2), (3,2,1), or (2,1,3);
    else return 0.

    This is used in the multiplication of Pauli matrices.

    Examples
    ========

    >>> from sympy.physics.paulialgebra import epsilon
    >>> epsilon(1, 2, 3)
    1
    >>> epsilon(1, 3, 2)
    -1
    """
    if (i, j, k) in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]:
        return 1
    elif (i, j, k) in [(1, 3, 2), (3, 2, 1), (2, 1, 3)]:
        return -1
    else:
        return 0


class Pauli(Symbol):
    """The class representing algebraic properties of Pauli matrices

    If the left multiplication of symbol or number with Pauli matrix is needed,
    please use parentheses  to separate Pauli and symbolic multiplication
    (for example: 2*I*(Pauli(3)*Pauli(2)))

    Another variant is to use evaluate_pauli_product function to evaluate
    the product of Pauli matrices and other symbols (with commutative
    multiply rules)

    See Also
    =======
    evaluate_pauli_product

    Examples
    ========

    >>> from sympy.physics.paulialgebra import Pauli
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

    >>> from sympy import I
    >>> I*(Pauli(2)*Pauli(3))
    -sigma1

    >>> from sympy.physics.paulialgebra import evaluate_pauli_product
    >>> f = I*Pauli(2)*Pauli(3)
    >>> f
    I*sigma2*sigma3
    >>> evaluate_pauli_product(f)
    -sigma1

    """

    __slots__ = ["i"]

    def __new__(cls, i):
        if not i in [1, 2, 3]:
            raise IndexError("Invalid Pauli index")
        obj = Symbol.__new__(cls, "sigma%d" % i, commutative=False)
        obj.i = i
        return obj

    def __getnewargs__(self):
        return (self.i,)

    # FIXME don't work for -I*Pauli(2)*Pauli(3)
    def __mul__(self, other):
        if isinstance(other, Pauli):
            j = self.i
            k = other.i
            return delta(j, k) \
                + I*epsilon(j, k, 1)*Pauli(1) \
                + I*epsilon(j, k, 2)*Pauli(2) \
                + I*epsilon(j, k, 3)*Pauli(3)
        return super(Pauli, self).__mul__(other)

    def _eval_power(b, e):
        if e.is_Integer and e.is_positive:
            return super(Pauli, b).__pow__(int(e) % 2)


def evaluate_pauli_product(arg):
    '''Help function to evaluate Pauli matrices product
    with symbolic objects

    Parameters
    ==========

    arg: symbolic expression that contains Paulimatrices

    Examples
    ========

    >>> from sympy.physics.paulialgebra import Pauli, evaluate_pauli_product
    >>> from sympy import I
    >>> evaluate_pauli_product(I*Pauli(1)*Pauli(2))
    -sigma3

    >>> from sympy.abc import x,y
    >>> evaluate_pauli_product(x**2*Pauli(2)*Pauli(1))
    -I*x**2*sigma3
    '''
    tmp = arg.as_coeff_mul()
    sigma_product = 1
    com_product = 1
    for el in tmp[1]:
        if isinstance(el, Pauli):
            sigma_product *= el
        else:
            com_product *= el
    return (tmp[0]*sigma_product*com_product)
