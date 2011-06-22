import operator
from functools import reduce
from sympy.core.function import Function
from sympy.core import sympify, S, Integer

###############################################################################
###################### Kronecker Delta, Levi-Civita etc. ######################
###############################################################################

class Dij(Function):
    """
    Represents the Kronecker Delta Function

    if i == j, Dij(i, j) = 1
    otherwise Dij(i, j) = 0
    where i, j are usually integers
    """
    nargs = (1, 2)

    @classmethod
    def eval(cls, i, j=0):
        i, j = map(sympify, (i, j))
        if i == j:
            return S.One
        elif i.is_number and j.is_number:
            return S.Zero

class Eijk(Function):
    """
    Represents the Levi-Civita symbol (antisymmetric symbol)
    """
    nargs = 3

    @classmethod
    def eval(cls, i, j, k):
        i, j, k = map(sympify, (i, j, k))
        if (i,j,k) in [(1,2,3), (2,3,1), (3,1,2)]:
            return S.One
        elif (i,j,k) in [(1,3,2), (3,2,1), (2,1,3)]:
            return S.NegativeOne
        elif i==j or j==k or k==i:
            return S.Zero

def prod(a):
    return reduce(operator.mul, a, 1)

def eval_levicivita(*args):
    """Evaluate Levi-Civita symbol."""
    from sympy import factorial
    n = len(args)
    return prod(
        prod(args[j] - args[i] for j in xrange(i + 1, n))
        / factorial(i) for i in xrange(n))
    # converting factorial(i) to int is slightly faster

class LeviCivita(Function):
    """Represent the Levi-Civita symbol.

    For even permutations of indices it returns 1, for odd permutations -1, and
    for everything else (a repeated index) it returns 0.

    Thus it represents an alternating pseudotensor.

    >>> from sympy import LeviCivita, symbols
    >>> LeviCivita(1,2,3)
    1
    >>> LeviCivita(1,3,2)
    -1
    >>> LeviCivita(1,2,2)
    0
    >>> i,j,k = symbols('i j k')
    >>> LeviCivita(i,j,k)
    LeviCivita(i, j, k)
    >>> LeviCivita(i,j,i)
    0
    """
    @classmethod
    def eval(cls, *args):
        if all(isinstance(a, (int, Integer)) for a in args):
            return eval_levicivita(*args)
        if len(set(args)) < len(args):
            return S.Zero

    def doit(self):
        return eval_levicivita(*self.args)

