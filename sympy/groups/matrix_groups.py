from __future__ import print_function
from sympy.core import Basic
from sympy import Matrix
from sympy.ntheory.primetest import isprime
from math import sin, cos, pi
import random

def spl_linear_group(dim = 2):
    """
    Returns a Special Linear Group with integer entries and
    of specified dimension.
    """

    if dim < 1:
        raise ValueError("Size must be greater than or equal to 1")
    elif dim is 1:
        return [1]

    _spl_linear_group = SplLinearGroup(dim)
    return _spl_linear_group

def spl_orthogonal_group(dim = 2, x = None):
    """
    Returns a Special Orthogonal Group of 2x2 matrices with entries as
    [  cos(x) sin(x)]
    [ -sin(x) cos(x)]
    where x is a real number.
    """

    if dim < 1:
        raise ValueError("Size must be greater than or equal to 1")
    elif dim is 1:
        return [1]
    _spl_orthogonal_group = SplOrthogonalGroup(dim, x)
    return _spl_orthogonal_group

def finite_matrix_group(dim = 2, p = 2):
    """
    Returns a dim x dim matrix with entries as integer modulo p, where
    p is a Prime number.
    """
    if dim < 1:
        raise ValueError("Size must be greater than or equal to 1")
    if not (isprime(p)):
        raise ValueError("p must be a Prime Number")

    _finite_matrix_group = FiniteMatrixGroup(dim, p)
    return _finite_matrix_group

class SplLinearGroup(Basic):
    def __new__(cls, dim):
        """
        Generates two dim x dim matrices, one upper triangular
        and other lower triangular with diagonal entries as 1
        and other entries as random integers. The resultant of
        their multiplication is returned.
        """

        def f_1(i, j):
            if i == j:
                return 1
            elif i < j:
                return random.randint(1, 10)
            else:
                return 0
        def f_2(i, j):
            if i == j:
                return 1
            elif i > j:
                return random.randint(1, 10)
            else:
                return 0

        M1 = Matrix(dim, dim, f_1)
        M2 = Matrix(dim, dim, f_2)
        return M1*M2

class SplOrthogonalGroup(Basic):
    def __new__(cls, dim, x = None):
        if x is None:
            random.seed()
            x = random.uniform(0, 2*pi)

        M = Matrix(([cos(x), sin(x)], [-sin(x), cos(x)]))
        return M

class FiniteMatrixGroup(Basic):
    def __new__(cls, dim, p):
        random.seed()
        _arr = [[] for i in range(dim)]
        M = Matrix(dim, dim, lambda i, j: random.randint(1, 100) % p)
        return M
