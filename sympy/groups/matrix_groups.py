from __future__ import print_function
from sympy.core import Basic
from math import sin, cos, pi, floor
import random

def MatrixMultiply(a, b):
    """Returns the product of two matrices"""

    _arr = [[] for i in range(len(a))]
    for i in range(len(a)):
        for j in range(len(b[0])):
            c = 0
            for k in range(len(b)):
                c = c + (a[i][k] * b[k][j])
            _arr[i].append(c)
    return _arr

def isPrime(p):
    if p is 2:
        return True
    elif p%2 is 0:
        return False
    for i in range(3, p/2 + 2, 2):
        if (p % i) is 0:
            return False
    return True

def spl_linear_group(dim = 2, min_val = 1, max_val = 10):
    """
    Returns a Special Linear Group with integer entries and
    of specified dimension.
    """

    if dim < 1:
        raise ValueError("Size must be greater than or equal to 1")
    elif dim is 1:
        return [1]

    _spl_linear_group = SplLinearGroup(dim, min_val, max_val)
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
    if not (isPrime(p)):
        raise ValueError("p must be a Prime Number")

    _finite_matrix_group = FiniteMatrixGroup(dim, p)
    return _finite_matrix_group


class SplLinearGroup(Basic):
    def __new__(cls, dim, _min, _max):
        """
        Generates two dim x dim matrices, one upper triangular
        and other lower triangular with diagonal entries as 1
        and other entries as random integers. The resultant of
        their multiplication is returned.
        """
        random.seed()
        a = [[] for i in range(dim)]
        b = [[] for i in range(dim)]

        for i in range(dim):
            for j in range(dim):
                if i is j:
                    a[i].append(1)
                    b[i].append(1)
                elif(i < j):
                    a[i].append(random.randint(_min, _max))
                    b[i].append(0)
                else:
                    a[i].append(0)
                    b[i].append(random.randint(_min, _max))

        _arr = MatrixMultiply(a, b)
        return _arr

class SplOrthogonalGroup(Basic):
    def __new__(cls, dim, x = None):
        _arr = [[] for i in range(dim)]
        if x is None:
            random.seed()
            x = random.uniform(0, 2*pi)

        _arr[0].append(cos(x))
        _arr[0].append(sin(x))
        _arr[1].append(-sin(x))
        _arr[1].append(cos(x))

        return _arr

class FiniteMatrixGroup(Basic):
    def __new__(cls, dim, p):
        random.seed()
        _arr = [[] for i in range(dim)]

        for i in range(dim):
            for j in range(dim):
                x = random.randint(1, 100)
                _arr[i].append(x % p)
        return _arr
