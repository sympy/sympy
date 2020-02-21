"""
The Schur number S(k) is the largest integer n for which the interval [1,n]
can be partitioned into k sum-free sets.(http://mathworld.wolfram.com/SchurNumber.html)
"""
import math
from sympy.core import S
from sympy.core.function import Function
from sympy.core.numbers import Integer


class schur_number_lower_bound(Function):
    """
    Gives the value for 1/2 * (3**n  - 1)
    """

    @classmethod
    def eval(cls, k):
        expr = (3**k - 1)/2

        if k.is_Number:
            if k is S.Infinity:
                return S.Infinity
            if k.is_zero:
                return 0
            if not k.is_Integer or k.is_negative:
                raise ValueError("k should be a positive integer")

        return Integer(expr)


class schur_number_subsets_lower_bound(Function):
    """
    This function returns a lower bound to Schur's Number
    """

    @classmethod
    def eval(cls, n):
        n = int(n)

        if n is S.Infinity:
            return S.Infinity
        if n <= 0:
            raise ValueError("n must be a positive integer.")
        elif n <= 3:
            min_k = 1
        else:
            min_k = math.ceil(math.log(2*n + 1, 3))

        return Integer(min_k)


def schur_partition(n):
    """

    This function returns the partition in the minimum number of sum-free subsets
    according to the lower bound given by the Schur Number.

    Parameters
    ==========

    n: a number
        n is the upper limit of the range [1, n] for which we need to find and
        return the minimum number of free subsets according to the lower bound
        of schur number

    Returns
    =======

    List of lists
        List of the minimum number of sum-free subsets

    Notes
    =====

    It is possible for some n to make the partition into less
    subsets since the only known Schur numbers are:
    S(1) = 1, S(2) = 4 , S(3) = 13, S(4) = 44.
    e.g for n = 44 the lower bound from the function above is 5 subsets but it has been proven
    that can be done with 4 subsets.

    Examples
    ========

    For n = 1, 2, 3 the answer is the set itself

    >>> from sympy.combinatorics.schur_number import schur_partition
    >>> schur_partition(2)
    [[1, 2]]

    For n > 3, the answer is the minimum number of sum-free subsets:

    >>> schur_partition(5)
    [[3, 2], [5], [1, 4]]

    >>> schur_partition(8)
    [[3, 2], [6, 5, 8], [1, 4, 7]]
    """
    n = int(n)
    number_of_subsets = schur_number_subsets_lower_bound(n)
    if n == 1:
        sum_free_subsets = [[1]]
    elif n == 2:
        sum_free_subsets = [[1, 2]]
    elif n == 3:
        sum_free_subsets = [[1, 2, 3]]
    else:
        sum_free_subsets = [[1, 4], [2, 3]]

    while len(sum_free_subsets) < number_of_subsets:
        sum_free_subsets = _generate_next_list(sum_free_subsets, n)
        missed_elements = [3*k + 1 for k in range(len(sum_free_subsets), (n-1)//3 + 1)]
        sum_free_subsets[-1] += missed_elements

    return sum_free_subsets


def _generate_next_list(current_list, n):
    new_list = []

    for item in current_list:
        temp_1 = [number*3 for number in item if number*3 <= n]
        temp_2 = [number*3 - 1 for number in item if number*3 - 1 <= n]
        new_item = temp_1 + temp_2
        new_list.append(new_item)

    last_list = [3*k + 1 for k in range(0, len(current_list)+1) if 3*k + 1 <= n]
    new_list.append(last_list)
    current_list = new_list

    return current_list
