"""
The Schur number S(k) is the largest integer n for which the interval [1,n]
can be partitioned into k sum-free sets.(http://mathworld.wolfram.com/SchurNumber.html)
"""
import math
from sympy import Symbol


def schur_number_lower_bound(n):
    """
    This function returns a lower bound to Schur's Number
    """
    n = int(n)
    if n <= 0:
        raise ValueError("n must be a positive integer.")
    elif n <= 3:
        min_k = 1
    else:
        min_k = math.ceil(math.log(2*n + 1,3))

    return min_k

def schur_partition(n):
    """
    This function makes the partition in the minimum number of sum-free subsets
    according to the lower bound given by the Schur Number.
    The partition is returned in a list of lists.
    *** Note that it is possible for some n to make the partition into less
    subsets since the only known Schur numbers are:
    S(1) = 1, S(2) = 4 , S(3) = 13, S(4) = 44.
    e.g for n = 44 the lower bound from the function above is 5 subsets but it has been proven
    that can be done with 4 subsets.
    """
    n = int(n)
    number_of_subsets = schur_number_lower_bound(n)
    if  n == 1:
        sum_free_subsets = [[1]]
    elif n == 2:
        sum_free_subsets = [[1,2]]
    elif n == 3:
        sum_free_subsets = [[1,2,3]]
    else:
        sum_free_subsets = [[1,4],[2,3]]

    while len(sum_free_subsets) < number_of_subsets  :
            sum_free_subsets = _generate_next_list(sum_free_subsets,n)
            missed_elements = [3*k + 1 for k in range(len(sum_free_subsets),(n-1)//3 + 1)]
            sum_free_subsets[-1] += missed_elements

    return sum_free_subsets

def _generate_next_list(current_list,n):
    new_list = []
    for item in current_list:
        new_item = []
        temp_1 = [number*3 for number in item if number*3 <= n]
        temp_2 = [number*3 - 1 for number in item if number*3 - 1 <= n]
        new_item = temp_1 + temp_2
        new_list.append(new_item)
    last_list = [3*k + 1 for k in range(0, len(current_list)+1) if 3*k +1 <= n]
    new_list.append(last_list)
    current_list = new_list

    return current_list

def Schur_number_lower_bound():
    """
    This function returns the inequality to be solved symbolically
    in order to find the Schur Number
    """
    n = Symbol("n")
    k = Symbol("k")
    return k >= 1/2*(3**n-1)
