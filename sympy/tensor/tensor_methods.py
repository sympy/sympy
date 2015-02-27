# -*- coding: utf-8 -*-

from itertools import *
from sympy.tensor.arraypy import Arraypy, Tensor
from random import randint
from sympy.functions.combinatorial.factorials import factorial

def symmetric(in_arr):
    """
    Creates the symmetric form of input tensor.
    Input: Arraypy or Tensor with equal axes (array shapes).
    Output: symmetric array. Output type - Arraypy or Tensor, depends of input
            
    Examples
    ========
    
    >>> a = list2arraypy(range(9), (3,3))
    >>> b = symmetric(a)
    >>> print (b)
    0.0 2.0 4.0
    2.0 4.0 6.0
    4.0 6.0 8.0
    """
    if not isinstance(in_arr, Arraypy):
        raise TypeError('Input must be Arraypy or Tensor type')

    flag = 0
    for j in range(0, in_arr.rank):
        if (in_arr.shape[0] != in_arr.shape[j]):
            raise ValueError('Different size of arrays axes')

    # forming list of tuples for Arraypy constructor of type a = Arraypy( [(a,
    # b), (c, d), ... , (y, z)] )
    arg = [(in_arr.start_index[i], in_arr.end_index[i])
           for i in range(in_arr.rank)]

    # if in_arr Tensor, then res_arr will be Tensor, else it will be Arraypy
    if isinstance(in_arr, Tensor):
        res_arr = Tensor(Arraypy(arg), in_arr.ind_char)
    else:
        res_arr = Arraypy(arg)

    index = [in_arr.start_index[i] for i in range(in_arr.rank)]

    for i in range(len(in_arr)):
        perm = list(permutations(index))
        for temp_index in perm:
            res_arr[tuple(index)] += in_arr[tuple(temp_index)]
        if isinstance(res_arr[tuple(index)], int):
            res_arr[tuple(index)] = float(res_arr[tuple(index)])
        res_arr[tuple(index)] /= factorial(in_arr.rank)

        index = in_arr.next_index(index)

    return res_arr


def asymmetric(in_arr):
    """
    Creates the asymmetric form of input tensor.
    Input: Arraypy or Tensor with equal axes (array shapes).
    Output: asymmetric array. Output type - Arraypy or Tensor, depends of input
            
    Examples
    ========
    
    >>> a = list2arraypy(range(9), (3,3))
    >>> b = asymmetric(a)
    >>> print (b)
    0.0 -1.0 -2.0
    1.0 0.0 -1.0
    2.0 1.0 0.0
    """
    if not isinstance(in_arr, Arraypy):
        raise TypeError('Input must be Arraypy or Tensor type')

    flag = 0
    for j in range(in_arr.rank):
        if (in_arr.shape[0] != in_arr.shape[j]):
            raise ValueError('Different size of arrays axes')

    # forming list of tuples for Arraypy constructor of type a = Arraypy( [(a,
    # b), (c, d), ... , (y, z)] )
    arg = [(in_arr.start_index[i], in_arr.end_index[i])
           for i in range(in_arr.rank)]

    # if in_arr Tensor, then res_arr will be Tensor, else it will be Arraypy
    if isinstance(in_arr, Tensor):
        res_arr = Tensor(Arraypy(arg), in_arr.ind_char)
    else:
        res_arr = Arraypy(arg)

    signs = [0 for i in range(factorial(in_arr.rank))]
    temp_i = 0
    for p in permutations(range(in_arr.rank)):
        signs[temp_i] = perm_parity(list(p))
        temp_i += 1

    index = [in_arr.start_index[i] for i in range(in_arr.rank)]

    for i in range(len(in_arr)):
        perm = list(permutations(index))
        perm_number = 0
        for temp_index in perm:
            res_arr[tuple(index)] += signs[perm_number] * \
                in_arr[tuple(temp_index)]
            perm_number += 1
        if isinstance(res_arr[tuple(index)], int):
            res_arr[tuple(index)] = float(res_arr[tuple(index)])
        res_arr[tuple(index)] /= factorial(in_arr.rank)

        index = in_arr.next_index(index)

    return res_arr


def perm_parity(lst):
    '''\
    THANKS TO Paddy McCarthy FROM http://code.activestate.com/ FOR THIS FUNCTION!
    Given a permutation of the digits 0..N in order as a list,
    returns its parity (or sign): +1 for even parity; -1 for odd.
            
    Examples
    ========
    
    >>> signs=zeros(6)
    >>> temp_i=0
    >>> for p in permutations(range(3)):
            signs[temp_i]=perm_parity(list(p))
            print(signs[temp_i], p)
            temp_i+=1

   (1.0, (0, 1, 2))
   (-1.0, (0, 2, 1))
   (-1.0, (1, 0, 2))
   (1.0, (1, 2, 0))
   (1.0, (2, 0, 1))
   (-1.0, (2, 1, 0))
    '''
    parity = 1
    for i in range(0, len(lst) - 1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i, len(lst)), key=lst.__getitem__)
            lst[i], lst[mn] = lst[mn], lst[i]
    return parity
#====================================================
