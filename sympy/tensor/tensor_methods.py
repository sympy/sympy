# -*- coding: utf-8 -*-

from itertools import permutations
from sympy.tensor.arraypy import Arraypy, Tensor, copy
from random import randint
from sympy.functions.combinatorial.factorials import factorial
from sympy import expand


def symmetric(in_arr):
    """
    Creates the symmetric form of input tensor.
    Input: Arraypy or Tensor with equal axes (array shapes).
    Output: symmetric array. Output type - Arraypy or Tensor, depends of input

    Examples
    ========

    >>> from sympy.tensor.arraypy import list2arraypy
    >>> from sympy.tensor.tensor_methods import symmetric
    >>> a = list2arraypy(list(range(9)), (3,3))
    >>> b = symmetric(a)
    >>> print (b)
    0  2.00000000000000  4.00000000000000
    2.00000000000000  4.00000000000000  6.00000000000000
    4.00000000000000  6.00000000000000  8.00000000000000
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

    >>> from sympy.tensor.arraypy import list2arraypy
    >>> from sympy.tensor.tensor_methods import asymmetric
    >>> a = list2arraypy(list(range(9)), (3,3))
    >>> b = asymmetric(a)
    >>> print (b)
    0  -1.00000000000000  -2.00000000000000
    1.00000000000000  0  -1.00000000000000
    2.00000000000000  1.00000000000000  0
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


def tensor_product(first_tensor, second_tensor):
    """Returns tensor product. Rank of resulted tensor is a summary rank of two
    tensors.

    Examples
    ========

    >>> from sympy.tensor.arraypy import Arraypy, Tensor
    >>> from sympy.tensor.tensor_methods import tensor_product
    >>> a = Tensor( Arraypy ('1..2', 'X'), 1)
    >>> b = Tensor( Arraypy ('1..2', 'Y'), -1)
    >>> c = tensor_product(a, b)
    >>> print (c)
    X[1]*Y[1] X[1]*Y[2]
    X[2]*Y[1] X[2]*Y[2]

    >>> print(c.start_index)
    (1, 1)
    >>> print(c.end_index)
    (2, 2)
    >>> print(c.rank)
    2
    >>> print(c.ind_char)
    (1, -1)

    """

    if not isinstance(first_tensor, Tensor):
        raise TypeError('first attribute must be Tensor')
    if not isinstance(second_tensor, Tensor):
        raise TypeError('second attribute must be Tensor')

    # forming list of tuples for Arraypy constructor of type a = Arraypy(
    # [(a, b), (c, d), ... , (y, z)] )
    arg = [(first_tensor.start_index[i], first_tensor.end_index[i])
           for i in range(first_tensor._rank)]
    arg = arg + [(second_tensor.start_index[i], second_tensor.end_index[i])
                 for i in range(second_tensor._rank)]

    # index charater of resulted tensor will be a concatination of two
    # index characters
    res = Tensor(Arraypy(arg), first_tensor.ind_char + second_tensor.ind_char)

    # start indexes of: current tensor and second_tensor tensor
    cur_idx = first_tensor.start_index
    second_tensor_idx = second_tensor.start_index

    # loop over current tensor
    for i in range(len(first_tensor)):
        # loop over second_tensor tensor
        for j in range(len(second_tensor)):
            res[cur_idx + second_tensor_idx] = first_tensor[cur_idx] * \
                second_tensor[second_tensor_idx]
            # second_tensor tensor next index
            second_tensor_idx = second_tensor.next_index(second_tensor_idx)
        # current tensor next index
        cur_idx = first_tensor.next_index(cur_idx)

    return res


def wedge(first_tensor, second_tensor):
    """
    Finds outer product (wedge) of two tensor arguments.
    The algoritm is too find the asymmetric form of tensor product of
    two arguments. The resulted array is multiplied on coefficient which is
    coeff = factorial(p+s)/factorial(p)*factorial(s)

    Examples
    ========

    """
    if not isinstance(first_tensor, Tensor):
        raise TypeError('Input must be of Tensor type')
    if not isinstance(second_tensor, Tensor):
        raise TypeError('Input must be of Tensor type')

    p = len(first_tensor)
    s = len(second_tensor)

    coeff = factorial(p + s) / factorial(p) * factorial(s)
    return coeff * asymmetric(tensor_product(first_tensor, second_tensor))


def lower_index(tensor, metric_tensor, *index_numbers_to_low):
    """Lowering one upper index of the tensor.
    The index count starts from 1.

    Examples
    ========
    >>> from sympy import symbols, sin
    >>> from sympy.tensor.arraypy import list2tensor
    >>> from sympy.tensor.tensor_methods import lower_index
    >>> x, y, z, w, r, phi = symbols('x y z w r phi')
    >>> A = list2tensor([1, 0, 0, 0, r**2, 0, 0, 0, (r**2)*sin(phi)], (3, 3) ,(-1, -1))
    >>> print(A)
    1  0  0
    0  r**2  0
    0  0  r**2*sin(phi)

    >>> T = list2tensor([w, x, 0, y, z, 0, 0, y**2, x*y*w], (3, 3) ,(1, -1))
    >>> print(T)
    w  x  0
    y  z  0
    0  y**2  w*x*y

    >>> S1 = lower_index(T, A, 1)
    >>> print(S1)
    w  x  0
    r**2*y  r**2*z  0
    0  r**2*y**2*sin(phi)  r**2*w*x*y*sin(phi)

    >>> print(S1.ind_char)
    (-1, -1)

    """
    index_numbers_to_low = list(index_numbers_to_low)
    for i in range(len(index_numbers_to_low)):
        index_numbers_to_low[i] -= 1

    if not isinstance(tensor, Tensor):
        raise TypeError('Input tensor must be of Tensor type')

    if not isinstance(metric_tensor, Tensor):
        raise TypeError('Metric tensor must be of Tensor type')

    if not metric_tensor.rank == 2:
        raise ValueError('Metric tensor rank must be equal to 2')

    if not metric_tensor.ind_char == (-1, -1):
        raise ValueError('Metric tensor must be covariant')

    for index_number in index_numbers_to_low:
        if tensor.ind_char[index_number] == -1:
            raise ValueError('Index number should point on upper index')

    # forming list of tuples for Arraypy constructor of type a = Arraypy( [(a,
    # b), (c, d), ... , (y, z)] )
    arg = [(tensor.start_index[i], tensor.end_index[i])
           for i in range(tensor.rank)]
    new_ind_char = [i for i in tensor.ind_char]
    for i in index_numbers_to_low:
        new_ind_char[i] = -1
    result_tensor = Tensor(Arraypy(arg), new_ind_char)

    for index_number in index_numbers_to_low:
        # loop over all tensor elements
        for cur_index in tensor.index_list:
            # loop over dimension pointed in index_number_to_low
            for j in range(
                    tensor.start_index[index_number],
                    tensor.end_index[index_number] + 1):
                # forming indexes
                metric_tensor_index = (j, cur_index[index_number])
                temp_index = [i for i in cur_index]
                temp_index[index_number] = j
                result_tensor[
                    cur_index] += tensor[tuple(temp_index)] * metric_tensor[metric_tensor_index]
        tensor = copy(result_tensor)

    return result_tensor


def raise_index(tensor, metric_tensor, *index_numbers_to_raise):
    """Raising one lower index of the tensor.
    The index count starts from 1.

    Examples
    ========
    >>> from sympy import symbols, sin
    >>> from sympy.tensor.arraypy import list2tensor
    >>> from sympy.tensor.tensor_methods import raise_index
    >>> x, y, z, w, r, phi = symbols('x y z w r phi')
    >>> A = list2tensor([1, 0, 0, 0, r**2, 0, 0, 0, (r**2)*sin(phi)], (3, 3) ,(-1, -1))
    >>> print(A)
    1  0  0
    0  r**2  0
    0  0  r**2*sin(phi)

    >>> T = list2tensor([w, x, 0, y, z, 0, 0, y**2, x*y*w], (3, 3) ,(1, -1))
    >>> print(T)
    w  x  0
    y  z  0
    0  y**2  w*x*y

    >>> S = raise_index(T, A, 2)
    >>> print(S)
    w  x/r**2  0
    y  z/r**2  0
    0  y**2/r**2  w*x*y/(r**2*sin(phi))

    >>> print(S.ind_char)
    (1, 1)

    """
    index_numbers_to_raise = list(index_numbers_to_raise)
    for i in range(len(index_numbers_to_raise)):
        index_numbers_to_raise[i] -= 1

    if not isinstance(tensor, Tensor):
        raise TypeError('Input tensor must be of Tensor type')

    if not isinstance(metric_tensor, Tensor):
        raise TypeError('Metric tensor must be of Tensor type')

    if not metric_tensor.rank == 2:
        raise ValueError('Metric tensor rank must be equal to 2')

    if not metric_tensor.ind_char == (-1, -1):
        raise ValueError('Metric tensor must be covariant')

    for index_number in index_numbers_to_raise:
        if tensor.ind_char[index_number] == 1:
            raise ValueError('Index number should point on lower index')

    # forming list of tuples for Arraypy constructor of type a = Arraypy( [(a,
    # b), (c, d), ... , (y, z)] )
    arg = [(tensor.start_index[i], tensor.end_index[i])
           for i in range(tensor.rank)]
    new_ind_char = [i for i in tensor.ind_char]
    for i in index_numbers_to_raise:
        new_ind_char[i] = 1
    result_tensor = Tensor(Arraypy(arg), new_ind_char)

    # finding contravariant reversed matrix
    metric_tensor = metric_tensor.to_matrix().inv()

    for index_number in index_numbers_to_raise:
        # loop over all tensor elements
        for cur_index in tensor.index_list:
            # loop over dimension pointed in index_number_to_low
            for j in range(
                    tensor.start_index[index_number],
                    tensor.end_index[index_number] + 1):
                # forming indexes
                metric_tensor_index = (j, cur_index[index_number])
                temp_index = [i for i in cur_index]
                temp_index[index_number] = j
                result_tensor[
                    cur_index] += tensor[tuple(temp_index)] * metric_tensor[metric_tensor_index]
        tensor = copy(result_tensor)

    return result_tensor


def change_basis(tensor, transformation_matrix, old_to_new=True):
    """
    WORK IN PROGRESS
    
    Changes basis of input Tensor with help of transformation matrix from old
    to new (if corresponding argument is True) or from new to old (if False).
    
    """
    if not isinstance(tensor, Tensor):
        raise TypeError('First argument must be of Tensor type')
    if not isinstance(transformation_matrix, Tensor):
        raise TypeError(
            'Transformation matrix must be of Tensor type')
    if transformation_matrix.rank != 2:
        raise ValueError('Transformation matrix must be of rank == 2')
    if transformation_matrix.ind_char != (1, -1):
        raise ValueError('Transformation matrix valency must be (1, -1)')

    inv_transformation_matrix = transformation_matrix.to_matrix().inv()
    # for testing reasons
    # inv_transformation_matrix = Tensor(Arraypy((2,2), 'B'), (-1, 1))

    # forming list of tuples for Arraypy constructor of type a = Arraypy( [(a,
    # b), (c, d), ... , (y, z)] )
    arg = [(tensor.start_index[i], tensor.end_index[i])
           for i in range(tensor.rank)]
    temp_tensor = Tensor(Arraypy(arg), tensor.ind_char)
    result_tensor = Tensor(Arraypy(arg), tensor.ind_char)

    # lists, that represents where is upper index and where is low
    upper_idx_numbers = [i for i in range(len(tensor.ind_char)) \
                         if tensor.ind_char[i] == 1]
    low_idx_numbers = [i for i in range(len(tensor.ind_char)) \
                         if tensor.ind_char[i] == -1]

    # summ over upper indicies
    for idx in tensor.index_list:
        for i in upper_idx_numbers:
            for j in range(tensor.start_index[i], tensor.end_index[i] + 1):
                    # forming index for multiplied tensor element
                    temp_idx = []
                    for k in range(len(idx)):
                        if k == i:
                            temp_idx.append(j)
                        else:
                            temp_idx.append(idx[k])
                    temp_tensor[idx] += transformation_matrix[(idx[i], j)] * tensor[tuple(temp_idx)]

    # summ over low indicies
    for idx in tensor.index_list:
        for i in low_idx_numbers:
            for j in range(tensor.start_index[i], tensor.end_index[i] + 1):
                # forming index for multiplied tensor element
                temp_idx = []
                for k in range(len(idx)):
                    if k == i:
                        temp_idx.append(j)
                    else:
                        temp_idx.append(idx[k])
                result_tensor[idx] += temp_tensor[tuple(temp_idx)] * inv_transformation_matrix[(idx[i], j)]
    
    # expanding for more clear view. Should I do this?
    for idx in result_tensor.index_list:
        result_tensor[idx] = expand(result_tensor[idx])
        
    return result_tensor


def perm_parity(lst):
    """\
    THANKS TO Paddy McCarthy FROM http://code.activestate.com/ FOR THIS FUNCTION!
    Given a permutation of the digits 0..N in order as a list,
    returns its parity (or sign): +1 for even parity; -1 for odd.

    Examples
    ========

    >>> from sympy.matrices import zeros
    >>> from itertools import permutations
    >>> from sympy.tensor.tensor_methods import perm_parity
    >>> signs=zeros(6)
    >>> temp_i=0
    >>> for p in permutations(range(3)):
    ...     signs[temp_i]=perm_parity(list(p))
    ...     print(str(signs[temp_i]) + ' ' + str(p))
    ...     temp_i+=1
    1 (0, 1, 2)
    -1 (0, 2, 1)
    -1 (1, 0, 2)
    1 (1, 2, 0)
    1 (2, 0, 1)
    -1 (2, 1, 0)
    """
    parity = 1
    for i in range(0, len(lst) - 1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i, len(lst)), key=lst.__getitem__)
            lst[i], lst[mn] = lst[mn], lst[i]
    return parity


def is_symmetric(array):
    """Check if array or tensor is already symmetric."""
    return array == symmetric(array)


def is_asymmetric(array):
    """Check if array or tensor is already asymmetric."""
    return array == asymmetric(array)
