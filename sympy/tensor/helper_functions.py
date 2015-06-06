# -*- coding: utf-8 -*-

from sympy.tensor.tensor_methods import is_symmetric, is_asymmetric
from sympy.tensor.arraypy import Arraypy, TensorArray
from sympy.matrices import Matrix


def check_vector_of_arguments(args):
    """The function contains checks for a one-dimensional list of arguments

    """

    if not isinstance(args, (list, TensorArray, Arraypy)):
        raise TypeError(
            "The type of vector of arguments must be list, TensorArray or \
            Arraypy")
    if isinstance(args, (TensorArray, Arraypy)):
        if len(args.shape) != 1:
            raise ValueError("The dimension of vector of arguments must be 1")
        if isinstance(args, TensorArray):
            if args.type_pq != (1, 0):
                raise ValueError(
                    "The valency of vector of arguments must be (+1)")


def check_metric_tensor(g):
    """The function contains checks for a metric tensor."""

    if not isinstance(g, (TensorArray, Matrix, Arraypy)):
        raise ValueError(
            "Type of metric tensor must be Matrix, TensorArray or Arraypy")
    if isinstance(g, (TensorArray, Arraypy)):
        if not is_symmetric(g):
            raise ValueError("The metric is not symmetric")
        if isinstance(g, (TensorArray)) and g.type_pq != (0, 2):
            raise ValueError(
                "The indices of metric tensor must be (-1,-1)")
    if isinstance(g, Matrix):
        if not g.is_symmetric:
            raise ValueError("The metric is not symmetric")
        
        
def check_the_christoffel_symbols_2(ch_2):
    """The function contains checks for a christoffel symbols of second kind."""
    if not isinstance(ch_2, (Arraypy, TensorArray)):
        raise TypeError(
            'The type of Christoffel symbol of second kind must be Arraypy \
            or TensorArray')
    if isinstance(ch_2, TensorArray) and ch_2.type_pq != (1, 2):
            raise ValueError(
                'The valence of Christoffel symbol of second \
                 kind must be (1, -1, -1)')


def check_the_vector_field(X):
    """The function contains checks for a vector field."""

    if not isinstance(X, (list, TensorArray, Arraypy)):
        raise ValueError(
            "The type of vector field must be list, TensorArray or Arraypy")
    if isinstance(X, (TensorArray, Arraypy)):
        if len(X.shape) != 1:
            raise ValueError("The dimension of vector field must be 1")
        if isinstance(X, TensorArray):
            if X.type_pq != (1, 0):
                raise ValueError("The valency of vector field must be (+1)")


def delete_index_from_list(_data, index):
    """The function takes the arguments '_data', which can be list or tuple. It
    returns tuple after remove the element standing in the position 'index'.

    Examples:
    =========

    >>> from sympy.tensor.helper_functions import delete_index_from_list
    >>> print(delete_index_from_list([10, 20, 30, 40, 50, 60], 3))
    (10, 20, 30, 50, 60)
    >>> print(delete_index_from_list((10, 20, 30, 40, 50, 60), 0))
    (20, 30, 40, 50, 60)

    """
    res = []
    for i in range(len(_data)):
        if i != index:
            res.append((_data[i]))
    return (tuple(res))


def replace_index_to_k(_data, index, k):
    """The function takes the arguments '_data' it's list or tuple,
    index-it's integer and k-any element.
    It returns tuple after replaces the item "index" on the element "k".

    Examples:
    =========

    >>> from sympy.tensor.helper_functions import replace_index_to_k
    >>> print(replace_index_to_k([10, 20, 30, 40], 0, 100))
    (100, 20, 30, 40)
    >>> print(replace_index_to_k((10, 20, 30, 40), 2, 100))
    (10, 20, 100, 40)

    """
    output = []
    for i in range(len(_data)):
        if i != index:
            output.append((_data[i]))
        else:
            output.append(k)
    return (tuple(output))


def sign_permutations(lst):
    """Return the sign of index.

    Examales:
    =========

    >>> from sympy.tensor.helper_functions import sign_permutations
    >>> print(sign_permutations([0, 1, 0]))
    0
    >>> print(sign_permutations([1, 1]))
    0
    >>> print(sign_permutations([3, 7, 10, 5]))
    1
    >>> print(sign_permutations([10, 5, 8, 11, 6]))
    -1

    """
    list_sort = sorted(lst)
    if (list_sort) == lst:
        if len(set(lst)) != len(lst):
            parity = 0
        else:
            parity = 1
    else:
        if len(set(lst)) != len(lst):
            parity = 0
        else:
            p = 0
            while len(lst) != 1:
                pos = lst.index(min(lst))
                p += pos
                lst.remove(min(lst))
            if p % 2 == 0:
                parity = 1
            else:
                parity = -1
    return parity
