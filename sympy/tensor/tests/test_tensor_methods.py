# -*- coding: utf-8 -*-

from sympy.tensor.arraypy import Arraypy, TensorArray, list2arraypy, \
    list2tensor, matrix2arraypy, matrix2tensor
from sympy.tensor.tensor_methods import symmetric, asymmetric, is_symmetric, \
    is_asymmetric, change_basis, tensor_product, wedge, raise_index,\
    lower_index
from sympy import Symbol, symbols, sin
from copy import copy

arr = list2arraypy(list(range(9)), (3, 3))
# 0 1 2
# 3 4 5
# 6 7 8

brr = list2arraypy([1 for i in range(9)], (3, 3))
# 1 1 1
# 1 1 1
# 1 1 1


def test_symmetric():
    # Arraypy
    sym_arr = symmetric(arr)
    # 0.0 2.0 4.0
    # 2.0 4.0 6.0
    # 4.0 6.0 8.0
    assert sym_arr[0, 1] == sym_arr[1, 0]
    assert sym_arr[0, 2] == sym_arr[1, 1] == sym_arr[2, 0]
    assert sym_arr[2, 1] == sym_arr[1, 2]

    # TensorArray
    tensor = arr.to_tensor((1, 1))
    sym_tensor = asymmetric(tensor)

    assert sym_tensor[0, 1] == -sym_tensor[1, 0]
    assert sym_tensor[0, 2] == -sym_tensor[2, 0]
    assert sym_tensor[2, 1] == -sym_tensor[1, 2]


def test_assymteric():
    # Arraypy
    asym_arr = asymmetric(arr)
    # 0.0 -1.0 -2.0
    # 1.0 0.0 -1.0
    # 2.0 1.0 0.0
    assert asym_arr[0, 1] == -asym_arr[1, 0]
    assert asym_arr[0, 2] == -asym_arr[2, 0]
    assert asym_arr[2, 1] == -asym_arr[1, 2]

    # TensorArray
    tensor = arr.to_tensor((1, 1))
    asym_tensor = asymmetric(tensor)

    assert asym_tensor[0, 1] == -asym_tensor[1, 0]
    assert asym_tensor[0, 2] == -asym_tensor[2, 0]
    assert asym_tensor[2, 1] == -asym_arr[1, 2]


def test_input_and_output_arguments():
    sym_arr = symmetric(arr)
    asym_arr = asymmetric(arr)
    assert isinstance(sym_arr, Arraypy)
    assert isinstance(asym_arr, Arraypy)

    tensor = arr.to_tensor((1, 1))
    sym_tensor = symmetric(tensor)
    asym_tensor = asymmetric(tensor)
    assert isinstance(sym_tensor, TensorArray)
    assert isinstance(asym_tensor, TensorArray)


def test_is_asymmetric():
    a = list2arraypy([0, -3, 5,
                      3, 0, 0,
                      -5, 0, 0], (3, 3))
    assert is_asymmetric(a)
    assert not is_symmetric(a)


def test_is_symmetric():
    a = list2arraypy([0, 3, 5,
                      3, 0, 0,
                      5, 0, 0], (3, 3))
    assert not is_asymmetric(a)
    assert is_symmetric(a)


def test_tensor_product():
    t_arr = arr.to_tensor((1, -1))
    t_brr = brr.to_tensor((-1, 1))

    t_prod = tensor_product(t_arr, t_brr)
    assert t_prod.shape == t_arr.shape + t_brr.shape == (3, 3, 3, 3)
    assert t_prod.ind_char == t_arr.ind_char + t_brr.ind_char == (1, -1, -1, 1)
    assert t_prod[1, 1, 2, 2] == t_arr[1, 1] * t_brr[2, 2]


def test_raise_index():
    x, y, z, w, r, phi = symbols('x y z w r phi')
    A = list2tensor([1, 0, 0,
                     0, r**2, 0,
                     0, 0, (r**2) * sin(phi)], (3, 3), (-1, -1))
    original_A = copy(A)

    T = list2tensor([w, x, 0,
                     y, z, 0,
                     0, y**2, x * y * w], (3, 3), (-1, -1))
    original_T = copy(T)

    S1 = raise_index(T, A, 1)
    assert S1.ind_char == (1, -1)

    # check that original A and T not changed
    assert original_A == A
    assert original_T == T

    S2 = raise_index(T, A, 2)
    # w  x/r**2  0
    # y  z/r**2  0
    # 0  y**2/r**2  w*x*y/(r**2*sin(phi))
    assert S2.ind_char == (-1, 1)
    assert S2[0, 0] == w
    assert S2[0, 1] == x / r**2
    assert S2[0, 2] == S2[1, 2] == S2[2, 0] == 0
    assert S2[1, 0] == y
    assert S2[1, 1] == z / r**2
    assert S2[2, 1] == y**2 / r**2
    assert S2[2, 2] == w * x * y / (r**2 * sin(phi))

    S3 = raise_index(T, A, 1, 2)
    S4 = raise_index(T, A, *[1, 2])
    assert S3 == S4
    assert S3.ind_char == S4.ind_char == (1, 1)


def test_lower_index():
    x, y, z, w, r, phi = symbols('x y z w r phi')
    A = list2tensor([1, 0, 0,
                     0, r**2, 0,
                     0, 0, (r**2) * sin(phi)], (3, 3), (-1, -1))
    original_A = copy(A)

    T = list2tensor([w, x, 0,
                     y, z, 0,
                     0, y**2, x * y * w], (3, 3), (1, 1))
    original_T = copy(T)

    S1 = lower_index(T, A, 1)
    assert S1.ind_char == (-1, 1)

    assert original_A == A
    assert original_T == T

    # w  x  0
    # r**2*y  r**2*z  0
    # 0  r**2*y**2*sin(phi)  r**2*w*x*y*sin(phi)
    assert S1[0, 0] == w
    assert S1[0, 1] == x
    assert S1[0, 2] == S1[1, 2] == S1[2, 0] == 0
    assert S1[1, 0] == r**2 * y
    assert S1[1, 1] == r**2 * z
    assert S1[2, 1] == r**2 * y**2 * sin(phi)
    assert S1[2, 2] == r**2 * w * x * y * sin(phi)

    S2 = lower_index(T, A, 2)
    assert S2.ind_char == (1, -1)

    S3 = lower_index(T, A, 1, 2)
    S4 = lower_index(T, A, *[1, 2])
    assert S3 == S4
    assert S3.ind_char == S4.ind_char == (-1, -1)
