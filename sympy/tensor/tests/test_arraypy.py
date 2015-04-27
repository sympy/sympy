# -*- coding: utf-8 -*-

from sympy.tensor.arraypy import Arraypy, Tensor, list2arraypy, \
    list2tensor, matrix2arraypy, matrix2tensor
from sympy.tensor.tensor_methods import tensor_product
from sympy import Symbol, symbols
from sympy.matrices import Matrix, MatrixSymbol
import sys


def test_arraypy_initiation():
    arr_with_one_element = Arraypy()
    assert len(arr_with_one_element) == 1
    assert arr_with_one_element[0] == 0
    assert arr_with_one_element.rank == 1

    arr_with_symbol_element = Arraypy('Py')
    assert len(arr_with_symbol_element) == 1
    assert arr_with_symbol_element[0] == Symbol('Py[0]')
    assert arr_with_symbol_element.rank == 1

    vector_length = 5
    vector = Arraypy(vector_length)
    assert len(vector) == vector_length
    assert vector.shape == (vector_length,)
    assert vector.start_index == (0,)
    assert vector.end_index == (vector_length - 1,)
    assert vector.rank == 1

    array_shape = (3, 3, 3, 3)
    n_dim_array = Arraypy(array_shape)
    assert len(n_dim_array) == 3 * 3 * 3 * 3
    assert n_dim_array.shape == array_shape
    assert n_dim_array.start_index == (0, 0, 0, 0)
    assert n_dim_array.end_index == (2, 2, 2, 2)
    assert n_dim_array.rank == 4

    sparse_array = Arraypy(array_shape, 'sparse')
    assert sparse_array._sparse is True
    assert len(sparse_array._output) == 0
    assert len(sparse_array) == 3 * 3 * 3 * 3
    assert n_dim_array.shape == array_shape
    assert n_dim_array.start_index == (0, 0, 0, 0)
    assert n_dim_array.end_index == (2, 2, 2, 2)
    assert n_dim_array.rank == 4

    arr_with_ranged_index = Arraypy('1..3, 2 .. 4,   3..5')
    assert arr_with_ranged_index.shape == (3, 3, 3)
    assert len(arr_with_ranged_index) == 3 * 3 * 3
    assert arr_with_ranged_index.start_index == (1, 2, 3)
    assert arr_with_ranged_index.end_index == (3, 4, 5)
    assert arr_with_ranged_index.rank == 3

    combined_arg = [2, 3, 1]
    array_with_combined_arg = Arraypy(combined_arg)
    assert len(array_with_combined_arg) == 3 * 3
    assert array_with_combined_arg.shape == (3, 3)
    assert array_with_combined_arg.start_index == (1, 1)
    assert array_with_combined_arg.end_index == (3, 3)
    assert array_with_combined_arg.rank == 2

    shape = (3, 3)
    array_with_many_arg = Arraypy(shape, 'X', 'sparse')
    assert len(array_with_many_arg) == 3 * 3
    assert array_with_many_arg.shape == shape
    assert array_with_many_arg._sparse is True
    assert array_with_many_arg[0, 0] == Symbol('X[0, 0]')
    assert array_with_many_arg.rank == 2


def test_reshape():
    array = Arraypy(50)
    assert array.shape == (50,)
    assert array.rank == 1

    array = array.reshape((5, 5, 2))
    assert array.shape == (5, 5, 2)
    assert array.rank == 3
    assert len(array) == 50


def test_next_index():
    array = Arraypy((2, 2), 'Py')
    assert array.next_index((0, 0)) == (0, 1)
    assert array.next_index((0, 1)) == (1, 0)
    assert array.next_index((1, 0)) == (1, 1)
    assert array.next_index((1, 1)) == (0, 0)


def test_iterator():
    if sys.version_info[0] >= 3:
        array = list2arraypy([0, 1, 2, 3], (2, 2))
        j = 0
        for i in array:
            assert i == j
            j += 1
    
        array = array.reshape(4)
        j = 0
        for i in array:
            assert i == j
            j += 1


def test_sparse():
    sparse_array = Arraypy((2, 2), 'sparse')
    assert len(sparse_array) == 2 * 2
    # dictionary where all data is
    assert len(sparse_array._output) == 0
    # it's empty, even thought Arraypy knows that 'empty' data is zero
    
    if sys.version_info[0] >= 3:
        for i in sparse_array:
            assert i == 0
    else:
        idx = sparse_array.start_index
        for i in len(sparse_array.start_index):
            assert sparse_array[idx] == 0
            idx = sparse_array.next_index(idx)

    sparse_array[0, 0] = 123
    assert len(sparse_array._output) == 1
    assert sparse_array[0, 0] == 123

    # when element in sparse array become zero it will disappear from
    # dictionary
    sparse_array[0, 0] = 0
    assert len(sparse_array._output) == 0
    assert sparse_array[0, 0] == 0


def test_calculation():
    # Arraypy
    list_of_ones = [1 for i in range(9)]
    list_of_nines = [9 for i in range(9)]
    shape = (3, 3)

    a = list2arraypy(list_of_ones, shape)
    b = list2arraypy(list_of_nines, shape)
    
    if sys.version_info[0] >= 3:
        c = a + b
        for i in c:
            assert i == 10
    
        c = b - a
        for i in c:
            assert i == 8
    
    else:        
        c = a + b
        idx = c.start_index
        for i in len(c):
            assert c[idx] == 10
            idx = c.next_index(idx)
            
        idx = c.start_index
        c = a + b
        for i in len(c):
            assert c[idx] == 8
            idx = c.next_index(idx)        

    # Tensor
    x0, x1, y0, y1 = symbols('X[0], X[1], Y[0], Y[1]')
    tensor1 = Tensor(Arraypy(2, 'X'), 1)
    tensor2 = Tensor(Arraypy(2, 'Y'), -1)
    assert tensor1.rank == tensor2.rank == 1

    res_tensor = tensor_product(tensor1, tensor2)
    assert len(res_tensor) == 4
    assert res_tensor.rank == 2
    assert res_tensor[0, 0] == x0 * y0
    assert res_tensor[0, 1] == x0 * y1
    assert res_tensor[1, 0] == x1 * y0
    assert res_tensor[1, 1] == x1 * y1


def test_tensor_contract():
    tensor = list2tensor(list(range(9)), (3, 3), (1, -1))
    # 0 1 2
    # 3 4 5
    # 6 7 8
    assert tensor.rank == 2

    # contract of matrix will be summ of diagonal elements: 0 + 4 + 8 == 12
    res_tensor = tensor.contract(1, 2)
    assert res_tensor.rank == 1
    assert len(res_tensor) == 1
    assert res_tensor[0] == 12


def test_arraypy_converting():
    arr_arraypy = list2arraypy([1, 2, 3, 4], (2, 2))
    arr_list = arr_arraypy.to_list()
    assert (isinstance(arr_list, list))

    arr_tensor = arr_arraypy.to_tensor((1, -1))
    assert (isinstance(arr_tensor, Tensor))

    arr_matrix = arr_arraypy.to_matrix()
    assert (isinstance(arr_matrix, Matrix))

    idx = (0, 0)
    for i in range(len(arr_arraypy)):
        assert arr_arraypy[idx] == arr_tensor[idx] == arr_matrix[idx]
        idx = arr_arraypy.next_index(idx)


def test_tensor_initiation():
    base_shape = (2, 2, 2)
    tensor_base = Arraypy(base_shape)
    index_character = (1, -1, 1)
    tensor = Tensor(tensor_base, index_character)

    assert tensor.ind_char == index_character
    assert tensor.shape == base_shape
    assert tensor.type_pq == (2, 1)
    assert tensor.rank == 3


def test_tensor_converting():
    arr_tensor = Tensor(Arraypy((2, 2)), (1, 1))

    arr_list = arr_tensor.to_list()
    assert (isinstance(arr_list, list))

    arr_arraypy = arr_tensor.to_arraypy()
    assert (isinstance(arr_arraypy, Arraypy))

    arr_matrix = arr_tensor.to_matrix()
    assert (isinstance(arr_matrix, Matrix))


def test_converting_functions():
    arr_list = [1, 2, 3, 4]
    arr_matrix = Matrix(((1, 2), (3, 4)))

    # list
    arr_arraypy = list2arraypy(arr_list, (2, 2))
    assert (isinstance(arr_arraypy, Arraypy))

    arr_tensor = list2tensor(arr_list, (2, 2), (-1, -1))
    assert (isinstance(arr_tensor, Tensor))

    # Matrix
    arr_arraypy = matrix2arraypy(arr_matrix)
    assert (isinstance(arr_arraypy, Arraypy))

    arr_tensor = matrix2tensor(arr_matrix, (-1, -1))
    assert (isinstance(arr_tensor, Tensor))


def test_equality():
    first_list = [1, 2, 3, 4]
    second_list = [1, 2, 3, 4]
    third_list = [4, 3, 2, 1]
    assert first_list == second_list
    assert first_list != third_list

    first_arraypy = list2arraypy(first_list, (2, 2))
    second_arraypy = list2arraypy(second_list, (2, 2))
    third_arraypy = list2arraypy(third_list, (2, 2))
    fourth_arraypy = list2arraypy(first_list, 4)

    assert first_arraypy == second_arraypy
    second_arraypy[0, 0] = 0
    assert first_arraypy != second_arraypy
    assert first_arraypy != third_arraypy
    assert first_arraypy != fourth_arraypy

    first_tensor = list2tensor(first_list, (2, 2), (1, 1))
    second_tensor = list2tensor(second_list, (2, 2), (1, 1))
    third_tensor = list2tensor(second_list, (2, 2), (-1, 1))

    assert first_tensor == second_tensor
    assert first_tensor != third_tensor
