from sympy.tensor.arraypy import Arraypy, Tensor, list2arraypy, \
     list2tensor, matrix2arraypy, matrix2tensor
from sympy.tensor.tensor_methods import symmetric, asymmetric
from sympy import Symbol, symbols

arr = list2arraypy(list(range(9)), (3,3))
# 0 1 2
# 3 4 5
# 6 7 8

def test_symmetric():
    # Arraypy
    sym_arr = symmetric(arr)
    # 0.0 2.0 4.0
    # 2.0 4.0 6.0
    # 4.0 6.0 8.0
    assert sym_arr[0, 1] == sym_arr[1, 0]
    assert sym_arr[0, 2] == sym_arr[1, 1] == sym_arr[2, 0]
    assert sym_arr[2, 1] == sym_arr[1, 2]
    
    # Tensor
    tensor = arr.to_tensor((1,1))
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
    
    # Tensor
    tensor = arr.to_tensor((1,1))
    asym_tensor = asymmetric(tensor)
    
    assert asym_tensor[0, 1] == -asym_tensor[1, 0]
    assert asym_tensor[0, 2] == -asym_tensor[2, 0]
    assert asym_tensor[2, 1] == -asym_arr[1, 2]     

def test_input_and_output_arguments():
    sym_arr = symmetric(arr)
    asym_arr = asymmetric(arr)
    assert isinstance(sym_arr, Arraypy)
    assert isinstance(asym_arr, Arraypy)
    
    tensor = arr.to_tensor((1,1))
    sym_tensor = symmetric(tensor)
    asym_tensor = asymmetric(tensor)
    assert isinstance(sym_tensor, Tensor)
    assert isinstance(asym_tensor, Tensor)