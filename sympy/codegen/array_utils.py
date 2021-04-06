from sympy.tensor.array.expressions.array_expressions import ArrayTensorProduct, ArrayContraction, ArrayAdd, \
    ArrayDiagonal
from sympy.tensor.array.expressions.conv_array_to_matrix import convert_array_to_matrix
from sympy.tensor.array.expressions.conv_indexed_to_array import convert_indexed_to_array
from sympy.tensor.array.expressions.conv_matrix_to_array import convert_matrix_to_array
from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(
    feature="Array expressions inside the codegen module",
    useinstead="Experimental module in sympy.tensor.array.expressions",
    issue=20996,
    deprecated_since_version="1.8"
).warn()


CodegenArrayTensorProduct = ArrayTensorProduct
CodegenArrayContraction = ArrayContraction
CodegenArrayElementwiseAdd = ArrayAdd
CodegenArrayDiagonal = ArrayDiagonal

recognize_matrix_expression = convert_array_to_matrix
parse_matrix_expression = convert_matrix_to_array
parse_indexed_expression = convert_indexed_to_array
