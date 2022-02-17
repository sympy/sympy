from sympy.tensor.array.expressions.array_expressions import ArrayTensorProduct, ArrayContraction, ArrayAdd, \
    ArrayDiagonal
from sympy.tensor.array.expressions.conv_array_to_matrix import convert_array_to_matrix
from sympy.tensor.array.expressions.conv_indexed_to_array import convert_indexed_to_array
from sympy.tensor.array.expressions.conv_matrix_to_array import convert_matrix_to_array
from sympy.utilities.exceptions import sympy_deprecation_warning

sympy_deprecation_warning(
    """
    Array expressions inside the codegen module are deprecated. Use the
    experimental module in sympy.tensor.array.expressions instead
    """,
    deprecated_since_version="1.8",
    active_deprecations_target='deprecated-arrayexpressions',
)


CodegenArrayTensorProduct = ArrayTensorProduct
CodegenArrayContraction = ArrayContraction
CodegenArrayElementwiseAdd = ArrayAdd
CodegenArrayDiagonal = ArrayDiagonal

recognize_matrix_expression = convert_array_to_matrix
parse_matrix_expression = convert_matrix_to_array
parse_indexed_expression = convert_indexed_to_array
