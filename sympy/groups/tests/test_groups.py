import sympy
from sympy.matrices import Matrix
from sympy.combinatorics.matrix_groups import spl_linear_group
from sympy.combinatorics.matrix_groups import spl_orthogonal_group
from sympy.combinatorics.matrix_groups import finite_matrix_group

def test_Groups():
    m = sympy.Matrix(spl_linear_group(2))
    assert m.det() == 1
    m = sympy.Matrix(spl_orthogonal_group(2, 0.1))
    assert m.det() == 1
    mat = sympy.Matrix(finite_matrix_group(2, 3))
    assert max(mat) < 3
    mat = sympy.Matrix(finite_matrix_group(4, 11))
    assert max(mat) < 11
