from sympy import ImmutableMatrix as Matrix
from sympy.vector import CoordSys3D, matrix_to_vector

def test_matrix_to_vector_float_zero():
    C = CoordSys3D('C')
    m = Matrix([0.0, 2, 3])
    v = matrix_to_vector(m, C)
    assert v.to_matrix(C).equals(m)

