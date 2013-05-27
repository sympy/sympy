from sympy.matrices import Matrix
from sympy.tensor.vtensor import VTensorIndexType, vtensorhead
from sympy.tensor.tensor import tensor_indices
from sympy.core.symbol import symbols
from sympy.functions.elementary.miscellaneous import sqrt
from sympy import eye


minkowski = Matrix((
    (1, 0, 0, 0),
    (0, -1, 0, 0),
    (0, 0, -1, 0),
    (0, 0, 0, -1),
))
Lorentz = VTensorIndexType('Lorentz', minkowski)
i0, i1, i2, i3, i4 = tensor_indices('i0:5', Lorentz)
E, px, py, pz = symbols('E px py pz')
A = vtensorhead('A', [Lorentz], [[1]], values=[E, px, py, pz])
B = vtensorhead('B', [Lorentz], [[1]], range(4), 'Gcomm')
AB = vtensorhead("AB", [Lorentz] * 2, [[1], [1]], values=minkowski)
ba_matrix = Matrix((
    (1, 2, 3, 4),
    (5, 6, 7, 8),
    (9, 0, -1, -2),
    (-3, -4, -5, -6),
))
BA = vtensorhead("BA", [Lorentz] * 2, [[1], [1]], values=ba_matrix)
# Let's test the diagonal metric, with inverted Minkowski metric:
LorentzD = VTensorIndexType('LorentzD', [-1, 1, 1, 1])
mu0, mu1, mu2 = tensor_indices('mu0:3', LorentzD)
C = vtensorhead('C', [LorentzD], [[1]], values=[E, px, py, pz])


def test_vtensor_iter():
    # iteration on VTensorHead
    assert list(A) == [E, px, py, pz]
    assert list(ba_matrix) == list(BA)

    # iteration on VTensMul
    assert list(A(i1)) == [E, px, py, pz]
    assert list(BA(i1, i2)) == list(ba_matrix)
    assert list(3 * BA(i1, i2)) == [3 * i for i in list(ba_matrix)]
    assert list(-5 * BA(i1, i2)) == [-5 * i for i in list(ba_matrix)]

    # iteration on VTensAdd
    assert list(A(i1) + A(i1)) == [2*E, 2*px, 2*py, 2*pz]
    assert list(BA(i1, i2) - BA(i1, i2)) == [0] * len(list(ba_matrix))
    assert list(BA(i1, i2) - 2 * BA(i1, i2)) == [-i for i in list(ba_matrix)]


def test_vtensor_covariant_contravariant_elements():
    # TODO: handle metric with covariant/contravariant indices!
    assert A(-i0)[0] == A(i0)[0]
    assert A(-i0)[1] == -A(i0)[1]

    assert AB(i0, i1)[1, 1] == -1
    assert AB(i0, -i1)[1, 1] == 1
    assert AB(-i0, -i1)[1, 1] == -1
    assert AB(-i0, i1)[1, 1] == 1


def test_vtensor_get_matrix():
    matab = AB(i0, i1).get_matrix()
    assert matab == Matrix([
                            [1,  0,  0,  0],
                            [0, -1,  0,  0],
                            [0,  0, -1,  0],
                            [0,  0,  0, -1],
                            ])
    # when alternating contravariant/covariant with [1, -1, -1, -1] metric
    # it becomes the identity matrix:
    assert AB(i0, -i1).get_matrix() == eye(4)

    # covariant and contravariant forms:
    assert A(i0).get_matrix() == Matrix([E, px, py, pz])
    assert A(-i0).get_matrix() == Matrix([E, -px, -py, -pz])
    # VTensorHead should also be endowed with a get_matrix method?
    # VTensMul

test_vtensor_get_matrix()


def test_vtensor_contraction():
    # TODO: decide how contract indices should work.
    assert A(i0) * A(-i0) == E ** 2 - px ** 2 - py ** 2 - pz ** 2
    assert A(i0) * A(-i0) == A ** 2
    assert A(i0) * A(-i0) == A(i0) ** 2
    assert (A(i0) * B(-i0)) == -px - 2 * py - 3 * pz

    for i in range(4):
        for j in range(4):
            assert (A(i0) * B(-i1))[i, j] == [E, px, py, pz][i] * [0, -1, -2, -3][j]

    # test contraction on the alternative Minkowski metric: [-1, 1, 1, 1]
    assert C(mu0) * C(-mu0) == -E ** 2 + px ** 2 + py ** 2 + pz ** 2

    contrexp = A(i0) * AB(i1, -i0)
    assert A(i0).rank == 1
    assert AB(i1, -i0).rank == 2
    assert contrexp.rank == 1
    for i in range(4):
        assert contrexp[i] == [E, px, py, pz][i]


def test_vtensor_self_contraction():
    AB(i0, -i0)
test_vtensor_self_contraction()


def test_vtensor_pow():
    assert C ** 2 == -E ** 2 + px ** 2 + py ** 2 + pz ** 2
    assert C ** 1 == sqrt(-E ** 2 + px ** 2 + py ** 2 + pz ** 2)


def test_vtensor_expressions():
    x1, x2, x3 = symbols('x1:4')

    # test coefficient in contraction:
    rank2coeff = x1 * A(i3) * B(i2)
    assert rank2coeff[1, 1] == x1 * px
    assert rank2coeff[3, 3] == 3 * pz * x1
    coeff_expr = (x1 * A(i4)) * (B(-i4) / x2)

    assert coeff_expr == -px * x1 / x2 - 2 * py * x1 / x2 - 3 * pz * x1 / x2

    add_expr = A(i0) + B(i0)

    assert add_expr[0] == E
    assert add_expr[1] == px + 1
    assert add_expr[2] == py + 2
    assert add_expr[3] == pz + 3

    sub_expr = A(i0) - B(i0)

    assert sub_expr[0] == E
    assert sub_expr[1] == px - 1
    assert sub_expr[2] == py - 2
    assert sub_expr[3] == pz - 3

    assert add_expr * B(-i0) == -px - 2 * py - 3 * pz - 14

    expr1 = x1 * A(i0) + x2 * B(i0)
    expr2 = expr1 * B(i1) * (-4)
    expr3 = expr2 + 3 * x3 * AB(i0, i1)
    expr4 = expr3 / 2
    assert expr4 * 2 == expr3
    expr5 = expr4 * BA(-i1, -i0)

    # TODO: check numerically that this value is correct:
    assert expr5 == (-2 * x1 * (-20 * E + 44 * px - 8 * py - 32 * pz) + 136 * x2 + 3 * x3).expand()
    # test contraction in more steps, verify that _vmetric is still there.


def test_noncommuting_components():
    euclid = VTensorIndexType('Euclidean', [1, 1])
    i1, i2, i3 = tensor_indices('i1:4', euclid)

    a, b, c, d = symbols('a b c d', commutative=False)
    V1 = vtensorhead('V1', [euclid] * 2, [[1], [1]], values=[[a, b], (c, d)])
    V2 = vtensorhead('V2', [euclid] * 2, [[1], [1]], values=[[a, c], [b, d]])

    vtp = V1(i1, i2) * V2(-i2, -i1)

    assert vtp == a ** 2 + b * c + c * b + d ** 2
    assert vtp != a**2 + 2*b*c + d**2

    # TODO: test non-commutative scalar coefficients in Mul, Add, subtractions...
    Vc = b * V1(i1, -i1)
    assert Vc.expand() == b * a + b * d


def test_vtensor_tensor_index_type_metric():
    # TODO: metric rank mismatch.
    # TODO: metric not the same as
    pass
