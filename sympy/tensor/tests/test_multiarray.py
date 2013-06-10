from sympy.core.numbers import Rational
from sympy.core.symbol import symbols
from sympy.matrices import Matrix
from sympy.tensor.multiarray import MultiArray, multiempty, multiones


# Define general variables:
array1 = [1, 2, 3, 4]
M = MultiArray.create(array1)
array2 = [[i * j for j in range(20)] for i in range(10)]
M2 = MultiArray.create(array2)
M3 = MultiArray.create(lambda x: x[0] ** 3 + 7, (5,))
mat = Matrix((
    (1, 2, 3, -1),
    (4, 5, 6, -2),
    (7, 8, 9, -3),
))
M4 = MultiArray.create(mat)
M5 = MultiArray.create(lambda x: x[0] + x[1] + x[2] + x[3], (5, 6, 7, 8,))
ex_sq_mat = Matrix([[i * j for j in range(1, 4)] for i in range(1, 4)])
M6 = MultiArray.from_matrix(ex_sq_mat)
M7 = MultiArray.create([[[i * j * k for k in range(1, 4)] for j in range(1, 4)] for i in range(1, 4)])
# rank-0 tests, creation:
M0 = MultiArray.create(lambda x: 10, [])


def test_multiarray_stuff():
    for i in xrange(len(array1)):
        assert M[i] == array1[i]
    assert M[:] == M
    assert list(M) == [1, 2, 3, 4]


def test_multiarray_get_matrix():
    # * a Matrix (add from Matrix and get_matrix)
    # * get a matrix when only two indices remain.
    # * get a matrix when a single index remains.

    # rank 0
    assert M0.get_matrix() == Matrix([10])

    # rank 1
    M_matrix = M.get_matrix()
    for i, el in enumerate([1, 2, 3, 4]):
        assert M_matrix[i, 0] == el
        assert M_matrix[i] == el

    # rank 2
    M2_matrix = M2.get_matrix()

    assert len(list(M2)) == 200

    for i in xrange(10):
        # assert M2[i]._data_list == array2[i]
        for j in xrange(20):
            assert M2[i, j] == array2[i][j]
            # assert M2[i][j] == array2[i][j]
            assert M2_matrix[i, j] == array2[i][j]

    M4_matrix = M4.get_matrix()

    for i in xrange(mat.shape[0]):
        for j in xrange(mat.shape[1]):
            assert M4[i, j] == mat[i, j]
            # assert M4[i][j] == mat[i, j]
            assert M4_matrix[i, j] == mat[i, j]


def test_multiarray_slice():
    M2_first_sliced = M2[:, 2]
    for i in range(10):
        assert M2_first_sliced[i] == M2[i, 2]

    for i in range(5):
        assert M3[i] == i ** 3 + 7

    # test slices:
    M5slice1 = M5[2, 3:5, 2:, 7]
    M5slice2 = M5[2, 3:-1, :, 7]
    M5slice3 = M5[2, -3:-1, -5:, 7]
    for i in xrange(2):
        for j in xrange(5):
                assert M5slice1[i, j] == M5[2, 3 + i, 2 + j, 7]
                assert M5slice2[i, j] == M5[2, 3 + i, j, 7]
                assert M5slice3[i, j] == M5[2, 3 + i, 2 + j, 7]
    assert M5.rank == 4

    assert M5slice1.rank == 2
    assert M5slice2.rank == 2
    assert M5slice3.rank == 2


def test_multiarray_creation():
    M5copy = MultiArray(M5)

    for i in xrange(5):
        for j in xrange(6):
            for k in xrange(7):
                for l in xrange(8):
                    assert M5copy[i, j, k, l] == M5[i, j, k, l]
                    assert M5[i, j, k, l] == i + j + k + l

    # TODO: create from
    # * ranges of sub MultiArray?
    # * fail on malformed lists, like [[1,2],3]
    # * fail on excess indices.
    # * way to get arbitrary subMultiArray.
    # * slices...
    # * MultiArray(function, dimensions)
    #   -> mismatch between func and dims (two kinds).
    # * negative index access
    # * slice with step

    # * test tensor_product
    # * test contract_positions
    #   if dimensions do not match.

    # * get_matrix when rank == 0
    # * __str__ when rank == 0
    # * __str__ when rank == 1 or 2.
    # rank 0 should be 0, not 1
    # TODO: contract() or trace() for Einstein summation?


def test_multiarray_index_contraction():
    try:
        M5.contract_positions(1, 2)
    except ValueError:
        pass
    else:
        raise ValueError("This contraction should have raised a ValueError.")

    M6_contracted = M6.contract_positions(0, 1)

    assert M6.rank == 2
    assert M6_contracted.rank == 0
    assert M6_contracted[0] == ex_sq_mat.trace()
    # assert M6_cont... == trace?

    M7c1 = M7.contract_positions(0, 1)

    assert M7c1.rank == 1

    for i in range(3):
        assert M7c1[i] == 14 * (i + 1)
    try:
        M7.contract_positions(1, 1)
    except ValueError:
        pass
    else:
        raise ValueError("Expression did not raise an exception.")


def test_multiarray_symmetry():
    # test symmetries
    return
    # TODO: this is still unimplemented, remove this test?
    S1 = MultiArray.create([[1, 2], [2, 1]], sym=[[1, 1]])
    A1 = MultiArray.create([[0, 3], [-3, 0]], sym=[[2]])

    S2 = MultiArray.create([[1, 2, 3], [2, 4, 5], [3, 5, 6]], sym=[[1, 1]])
    A2 = MultiArray.create([[0, 1, 2], [-1, 0, 3], [-2, -3, 0]], sym=[[2]])


def test_multiarray_rank_0():
    assert M0.rank == 0
    assert M0[()] == 10

# * test multiarray contract all positions
#   -> give rank == 0
#   -> M0[()]


def test_multiarray_inequalities():
    # ASSERT INEQUALITIES:
    assert M != M2
    assert M2 != M3
    assert M3 != M4
    assert M4 != M5
    assert M5 != M6


def test_multiarray_tensor_product():
    # test tensor product:
    M14 = M.tensor_product(M4)

    for i in xrange(M.dimensions[0]):
        for j in xrange(M4.dimensions[0]):
            for k in xrange(M4.dimensions[1]):
                assert M14[i, j, k] == M[i] * M4[j, k]


def test_multiarray_mixed_expressions_with_scalars():
    # multiplication of MultiArray with scalar:
    x1, x2, x3 = symbols('x1 x2 x3')
    assert M[0] != 0  # otherwise the following tests are useless
    assert (M * 3)[0] == M[0] * 3
    assert (4 * M)[0] == 4 * M[0]

    assert (x1 * M)[0] == x1 * M[0]
    assert (M * x1)[0] == M[0] * x1

    exprm2 = x1 * M2 * 3 * x2 * (-4)

    assert exprm2[1, 1] == -12 * x1 * x2
    assert exprm2[8, 8] == -768 * x1 * x2
    assert exprm2 == x1 * M2 * 3 * x2 * (-4)

    exprm3 = M / 4

    assert exprm3[0] == Rational(1, 4)
    assert exprm3 == M / 4
    assert exprm3 * 4 == M
    assert 4 * exprm3 == M


def test_multiarray_add_radd_sub_rsub():
    traM = M + 3
    for i in range(4):
        assert traM[i] == M[i] + 3

    tramM = M - 2
    for i in range(4):
        assert tramM[i] == M[i] - 2

    rtraM = 7 + M
    for i in range(4):
        assert rtraM[i] == M[i] + 7
    assert rtraM - 7 == M

    rtramM = 5 - M
    for i in range(4):
        assert rtramM[i] == 5 - M[i]
    assert rtramM + M == 5 * multiones(1, 4)

    # excepted to fail: 4 / M TODO
    try:
        div4m = 4 / M
    except Exception:
        pass
    else:
        raise ValueError("division by MultiArray did not raise an Exception!")

    x1, x2, x3, x4 = symbols('x1 x2 x3 x4')
    tM = M + x1
    assert tM - x1 == M
    for i in range(4):
        assert tM[i] == M[i] + x1

    tmM = M - x2
    assert tmM + x2 == M
    for i in range(4):
        assert tmM[i] == M[i] - x2

    rtM = x3 - M
    assert rtM + M == x3 * multiones(1, 4)
    for i in range(4):
        assert rtM[i] == x3 - M[i]

    rtmM = M + x1
    assert rtmM - x1 == M
    for i in range(4):
        assert rtmM[i] == M[i] + x1

    # MISSING, TODO: addition and subtraction of MultiArray with scalars, numbers and symbols.


def test_multiarray_applyfunc():
    newm = M.applyfunc(lambda x: x**4 + 17)
    for i in xrange(newm.dimensions[0]):
        assert newm[i] == M[i]**4 + 17

    newm2 = M2.applyfunc(lambda x: x**3 - 97)
    for i in xrange(newm2.dimensions[0]):
        for j in xrange(newm2.dimensions[1]):
            assert newm2[i, j] == M2[i, j]**3 - 97


def test_multiarray_add_sub_with_noncommutative():
    x1, x2, x3 = symbols('x1, x2, x3')
    nx1, nx2, nx3 = symbols('nx1, nx2, nx3', commutative=False)

    assert M * x1 == x1 * M
    assert nx2 * M * nx1 != nx1 * M * nx2

    ncM = MultiArray.create([[nx1, nx2], [nx3, x3]])

    assert nx1 * ncM != ncM * nx1
    assert x1 * ncM == ncM * x1


def test_self_extract():
    pass
    Mextr = [M.self_extract(0, _) for _ in range(4)]
    for i in range(4):
        assert Mextr[i] == M[i]
    M2extr1 = [M2.self_extract(0, _) for _ in range(10)]
    M2extr2 = [M2.self_extract(1, _) for _ in range(20)]
    for i in range(10):
        for j in range(20):
            assert M2extr1[i][j] == M2[i, j]
            assert M2extr2[j][i] == M2[i, j]


def test_get_indexwise_linear_transformation():
    pass
    m6eq1 = M6.get_indexwise_linear_transformation(None, None)
    m6eq2 = M6.get_indexwise_linear_transformation(None, MultiArray.create([1]*3))
    m6eq3 = M6.get_indexwise_linear_transformation(MultiArray.create([1]*3), None)

    m6extr1 = M6.get_indexwise_linear_transformation(MultiArray.create([-1]*3), None)
    m6extr2 = M6.get_indexwise_linear_transformation(None, MultiArray.create([-1]*3))
    m6extr3 = M6.get_indexwise_linear_transformation(MultiArray.create([-1]*3), MultiArray.create([-1]*3))

    m6extr4 = M6.get_indexwise_linear_transformation(None, MultiArray.create([_**2 - 104 for _ in range(3)]))
    m6extr5 = M6.get_indexwise_linear_transformation(MultiArray.create([_**3 + 93 for _ in range(3)]), None)

    for i in range(3):
        for j in range(3):
            assert m6eq1[i, j] == M6[i, j]
            assert m6eq2[i, j] == M6[i, j]
            assert m6eq3[i, j] == M6[i, j]
            assert m6extr1[i, j] == -M6[i, j]
            assert m6extr2[i, j] == -M6[i, j]
            assert m6extr3[i, j] == M6[i, j]
            assert m6extr4[i, j] == (j**2 - 104) * M6[i, j]
            assert m6extr5[i, j] == (i**3 + 93) * M6[i, j]


def test_multiempty():
    E1 = multiempty(5)
    E2 = multiempty(4, 5)
    E3 = multiempty(7, 6, 5, 4, 3, 2, 1)
    for i in xrange(5):
        assert E1[i] == 0
    for i in xrange(4):
        for j in xrange(5):
                assert E2[i, j] == 0
    for i0 in xrange(7):
        for i1 in xrange(6):
            for i2 in xrange(5):
                for i3 in xrange(4):
                    for i4 in xrange(3):
                        for i5 in xrange(2):
                            for i6 in xrange(1):
                                assert E3[i0, i1, i2, i3, i4, i5, i6] == 0
