from sympy.tensor.array.sp_utils import LilSparseArray, CooSparseArray, CsrSparseArray


a = [[1, 0, 1, 0],
     [0, 0, 0, 1],
     [1, 0, 0, 1],
     [0, 0, 1, 0]]

def test_sparse_array_format():
    for sp_format in [LilSparseArray, CooSparseArray, CsrSparseArray]:
        A = sp_format(a)
        assert len(A) == 16
        assert A.shape == (4, 4)
        assert A.rank() == 2
        assert A.tolist() == [[1, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 1], [0, 0, 1, 0]]
        assert (A + A).tolist() == [[2, 0, 2, 0], [0, 0, 0, 2], [2, 0, 0, 2], [0, 0, 2, 0]]

def test_lil_format():
    A = LilSparseArray(a)
    assert A._data == [[1, 1], [1], [1, 1], [1]]
    assert A._rows == [[0, 2], [3], [0, 3], [2]]
    B = LilSparseArray(([[1, 1], [1], [1, 1], [1]], [[0, 2], [3], [0, 3], [2]]), (4, 4))
    assert A.tolist() == B.tolist()

    C = B.tocsr()
    C.tolist() == CsrSparseArray(a).tolist()
    

def test_coo_format():
    A = CooSparseArray(a)
    assert A._data == [1, 1, 1, 1, 1, 1]
    assert A._coor == [[0, 0, 1, 2, 2, 3], [0, 2, 3, 0, 3, 2]]
    B = CooSparseArray(([1, 1, 1, 1, 1, 1], [[0, 0, 1, 2, 2, 3], [0, 2, 3, 0, 3, 2]]), (4, 4))
    assert A.tolist() == B.tolist()

    C = B.tocsr()
    C.tolist() == CsrSparseArray(a).tolist()


def test_csr_format():
    A = CsrSparseArray(a)
    assert A._data == [1, 1, 1, 1, 1, 1]
    assert A._col_ind == [0, 2, 3, 0, 3, 2]
    assert A._row_ptr == [0, 2, 3, 5, 7]
    B = CsrSparseArray(([1, 1, 1, 1, 1, 1], [0, 2, 3, 0, 3, 2], [0, 2, 3, 5, 7]), (4, 4))
    assert A.tolist() == B.tolist()
