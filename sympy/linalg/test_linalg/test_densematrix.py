from sympy.matrices.lilmatrix import LILMatrix

# Search and replace matrix to matrix

matrix = DenseMatrix

def test_setelement():
    '''__setitem__, element assignment'''
    A = matrix([[1,0,3],
                [0,0,1],
                [0,0,0]])
    A[2, 2] = 1
    assert A == matrix([[1,0,3],
                        [0,0,1],
                        [0,0,1]])
    A[2, 2] = 0
    A[1, 2] = 0
    assert A == matrix([[1,0,3],
                        [0,0,0],
                        [0,0,0]])
    A[0, 2] = 10
    A[2, 2] = 0
    assert A == matrix([[1,0,10],
                        [0,0,0],
                        [0,0,0]])
    A = matrix([[0,0,0],
                [0,0,0],
                [0,0,0]])
    A[1, 2] = 1
    assert A == matrix([[0,0,0],
                        [0,0,1],
                        [0,0,0]])
    A = matrix([[1,2,3],
                [4,5,6],
                [7,8,9]])
    A[1, 1] = 50
    A[0, 0] = 0
    assert A == matrix([[0,2,3],
                        [4,50,6],
                        [7,8,9]])

def test_getelement():
    '''__getitem__, element access'''
    A = matrix([[1,0,3],
                [0,0,1],
                [0,0,0]])
    assert A[2, 2] == 0
    assert A[0, 0] == 1
    assert A[0, 2] == 3
    assert A[1, 1] == 0
    assert A[1, 2] == 1
    assert A[2, 1] == 0
    A = matrix([[1,2,3],
                [4,5,6],
                [7,8,9]])
    assert A[0, 0] == 1
    assert A[1, 2] == 6
    assert A[2, 0] == 7


def test_submatrix():
    '''
    The A[:, :] syntax will no longer work in the low-level matrix.
    function matrix.get_submatrix(rlo, rhi, clo, chi) has been provided
    whcih returns the correspoding submatrix.
    '''
    A = matrix([[1,0,3],
                [0,0,1],
                [0,0,0]])
    A.submatrix(0, 2, 0, 2) == matrix([[1,0],[0,0]])
    A.submatrix(0, 1, 1, 2) == matrix([[0]])
    A.submatrix(0, 3, 0, 3) == A
    A.submatrix(1, 3, 1, 3) == matrix([[0, 1], [0, 0]])
    
    A = matrix([[1,2,3],[4,5,6],[7,8,9]])
    A.submatrix(0, 2, 0, 2) == matrix([[1,2],[4,5]])
    A.submatrix(0, 1, 1, 2) == matrix([[2]])
    A.submatrix(0, 3, 0, 3) == A
    A.submatrix(1, 3, 1, 3) == matrix([[5, 6], [8, 9]])
    
def test_setsubmatrix():
    # Is a function of 5 arguments acceptable ?
    pass

def test_addition():
    A = matrix([[1,0,3],[0,0,1],[0,0,0]])
    B = matrix([[0,1,1],[4,7,-1],[1,0,-1]])
    assert A + B == matrix([[1,1,4],[4,7,0],[1,0,-1]])
    C = matrix([[-1,0,-3],[0,0,-1],[0,0,0]])
    D = A + C
    assert D == matrix([[0,0,0],[0,0,0],[0,0,0]])
    print D, D.mat
    assert D.mat == [[],[],[]]
    # Add more ?

def test_subtraction():
    A = matrix([[0,1,1],[4,7,-1],[1,0,-1]])
    assert A - A == matrix([[0,0,0],[0,0,0],[0,0,0]])
    (A - A).mat == [[],[],[]]
    B = matrix([[-1,0,-3],[0,0,-1],[0,0,0]])
    assert A - B == matrix([[1,1,4],[4,7,0],[1,0,-1]])

def test_scalar_multiplication():
    A = matrix([[0,1,1],[4,7,-1],[1,0,-1]])
    assert 5 * A == matrix([[0,5,5],[20,35,-5],[5,0,-5]])
    assert 0 * A == matrix([[0,0,0],[0,0,0],[0,0,0]])

def test_matrix_multiplication():
    A = matrix([[0,1,1],[4,7,-1],[1,0,-1]])
    B = matrix([[-1,0,-3],[0,0,-1],[0,0,0]])
    assert A * B == matrix([[0,0,-1],[-4,0,-19],[-1,0,-3]])
    C = matrix([[7,-1,8],[-3,1,-4],[7,-1,4]])
    assert A * C == 4 * matrix([[1,0,0],[0,1,0],[0,0,1]])

def test_LU():
    A = matrix([[1,2,0],[3,6,-1],[1,2,1]])
    L, U = LUdecomp(A)
    assert L.is_lower()
    assert U.is_upper()
    assert L * U == A
    assert L == matrix([[1,0,0],[3,1,0],[1,-1,1]])
    assert U == matrix([[1,2,0],[0,0,-1],[0,0,0]])

def test_cholesky():
    A = matrix([[25,15,-5],[15,18,0],[-5,0,11]])
    L = cholesky(A)
    assert L.is_lower()
    assert L * L.T == A
    assert L == matrix([[5,3,-1],[0,3,1],[0,0,3]])

def test_LDL():
    A = matrix([[25,15,-5],[15,18,0],[-5,0,11]])
    L, D = LDLdecomp(A)
    assert L.is_lower()
    assert D.is_diagonal()
    assert L * D * L.T == A
    assert L == matrix([[1,0,0],[3/5,1,0],[-1/5,1/3,1]])
    assert D == matrix([[25,0,0],[0,9,0],[0,0,9]])

def test_det_GE():
    A = matrix([[3,1,8],[2,-5,4],[-1,6,-2]])
    assert A.det_GE() == 14
    A = matrix([[1,2,3],[4,5,6],[7,8,9]])
    assert A.det_GE() == 0
    A = matrix([[1,5,3],[6,8,3],[3,5,6]])
    assert A.det_GE() == -84

def test_LUsolve():
    A = matrix([[3,1,8],[2,-5,4],[-1,6,-2]])
    B = matrix([14,21,7])
    X = LUsolve(A, B)
    assert A * X == B
    assert X == matrix([83,5,-30])

def test_cholesky_solve():
    A = matrix([[3,1,8],[2,-5,4],[-1,6,-2]])
    B = matrix([14,21,7])
    X = cholesky_solve(A, B)
    assert A * X == B
    assert X == matrix([83,5,-30])

def test_LDLsolve():
    A = matrix([[3,1,8],[2,-5,4],[-1,6,-2]])
    B = matrix([14,21,7])
    X = LDLsolve(A, B)
    assert A * X == B
    assert X == matrix([83,5,-30])

# Merge the above three functions ?
    
