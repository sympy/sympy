
from sympy.core import *
from sympy.matrices import *

def test_multiplication():
    a=Matrix((
        (1, 2),
        (3, 1),
        (0, 6),
        ))

    b = Matrix ((
        (1, 2),
        (3, 0),
        ))

    c= a*b
    assert c[0,0]==7
    assert c[0,1]==2
    assert c[1,0]==6
    assert c[1,1]==6
    assert c[2,0]==18
    assert c[2,1]==0

    x = Symbol("x")

    c = b * Symbol("x")
    assert isinstance(c,Matrix)
    assert c[0,0] == x
    assert c[0,1] == 2*x
    assert c[1,0] == 3*x
    assert c[1,1] == 0

    c = 5 * b
    assert isinstance(c,Matrix)
    assert c[0,0] == 5
    assert c[0,1] == 2*5
    assert c[1,0] == 3*5
    assert c[1,1] == 0

def test_power():
    A = Matrix([[2,3],[4,5]])
    assert (A**5)[:] == [6140, 8097, 10796, 14237]
    A = Matrix([[2, 1, 3],[4,2, 4], [6,12, 1]])
    assert (A**3)[:] == [290, 262, 251, 448, 440, 368, 702, 954, 433]

def test_creation():
    x = Symbol("x")
    a = Matrix([x, 0], [0, 0])
    m = a
    assert m.cols == m.lines
    assert m.cols == 2
    assert m[:] == [x,0,0,0]
    b = Matrix(2,2, [x, 0, 0, 0])
    m = b
    assert m.cols == m.lines
    assert m.cols == 2
    assert m[:] == [x,0,0,0]

    assert a == b

def test_determinant():
    x, y = Symbol('x'), Symbol('y')

    assert Matrix([ [1] ]).det() == 1

    assert Matrix(( (-3,  2),
                    ( 8, -5) )).det() == -1

    assert Matrix(( (x,   1),
                    (y, 2*y) )).det() == 2*x*y-y

    assert Matrix(( (1, 1, 1),
                    (1, 2, 3),
                    (1, 3, 6) )).det() == 1

    assert Matrix(( ( 3, -2,  0, 5),
                    (-2,  1, -2, 2),
                    ( 0, -2,  5, 0),
                    ( 5,  0,  3, 4) )).det() == -289

    assert Matrix(( ( 1,  2,  3,  4),
                    ( 5,  6,  7,  8),
                    ( 9, 10, 11, 12),
                    (13, 14, 15, 16) )).det() == 0

    assert Matrix(( (3, 2, 0, 0, 0),
                    (0, 3, 2, 0, 0),
                    (0, 0, 3, 2, 0),
                    (0, 0, 0, 3, 2),
                    (2, 0, 0, 0, 3) )).det() == 275

    assert Matrix(( (1, 0,  1,  2, 12),
                    (2, 0,  1,  1,  4),
                    (2, 1,  1, -1,  3),
                    (3, 2, -1,  1,  8),
                    (1, 1,  1,  0,  6) )).det() == -55

    assert Matrix(( (-5,  2,  3,  4,  5),
                    ( 1, -4,  3,  4,  5),
                    ( 1,  2, -3,  4,  5),
                    ( 1,  2,  3, -2,  5),
                    ( 1,  2,  3,  4, -1) )).det() == 11664

    assert Matrix(( ( 2,  7, -1, 3, 2),
                    ( 0,  0,  1, 0, 1),
                    (-2,  0,  7, 0, 2),
                    (-3, -2,  4, 5, 3),
                    ( 1,  0,  0, 0, 1) )).det() == 123

def test_submatrix():
    m0 = eye(4)
    assert m0[0:3, 0:3] == eye(3)
    assert m0[2:4, 0:2] == zero(2)

    m1 = Matrix(3,3, lambda i,j: i+j)
    assert m1[0,:] == Matrix(1,3,(0,1,2))
    assert m1[1:3, 1] == Matrix(2,1,(2,3))

    m2 = Matrix([0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15])
    assert m2[:,-1] == Matrix(4,1,[3,7,11,15])
    assert m2[-2:,:] == Matrix([[8,9,10,11],[12,13,14,15]])

def test_submatrix_assignment():
    m = zero(4)
    m[2:4, 2:4] = eye(2)
    assert m == Matrix((0,0,0,0),
                        (0,0,0,0),
                        (0,0,1,0),
                        (0,0,0,1))
    m[0:2, 0:2] = eye(2)
    assert m == eye(4)
    m[:,0] = Matrix(4,1,(1,2,3,4))
    assert m == Matrix((1,0,0,0),
                        (2,1,0,0),
                        (3,0,1,0),
                        (4,0,0,1))
    m[:,:] = zero(4)
    assert m == zero(4)
    m[:,:] = ((1,2,3,4),(5,6,7,8),(9, 10, 11, 12),(13,14,15,16))
    assert m == Matrix(((1,2,3,4),
                        (5,6,7,8),
                        (9, 10, 11, 12),
                        (13,14,15,16)))
    m[0:2, 0] = [0,0]
    assert m == Matrix(((0,2,3,4),
                        (0,6,7,8),
                        (9, 10, 11, 12),
                        (13,14,15,16)))

def test_reshape():
    m0 = eye(3)
    assert m0.reshape(1,9) == Matrix(1,9,(1,0,0,0,1,0,0,0,1))
    m1 = Matrix(3,4, lambda i,j: i+j)
    assert m1.reshape(4,3) == Matrix((0,1,2), (3,1,2), (3,4,2), (3,4,5))
    assert m1.reshape(2,6) == Matrix((0,1,2,3,1,2), (3,4,2,3,4,5))

def test_applyfunc():
    m0 = eye(3)
    assert m0.applyfunc(lambda x:2*x) == eye(3)*2
    assert m0.applyfunc(lambda x: 0 ) == zero(3)

def test_random():
    M = randMatrix(3,3)
    M = randMatrix(3,3,seed=3)
    M = randMatrix(3,4,0,150)

def test_LUdecomp():
    testmat = Matrix([[0,2,5,3],
                      [3,3,7,4],
                      [8,4,0,2],
                      [-2,6,3,4]])
    L,U,p = testmat.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-testmat == zero(4)

    testmat = Matrix([[6,-2,7,4],
                      [0,3,6,7],
                      [1,-2,7,4],
                      [-9,2,6,3]])
    L,U,p = testmat.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-testmat == zero(4)

    x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
    M = Matrix((1, x, 1), (2, y, 0), (y, 0, z))
    L, U, p = M.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-M == zero(3)

    # test FF LUdecomp
    M = Matrix([[1, 3, 3],
                [3, 2, 6],
                [3, 2, 2]])
    P, L, Dee, U = M.LUdecompositionFF()
    assert P*M == L*Dee.inv()*U

    M = Matrix([[1, 2, 3, 4],
                 [3, -1, 2, 3],
                 [3, 1, 3, -2],
                 [6, -1, 0, 2]])
    P, L, Dee, U = M.LUdecompositionFF()
    assert P*M == L*Dee.inv()*U



def test_LUsolve():
    A = Matrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    x = Matrix(3,1,[3,7,5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x
    A = Matrix([[0,-1,2],
                [5,10,7],
                [8,3,4]])
    x = Matrix(3,1,[-1,2,5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x

def test_inverse():
    A = eye(4)
    assert A.inv() == eye(4)
    A = Matrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    Ainv = A.inv()
    assert A*Ainv == eye(3)
    Ainv2 = A.inv("LU")
    assert Ainv == Ainv2

def test_util():
    v1 = Matrix(1,3,[1,2,3])
    v2 = Matrix(1,3,[3,4,5])
    assert v1.cross(v2) == Matrix(1,3,[-2,4,-2])
    assert v1.norm() == sqrt(14)
    # cofactor
    assert eye(3) == eye(3).cofactorMatrix()
    test = Matrix([[1,3,2],[2,6,3],[2,3,6]])
    assert test.cofactorMatrix() == Matrix([[27,-6,-6],[-12,2,3],[-3,1,0]])
    test = Matrix([[1,2,3],[4,5,6],[7,8,9]])
    assert test.cofactorMatrix() == Matrix([[-3,6,-3],[6,-12,6],[-3,6,-3]])

def test_jacobian_hessian():
    x = Symbol('x')
    y = Symbol('y')
    L = Matrix(1,2,[x**2*y, 2*y**2 + x*y])
    syms = [x,y]
    assert L.jacobian(syms) == Matrix([[2*x*y, x**2],[y, 4*y+x]])

    L = Matrix(1,2,[x, x**2*y**3])
    assert L.jacobian(syms) == Matrix([[1, 0], [2*x*y**3, x**2*3*y**2]])

    f = x**2*y
    syms = [x,y]
    assert hessian(f, syms) == Matrix([[2*y, 2*x], [2*x, 0]])

    f = x**2*y**3
    assert hessian(f, syms) == Matrix([[2*y**3, 6*x*y**2],[6*x*y**2, 6*x**2*y]])


def test_QR():
    A = Matrix([[1,2],[2,3]])
    Q, S = A.QRdecomposition()
    R = Rational
    assert Q == Matrix([[5**R(-1,2), (R(2)/5)*(R(1)/5)**R(-1,2)], [2*5**R(-1,2), (-R(1)/5)*(R(1)/5)**R(-1,2)]])
    assert S == Matrix([[5**R(1,2), 8*5**R(-1,2)], [0, (R(1)/5)**R(1,2)]])
    assert Q*S == A
    assert Q.T * Q == eye(2)

    A = Matrix([[1,1,1],[1,1,3],[2,3,4]])
    Q, R = A.QRdecomposition()
    assert Q*Q.T == eye(Q.lines)
    assert R.is_upper()
    assert A == Q*R

def test_nullspace():
    # first test reduced row-ech form
    R = Rational

    M = Matrix([[5,7,2,1],
               [1,6,2,-1]])
    out, tmp = M.rref()
    assert out == Matrix([[1,0,-R(2)/23,R(13)/23],
                              [0,1,R(8)/23, R(-6)/23]])

    M = Matrix([[1,3,0,2,6,3,1],
                [-2,-6,0,-2,-8,3,1],
                [3,9,0,0,6,6,2],
                [-1,-3,0,1,0,9,3]])
    out, tmp = M.rref()
    assert out == Matrix([[1,3,0,0,2,0,0],
                               [0,0,0,1,2,0,0],
                               [0,0,0,0,0,1,R(1)/3],
                               [0,0,0,0,0,0,0]])
    # now check the vectors
    basis = M.nullspace()
    assert basis[0] == Matrix([-3,1,0,0,0,0,0])
    assert basis[1] == Matrix([0,0,1,0,0,0,0])
    assert basis[2] == Matrix([-2,0,0,-2,1,0,0])
    assert basis[3] == Matrix([0,0,0,0,0,R(-1)/3, 1])

def test_eigen():
    # test charpoly
    x = Symbol('x')
    y = Symbol('y')
    eye3 = eye(3)
    assert eye3.charpoly(x) == (1-x)**3
    assert eye3.charpoly(y) == (1-y)**3
    # test values
    M = Matrix([(0,1,-1),
                (1,1,0),
                (-1,0,1) ])
    vals = M.eigenvals()
    vals.sort()
    assert vals == [[-1,1],[1,1],[2,1]]
    # test vectors
    R = Rational
    M = Matrix([ [1,0,0],
                 [0,1,0],
                 [0,0,1]])
    assert M.eigenvects() == [[1, 3, [Matrix([1,0,0]), Matrix([0,1,0]), Matrix([0,0,1])]]]
    M = Matrix([ [5,0,2],
                 [3,2,0],
                 [0,0,1]])
    assert M.eigenvects() == [[1, 1, [Matrix([R(-1)/2,R(3)/2,1])]],
                              [2, 1, [Matrix([0,1,0])]],
                              [5, 1, [Matrix([1,1,0])]]]

def test_sparse_matrix():
    return
    def eye(n):
        tmp = SMatrix(n,n,lambda i,j:0)
        for i in range(tmp.lines):
            tmp[i,i] = 1
        return tmp
    def zero(n):
        return SMatrix(n,n,lambda i,j:0)

    # test_multiplication
    a=SMatrix((
        (1, 2),
        (3, 1),
        (0, 6),
        ))

    b = SMatrix ((
        (1, 2),
        (3, 0),
        ))

    c= a*b
    assert c[0,0]==7
    assert c[0,1]==2
    assert c[1,0]==6
    assert c[1,1]==6
    assert c[2,0]==18
    assert c[2,1]==0

    x = Symbol("x")

    c = b * Symbol("x")
    assert isinstance(c,SMatrix)
    assert c[0,0] == x
    assert c[0,1] == 2*x
    assert c[1,0] == 3*x
    assert c[1,1] == 0

    c = 5 * b
    assert isinstance(c,SMatrix)
    assert c[0,0] == 5
    assert c[0,1] == 2*5
    assert c[1,0] == 3*5
    assert c[1,1] == 0

    #test_power
    A = SMatrix([[2,3],[4,5]])
    assert (A**5)[:] == [6140, 8097, 10796, 14237]
    A = SMatrix([[2, 1, 3],[4,2, 4], [6,12, 1]])
    assert (A**3)[:] == [290, 262, 251, 448, 440, 368, 702, 954, 433]


    # test_creation
    x = Symbol("x")
    a = SMatrix([x, 0], [0, 0])
    m = a
    assert m.cols == m.lines
    assert m.cols == 2
    assert m[:] == [x,0,0,0]
    b = SMatrix(2,2, [x, 0, 0, 0])
    m = b
    assert m.cols == m.lines
    assert m.cols == 2
    assert m[:] == [x,0,0,0]

    assert a == b

    # test_determinant
    x, y = Symbol('x'), Symbol('y')

    assert SMatrix([ [1] ]).det() == 1

    assert SMatrix(( (-3,  2),
                    ( 8, -5) )).det() == -1

    assert SMatrix(( (x,   1),
                    (y, 2*y) )).det() == 2*x*y-y

    assert SMatrix(( (1, 1, 1),
                    (1, 2, 3),
                    (1, 3, 6) )).det() == 1

    assert SMatrix(( ( 3, -2,  0, 5),
                    (-2,  1, -2, 2),
                    ( 0, -2,  5, 0),
                    ( 5,  0,  3, 4) )).det() == -289

    assert SMatrix(( ( 1,  2,  3,  4),
                    ( 5,  6,  7,  8),
                    ( 9, 10, 11, 12),
                    (13, 14, 15, 16) )).det() == 0

    assert SMatrix(( (3, 2, 0, 0, 0),
                    (0, 3, 2, 0, 0),
                    (0, 0, 3, 2, 0),
                    (0, 0, 0, 3, 2),
                    (2, 0, 0, 0, 3) )).det() == 275

    assert SMatrix(( (1, 0,  1,  2, 12),
                    (2, 0,  1,  1,  4),
                    (2, 1,  1, -1,  3),
                    (3, 2, -1,  1,  8),
                    (1, 1,  1,  0,  6) )).det() == -55

    assert SMatrix(( (-5,  2,  3,  4,  5),
                    ( 1, -4,  3,  4,  5),
                    ( 1,  2, -3,  4,  5),
                    ( 1,  2,  3, -2,  5),
                    ( 1,  2,  3,  4, -1) )).det() == 11664

    assert SMatrix(( ( 2,  7, -1, 3, 2),
                    ( 0,  0,  1, 0, 1),
                    (-2,  0,  7, 0, 2),
                    (-3, -2,  4, 5, 3),
                    ( 1,  0,  0, 0, 1) )).det() == 123

    # test_submatrix
    m0 = eye(4)
    assert m0[0:3, 0:3] == eye(3)
    assert m0[2:4, 0:2] == zero(2)

    m1 = SMatrix(3,3, lambda i,j: i+j)
    assert m1[0,:] == SMatrix(1,3,(0,1,2))
    assert m1[1:3, 1] == SMatrix(2,1,(2,3))

    m2 = SMatrix([0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15])
    assert m2[:,-1] == SMatrix(4,1,[3,7,11,15])
    assert m2[-2:,:] == SMatrix([[8,9,10,11],[12,13,14,15]])

    # test_submatrix_assignment
    m = zero(4)
    m[2:4, 2:4] = eye(2)
    assert m == SMatrix((0,0,0,0),
                        (0,0,0,0),
                        (0,0,1,0),
                        (0,0,0,1))
    m[0:2, 0:2] = eye(2)
    assert m == eye(4)
    m[:,0] = SMatrix(4,1,(1,2,3,4))
    assert m == SMatrix((1,0,0,0),
                        (2,1,0,0),
                        (3,0,1,0),
                        (4,0,0,1))
    m[:,:] = zero(4)
    assert m == zero(4)
    m[:,:] = ((1,2,3,4),(5,6,7,8),(9, 10, 11, 12),(13,14,15,16))
    assert m == SMatrix(((1,2,3,4),
                        (5,6,7,8),
                        (9, 10, 11, 12),
                        (13,14,15,16)))
    m[0:2, 0] = [0,0]
    assert m == SMatrix(((0,2,3,4),
                        (0,6,7,8),
                        (9, 10, 11, 12),
                        (13,14,15,16)))

    # test_reshape
    m0 = eye(3)
    assert m0.reshape(1,9) == SMatrix(1,9,(1,0,0,0,1,0,0,0,1))
    m1 = SMatrix(3,4, lambda i,j: i+j)
    assert m1.reshape(4,3) == SMatrix((0,1,2), (3,1,2), (3,4,2), (3,4,5))
    assert m1.reshape(2,6) == SMatrix((0,1,2,3,1,2), (3,4,2,3,4,5))

    # test_applyfunc
    m0 = eye(3)
    assert m0.applyfunc(lambda x:2*x) == eye(3)*2
    assert m0.applyfunc(lambda x: 0 ) == zero(3)

    # test_LUdecomp
    testmat = SMatrix([[0,2,5,3],
                      [3,3,7,4],
                      [8,4,0,2],
                      [-2,6,3,4]])
    L,U,p = testmat.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-testmat == zero(4)

    testmat = SMatrix([[6,-2,7,4],
                      [0,3,6,7],
                      [1,-2,7,4],
                      [-9,2,6,3]])
    L,U,p = testmat.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-testmat == zero(4)

    x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
    M = Matrix((1, x, 1), (2, y, 0), (y, 0, z))
    L, U, p = M.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-M == zero(3)

    # test_LUsolve
    A = SMatrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    x = SMatrix(3,1,[3,7,5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x
    A = SMatrix([[0,-1,2],
                [5,10,7],
                [8,3,4]])
    x = SMatrix(3,1,[-1,2,5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x

    # test_inverse
    A = eye(4)
    assert A.inv() == eye(4)
    A = SMatrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    Ainv = A.inv()
    assert A*Ainv == eye(3)
    Ainv2 = A.inv("LU")
    assert Ainv == Ainv2

    # test_cross
    v1 = Matrix(1,3,[1,2,3])
    v2 = Matrix(1,3,[3,4,5])
    assert v1.cross(v2) == Matrix(1,3,[-2,4,-2])
    assert v1.norm(v1) == 14

    # test_cofactor
    assert eye(3) == eye(3).cofactorMatrix()
    test = SMatrix([[1,3,2],[2,6,3],[2,3,6]])
    assert test.cofactorMatrix() == SMatrix([[27,-6,-6],[-12,2,3],[-3,1,0]])
    test = SMatrix([[1,2,3],[4,5,6],[7,8,9]])
    assert test.cofactorMatrix() == SMatrix([[-3,6,-3],[6,-12,6],[-3,6,-3]])

    # test_jacobian
    x = Symbol('x')
    y = Symbol('y')
    L = SMatrix(1,2,[x**2*y, 2*y**2 + x*y])
    syms = [x,y]
    assert L.jacobian(syms) == Matrix([[2*x*y, x**2],[y, 4*y+x]])

    L = SMatrix(1,2,[x, x**2*y**3])
    assert L.jacobian(syms) == SMatrix([[1, 0], [2*x*y**3, x**2*3*y**2]])

    # test_QR
    A = Matrix([[1,2],[2,3]])
    Q, S = A.QRdecomposition()
    R = Rational
    assert Q == Matrix([[5**R(-1,2), (R(2)/5)*(R(1)/5)**R(-1,2)], [2*5**R(-1,2), (-R(1)/5)*(R(1)/5)**R(-1,2)]])
    assert S == Matrix([[5**R(1,2), 8*5**R(-1,2)], [0, (R(1)/5)**R(1,2)]])
    assert Q*S == A
    assert Q.T * Q == eye(2)

    # test nullspace
    # first test reduced row-ech form
    R = Rational

    M = Matrix([[5,7,2,1],
               [1,6,2,-1]])
    out, tmp = M.rref()
    assert out == Matrix([[1,0,-R(2)/23,R(13)/23],
                              [0,1,R(8)/23, R(-6)/23]])

    M = Matrix([[1,3,0,2,6,3,1],
                [-2,-6,0,-2,-8,3,1],
                [3,9,0,0,6,6,2],
                [-1,-3,0,1,0,9,3]])
    out, tmp = M.rref()
    assert out == Matrix([[1,3,0,0,2,0,0],
                               [0,0,0,1,2,0,0],
                               [0,0,0,0,0,1,R(1)/3],
                               [0,0,0,0,0,0,0]])
    # now check the vectors
    basis = M.nullspace()
    assert basis[0] == Matrix([[-3,1,0,0,0,0,0]])
    assert basis[1] == Matrix([[0,0,1,0,0,0,0]])
    assert basis[2] == Matrix([[-2,0,0,-2,1,0,0]])
    assert basis[3] == Matrix([[0,0,0,0,0,R(-1)/3, 1]])


    # test eigen
    x = Symbol('x')
    y = Symbol('y')
    eye3 = eye(3)
    assert eye3.charpoly(x) == (1-x)**3
    assert eye3.charpoly(y) == (1-y)**3
    # test values
    M = Matrix([(0,1,-1),
                (1,1,0),
                (-1,0,1) ])
    vals = M.eigenvals()
    vals.sort()
    assert vals == [-1, 1, 2]

    R = Rational
    M = Matrix([ [1,0,0],
                 [0,1,0],
                 [0,0,1]])
    assert M.eigenvects() == [[1, 3, [Matrix(1,3,[1,0,0]), Matrix(1,3,[0,1,0]), Matrix(1,3,[0,0,1])]]]
    M = Matrix([ [5,0,2],
                 [3,2,0],
                 [0,0,1]])
    assert M.eigenvects() == [[1, 1, [Matrix(1,3,[R(-1)/2,R(3)/2,1])]],
                              [2, 1, [Matrix(1,3,[0,1,0])]],
                              [5, 1, [Matrix(1,3,[1,1,0])]]]

def test_subs():
    x = Symbol('x')
    assert Matrix([[1,x],[x,4]]).subs(x, 5) == Matrix([[1,5],[5,4]])
