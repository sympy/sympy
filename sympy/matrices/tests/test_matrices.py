from sympy import symbols, Matrix, eye, I, Symbol, Rational, wronskian, cos, \
        sin, exp, hessian, sqrt, zeros, ones, randMatrix, Poly, S, pi, \
        integrate, oo, raises, trigsimp, Integer, block_diag
from sympy.matrices.matrices import ShapeError, MatrixError
from sympy.printing import srepr
from sympy.utilities.pytest import XFAIL

def test_division():
    x, y, z = symbols('x','y','z')
    v = Matrix(1,2,[x, y])
    assert v.__div__(z) == Matrix(1,2,[x/z, y/z])
    assert v.__truediv__(z) == Matrix(1,2,[x/z, y/z])
    assert v/z == Matrix(1,2,[x/z, y/z])

def test_sum():
    x, y, z = symbols('xyz')
    m = Matrix([[1,2,3],[x,y,x],[2*y,-50,z*x]])
    assert m+m == Matrix([[2,4,6],[2*x,2*y,2*x],[4*y,-100,2*z*x]])

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

    c2 = x * b
    assert c == c2

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
    assert A**0 == eye(3)
    assert A**1 == A
    assert (Matrix([[2]]) ** 100)[0,0] == 2**100
    assert eye(2)**10000000 == eye(2)
    assert Matrix([[1, 2], [3, 4]])**Integer(2) == Matrix([[7, 10], [15, 22]])

def test_creation():
    raises(MatrixError, 'Matrix(5,5,range(20))')

    x = Symbol("x")
    a = Matrix([[x, 0], [0, 0]])
    m = a
    assert m.cols == m.rows
    assert m.cols == 2
    assert m[:] == [x,0,0,0]
    b = Matrix(2,2, [x, 0, 0, 0])
    m = b
    assert m.cols == m.rows
    assert m.cols == 2
    assert m[:] == [x,0,0,0]

    assert a == b

    assert Matrix(b) == b

def test_tolist():
    x, y, z = symbols('xyz')
    lst = [[S.One,S.Half,x*y,S.Zero],[x,y,z,x**2],[y,-S.One,z*x,3]]
    m = Matrix(lst)
    assert m.tolist() == lst

def test_determinant():
    x, y, z = Symbol('x'), Symbol('y'), Symbol('z')

    M = Matrix((1,))

    assert M.det(method="bareis") == 1
    assert M.det(method="berkowitz") == 1

    M = Matrix(( (-3,  2),
                 ( 8, -5) ))

    assert M.det(method="bareis") == -1
    assert M.det(method="berkowitz") == -1

    M = Matrix(( (x,   1),
                 (y, 2*y) ))

    assert M.det(method="bareis") == 2*x*y-y
    assert M.det(method="berkowitz") == 2*x*y-y

    M = Matrix(( (1, 1, 1),
                 (1, 2, 3),
                 (1, 3, 6) ))

    assert M.det(method="bareis") == 1
    assert M.det(method="berkowitz") == 1

    M = Matrix(( ( 3, -2,  0, 5),
                 (-2,  1, -2, 2),
                 ( 0, -2,  5, 0),
                 ( 5,  0,  3, 4) ))

    assert M.det(method="bareis") == -289
    assert M.det(method="berkowitz") == -289

    M = Matrix(( ( 1,  2,  3,  4),
                 ( 5,  6,  7,  8),
                 ( 9, 10, 11, 12),
                 (13, 14, 15, 16) ))

    assert M.det(method="bareis") == 0
    assert M.det(method="berkowitz") == 0

    M = Matrix(( (3, 2, 0, 0, 0),
                 (0, 3, 2, 0, 0),
                 (0, 0, 3, 2, 0),
                 (0, 0, 0, 3, 2),
                 (2, 0, 0, 0, 3) ))

    assert M.det(method="bareis") == 275
    assert M.det(method="berkowitz") == 275

    M = Matrix(( (1, 0,  1,  2, 12),
                 (2, 0,  1,  1,  4),
                 (2, 1,  1, -1,  3),
                 (3, 2, -1,  1,  8),
                 (1, 1,  1,  0,  6) ))

    assert M.det(method="bareis") == -55
    assert M.det(method="berkowitz") == -55

    M = Matrix(( (-5,  2,  3,  4,  5),
                 ( 1, -4,  3,  4,  5),
                 ( 1,  2, -3,  4,  5),
                 ( 1,  2,  3, -2,  5),
                 ( 1,  2,  3,  4, -1) ))

    assert M.det(method="bareis") == 11664
    assert M.det(method="berkowitz") == 11664

    M = Matrix(( ( 2,  7, -1, 3, 2),
                 ( 0,  0,  1, 0, 1),
                 (-2,  0,  7, 0, 2),
                 (-3, -2,  4, 5, 3),
                 ( 1,  0,  0, 0, 1) ))

    assert M.det(method="bareis") == 123
    assert M.det(method="berkowitz") == 123

    M = Matrix(( (x,y,z),
                 (1,0,0),
                 (y,z,x) ))

    assert M.det(method="bareis") == z**2 - x*y
    assert M.det(method="berkowitz") == z**2 - x*y

def test_submatrix():
    m0 = eye(4)
    assert m0[0:3, 0:3] == eye(3)
    assert m0[2:4, 0:2] == zeros(2)

    m1 = Matrix(3,3, lambda i,j: i+j)
    assert m1[0,:] == Matrix(1,3,(0,1,2))
    assert m1[1:3, 1] == Matrix(2,1,(2,3))

    m2 = Matrix([[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]])
    assert m2[:,-1] == Matrix(4,1,[3,7,11,15])
    assert m2[-2:,:] == Matrix([[8,9,10,11],[12,13,14,15]])

def test_submatrix_assignment():
    m = zeros(4)
    m[2:4, 2:4] = eye(2)
    assert m == Matrix(((0,0,0,0),
                        (0,0,0,0),
                        (0,0,1,0),
                        (0,0,0,1)))
    m[0:2, 0:2] = eye(2)
    assert m == eye(4)
    m[:,0] = Matrix(4,1,(1,2,3,4))
    assert m == Matrix(((1,0,0,0),
                        (2,1,0,0),
                        (3,0,1,0),
                        (4,0,0,1)))
    m[:,:] = zeros(4)
    assert m == zeros(4)
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
    assert m1.reshape(4,3) == Matrix(((0,1,2), (3,1,2), (3,4,2), (3,4,5)))
    assert m1.reshape(2,6) == Matrix(((0,1,2,3,1,2), (3,4,2,3,4,5)))

def test_applyfunc():
    m0 = eye(3)
    assert m0.applyfunc(lambda x:2*x) == eye(3)*2
    assert m0.applyfunc(lambda x: 0 ) == zeros(3)

def test_expand():
    x,y = symbols('xy')
    m0 = Matrix([[x*(x+y),2],[((x+y)*y)*x,x*(y+x*(x+y))]])
    # Test if expand() returns a matrix
    m1 = m0.expand()
    assert m1 == Matrix([[x*y+x**2,2],[x*y**2+y*x**2,x*y+y*x**2+x**3]])

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
    assert (L*U).permuteBkwd(p)-testmat == zeros(4)

    testmat = Matrix([[6,-2,7,4],
                      [0,3,6,7],
                      [1,-2,7,4],
                      [-9,2,6,3]])
    L,U,p = testmat.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-testmat == zeros(4)

    x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
    M = Matrix(((1, x, 1), (2, y, 0), (y, 0, z)))
    L, U, p = M.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-M == zeros(3)

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
    assert A.inv("LU") == eye(4)
    assert A.inv("ADJ") == eye(4)
    A = Matrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    Ainv = A.inv()
    assert A*Ainv == eye(3)
    assert A.inv("LU") == Ainv
    assert A.inv("ADJ") == Ainv

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
    assert Q*Q.T == eye(Q.rows)
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

    M = Matrix([[-5,-1, 4,-3,-1],
                [ 1,-1,-1, 1, 0],
                [-1, 0, 0, 0, 0],
                [ 4, 1,-4, 3, 1],
                [-2, 0, 2,-2,-1]])
    assert M*M.nullspace()[0] == Matrix(5,1,[0]*5)

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

def test_wronskian():
    x = Symbol('x')
    assert wronskian([cos(x), sin(x)], x) == cos(x)**2 + sin(x)**2
    assert wronskian([exp(x), exp(2*x)], x) == exp(3*x)
    assert wronskian([exp(x), x], x) == exp(x) - x*exp(x)
    assert wronskian([1, x, x**2], x) == 2
    w1 = -6*exp(x)*sin(x)*x + 6*cos(x)*exp(x)*x**2 - 6*exp(x)*cos(x)*x - \
        exp(x)*cos(x)*x**3 + exp(x)*sin(x)*x**3
    assert wronskian([exp(x), cos(x), x**3], x).expand() == w1
    assert wronskian([exp(x), cos(x), x**3], x, method='berkowitz').expand() == w1
    w2 = -x**3*cos(x)**2 - x**3*sin(x)**2 - 6*x*cos(x)**2 - 6*x*sin(x)**2
    assert wronskian([sin(x), cos(x), x**3], x).expand() == w2
    assert wronskian([sin(x), cos(x), x**3], x, \
        method='berkowitz').expand() == w2


def canonicalize(v):
    """
    Takes the output of eigenvects() and makes it canonical, so that we can
    compare it across platforms.

    It converts Matrices to lists, and uses set() to list the outer list in a
    platform independent way.
    """
    def c(x):
        a, b, c = x
        return (S(a), S(b), tuple(c[0]))
    return tuple(set([c(x) for x in v]))

def test_eigen():
    x,y = symbols('xy')

    R = Rational

    assert eye(3).charpoly(x) == Poly((x-1)**3, x)
    assert eye(3).charpoly(y) == Poly((y-1)**3, y)

    M = Matrix([[1,0,0],
                [0,1,0],
                [0,0,1]])

    assert M.eigenvals() == {S.One: 3}

    assert canonicalize(M.eigenvects()) == canonicalize(
        [(1, 3, [Matrix([1,0,0]),
                 Matrix([0,1,0]),
                 Matrix([0,0,1])])])

    M = Matrix([[0,1,1],
                [1,0,0],
                [1,1,1]])

    assert M.eigenvals() == {2*S.One: 1, -S.One: 1, S.Zero: 1}

    assert canonicalize(M.eigenvects()) == canonicalize(
        [( 2, 1, [Matrix([R(2,3), R(1,3), 1])]),
         (-1, 1, [Matrix([-1, 1, 0])]),
         ( 0, 1, [Matrix([ 0,-1, 1])])])

    eps = Symbol('eps',real=True)

    M = Matrix([[abs(eps), I*eps    ],
               [-I*eps,   abs(eps) ]])

    assert canonicalize(M.eigenvects()) == canonicalize(
        [( 2*abs(eps), 1, [ Matrix([[I*eps/abs(eps)],[1]]) ] ),
         ( 0, 1, [Matrix([[-I*eps/abs(eps)],[1]])]) ])

def test_sparse_matrix():
    return
    def eye(n):
        tmp = SMatrix(n,n,lambda i,j:0)
        for i in range(tmp.rows):
            tmp[i,i] = 1
        return tmp
    def zeros(n):
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
    assert m.cols == m.rows
    assert m.cols == 2
    assert m[:] == [x,0,0,0]
    b = SMatrix(2,2, [x, 0, 0, 0])
    m = b
    assert m.cols == m.rows
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
    assert m0[2:4, 0:2] == zeros(2)

    m1 = SMatrix(3,3, lambda i,j: i+j)
    assert m1[0,:] == SMatrix(1,3,(0,1,2))
    assert m1[1:3, 1] == SMatrix(2,1,(2,3))

    m2 = SMatrix([0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15])
    assert m2[:,-1] == SMatrix(4,1,[3,7,11,15])
    assert m2[-2:,:] == SMatrix([[8,9,10,11],[12,13,14,15]])

    # test_submatrix_assignment
    m = zeros(4)
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
    m[:,:] = zeros(4)
    assert m == zeros(4)
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
    assert m0.applyfunc(lambda x: 0 ) == zeros(3)

    # test_LUdecomp
    testmat = SMatrix([[0,2,5,3],
                      [3,3,7,4],
                      [8,4,0,2],
                      [-2,6,3,4]])
    L,U,p = testmat.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-testmat == zeros(4)

    testmat = SMatrix([[6,-2,7,4],
                      [0,3,6,7],
                      [1,-2,7,4],
                      [-9,2,6,3]])
    L,U,p = testmat.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-testmat == zeros(4)

    x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
    M = Matrix(((1, x, 1), (2, y, 0), (y, 0, z)))
    L, U, p = M.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-M == zeros(3)

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
    assert A.inv("LU") == eye(4)
    assert A.inv("ADJ") == eye(4)
    A = SMatrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    Ainv = A.inv()
    assert A*Ainv == eye(3)
    assert A.inv("LU") == Ainv
    assert A.inv("ADJ") == Ainv

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

    assert M.zeros((3, 5)) == SMatrix(3, 5, {})

def test_subs():
    x = Symbol('x')
    assert Matrix([[1,x],[x,4]]).subs(x, 5) == Matrix([[1,5],[5,4]])
    y = Symbol('y')
    assert Matrix([[x,2],[x+y,4]]).subs([[x,-1],[y,-2]]) == Matrix([[-1,2],[-3,4]])
    assert Matrix([[x,2],[x+y,4]]).subs([(x,-1),(y,-2)]) == Matrix([[-1,2],[-3,4]])
    assert Matrix([[x,2],[x+y,4]]).subs({x:-1,y:-2}) == Matrix([[-1,2],[-3,4]])

def test_simplify():
    x,y,f,n = symbols('xyfn')
    M = Matrix([ [    1/x + 1/y,               (x+x*y)/ x             ],
                 [(f(x)+y*f(x))/f(x), (2 * (1/n - cos(n * pi)/n))/ pi ]
                 ])
    M.simplify()
    assert M ==  Matrix([[(x+y)/(x*y),               1 + y           ],
                         [   1 + y,       (2 - 2*cos(pi*n))/ (pi*n)   ]])

def test_transpose():
    M = Matrix([[1,2,3,4,5,6,7,8,9,0],
                [1,2,3,4,5,6,7,8,9,0]])
    assert M.T == Matrix( [ [1,1],
                            [2,2],
                            [3,3],
                            [4,4],
                            [5,5],
                            [6,6],
                            [7,7],
                            [8,8],
                            [9,9],
                            [0,0] ])
    assert M.T.T == M
    assert M.T == M.transpose()

def test_conjugate():
    M = Matrix([ [0,I,5],
                 [1,2,0]])

    assert M.T == Matrix([ [0,1],
                           [I,2],
                           [5,0]])

    assert M.C == Matrix([ [0,-I,5],
                           [1,2,0]])
    assert M.C == M.conjugate()

    assert M.H == M.T.C
    assert M.H == Matrix([ [0,1],
                           [-I,2],
                           [5,0]])

def test_conj_dirac():
    raises(ShapeError, "eye(3).D")

    M = Matrix([ [1,I,I,I],
                 [0,1,I,I],
                 [0,0,1,I],
                 [0,0,0,1] ])

    assert M.D == Matrix([
                 [1,0,0,0],
                 [-I,1,0,0],
                 [-I,-I,-1,0],
                 [-I,-I,I,-1] ])


def test_trace():
    M = Matrix([[1,0,0],
                [0,5,0],
                [0,0,8]])
    assert M.trace() == 14

def test_shape():
    x, y = symbols("xy")
    M = Matrix([[x,0,0],
                [0,y,0]])
    assert M.shape == (2, 3)

def test_col_row():
    x, y = symbols("xy")
    M = Matrix([[x,0,0],
                [0,y,0]])
    M.row(1,lambda x,i: x+i+1)
    assert M == Matrix([[x,0,0],
                        [1,y+2,3]])
    M.col(0,lambda x,i: x+y**i)
    assert M == Matrix([[x+1,0,0],
                        [1+y,y+2,3]])

def test_issue851():
    m = Matrix([1, 2, 3])
    a = Matrix([1, 2, 3])
    b = Matrix([2, 2, 3])
    assert not (m in [])
    assert not (m in [1])
    assert m != 1
    assert m == a
    assert m != b

def test_issue882():
    class Int1(object):
        def __int__(self):
            return 1
    class Int2(object):
        def __int__(self):
            return 2
    class Index1(object):
        def __index__(self):
            return 1
    class Index2(object):
        def __index__(self):
            return 2
    int1 = Int1()
    int2 = Int2()
    index1 = Index1()
    index2 = Index2()

    m = Matrix([1, 2, 3])
    assert m[int2] == 3
    assert m[index2] == 3

    m[int2] = 4
    assert m[2] == 4
    m[index2] = 5
    assert m[2] == 5

    m = Matrix([[1, 2, 3], [4, 5, 6]])
    assert m[int1,int2] == 6
    assert m[index1,int2] == 6
    assert m[int1,index2] == 6
    assert m[index1,index2] == 6
    assert m[1,int2] == 6
    assert m[1,index2] == 6
    assert m[int1,2] == 6
    assert m[index1,2] == 6

    m[int1,int2] = 1
    assert m[1, 2] == 1
    m[index1,int2] = 2
    assert m[1, 2] == 2
    m[int1,index2] = 3
    assert m[1, 2] == 3
    m[index1,index2] = 4
    assert m[1, 2] == 4
    m[1,int2] = 5
    assert m[1, 2] == 5
    m[1,index2] = 6
    assert m[1, 2] == 6
    m[int1,2] = 7
    assert m[1, 2] == 7
    m[index1,2] = 8
    assert m[1, 2] == 8

def test_evalf():
    a = Matrix([sqrt(5), 6])
    assert abs(a.evalf()[0] - a[0].evalf()) < 1e-10
    assert abs(a.evalf()[1] - a[1].evalf()) < 1e-10

def test_is_symbolic():
    x = Symbol('x')
    a = Matrix([[x,x],[x,x]])
    assert a.is_symbolic() == True
    a = Matrix([[1,2],[3,4]])
    assert a.is_symbolic() == False
    a = Matrix([[1,x],[3,4]])
    assert a.is_symbolic() == True

def test_zeros_ones_fill():
    n, m = 3, 5

    a = zeros( (n, m) )
    a.fill( 5 )

    b = 5 * ones( (n, m) )

    assert a == b
    assert a.rows == b.rows == 3
    assert a.cols == b.cols == 5
    assert a.shape == b.shape == (3, 5)

def test_issue650():
    x, y = symbols('x','y')
    a = Matrix([[x**2, x*y],[x*sin(y), x*cos(y)]])
    assert a.diff(x) == Matrix([[2*x, y],[sin(y), cos(y)]])
    assert Matrix([[x, -x, x**2],[exp(x),1/x-exp(-x), x+1/x]]).limit(x, oo) == Matrix([[oo, -oo, oo],[oo, 0, oo]])
    assert Matrix([[(exp(x)-1)/x, 2*x + y*x, x**x ],
                    [1/x, abs(x) , abs(sin(x+1))]]).limit(x, 0) == Matrix([[1, 0, 1],[oo, 0, sin(1)]])
    assert a.integrate(x) == Matrix([[Rational(1,3)*x**3, y*x**2/2],[x**2*sin(y)/2, x**2*cos(y)/2]])

def test_inv_iszerofunc():
    A = eye(4)
    A.col_swap(0,1)
    for method in "GE", "LU":
        assert A.inv(method, iszerofunc=lambda x: x==0) == A.inv("ADJ")

def test_jacobian_metrics():
    rho, phi = symbols("rho phi")
    X = Matrix([rho*cos(phi), rho*sin(phi)])
    Y = Matrix([rho, phi])
    J = X.jacobian(Y)
    assert J == X.jacobian(Y.T)
    assert J == (X.T).jacobian(Y)
    assert J == (X.T).jacobian(Y.T)
    g = J.T*eye(J.shape[0])*J
    g = g.applyfunc(trigsimp)
    assert g == Matrix([[1, 0], [0, rho**2]])

def test_jacobian2():
    rho, phi = symbols("rho phi")
    X = Matrix([rho*cos(phi), rho*sin(phi), rho**2])
    Y = Matrix([rho, phi])
    J = Matrix([
            [cos(phi), -rho*sin(phi)],
            [sin(phi),  rho*cos(phi)],
            [   2*rho,             0],
        ])
    assert X.jacobian(Y) == J

def test_issue1465():
    x, y, z = symbols('x', 'y', 'z')
    X = Matrix([exp(x + y + z), exp(x + y + z), exp(x + y + z)])
    Y = Matrix([x, y, z])
    for i in range(1, 3):
        for j in range(1, 3):
            X_slice = X[:i,:]
            Y_slice = Y[:j,:]
            J = X_slice.jacobian(Y_slice)
            assert J.rows == i
            assert J.cols == j
            for k in range(j):
                assert J[:,k] == X_slice

def test_nonvectorJacobian():
    x, y, z = symbols('x', 'y', 'z')
    X = Matrix([ [exp(x + y + z), exp(x + y + z)],
                 [exp(x + y + z), exp(x + y + z)] ])
    Y = Matrix([x, y, z])
    raises(TypeError, 'X.jacobian(Y)')
    X = X[0,:]
    Y = Matrix([ [x, y], [x,z] ])
    raises(TypeError, 'X.jacobian(Y)')

def test_vec():
    m = Matrix([ [1,3], [2,4] ])
    m_vec = m.vec()
    assert m_vec.cols == 1
    for i in xrange(4):
        assert m_vec[i] == i + 1

def test_vech():
    m = Matrix([ [1,2], [2,3] ])
    m_vech = m.vech()
    assert m_vech.cols == 1
    for i in xrange(3):
        assert m_vech[i] == i + 1
    m_vech = m.vech(diagonal = False)
    assert m_vech[0] == 2

def test_vech_TypeError():
    m = Matrix([ [1,3] ])
    raises(TypeError, 'm.vech()')
    m = Matrix([ [1,3], [2,4] ])
    raises(TypeError, 'm.vech()')

def test_block_diag1():
    x, y, z = symbols("x y z")
    a = Matrix([[1, 2], [2, 3]])
    b = Matrix([[3, x], [y, 3]])
    c = Matrix([[3, x, 3], [y, 3, z], [x, y, z]])
    assert block_diag([a, b, b]) == Matrix([
            [1, 2, 0, 0, 0, 0],
            [2, 3, 0, 0, 0, 0],
            [0, 0, 3, x, 0, 0],
            [0, 0, y, 3, 0, 0],
            [0, 0, 0, 0, 3, x],
            [0, 0, 0, 0, y, 3],
            ])
    assert block_diag([a, b, c]) == Matrix([
            [1, 2, 0, 0, 0, 0, 0],
            [2, 3, 0, 0, 0, 0, 0],
            [0, 0, 3, x, 0, 0, 0],
            [0, 0, y, 3, 0, 0, 0],
            [0, 0, 0, 0, 3, x, 3],
            [0, 0, 0, 0, y, 3, z],
            [0, 0, 0, 0, x, y, z],
            ])
    assert block_diag([a, c, b]) == Matrix([
            [1, 2, 0, 0, 0, 0, 0],
            [2, 3, 0, 0, 0, 0, 0],
            [0, 0, 3, x, 3, 0, 0],
            [0, 0, y, 3, z, 0, 0],
            [0, 0, x, y, z, 0, 0],
            [0, 0, 0, 0, 0, 3, x],
            [0, 0, 0, 0, 0, y, 3],
            ])

def test_get_diag_blocks1():
    x, y, z = symbols("x y z")
    a = Matrix([[1, 2], [2, 3]])
    b = Matrix([[3, x], [y, 3]])
    c = Matrix([[3, x, 3], [y, 3, z], [x, y, z]])
    a.get_diag_blocks() == [a]
    b.get_diag_blocks() == [b]
    c.get_diag_blocks() == [c]

def test_get_diag_blocks2():
    x, y, z = symbols("x y z")
    a = Matrix([[1, 2], [2, 3]])
    b = Matrix([[3, x], [y, 3]])
    c = Matrix([[3, x, 3], [y, 3, z], [x, y, z]])
    assert block_diag([a, b, b]).get_diag_blocks() == [a, b, b]
    assert block_diag([a, b, c]).get_diag_blocks() == [a, b, c]
    assert block_diag([a, c, b]).get_diag_blocks() == [a, c, b]
    assert block_diag([c, c, b]).get_diag_blocks() == [c, c, b]

def test_inv_block():
    x, y, z = symbols("x y z")
    a = Matrix([[1, 2], [2, 3]])
    b = Matrix([[3, x], [y, 3]])
    c = Matrix([[3, x, 3], [y, 3, z], [x, y, z]])
    A = block_diag([a, b, b])
    assert A.inv(try_block_diag=True) == block_diag([a.inv(), b.inv(), b.inv()])
    A = block_diag([a, b, c])
    assert A.inv(try_block_diag=True) == block_diag([a.inv(), b.inv(), c.inv()])
    A = block_diag([a, c, b])
    assert A.inv(try_block_diag=True) == block_diag([a.inv(), c.inv(), b.inv()])
    A = block_diag([a, a, b, a, c, a])
    assert A.inv(try_block_diag=True) == block_diag([
        a.inv(), a.inv(), b.inv(), a.inv(), c.inv(), a.inv()])
    assert A.inv(try_block_diag=True, method="ADJ") == block_diag([
        a.inv(method="ADJ"), a.inv(method="ADJ"), b.inv(method="ADJ"),
        a.inv(method="ADJ"), c.inv(method="ADJ"), a.inv(method="ADJ")])
