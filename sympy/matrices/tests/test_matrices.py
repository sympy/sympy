from sympy import (symbols, Matrix, SparseMatrix, eye, I, Symbol, Rational,
    Float, wronskian, cos, sin, exp, hessian, sqrt, zeros, ones, randMatrix,
    Poly, S, pi, E, I, oo, trigsimp, Integer, block_diag, N, zeros, sympify,
    Pow, simplify, Min, Max, Abs)
from sympy.matrices.matrices import (ShapeError, MatrixError,
    matrix_multiply_elementwise, diag,

    SparseMatrix, SparseMatrix, NonSquareMatrixError, _dims_to_nm,
    matrix_multiply_elementwise)
from sympy.utilities.pytest import raises
#from sympy.functions.elementary.miscellaneous import Max, Min
#from sympy.functions.elementary.miscellaneous import Max, Min

def test_division():
    x, y, z = symbols('x y z')
    v = Matrix(1,2,[x, y])
    assert v.__div__(z) == Matrix(1,2,[x/z, y/z])
    assert v.__truediv__(z) == Matrix(1,2,[x/z, y/z])
    assert v/z == Matrix(1,2,[x/z, y/z])

def test_sum():
    x, y, z = symbols('x y z')
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

    h = matrix_multiply_elementwise(a, c)
    assert h == a.multiply_elementwise(c)
    assert h[0,0]==7
    assert h[0,1]==4
    assert h[1,0]==18
    assert h[1,1]==6
    assert h[2,0]==0
    assert h[2,1]==0
    raises(ShapeError, 'matrix_multiply_elementwise(a, b)')

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

    A = Matrix([[33, 24], [48, 57]])
    assert (A**(S(1)/2))[:] == [5, 2, 4, 7]
    A = Matrix([[0, 4], [-1, 5]])
    assert (A**(S(1)/2))**2 == A

def test_creation():
    raises(ValueError, 'Matrix(5, 5, range(20))')

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

    c = Matrix((
          Matrix((
            (1, 2, 3),
            (4, 5, 6)
          )),
          (7, 8, 9)
    ))
    assert c.cols == 3
    assert c.rows == 3
    assert c[:] == [1,2,3,4,5,6,7,8,9]

def test_tolist():
    x, y, z = symbols('x y z')
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

def test_extract():
    m = Matrix(4, 3, lambda i, j: i*3 + j)
    assert m.extract([0,1,3],[0,1]) == Matrix(3,2,[0,1,3,4,9,10])
    assert m.extract([0,3],[0,0,2]) == Matrix(2,3,[0,0,2,9,9,11])
    assert m.extract(range(4),range(3)) == m
    raises(IndexError, 'm.extract([4], [0])')
    raises(IndexError, 'm.extract([0], [3])')

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
    x,y = symbols('x y')
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

    mL = Matrix((
      (1,0,0),
      (2,3,0),
    ))
    assert mL.is_lower() == True
    assert mL.is_upper() == False
    mU = Matrix((
      (1,2,3),
      (0,4,5),
    ))
    assert mU.is_lower() == False
    assert mU.is_upper() == True

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

    M = Matrix([[0, 0, 1],
                 [2,3,0],
                 [3, 1, 4]])
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

def test_QRsolve():
    A = Matrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    x = Matrix(3,1,[3,7,5])
    b = A*x
    soln = A.QRsolve(b)
    assert soln == x
    x = Matrix([[1,2],[3,4],[5,6]])
    b = A*x
    soln = A.QRsolve(b)
    assert soln == x

    A = Matrix([[0,-1,2],
                [5,10,7],
                [8,3,4]])
    x = Matrix(3,1,[-1,2,5])
    b = A*x
    soln = A.QRsolve(b)
    assert soln == x
    x = Matrix([[7,8],[9,10],[11,12]])
    b = A*x
    soln = A.QRsolve(b)
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
    assert Q.T * Q == eye(Q.cols)
    assert R.is_upper()
    assert A == Q*R

def test_QR_non_square():
    A = Matrix([[9,0,26],[12,0,-7],[0,4,4],[0,-3,-3]])
    Q, R = A.QRdecomposition()
    assert Q.T * Q == eye(Q.cols)
    assert R.is_upper()
    assert A == Q*R

    A = Matrix([[1,-1,4],[1,4,-2],[1,4,2],[1,-1,0]])
    Q, R = A.QRdecomposition()
    assert Q.T * Q == eye(Q.cols)
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

    # issue 1698; just see that we can do it when rows > cols
    M = Matrix([[1,2],[2,4],[3,6]])
    assert M.nullspace()

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
    x,y = symbols('x y')

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

    M = Matrix([ [1, -1],
                 [1,  3]])
    assert canonicalize(M.eigenvects()) == canonicalize(
        [[2, 2, [Matrix(1,2,[-1,1])]]])

    M = Matrix([ [1, 2, 3], [4, 5, 6], [7, 8, 9] ])
    a=R(15,2)
    b=3*33**R(1,2)
    c=R(13,2)
    d=(R(33,8) + 3*b/8)
    e=(R(33,8) - 3*b/8)
    def NS(e, n):
        return str(N(e, n))
    r = [
        (a - b/2, 1, [Matrix([(12 + 24/(c - b/2))/((c - b/2)*e) + 3/(c - b/2),
                              (6 + 12/(c - b/2))/e,1])]),
        (      0, 1, [Matrix([1,-2,1])]),
        (a + b/2, 1, [Matrix([(12 + 24/(c + b/2))/((c + b/2)*d) + 3/(c + b/2),
                              (6 + 12/(c + b/2))/d,1])]),
        ]
    r1 = [(NS(r[i][0],2),NS(r[i][1],2),[NS(j,2) for j in r[i][2][0]]) for i in range(len(r))]
    r = M.eigenvects()
    r2=[(NS(r[i][0],2),NS(r[i][1],2),[NS(j,2) for j in r[i][2][0]]) for i in range(len(r))]
    assert sorted(r1) == sorted(r2)

    eps = Symbol('eps',real=True)

    M = Matrix([[abs(eps), I*eps    ],
               [-I*eps,   abs(eps) ]])

    assert canonicalize(M.eigenvects()) == canonicalize(
        [( 2*abs(eps), 1, [ Matrix([[I*eps/abs(eps)],[1]]) ] ),
         ( 0, 1, [Matrix([[-I*eps/abs(eps)],[1]])]) ])

def test_sparse_matrix():
    return
    def eye(n):
        tmp = SparseMatrix(n,n,lambda i,j:0)
        for i in range(tmp.rows):
            tmp[i,i] = 1
        return tmp
    def zeros(n):
        return SparseMatrix(n,n,lambda i,j:0)

    # test element assignment
    a = SparseMatrix((
        (1, 0),
        (0, 1)
    ))
    a[0, 0] = 2
    assert a == SparseMatrix((
        (2, 0),
        (0, 1)
    ))
    a[1, 0] = 5
    assert a == SparseMatrix((
        (2, 0),
        (5, 1)
    ))
    a[1, 1] = 0
    assert a == SparseMatrix((
        (2, 0),
        (5, 0)
    ))
    assert a.mat == {(0, 0): 2, (1, 0): 5}

    # test_multiplication
    a=SparseMatrix((
        (1, 2),
        (3, 1),
        (0, 6),
        ))

    b = SparseMatrix ((
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
    assert isinstance(c,SparseMatrix)
    assert c[0,0] == x
    assert c[0,1] == 2*x
    assert c[1,0] == 3*x
    assert c[1,1] == 0

    c = 5 * b
    assert isinstance(c,SparseMatrix)
    assert c[0,0] == 5
    assert c[0,1] == 2*5
    assert c[1,0] == 3*5
    assert c[1,1] == 0

    #test_power
    A = SparseMatrix([[2,3],[4,5]])
    assert (A**5)[:] == [6140, 8097, 10796, 14237]
    A = SparseMatrix([[2, 1, 3],[4,2, 4], [6,12, 1]])
    assert (A**3)[:] == [290, 262, 251, 448, 440, 368, 702, 954, 433]


    # test_creation
    x = Symbol("x")
    a = SparseMatrix([x, 0], [0, 0])
    m = a
    assert m.cols == m.rows
    assert m.cols == 2
    assert m[:] == [x,0,0,0]
    b = SparseMatrix(2,2, [x, 0, 0, 0])
    m = b
    assert m.cols == m.rows
    assert m.cols == 2
    assert m[:] == [x,0,0,0]

    assert a == b

    # test_determinant
    x, y = Symbol('x'), Symbol('y')

    assert SparseMatrix([ [1] ]).det() == 1

    assert SparseMatrix(( (-3,  2),
                    ( 8, -5) )).det() == -1

    assert SparseMatrix(( (x,   1),
                    (y, 2*y) )).det() == 2*x*y-y

    assert SparseMatrix(( (1, 1, 1),
                    (1, 2, 3),
                    (1, 3, 6) )).det() == 1

    assert SparseMatrix(( ( 3, -2,  0, 5),
                    (-2,  1, -2, 2),
                    ( 0, -2,  5, 0),
                    ( 5,  0,  3, 4) )).det() == -289

    assert SparseMatrix(( ( 1,  2,  3,  4),
                    ( 5,  6,  7,  8),
                    ( 9, 10, 11, 12),
                    (13, 14, 15, 16) )).det() == 0

    assert SparseMatrix(( (3, 2, 0, 0, 0),
                    (0, 3, 2, 0, 0),
                    (0, 0, 3, 2, 0),
                    (0, 0, 0, 3, 2),
                    (2, 0, 0, 0, 3) )).det() == 275

    assert SparseMatrix(( (1, 0,  1,  2, 12),
                    (2, 0,  1,  1,  4),
                    (2, 1,  1, -1,  3),
                    (3, 2, -1,  1,  8),
                    (1, 1,  1,  0,  6) )).det() == -55

    assert SparseMatrix(( (-5,  2,  3,  4,  5),
                    ( 1, -4,  3,  4,  5),
                    ( 1,  2, -3,  4,  5),
                    ( 1,  2,  3, -2,  5),
                    ( 1,  2,  3,  4, -1) )).det() == 11664

    assert SparseMatrix(( ( 2,  7, -1, 3, 2),
                    ( 0,  0,  1, 0, 1),
                    (-2,  0,  7, 0, 2),
                    (-3, -2,  4, 5, 3),
                    ( 1,  0,  0, 0, 1) )).det() == 123

    # test_submatrix
    m0 = eye(4)
    assert m0[0:3, 0:3] == eye(3)
    assert m0[2:4, 0:2] == zeros(2)

    m1 = SparseMatrix(3,3, lambda i,j: i+j)
    assert m1[0,:] == SparseMatrix(1,3,(0,1,2))
    assert m1[1:3, 1] == SparseMatrix(2,1,(2,3))

    m2 = SparseMatrix([0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15])
    assert m2[:,-1] == SparseMatrix(4,1,[3,7,11,15])
    assert m2[-2:,:] == SparseMatrix([[8,9,10,11],[12,13,14,15]])

    # test_submatrix_assignment
    m = zeros(4)
    m[2:4, 2:4] = eye(2)
    assert m == SparseMatrix((0,0,0,0),
                        (0,0,0,0),
                        (0,0,1,0),
                        (0,0,0,1))
    m[0:2, 0:2] = eye(2)
    assert m == eye(4)
    m[:,0] = SparseMatrix(4,1,(1,2,3,4))
    assert m == SparseMatrix((1,0,0,0),
                        (2,1,0,0),
                        (3,0,1,0),
                        (4,0,0,1))
    m[:,:] = zeros(4)
    assert m == zeros(4)
    m[:,:] = ((1,2,3,4),(5,6,7,8),(9, 10, 11, 12),(13,14,15,16))
    assert m == SparseMatrix(((1,2,3,4),
                        (5,6,7,8),
                        (9, 10, 11, 12),
                        (13,14,15,16)))
    m[0:2, 0] = [0,0]
    assert m == SparseMatrix(((0,2,3,4),
                        (0,6,7,8),
                        (9, 10, 11, 12),
                        (13,14,15,16)))

    # test_reshape
    m0 = eye(3)
    assert m0.reshape(1,9) == SparseMatrix(1,9,(1,0,0,0,1,0,0,0,1))
    m1 = SparseMatrix(3,4, lambda i,j: i+j)
    assert m1.reshape(4,3) == SparseMatrix((0,1,2), (3,1,2), (3,4,2), (3,4,5))
    assert m1.reshape(2,6) == SparseMatrix((0,1,2,3,1,2), (3,4,2,3,4,5))

    # test_applyfunc
    m0 = eye(3)
    assert m0.applyfunc(lambda x:2*x) == eye(3)*2
    assert m0.applyfunc(lambda x: 0 ) == zeros(3)

    # test_LUdecomp
    testmat = SparseMatrix([[0,2,5,3],
                      [3,3,7,4],
                      [8,4,0,2],
                      [-2,6,3,4]])
    L,U,p = testmat.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-testmat == zeros(4)

    testmat = SparseMatrix([[6,-2,7,4],
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
    A = SparseMatrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    x = SparseMatrix(3,1,[3,7,5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x
    A = SparseMatrix([[0,-1,2],
                [5,10,7],
                [8,3,4]])
    x = SparseMatrix(3,1,[-1,2,5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x

    # test_inverse
    A = eye(4)
    assert A.inv() == eye(4)
    assert A.inv("LU") == eye(4)
    assert A.inv("ADJ") == eye(4)
    A = SparseMatrix([[2,3,5],
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
    test = SparseMatrix([[1,3,2],[2,6,3],[2,3,6]])
    assert test.cofactorMatrix() == SparseMatrix([[27,-6,-6],[-12,2,3],[-3,1,0]])
    test = SparseMatrix([[1,2,3],[4,5,6],[7,8,9]])
    assert test.cofactorMatrix() == SparseMatrix([[-3,6,-3],[6,-12,6],[-3,6,-3]])

    # test_jacobian
    x = Symbol('x')
    y = Symbol('y')
    L = SparseMatrix(1,2,[x**2*y, 2*y**2 + x*y])
    syms = [x,y]
    assert L.jacobian(syms) == Matrix([[2*x*y, x**2],[y, 4*y+x]])

    L = SparseMatrix(1,2,[x, x**2*y**3])
    assert L.jacobian(syms) == SparseMatrix([[1, 0], [2*x*y**3, x**2*3*y**2]])

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

    assert M.zeros((3, 5)) == SparseMatrix(3, 5, {})

def test_subs():
    x = Symbol('x')
    assert Matrix([[1,x],[x,4]]).subs(x, 5) == Matrix([[1,5],[5,4]])
    y = Symbol('y')
    assert Matrix([[x,2],[x+y,4]]).subs([[x,-1],[y,-2]]) == Matrix([[-1,2],[-3,4]])
    assert Matrix([[x,2],[x+y,4]]).subs([(x,-1),(y,-2)]) == Matrix([[-1,2],[-3,4]])
    assert Matrix([[x,2],[x+y,4]]).subs({x:-1,y:-2}) == Matrix([[-1,2],[-3,4]])

def test_simplify():
    x,y,f,n = symbols('x y f n')
    M = Matrix([ [    1/x + 1/y,               (x + x*y)/ x           ],
                 [(f(x) + y*f(x))/f(x), 2 * (1/n - cos(n * pi)/n)/ pi ]
                 ])
    M.simplify()
    assert M ==  Matrix([[(x + y)/(x * y),                 1 + y       ],
                         [   1 + y,       2*((1 - 1*cos(pi*n))/(pi*n)) ]])
    M = Matrix([[(1 + x)**2]])
    M.simplify()
    assert M == Matrix([[(1 + x)**2]])
    M.simplify(ratio=oo)
    assert M == Matrix([[1 + 2*x + x**2]])

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
    x, y = symbols("x y")
    M = Matrix([[x,0,0],
                [0,y,0]])
    assert M.shape == (2, 3)

def test_col_row():
    x, y = symbols("x y")
    M = Matrix([[x,0,0],
                [0,y,0]])
    M.row(1,lambda r, j: r+j+1)
    assert M == Matrix([[x,0,0],
                        [1,y+2,3]])
    M.col(0,lambda c, j: c+y**j)
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
    class Index1(object):
        def __index__(self):
            return 1
    class Index2(object):
        def __index__(self):
            return 2
    index1 = Index1()
    index2 = Index2()

    m = Matrix([1, 2, 3])

    assert m[index2] == 3

    m[index2] = 5
    assert m[2] == 5

    m = Matrix([[1, 2, 3], [4, 5, 6]])
    assert m[index1,index2] == 6
    assert m[1,index2] == 6
    assert m[index1,2] == 6

    m[index1,index2] = 4
    assert m[1, 2] == 4
    m[1,index2] = 6
    assert m[1, 2] == 6
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
    a = Matrix([[1,2,3,4],[5,6,7,8]])
    assert a.is_symbolic() == False
    a = Matrix([[1,2,3,4],[5,6,x,8]])
    assert a.is_symbolic() == True
    a = Matrix([[1,x,3]])
    assert a.is_symbolic() == True
    a = Matrix([[1,2,3]])
    assert a.is_symbolic() == False
    a = Matrix([[1],[x],[3]])
    assert a.is_symbolic() == True
    a = Matrix([[1],[2],[3]])
    assert a.is_symbolic() == False

def test_is_upper():
    a = Matrix([[1,2,3]])
    assert a.is_upper() == True
    a = Matrix([[1],[2],[3]])
    assert a.is_upper() == False

def test_is_lower():
    a = Matrix([[1,2,3]])
    assert a.is_lower() == False
    a = Matrix([[1],[2],[3]])
    assert a.is_lower() == True

def test_is_nilpotent():
    a = Matrix(4, 4, [0,2,1,6,0,0,1,2,0,0,0,3,0,0,0,0])
    assert a.is_nilpotent()
    a = Matrix([[1,0],[0,1]])
    assert not a.is_nilpotent()

def test_zeros_ones_fill():
    n, m = 3, 5

    a = zeros( (n, m) )
    a.fill( 5 )

    b = 5 * ones( (n, m) )

    assert a == b
    assert a.rows == b.rows == 3
    assert a.cols == b.cols == 5
    assert a.shape == b.shape == (3, 5)

def test_empty_zeros():
    a = zeros(0)
    assert a == Matrix()
    a = zeros([0, 2])
    assert a.rows == 0
    assert a.cols == 2
    a = zeros([2, 0])
    assert a.rows == 2
    assert a.cols == 0

def test_issue650():
    x, y = symbols('x y')
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
    x, y, z = symbols('x y z')
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
    x, y, z = symbols('x y z')
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
    m_vech = m.vech(diagonal=False)
    assert m_vech[0] == 2
    x,y = symbols('x,y')
    m = Matrix([ [1, x*(x+y)], [y*x+x**2, 1] ])
    m_vech = m.vech(diagonal=False)
    assert m_vech[0] == x*(x + y)
    x,y = symbols('x,y')
    m = Matrix([ [1, x*(x+y)], [y*x, 1] ])
    m_vech = m.vech(diagonal=False, check_symmetry=False)
    assert m_vech[0] == y*x

def test_vech_errors():
    m = Matrix([ [1,3] ])
    raises(ShapeError, 'm.vech()')
    m = Matrix([ [1,3], [2,4] ])
    raises(ValueError, 'm.vech()')

def test_diag():
    x, y, z = symbols("x y z")
    a = Matrix([[1, 2], [2, 3]])
    b = Matrix([[3, x], [y, 3]])
    c = Matrix([[3, x, 3], [y, 3, z], [x, y, z]])
    assert diag(a, b, b) == Matrix([
            [1, 2, 0, 0, 0, 0],
            [2, 3, 0, 0, 0, 0],
            [0, 0, 3, x, 0, 0],
            [0, 0, y, 3, 0, 0],
            [0, 0, 0, 0, 3, x],
            [0, 0, 0, 0, y, 3],
            ])
    assert diag(a, b, c) == Matrix([
            [1, 2, 0, 0, 0, 0, 0],
            [2, 3, 0, 0, 0, 0, 0],
            [0, 0, 3, x, 0, 0, 0],
            [0, 0, y, 3, 0, 0, 0],
            [0, 0, 0, 0, 3, x, 3],
            [0, 0, 0, 0, y, 3, z],
            [0, 0, 0, 0, x, y, z],
            ])
    assert diag(a, c, b) == Matrix([
            [1, 2, 0, 0, 0, 0, 0],
            [2, 3, 0, 0, 0, 0, 0],
            [0, 0, 3, x, 3, 0, 0],
            [0, 0, y, 3, z, 0, 0],
            [0, 0, x, y, z, 0, 0],
            [0, 0, 0, 0, 0, 3, x],
            [0, 0, 0, 0, 0, y, 3],
            ])
    a = Matrix([x, y, z])
    b = Matrix([[1, 2], [3, 4]])
    c = Matrix([[5, 6]])
    assert diag(a, 7, b, c) == Matrix([
            [x, 0, 0, 0, 0, 0],
            [y, 0, 0, 0, 0, 0],
            [z, 0, 0, 0, 0, 0],
            [0, 7, 0, 0, 0, 0],
            [0, 0, 1, 2, 0, 0],
            [0, 0, 3, 4, 0, 0],
            [0, 0, 0, 0, 5, 6],
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
    assert diag(a, b, b).get_diag_blocks() == [a, b, b]
    assert diag(a, b, c).get_diag_blocks() == [a, b, c]
    assert diag(a, c, b).get_diag_blocks() == [a, c, b]
    assert diag(c, c, b).get_diag_blocks() == [c, c, b]

def test_inv_block():
    x, y, z = symbols("x y z")
    a = Matrix([[1, 2], [2, 3]])
    b = Matrix([[3, x], [y, 3]])
    c = Matrix([[3, x, 3], [y, 3, z], [x, y, z]])
    A = diag(a, b, b)
    assert A.inv(try_block_diag=True) == diag(a.inv(), b.inv(), b.inv())
    A = diag(a, b, c)
    assert A.inv(try_block_diag=True) == diag(a.inv(), b.inv(), c.inv())
    A = diag(a, c, b)
    assert A.inv(try_block_diag=True) == diag(a.inv(), c.inv(), b.inv())
    A = diag(a, a, b, a, c, a)
    assert A.inv(try_block_diag=True) == diag(
        a.inv(), a.inv(), b.inv(), a.inv(), c.inv(), a.inv())
    assert A.inv(try_block_diag=True, method="ADJ") == diag(
        a.inv(method="ADJ"), a.inv(method="ADJ"), b.inv(method="ADJ"),
        a.inv(method="ADJ"), c.inv(method="ADJ"), a.inv(method="ADJ"))

def test_creation_args():
    """
    Check that matrix dimensions can be specified using any reasonable type
    (see issue 1515).
    """
    raises(ValueError, 'zeros((3, -1))')
    raises(ValueError, 'zeros((1, 2, 3, 4))')
    assert zeros(3L) == zeros(3)
    assert zeros(Integer(3)) == zeros(3)
    assert zeros(3.) == zeros(3)
    assert eye(3L) == eye(3)
    assert eye(Integer(3)) == eye(3)
    assert eye(3.) == eye(3)
    assert ones((3L, Integer(4))) == ones((3, 4))
    raises(TypeError, 'Matrix(1, 2)')

def test_diagonal_symmetrical():
    m = Matrix(2,2,[0, 1, 1, 0])
    assert not m.is_diagonal()
    assert m.is_symmetric()
    assert m.is_symmetric(simplify=False)

    m = Matrix(2,2,[1, 0, 0, 1])
    assert m.is_diagonal()

    m = diag(1, 2, 3)
    assert m.is_diagonal()
    assert m.is_symmetric()

    m = Matrix(3,3,[1, 0, 0, 0, 2, 0, 0, 0, 3])
    assert m == diag(1, 2, 3)

    m = Matrix(2, 3, [0, 0, 0, 0, 0, 0])
    assert not m.is_symmetric()
    assert m.is_diagonal()

    m = Matrix(((5, 0), (0, 6), (0, 0)))
    assert m.is_diagonal()

    m = Matrix(((5, 0, 0), (0, 6, 0)))
    assert m.is_diagonal()

    x, y = symbols('x y')
    m = Matrix(3,3,[1, x**2 + 2*x + 1, y, (x + 1)**2 , 2, 0, y, 0, 3])
    assert m.is_symmetric()
    assert not m.is_symmetric(simplify=False)
    assert m.expand().is_symmetric(simplify=False)


def test_diagonalization():
    x, y, z = symbols('x y z')
    m = Matrix(3,2,[-3, 1, -3, 20, 3, 10])
    assert not m.is_diagonalizable()
    assert not m.is_symmetric()
    raises(NonSquareMatrixError, 'm.diagonalize()')

    # diagonalizable
    m = diag(1, 2, 3)
    (P, D) = m.diagonalize()
    assert P == eye(3)
    assert D == m

    m = Matrix(2,2,[0, 1, 1, 0])
    assert m.is_symmetric()
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D

    m = Matrix(2,2,[1, 0, 0, 3])
    assert m.is_symmetric()
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D
    assert P == eye(2)
    assert D == m

    m = Matrix(2,2,[1, 1, 0, 0])
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D

    m = Matrix(3,3,[1, 2, 0, 0, 3, 0, 2, -4, 2])
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D

    m = Matrix(2,2,[1, 0, 0, 0])
    assert m.is_diagonal()
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D
    assert P == eye(2)

    # diagonalizable, complex only
    m = Matrix(2,2,[0, 1, -1, 0])
    assert not m.is_diagonalizable(True)
    raises(MatrixError, '(D, P) = m.diagonalize(True)')
    assert m.is_diagonalizable()
    (P, D) = m.diagonalize()
    assert P.inv() * m * P == D

    m = Matrix(2,2,[1, 0, 0, I])
    raises(NotImplementedError, 'm.is_diagonalizable(True)')
    # !!! bug because of eigenvects() or roots(x**2 + (-1 - I)*x + I, x)
    # see issue 2193
    # assert not m.is_diagonalizable(True)
    # raises(MatrixError, '(P, D) = m.diagonalize(True)')
    # (P, D) = m.diagonalize(True)

    # not diagonalizable
    m = Matrix(2,2,[0, 1, 0, 0])
    assert not m.is_diagonalizable()
    raises(MatrixError, '(D, P) = m.diagonalize()')

    m = Matrix(3,3,[-3, 1, -3, 20, 3, 10, 2, -2, 4])
    assert not m.is_diagonalizable()
    raises(MatrixError, '(D, P) = m.diagonalize()')

    # symbolic
    a, b, c, d = symbols('a b c d')
    m = Matrix(2,2,[a, c, c, b])
    assert m.is_symmetric()
    assert m.is_diagonalizable()

def test_jordan_form():

    m = Matrix(3,2,[-3, 1, -3, 20, 3, 10])
    raises(NonSquareMatrixError, 'm.jordan_form()')

    # diagonalizable
    m = Matrix(3, 3, [7, -12, 6, 10, -19, 10, 12, -24, 13])
    Jmust = Matrix(3, 3, [1, 0, 0, 0, 1, 0, 0, 0, -1])
    (P, J) = m.jordan_form()
    assert Jmust == J
    assert Jmust == m.diagonalize()[1]

    #m = Matrix(3, 3, [0, 6, 3, 1, 3, 1, -2, 2, 1])
    #m.jordan_form() # very long
    # m.jordan_form() #

    # diagonalizable, complex only

    # Jordan cells
    # complexity: one of eigenvalues is zero
    m = Matrix(3, 3, [0, 1, 0, -4, 4, 0, -2, 1, 2])
    Jmust = Matrix(3, 3, [2, 0, 0, 0, 2, 1, 0, 0, 2])
    assert Jmust == m.jordan_form()[1]
    (P, Jcells) = m.jordan_cells()
    assert Jcells[0] == Matrix(1, 1, [2])
    assert Jcells[1] == Matrix(2, 2, [2, 1, 0, 2])

    #complexity: all of eigenvalues are equal
    m = Matrix(3, 3, [2, 6, -15, 1, 1, -5, 1, 2, -6])
    Jmust = Matrix(3, 3, [-1, 0, 0, 0, -1, 1, 0, 0, -1])
    (P, J) = m.jordan_form()
    assert Jmust == J

    #complexity: two of eigenvalues are zero
    m = Matrix(3, 3, [4, -5, 2, 5, -7, 3, 6, -9, 4])
    Jmust = Matrix(3, 3, [1, 0, 0, 0, 0, 1, 0, 0, 0])
    (P, J) = m.jordan_form()
    assert Jmust == J

    m = Matrix(4, 4, [6, 5, -2, -3, -3, -1, 3, 3, 2, 1, -2, -3, -1, 1, 5, 5])
    Jmust = Matrix(4, 4, [2, 1, 0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 0, 0, 0, 2])
    (P, J) = m.jordan_form()
    assert Jmust == J

    m = Matrix(4, 4, [6, 2, -8, -6, -3, 2, 9, 6, 2, -2, -8, -6, -1, 0, 3, 4])
    Jmust = Matrix(4, 4, [2, 0, 0, 0, 0, 2, 1, 0, 0, 0, 2, 0, 0, 0, 0, -2])
    (P, J) = m.jordan_form()
    assert Jmust == J

    m = Matrix(4, 4, [5, 4, 2, 1, 0, 1, -1, -1, -1, -1, 3, 0, 1, 1, -1, 2])
    assert not m.is_diagonalizable()
    Jmust = Matrix(4, 4, [1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 4, 1, 0, 0, 0, 4])
    (P, J) = m.jordan_form()
    assert Jmust == J

def test_Matrix_berkowitz_charpoly():
    x, UA, K_i, K_w = symbols('x UA K_i K_w')

    A = Matrix([[-K_i - UA + K_i**2/(K_i + K_w),       K_i*K_w/(K_i + K_w)],
                [           K_i*K_w/(K_i + K_w), -K_w + K_w**2/(K_i + K_w)]])

    assert A.berkowitz_charpoly(x) == \
        Poly(x**2 + (K_i*UA + K_w*UA + 2*K_i*K_w)/(K_i + K_w)*x + K_i*K_w*UA/(K_i + K_w), x, domain='ZZ(K_i,K_w,UA)')

def test_exp():
    m = Matrix([[3,4],[0,-2]])
    assert m.exp() == Matrix([[exp(3),-4*exp(-2)/5 + 4*exp(3)/5],[0,exp(-2)]])

    m = Matrix([[1,0],[0,1]])
    assert m.exp() == Matrix([[E,0],[0,E]])

def test_SparseMatrix_transpose():
    assert SparseMatrix((1,2),(3,4)).transpose() == SparseMatrix((1,3),(2,4))

def test_SparseMatrix_CL_RL():
    assert SparseMatrix((1,2),(3,4)).row_list() == [(0, 0, 1), (0, 1, 2), (1, 0, 3), (1, 1, 4)]
    assert SparseMatrix((1,2),(3,4)).col_list() == [(0, 0, 1), (1, 0, 3), (0, 1, 2), (1, 1, 4)]

def test_SparseMatrix_add():
    assert SparseMatrix(((1,0), (0,1))) + SparseMatrix(((0,1), (1,0))) == SparseMatrix(((1,1), (1,1)))
    a = SparseMatrix(100, 100, lambda i, j : int(j != 0 and i % j == 0))
    b = SparseMatrix(100, 100, lambda i, j : int(i != 0 and j % i == 0))
    assert (len(a.mat) + len(b.mat) - len((a+b).mat) > 0)

def test_has():
    x, y, z = symbols('x,y,z')
    A = Matrix(((x,y),(2,3)))
    assert A.has(x)
    assert not A.has(z)
    assert A.has(Symbol)

    A = A.subs(x,2)
    assert not A.has(x)

def test_errors():
    # Note, some errors not tested.  See 'XXX' in code.
    raises(ValueError, "_dims_to_nm([1, 0, 2])")
    raises(ValueError, "Matrix([[1, 2], [1]])")
    raises(ShapeError, "Matrix([[1, 2], [3, 4]]).copyin_matrix([1, 0], Matrix([1, 2]))")
    raises(TypeError, "Matrix([[1, 2], [3, 4]]).copyin_list([0, 1], set([]))")
    raises(NonSquareMatrixError, "Matrix([[1, 2, 3], [2, 3, 0]]).inv()")
    raises(ShapeError, "Matrix(1, 2, [1, 2]).row_join(Matrix([[1, 2], [3, 4]]))")
    raises(ShapeError, "Matrix([1, 2]).col_join(Matrix([[1, 2], [3, 4]]))")
    raises(ShapeError, "Matrix([1]).row_insert(1, Matrix([[1, 2], [3, 4]]))")
    raises(ShapeError, "Matrix([1]).col_insert(1, Matrix([[1, 2], [3, 4]]))")
    raises(NonSquareMatrixError, "Matrix([1, 2]).trace()")
    raises(TypeError, "SparseMatrix([[1, 2], [3, 4]]).submatrix([1, 1])")
    raises(TypeError, "Matrix([1]).applyfunc(1)")
    raises(ShapeError, "Matrix([1]).LUsolve(Matrix([[1, 2], [3, 4]]))")
    raises(MatrixError, "Matrix([[1,2,3],[4,5,6],[7,8,9]]).QRdecomposition()")
    raises(NonSquareMatrixError, "Matrix([1, 2]).LUdecomposition_Simple()")
    raises(ValueError, "Matrix([[1, 2], [3, 4]]).minorEntry(4, 5)")
    raises(ValueError, "Matrix([[1, 2], [3, 4]]).minorMatrix(4, 5)")
    raises(TypeError, "Matrix([1, 2, 3]).cross(1)")
    raises(TypeError, "Matrix([1, 2, 3]).dot(1)")
    raises(ShapeError, "Matrix([1, 2, 3]).dot(Matrix([1, 2]))")
    raises(NotImplementedError, "Matrix([[0,1,2],[0,0,-1], [0,0,0]]).exp()")
    raises(NonSquareMatrixError, "Matrix([1, 2, 3]).exp()")
    raises(ShapeError, "Matrix([[1, 2], [3, 4]]).normalized()")
    raises(NonSquareMatrixError, "Matrix([1, 2]).inverse_GE()")
    raises(ValueError, "Matrix([[1, 2], [1, 2]]).inverse_GE()")
    raises(NonSquareMatrixError, "Matrix([1, 2]).inverse_ADJ()")
    raises(ValueError, "Matrix([[1, 2], [1, 2]]).inverse_ADJ()")
    raises(NonSquareMatrixError, "Matrix([1,2]).is_nilpotent()")
    raises(ValueError, "hessian(Matrix([[1, 2], [3, 4]]), Matrix([[1, 2], [2, 1]]))")
    raises(ValueError, "hessian(Matrix([[1, 2], [3, 4]]), [])")
    raises(TypeError, "SparseMatrix(1.4, 2, lambda i, j: 0)")
    raises(ValueError, "SparseMatrix([1, 2, 3], [1, 2])")
    raises(ValueError, "SparseMatrix([[1, 2], [3, 4]])[(1, 2, 3)]")
    raises(ValueError, "SparseMatrix([[1, 2], [3, 4]]).rowdecomp(5)")
    raises(ValueError, "SparseMatrix([[1, 2], [3, 4]])[1, 2, 3] = 4")
    raises(TypeError, "SparseMatrix([[1, 2], [3, 4]]).copyin_list([0, 1], set([]))")
    raises(TypeError, "SparseMatrix([[1, 2], [3, 4]]).submatrix((1, 2))")
    raises(TypeError, "SparseMatrix([1, 2, 3]).cross(1)")
    raises(ValueError, "Matrix([[5, 10, 7],[0, -1, 2],[8,  3, 4]]).LUdecomposition_Simple(iszerofunc=lambda x:abs(x)<=4)")
    raises(NotImplementedError, "Matrix([[1, 0],[1, 1]])**(S(1)/2)")
    raises(NotImplementedError, "Matrix([[1, 2, 3],[4, 5, 6],[7,  8, 9]])**(0.5)")

def test_len():
    assert len(Matrix()) == 0
    assert len(Matrix([[1, 2]])) == len(Matrix([[1], [2]])) == 2
    assert len(Matrix(0, 2, lambda i, j: 0)) == len(Matrix(2, 0, lambda i, j: 0)) == 0
    assert len(Matrix([[0, 1, 2], [3, 4, 5]])) == 6
    assert Matrix([1])
    assert not Matrix()

def test_integrate():
    x, y = symbols('x,y')
    A = Matrix(((1,4,x),(y,2,4),(10,5,x**2)))
    assert A.integrate(x) == Matrix(((x, 4*x, x**2/2), (x*y, 2*x, 4*x), (10*x, 5*x, x**3/3)))
    assert A.integrate(y) == Matrix(((y, 4*y, x*y),(y**2/2, 2*y, 4*y), (10*y, 5*y, y*x**2)))

def test_limit():
    x, y = symbols('x,y')
    A = Matrix(((1,4,sin(x)/x),(y,2,4),(10,5,x**2+1)))
    assert A.limit(x,0) == Matrix(((1,4,1),(y,2,4),(10,5,1)))

def test_diff():
    x, y = symbols('x,y')
    A = Matrix(((1,4,x),(y,2,4),(10,5,x**2+1)))
    assert A.diff(x) == Matrix(((0,0,1),(0,0,0),(0,0,2*x)))
    assert A.diff(y) == Matrix(((0,0,0),(1,0,0),(0,0,0)))

def test_getattr():
    x, y = symbols('x,y')
    A = Matrix(((1,4,x),(y,2,4),(10,5,x**2+1)))
    raises (AttributeError, 'A.nonexistantattribute')

def test_hessenberg():
    A = Matrix([[3, 4, 1],[2, 4 ,5],[0, 1, 2]])
    assert A.is_upper_hessenberg()
    assert A.transpose().is_lower_hessenberg()

    A = Matrix([[3, 4, 1],[2, 4 ,5],[3, 1, 2]])
    assert not A.is_upper_hessenberg()

def test_cholesky():
    A = Matrix(((25,15,-5),(15,18,0),(-5,0,11)))
    assert A.cholesky() * A.cholesky().T == A
    assert A.cholesky().is_lower()
    assert A.cholesky() == Matrix([[5, 0, 0], [3, 3, 0], [-1, 1, 3]])

def test_LDLdecomposition():
    A = Matrix(((25,15,-5), (15,18,0), (-5,0,11)))
    L, D = A.LDLdecomposition()
    assert L * D * L.T == A
    assert L.is_lower()
    assert L == Matrix([[1, 0, 0], [ S(3)/5, 1, 0], [S(-1)/5, S(1)/3, 1]])
    assert D.is_diagonal()
    assert D == Matrix([[25, 0, 0], [0, 9, 0], [0, 0, 9]])

def test_cholesky_solve():
    A = Matrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    x = Matrix(3,1,[3,7,5])
    b = A*x
    soln = A.cholesky_solve(b)
    assert soln == x
    A = Matrix([[0,-1,2],
                [5,10,7],
                [8,3,4]])
    x = Matrix(3,1,[-1,2,5])
    b = A*x
    soln = A.cholesky_solve(b)
    assert soln == x

def test_LDLsolve():
    A = Matrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    x = Matrix(3,1,[3,7,5])
    b = A*x
    soln = A.LDLsolve(b)
    assert soln == x
    A = Matrix([[0,-1,2],
                [5,10,7],
                [8,3,4]])
    x = Matrix(3,1,[-1,2,5])
    b = A*x
    soln = A.LDLsolve(b)
    assert soln == x

def test_matrix_norm():
    # Vector Tests
    # Test columns and symbols
    x = Symbol('x', real=True)
    v = Matrix([cos(x), sin(x)])
    assert trigsimp(v.norm(2)) == 1
    assert v.norm(10) == Pow(cos(x)**10 + sin(x)**10, S(1)/10)

    # Test Rows
    y = Matrix([[5, Rational(3,2)]])
    assert y.norm() == Pow(25 + Rational(9,4),S(1)/2)
    assert y.norm(oo) == max(y.mat)
    assert y.norm(-oo) == min(y.mat)

    # Matrix Tests
    # Intuitive test
    A = Matrix([[1,1], [1,1]])
    assert A.norm(2)==2
    assert A.norm(-2)==0
    assert A.norm('frobenius')==2
    assert eye(10).norm(2)==eye(10).norm(-2)==1

    # Test with Symbols and more complex entries
    y = Symbol('y')
    A = Matrix([[3,y,y],[x,S(1)/2, -pi]])
    assert (A.norm('fro')
           == (S(37)/4 + 2*abs(y)**2 + pi**2 + x**2)**(S(1)/2))

    # Check non-square
    A = Matrix([[1,2,-3],[4,5,Rational(13,2)]])
    assert A.norm(2) == sympify('(389/8 + 78665**(1/2)/8)**(1/2)')
    assert A.norm(-2) == S(0)
    assert A.norm('frobenius') == 389**Rational(1,2)/2

    # Test properties of matrix norms
    # http://en.wikipedia.org/wiki/Matrix_norm#Definition
    # Two matrices
    A = Matrix([[1,2],[3,4]])
    B = Matrix([[5,5],[-2,2]])
    C = Matrix([[0,-I],[I,0]])
    D = Matrix([[1,0],[0,-1]])
    L = [A,B,C,D]
    alpha = Symbol('alpha', real=True)

    for order in ['fro', 2, -2]:
        # Zero Check
        assert zeros(3).norm(order) == S(0)
        # Check Triangle Inequality for all Pairs of Matrices
        for X in L:
            for Y in L:
                assert X.norm(order)+Y.norm(order) >= (X+Y).norm(order)
        # Scalar multiplication linearity
        for M in [A,B,C,D]:
            if order in [2,-2]:
                # Abs is causing tests to fail when Abs(alpha) is inside a Max
                # or Min. The tests produce mathematically true statements that
                # are too complex to be simplified well.
                continue;
            try:
                assert ((alpha*M).norm(order) ==
                        abs(alpha) * M.norm(order))
            except NotImplementedError:
                pass; # Some Norms fail on symbolic matrices due to Max issue

    # Test Properties of Vector Norms
    # http://en.wikipedia.org/wiki/Vector_norm
    # Two column vectors
    a = Matrix([1,1-1*I,-3])
    b = Matrix([S(1)/2, 1*I, 1])
    c = Matrix([-1,-1,-1])
    d = Matrix([3, 2, I])
    e = Matrix([Integer(1e2),Rational(1,1e2),1])
    L = [a,b,c,d,e]
    alpha = Symbol('alpha', real=True)

    for order in [1,2,-1, -2, S.Infinity, S.NegativeInfinity, pi]:
        # Zero Check
        assert Matrix([0,0,0]).norm(order) == S(0)
        # Triangle inequality on all pairs
        if order >= 1: # Triangle InEq holds only for these norms
            for v in L:
                for w in L:
                    assert v.norm(order)+w.norm(order) >= (v+w).norm(order)
        # Linear to scalar multiplication
        if order in [1,2, -1, -2, S.Infinity, S.NegativeInfinity]:
            for vec in L:
                try:
                    assert simplify(  (alpha*v).norm(order) -
                            (abs(alpha) * v.norm(order))  ) == 0
                except NotImplementedError:
                    pass; # Some Norms fail on symbolics due to Max issue


def test_singular_values():
    x = Symbol('x', real=True)

    A = Matrix([[0,1*I],[2,0]])
    assert A.singular_values() == [2,1]

    A = eye(3); A[1,1] = x; A[2,2] = 5
    vals = A.singular_values();
    assert 1 in vals and 5 in vals and abs(x) in vals

    A = Matrix([[sin(x), cos(x)],[-cos(x), sin(x)]])
    vals = [sv.trigsimp() for sv in A.singular_values()]
    assert vals == [S(1), S(1)]

def test_condition_number():
    x = Symbol('x', real=True)
    A = eye(3);
    A[0,0] = 10;
    A[2,2] = S(1)/10;
    assert A.condition_number() == 100

    A[1,1] = x
    assert A.condition_number() == Max(10, Abs(x)) / Min(S(1)/10 , Abs(x))

    M = Matrix([[cos(x), sin(x)], [-sin(x), cos(x)]])
    Mc = M.condition_number()
    assert all(Float(1.).epsilon_eq(Mc.subs(x, val).evalf()) for val in \
            [Rational(1,5), Rational(1, 2), Rational(1, 10), pi/2, pi, 7*pi/4 ])

def test_len():
    assert len(Matrix()) == 0
    assert len(Matrix([[1, 2]])) == len(Matrix([[1], [2]])) == 2
    assert len(Matrix(0, 2, lambda i, j: 0)) == len(Matrix(2, 0, lambda i, j: 0)) == 0
    assert len(Matrix([[0, 1, 2], [3, 4, 5]])) == 6
    assert Matrix([1])
    assert not Matrix()

def test_equality():
    A = Matrix(((1,2,3),(4,5,6),(7,8,9)))
    B = Matrix(((9,8,7),(6,5,4),(3,2,1)))
    assert A == A[:, :]
    assert not A != A[:, :]
    assert not A == B
    assert A != B
    assert A != 10
    assert not A == 10

    # A SparseMatrix can be equal to a Matrix
    C = SparseMatrix(((1,0,0),(0,1,0),(0,0,1)))
    D = Matrix(((1,0,0),(0,1,0),(0,0,1)))
    assert C == D
    assert not C != D

