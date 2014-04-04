Matrices
========

Creating matrices
-----------------

Basic methods
.............

Matrices in mpmath are implemented using dictionaries. Only non-zero values are
stored, so it is cheap to represent sparse matrices.

The most basic way to create one is to use the ``matrix`` class directly. You
can create an empty matrix specifying the dimensions::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = False
    >>> matrix(2)
    matrix(
    [['0.0', '0.0'],
     ['0.0', '0.0']])
    >>> matrix(2, 3)
    matrix(
    [['0.0', '0.0', '0.0'],
     ['0.0', '0.0', '0.0']])

Calling ``matrix`` with one dimension will create a square matrix.

To access the dimensions of a matrix, use the ``rows`` or ``cols`` keyword::

    >>> A = matrix(3, 2)
    >>> A
    matrix(
    [['0.0', '0.0'],
     ['0.0', '0.0'],
     ['0.0', '0.0']])
    >>> A.rows
    3
    >>> A.cols
    2

You can also change the dimension of an existing matrix. This will set the
new elements to 0. If the new dimension is smaller than before, the
concerning elements are discarded::

    >>> A.rows = 2
    >>> A
    matrix(
    [['0.0', '0.0'],
     ['0.0', '0.0']])

Internally ``convert`` is applied every time an element is set. This is
done using the syntax A[row,column], counting from 0::

    >>> A = matrix(2)
    >>> A[1,1] = 1 + 1j
    >>> print A
    [0.0           0.0]
    [0.0  (1.0 + 1.0j)]

A more comfortable way to create a matrix lets you use nested lists::

    >>> matrix([[1, 2], [3, 4]])
    matrix(
    [['1.0', '2.0'],
     ['3.0', '4.0']])

Advanced methods
................

Convenient functions are available for creating various standard matrices::

    >>> zeros(2)
    matrix(
    [['0.0', '0.0'],
     ['0.0', '0.0']])
    >>> ones(2)
    matrix(
    [['1.0', '1.0'],
     ['1.0', '1.0']])
    >>> diag([1, 2, 3]) # diagonal matrix
    matrix(
    [['1.0', '0.0', '0.0'],
     ['0.0', '2.0', '0.0'],
     ['0.0', '0.0', '3.0']])
    >>> eye(2) # identity matrix
    matrix(
    [['1.0', '0.0'],
     ['0.0', '1.0']])

You can even create random matrices::

    >>> randmatrix(2) # doctest:+SKIP
    matrix(
    [['0.53491598236191806', '0.57195669543302752'],
     ['0.85589992269513615', '0.82444367501382143']])

Vectors
.......

Vectors may also be represented by the ``matrix`` class (with rows = 1 or cols = 1).
For vectors there are some things which make life easier. A column vector can
be created using a flat list, a row vectors using an almost flat nested list::

    >>> matrix([1, 2, 3])
    matrix(
    [['1.0'],
     ['2.0'],
     ['3.0']])
    >>> matrix([[1, 2, 3]])
    matrix(
    [['1.0', '2.0', '3.0']])

Optionally vectors can be accessed like lists, using only a single index::

    >>> x = matrix([1, 2, 3])
    >>> x[1]
    mpf('2.0')
    >>> x[1,0]
    mpf('2.0')

Other
.....

Like you probably expected, matrices can be printed::

    >>> print randmatrix(3) # doctest:+SKIP
    [ 0.782963853573023  0.802057689719883  0.427895717335467]
    [0.0541876859348597  0.708243266653103  0.615134039977379]
    [ 0.856151514955773  0.544759264818486  0.686210904770947]

Use ``nstr`` or ``nprint`` to specify the number of digits to print::

    >>> nprint(randmatrix(5), 3) # doctest:+SKIP
    [2.07e-1  1.66e-1  5.06e-1  1.89e-1  8.29e-1]
    [6.62e-1  6.55e-1  4.47e-1  4.82e-1  2.06e-2]
    [4.33e-1  7.75e-1  6.93e-2  2.86e-1  5.71e-1]
    [1.01e-1  2.53e-1  6.13e-1  3.32e-1  2.59e-1]
    [1.56e-1  7.27e-2  6.05e-1  6.67e-2  2.79e-1]

As matrices are mutable, you will need to copy them sometimes::

    >>> A = matrix(2)
    >>> A
    matrix(
    [['0.0', '0.0'],
     ['0.0', '0.0']])
    >>> B = A.copy()
    >>> B[0,0] = 1
    >>> B
    matrix(
    [['1.0', '0.0'],
     ['0.0', '0.0']])
    >>> A
    matrix(
    [['0.0', '0.0'],
     ['0.0', '0.0']])

Finally, it is possible to convert a matrix to a nested list. This is very useful,
as most Python libraries involving matrices or arrays (namely NumPy or SymPy)
support this format::

    >>> B.tolist()
    [[mpf('1.0'), mpf('0.0')], [mpf('0.0'), mpf('0.0')]]


Matrix operations
-----------------

You can add and substract matrices of compatible dimensions::

    >>> A = matrix([[1, 2], [3, 4]])
    >>> B = matrix([[-2, 4], [5, 9]])
    >>> A + B
    matrix(
    [['-1.0', '6.0'],
     ['8.0', '13.0']])
    >>> A - B
    matrix(
    [['3.0', '-2.0'],
     ['-2.0', '-5.0']])
    >>> A + ones(3) # doctest:+ELLIPSIS
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "...", line 238, in __add__
        raise ValueError('incompatible dimensions for addition')
    ValueError: incompatible dimensions for addition

It is possible to multiply or add matrices and scalars. In the latter case the
operation will be done element-wise::

    >>> A * 2
    matrix(
    [['2.0', '4.0'],
     ['6.0', '8.0']])
    >>> A / 4
    matrix(
    [['0.25', '0.5'],
     ['0.75', '1.0']])
    >>> A - 1
    matrix(
    [['0.0', '1.0'],
     ['2.0', '3.0']])

Of course you can perform matrix multiplication, if the dimensions are
compatible::

    >>> A * B
    matrix(
    [['8.0', '22.0'],
     ['14.0', '48.0']])
    >>> matrix([[1, 2, 3]]) * matrix([[-6], [7], [-2]])
    matrix(
    [['2.0']])

You can raise powers of square matrices::

    >>> A**2
    matrix(
    [['7.0', '10.0'],
     ['15.0', '22.0']])

Negative powers will calculate the inverse::

    >>> A**-1
    matrix(
    [['-2.0', '1.0'],
     ['1.5', '-0.5']])
    >>> nprint(A * A**-1, 3)
    [      1.0  1.08e-19]
    [-2.17e-19       1.0]

Matrix transposition is straightforward::

    >>> A = ones(2, 3)
    >>> A
    matrix(
    [['1.0', '1.0', '1.0'],
     ['1.0', '1.0', '1.0']])
    >>> A.T
    matrix(
    [['1.0', '1.0'],
     ['1.0', '1.0'],
     ['1.0', '1.0']])


Norms
.....

Sometimes you need to know how "large" a matrix or vector is. Due to their
multidimensional nature it's not possible to compare them, but there are
several functions to map a matrix or a vector to a positive real number, the
so called norms.

.. autofunction :: mpmath.norm

.. autofunction :: mpmath.mnorm


Linear algebra
--------------

Decompositions
..............

.. autofunction :: mpmath.cholesky


Linear equations
................

Basic linear algebra is implemented; you can for example solve the linear
equation system::

      x + 2*y = -10
    3*x + 4*y =  10

using ``lu_solve``::

    >>> A = matrix([[1, 2], [3, 4]])
    >>> b = matrix([-10, 10])
    >>> x = lu_solve(A, b)
    >>> x
    matrix(
    [['30.0'],
     ['-20.0']])

If you don't trust the result, use ``residual`` to calculate the residual ||A*x-b||::

    >>> residual(A, x, b)
    matrix(
    [['3.46944695195361e-18'],
     ['3.46944695195361e-18']])
    >>> str(eps)
    '2.22044604925031e-16'

As you can see, the solution is quite accurate. The error is caused by the
inaccuracy of the internal floating point arithmetic. Though, it's even smaller
than the current machine epsilon, which basically means you can trust the
result.

If you need more speed, use NumPy, or use ``fp`` instead ``mp`` matrices
and methods::

    >>> A = fp.matrix([[1, 2], [3, 4]])
    >>> b = fp.matrix([-10, 10])
    >>> fp.lu_solve(A, b)
    matrix(
    [['30.0'],
     ['-20.0']])

``lu_solve`` accepts overdetermined systems. It is usually not possible to solve
such systems, so the residual is minimized instead. Internally this is done
using Cholesky decomposition to compute a least squares approximation. This means
that that ``lu_solve`` will square the errors. If you can't afford this, use
``qr_solve`` instead. It is twice as slow but more accurate, and it calculates
the residual automatically.


Matrix factorization
....................

The function ``lu`` computes an explicit LU factorization of a matrix::

    >>> P, L, U = lu(matrix([[0,2,3],[4,5,6],[7,8,9]]))
    >>> print P
    [0.0  0.0  1.0]
    [1.0  0.0  0.0]
    [0.0  1.0  0.0]
    >>> print L
    [              1.0                0.0  0.0]
    [              0.0                1.0  0.0]
    [0.571428571428571  0.214285714285714  1.0]
    >>> print U
    [7.0  8.0                9.0]
    [0.0  2.0                3.0]
    [0.0  0.0  0.214285714285714]
    >>> print P.T*L*U
    [0.0  2.0  3.0]
    [4.0  5.0  6.0]
    [7.0  8.0  9.0]

The function ``qr`` computes a QR factorization of a matrix::

    >>> A = matrix([[1, 2], [3, 4], [1, 1]])
    >>> Q, R = qr(A)
    >>> print Q
    [-0.301511344577764   0.861640436855329   0.408248290463863]
    [-0.904534033733291  -0.123091490979333  -0.408248290463863]
    [-0.301511344577764  -0.492365963917331   0.816496580927726]
    >>> print R
    [-3.3166247903554  -4.52267016866645]
    [             0.0  0.738548945875996]
    [             0.0                0.0]
    >>> print Q * R
    [1.0  2.0]
    [3.0  4.0]
    [1.0  1.0]
    >>> print chop(Q.T * Q)
    [1.0  0.0  0.0]
    [0.0  1.0  0.0]
    [0.0  0.0  1.0]


The singular value decomposition
................................

The routines ``svd_r`` and ``svd_c`` compute the singular value decomposition
of a real or complex matrix A. ``svd`` is an unified interface calling
either ``svd_r`` or ``svd_c`` depending on whether *A* is real or complex.

Given *A*, two orthogonal (*A* real) or unitary (*A* complex) matrices *U* and *V*
are calculated such that

.. math ::

       A = U S V, \quad U' U = 1, \quad V V' = 1

where *S* is a suitable shaped matrix whose off-diagonal elements are zero.
Here ' denotes the hermitian transpose (i.e. transposition and complex
conjugation). The diagonal elements of *S* are the singular values of *A*,
i.e. the square roots of the eigenvalues of `A' A` or `A A'`.

Examples::

   >>> from mpmath import mp
   >>> A = mp.matrix([[2, -2, -1], [3, 4, -2], [-2, -2, 0]])
   >>> S = mp.svd_r(A, compute_uv = False)
   >>> print S
   [6.0]
   [3.0]
   [1.0]
   >>> U, S, V = mp.svd_r(A)
   >>> print mp.chop(A - U * mp.diag(S) * V)
   [0.0  0.0  0.0]
   [0.0  0.0  0.0]
   [0.0  0.0  0.0]


The Schur decomposition
.......................

This routine computes the Schur decomposition of a square matrix *A*.
Given *A*, a unitary matrix *Q* is determined such that

.. math ::

      Q' A Q = R, \quad Q' Q = Q Q' = 1

where *R* is an upper right triangular matrix. Here ' denotes the
hermitian transpose (i.e. transposition and conjugation).

Examples::

    >>> from mpmath import mp
    >>> A = mp.matrix([[3, -1, 2], [2, 5, -5], [-2, -3, 7]])
    >>> Q, R = mp.schur(A)
    >>> mp.nprint(R, 3) # doctest:+SKIP
    [2.0  0.417  -2.53]
    [0.0    4.0  -4.74]
    [0.0    0.0    9.0]
    >>> print(mp.chop(A - Q * R * Q.transpose_conj()))
    [0.0  0.0  0.0]
    [0.0  0.0  0.0]
    [0.0  0.0  0.0]


The eigenvalue problem
......................

The routine ``eig`` solves the (ordinary) eigenvalue problem for a real or complex
square matrix *A*. Given *A*, a vector *E* and matrices *ER* and *EL* are calculated such that

.. code ::

              A ER[:,i] =         E[i] ER[:,i]
      EL[i,:] A         = EL[i,:] E[i]

*E* contains the eigenvalues of *A*. The columns of *ER* contain the right eigenvectors
of *A* whereas the rows of *EL* contain the left eigenvectors.


Examples::

    >>> from mpmath import mp
    >>> A = mp.matrix([[3, -1, 2], [2, 5, -5], [-2, -3, 7]])
    >>> E, ER = mp.eig(A)
    >>> print(mp.chop(A * ER[:,0] - E[0] * ER[:,0]))
    [0.0]
    [0.0]
    [0.0]
    >>> E, EL, ER = mp.eig(A,left = True, right = True)
    >>> E, EL, ER = mp.eig_sort(E, EL, ER)
    >>> mp.nprint(E)
    [2.0, 4.0, 9.0]
    >>> print(mp.chop(A * ER[:,0] - E[0] * ER[:,0]))
    [0.0]
    [0.0]
    [0.0]
    >>> print(mp.chop( EL[0,:] * A - EL[0,:] * E[0]))
    [0.0  0.0  0.0]


The symmetric eigenvalue problem
................................

The routines ``eigsy`` and ``eighe`` solve the (ordinary) eigenvalue problem
for a real symmetric or complex hermitian square matrix *A*.
``eigh`` is an unified interface for this two functions calling either
``eigsy`` or ``eighe`` depending on whether *A* is real or complex.

Given *A*, an orthogonal (*A* real) or unitary matrix *Q* (*A* complex) is
calculated which diagonalizes A:

.. math ::

        Q' A Q = \operatorname{diag}(E), \quad Q Q' = Q' Q = 1

Here diag(*E*) a is diagonal matrix whose diagonal is *E*.
' denotes the hermitian transpose (i.e. ordinary transposition and
complex conjugation).

The columns of *Q* are the eigenvectors of *A* and *E* contains the eigenvalues:

.. code ::

        A Q[:,i] = E[i] Q[:,i]

Examples::

    >>> from mpmath import mp
    >>> A = mp.matrix([[3, 2], [2, 0]])
    >>> E = mp.eigsy(A, eigvals_only = True)
    >>> print E
    [-1.0]
    [ 4.0]
    >>> A = mp.matrix([[1, 2], [2, 3]])
    >>> E, Q = mp.eigsy(A)                     # alternative: E, Q = mp.eigh(A)
    >>> print mp.chop(A * Q[:,0] - E[0] * Q[:,0])
    [0.0]
    [0.0]
    >>> A = mp.matrix([[1, 2 + 5j], [2 - 5j, 3]])
    >>> E, Q = mp.eighe(A)                     # alternative: E, Q = mp.eigh(A)
    >>> print mp.chop(A * Q[:,0] - E[0] * Q[:,0])
    [0.0]
    [0.0]


Interval and double-precision matrices
--------------------------------------

The ``iv.matrix`` and ``fp.matrix`` classes convert inputs
to intervals and Python floating-point numbers respectively.

Interval matrices can be used to perform linear algebra operations
with rigorous error tracking::

    >>> a = iv.matrix([['0.1','0.3','1.0'],
    ...                ['7.1','5.5','4.8'],
    ...                ['3.2','4.4','5.6']])
    >>>
    >>> b = iv.matrix(['4','0.6','0.5'])
    >>> c = iv.lu_solve(a, b)
    >>> print c
    [  [5.2582327113062393041, 5.2582327113062749951]]
    [[-13.155049396267856583, -13.155049396267821167]]
    [  [7.4206915477497212555, 7.4206915477497310922]]
    >>> print a*c
    [  [3.9999999999999866773, 4.0000000000000133227]]
    [[0.59999999999972430942, 0.60000000000027142733]]
    [[0.49999999999982236432, 0.50000000000018474111]]

Matrix functions
----------------

.. autofunction :: mpmath.expm
.. autofunction :: mpmath.cosm
.. autofunction :: mpmath.sinm
.. autofunction :: mpmath.sqrtm
.. autofunction :: mpmath.logm
.. autofunction :: mpmath.powm
