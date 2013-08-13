Matrices (linear algebra)
=========================

.. module:: sympy.matrices.matrices

Creating Matrices
-----------------

The linear algebra module is designed to be as simple as possible. First, we
import and declare our first ``Matrix`` object:

    >>> from sympy.interactive.printing import init_printing
    >>> init_printing(use_unicode=False, wrap_line=False, no_global=True)
    >>> from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt
    >>> M = Matrix([[1,0,0], [0,0,0]]); M
    [1  0  0]
    [       ]
    [0  0  0]
    >>> Matrix([M, (0, 0, -1)])
    [1  0  0 ]
    [        ]
    [0  0  0 ]
    [        ]
    [0  0  -1]
    >>> Matrix([[1, 2, 3]])
    [1 2 3]
    >>> Matrix([1, 2, 3])
    [1]
    [ ]
    [2]
    [ ]
    [3]

In addition to creating a matrix from a list of appropriately-sized lists
and/or matrices, SymPy also supports more advanced methods of matrix creation
including a single list of values and dimension inputs:

    >>> Matrix(2, 3, [1, 2, 3, 4, 5, 6])
    [1  2  3]
    [       ]
    [4  5  6]

More interesting (and useful), is the ability to use a 2-variable function
(or ``lambda``) to create a matrix. Here we create an indicator function which
is 1 on the diagonal and then use it to make the identity matrix:

    >>> def f(i,j):
    ...     if i == j:
    ...         return 1
    ...     else:
    ...         return 0
    ...
    >>> Matrix(4, 4, f)
    [1  0  0  0]
    [          ]
    [0  1  0  0]
    [          ]
    [0  0  1  0]
    [          ]
    [0  0  0  1]

Finally let's use ``lambda`` to create a 1-line matrix with 1's in the even
permutation entries:

    >>> Matrix(3, 4, lambda i,j: 1 - (i+j) % 2)
    [1  0  1  0]
    [          ]
    [0  1  0  1]
    [          ]
    [1  0  1  0]

There are also a couple of special constructors for quick matrix construction:
``eye`` is the identity matrix, ``zeros`` and ``ones`` for matrices of all
zeros and ones, respectively, and ``diag`` to put matrices or elements along
the diagonal:

    >>> eye(4)
    [1  0  0  0]
    [          ]
    [0  1  0  0]
    [          ]
    [0  0  1  0]
    [          ]
    [0  0  0  1]
    >>> zeros(2)
    [0  0]
    [    ]
    [0  0]
    >>> zeros(2, 5)
    [0  0  0  0  0]
    [             ]
    [0  0  0  0  0]
    >>> ones(3)
    [1  1  1]
    [       ]
    [1  1  1]
    [       ]
    [1  1  1]
    >>> ones(1, 3)
    [1  1  1]
    >>> diag(1, Matrix([[1, 2], [3, 4]]))
    [1  0  0]
    [       ]
    [0  1  2]
    [       ]
    [0  3  4]


Basic Manipulation
------------------

While learning to work with matrices, let's choose one where the entries are
readily identifiable. One useful thing to know is that while matrices are
2-dimensional, the storage is not and so it is allowable - though one should be
careful - to access the entries as if they were a 1-d list.

    >>> M = Matrix(2, 3, [1, 2, 3, 4, 5, 6])
    >>> M[4]
    5

Now, the more standard entry access is a pair of indices which will always
return the value at the corresponding row and column of the matrix:

    >>> M[1, 2]
    6
    >>> M[0, 0]
    1
    >>> M[1, 1]
    5

Since this is Python we're also able to slice submatrices; slices always
give a matrix in return, even if the dimension is 1 x 1::

    >>> M[0:2, 0:2]
    [1  2]
    [    ]
    [4  5]
    >>> M[2:2, 2]
    []
    >>> M[:, 2]
    [3]
    [ ]
    [6]
    >>> M[:1, 2]
    [3]

In the second example above notice that the slice 2:2 gives an empty range. Note
also (in keeping with 0-based indexing of Python) the first row/column is 0.

You cannot access rows or columns that are not present unless they
are in a slice:

    >>> M[:, 10] # the 10-th column (not there)
    Traceback (most recent call last):
    ...
    IndexError: Index out of range: a[[0, 10]]
    >>> M[:, 10:11] # the 10-th column (if there)
    []
    >>> M[:, :10] # all columns up to the 10-th
    [1  2  3]
    [       ]
    [4  5  6]

Slicing an empty matrix works as long as you use a slice for the coordinate
that has no size:

    >>> Matrix(0, 3, [])[:, 1]
    []

Slicing gives a copy of what is sliced, so modifications of one object
do not affect the other:

    >>> M2 = M[:, :]
    >>> M2[0, 0] = 100
    >>> M[0, 0] == 100
    False

Notice that changing ``M2`` didn't change ``M``. Since we can slice, we can also assign
entries:

    >>> M = Matrix(([1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]))
    >>> M
    [1   2   3   4 ]
    [              ]
    [5   6   7   8 ]
    [              ]
    [9   10  11  12]
    [              ]
    [13  14  15  16]
    >>> M[2,2] = M[0,3] = 0
    >>> M
    [1   2   3   0 ]
    [              ]
    [5   6   7   8 ]
    [              ]
    [9   10  0   12]
    [              ]
    [13  14  15  16]

as well as assign slices:

    >>> M = Matrix(([1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]))
    >>> M[2:,2:] = Matrix(2,2,lambda i,j: 0)
    >>> M
    [1   2   3  4]
    [            ]
    [5   6   7  8]
    [            ]
    [9   10  0  0]
    [            ]
    [13  14  0  0]

All the standard arithmetic operations are supported:

    >>> M = Matrix(([1,2,3],[4,5,6],[7,8,9]))
    >>> M - M
    [0  0  0]
    [       ]
    [0  0  0]
    [       ]
    [0  0  0]
    >>> M + M
    [2   4   6 ]
    [          ]
    [8   10  12]
    [          ]
    [14  16  18]
    >>> M * M
    [30   36   42 ]
    [             ]
    [66   81   96 ]
    [             ]
    [102  126  150]
    >>> M2 = Matrix(3,1,[1,5,0])
    >>> M*M2
    [11]
    [  ]
    [29]
    [  ]
    [47]
    >>> M**2
    [30   36   42 ]
    [             ]
    [66   81   96 ]
    [             ]
    [102  126  150]

As well as some useful vector operations:

    >>> M.row_del(0)
    >>> M
    [4  5  6]
    [       ]
    [7  8  9]
    >>> M.col_del(1)
    >>> M
    [4  6]
    [    ]
    [7  9]
    >>> v1 = Matrix([1,2,3])
    >>> v2 = Matrix([4,5,6])
    >>> v3 = v1.cross(v2)
    >>> v1.dot(v2)
    32
    >>> v2.dot(v3)
    0
    >>> v1.dot(v3)
    0

Recall that the ``row_del()`` and ``col_del()`` operations don't return a value - they
simply change the matrix object. We can also ''glue'' together matrices of the
appropriate size:

    >>> M1 = eye(3)
    >>> M2 = zeros(3, 4)
    >>> M1.row_join(M2)
    [1  0  0  0  0  0  0]
    [                   ]
    [0  1  0  0  0  0  0]
    [                   ]
    [0  0  1  0  0  0  0]
    >>> M3 = zeros(4, 3)
    >>> M1.col_join(M3)
    [1  0  0]
    [       ]
    [0  1  0]
    [       ]
    [0  0  1]
    [       ]
    [0  0  0]
    [       ]
    [0  0  0]
    [       ]
    [0  0  0]
    [       ]
    [0  0  0]


Operations on entries
---------------------

We are not restricted to having multiplication between two matrices:

    >>> M = eye(3)
    >>> 2*M
    [2  0  0]
    [       ]
    [0  2  0]
    [       ]
    [0  0  2]
    >>> 3*M
    [3  0  0]
    [       ]
    [0  3  0]
    [       ]
    [0  0  3]

but we can also apply functions to our matrix entries using ``applyfunc()``. Here we'll declare a function that double any input number. Then we apply it to the 3x3 identity matrix:

    >>> f = lambda x: 2*x
    >>> eye(3).applyfunc(f)
    [2  0  0]
    [       ]
    [0  2  0]
    [       ]
    [0  0  2]

One more useful matrix-wide entry application function is the substitution function. Let's declare a matrix with symbolic entries then substitute a value. Remember we can substitute anything - even another symbol!:

    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> M = eye(3) * x
    >>> M
    [x  0  0]
    [       ]
    [0  x  0]
    [       ]
    [0  0  x]
    >>> M.subs(x, 4)
    [4  0  0]
    [       ]
    [0  4  0]
    [       ]
    [0  0  4]
    >>> y = Symbol('y')
    >>> M.subs(x, y)
    [y  0  0]
    [       ]
    [0  y  0]
    [       ]
    [0  0  y]


Linear algebra
--------------

Now that we have the basics out of the way, let's see what we can do with the
actual matrices. Of course, one of the first things that comes to mind is the
determinant:

    >>> M = Matrix(( [1, 2, 3], [3, 6, 2], [2, 0, 1] ))
    >>> M.det()
    -28
    >>> M2 = eye(3)
    >>> M2.det()
    1
    >>> M3 = Matrix(( [1, 0, 0], [1, 0, 0], [1, 0, 0] ))
    >>> M3.det()
    0

Another common operation is the inverse: In SymPy, this is computed by Gaussian
elimination by default (for dense matrices) but we can specify it be done by `LU`
decomposition as well:

    >>> M2.inv()
    [1  0  0]
    [       ]
    [0  1  0]
    [       ]
    [0  0  1]
    >>> M2.inv(method="LU")
    [1  0  0]
    [       ]
    [0  1  0]
    [       ]
    [0  0  1]
    >>> M.inv(method="LU")
    [-3/14  1/14  1/2 ]
    [                 ]
    [-1/28  5/28  -1/4]
    [                 ]
    [ 3/7   -1/7   0  ]
    >>> M * M.inv(method="LU")
    [1  0  0]
    [       ]
    [0  1  0]
    [       ]
    [0  0  1]

We can perform a `QR` factorization which is handy for solving systems:

    >>> A = Matrix([[1,1,1],[1,1,3],[2,3,4]])
    >>> Q, R = A.QRdecomposition()
    >>> Q
    [  ___     ___      ___ ]
    [\/ 6   -\/ 3    -\/ 2  ]
    [-----  -------  -------]
    [  6       3        2   ]
    [                       ]
    [  ___     ___      ___ ]
    [\/ 6   -\/ 3     \/ 2  ]
    [-----  -------   ----- ]
    [  6       3        2   ]
    [                       ]
    [  ___     ___          ]
    [\/ 6    \/ 3           ]
    [-----   -----      0   ]
    [  3       3            ]
    >>> R
    [           ___         ]
    [  ___  4*\/ 6       ___]
    [\/ 6   -------  2*\/ 6 ]
    [          3            ]
    [                       ]
    [          ___          ]
    [        \/ 3           ]
    [  0     -----      0   ]
    [          3            ]
    [                       ]
    [                   ___ ]
    [  0       0      \/ 2  ]
    >>> Q*R
    [1  1  1]
    [       ]
    [1  1  3]
    [       ]
    [2  3  4]


In addition to the solvers in the ``solver.py`` file, we can solve the system Ax=b
by passing the b vector to the matrix A's LUsolve function. Here we'll cheat a
little choose A and x then multiply to get b. Then we can solve for x and check
that it's correct:

    >>> A = Matrix([ [2, 3, 5], [3, 6, 2], [8, 3, 6] ])
    >>> x = Matrix(3,1,[3,7,5])
    >>> b = A*x
    >>> soln = A.LUsolve(b)
    >>> soln
    [3]
    [ ]
    [7]
    [ ]
    [5]

There's also a nice Gram-Schmidt orthogonalizer which will take a set of
vectors and orthogonalize then with respect to another another. There is an
optional argument which specifies whether or not the output should also be
normalized, it defaults to ``False``. Let's take some vectors and orthogonalize
them - one normalized and one not:

    >>> L = [Matrix([2,3,5]), Matrix([3,6,2]), Matrix([8,3,6])]
    >>> out1 = GramSchmidt(L)
    >>> out2 = GramSchmidt(L, True)

Let's take a look at the vectors:

    >>> for i in out1:
    ...     print(i)
    ...
    Matrix([[2], [3], [5]])
    Matrix([[23/19], [63/19], [-47/19]])
    Matrix([[1692/353], [-1551/706], [-423/706]])
    >>> for i in out2:
    ...      print(i)
    ...
    Matrix([[sqrt(38)/19], [3*sqrt(38)/38], [5*sqrt(38)/38]])
    Matrix([[23*sqrt(6707)/6707], [63*sqrt(6707)/6707], [-47*sqrt(6707)/6707]])
    Matrix([[12*sqrt(706)/353], [-11*sqrt(706)/706], [-3*sqrt(706)/706]])

We can spot-check their orthogonality with dot() and their normality with
norm():

    >>> out1[0].dot(out1[1])
    0
    >>> out1[0].dot(out1[2])
    0
    >>> out1[1].dot(out1[2])
    0
    >>> out2[0].norm()
    1
    >>> out2[1].norm()
    1
    >>> out2[2].norm()
    1

So there is quite a bit that can be done with the module including eigenvalues,
eigenvectors, nullspace calculation, cofactor expansion tools, and so on. From
here one might want to look over the ``matrices.py`` file for all functionality.

MatrixBase Class Reference
--------------------------
.. autoclass:: MatrixBase
   :members:

Matrix Exceptions Reference
---------------------------

.. autoclass:: MatrixError

.. autoclass:: ShapeError

.. autoclass:: NonSquareMatrixError


Matrix Functions Reference
--------------------------

.. autofunction:: classof

.. autofunction:: sympy.matrices.dense.matrix_multiply_elementwise

.. autofunction:: sympy.matrices.dense.zeros

.. autofunction:: sympy.matrices.dense.ones

.. autofunction:: sympy.matrices.dense.eye

.. autofunction:: sympy.matrices.dense.diag

.. autofunction:: sympy.matrices.dense.jordan_cell

.. autofunction:: sympy.matrices.dense.hessian

.. autofunction:: sympy.matrices.dense.GramSchmidt

.. autofunction:: sympy.matrices.dense.wronskian

.. autofunction:: sympy.matrices.dense.casoratian

.. autofunction:: sympy.matrices.dense.randMatrix

Numpy Utility Functions Reference
---------------------------------

.. autofunction:: sympy.matrices.dense.list2numpy

.. autofunction:: sympy.matrices.dense.matrix2numpy

.. autofunction:: sympy.matrices.dense.symarray

.. autofunction:: sympy.matrices.dense.rot_axis1

.. autofunction:: sympy.matrices.dense.rot_axis2

.. autofunction:: sympy.matrices.dense.rot_axis3

.. autofunction:: a2idx
