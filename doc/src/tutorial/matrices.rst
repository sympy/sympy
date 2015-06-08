==========
 Matrices
==========

    >>> from sympy import *
    >>> init_printing(use_unicode=True)

To make a matrix in SymPy, use the ``Matrix`` object.  A matrix is constructed
by providing a list of row vectors that make up the matrix.  For example,
to construct the matrix

.. math::


   \left[\begin{array}{cc}1 & -1\\3 & 4\\0 & 2\end{array}\right]

use

    >>> Matrix([[1, -1], [3, 4], [0, 2]])
    ⎡1  -1⎤
    ⎢     ⎥
    ⎢3  4 ⎥
    ⎢     ⎥
    ⎣0  2 ⎦

To make it easy to make column vectors, a list of elements is considered to be
a column vector.

    >>> Matrix([1, 2, 3])
    ⎡1⎤
    ⎢ ⎥
    ⎢2⎥
    ⎢ ⎥
    ⎣3⎦

Matrices are manipulated just like any other object in SymPy or Python.

    >>> M = Matrix([[1, 2, 3], [3, 2, 1]])
    >>> N = Matrix([0, 1, 1])
    >>> M*N
    ⎡5⎤
    ⎢ ⎥
    ⎣3⎦

One important thing to note about SymPy matrices is that, unlike every other
object in SymPy, they are mutable.  This means that they can be modified in
place, as we will see below.  The downside to this is that ``Matrix`` cannot
be used in places that require immutability, such as inside other SymPy
expressions or as keys to dictionaries.  If you need an immutable version of
``Matrix``, use ``ImmutableMatrix``.

Basic Operations
================

Shape
-----

Here are some basic operations on ``Matrix``.  To get the shape of a matrix
use ``shape``

    >>> M = Matrix([[1, 2, 3], [-2, 0, 4]])
    >>> M
    ⎡1   2  3⎤
    ⎢        ⎥
    ⎣-2  0  4⎦
    >>> M.shape
    (2, 3)

Accessing Rows and Columns
--------------------------

To get an individual row or column of a matrix, use ``row`` or ``col``.  For
example, ``M.row(0)`` will get the first row. ``M.col(-1)`` will get the last
column.

    >>> M.row(0)
    [1  2  3]
    >>> M.col(-1)
    ⎡3⎤
    ⎢ ⎥
    ⎣4⎦

Deleting and Inserting Rows and Columns
---------------------------------------

To delete a row or column, use ``row_del`` or ``col_del``.  These operations
will modify the Matrix **in place**.

    >>> M.col_del(0)
    >>> M
    ⎡2  3⎤
    ⎢    ⎥
    ⎣0  4⎦
    >>> M.row_del(1)
    >>> M
    [2  3]

.. TODO: This is a mess. See issue 6992.

To insert rows or columns, use ``row_insert`` or ``col_insert``.  These
operations will modify the Matrix **in place**.

    >>> M
    [2  3]
    >>> M.row_insert(1, Matrix([[0, 4]]))
    >>> M
    ⎡2  3⎤
    ⎢    ⎥
    ⎣0  4⎦
    >>> M.col_insert(0, Matrix([1, -2]))
    >>> M
    ⎡1   2  3⎤
    ⎢        ⎥
    ⎣-2  0  4⎦

Unless explicitly stated, the methods mentioned below do not operate in
place. In general, a method that does not operate in place will return a new
``Matrix`` and a method that does operate in place will return ``None``.

Basic Methods
=============

As noted above, simple operations like addition and multiplication are done
just by using ``+``, ``*``, and ``**``.  To find the inverse of a matrix, just
raise it to the ``-1`` power.

    >>> M = Matrix([[1, 3], [-2, 3]])
    >>> N = Matrix([[0, 3], [0, 7]])
    >>> M + N
    ⎡1   6 ⎤
    ⎢      ⎥
    ⎣-2  10⎦
    >>> M*N
    ⎡0  24⎤
    ⎢     ⎥
    ⎣0  15⎦
    >>> 3*M
    ⎡3   9⎤
    ⎢     ⎥
    ⎣-6  9⎦
    >>> M**2
    ⎡-5  12⎤
    ⎢      ⎥
    ⎣-8  3 ⎦
    >>> M**-1
    ⎡1/3  -1/3⎤
    ⎢         ⎥
    ⎣2/9  1/9 ⎦
    >>> N**-1
    Traceback (most recent call last):
    ...
    ValueError: Matrix det == 0; not invertible.

To take the transpose of a Matrix, use ``T``.

    >>> M = Matrix([[1, 2, 3], [4, 5, 6]])
    >>> M
    ⎡1  2  3⎤
    ⎢       ⎥
    ⎣4  5  6⎦
    >>> M.T
    ⎡1  4⎤
    ⎢    ⎥
    ⎢2  5⎥
    ⎢    ⎥
    ⎣3  6⎦

Matrix Constructors
===================

Several constructors exist for creating common matrices.  To create an
identity matrix, use ``eye``.  ``eye(n)`` will create an `n\times n` identity matrix.

    >>> eye(3)
    ⎡1  0  0⎤
    ⎢       ⎥
    ⎢0  1  0⎥
    ⎢       ⎥
    ⎣0  0  1⎦
    >>> eye(4)
    ⎡1  0  0  0⎤
    ⎢          ⎥
    ⎢0  1  0  0⎥
    ⎢          ⎥
    ⎢0  0  1  0⎥
    ⎢          ⎥
    ⎣0  0  0  1⎦

To create a matrix of all zeros, use ``zeros``.  ``zeros(n, m)`` creates an
`n\times m` matrix of `0`\ s.

    >>> zeros(2, 3)
    ⎡0  0  0⎤
    ⎢       ⎥
    ⎣0  0  0⎦

Similarly, ``ones`` creates a matrix of ones.

    >>> ones(3, 2)
    ⎡1  1⎤
    ⎢    ⎥
    ⎢1  1⎥
    ⎢    ⎥
    ⎣1  1⎦

To create diagonal matrices, use ``diag``.  The arguments to ``diag`` can be
either numbers or matrices.  A number is interpreted as a `1\times 1`
matrix. The matrices are stacked diagonally.  The remaining elements are
filled with `0`\ s.

    >>> diag(1, 2, 3)
    ⎡1  0  0⎤
    ⎢       ⎥
    ⎢0  2  0⎥
    ⎢       ⎥
    ⎣0  0  3⎦
    >>> diag(-1, ones(2, 2), Matrix([5, 7, 5]))
    ⎡-1  0  0  0⎤
    ⎢           ⎥
    ⎢0   1  1  0⎥
    ⎢           ⎥
    ⎢0   1  1  0⎥
    ⎢           ⎥
    ⎢0   0  0  5⎥
    ⎢           ⎥
    ⎢0   0  0  7⎥
    ⎢           ⎥
    ⎣0   0  0  5⎦

Advanced Methods
================

Determinant
-----------

To compute the determinant of a matrix, use ``det``.

    >>> M = Matrix([[1, 0, 1], [2, -1, 3], [4, 3, 2]])
    >>> M
    ⎡1  0   1⎤
    ⎢        ⎥
    ⎢2  -1  3⎥
    ⎢        ⎥
    ⎣4  3   2⎦
    >>> M.det()
    -1

RREF
----

To put a matrix into reduced row echelon form, use ``rref``.  ``rref`` returns
a tuple of two elements. The first is the reduced row echelon form, and the
second is a list of indices of the pivot columns.

    >>> M = Matrix([[1, 0, 1, 3], [2, 3, 4, 7], [-1, -3, -3, -4]])
    >>> M
    ⎡1   0   1   3 ⎤
    ⎢              ⎥
    ⎢2   3   4   7 ⎥
    ⎢              ⎥
    ⎣-1  -3  -3  -4⎦
    >>> M.rref()
    ⎛⎡1  0   1    3 ⎤, [0, 1]⎞
    ⎜⎢              ⎥        ⎟
    ⎜⎢0  1  2/3  1/3⎥        ⎟
    ⎜⎢              ⎥        ⎟
    ⎝⎣0  0   0    0 ⎦        ⎠

.. Note:: The first element of the tuple returned by ``rref`` is of type
   ``Matrix``. The second is of type ``list``.

Nullspace
---------

To find the nullspace of a matrix, use ``nullspace``. ``nullspace`` returns a
``list`` of column vectors that span the nullspace of the matrix.

    >>> M = Matrix([[1, 2, 3, 0, 0], [4, 10, 0, 0, 1]])
    >>> M
    ⎡1  2   3  0  0⎤
    ⎢              ⎥
    ⎣4  10  0  0  1⎦
    >>> M.nullspace()
    ⎡⎡-15⎤, ⎡0⎤, ⎡ 1  ⎤⎤
    ⎢⎢   ⎥  ⎢ ⎥  ⎢    ⎥⎥
    ⎢⎢ 6 ⎥  ⎢0⎥  ⎢-1/2⎥⎥
    ⎢⎢   ⎥  ⎢ ⎥  ⎢    ⎥⎥
    ⎢⎢ 1 ⎥  ⎢0⎥  ⎢ 0  ⎥⎥
    ⎢⎢   ⎥  ⎢ ⎥  ⎢    ⎥⎥
    ⎢⎢ 0 ⎥  ⎢1⎥  ⎢ 0  ⎥⎥
    ⎢⎢   ⎥  ⎢ ⎥  ⎢    ⎥⎥
    ⎣⎣ 0 ⎦  ⎣0⎦  ⎣ 1  ⎦⎦

Columnspace
-----------

To find the columnspace of a matrix, use ``columnspace``. ``columnspace`` returns a
``list`` of column vectors that span the columnspace of the matrix.

    >>> M = Matrix([[1, 1, 2], [2 ,1 , 3], [3 , 1, 4]])
    >>> M
    ⎡1  1  2⎤
    ⎢       ⎥
    ⎢2  1  3⎥
    ⎢       ⎥
    ⎣3  1  4⎦
    >>> M.columnspace()
    ⎡⎡1⎤, ⎡1⎤⎤
    ⎢⎢ ⎥  ⎢ ⎥⎥
    ⎢⎢2⎥  ⎢1⎥⎥
    ⎢⎢ ⎥  ⎢ ⎥⎥
    ⎣⎣3⎦  ⎣1⎦⎦

Eigenvalues, Eigenvectors, and Diagonalization
----------------------------------------------

To find the eigenvalues of a matrix, use ``eigenvals``.  ``eigenvals``
returns a dictionary of ``eigenvalue:algebraic multiplicity`` pairs (similar to the
output of :ref:`roots <tutorial-roots>`).

    >>> M = Matrix([[3, -2,  4, -2], [5,  3, -3, -2], [5, -2,  2, -2], [5, -2, -3,  3]])
    >>> M
    ⎡3  -2  4   -2⎤
    ⎢             ⎥
    ⎢5  3   -3  -2⎥
    ⎢             ⎥
    ⎢5  -2  2   -2⎥
    ⎢             ⎥
    ⎣5  -2  -3  3 ⎦
    >>> M.eigenvals()
    {-2: 1, 3: 1, 5: 2}

This means that ``M`` has eigenvalues -2, 3, and 5, and that the
eigenvalues -2 and 3 have algebraic multiplicity 1 and that the eigenvalue 5
has algebraic multiplicity 2.

To find the eigenvectors of a matrix, use ``eigenvects``.  ``eigenvects``
returns a list of tuples of the form ``(eigenvalue:algebraic multiplicity,
[eigenvectors])``.

    >>> M.eigenvects()
    ⎡⎛-2, 1, ⎡⎡0⎤⎤⎞, ⎛3, 1, ⎡⎡1⎤⎤⎞, ⎛5, 2, ⎡⎡1⎤, ⎡0 ⎤⎤⎞⎤
    ⎢⎜       ⎢⎢ ⎥⎥⎟  ⎜      ⎢⎢ ⎥⎥⎟  ⎜      ⎢⎢ ⎥  ⎢  ⎥⎥⎟⎥
    ⎢⎜       ⎢⎢1⎥⎥⎟  ⎜      ⎢⎢1⎥⎥⎟  ⎜      ⎢⎢1⎥  ⎢-1⎥⎥⎟⎥
    ⎢⎜       ⎢⎢ ⎥⎥⎟  ⎜      ⎢⎢ ⎥⎥⎟  ⎜      ⎢⎢ ⎥  ⎢  ⎥⎥⎟⎥
    ⎢⎜       ⎢⎢1⎥⎥⎟  ⎜      ⎢⎢1⎥⎥⎟  ⎜      ⎢⎢1⎥  ⎢0 ⎥⎥⎟⎥
    ⎢⎜       ⎢⎢ ⎥⎥⎟  ⎜      ⎢⎢ ⎥⎥⎟  ⎜      ⎢⎢ ⎥  ⎢  ⎥⎥⎟⎥
    ⎣⎝       ⎣⎣1⎦⎦⎠  ⎝      ⎣⎣1⎦⎦⎠  ⎝      ⎣⎣0⎦  ⎣1 ⎦⎦⎠⎦

This shows us that, for example, the eigenvalue 5 also has geometric
multiplicity 2, because it has two eigenvectors.  Because the algebraic and
geometric multiplicities are the same for all the eigenvalues, ``M`` is
diagonalizable.

To diagonalize a matrix, use ``diagonalize``. ``diagonalize`` returns a tuple
`(P, D)`, where `D` is diagonal and `M = PDP^{-1}`.

    >>> P, D = M.diagonalize()
    >>> P
    ⎡0  1  1  0 ⎤
    ⎢           ⎥
    ⎢1  1  1  -1⎥
    ⎢           ⎥
    ⎢1  1  1  0 ⎥
    ⎢           ⎥
    ⎣1  1  0  1 ⎦
    >>> D
    ⎡-2  0  0  0⎤
    ⎢           ⎥
    ⎢0   3  0  0⎥
    ⎢           ⎥
    ⎢0   0  5  0⎥
    ⎢           ⎥
    ⎣0   0  0  5⎦
    >>> P*D*P**-1
    ⎡3  -2  4   -2⎤
    ⎢             ⎥
    ⎢5  3   -3  -2⎥
    ⎢             ⎥
    ⎢5  -2  2   -2⎥
    ⎢             ⎥
    ⎣5  -2  -3  3 ⎦
    >>> P*D*P**-1 == M
    True

.. sidebar:: Quick Tip

   ``lambda`` is a reserved keyword in Python, so to create a Symbol called
   `\lambda`, while using the same names for SymPy Symbols and Python
   variables, use ``lamda`` (without the ``b``).  It will still pretty print
   as `\lambda`.

Note that since ``eigenvects`` also includes the eigenvalues, you should use
it instead of ``eigenvals`` if you also want the eigenvectors. However, as
computing the eigenvectors may often be costly, ``eigenvals`` should be
preferred if you only wish to find the eigenvalues.

If all you want is the characteristic polynomial, use ``charpoly``.  This is
more efficient than ``eigenvals``, because sometimes symbolic roots can be
expensive to calculate.

    >>> lamda = symbols('lamda')
    >>> p = M.charpoly(lamda)
    >>> factor(p)
           2
    (λ - 5) ⋅(λ - 3)⋅(λ + 2)

.. TODO: Add an example for ``jordan_form``, once it is fully implemented.
