Immutable Matrices
==================

.. module:: sympy.matrices.immutable_matrix

The standard ``Matrix`` class in SymPy is mutable. This is important for
performance reasons but means that standard matrices can not interact well with
the rest of SymPy. This is because the ``Basic`` object, from which most SymPy
classes inherit, is immutable.

To bridge this tension between performance/mutability and safety/immutability
is the ``ImmutableMatrix`` class. Immutable matrices can do almost everything
that normal matrices can do but they inherit from Basic and can thus interact
more naturally with the rest of SymPy. 

You can turn any Matrix-like object into an ``ImmutableMatrix`` by calling the
constructor 

    >>> M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> M[1, 1] = 0
    >>> IM = ImmutableMatrix(M)
    >>> IM
    [1, 2, 3]
    [4, 0, 6]
    [7, 8, 9]
    
    >>> IM[1, 1] = 5
    TypeError: Can not set values in Immutable Matrix

Matrix Expressions
------------------

ImmutableMatrix also inherits from MatrixExpr, allowing it to interact freely
with SymPy's Matrix Expression module.
