Matrix Expressions
==================

.. module:: sympy.matrices.expressions

The Matrix expression module allows users to write down statements like

    >>> from sympy import MatrixSymbol, Matrix
    >>> X = MatrixSymbol('X', 3, 3)
    >>> Y = MatrixSymbol('Y', 3, 3)
    >>> (X.T*X).I*Y
    X**-1*X'**-1*Y

    >>> Matrix(X)
    Matrix([
    [X[0, 0], X[0, 1], X[0, 2]],
    [X[1, 0], X[1, 1], X[1, 2]],
    [X[2, 0], X[2, 1], X[2, 2]]])

    >>> (X*Y)[1, 2]
    X[1, 0]*Y[0, 2] + X[1, 1]*Y[1, 2] + X[1, 2]*Y[2, 2]

where ``X`` and ``Y`` are :class:`MatrixSymbol`'s rather than scalar symbols.

Matrix Expressions Core Reference
---------------------------------
.. autoclass:: MatrixExpr
   :members:
.. autoclass:: MatrixSymbol
   :members:
.. autoclass:: MatAdd
   :members:
.. autoclass:: MatMul
   :members:
.. autoclass:: MatPow
   :members:
.. autoclass:: Inverse
   :members:
.. autoclass:: Transpose
   :members:
.. autoclass:: Trace
   :members:
.. autoclass:: FunctionMatrix
   :members:
.. autoclass:: Identity
   :members:
.. autoclass:: ZeroMatrix
   :members:

Block Matrices
--------------

Block matrices allow you to construct larger matrices out of smaller
sub-blocks. They can work with :class:`MatrixExpr` or
:class:`ImmutableMatrix` objects.

.. module:: sympy.matrices.expressions.blockmatrix

.. autoclass:: BlockMatrix
   :members:
.. autoclass:: BlockDiagMatrix
   :members:
.. autofunction:: block_collapse
