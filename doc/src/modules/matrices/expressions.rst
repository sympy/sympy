Matrix Expressions
==================

.. module:: sympy.matrices.expressions

The Matrix expression module allow users to write down statements like 

    >>> X = MatrixSymbol('X', 3, 3)
    >>> Y = MatrixSymbol('Y', 3, 3)
    >>> (X.T*X).I*Y
    X^-1⋅X'^-1⋅Y

    >>> Matrix(X)
    [X_00, X_01, X_02]
    [X_10, X_11, X_12]
    [X_20, X_21, X_22]

    >>> (X*Y)[1, 2]
    X_10*Y_02 + X_11*Y_12 + X_12*Y_22

where `X` and `Y` are matrix symbols rather than scalar symbols.

Core Class Reference
--------------------
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

.. module:: sympy.matrices.expressions.blockmatrix

Block matrices allow you to construct larger matrices out of smaller
sub-blocks. They can work with MatrixExprs or ImmutableMatrices

.. autoclass:: BlockMatrix 
   :members:
.. autoclass:: BlockDiagMatrix 
   :members:
.. autofunction:: block_collapse

