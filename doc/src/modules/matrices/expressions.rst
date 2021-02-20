Matrix Expressions
==================

.. module:: sympy.matrices.expressions

The Matrix expression module allows users to write down statements like

    >>> from sympy import MatrixSymbol, Matrix
    >>> X = MatrixSymbol('X', 3, 3)
    >>> Y = MatrixSymbol('Y', 3, 3)
    >>> (X.T*X).I*Y
    X**(-1)*X.T**(-1)*Y

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
   :undoc-members:
   :private-members:
.. autoclass:: MatrixSymbol
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: MatAdd
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: MatMul
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: MatPow
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: HadamardProduct
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: HadamardPower
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: Inverse
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: Transpose
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: Trace
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: FunctionMatrix
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: PermutationMatrix
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: MatrixPermute
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: Identity
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: ZeroMatrix
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: CompanionMatrix
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: MatrixSet
   :members:
   :undoc-members:
   :private-members:

Block Matrices
--------------

Block matrices allow you to construct larger matrices out of smaller
sub-blocks. They can work with :class:`MatrixExpr` or
:obj:`~.ImmutableMatrix` objects.

.. module:: sympy.matrices.expressions.blockmatrix

.. autoclass:: BlockMatrix
   :members:
   :undoc-members:
   :private-members:
.. autoclass:: BlockDiagMatrix
   :members:
   :undoc-members:
   :private-members:
.. autofunction:: block_collapse
