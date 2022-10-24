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

Matrix expression derivatives are supported. The derivative of a matrix by another matrix
is generally a 4-dimensional array, but if some dimensions are trivial or diagonal,
the derivation algorithm will try to express the result as a matrix expression:

    >>> a = MatrixSymbol("a", 3, 1)
    >>> b = MatrixSymbol("b", 3, 1)
    >>> (a.T*X**2*b).diff(X)
    a*b.T*X.T + X.T*a*b.T

    >>> X.diff(X)
    PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))

The last output is an array expression, as the returned symbol
is 4-dimensional.

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
.. autofunction:: hadamard_product
.. autoclass:: HadamardProduct
   :members:
.. autoclass:: HadamardPower
   :members:
.. autoclass:: Inverse
   :members:
.. autoclass:: Transpose
   :members:
.. autoclass:: Trace
   :members:
.. autoclass:: FunctionMatrix
   :members:
.. autoclass:: PermutationMatrix
   :members:
.. autoclass:: MatrixPermute
   :members:
.. autoclass:: Identity
   :members:
.. autoclass:: ZeroMatrix
   :members:
.. autoclass:: CompanionMatrix
   :members:
.. autoclass:: MatrixSet
   :members:

Block Matrices
--------------

Block matrices allow you to construct larger matrices out of smaller
sub-blocks. They can work with :class:`MatrixExpr` or
:obj:`~.ImmutableMatrix` objects.

.. module:: sympy.matrices.expressions.blockmatrix

.. autoclass:: BlockMatrix
   :members:
.. autoclass:: BlockDiagMatrix
   :members:
.. autofunction:: block_collapse
