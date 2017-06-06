Immutable Matrices
==================

.. module:: sympy

The standard :class:`Matrix` class in SymPy is mutable. This is important for
performance reasons but means that standard matrices can not interact well with
the rest of SymPy. This is because the :class:`Basic` object, from which most
SymPy classes inherit, is immutable.

The mission of the :class:`ImmutableMatrix` class is to bridge the tension
between performance/mutability and safety/immutability. Immutable matrices can
do almost everything that normal matrices can do but they inherit from
:class:`Basic` and can thus interact more naturally with the rest of SymPy.
:class:`ImmutableMatrix` also inherits from :class:`MatrixExpr`, allowing it to
interact freely with SymPy's Matrix Expression module.

You can turn any Matrix-like object into an :class:`ImmutableMatrix` by calling
the constructor

    >>> from sympy import Matrix, ImmutableMatrix
    >>> M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> M[1, 1] = 0
    >>> IM = ImmutableMatrix(M)
    >>> IM
    Matrix([
    [1, 2, 3],
    [4, 0, 6],
    [7, 8, 9]])
    >>> IM[1, 1] = 5
    Traceback (most recent call last):
    ...
    TypeError: Can not set values in Immutable Matrix. Use Matrix instead.


ImmutableMatrix Class Reference
-------------------------------

.. module:: sympy.matrices.immutable

.. autoclass:: ImmutableMatrix
   :members:
