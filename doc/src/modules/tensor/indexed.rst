===============
Indexed Objects
===============

The classes ``IndexedBase``, ``Indexed``, and ``Idx`` represent a
matrix element ``M[i, j]`` as in the following diagram::

   1) The Indexed class represents the entire indexed object.
              |
           ___|___
          '       '
           M[i, j]
          /   \__\______
          |             |
          |             |
          |     2) The Idx class represents indices; each Idx can
          |        optionally contain information about its range.
          |
    3) IndexedBase represents the 'stem' of an Indexed object, here M.

There can be any number of indices on an ``Indexed`` object.  No symmetry
properties are implemented in these objects. Implicit contraction of
repeated indices is supported via the ``EinsteinSum`` function.

Note that the support for complicated (i.e. non-atomic) integer
expressions as indices is limited.  (This should be improved in
future releases.)

Basic Examples
==============

To express the above matrix element example:

>>> from sympy import symbols, IndexedBase, Idx
>>> M = IndexedBase('M')
>>> i, j = symbols('i j', cls=Idx)
>>> M[i, j]
M[i, j]

Repeated indices in a monomial imply summation when enclosed by the
``EinsteinSum`` function. To express a matrix-vector product in terms of
``Indexed`` objects:

>>> from sympy import EinsteinSum
>>> x = IndexedBase('x')
>>> EinsteinSum(M[i, j] * x[j])
EinsteinSum(x[j]*M[i, j])

``EinsteinSum`` expressions can be differentiated:

>>> EinsteinSum(x[j] * x[j]).diff(x[i])
EinsteinSum(2*x[i])

They can also be converted into Python functions accepting NumPy arrays:

>>> func = EinsteinSum(M[i, j] * x[j]).numpify([M, x], [i])  # doctest: +SKIP
>>> func([[1, -1], [0, 1]], [1, 2])  # doctest: +SKIP
[-1. 2.]

Reference
=========

Indexed Objects
---------------

.. module:: sympy.tensor.indexed

.. autoclass:: Indexed
   :members:

.. autoclass:: IndexedBase
   :members:

.. autoclass:: Idx
   :members:

.. autoclass:: DeltaIndexedBase
   :members:

.. module:: sympy.tensor.indexed_sums

.. autofunction:: get_indices

Implicit Sums (Einstein Notation)
---------------------------------

.. _tensor_sums:
.. autoclass:: EinsteinSum
   :members:
