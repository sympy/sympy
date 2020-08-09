===
Map
===

Introduction
------------

This module implements mathematical maps as SymPy object.

Basic Map
---------

.. module:: sympy.map.map

Set of functions
^^^^^^^^^^^^^^^^

.. autoclass:: FunctionSet
   :members:

For general use, ``function_set`` is provided which is a set of every functions.

Maps
^^^^

.. autoclass:: Map
   :members:

.. autoclass:: UndefinedMap
   :members:

.. autoclass:: RestrictedMap
   :members:

.. autoclass:: InverseMap
   :members:

.. autoclass:: IdentityMap
   :members:

.. autoclass:: ConstantMap
   :members:

Applied map
^^^^^^^^^^^

.. autoclass:: AppliedMap
   :members:

.. autofunction:: isappliedmap

Binary operators
----------------

.. module:: sympy.map.operator

Operators
^^^^^^^^^

.. autoclass:: BinaryOperator
   :members:

.. autoclass:: LeftDivisionOperator
   :members:

.. autoclass:: RightDivisionOperator
   :members:

Derived operators
^^^^^^^^^^^^^^^^^

.. autoclass:: InverseOperator
   :members:

.. autoclass:: ExponentOperator
   :members:

Inverse element and exponent element
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: InverseElement
   :members:

.. autoclass:: ExponentElement
   :members:

Addition operators
------------------

.. module:: sympy.map.add

Operators for scalar and vector addition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: AdditionOperator
   :members:

.. autoclass:: NumericAdditionOperator
   :members:

.. autoclass:: VectorAdditionOperator
   :members:

Result of addition
^^^^^^^^^^^^^^^^^^

.. autoclass:: Addition
   :members:

Multiplication operators
------------------------

.. module:: sympy.map.mul

Operators for scalar and vector multiplication
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: MultiplicationOperator
   :members:

.. autoclass:: NumericMultiplicationOperator
   :members:

.. autoclass:: ScalarMultiplicationOperator
   :members:

.. autoclass:: VectorMultiplicationOperator
   :members:

Result of multiplication
^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: Multiplication
   :members:

Function composition operators
------------------------------

.. module:: sympy.map.composite

Operators for function composition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: CompositionOperator
   :members:

For general use, ``composite_op`` is provided which is
a composition operator between any function.

.. autoclass:: IterationOperator
   :members:

Result of function composition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: CompositeMap
   :members:

.. autoclass:: IteratedMap
   :members:
