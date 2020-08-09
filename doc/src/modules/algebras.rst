========
Algebras
========

Introduction
------------

The Algebras module for SymPy provides support for abstract algebra, and
basic algebraic operations on Quaternions. Interface between these two is
not implemented yet.

Abstract algebra
----------------

.. module:: sympy.algebras.abstract

Basic algebraic structure
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: AlgebraicStructure
   :members:

Group-like structures
^^^^^^^^^^^^^^^^^^^^^

.. module:: sympy.algebras.abstract.group

.. autoclass:: Magma
   :members:

.. autoclass:: Semigroup
   :members:

.. autoclass:: Monoid
   :members:

.. autoclass:: Quasigroup
   :members:

.. autoclass:: Loop
   :members:

.. autoclass:: Group
   :members:

Ring-like structures
^^^^^^^^^^^^^^^^^^^^

.. module:: sympy.algebras.abstract.ring

.. autoclass:: Ring
   :members:

.. autoclass:: Field
   :members:

Ring-like structures
^^^^^^^^^^^^^^^^^^^^

.. module:: sympy.algebras.abstract.module

.. autoclass:: Module
   :members:

.. autoclass:: VectorSpace
   :members:

Algebra-like structures
^^^^^^^^^^^^^^^^^^^^^^^

.. module:: sympy.algebras.abstract.algebra

.. autoclass:: Algebra
   :members:

Operators to construct function space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. module:: sympy.algebras.abstract.module.functionspace

.. autoclass:: FunctionAdditionOperator
   :members:

.. autoclass:: FunctionAddition
   :members:

.. autoclass:: FunctionScalarMultiplicationOperator
   :members:

.. autoclass:: FunctionVectorMultiplicationOperator
   :members:

.. autoclass:: FunctionMultiplication
   :members:

.. autoclass:: FunctionExponent
   :members:

.. autoclass:: ReciprocalFunction
   :members:

Quaternion Reference
--------------------

.. automodule:: sympy.algebras

This section lists the classes implemented by the Algebras module.

.. autoclass:: Quaternion
   :members:
