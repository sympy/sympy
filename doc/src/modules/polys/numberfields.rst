.. _polys-numberfields:

=============
Number Fields
=============

Support for algebraic number theory, including calculation of
minimal polynomials, and solution of
:ref:`the subfield problem<SubfieldProblem>` and allied problems.

Algebraic Numbers
=================

.. currentmodule:: sympy.polys.numberfields
.. autoclass:: AlgebraicNumber
   :members:


Minimal Polynomials
===================

.. currentmodule:: sympy.polys.numberfields.minpoly

.. autofunction:: minimal_polynomial
.. autofunction:: minpoly


The Subfield Problem
====================
.. _SubfieldProblem:

.. automodule:: sympy.polys.numberfields.subfield

.. autofunction:: field_isomorphism

.. autofunction:: primitive_element

.. autofunction:: to_number_field
