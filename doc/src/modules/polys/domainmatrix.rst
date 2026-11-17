.. _polys-domainmatrix:

===============================================
Introducing the domainmatrix of the poly module
===============================================

This page introduces the idea behind domainmatrix which is used in SymPy's
:mod:`sympy.polys` module. This is a relatively advanced topic so for a better understanding
it is recommended to read about :py:class:`~.Domain` and :py:class:`~.DDM` along with
:mod:`sympy.matrices` module.

What is domainmatrix?
=====================

It is way of associating Matrix with :py:class:`~.Domain`.

A domainmatrix represents a matrix with elements that are in a particular
Domain. Each domainmatrix internally wraps a DDM which is used for the lower-level operations.
The idea is that the domainmatrix class provides the convenience routines for converting
between Expr and the poly domains as well as unifying matrices with different domains.

In general, we represent a matrix without concerning about the :py:class:`~.Domain` as:
   >>> from sympy import Matrix
   >>> from sympy.polys.matrices import DomainMatrix
   >>> A = Matrix([
   ... [1, 2],
   ... [3, 4]])
   >>> A
   Matrix([
   [1, 2],
   [3, 4]])

.. currentmodule:: sympy.polys.matrices

.. automodule:: sympy.polys.matrices.domainmatrix
   :members:

.. automodule:: sympy.polys.matrices.ddm
   :members:

.. automodule:: sympy.polys.matrices.dense
   :members:

.. automodule:: sympy.polys.matrices._typing
   :members:

.. automodule:: sympy.polys.matrices.sdm
   :members:

.. automodule:: sympy.polys.matrices._dfm
   :members:

.. currentmodule:: sympy.polys.matrices.normalforms

.. autofunction:: smith_normal_form
.. autofunction:: hermite_normal_form
