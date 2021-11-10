.. _reference:

====================
SymPy API Reference
====================

.. module:: sympy

This page has a detailed description of the SymPy API, categorized as ``public`` and
``internal`` for public and private modules respectively. It discusses 
what the different functions and modules do, their arguments and return values.
For complete documentation including how to get started, 
see `documentation overview <https://docs.sympy.org>`_.

The public API documentation is further subdivided into eight categories for easy
navigation. The documentation for the public API is manually implemented using 
``autoapi`` directives while ``autosummary`` is used for auto generation of
the internal API documentation.

Public API Reference
=====================

This section contains a summary of the public SymPy modules, that is to say
functions and methods. All functions and objects implemented in the ``sympy``
core subpackage are documented below.

..  rst-class:: column

:ref:`Basics <basic_modules>`
-------------------------------------------------------

   Contains method docstrings for the basic modules. Subcategories include: *absolute basics*, 
   *manipulation*, *assumptions*, *functions*, *simplification*, *calculus*, *solvers*, and some 
   other subcategories.


..  rst-class:: column

:ref:`Code Generation <codegen_module>`
-------------------------------------------------------

   Contains a description of methods for the generation of compilable and executable code.


..  rst-class:: column

:ref:`Logic <logic>`
-------------------------------------------------------

   Contains method docstrings for the *logic* and *sets* modules.

..  rst-class:: column

:ref:`Matrices <matrices_modules>`
-------------------------------------------------------

   Contains method docstrings for the matrices, tensor and vector modules.

..  rst-class:: column

:ref:`Number Theory <numtheory_module>`
-------------------------------------------------------

   Contains doc strings for the Number theory methods.

..  rst-class:: column

:ref:`Physics <physics-docs>`
-------------------------------------------------------

   Contains doc strings for Physics methods.

..  rst-class:: column

:ref:`Utilities <utilities>`
-------------------------------------------------------

   Contains docstrings for methods of several utility modules. Subcategories
   include: *Interactive*, *Parsing*, *Printing*, *Testing*, *Utilities*.


..  rst-class:: column

:ref:`Topics <topics>`
-------------------------------------------------------

   Contains method docstrings for several modules. Subcategories include : *Plotting*, 
   *Polynomials*, *Geometry*, *Category Theory*, *Cryptography*, *Differential*, *Holonomic*, 
   *Lie Algebra*, and *Stats*.


.. toctree::
   :hidden:
   :maxdepth: 2

   public/basics/index.rst
   public/codegeneration/index.rst
   public/logic/index.rst
   public/matrices/index.rst
   public/numbertheory/index.rst
   public/physics/index.rst
   public/utilities/index.rst
   public/topics/index.rst
