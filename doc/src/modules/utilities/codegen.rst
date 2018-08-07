.. _codegen_API:

=======
Codegen
=======

This module provides functionality to generate directly compilable code from
SymPy expressions.  The ``codegen`` function is the user interface to the code
generation functionality in SymPy.  Some details of the implementation is given
below for advanced users that may want to use the framework directly.

.. note:: The ``codegen`` callable is not in the sympy namespace automatically,
   to use it you must first execute

   >>> from sympy.utilities.codegen import codegen

Implementation Details
======================

Here we present the most important pieces of the internal structure, as
advanced users may want to use it directly, for instance by subclassing a code
generator for a specialized application.  **It is very likely that you would
prefer to use the codegen() function documented above.**

Basic assumptions:

* A generic Routine data structure describes the routine that must be translated
  into C/Fortran/... code. This data structure covers all features present in
  one or more of the supported languages.

* Descendants from the CodeGen class transform multiple Routine instances into
  compilable code. Each derived class translates into a specific language.

* In many cases, one wants a simple workflow. The friendly functions in the last
  part are a simple api on top of the Routine/CodeGen stuff. They are easier to
  use, but are less powerful.

Routine
=======

The Routine class is a very important piece of the codegen module. Viewing the
codegen utility as a translator of mathematical expressions into a set of
statements in a programming language, the Routine instances are responsible for
extracting and storing information about how the math can be encapsulated in a
function call.  Thus, it is the Routine constructor that decides what arguments
the routine will need and if there should be a return value.

API Reference
=============

.. automodule:: sympy.utilities.codegen
   :members:
