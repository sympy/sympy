.. _assumptions_module:

===========
Assumptions
===========

.. automodule:: sympy.assumptions


Predicate
=========

.. autoclass:: sympy.assumptions.assume::Predicate
   :members:
   :noindex:

.. autoclass:: sympy.assumptions.assume::AppliedPredicate
   :members:
   :noindex:


Querying
========

Queries are used to ask information about expressions. Main method for this
is ``ask()``:

.. autofunction:: sympy.assumptions.ask::ask
   :noindex:

``ask``'s optional second argument should be a boolean expression involving
assumptions about objects in *expr*. Valid values include:

    * ``Q.integer(x)``
    * ``Q.positive(x)``
    * ``Q.integer(x) & Q.positive(x)``
    * etc.

``Q`` is an object holding known predicates.

See documentation for the logic module for a complete list of valid boolean
expressions.

You can also define a context so you don't have to pass that argument
each time to function ``ask()``. This is done by using the assuming context manager
from module sympy.assumptions. ::

     >>> from sympy import *
     >>> x = Symbol('x')
     >>> y = Symbol('y')
     >>> facts = Q.positive(x), Q.positive(y)
     >>> with assuming(*facts):
     ...     print(ask(Q.positive(2*x + y)))
     True


Contents
========

.. toctree::
   :titlesonly:

   ask.rst
   assume.rst
   refine.rst
   predicates.rst


Performance improvements
========================

On queries that involve symbolic coefficients, logical inference is used. Work on
improving satisfiable function (sympy.logic.inference.satisfiable) should result
in notable speed improvements.

Logic inference used in one ask could be used to speed up further queries, and
current system does not take advantage of this. For example, a truth maintenance
system (https://en.wikipedia.org/wiki/Truth_maintenance_system) could be implemented.

Misc
====

You can find more examples in the form of tests in the directory
``sympy/assumptions/tests/``
