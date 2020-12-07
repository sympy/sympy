.. _assumptions_module:

===========
Assumptions
===========

.. automodule:: sympy.assumptions

Contents
========

.. toctree::
    :maxdepth: 3

    ask.rst
    assume.rst
    refine.rst
    handlers/index.rst


Predicates
==========


.. autoclass:: sympy.assumptions.assume::Predicate
   :noindex:

.. autoclass:: sympy.assumptions.assume::AppliedPredicate
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


Design
======

Each time ``ask()`` is called, the appropriate handler for the current predicate is called.
The handler is always an instance of ``sympy.multipledispatch.Dispatcher``.
Class signature for its arguments is registered to the handler. For example, a handler for the key
'positive' would look like this::

    AskPositiveHandler = Dispatcher('positive')

    @AskPositiveHandler.register(Mul)
    def Mul_handler(expr, assumptions):
        # check the arguments of expr and return True, False or None accordingly
        ...

    @AskPositiveHandler.register(Add)
    def Add_handler(expr, assumptions):
        ...

The function ``Mul_handler()`` is called when ``expr`` is an instance of :class:`~.Mul()`, and ``Add_handler()``
is called when ``expr`` is an instance of :class:`~.Add()`.


Extensibility
=============

You can define new queries or support new types by subclassing ``Predicate`` and registering the instance
to ``Q``. Supporting new types for the handler is done by dispatching.

In the following example, we will create a new predicate which checks if the argument is mersenne number [1]_.

.. parsed-literal::

    >>> from sympy.assumptions import Predicate, Q
    >>> class MersennePredicate(Predicate):
    ...     """Return True if argument is (2**n)-1 pattern."""
    >>> Q.mersenne = MersennePredicate("mersenne")
    >>> @Q.mersenne.register(Integer)
    ... def _(expr, assumptions):
    ...     from sympy import log
    ...     if ask(Q.integer(log(expr + 1, 2))):
    ...             return True
    >>> ask(Q.mersenne(7))
    True

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

You can find more examples in the in the form of test under directory
sympy/assumptions/tests/

.. [1] https://en.wikipedia.org/wiki/Mersenne_prime
