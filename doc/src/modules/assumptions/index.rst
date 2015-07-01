==================
Assumptions module
==================

.. automodule:: sympy.assumptions

Contents
========

.. toctree::
    :maxdepth: 3

    ask.rst
    assume.rst
    refine.rst
    handlers/index.rst

Queries are used to ask information about expressions. Main method for this
is ask():

.. autofunction:: sympy.assumptions.ask.ask
   :noindex:

Querying
========

ask's optional second argument should be a boolean expression involving
assumptions about objects in expr. Valid values include:

    * Q.integer(x)
    * Q.positive(x)
    * Q.integer(x) & Q.positive(x)
    * etc.

Q is an object holding known predicates.

See documentation for the logic module for a complete list of valid boolean
expressions.

You can also define a context so you don't have to pass that argument
each time to function ask(). This is done by using the assuming context manager
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

Each time ask is called, the appropriate Handler for the current key is called. This is
always a subclass of sympy.assumptions.AskHandler. It's classmethods have the name's of the classes
it supports. For example, a (simplified) AskHandler for the ask 'positive' would
look like this::

    class AskPositiveHandler(CommonHandler):

        def Mul(self):
            # return True if all argument's in self.expr.args are positive
            ...

        def Add(self):
            for arg in self.expr.args:
                if not ask(arg, positive, self.assumptions):
                    break
            else:
                # if all argument's are positive
                return True
        ...

The .Mul() method is called when self.expr is an instance of Mul, the Add method
would be called when self.expr is an instance of Add and so on.


Extensibility
=============

You can define new queries or support new types by subclassing sympy.assumptions.AskHandler
 and registering that handler for a particular key by calling register_handler:

.. autofunction:: sympy.assumptions.ask.register_handler
                  :noindex:

You can undo this operation by calling remove_handler.

.. autofunction:: sympy.assumptions.ask.remove_handler
                  :noindex:

You can support new types [1]_ by adding a handler to an existing key. In the
following example, we will create a new type MyType and extend the key 'prime'
to accept this type (and return True)

.. parsed-literal::

    >>> from sympy.core import Basic
    >>> from sympy.assumptions import register_handler
    >>> from sympy.assumptions.handlers import AskHandler
    >>> class MyType(Basic):
    ...     pass
    >>> class MyAskHandler(AskHandler):
    ...     @staticmethod
    ...     def MyType(expr, assumptions):
    ...         return True
    >>> a = MyType()
    >>> register_handler('prime', MyAskHandler)
    >>> ask(Q.prime(a))
    True


Performance improvements
========================

On queries that involve symbolic coefficients, logical inference is used. Work on
improving satisfiable function (sympy.logic.inference.satisfiable) should result
in notable speed improvements.

Logic inference used in one ask could be used to speed up further queries, and
current system does not take advantage of this. For example, a truth maintenance
system (http://en.wikipedia.org/wiki/Truth_maintenance_system) could be implemented.

Misc
====

You can find more examples in the in the form of test under directory
sympy/assumptions/tests/

.. [1] New type must inherit from Basic, otherwise an exception will be raised.
   This is a bug and should be fixed.
