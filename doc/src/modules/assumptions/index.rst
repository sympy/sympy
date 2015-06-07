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

Q is a class in sympy.assumptions holding known predicates.

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

Supported predicates
====================

bounded
-------

Test that a function is bounded with respect to its variables. For example,
sin(x) is a bounded functions, but exp(x) is not.

Examples::

    >>> from sympy import *
    >>> x = Symbol('x')
    >>> ask(Q.finite(exp(x)), ~Q.finite(x))
    False
    >>> ask(Q.finite(exp(x)) , Q.finite(x))
    True
    >>> ask(Q.finite(sin(x)), ~Q.finite(x))
    True


commutative
-----------

Test that objects are commutative. By default, symbols in SymPy are considered
commutative except otherwise stated.

Examples::

    >>> from sympy import *
    >>> x, y = symbols('x,y')
    >>> ask(Q.commutative(x))
    True
    >>> ask(Q.commutative(x), ~Q.commutative(x))
    False
    >>> ask(Q.commutative(x*y), ~Q.commutative(x))
    False


complex
-------

Test that expression belongs to the field of complex numbers.

Examples::

    >>> from sympy import *
    >>> ask(Q.complex(2))
    True
    >>> ask(Q.complex(I))
    True
    >>> x, y = symbols('x,y')
    >>> ask(Q.complex(x+I*y), Q.real(x) & Q.real(y))
    True


even
----

Test that expression represents an even number, that is, an number that
can be written in the form 2*n, n integer.

Examples::

    >>> from sympy import *
    >>> ask(Q.even(2))
    True
    >>> n = Symbol('n')
    >>> ask(Q.even(2*n), Q.integer(n))
    True


extended_real
-------------

Test that an expression belongs to the field of extended real numbers, that is, real
numbers union {Infinity, -Infinity}.

Examples::

    >>> from sympy import *
    >>> ask(Q.extended_real(oo))
    True
    >>> ask(Q.extended_real(2))
    True
    >>> ask(Q.extended_real(x), Q.real(x))
    True


imaginary
---------

Test that an expression belongs to the set of imaginary numbers, that is,
 it can be written as x*I, where x is real and I is the imaginary unit.

Examples::

    >>> from sympy import *
    >>> ask(Q.imaginary(2*I))
    True
    >>> x = Symbol('x')
    >>> ask(Q.imaginary(x*I), Q.real(x))
    True


infinitesimal
-------------

Test that an expression is equivalent to an infinitesimal number.

Examples::

    >>> from sympy import *
    >>> ask(Q.infinitesimal(1/oo))
    True
    >>> x, y = symbols('x,y')
    >>> ask(Q.infinitesimal(2*x), Q.infinitesimal(x))
    True
    >>> ask(Q.infinitesimal(x*y), Q.infinitesimal(x) & Q.finite(y))
    True


integer
-------

Test that an expression belongs to the set of integer numbers.

Examples::

    >>> from sympy import *
    >>> ask(Q.integer(2))
    True
    >>> ask(Q.integer(sqrt(2)))
    False
    >>> x = Symbol('x')
    >>> ask(Q.integer(x/2), Q.even(x))
    True


irrational
----------

Test that an expression represents an irrational number.

Examples::

     >>> from sympy import *
     >>> ask(Q.irrational(pi))
     True
     >>> ask(Q.irrational(sqrt(2)))
     True
     >>> ask(Q.irrational(x*sqrt(2)), Q.rational(x))
     True


rational
--------

Test that an expression represents a rational number.

Examples::

     >>> from sympy import *
     >>> ask(Q.rational(Rational(3, 4)))
     True
     >>> x, y = symbols('x,y')
     >>> ask(Q.rational(x/2), Q.integer(x))
     True
     >>> ask(Q.rational(x/y), Q.integer(x) & Q.integer(y))
     True


negative
--------

Test that an expression is less (strict) than zero.

Examples::

     >>> from sympy import *
     >>> ask(Q.negative(0.3))
     False
     >>> x = Symbol('x')
     >>> ask(Q.negative(-x), Q.positive(x))
     True

Remarks
^^^^^^^
negative numbers are defined as real numbers that are not zero nor positive, so
complex numbers (with nontrivial imaginary coefficients) will return False
for this predicate. The same applies to Q.positive.


positive
--------

Test that a given expression is greater (strict) than zero.

Examples::

     >>> from sympy import *
     >>> ask(Q.positive(0.3))
     True
     >>> x = Symbol('x')
     >>> ask(Q.positive(-x), Q.negative(x))
     True

Remarks
^^^^^^^
see Remarks for negative


prime
-----

Test that an expression represents a prime number.

Examples::

    >>> from sympy import *
    >>> ask(Q.prime(13))
    True

Remarks: Use sympy.ntheory.isprime to test numeric values efficiently.


real
----

Test that an expression belongs to the field of real numbers.

Examples::

    >>> from sympy import *
    >>> ask(Q.real(sqrt(2)))
    True
    >>> x, y = symbols('x,y')
    >>> ask(Q.real(x*y), Q.real(x) & Q.real(y))
    True


odd
---

Test that an expression represents an odd number.

Examples::

    >>> from sympy import *
    >>> ask(Q.odd(3))
    True
    >>> n = Symbol('n')
    >>> ask(Q.odd(2*n + 1), Q.integer(n))
    True


nonzero
-------

Test that an expression is not zero.

Examples::

     >>> from sympy import *
     >>> x = Symbol('x')
     >>> ask(Q.nonzero(x), Q.positive(x) | Q.negative(x))
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
