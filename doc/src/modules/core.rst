SymPy Core
==========

sympify
-------
.. module:: sympy.core.sympify

sympify
^^^^^^^
.. autofunction:: sympify

assumptions
-----------

.. automodule:: sympy.core.assumptions

cache
-------
.. module:: sympy.core.cache

cacheit
^^^^^^^
.. autofunction:: cacheit

basic
-----
.. module:: sympy.core.basic

Basic
^^^^^
.. autoclass:: Basic
   :members:

Atom
^^^^
.. autoclass:: Atom
   :members:

core
----
.. module:: sympy.core.core

singleton
---------
.. module:: sympy.core.singleton

S
^
.. autoclass:: Singleton
   :members:

.. note:: Some important points.

   1.Normally the first invocation of this class creates the unique instance of the class and
   subsequent invocations would simply return a reference to the instance created earlier.
   ``S.One`` denotes a singleton object in SymPy having value ``Integer(1)``.The ``S`` object makes it so
   that only one instance of ``Integer(1)`` is ever created. The reason that these are stored in this way
   is because of their frequent use throughout the codebase, which makes it more efficient.

   2. Most Integers -- even identical ones like ``Integer(1)`` and ``Integer(1)`` -- get stored
   at arbitrary memory locations and tests for equality between them and a value must be
   done with the ``==`` operator: ``x == 1``. By having a single instance, SymPy enables two optimizations.
   First, it saves memory. Second, you can compare against these objects using is comparison, like ``x is S.One``. Because only one instance can ever exist ``Integer(1)`` is always the same as ``S.One``.

   * There are several objects in SymPy which are implemented as singletons. Some of them are ``S.NegativeOne``, ``S.Zero``, ``S.NegativeInfinity``, ``S.EmptySet``, ``S.Infinity`` and others can be seen by running ``dir(S)``. Ignore the ones that start with ``__``  those are Python internal methods.


expr
----
.. module:: sympy.core.expr

Expr
----
.. autoclass:: Expr
   :members:

AtomicExpr
----------
.. autoclass:: AtomicExpr
   :members:

symbol
------
.. module:: sympy.core.symbol

Symbol
^^^^^^
.. autoclass:: Symbol
   :members:

Wild
^^^^
.. autoclass:: Wild
   :members:

Dummy
^^^^^
.. autoclass:: Dummy
   :members:

symbols
^^^^^^^
.. autofunction:: symbols

var
^^^
.. autofunction:: var

numbers
-------
.. module:: sympy.core.numbers

Number
^^^^^^
.. autoclass:: Number
   :members:

Float
^^^^^
.. autoclass:: Float
   :members:

Rational
^^^^^^^^
.. autoclass:: Rational
   :members:

Integer
^^^^^^^
.. autoclass:: Integer
   :members:

NumberSymbol
^^^^^^^^^^^^
.. autoclass:: NumberSymbol
   :members:

RealNumber
^^^^^^^^^^
.. autoclass:: RealNumber
   :members:

igcd
^^^^
.. autofunction:: igcd

ilcm
^^^^
.. autofunction:: ilcm

seterr
^^^^^^
.. autofunction:: seterr

Zero
^^^^

.. autoclass:: Zero
   :members:

One
^^^

.. autoclass:: One
   :members:

NegativeOne
^^^^^^^^^^^

.. autoclass:: NegativeOne
   :members:

Half
^^^^

.. autoclass:: Half
   :members:

NaN
^^^

.. autoclass:: NaN
   :members:

Infinity
^^^^^^^^

.. autoclass:: Infinity
   :members:

NegativeInfinity
^^^^^^^^^^^^^^^^

.. autoclass:: NegativeInfinity
   :members:

ComplexInfinity
^^^^^^^^^^^^^^^

.. autoclass:: ComplexInfinity
   :members:

Exp1
^^^^

.. autoclass:: Exp1
   :members:

ImaginaryUnit
^^^^^^^^^^^^^

.. autoclass:: ImaginaryUnit
   :members:

Pi
^^

.. autoclass:: Pi
   :members:

EulerGamma
^^^^^^^^^^

.. autoclass:: EulerGamma
   :members:

Catalan
^^^^^^^

.. autoclass:: Catalan
   :members:

GoldenRatio
^^^^^^^^^^^

.. autoclass:: GoldenRatio
   :members:

power
-----
.. module:: sympy.core.power

Pow
^^^
.. autoclass:: Pow
   :members:

integer_nthroot
^^^^^^^^^^^^^^^
.. autofunction:: integer_nthroot

mul
---
.. module:: sympy.core.mul

Mul
^^^
.. autoclass:: Mul
   :members:

prod
^^^^
.. autofunction:: prod

add
---
.. module:: sympy.core.add

Add
^^^
.. autoclass:: Add
   :members:

mod
---
.. module:: sympy.core.mod

Mod
^^^
.. autoclass:: Mod
   :members:

relational
----------
.. module:: sympy.core.relational

Rel
^^^
.. autoclass:: Rel
   :members:

Eq
^^
.. autoclass:: Eq
   :members:

Ne
^^
.. autoclass:: Ne
   :members:

Lt
^^
.. autoclass:: Lt
   :members:

Le
^^
.. autoclass:: Le
   :members:

Gt
^^
.. autoclass:: Gt
   :members:

Ge
^^
.. autoclass:: Ge
   :members:

Equality
^^^^^^^^
.. autoclass:: Equality
   :members:

GreaterThan
^^^^^^^^^^^
.. autoclass:: GreaterThan
   :members:

LessThan
^^^^^^^^
.. autoclass:: LessThan
   :members:

Unequality
^^^^^^^^^^
.. autoclass:: Unequality
   :members:

StrictGreaterThan
^^^^^^^^^^^^^^^^^
.. autoclass:: StrictGreaterThan
   :members:

StrictLessThan
^^^^^^^^^^^^^^
.. autoclass:: StrictLessThan
   :members:

multidimensional
----------------
.. module:: sympy.core.multidimensional

vectorize
^^^^^^^^^
.. autoclass:: vectorize
   :members:

function
--------
.. module:: sympy.core.function

Lambda
^^^^^^
.. autoclass:: Lambda
   :members:

WildFunction
^^^^^^^^^^^^
.. autoclass:: WildFunction
   :members:

Derivative
^^^^^^^^^^
.. autoclass:: Derivative
   :members:

diff
^^^^
.. autofunction:: diff

FunctionClass
^^^^^^^^^^^^^
.. autoclass:: FunctionClass
   :members:

Function
^^^^^^^^
.. autoclass:: Function
   :members:

.. note:: Not all functions are the same

   SymPy defines many functions (like ``cos`` and ``factorial``). It also
   allows the user to create generic functions which act as argument
   holders. Such functions are created just like symbols:

   >>> from sympy import Function, cos
   >>> from sympy.abc import x
   >>> f = Function('f')
   >>> f(2) + f(x)
   f(2) + f(x)

   If you want to see which functions appear in an expression you can use
   the atoms method:

   >>> e = (f(x) + cos(x) + 2)
   >>> e.atoms(Function)
   set([f(x), cos(x)])

   If you just want the function you defined, not SymPy functions, the
   thing to search for is AppliedUndef:

   >>> from sympy.core.function import AppliedUndef
   >>> e.atoms(AppliedUndef)
   set([f(x)])

Subs
^^^^
.. autoclass:: Subs
   :members:

expand
^^^^^^
.. autofunction:: expand

PoleError
^^^^^^^^^
.. autoclass:: PoleError
   :members:

count_ops
^^^^^^^^^
.. autofunction:: count_ops

expand_mul
^^^^^^^^^^
.. autofunction:: expand_mul

expand_log
^^^^^^^^^^
.. autofunction:: expand_log

expand_func
^^^^^^^^^^^
.. autofunction:: expand_func

expand_trig
^^^^^^^^^^^
.. autofunction:: expand_trig

expand_complex
^^^^^^^^^^^^^^
.. autofunction:: expand_complex

expand_multinomial
^^^^^^^^^^^^^^^^^^
.. autofunction:: expand_multinomial

expand_power_exp
^^^^^^^^^^^^^^^^
.. autofunction:: expand_power_exp

expand_power_base
^^^^^^^^^^^^^^^^^
.. autofunction:: expand_power_base

nfloat
^^^^^^
.. autofunction:: nfloat

evalf
-----
.. module:: sympy.core.evalf

PrecisionExhausted
^^^^^^^^^^^^^^^^^^
.. autoclass:: PrecisionExhausted
   :members:

N
^
.. autofunction:: N

containers
----------
.. module:: sympy.core.containers

Tuple
^^^^^
.. autoclass:: Tuple
   :members:

Dict
^^^^
.. autoclass:: Dict
   :members:

compatibility
-------------
.. module:: sympy.core.compatibility

iterable
^^^^^^^^
.. autofunction:: iterable

is_sequence
^^^^^^^^^^^
.. autofunction:: is_sequence

as_int
^^^^^^
.. autofunction:: as_int

exprtools
---------
.. module:: sympy.core.exprtools

gcd_terms
^^^^^^^^^
.. autofunction:: gcd_terms

factor_terms
^^^^^^^^^^^^
.. autofunction:: factor_terms
