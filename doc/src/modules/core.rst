.. _core_module:

====
Core
====

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
-----
.. module:: sympy.core.cache

cacheit
^^^^^^^
.. autofunction:: __cacheit

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

.. autoclass:: sympy.core.singleton.SingletonRegistry
   :members:

.. autoclass:: Singleton
   :members:

expr
----
.. module:: sympy.core.expr

Expr
^^^^
.. autoclass:: Expr
   :members:

UnevaluatedExpr
^^^^^^^^^^^^^^^
.. autoclass:: UnevaluatedExpr
   :members:

AtomicExpr
^^^^^^^^^^
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

AlgebraicNumber
^^^^^^^^^^^^^^^
.. autoclass:: AlgebraicNumber
   :members:

   .. automethod:: AlgebraicNumber.__new__

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

TribonacciConstant
^^^^^^^^^^^^^^^^^^

.. autoclass:: TribonacciConstant
   :members:

mod_inverse
^^^^^^^^^^^

.. autofunction:: mod_inverse

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

integer_log
^^^^^^^^^^^
.. autofunction:: integer_log

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
.. autoclass:: Relational
   :members:

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
   :private-members:

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
   {f(x), cos(x)}

   If you just want the function you defined, not SymPy functions, the
   thing to search for is AppliedUndef:

   >>> from sympy.core.function import AppliedUndef
   >>> e.atoms(AppliedUndef)
   {f(x)}

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

EvalfMixin
^^^^^^^^^^

.. autoclass:: EvalfMixin
   :members:

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

TupleKind
^^^^^^^^^
.. autoclass:: TupleKind
   :members:

Dict
^^^^
.. autoclass:: Dict
   :members:

exprtools
---------
.. module:: sympy.core.exprtools

gcd_terms
^^^^^^^^^
.. autofunction:: gcd_terms

factor_terms
^^^^^^^^^^^^
.. autofunction:: factor_terms

kind
----
.. module:: sympy.core.kind

Kind
^^^^
.. autoclass:: Kind
   :members:

NumberKind
^^^^^^^^^^
.. autoclass:: NumberKind
   :members:

UndefinedKind
^^^^^^^^^^^^^
.. autoclass:: UndefinedKind
   :members:

BooleanKind
^^^^^^^^^^^
.. autoclass:: BooleanKind
   :members:

Sorting
-------

default_sort_key
^^^^^^^^^^^^^^^^

.. autofunction:: sympy.core.sorting.default_sort_key

ordered
^^^^^^^

.. autofunction:: sympy.core.sorting.ordered

Random
------

.. automodule:: sympy.core.random

random_complex_number
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: random_complex_number

verify_numerically
^^^^^^^^^^^^^^^^^^
.. autofunction:: verify_numerically

test_derivative_numerically
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: test_derivative_numerically

_randrange
^^^^^^^^^^
.. autofunction:: _randrange

_randint
^^^^^^^^
.. autofunction:: _randint

Traversal
---------
.. module:: sympy.core.traversal

bottom_up
^^^^^^^^^
.. autofunction:: bottom_up

postorder_traversal
^^^^^^^^^^^^^^^^^^^
.. autofunction:: postorder_traversal

preorder_traversal
^^^^^^^^^^^^^^^^^^
.. autofunction:: preorder_traversal

use
^^^
.. autofunction:: use

walk
^^^^
.. autofunction:: walk
