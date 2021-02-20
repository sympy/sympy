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
   :undoc-members:
   :private-members:

Atom
^^^^
.. autoclass:: Atom
   :members:
   :undoc-members:
   :private-members:

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
   :undoc-members:
   :private-members:

.. autoclass:: Singleton
   :members:
   :undoc-members:
   :private-members:

expr
----
.. module:: sympy.core.expr

Expr
^^^^
.. autoclass:: Expr
   :members:
   :undoc-members:
   :private-members:

UnevaluatedExpr
^^^^^^^^^^^^^^^
.. autoclass:: UnevaluatedExpr
   :members:
   :undoc-members:
   :private-members:

AtomicExpr
^^^^^^^^^^
.. autoclass:: AtomicExpr
   :members:
   :undoc-members:
   :private-members:

symbol
------
.. module:: sympy.core.symbol

Symbol
^^^^^^
.. autoclass:: Symbol
   :members:
   :undoc-members:
   :private-members:

Wild
^^^^
.. autoclass:: Wild
   :members:
   :undoc-members:
   :private-members:

Dummy
^^^^^
.. autoclass:: Dummy
   :members:
   :undoc-members:
   :private-members:

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
   :undoc-members:
   :private-members:

Float
^^^^^
.. autoclass:: Float
   :members:
   :undoc-members:
   :private-members:

Rational
^^^^^^^^
.. autoclass:: Rational
   :members:
   :undoc-members:
   :private-members:

Integer
^^^^^^^
.. autoclass:: Integer
   :members:
   :undoc-members:
   :private-members:

NumberSymbol
^^^^^^^^^^^^
.. autoclass:: NumberSymbol
   :members:
   :undoc-members:
   :private-members:

RealNumber
^^^^^^^^^^
.. autoclass:: RealNumber
   :members:
   :undoc-members:
   :private-members:

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
   :undoc-members:
   :private-members:

One
^^^

.. autoclass:: One
   :members:
   :undoc-members:
   :private-members:

NegativeOne
^^^^^^^^^^^

.. autoclass:: NegativeOne
   :members:
   :undoc-members:
   :private-members:

Half
^^^^

.. autoclass:: Half
   :members:
   :undoc-members:
   :private-members:

NaN
^^^

.. autoclass:: NaN
   :members:
   :undoc-members:
   :private-members:

Infinity
^^^^^^^^

.. autoclass:: Infinity
   :members:
   :undoc-members:
   :private-members:

NegativeInfinity
^^^^^^^^^^^^^^^^

.. autoclass:: NegativeInfinity
   :members:
   :undoc-members:
   :private-members:

ComplexInfinity
^^^^^^^^^^^^^^^

.. autoclass:: ComplexInfinity
   :members:
   :undoc-members:
   :private-members:

Exp1
^^^^

.. autoclass:: Exp1
   :members:
   :undoc-members:
   :private-members:

ImaginaryUnit
^^^^^^^^^^^^^

.. autoclass:: ImaginaryUnit
   :members:
   :undoc-members:
   :private-members:

Pi
^^

.. autoclass:: Pi
   :members:
   :undoc-members:
   :private-members:

EulerGamma
^^^^^^^^^^

.. autoclass:: EulerGamma
   :members:
   :undoc-members:
   :private-members:

Catalan
^^^^^^^

.. autoclass:: Catalan
   :members:
   :undoc-members:
   :private-members:

GoldenRatio
^^^^^^^^^^^

.. autoclass:: GoldenRatio
   :members:
   :undoc-members:
   :private-members:

TribonacciConstant
^^^^^^^^^^^^^^^^^^

.. autoclass:: TribonacciConstant
   :members:
   :undoc-members:
   :private-members:

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
   :undoc-members:
   :private-members:

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
   :undoc-members:
   :private-members:

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
   :undoc-members:
   :private-members:

mod
---
.. module:: sympy.core.mod

Mod
^^^
.. autoclass:: Mod
   :members:
   :undoc-members:
   :private-members:

relational
----------
.. module:: sympy.core.relational

Rel
^^^
.. autoclass:: Relational
   :members:
   :undoc-members:
   :private-members:

.. autoclass:: Rel
   :members:
   :undoc-members:
   :private-members:

Eq
^^
.. autoclass:: Eq
   :members:
   :undoc-members:
   :private-members:

Ne
^^
.. autoclass:: Ne
   :members:
   :undoc-members:
   :private-members:

Lt
^^
.. autoclass:: Lt
   :members:
   :undoc-members:
   :private-members:

Le
^^
.. autoclass:: Le
   :members:
   :undoc-members:
   :private-members:

Gt
^^
.. autoclass:: Gt
   :members:
   :undoc-members:
   :private-members:

Ge
^^
.. autoclass:: Ge
   :members:
   :undoc-members:
   :private-members:

Equality
^^^^^^^^
.. autoclass:: Equality
   :members:
   :undoc-members:
   :private-members:

GreaterThan
^^^^^^^^^^^
.. autoclass:: GreaterThan
   :members:
   :undoc-members:
   :private-members:

LessThan
^^^^^^^^
.. autoclass:: LessThan
   :members:
   :undoc-members:
   :private-members:

Unequality
^^^^^^^^^^
.. autoclass:: Unequality
   :members:
   :undoc-members:
   :private-members:

StrictGreaterThan
^^^^^^^^^^^^^^^^^
.. autoclass:: StrictGreaterThan
   :members:
   :undoc-members:
   :private-members:

StrictLessThan
^^^^^^^^^^^^^^
.. autoclass:: StrictLessThan
   :members:
   :undoc-members:
   :private-members:

multidimensional
----------------
.. module:: sympy.core.multidimensional

vectorize
^^^^^^^^^
.. autoclass:: vectorize
   :members:
   :undoc-members:
   :private-members:

function
--------
.. module:: sympy.core.function

Lambda
^^^^^^
.. autoclass:: Lambda
   :members:
   :undoc-members:
   :private-members:

WildFunction
^^^^^^^^^^^^
.. autoclass:: WildFunction
   :members:
   :undoc-members:
   :private-members:

Derivative
^^^^^^^^^^
.. autoclass:: Derivative
   :members:
   :undoc-members:
   :private-members:

diff
^^^^
.. autofunction:: diff

FunctionClass
^^^^^^^^^^^^^
.. autoclass:: FunctionClass
   :members:
   :undoc-members:
   :private-members:

Function
^^^^^^^^
.. autoclass:: Function
   :members:
   :undoc-members:
   :private-members:

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
   :undoc-members:
   :private-members:

expand
^^^^^^
.. autofunction:: expand

PoleError
^^^^^^^^^
.. autoclass:: PoleError
   :members:
   :undoc-members:
   :private-members:

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
   :undoc-members:
   :private-members:

PrecisionExhausted
^^^^^^^^^^^^^^^^^^
.. autoclass:: PrecisionExhausted
   :members:
   :undoc-members:
   :private-members:

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
   :undoc-members:
   :private-members:

Dict
^^^^
.. autoclass:: Dict
   :members:
   :undoc-members:
   :private-members:

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

ordered
^^^^^^^

.. autofunction:: sympy.core.compatibility.ordered

kind
----
.. module:: sympy.core.kind

Kind
^^^^
.. autoclass:: Kind
   :members:
   :undoc-members:
   :private-members:

NumberKind
^^^^^^^^^^
.. autoclass:: NumberKind
   :members:
   :undoc-members:
   :private-members:

BooleanKind
^^^^^^^^^^^
.. autoclass:: BooleanKind
   :members:
   :undoc-members:
   :private-members:
