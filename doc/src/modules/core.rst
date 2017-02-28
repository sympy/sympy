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
.. autoclass:: sympy.core.singleton.SingletonRegistry
   :members:

expr
----
.. module:: sympy.core.expr

Expr
----
.. autoclass:: Expr
   :members:

UnevaluatedExpr
---------------

``UnevaluatedExpr()`` is a method provided by Sympy which lets the user to keep an expression unevaluated. By Unevaluated, it is meant that the final value of the expression is not printed and also operations are also avoided with unevaluated expressions. Here are some examples

      >>> from sympy import UnevaluatedExpr
      >>> uexpr = UnevaluatedExpr(S.One*5/7)*UnevaluatedExpr(S.One*3/4)
      >>> uexpr
      (5/7).(3/4)

      >>> from sympy import UnevaluatedExpr
      >>> from sympy.abc import a, b, x, y
      >>> x*UnevaluatedExpr(1/x)
      x.1/x

`` A point to be noted is that ``UnevaluatedExpr()`` can evaluate an expression which is given as argument`` . Example

      >>> from sympy import UnevaluatedExpr
      >>> uexpr = UnevaluatedExpr((2+5)/(3-1) + ((2**3)*3)*4)
      >>> uexpr
      99

`` Another point to be noted is that ``UnevaluatedExpr()`` gives different than what is given by ``Simpify(,evaluate=False)`` . Here are some examples

      >>> from sympy import *
      >>> x*sympify('1/x',evaluate=False)
      1

      >>> from sympy import *
      >>> uexpr = UnevaluatedExpr((2+5)/(3-1) + ((2**3)*3)*4)
      >>> simexpr = sympify('(2+5)/(3-1) + ((2**3)*3)*4', evaluate=False)
      >>> uexpr
      99
      >>> simexpr
       2+5        3
      ──── + 3.4.2
        -1+3

``UnevalutedExpr()`` sre supported by Sympy printers and can be used for printing the result in different forms. Example

      >>> from sympy import UnevaluatedExpr
      >>> uexpr = UnevaluatedExpr(S.One*5/7)*UnevaluatedExpr(S.One*3/4)
      >>> print(latex(uexpr))
      \frac{5}{7} \frac{3}{4}

      
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
