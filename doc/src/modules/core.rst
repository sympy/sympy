.. _core_module:

====
Core
====

sympify
-------
.. module:: sympy.core.sympify

.. autofunction:: sympify

assumptions
-----------

.. automodule:: sympy.core.assumptions

cache
-----
.. module:: sympy.core.cache

.. autofunction:: __cacheit

basic
-----
.. module:: sympy.core.basic

.. autoclass:: Basic
   :members:

.. autoclass:: Atom
   :members:

singleton
---------
.. module:: sympy.core.singleton


.. autoclass:: sympy.core.singleton.SingletonRegistry
   :members:

.. autoclass:: Singleton
   :members:

expr
----
.. module:: sympy.core.expr

.. autoclass:: Expr
   :members:

.. autoclass:: UnevaluatedExpr
   :members:

.. autoclass:: AtomicExpr
   :members:

symbol
------
.. module:: sympy.core.symbol

.. autoclass:: Symbol
   :members:

.. autoclass:: Wild
   :members:

.. autoclass:: Dummy
   :members:

.. autofunction:: symbols

.. autofunction:: var

intfunc
-------
.. module:: sympy.core.intfunc

.. autofunction:: num_digits

.. autofunction:: trailing

.. autofunction:: ilcm

.. autofunction:: igcd

.. autofunction:: igcd_lehmer

.. autofunction:: igcdex

.. autofunction:: isqrt

.. autofunction:: integer_nthroot

.. autofunction:: integer_log

.. autofunction:: mod_inverse

numbers
-------
.. module:: sympy.core.numbers

.. autoclass:: Number
   :members:

.. autoclass:: Float
   :members:

.. autoclass:: Rational
   :members:

.. autoclass:: Integer
   :members:

.. autoclass:: AlgebraicNumber
   :members:

   .. automethod:: AlgebraicNumber.__new__

.. autoclass:: NumberSymbol
   :members:

.. autoclass:: RealNumber
   :members:

.. autofunction:: seterr


.. autoclass:: Zero
   :members:


.. autoclass:: One
   :members:


.. autoclass:: NegativeOne
   :members:


.. autoclass:: Half
   :members:


.. autoclass:: NaN
   :members:


.. autoclass:: Infinity
   :members:


.. autoclass:: NegativeInfinity
   :members:


.. autoclass:: ComplexInfinity
   :members:


.. autoclass:: Exp1
   :members:


.. autoclass:: ImaginaryUnit
   :members:


.. autoclass:: Pi
   :members:


.. autoclass:: EulerGamma
   :members:


.. autoclass:: Catalan
   :members:


.. autoclass:: GoldenRatio
   :members:


.. autoclass:: TribonacciConstant
   :members:


.. autofunction:: mod_inverse

.. autofunction:: equal_valued


power
-----
.. module:: sympy.core.power

.. autoclass:: Pow
   :members:

mul
---
.. module:: sympy.core.mul

.. autoclass:: Mul
   :members:

.. autofunction:: prod

add
---
.. module:: sympy.core.add

.. autoclass:: Add
   :members:

mod
---
.. module:: sympy.core.mod

.. autoclass:: Mod
   :members:

relational
----------
.. module:: sympy.core.relational

.. autoclass:: Relational
   :members:

.. autoclass:: Rel
   :members:

.. autoclass:: Eq
   :members:

.. autoclass:: Ne
   :members:

.. autoclass:: Lt
   :members:

.. autoclass:: Le
   :members:

.. autoclass:: Gt
   :members:

.. autoclass:: Ge
   :members:

.. autoclass:: Equality
   :members:

.. autoclass:: GreaterThan
   :members:

.. autoclass:: LessThan
   :members:

.. autoclass:: Unequality
   :members:

.. autoclass:: StrictGreaterThan
   :members:

.. autoclass:: StrictLessThan
   :members:

multidimensional
----------------
.. module:: sympy.core.multidimensional

.. autoclass:: vectorize
   :members:

function
--------
.. module:: sympy.core.function

.. autoclass:: Lambda
   :members:

.. autoclass:: WildFunction
   :members:

.. autoclass:: Derivative
   :members:
   :private-members:

.. autofunction:: diff

.. autoclass:: FunctionClass
   :members:

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

.. autoclass:: Subs
   :members:

.. autofunction:: expand

.. autoclass:: PoleError
   :members:

.. autofunction:: count_ops

.. autofunction:: expand_mul

.. autofunction:: expand_log

.. autofunction:: expand_func

.. autofunction:: expand_trig

.. autofunction:: expand_complex

.. autofunction:: expand_multinomial

.. autofunction:: expand_power_exp

.. autofunction:: expand_power_base

.. autofunction:: nfloat

evalf
-----
.. module:: sympy.core.evalf


.. autoclass:: EvalfMixin
   :members:

.. autoclass:: PrecisionExhausted
   :members:

.. autofunction:: N

containers
----------
.. module:: sympy.core.containers

.. autoclass:: Tuple
   :members:

.. autoclass:: TupleKind
   :members:

.. autoclass:: Dict
   :members:

exprtools
---------
.. module:: sympy.core.exprtools

.. autofunction:: gcd_terms

.. autofunction:: factor_terms

kind
----
.. module:: sympy.core.kind

.. autoclass:: Kind
   :members:

.. autoclass:: NumberKind
   :members:

.. autoclass:: UndefinedKind
   :members:

.. autoclass:: BooleanKind
   :members:

Sorting
-------


.. autofunction:: sympy.core.sorting.default_sort_key


.. autofunction:: sympy.core.sorting.ordered

Random
------

.. automodule:: sympy.core.random

.. autofunction:: random_complex_number

.. autofunction:: verify_numerically

.. autofunction:: test_derivative_numerically

.. autofunction:: _randrange

.. autofunction:: _randint

Traversal
---------
.. module:: sympy.core.traversal

.. autofunction:: bottom_up

.. autofunction:: postorder_traversal

.. autofunction:: preorder_traversal

.. autofunction:: use

.. autofunction:: walk
