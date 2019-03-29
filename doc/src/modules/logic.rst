Logic Module
============

.. module:: sympy.logic

Introduction
------------

The logic module for SymPy allows to form and manipulate logic expressions
using symbolic and Boolean values.

Forming logical expressions
---------------------------

You can build Boolean expressions with the standard python operators ``&``
(:class:`And`), ``|`` (:class:`Or`), ``~`` (:class:`Not`)::

    >>> from sympy import *
    >>> x, y = symbols('x,y')
    >>> y | (x & y)
    y | (x & y)
    >>> x | y
    x | y
    >>> ~x
    ~x

You can also form implications with ``>>`` and ``<<``::

    >>> x >> y
    Implies(x, y)
    >>> x << y
    Implies(y, x)

Like most types in SymPy, Boolean expressions inherit from :class:`Basic`::

    >>> (y & x).subs({x: True, y: True})
    True
    >>> (x | y).atoms()
    {x, y}

The logic module also includes the following functions to derive boolean expressions
from their truth tables-

.. autofunction:: sympy.logic.boolalg.SOPform

.. autofunction:: sympy.logic.boolalg.POSform

Boolean functions
-----------------

.. autoclass:: sympy.logic.boolalg.BooleanTrue

.. autoclass:: sympy.logic.boolalg.BooleanFalse

.. autoclass:: sympy.logic.boolalg.And

.. autoclass:: sympy.logic.boolalg.Or

.. autoclass:: sympy.logic.boolalg.Not

.. autoclass:: sympy.logic.boolalg.Xor

.. autoclass:: sympy.logic.boolalg.Nand

.. autoclass:: sympy.logic.boolalg.Nor

.. autoclass:: sympy.logic.boolalg.Implies

.. autoclass:: sympy.logic.boolalg.Equivalent

.. autoclass:: sympy.logic.boolalg.ITE

The following functions can be used to handle Conjunctive and Disjunctive Normal
forms-

.. autofunction:: sympy.logic.boolalg.to_cnf

.. autofunction:: sympy.logic.boolalg.to_dnf

.. autofunction:: sympy.logic.boolalg.is_cnf

.. autofunction:: sympy.logic.boolalg.is_dnf

Simplification and equivalence-testing
--------------------------------------

.. autofunction:: sympy.logic.boolalg.simplify_logic

SymPy's simplify() function can also be used to simplify logic expressions to their
simplest forms.

.. autofunction:: sympy.logic.boolalg.bool_map

Inference
---------

.. module:: sympy.logic.inference

This module implements some inference routines in propositional logic.

The function satisfiable will test that a given Boolean expression is satisfiable,
that is, you can assign values to the variables to make the sentence `True`.

For example, the expression ``x & ~x`` is not satisfiable, since there are no
values for ``x`` that make this sentence ``True``. On the other hand, ``(x
| y) & (x | ~y) & (~x | y)`` is satisfiable with both ``x`` and ``y`` being
``True``.

    >>> from sympy.logic.inference import satisfiable
    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> satisfiable(x & ~x)
    False
    >>> satisfiable((x | y) & (x | ~y) & (~x | y))
    {x: True, y: True}

As you see, when a sentence is satisfiable, it returns a model that makes that
sentence ``True``. If it is not satisfiable it will return ``False``.

.. autofunction:: sympy.logic.inference.satisfiable

.. TODO: write about CNF file format
