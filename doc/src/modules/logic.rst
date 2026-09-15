=====
Logic
=====

.. module:: sympy.logic

Introduction
------------

The logic module for SymPy allows to form and manipulate logic expressions
using symbolic and Boolean values.

Forming logical expressions
---------------------------

You can build Boolean expressions with the standard python operators ``&``
(:class:`~.And`), ``|`` (:class:`~.Or`), ``~`` (:class:`~.Not`)::

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

Like most types in SymPy, Boolean expressions inherit from :class:`~.Basic`::

    >>> (y & x).subs({x: True, y: True})
    True
    >>> (x | y).atoms()
    {x, y}

The logic module also includes the following functions to derive boolean expressions
from their truth tables:

.. autofunction:: sympy.logic.boolalg::SOPform

.. autofunction:: sympy.logic.boolalg::POSform

.. autofunction:: sympy.logic.boolalg::ANFform

Boolean functions
-----------------

.. autoclass:: sympy.logic.boolalg::Boolean
   :members:

.. autoclass:: sympy.logic.boolalg::BooleanTrue
   :members:

.. autoclass:: sympy.logic.boolalg::BooleanFalse
   :members:

.. autoclass:: sympy.logic.boolalg::And
   :members:

.. autoclass:: sympy.logic.boolalg::Or
   :members:

.. autoclass:: sympy.logic.boolalg::Not
   :members:

.. autoclass:: sympy.logic.boolalg::Xor
   :members:

.. autoclass:: sympy.logic.boolalg::Nand
   :members:

.. autoclass:: sympy.logic.boolalg::Nor
   :members:

.. autoclass:: sympy.logic.boolalg::Xnor
   :members:

.. autoclass:: sympy.logic.boolalg::Implies
   :members:

.. autoclass:: sympy.logic.boolalg::Equivalent
   :members:

.. autoclass:: sympy.logic.boolalg::ITE
   :members:

.. autoclass:: sympy.logic.boolalg::Exclusive
   :members:

The following functions can be used to handle Algebraic, Conjunctive,
Disjunctive, and Negated Normal forms:

.. autofunction:: sympy.logic.boolalg::to_anf

.. autofunction:: sympy.logic.boolalg::to_cnf

.. autofunction:: sympy.logic.boolalg::to_dnf

.. autofunction:: sympy.logic.boolalg::to_nnf

.. autofunction:: sympy.logic.boolalg::is_anf

.. autofunction:: sympy.logic.boolalg::is_cnf

.. autofunction:: sympy.logic.boolalg::is_dnf

.. autofunction:: sympy.logic.boolalg::is_nnf

.. autofunction:: sympy.logic.boolalg::gateinputcount

Simplification and equivalence-testing
--------------------------------------

.. autofunction:: sympy.logic.boolalg::simplify_logic

SymPy's :py:func:`~.simplify` function can also be used to simplify logic expressions to their
simplest forms.

.. autofunction:: sympy.logic.boolalg::bool_map

Manipulating expressions
------------------------

The following functions can be used to manipulate Boolean expressions:

.. autofunction:: sympy.logic.boolalg::distribute_and_over_or

.. autofunction:: sympy.logic.boolalg::distribute_or_over_and

.. autofunction:: sympy.logic.boolalg::distribute_xor_over_and

.. autofunction:: sympy.logic.boolalg::eliminate_implications

Truth tables and related functions
----------------------------------

It is possible to create a truth table for a Boolean function with

.. autofunction:: sympy.logic.boolalg::truth_table

For mapping between integer representations of truth table positions, lists of
zeros and ones and symbols, the following functions can be used:

.. autofunction:: sympy.logic.boolalg::integer_to_term

.. autofunction:: sympy.logic.boolalg::term_to_integer

.. autofunction:: sympy.logic.boolalg::bool_maxterm

.. autofunction:: sympy.logic.boolalg::bool_minterm

.. autofunction:: sympy.logic.boolalg::bool_monomial

.. autofunction:: sympy.logic.boolalg::anf_coeffs

.. autofunction:: sympy.logic.boolalg::to_int_repr


Inference
---------

.. module:: sympy.logic.inference

This module implements some inference routines in propositional logic.

The function satisfiable will test that a given Boolean expression is satisfiable,
that is, you can assign values to the variables to make the sentence ``True``.

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

.. autofunction:: sympy.logic.inference::satisfiable

.. TODO: write about CNF file format
