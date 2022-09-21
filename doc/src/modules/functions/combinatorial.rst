.. _combinatorial-functions:

Combinatorial
=============

This module implements various combinatorial functions.

bell
----

.. autoclass:: sympy.functions.combinatorial.numbers.bell
   :members:

bernoulli
---------

.. autoclass:: sympy.functions.combinatorial.numbers.bernoulli
   :members:

binomial
--------

.. autoclass:: sympy.functions.combinatorial.factorials.binomial
   :members:

catalan
-------

.. autoclass:: sympy.functions.combinatorial.numbers.catalan
   :members:


euler
-----

.. autoclass:: sympy.functions.combinatorial.numbers.euler
   :members:


factorial
---------

.. autoclass:: sympy.functions.combinatorial.factorials.factorial
   :members:

subfactorial
------------

.. autoclass:: sympy.functions.combinatorial.factorials.subfactorial
   :members:

factorial2 / double factorial
-----------------------------

.. autoclass:: sympy.functions.combinatorial.factorials.factorial2
   :members:


FallingFactorial
----------------

.. autoclass:: sympy.functions.combinatorial.factorials.FallingFactorial
   :members:

fibonacci
---------

.. autoclass:: sympy.functions.combinatorial.numbers.fibonacci
   :members:

tribonacci
----------

.. autoclass:: sympy.functions.combinatorial.numbers.tribonacci
   :members:


harmonic
--------

.. autoclass:: sympy.functions.combinatorial.numbers.harmonic
   :members:

lucas
-----

.. autoclass:: sympy.functions.combinatorial.numbers.lucas
   :members:

genocchi
--------

.. autoclass:: sympy.functions.combinatorial.numbers.genocchi
   :members:

andre
-----

.. autoclass:: sympy.functions.combinatorial.numbers.andre
   :members:

partition
---------

.. autoclass:: sympy.functions.combinatorial.numbers.partition
   :members:

MultiFactorial
--------------

.. autoclass:: sympy.functions.combinatorial.factorials.MultiFactorial
   :members:


RisingFactorial
---------------

.. autoclass:: sympy.functions.combinatorial.factorials.RisingFactorial
   :members:

stirling
--------

.. autofunction:: sympy.functions.combinatorial.numbers.stirling

Enumeration
===========

Three functions are available. Each of them attempts to efficiently compute
a given combinatorial quantity for a given set or multiset which can be
entered as an integer, sequence or multiset (dictionary with
elements as keys and multiplicities as values). The ``k`` parameter indicates
the number of elements to pick (or the number of partitions to make). When
``k`` is None, the sum of the enumeration for all ``k`` (from 0 through the
number of items represented by ``n``) is returned. A ``replacement`` parameter
is recognized for combinations and permutations; this indicates that any item
may appear with multiplicity as high as the number of items in the original
set.

>>> from sympy.functions.combinatorial.numbers import nC, nP, nT
>>> items = 'baby'

nC
--

.. autofunction:: sympy.functions.combinatorial.numbers.nC

nP
--

.. autofunction:: sympy.functions.combinatorial.numbers.nP

nT
--

.. autofunction:: sympy.functions.combinatorial.numbers.nT
