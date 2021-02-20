Combinatorial
=============

This module implements various combinatorial functions.

bell
----

.. autoclass:: sympy.functions.combinatorial.numbers.bell
   :members:
   :undoc-members:
   :private-members:

bernoulli
---------

.. autoclass:: sympy.functions.combinatorial.numbers.bernoulli
   :members:
   :undoc-members:
   :private-members:

binomial
--------

.. autoclass:: sympy.functions.combinatorial.factorials.binomial
   :members:
   :undoc-members:
   :private-members:

catalan
-------

.. autoclass:: sympy.functions.combinatorial.numbers.catalan
   :members:
   :undoc-members:
   :private-members:


euler
-----

.. autoclass:: sympy.functions.combinatorial.numbers.euler
   :members:
   :undoc-members:
   :private-members:


factorial
---------

.. autoclass:: sympy.functions.combinatorial.factorials.factorial
   :members:
   :undoc-members:
   :private-members:

subfactorial
------------

.. autoclass:: sympy.functions.combinatorial.factorials.subfactorial
   :members:
   :undoc-members:
   :private-members:

factorial2 / double factorial
-----------------------------

.. autoclass:: sympy.functions.combinatorial.factorials.factorial2
   :members:
   :undoc-members:
   :private-members:


FallingFactorial
----------------

.. autoclass:: sympy.functions.combinatorial.factorials.FallingFactorial
   :members:
   :undoc-members:
   :private-members:

fibonacci
---------

.. autoclass:: sympy.functions.combinatorial.numbers.fibonacci
   :members:
   :undoc-members:
   :private-members:

tribonacci
----------

.. autoclass:: sympy.functions.combinatorial.numbers.tribonacci
   :members:
   :undoc-members:
   :private-members:


harmonic
--------

.. autoclass:: sympy.functions.combinatorial.numbers.harmonic
   :members:
   :undoc-members:
   :private-members:

lucas
-----

.. autoclass:: sympy.functions.combinatorial.numbers.lucas
   :members:
   :undoc-members:
   :private-members:

genocchi
--------

.. autoclass:: sympy.functions.combinatorial.numbers.genocchi
   :members:
   :undoc-members:
   :private-members:

partition
---------

.. autoclass:: sympy.functions.combinatorial.numbers.partition
   :members:
   :undoc-members:
   :private-members:

MultiFactorial
--------------

.. autoclass:: sympy.functions.combinatorial.factorials.MultiFactorial
   :members:
   :undoc-members:
   :private-members:


RisingFactorial
---------------

.. autoclass:: sympy.functions.combinatorial.factorials.RisingFactorial
   :members:
   :undoc-members:
   :private-members:

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
