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

harmonic
--------

.. autoclass:: sympy.functions.combinatorial.numbers.harmonic
   :members:


lucas
-----

.. autoclass:: sympy.functions.combinatorial.numbers.lucas
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

Calculate the number of combinations of length ``k``.

>>> [nC(items, k) for k in range(len(items) + 1)], nC(items)
([1, 3, 4, 3, 1], 12)
>>> nC('aaa', 2)
1
>>> nC('abc', 2)
3
>>> nC(3, 2)
3

nP
--

Calculate the number of permutations of length ``k``.

>>> [nP(items, k) for k in range(len(items) + 1)], nP(items)
([1, 3, 7, 12, 12], 35)
>>> nC('aaa', 2)
1
>>> nC('abc', 2)
3
>>> nC(3, 2)
3

nT
--

Calculate the number of partitions that have ``k`` parts.

>>> [nT(items, k) for k in range(len(items) + 1)], nT(items)
([0, 1, 5, 4, 1], 11)
>>> nT('aaa', 2)
1
>>> nT('abc', 2)
3
>>> nT(3, 2)
1

Note that the integer for ``n`` indicates *identical* items for ``nT`` but
indicates ``n`` *different* items for ``nC`` and ``nP``.
