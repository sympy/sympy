.. _combinatorial-functions:

Combinatorial
=============

This module implements various combinatorial functions.

.. autoclass:: sympy.functions.combinatorial.numbers.bell
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.bernoulli
   :members:

.. autoclass:: sympy.functions.combinatorial.factorials.binomial
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.catalan
   :members:


.. autoclass:: sympy.functions.combinatorial.numbers.euler
   :members:


.. autoclass:: sympy.functions.combinatorial.factorials.factorial
   :members:

.. autoclass:: sympy.functions.combinatorial.factorials.subfactorial
   :members:

.. autoclass:: sympy.functions.combinatorial.factorials.factorial2
   :members:


.. autoclass:: sympy.functions.combinatorial.factorials.FallingFactorial
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.fibonacci
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.tribonacci
   :members:


.. autoclass:: sympy.functions.combinatorial.numbers.harmonic
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.lucas
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.genocchi
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.andre
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.partition
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.divisor_sigma
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.udivisor_sigma
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.legendre_symbol
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.jacobi_symbol
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.kronecker_symbol
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.mobius
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.primenu
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.primeomega
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.totient
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.reduced_totient
   :members:

.. autoclass:: sympy.functions.combinatorial.numbers.primepi
   :members:

.. autoclass:: sympy.functions.combinatorial.factorials.MultiFactorial
   :members:


.. autoclass:: sympy.functions.combinatorial.factorials.RisingFactorial
   :members:

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

.. autofunction:: sympy.functions.combinatorial.numbers.nC

.. autofunction:: sympy.functions.combinatorial.numbers.nP

.. autofunction:: sympy.functions.combinatorial.numbers.nT
