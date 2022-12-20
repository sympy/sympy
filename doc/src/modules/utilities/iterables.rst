=========
Iterables
=========


.. automodule:: sympy.utilities.iterables
   :members:


variations
----------

variations(seq, n) Returns all the variations of the list of size n.

Has an optional third argument. Must be a boolean value and makes the method
return the variations with repetition if set to True, or the variations
without repetition if set to False.

Examples::
    >>> from sympy.utilities.iterables import variations
    >>> list(variations([1,2,3], 2))
    [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]
    >>> list(variations([1,2,3], 2, True))
    [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]


partitions
----------

Although the combinatorics module contains Partition and IntegerPartition
classes for investigation and manipulation of partitions, there are a few
functions to generate partitions that can be used as low-level tools for
routines:  ``partitions`` and ``multiset_partitions``. The former gives
integer partitions, and the latter gives enumerated partitions of elements.
There is also a routine ``kbins`` that will give a variety of permutations
of partitions. And to obtain partitions as a list instead of a dictionary, there
is the ``ordered_partition`` function which is quite fast. Finally, to simply
obtain a count of the number of partitions without enumerating them, there
is the ``nT`` function.

See Also
========

sympy.utilities.iterables.ordered_partitions, sympy.functions.combinatorial.numbers.nT

partitions::

    >>> from sympy.utilities.iterables import partitions
    >>> [p.copy() for s, p in partitions(7, m=2, size=True) if s == 2]
    [{1: 1, 6: 1}, {2: 1, 5: 1}, {3: 1, 4: 1}]

multiset_partitions::

    >>> from sympy.utilities.iterables import multiset_partitions
    >>> [p for p in multiset_partitions(3, 2)]
    [[[0, 1], [2]], [[0, 2], [1]], [[0], [1, 2]]]
    >>> [p for p in multiset_partitions([1, 1, 1, 2], 2)]
    [[[1, 1, 1], [2]], [[1, 1, 2], [1]], [[1, 1], [1, 2]]]

kbins::

    >>> from sympy.utilities.iterables import kbins
    >>> def show(k):
    ...     rv = []
    ...     for p in k:
    ...         rv.append(','.join([''.join(j) for j in p]))
    ...     return sorted(rv)
    ...
    >>> show(kbins("ABCD", 2))
    ['A,BCD', 'AB,CD', 'ABC,D']
    >>> show(kbins("ABC", 2))
    ['A,BC', 'AB,C']
    >>> show(kbins("ABC", 2, ordered=0))  # same as multiset_partitions
    ['A,BC', 'AB,C', 'AC,B']
    >>> show(kbins("ABC", 2, ordered=1))
    ['A,BC', 'A,CB',
     'B,AC', 'B,CA',
     'C,AB', 'C,BA']
    >>> show(kbins("ABC", 2, ordered=10))
    ['A,BC', 'AB,C', 'AC,B',
     'B,AC', 'BC,A',
     'C,AB']
    >>> show(kbins("ABC", 2, ordered=11))
    ['A,BC', 'A,CB', 'AB,C', 'AC,B',
     'B,AC', 'B,CA', 'BA,C', 'BC,A',
     'C,AB', 'C,BA', 'CA,B', 'CB,A']
