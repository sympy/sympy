===========
Enumerative
===========
.. module:: sympy.utilities.enumerative

This module includes functions and classes for enumerating and
counting multiset partitions.

.. autofunction:: multiset_partitions_taocp
.. autofunction:: factoring_visitor
.. autofunction:: list_visitor

The approach of the function ``multiset_partitions_taocp`` is extended
and generalized by the class ``MultisetPartitionTraverser``.

.. autoclass:: MultisetPartitionTraverser
   :members: count_partitions,
             enum_all,
             enum_large,
             enum_range,
             enum_small
