===========
Enumerative
===========
.. module:: sympy.utilities.enumerative

This module includes functions and classes for enumerating and
counting multiset partitions.

.. autofunction:: multiset_partitions_taocp
.. autofunction:: factoring_visitor
.. autofunction:: list_visitor

And also a class, which reimplements and extends the basic algorithm

.. autoclass:: MultisetPartitionTraverser
   :members: count_partitions,
             enum_all,
             enum_large,
             enum_range,
             enum_small

And another class, which is really only useful internally

.. autoclass:: PartComponent
   :members:



