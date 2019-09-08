====
Sets
====

.. automodule:: sympy.sets.sets

Set
^^^
.. autoclass:: Set
   :members:

.. autofunction:: imageset

Elementary Sets
---------------

Interval
^^^^^^^^
.. autoclass:: Interval
   :members:

FiniteSet
^^^^^^^^^
.. autoclass:: FiniteSet
   :members:

Compound Sets
-------------
Union
^^^^^
.. autoclass:: Union
   :members:

Intersection
^^^^^^^^^^^^
.. autoclass:: Intersection
   :members:

ProductSet
^^^^^^^^^^
.. autoclass:: ProductSet
   :members:

Complement
^^^^^^^^^^
.. autoclass:: Complement
   :members:

Singleton Sets
--------------

EmptySet
^^^^^^^^
.. autoclass:: EmptySet
   :members:

UniversalSet
^^^^^^^^^^^^
.. autoclass:: UniversalSet
   :members:

Special Sets
------------
.. automodule:: sympy.sets.fancysets

Naturals
^^^^^^^^
.. autoclass:: Naturals
   :members:

Naturals0
^^^^^^^^^
.. autoclass:: Naturals0
   :members:

Integers
^^^^^^^^
.. autoclass:: Integers
   :members:

ImageSet
^^^^^^^^
.. autoclass:: ImageSet
   :members:

Range
^^^^^
.. autoclass:: Range
   :members:

ComplexRegion
^^^^^^^^^^^^^
.. autoclass:: ComplexRegion
   :members:

.. autofunction:: normalize_theta_set

Iteration over sets
^^^^^^^^^^^^^^^^^^^

For set unions, `\{a, b\} \cup \{x, y\}` may have no problem to be
comprehended like `\{a, b, x, y\}` regardless of the distinctiveness of
the elements, however, for set intersections, assuming that
`\{a, b\} \cap \{x, y\}` is `\varnothing` or `\{a, b \}` would not
always be valid.

So, we have designed ``__iter__`` method for composite sets to sample an
element from a enumeratible set and testing its membership for an
another set, until an undecidable `\in` appears, which would raise a
``TypeError``.

There are some reasons to implement like this, even if it breaks the
consistency with how the python set iterator works.
We keep in mind that sympy set comprehension like ``FiniteSet(*s)`` from
a existing sympy sets could be a common usage, and at least try to give
a consistent result with ``s`` simplified in a symbolic way,
than giving a mathematically wrong iteration.
