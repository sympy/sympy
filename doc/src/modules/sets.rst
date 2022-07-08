.. _sets-module:

====
Sets
====

Basic Sets
----------

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

.. module:: sympy.sets.sets
    :noindex:

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

SymmetricDifference
^^^^^^^^^^^^^^^^^^^
.. autoclass:: SymmetricDifference
   :members:

DisjointUnion
^^^^^^^^^^^^^
.. autoclass:: DisjointUnion
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

Rationals
^^^^^^^^^
.. autoclass:: Rationals
   :members:

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


Reals
^^^^^
.. autoclass:: Reals
   :members:

Complexes
^^^^^^^^^
.. autoclass:: Complexes
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

.. autoclass:: CartesianComplexRegion
   :members:

.. autoclass:: PolarComplexRegion
   :members:

.. autofunction:: normalize_theta_set

Power sets
----------

.. automodule:: sympy.sets.powerset

PowerSet
^^^^^^^^
.. autoclass:: PowerSet
   :members:

Iteration over sets
^^^^^^^^^^^^^^^^^^^

For set unions, `\{a, b\} \cup \{x, y\}` can be treated as
`\{a, b, x, y\}` for iteration regardless of the distinctiveness of
the elements, however, for set intersections, assuming that
`\{a, b\} \cap \{x, y\}` is `\varnothing` or `\{a, b \}` would not
always be valid, since some of `a`, `b`, `x` or `y` may or may not be
the elements of the intersection.

Iterating over the elements of a set involving intersection, complement,
or symmetric difference yields (possibly duplicate) elements of the set
provided that all elements are known to be the elements of the set.
If any element cannot be determined to be a member of a set then the
iteration gives ``TypeError``.
This happens in the same cases where ``x in y`` would give an error.

There are some reasons to implement like this, even if it breaks the
consistency with how the python set iterator works.
We keep in mind that sympy set comprehension like ``FiniteSet(*s)`` from
a existing sympy sets could be a common usage.
And this approach would make ``FiniteSet(*s)`` to be consistent with any
symbolic set processing methods like ``FiniteSet(*simplify(s))``.

Condition Sets
--------------

.. automodule:: sympy.sets.conditionset

ConditionSet
^^^^^^^^^^^^

.. autoclass:: ConditionSet
   :members:

Relations on sets
^^^^^^^^^^^^^^^^^

.. autoclass:: Contains
   :members:

SetKind
--------------

.. autoclass:: SetKind
   :members:
