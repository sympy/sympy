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
   :undoc-members:
   :private-members:

.. autofunction:: imageset

Elementary Sets
---------------

Interval
^^^^^^^^
.. autoclass:: Interval
   :members:
   :undoc-members:
   :private-members:

FiniteSet
^^^^^^^^^
.. autoclass:: FiniteSet
   :members:
   :undoc-members:
   :private-members:

Compound Sets
-------------

.. module:: sympy.sets.sets
    :noindex:

Union
^^^^^
.. autoclass:: Union
   :members:
   :undoc-members:
   :private-members:

Intersection
^^^^^^^^^^^^
.. autoclass:: Intersection
   :members:
   :undoc-members:
   :private-members:

ProductSet
^^^^^^^^^^
.. autoclass:: ProductSet
   :members:
   :undoc-members:
   :private-members:

Complement
^^^^^^^^^^
.. autoclass:: Complement
   :members:
   :undoc-members:
   :private-members:

SymmetricDifference
^^^^^^^^^^^^^^^^^^^
.. autoclass:: SymmetricDifference
   :members:
   :undoc-members:
   :private-members:

Singleton Sets
--------------

EmptySet
^^^^^^^^
.. autoclass:: EmptySet
   :members:
   :undoc-members:
   :private-members:

UniversalSet
^^^^^^^^^^^^
.. autoclass:: UniversalSet
   :members:
   :undoc-members:
   :private-members:

Special Sets
------------
.. automodule:: sympy.sets.fancysets

Naturals
^^^^^^^^
.. autoclass:: Naturals
   :members:
   :undoc-members:
   :private-members:

Naturals0
^^^^^^^^^
.. autoclass:: Naturals0
   :members:
   :undoc-members:
   :private-members:

Integers
^^^^^^^^
.. autoclass:: Integers
   :members:
   :undoc-members:
   :private-members:


Reals
^^^^^
.. autoclass:: Reals
   :members:
   :undoc-members:
   :private-members:

Complexes
^^^^^^^^^
.. autoclass:: Complexes
   :members:
   :undoc-members:
   :private-members:

ImageSet
^^^^^^^^
.. autoclass:: ImageSet
   :members:
   :undoc-members:
   :private-members:

Range
^^^^^
.. autoclass:: Range
   :members:
   :undoc-members:
   :private-members:

ComplexRegion
^^^^^^^^^^^^^
.. autoclass:: ComplexRegion
   :members:
   :undoc-members:
   :private-members:

.. autoclass:: CartesianComplexRegion
   :members:
   :undoc-members:
   :private-members:

.. autoclass:: PolarComplexRegion
   :members:
   :undoc-members:
   :private-members:

.. autofunction:: normalize_theta_set

Power sets
----------

.. automodule:: sympy.sets.powerset

PowerSet
^^^^^^^^
.. autoclass:: PowerSet
   :members:
   :undoc-members:
   :private-members:

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
   :undoc-members:
   :private-members:
