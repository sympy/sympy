Introduction
============

This module presents the functionality designed for computing with finitely-
presented groups (fp-groups for short). The name of the corresponding SymPy
object is ``FpGroup``. The functions or classes described here are studied
under **computational group theory**.

Overview of Facilities
----------------------

The facilities provided for fp-groups fall into a number of natural groupings

* The construction of fp-groups using a free group and a list of words in
  generators of that free group.

* Index determination using the famous Todd-Coxeter procedure.

* The construction of all subgroups having index less than some (small)
  specified positive integer, using the *Low-Index Subgroups* algorithm.

For a description of fundamental algorithms of finitely presented groups
we often make use of *Handbook of Computational Group Theory*.

The Construction of Finitely Presented Groups
---------------------------------------------

Finitely presented groups are construced by factoring a free group by a
set of relators. The set of relators is taken in as a list of words in
generators of free group in SymPy, using a list provides ordering to the
relators. If the list of relators is empty, the assosciated free group is
returned.

For example:

    >>> F, x, y = free_group("x, y")
    >>> f = FpGroup(F, [x**2, y**3, (x*y)**4])
    >>> f
    <fp group on the generators (x, y)>

Free Groups and Words
=====================

Construction of a Free Group
----------------------------

``free_group("gen0, gen1, \ldots gen_(n-1)")`` constructs a free group ``F`` on ``n``
generators, where ``n`` is a positive integer.
The `i`-th generator of `F` may be obtained using the method `.generators[i]`, `i = 0, \ldots n-1`.

Ex

F, x, y = free_group("x, y")

creates a free group of rank 2 and assigns the variables ``x`` and ``y`` to the two
generators.

F = vfree_group("x, y")

creates a free group ``F[0]`` of rank 2, and a tuple of generators ``F[1]``.


Construction of words
---------------------

Can be called with different public versions of ``FreeGroup`` i.e ``free_group``
, ``vfree_group`` and ``xfree_group``.

Methods for manipulating words
------------------------------

This section describes some basic access functions for words. These operations apply
to both, free groups and arbitrary fp-groups.

copy
return a

is_identity

array_form

letter_form

index

ext_rep

contains

eliminate_word

len

``in``

``len(w)``
The lenth of word ``w``.


Comparison of words
-------------------

Coset Enumeration: The Todd-Coxeter Algorithm
=============================================

Subgroups of Finite Index: Low Index Subgroups algorithm
========================================================

Bibliography
============

[CDHW73]
    John J. Cannon, Lucien A. Dimino, George Havas, and Jane M. Watson.
    Implementation and analysis of the Todd-Coxeter algorithm. Math. Comp., 27:463–
    490, 1973.

[Ho05]
    Derek F. Holt,
    Handbook of Computational Group Theory.
    In the series 'Discrete Mathematics and its Applications',
