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
relators. If the list of relators is empty, the associated free group is
returned.

Example of construction of a finitely-presented group.
The symmetric group of degree 4 may be represented as a two generator group
with presentation `< a, b | a^2, b^3, (a*b)^4 >`. Giving the relations as a
list of relators, the presentation of would be specified as:

>>> F, a, b = free_group("a, b")
>>> G = FpGroup(F, [a**2, b**3, (a*b)**4])
>>> G
<free group on the generators (a, b)>

Currently groups with relators having presentation like
`< r, s | r^2, s^2, t^2, rst = str = trs >` will have to be specified as:

>>> F, r, s = free_group("r, s")
>>> G = FpGroup(F, [r**2, s**2, t**2, r*s*t*r**-1*t**-1*s**-1, s*t*r*s**-1*r**-1*t**-1])

obviously this is not a unique way to make such a group, but the point is that
in case of equality with non-identity the user has to manually do that.

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

This section describes the use of coset enumeration techniques in SymPy. The
algorithm used for coset enumeration procedure is Todd-Coxeter algorithm and
is developed in Sympy using [Ho05] and [CDHW73]. The reader should consult
[CDHW73] and [Hav91] for a general description of the algorithm.

For experienced users we have a number of parameters that can be used to
manipulate the algorithm, like

``coset_table_max_limit``: manipulate the maximum number of cosets allowed in
coset enumeration.

``max_stack_size``
maximum size of ``deduction_stack`` above to equal to which the stack is
emptied.

Compression and Standardization
-------------------------------


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
    https://www.crcpress.com/Handbook-of-Computational-Group-Theory/Holt-Eick-OBrien/p/book/9781584883722

    A practical method for enumerating cosets of a finite abstract group
