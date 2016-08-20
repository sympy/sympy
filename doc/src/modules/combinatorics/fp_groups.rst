Finitely Presented Groups
=========================

Introduction
------------

This module presents the functionality designed for computing with finitely-
presented groups (fp-groups for short). The name of the corresponding SymPy
object is ``FpGroup``. The functions or classes described here are studied
under **computational group theory**. All code examples assume:

>>> from sympy.combinatorics.free_groups import free_group, vfree_group, xfree_group
>>> from sympy.combinatorics.fp_groups import FpGroup, CosetTable, coset_enumeration_r

Overview of Facilities
``````````````````````

The facilities provided for fp-groups fall into a number of natural groupings

* The construction of fp-groups using a free group and a list of words in
  generators of that free group.

* Index determination using the famous Todd-Coxeter procedure.

* The construction of all subgroups having index less than some (small)
  specified positive integer, using the *Low-Index Subgroups* algorithm.

* Algorithms for computing presentations of a subgroup of finite index
  in a group defined by finite presentation.

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
with presentation `< a, b \mid a^2, b^3, (a*b)^4 >`. Giving the relations as a
list of relators, the presentation of would be specified as:

>>> F, a, b = free_group("a, b")
>>> G = FpGroup(F, [a**2, b**3, (a*b)**4])
>>> G
<fp group on the generators (a, b)>

Currently groups with relators having presentation like
`< r, s, t \mid r^2, s^2, t^2, rst = str = trs >` will have to be specified as:

>>> F, r, s, t = free_group("r, s, t")
>>> G = FpGroup(F, [r**2, s**2, t**2, r*s*t*r**-1*t**-1*s**-1, s*t*r*s**-1*r**-1*t**-1])

obviously this is not a unique way to make such a group, but the point is that
in case of equality with non-identity the user has to manually do that.

Free Groups and Words
---------------------

Construction of a Free Group
````````````````````````````

``free_group("gen0, gen1, \ldots gen_(n-1)")`` constructs a free group ``F`` on ``n``
generators, where ``n`` is a positive integer.
The `i`-th generator of `F` may be obtained using the method ``.generators[i]``, `i = 0, \ldots n-1`.

Ex

>>> F, x, y = free_group("x, y")

creates a free group of rank 2 and assigns the variables ``x`` and ``y`` to the two
generators.

>>> F = vfree_group("x, y")

creates a free group ``F[0]`` of rank 2, and a tuple of generators ``F[1]``.


Construction of words
`````````````````````

Can be called with different public versions of ``FreeGroup`` i.e
``free_group``, ``vfree_group`` and ``xfree_group``.

Methods for manipulating Words
``````````````````````````````

This section describes some basic access functions for *words*. These operations apply
to words of both ``FreeGroup`` as well as that of arbitrary ``FpGroup``.

``copy``
    return a

``is_identity``

``array_form``
    Returns one of the internal use representation of ``word``.

``letter_form``

index

ext_rep

contains

eliminate_word

``in``

``len(w)``
    The lenth of word ``w``.


Comparison of words
```````````````````

Coset Enumeration: The Todd-Coxeter Algorithm
---------------------------------------------

This section describes the use of coset enumeration techniques in SymPy. The
algorithm used for coset enumeration procedure is Todd-Coxeter algorithm and
is developed in Sympy using [Ho05] and [CDHW73]. The reader should consult
[CDHW73] and [Hav91] for a general description of the algorithm.

We have two strategies of coset enumeration *relator-based* and
*coset-table based* and the two have been implemented as
``coset_enumeration_r``, ``coset_enumeration_c`` respectively. The two
strategies differ in the way they make new definitions for the cosets.

Though from the user point of view it is suggested to rather use the
``.coset_enumeration`` method of ``FpGroup`` and specify the ``strategy``
argument

``strategy``:
    (default="relator_based") specifies the strategy of coset
    enumeration to be used, possible values are *"relator_based"* or
    *"coset_table_based"*.

CosetTable
``````````

Class used to manipulate the information regarding the coset enumeration of
the finitely presented group ``G`` on the cosets of the subgroup ``H``.

Basically a *coset table* ``CosetTable(G,H)``, is the permutation representation
of the finitely presented group on the cosets of a subgroup. Most of the set
theoretic and group functions use the regular representation of ``G``, i.e.,
the coset table of ``G`` over the trivial subgroup.

The actual mathematical coset table is obtained using ``.table`` attribute and
is a list of lists. For each generator ``g`` of ``G`` it contains a coloumn and
the next coloumn corresponds to ``g**-1`` and so on for other generators, so in
total it has ``2*G.rank()`` coloumns. Each coloumn is simply a list of integers.
If ``l`` is the generator list for the generator `g` and if ``l[i] = j`` then
generator ``g`` takes the coset `i` to the coset `j` by multiplication from the
right.

For finitely presented groups, a coset table is computed by a Todd-Coxeter
coset enumeration. Note that you may influence the performance of that
enumeration by changing the values of the variable
``CosetTable.coset_table_max_limit``.

Attributes of CosetTable
````````````````````````

For ``CosetTable(G, H)`` where ``G`` is the group and ``H`` is the subgroup.

* ``n``: A non-negative integer, non-mutable attriburte, dependently
  calculated as the maximum among the live-cosets (i.e `\Omega`).

* ``table``: A list of lists, mutable attribute, mathematically represents the
  coset table.

* ``omega``: A list, dependent on the internal attribute ``p``. `\Omega`
  represents the list of live-cosets. A *standard* coset-table has its
  `\Omega = [0, 1, \ldots, index-1]` where `index` is the index of subgroup
  ``H`` in ``G``.

For experienced users we have a number of parameters that can be used to
manipulate the algorithm, like

* ``coset_table_max_limit`` (default value = `4096000`): manipulate the maximum
  number of cosets allowed in coset enumeration, i.e the number of rows allowed
  in coset table. A coset enumeration will not finish if the subgroup does not
  have finite index, and even if it has it may take many more intermediate
  cosets than the actual index of the subgroup is. To avoid a coset enumeration
  "running away" therefore SymPy has a "safety stop" built-in. This is
  controlled by this variable. For example:

  >>> CosetTable.coset_table_max_limit = 50
  >>> F, a, b = free_group("a, b")
  >>> Cox = FpGroup(F, [a**6, b**6, (a*b)**2, (a**2*b**2)**2, (a**3*b**3)**5])
  >>> C_r = coset_enumeration_r(Cox, [a])
  Traceback (most recent call last):
    ...
  ValueError: the coset enumeration has defined more than 50 cosets


* ``max_stack_size`` (default value = `500`): manipulate the maximum size of
  ``deduction_stack`` above or equal to which the stack is emptied.

Compression and Standardization
```````````````````````````````

For any two entries `i, j` with `i < j` in coset table, the first
occurrence of `i` in a coset table precedes the first occurrence of `j` with
respect to the usual row-wise ordering of the table entries. We call such a
table a standard coset table. To standardize a ``CosetTable`` we use the
``.standardize`` method.

**Note** the method alters the given table, it does not create a copy.

Subgroups of Finite Index
-------------------------

The functionality in this section are concerned with the construction of
subgroups of finite index. We describe a method for computing all subgroups
whose index does not exceed some (modest) integer bound.

Low Index Subgroups
```````````````````

``low_index_subgroups(G, N)``: Given a finitely presented group ``G`` (can be a free
group), and ``N`` a positive integer, determine the conjugacy classes of
subgroups of ``G`` whose indices is less than or equal to ``N``.

Bibliography
------------

[CDHW73]
    John J. Cannon, Lucien A. Dimino, George Havas, and Jane M. Watson.
    Implementation and analysis of the Todd-Coxeter algorithm. Math. Comp., 27:463–
    490, 1973.

[Ho05]
    Derek F. Holt,
    Handbook of Computational Group Theory.
    In the series 'Discrete Mathematics and its Applications',
    `Chapman & Hall/CRC 2005, xvi + 514 p <https://www.crcpress.com/Handbook-of-Computational-Group-Theory/Holt-Eick-OBrien/p/book/9781584883722>`_.

    A practical method for enumerating cosets of a finite abstract group
