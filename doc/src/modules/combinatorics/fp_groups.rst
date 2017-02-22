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

Finitely presented groups are constructed by factoring a free group by a
set of relators. The set of relators is taken in as a list of words in
generators of free group in SymPy, using a list provides ordering to the
relators. If the list of relators is empty, the associated free group is
returned.

Example of construction of a finitely-presented group.
The symmetric group of degree 4 may be represented as a two generator group
with presentation `\langle a, b \mid a^2, b^3, (ab)^4 \rangle`. Giving the relations as a
list of relators, group in SymPy would be specified as:

>>> F, a, b = free_group("a, b")
>>> G = FpGroup(F, [a**2, b**3, (a*b)**4])
>>> G
<fp group on the generators (a, b)>

Currently groups with relators having presentation like
`\langle r, s, t \mid r^2, s^2, t^2, rst = str = trs \rangle` will have to be specified as:

>>> F, r, s, t = free_group("r, s, t")
>>> G = FpGroup(F, [r**2, s**2, t**2, r*s*t*r**-1*t**-1*s**-1, s*t*r*s**-1*r**-1*t**-1])

Obviously this is not a unique way to make that particular group, but the point
is that in case of equality with non-identity the user has to manually do that.

Free Groups and Words
---------------------

Construction of a Free Group
````````````````````````````

``free_group("gen0, gen1, ..., gen_(n-1)")`` constructs a free group ``F`` on
``n`` generators, where ``n`` is a positive integer. The `i`-th generator of
`F` may be obtained using the method ``.generators[i]``, `i = 0, \ldots n-1`.

>>> F, x, y = free_group("x, y")

creates a free group ``F`` of rank 2 and assigns the variables ``x`` and ``y``
to the two generators.

>>> F = vfree_group("x, y")
>>> F
<free group on the generators (x, y)>

creates a free group ``F`` of rank 2, with tuple of generators ``F.generators``,
and inserts ``x`` and ``y`` as generators into the global namespace.

>>> F = xfree_group("x, y")
>>> F
(<free group on the generators (x, y)>, (x, y))
>>> x**2
x**2

creates a free groups ``F[0]`` of rank 2, with tuple of generators ``F[1]``.

Construction of words
`````````````````````

This section is applicable to words of ``FreeGroup`` as well as ``FpGroup``.
When we say *word* in SymPy, it actually means a `reduced word
<https://en.wikipedia.org/wiki/Word_(group_theory)#Reduced_words>`_ , since the
words are automatically reduced. Given a group ``G`` defined on `n` generators
`x_1, x_2, x_3, \ldots, x_n`, a word is constructed as
`s_1^{r_1}s_2^{r_2} \cdots s_k^{r_k}` where `s_i \in \{x_1, x_2, \ldots, x_n\}`
, `r_i \in \mathbb{Z}` for all `k`.

Each word can be constructed in a variety of ways, since after reduction they
may be equivalent.

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
argument.

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
is a list of lists. For each generator ``g`` of ``G`` it contains a column and
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
  `\Omega = \{0, 1, \ldots, index-1 \}` where `index` is the index of subgroup
  `H` in `G`.

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

``low_index_subgroups(G, N)``: Given a finitely presented group `G = \langle X \mid R \rangle`
(can be a free group), and ``N`` a positive integer, determine the conjugacy classes of
subgroups of ``G`` whose indices is less than or equal to ``N``.

For example to find all subgroups of `G = \langle a, b \mid a^2 = b^3 = (ab)^4 = 1 \rangle`
having index <= 4, can be found as follows:

>>> from sympy.combinatorics.fp_groups import low_index_subgroups
>>> F, a, b = free_group("a, b")
>>> G = FpGroup(F, [a**2, b**3, (a*b)**4])
>>> l = low_index_subgroups(G, 4)
>>> for coset_table in l:
...     print(coset_table.table)
...
[[0, 0, 0, 0]]
[[0, 0, 1, 2], [1, 1, 2, 0], [3, 3, 0, 1], [2, 2, 3, 3]]
[[0, 0, 1, 2], [2, 2, 2, 0], [1, 1, 0, 1]]
[[1, 1, 0, 0], [0, 0, 1, 1]]

This returns the coset tables of subgroups of satisfying the property
that index, `index`, of subgroup in group is `\le n`.

Constructing a presentation for a subgroup
------------------------------------------

In this section we discuss finding the presentation of a subgroup in a finitely
presentation group. While the *subgroup* is currently allowed as input only in
the form of a list of generators for the subgroup, you can expect the
functionality of a *coset table* as input for subgroup in the group in near
future.

There are two ways to construct a set of defining relations for subgroup from
those of ``G``. First is on a set of Schreier generators, known generally as
Reidemeister-Schreier algorithm or on the given list of generators of ``H``.

Reidemeister Schreier algorithm
```````````````````````````````

called using ``reidemeister_presentation(G, Y)`` where ``G`` is the group and
``Y`` is a list of generators for subgroup ``H`` whose presentation we want to
find.

>>> from sympy.combinatorics.fp_groups import reidemeister_presentation
>>> F, x, y = free_group("x, y")
>>> f = FpGroup(F, [x**3, y**5, (x*y)**2])
>>> H = [x*y, x**-1*y**-1*x*y*x]
>>> p1 = reidemeister_presentation(f, H)
>>> p1
((y_1, y_2), (y_1**2, y_2**3, y_2*y_1*y_2*y_1*y_2*y_1))

Bibliography
------------

.. [CDHW73] John J. Cannon, Lucien A. Dimino, George Havas, and Jane M. Watson.
    Implementation and analysis of the Todd-Coxeter algorithm. Math. Comp., 27:463–
    490, 1973.

.. [Ho05] Derek F. Holt,
    Handbook of Computational Group Theory.
    In the series 'Discrete Mathematics and its Applications',
    `Chapman & Hall/CRC 2005, xvi + 514 p <https://www.crcpress.com/Handbook-of-Computational-Group-Theory/Holt-Eick-OBrien/p/book/9781584883722>`_.

.. [Hav91] George Havas, Coset enumeration strategies.
    In Proceedings of the International Symposium on Symbolic and Algebraic Computation (ISSAC'91), Bonn 1991, pages 191--199. ACM Press, 1991.
