Polycyclic Groups
=================

Introduction
------------

This module presents the functionality designed for computing with
polycyclic groups(PcGroup for short). The name of the corresponding
SymPy object is ``PolycyclicGroup``. The functions or classes described
here are studied under **Computational Group Theory**.

Overview of functionalities
```````````````````````````

* The construction of PolycyclicGroup from a given PermutationGroup.

* Computation of polycyclic generating sequence(pcgs for short) and
  polycyclic series(pc_series).

* Computation of relative order for polycyclic series.

* Implementation of class Collector which can be treated as a base for polycylic groups.

* Implementation of polycyclic group presentation(pc_presentation for short).

* Computation of exponent vector, depth and leading exponent for a given element
  of a polycyclic group.

For a description of fundamental algorithms of polycyclic groups, we
often make use of *Handbook of Computational Group Theory*.


The Construction of Polycyclic Groups
-------------------------------------

Given a Permutation Group, A Polycyclic Group is constructed by computing the
corresponding polycylic generating sequence, polycyclic series and it's
relative order.

Attributes of PolycyclicGroup
`````````````````````````````

* ``pc_sequence`` : Polycyclic sequence is formed by collecting all the missing
  generators between the adjacent groups in the derived series of given
  permutation group.

* ``pc_series`` : Polycyclic series is formed by adding all the missing generators
  of ``der[i+1]`` in ``der[i]``, where ``der`` represents derived series.

* ``relative_order`` : A list, computed by the ratio of adjacent groups in pc_series.

* ``collector`` : By default, it is None. Collector class provides the polycyclic presentation.

>>> from sympy.combinatorics.named_groups import SymmetricGroup
>>> G = SymmetricGroup(4)
>>> PcGroup = G.polycyclic_group()
>>> len(PcGroup.pcgs)
4
>>> pc_series = PcGroup.pc_series
>>> pc_series[0].equals(G)  # use equals, not literal `==`
True
>>> gen = pc_series[len(pc_series) - 1].generators[0]
>>> gen.is_identity
True
>>> PcGroup.relative_order
[2, 3, 2, 2]


The Construction of Collector
-----------------------------

Collector is one of the attributes of class PolycyclicGroup.

Attributes of Collector
```````````````````````

Collector posses all the attributes of PolycyclicGroup, In addition there are
few more attributes which are defined below:

* ``free_group`` : free_group provides the mapping of polycyclic generating sequence with
  the free group elements.

* ``pc_presentation`` : Provides the presentation of polycyclic groups with the
  help of power and conjugate relators.

>>> from sympy.combinatorics.named_groups import SymmetricGroup
>>> G = SymmetricGroup(3)
>>> PcGroup = G.polycyclic_group()
>>> Collector = PcGroup.collector
>>> Collector.free_group
<free group on the generators (x0, x1)>
>>> Collector.pc_presentation
{x0**2: (), x1**3: (), x0**-1*x1*x0: x1**2}


Computation of Minimal Uncollected Subword
``````````````````````````````````````````

A word ``V`` defined on generators in the free_group of pc_group is a minimal
uncollected subword of the word ``W`` if ``V`` is a subword of ``W`` and it has one of
the following form:

* `v = {x_{i+1}}^{a_j}x_i`
* `v = {x_{i+1}}^{a_j}{x_i}^{-1}`
* `v = {x_i}^{a_j}`

`a_j \notin \{0, \ldots \mathrm{relative\_order}[j]-1\}`.

>>> from sympy.combinatorics.named_groups import SymmetricGroup
>>> from sympy.combinatorics.free_groups import free_group
>>> G = SymmetricGroup(4)
>>> PcGroup = G.polycyclic_group()
>>> collector = PcGroup.collector
>>> F, x1, x2 = free_group("x1, x2")
>>> word = x2**2*x1**7
>>> collector.minimal_uncollected_subword(word)
((x2, 2),)


Computation of Subword Index
````````````````````````````

For a given word and it's subword, subword_index computes the
starting and ending index of the subword in the word.

>>> from sympy.combinatorics.named_groups import SymmetricGroup
>>> from sympy.combinatorics.free_groups import free_group
>>> G = SymmetricGroup(4)
>>> PcGroup = G.polycyclic_group()
>>> collector = PcGroup.collector
>>> F, x1, x2 = free_group("x1, x2")
>>> word = x2**2*x1**7
>>> w = x2**2*x1
>>> collector.subword_index(word, w)
(0, 3)
>>> w = x1**7
>>> collector.subword_index(word, w)
(2, 9)


Computation of Collected Word
`````````````````````````````

A word ``W`` is called collected, if ``W`` `= {x_{i_1}}^{a_1} \ldots {x_{i_r}}^{a_r}`
with `i_1 < i_2< \ldots < i_r` and `a_j` is in `\{1 \ldots s_{j-1}\}`, where `s_j`
represents the respective relative order.

>>> from sympy.combinatorics.named_groups import SymmetricGroup
>>> from sympy.combinatorics.perm_groups import PermutationGroup
>>> from sympy.combinatorics.free_groups import free_group
>>> G = SymmetricGroup(4)
>>> PcGroup = G.polycyclic_group()
>>> collector = PcGroup.collector
>>> F, x0, x1, x2, x3 = free_group("x0, x1, x2, x3")
>>> word = x3*x2*x1*x0
>>> collected_word = collector.collected_word(word)
>>> free_to_perm = {}
>>> free_group = collector.free_group
>>> for sym, gen in zip(free_group.symbols, collector.pcgs):
...     free_to_perm[sym] = gen
>>> G1 = PermutationGroup()
>>> for w in word:
...     sym = w[0]
...     perm = free_to_perm[sym]
...     G1 = PermutationGroup([perm] + G1.generators)
>>> G2 = PermutationGroup()
>>> for w in collected_word:
...     sym = w[0]
...     perm = free_to_perm[sym]
...     G2 = PermutationGroup([perm] + G2.generators)

The two are not identical but they are equivalent:

>>> G1 == G2
False
>>> G1.equals(G2)
True


Computation of Polycyclic Presentation
--------------------------------------

The computation of presentation starts from the bottom of the pcgs and polycyclic series.
Storing all the previous generators from pcgs and then taking the last generator
as the generator which acts as a conjugator and conjugates all the previous
generators in the list.

To get a clear picture, start with an example of SymmetricGroup(4). For S(4) there are 4
generators in pcgs say `[x_0, x_1, x_2, x_3]` and the relative_order vector is [2, 3, 2, 2].
Starting from bottom of this sequence the presentation is computed in order as below.

using only `[x_3]` from ``pcgs`` and ``pc_series[4]`` compute:

* `x_3^2`

using only `[x_3]` from ``pcgs`` and ``pc_series[3]`` compute:

* `x_2^2`
* `x_2^{-1}x_3x_2`

using `[x_3, x_2]` from ``pcgs`` and ``pc_series[2]`` compute:

* `x_1^3`
* `x_1^{-1}x_3x_1`
* `x_1^{-1}x_2x_1`

using `[x_3, x_2, x_1]` from ``pcgs`` and ``pc_series[1]`` compute:

* `x_0^2`
* `x_0^{-1}x_3x_0`
* `x_0^{-1}x_2x_0`
* `x_0^{-1}x_1x_0`


One thing to note is same group can have different pcgs due to variying derived_series which,
results in different polycyclic presentations.

>>> from sympy.combinatorics.named_groups import SymmetricGroup
>>> from sympy.combinatorics.permutations import Permutation
>>> G = SymmetricGroup(4)
>>> PcGroup = G.polycyclic_group()
>>> collector = PcGroup.collector
>>> pcgs = PcGroup.pcgs
>>> len(pcgs)
4
>>> free_group = collector.free_group
>>> pc_resentation = collector.pc_presentation
>>> free_to_perm = {}
>>> for s, g in zip(free_group.symbols, pcgs):
...     free_to_perm[s] = g
>>> for k, v in pc_resentation.items():
...     k_array = k.array_form
...     if v != ():
...        v_array = v.array_form
...     lhs = Permutation()
...     for gen in k_array:
...         s = gen[0]
...         e = gen[1]
...         lhs = lhs*free_to_perm[s]**e
...     if v == ():
...         assert lhs.is_identity
...         continue
...     rhs = Permutation()
...     for gen in v_array:
...         s = gen[0]
...         e = gen[1]
...         rhs = rhs*free_to_perm[s]**e
...     assert lhs == rhs


Computation of Exponent Vector
------------------------------

Any generator of the polycyclic group can be represented with the help of it's
polycyclic generating sequence. Hence, the length of exponent vector is equal to
the length of the pcgs.

A given generator ``g`` of the polycyclic group, can be represented as
`g = x_1^{e_1} \ldots x_n^{e_n}`, where `x_i` represents polycyclic generators
and ``n`` is the number of generators in the free_group equal to the length of pcgs.

>>> from sympy.combinatorics.named_groups import SymmetricGroup
>>> from sympy.combinatorics.permutations import Permutation
>>> G = SymmetricGroup(4)
>>> PcGroup = G.polycyclic_group()
>>> collector = PcGroup.collector
>>> pcgs = PcGroup.pcgs
>>> collector.exponent_vector(G[0])
[1, 0, 0, 0]
>>> exp = collector.exponent_vector(G[1])
>>> g = Permutation()
>>> for i in range(len(exp)):
...     g = g*pcgs[i]**exp[i] if exp[i] else g
>>> assert g == G[1]


Depth of Polycyclic generator
`````````````````````````````

Depth of a given polycyclic generator is defined as the index of the first
non-zero entry in the exponent vector.

>>> from sympy.combinatorics.named_groups import SymmetricGroup
>>> G = SymmetricGroup(3)
>>> PcGroup = G.polycyclic_group()
>>> collector = PcGroup.collector
>>> collector.depth(G[0])
2
>>> collector.depth(G[1])
1


Computation of Leading Exponent
```````````````````````````````

Leading exponent represents the exponent of polycyclic generator at the above depth.

>>> from sympy.combinatorics.named_groups import SymmetricGroup
>>> G = SymmetricGroup(3)
>>> PcGroup = G.polycyclic_group()
>>> collector = PcGroup.collector
>>> collector.leading_exponent(G[1])
1


Bibliography
------------

.. [Ho05] Derek F. Holt,
    Handbook of Computational Group Theory.
    In the series 'Discrete Mathematics and its Applications',
    `Chapman & Hall/CRC 2005, xvi + 514 p <https://www.routledge.com/Handbook-of-Computational-Group-Theory/Holt-Eick-OBrien/p/book/9781584883722>`_.
