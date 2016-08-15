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
we make use of *Handbook of Computational Group Theory*.

The Construction of Finitely Presented Groups
---------------------------------------------

Finitely presented groups are construced by factoring a free group by a
set of relators.

Free Groups and word
====================

Construction of a Free Group
----------------------------

Construction of words
---------------------

Methods for manipulating words
------------------------------

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
