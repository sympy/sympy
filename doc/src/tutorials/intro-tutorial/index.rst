.. _intro-tutorial:

=====================
Introductory Tutorial
=====================

This tutorial aims to give an introduction to SymPy for someone who has not
used the library before.  Many features of SymPy will be introduced in this
tutorial, but they will not be exhaustive. In fact, virtually every
functionality shown in this tutorial will have more options or capabilities
than what will be shown.  The rest of the SymPy documentation serves as API
documentation, which extensively lists every feature and option of each
function.

These are the goals of this tutorial:

.. NB: This is mainly here for you, the person who is editing and adding to
   this tutorial. Try to keep these principles in mind.

- To give a guide, suitable for someone who has never used SymPy (but who has
  used Python and knows the necessary mathematics).

- To be written in a narrative format, which is both easy and fun to follow.
  It should read like a book.

- To give insightful examples and exercises, to help the reader learn and to
  make it entertaining to work through.

- To introduce concepts in a logical order.

.. In other words, don't try to get ahead of yourself.

- To use good practices and idioms, and avoid antipatterns.  Functions or
  methodologies that tend to lead to antipatterns are avoided. Features that
  are only useful to advanced users are not shown.

- To be consistent.  If there are multiple ways to do it, only the best way is
  shown.

.. For example, there are at least five different ways to create Symbols.
   ``symbols`` is the only one that is general and doesn't lead to
   antipatterns, so it is the only one used.

- To avoid unnecessary duplication, it is assumed that previous sections of
  the tutorial have already been read.

Feedback on this tutorial, or on SymPy in general is always welcome. Just
write to our `mailing list
<https://groups.google.com/forum/?fromgroups#!forum/sympy>`_.


**Content**

.. toctree::
   :maxdepth: 2

   preliminaries.rst
   intro.rst
   gotchas.rst
   features.rst
   next.rst
