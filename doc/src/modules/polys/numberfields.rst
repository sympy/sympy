.. _polys-numberfields:

.. currentmodule:: sympy.polys.numberfields

=============
Number Fields
=============

Introduction
============

  Like many other computations in algebraic number theory, the splitting of
  rational primes can be treated by *rational* methods only. This fact is very
  important if computation by automatic computing machinery is considered. Only
  the knowledge of the irreducible polynomial $f(x)$, a zero of which generates
  the field in question, is needed.

    Olga Taussky, 1953


Concepts like number fields and algebraic numbers are essential to our
understanding of algebraic number theory, but to the computer the subject is
all about polynomials: the ring $\mathbb{Q}[x]$ reduced modulo irreducible
polynomials $f(x) \in \mathbb{Q}[x]$. It thus finds a natural home under the
:py:mod:`~.polys` module in SymPy.

Various authors (such as Taussky, Zimmer, Pohst and Zassenhaus, or Cohen)
have articulated the main goals of computational algebraic number theory in
different ways, but invariably the list centers around a certain essential set
of tasks. As a goal for the ``numberfields`` module in SymPy, we may set the
following list, based on [Cohen93]_ (Sec. 4.9.3).

For a number field $K = \mathbb{Q}(\theta)$, whose ring of algebraic integers
is denoted $\mathbb{Z}_K$, compute:

1. an integral basis of $\mathbb{Z}_K$
2. the decomposition of rational primes in $\mathbb{Z}_K$
3. $\mathfrak{p}$-adic valuations for ideals and elements
4. the Galois group of the Galois closure of $K$
5. a system of fundamental units of $K$
6. the regulator $R(K)$
7. the class number
8. the structure of the class group $Cl(K)$
9. decide whether a given ideal is principal, and if so compute a generator.

As a foundation, and to support our basic ability to define and work with
number fields and algebraic numbers, we also set the following problems:

10. Given an algebraic number -- expressed by radicals and rational operations,
    or even as a special value of a transcendental function -- determine its
    minimal polynomial over $\mathbb{Q}$.
11. Given an irreducible polynomial, distinguish its roots as elements of
    $\mathbb{C}$ by computing isolating intervals.
12. The Primitive Element Problem: Given several algebraic numbers
    $\alpha_1, \ldots, \alpha_m$, compute a single algebraic number $\theta$
    such that $\mathbb{Q}(\alpha_1, \ldots, \alpha_m) = \mathbb{Q}(\theta)$.
13. The Field Isomorphism Problem: Decide whether two number fields
    $\mathbb{Q}(\alpha)$, $\mathbb{Q}(\beta)$ are isomorphic.
14. The Field Membership Problem: Given two algebraic numbers $\alpha$,
    $\beta$, decide whether $\alpha \in \mathbb{Q}(\beta)$, and if so write
    $\alpha = f(\beta)$ for some $f(x) \in \mathbb{Q}[x]$.

At present only a subset of the tasks enumerated in the lists above is yet
supported in SymPy, and if you are interested in expanding support, you are
encouraged to contribute!

At time of writing, the existing solutions to the above problems are found
in the following places:

=================================  ======================================
Task                               Implementation
=================================  ======================================
(1) integral basis                 :py:func:`~.round_two`
(2) prime decomposition            :py:func:`~.prime_decomp`
(3) $\mathfrak{p}$-adic valuation  :py:func:`~.prime_valuation`
(10) find min poly                 :py:func:`~.minpoly`
(11) isolating intervals           :py:func:`~.isolate`
(12) primitive element             :py:func:`~.primitive_element`
(13) field isomorphism             :py:func:`~.field_isomorphism`
(14) field membership              :py:func:`~.to_number_field`
=================================  ======================================


Internals
=========

Algebraic number fields
-----------------------

Algebraic number fields are represented in SymPy by the
:py:class:`~.AlgebraicField` class, which is a part of the
polynomial domains system.



Representing algebraic numbers
------------------------------

There are several different ways to represent algebraic numbers, and different
forms may be preferable for different computational tasks.
See [Cohen93]_ Section 4.2.



As number field elements
````````````````````````

In SymPy, an element of an :py:class:`~.AlgebraicField` is represented by
an :py:class:`~.AlgebraicNumber`.

.. autoclass:: AlgebraicNumber
   :members:


As points in the complex plane
``````````````````````````````

The minimal polynomial for an algebraic number determines it only up to
conjugacy; in other words, the number and all its conjugates are algebraically
indistinguishable.

In order to select a particular conjugate, we can identify
it with a point in the complex plane. This means identifying two real
intervals, one on the real axis and one on the imaginary axis, such that
exactly one root of the given irreducible polynomial over $\mathbb{Q}$ lies
within the complex region so defined.

Support for this in SymPy is currently limited. The :py:func:`~.isolate`
function is currently applicable only to algebraic numbers lying on the
real line.

.. autofunction:: isolate


As elements of finitely-generated modules
`````````````````````````````````````````

In computational algebraic number theory, finitely-generated
$\mathbb{Z}$-modules are of central importance. For example, every
order_ and every ideal_ is such a module.

In particular, the maximal order -- or `ring of integers`_ -- in a number field
is a finitely-generated $\mathbb{Z}$-module, whose generators form an
`integral basis`_ for the field.

Classes allowing us to represent such modules, and their elements, and
homomorphisms between them, are provided in the
:py:mod:`~.modules` module.


Modules
-------

.. automodule:: sympy.polys.numberfields.modules

Class Reference
```````````````
.. autoclass:: Module
   :members:

   .. automethod:: Module.__call__

.. autoclass:: PowerBasis
   :members:

   .. automethod:: PowerBasis.__init__

.. autoclass:: Submodule
   :members:

   .. automethod:: Submodule.__init__

.. autoclass:: ModuleElement
   :members:

   .. automethod:: ModuleElement.__init__
   .. automethod:: ModuleElement.__add__
   .. automethod:: ModuleElement.__mul__
   .. automethod:: ModuleElement.__mod__

.. autoclass:: PowerBasisElement
   :members:

.. autofunction:: make_mod_elt

.. autoclass:: ModuleHomomorphism
   :members:

   .. automethod:: ModuleHomomorphism.__init__

.. autoclass:: ModuleEndomorphism
   :members:

   .. automethod:: ModuleEndomorphism.__init__

.. autoclass:: InnerEndomorphism
   :members:

   .. automethod:: InnerEndomorphism.__init__

.. autoclass:: EndomorphismRing
   :members:

   .. automethod:: EndomorphismRing.__init__

.. autofunction:: find_min_poly


Utilities
---------

.. currentmodule:: sympy.polys.numberfields.utilities

.. autofunction:: is_rat
.. autofunction:: is_int
.. autofunction:: get_num_denom
.. autofunction:: extract_fundamental_discriminant

.. autoclass:: AlgIntPowers
   :members:

   .. automethod:: AlgIntPowers.__init__

.. autofunction:: coeff_search
.. autofunction:: supplement_a_subspace


Solving the Main Problems
=========================

Integral Basis
--------------

.. currentmodule:: sympy.polys.numberfields
.. autofunction:: round_two


Prime Decomposition
-------------------

.. autofunction:: prime_decomp
.. currentmodule:: sympy.polys.numberfields.primes
.. autoclass:: PrimeIdeal
   :members:

   .. automethod:: PrimeIdeal.__init__


p-adic Valuation
----------------

.. currentmodule:: sympy.polys.numberfields
.. autofunction:: prime_valuation


Finding Minimal Polynomials
---------------------------

.. autofunction:: minimal_polynomial
.. autofunction:: minpoly


The Primitive Element Problem
-----------------------------

.. autofunction:: primitive_element


The Field Isomorphism Problem
-----------------------------

.. autofunction:: field_isomorphism


The Field Membership Problem
----------------------------

.. autofunction:: to_number_field



.. _ideal: https://en.wikipedia.org/wiki/Ideal_(ring_theory)
.. _order: https://en.wikipedia.org/wiki/Order_(ring_theory)
.. _ring of integers: https://en.wikipedia.org/wiki/Ring_of_integers
.. _integral basis: https://en.wikipedia.org/wiki/Algebraic_number_field#Integral_basis

