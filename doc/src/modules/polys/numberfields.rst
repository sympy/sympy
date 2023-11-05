.. _polys-numberfields:

.. currentmodule:: sympy.polys.numberfields

=============
Number Fields
=============

Introduction
============

.. epigraph::
   Like many other computations in algebraic number theory, the splitting of
   rational primes can be treated by *rational* methods only. This fact is very
   important if computation by automatic computing machinery is considered. Only
   the knowledge of the irreducible polynomial $f(x)$, a zero of which generates
   the field in question, is needed.

   -- Olga Taussky, 1953


Concepts like number fields and algebraic numbers are essential to our
understanding of algebraic number theory, but to the computer the subject is
all about polynomials: the ring $\mathbb{Q}[x]$ reduced modulo irreducible
polynomials $f(x) \in \mathbb{Q}[x]$. It thus finds a natural home under the
:py:mod:`~.polys` module in SymPy.

Various authors (such as Taussky, Zimmer, Pohst and Zassenhaus, or Cohen)
have articulated the main goals of computational algebraic number theory in
different ways, but invariably the list centers around a certain essential set
of tasks. As a goal for the ``numberfields`` module in SymPy, we may set the
following list, based on [Cohen93]_, Sec. 4.9.3.

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
number fields and algebraic numbers, we also set the following problems,
following [Cohen93]_, Sec. 4.5.

10. Given an algebraic number -- expressed by radicals and rational operations,
    or even as a special value of a transcendental function -- determine its
    minimal polynomial over $\mathbb{Q}$.
11. The Subfield Problem: Given two number fields $\mathbb{Q}(\alpha)$,
    $\mathbb{Q}(\beta)$ via the minimal polynomials for their generators
    $\alpha$ and $\beta$, decide whether one field is isomorphic to a subfield
    of the other, and if so exhibit an embedding.
12. The Field Membership Problem: Given two algebraic numbers $\alpha$,
    $\beta$, decide whether $\alpha \in \mathbb{Q}(\beta)$, and if so write
    $\alpha = f(\beta)$ for some $f(x) \in \mathbb{Q}[x]$.
13. The Primitive Element Problem: Given several algebraic numbers
    $\alpha_1, \ldots, \alpha_m$, compute a single algebraic number $\theta$
    such that $\mathbb{Q}(\alpha_1, \ldots, \alpha_m) = \mathbb{Q}(\theta)$.

At present only a subset of the tasks enumerated above is yet supported in
SymPy, and if you are interested in expanding support, you are encouraged to
contribute! An excellent source, providing solutions to all the remaining
problems (as well as those already solved) is [Cohen93]_.

At time of writing, the existing solutions to the above problems are found
in the following places:

=================================  ======================================
Task                               Implementation
=================================  ======================================
(1) integral basis                 :py:func:`~.round_two`
(2) prime decomposition            :py:func:`~.prime_decomp`
(3) $\mathfrak{p}$-adic valuation  :py:func:`~.prime_valuation`
(4) Galois group                   :py:func:`~.galois_group`
(10) find minimal polynomial       :py:func:`~.minimal_polynomial`
(11) subfield                      :py:func:`~.field_isomorphism`
(12) field membership              :py:func:`~.to_number_field`
(13) primitive element             :py:func:`~.primitive_element`
=================================  ======================================


Solving the Main Problems
=========================

Integral Basis
--------------
.. _IntegralBasis:

.. currentmodule:: sympy.polys.numberfields.basis
.. autofunction:: round_two


Prime Decomposition
-------------------
.. _PrimeDecomposition:

.. currentmodule:: sympy.polys.numberfields.primes
.. autofunction:: prime_decomp
.. autoclass:: PrimeIdeal
   :members:

   .. automethod:: PrimeIdeal.__init__
   .. automethod:: PrimeIdeal.__add__
   .. automethod:: PrimeIdeal.__mul__


p-adic Valuation
----------------
.. _pAdicValuation:

.. currentmodule:: sympy.polys.numberfields.primes
.. autofunction:: prime_valuation


Galois Groups
-------------
.. _GaloisGroups:

.. currentmodule:: sympy.polys.numberfields.galoisgroups
.. autofunction:: galois_group


Finding Minimal Polynomials
---------------------------
.. _MinimalPolynomials:

.. currentmodule:: sympy.polys.numberfields.minpoly
.. autofunction:: minimal_polynomial
.. autofunction:: minpoly


The Subfield Problem
--------------------
.. _SubfieldProblem:

.. automodule:: sympy.polys.numberfields.subfield

.. autofunction:: field_isomorphism

.. autofunction:: primitive_element

.. autofunction:: to_number_field



Internals
=========

Algebraic number fields
-----------------------

Algebraic number fields are represented in SymPy by the
:py:class:`~.AlgebraicField` class, which is a part of
:ref:`the polynomial domains system<polys-domainsref>`.


Representing algebraic numbers
------------------------------

There are several different ways to represent algebraic numbers, and different
forms may be preferable for different computational tasks.
See [Cohen93]_, Sec. 4.2.


As number field elements
````````````````````````

In SymPy, there is a distinction between number and expression classes defined
in the :py:mod:`sympy.core.numbers` module on the one hand, and domains and
domain elements defined in the :py:mod:`~sympy.polys` module on the other.
This is explained in more detail :ref:`here<polys-domainsintro>`.

When it comes to algebraic numbers, the :py:mod:`sympy.core.numbers` module
offers the :py:class:`~.AlgebraicNumber` class, while the
:py:mod:`~sympy.polys` module offers the
:py:class:`~sympy.polys.polyclasses.ANP` class. This is the type of domain
elements belonging to the :py:class:`~.AlgebraicField` domain.


As elements of finitely-generated modules
`````````````````````````````````````````

In computational algebraic number theory, finitely-generated
$\mathbb{Z}$-modules are of central importance. For example, every
order_ and every ideal_ is such a module.

In particular, the maximal order -- or `ring of integers`_ -- in a number field
is a finitely-generated $\mathbb{Z}$-module, whose generators form an
`integral basis`_ for the field.

Classes allowing us to represent such modules, and their elements, are provided
in the :py:mod:`~.modules` module. Here, the :py:class:`~.ModuleElement` class
provides another way to represent algebraic numbers.


Finitely-generated modules
--------------------------

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

.. autofunction:: isolate



.. _ideal: https://en.wikipedia.org/wiki/Ideal_(ring_theory)
.. _order: https://en.wikipedia.org/wiki/Order_(ring_theory)
.. _ring of integers: https://en.wikipedia.org/wiki/Ring_of_integers
.. _integral basis: https://en.wikipedia.org/wiki/Algebraic_number_field#Integral_basis

