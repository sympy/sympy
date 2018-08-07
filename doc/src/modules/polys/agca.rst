.. _polys-agca:

========================================================
AGCA - Algebraic Geometry and Commutative Algebra Module
========================================================

Introduction
============

  Algebraic geometry is a mixture of the ideas of two Mediterranean
  cultures. It is the superposition of the Arab science of the lightening
  calculation of the solutions of equations over the Greek art of position
  and shape.
  This tapestry was originally woven on European soil and is still being refined
  under the influence of international fashion. Algebraic geometry studies the
  delicate balance between the geometrically plausible and the algebraically
  possible.  Whenever one side of this mathematical teeter-totter outweighs the
  other, one immediately loses interest and runs off in search of a more exciting
  amusement.

    George R. Kempf
    1944 -- 2002


Algebraic Geometry refers to the study of geometric problems via algebraic
methods (and sometimes vice versa). While this is a rather old topic,
algebraic geometry as understood today is very much a 20th century
development. Building on ideas of e.g. Riemann and Dedekind, it was realized
that there is an intimate connection between properties of the set of
solutions of a system of polynomial equations (called an algebraic variety)
and the behavior of the set of polynomial functions on that variety
(called the coordinate ring).

As in many geometric disciplines, we can distinguish between local and global
questions (and methods). Local investigations in algebraic geometry are
essentially equivalent to the study of certain rings, their ideals and modules.
This latter topic is also called commutative algebra. It is the basic local
toolset of algebraic geometers, in much the same way that differential analysis
is the local toolset of differential geometers.

A good conceptual introduction to commutative algebra is [Atiyah69]_. An
introduction more geared towards computations, and the work most of the
algorithms in this module are based on, is [Greuel2008]_.

This module aims to eventually allow expression and solution of both
local and global geometric problems, both in the classical case over a field
and in the more modern arithmetic cases. So far, however, there is no geometric
functionality at all. Currently the module only provides tools for computational
commutative algebra over fields.

All code examples assume::

    >>> from sympy import *
    >>> x, y, z = symbols('x,y,z')
    >>> init_printing(use_unicode=True, wrap_line=False, no_global=True)

Reference
=========

In this section we document the usage of the AGCA module. For convenience of
the reader, some definitions and examples/explanations are interspersed.

Base Rings
----------

Almost all computations in commutative algebra are relative to a "base ring".
(For example, when asking questions about an ideal, the base ring is the ring
the ideal is a subset of.) In principle all polys "domains" can be used as base
rings. However, useful functionality is only implemented for polynomial rings
over fields, and various localizations and quotients thereof.

As demonstrated in
the examples below, the most convenient method to create objects you are
interested in is to build them up from the ground field, and then use the
various methods to create new objects from old. For example, in order to
create the local ring of the nodal cubic `y^2 = x^3` at the origin, over
`\mathbb{Q}`, you do::

    >>> lr = QQ.old_poly_ring(x, y, order="ilex") / [y**2 - x**3]
    >>> lr
    ℚ[x, y, order=ilex]
    ───────────────────
        ╱   3    2╲
        ╲- x  + y ╱

Note how the python list notation can be used as a short cut to express ideals.
You can use the ``convert`` method to return ordinary sympy objects into
objects understood by the AGCA module (although in many cases this will be done
automatically -- for example the list was automatically turned into an ideal,
and in the process the symbols `x` and `y` were automatically converted into
other representations). For example::

    >>> X, Y = lr.convert(x), lr.convert(y) ; X
        ╱   3    2╲
    x + ╲- x  + y ╱

    >>> x**3 == y**2
    False

    >>> X**3 == Y**2
    True

When no localisation is needed, a more mathematical notation can be
used. For example, let us create the coordinate ring of three-dimensional
affine space `\mathbb{A}^3`::

    >>> ar = QQ.old_poly_ring(x, y, z); ar
    ℚ[x, y, z]

For more details, refer to the following class documentation. Note that
the base rings, being domains, are the main point of overlap between the
AGCA module and the rest of the polys module. All domains are documented
in detail in the polys reference, so we show here only an abridged version,
with the methods most pertinent to the AGCA module.

.. autoclass:: sympy.polys.domains.ring.Ring
   :members: free_module, ideal, quotient_ring
   :noindex:

.. autofunction:: sympy.polys.domains.polynomialring.PolynomialRing
   :noindex:

.. autoclass:: sympy.polys.domains.quotientring.QuotientRing
   :noindex:

Modules, Ideals and their Elementary Properties
-----------------------------------------------

Let `A` be a ring. An `A`-module is a set `M`, together with two binary
operations `+: M \times M \to M` and `\times: R \times M \to M` called
addition and scalar multiplication. These are required to satisfy certain
axioms, which can be found in e.g. [Atiyah69]_. In this way modules are
a direct generalisation of both vector spaces (`A` being a field) and abelian
groups (`A = \mathbb{Z}`). A *submodule* of the `A`-module `M` is a subset
`N \subset M`, such that the binary operations restrict to `N`, and `N` becomes
an `A`-module with these operations.

The ring `A` itself has a natural `A`-module structure where addition and
multiplication in the module coincide with addition and multiplication in
the ring. This `A`-module is also written as `A`. An `A`-submodule of `A`
is called an *ideal* of `A`. Ideals come up very naturally in algebraic
geometry. More general modules can be seen as a technically convenient "elbow
room" beyond talking only about ideals.

If `M`, `N` are `A`-modules,
then there is a natural (componentwise) `A`-module structure on `M \times N`.
Similarly there are `A`-module structures on cartesian products of more
components. (For the categorically inclined:
the cartesian product of finitely many `A`-modules, with this
`A`-module structure, is the finite biproduct in the category of all
`A`-modules. With infinitely many components, it is the direct product
(but the infinite direct sum has to be constructed differently).) As usual,
repeated product of the `A`-module `M` is denoted `M, M^2, M^3 \ldots`, or
`M^I` for arbitrary index sets `I`.

An `A`-module `M` is called *free* if it is isomorphic to the `A`-module
`A^I` for some (not necessarily finite) index set `I` (refer to the next
section for a definition of isomorphism). The cardinality of `I` is called
the *rank* of `M`; one may prove this is well-defined.
In general, the AGCA module only works with free modules of finite rank, and
other closely related modules. The easiest way to create modules is to use
member methods of the objects they are made up from. For example, let us create
a free module of rank 4 over the coordinate ring of `\mathbb{A}^2`
we created above, together with a submodule::

    >>> F = ar.free_module(4) ; F
              4
    ℚ[x, y, z]

    >>> S = F.submodule([1, x, x**2, x**3], [0, 1, 0, y]) ; S
    ╱⎡       2   3⎤              ╲
    ╲⎣1, x, x , x ⎦, [0, 1, 0, y]╱

Note how python lists can be used as a short-cut notation for module
elements (vectors). As usual, the ``convert`` method can be used to convert
sympy/python objects into the internal AGCA representation (see detailed
reference below).

Here is the detailed documentation of the classes for modules, free modules,
and submodules:

.. currentmodule:: sympy.polys.agca.modules

.. autoclass:: Module
   :members:

.. autoclass:: FreeModule
   :members:

.. autoclass:: SubModule
   :members:


Ideals are created very similarly to modules. For example, let's verify
that the nodal cubic is indeed singular at the origin::

    >>> I = lr.ideal(x, y)
    >>> I == lr.ideal(x)
    False

    >>> I == lr.ideal(y)
    False

We are using here the fact that a curve is non-singular at a point if and only
if the maximal ideal of the local ring is principal, and that in this case at
least one of `x` and `y` must be generators.

This is the detailed documentation of the class ideal. Please note that most
of the methods regarding properties of ideals (primality etc.) are not yet
implemented.

.. currentmodule:: sympy.polys.agca.ideals

.. autoclass:: Ideal
   :members:


If `M` is an `A`-module and `N` is an `A`-submodule, we can define two elements
`x` and `y` of `M` to be equivalent if `x - y \in N`. The set of equivalence
classes is written `M/N`, and has a natural `A`-module structure. This is
called the quotient module of `M` by `N`. If `K` is a submodule of `M`
containing `N`, then `K/N` is in a natural way a submodule of `M/N`. Such a
module is called a subquotient. Here is the documentation of quotient and
subquotient modules:

.. currentmodule:: sympy.polys.agca.modules

.. autoclass:: QuotientModule
   :members:

.. autoclass:: SubQuotientModule
   :members:

Module Homomorphisms and Syzygies
---------------------------------

Let `M` and `N` be `A`-modules. A mapping `f: M \to N` satisfying various
obvious properties (see [Atiyah69]_) is called an `A`-module homomorphism.
In this case `M` is called the *domain* and *N* the *codomain*. The
set `\{x \in M | f(x) = 0\}` is called the *kernel* `ker(f)`, whereas the
set `\{f(x) | x \in M\}` is called the *image* `im(f)`.
The kernel is a submodule of `M`, the image is a submodule of `N`.
The homomorphism `f` is injective if and only if `ker(f) = 0` and surjective
if and only if `im(f) = N`.
A bijective homomorphism is called an *isomorphism*. Equivalently, `ker(f) = 0`
and `im(f) = N`. (A related notion, which currently has no special name in
the AGCA module, is that of the *cokernel*, `coker(f) = N/im(f)`.)

Suppose now `M` is an `A`-module. `M` is called *finitely generated* if there
exists a surjective homomorphism `A^n \to M` for some `n`. If such a morphism
`f` is chosen, the images of the standard basis of `A^n` are called the
*generators* of `M`. The module `ker(f)` is called *syzygy module* with respect
to the generators. A module is called *finitely presented* if it is finitely
generated with a finitely generated syzygy module. The class of finitely
presented modules is essentially the largest class we can hope to be able to
meaningfully compute in.

It is an important theorem that, for all the rings we are considering,
all submodules of finitely generated modules are finitely generated, and hence
finitely generated and finitely presented modules are the same.

The notion of syzygies, while it may first seem rather abstract, is actually
very computational. This is because there exist (fairly easy) algorithms for
computing them, and more general questions (kernels, intersections, ...) are
often reduced to syzygy computation.

Let us say a few words about the definition of homomorphisms in the AGCA
module. Suppose first that `f : M \to N` is an arbitrary morphism of
`A`-modules. Then if `K` is a submodule of `M`, `f` naturally defines a new
homomorphism `g: K \to N` (via `g(x) = f(x)`), called the *restriction* of
`f` to `K`. If now `K` contained in the kernel of
`f`, then moreover `f` defines in a natural homomorphism `g: M/K \to N`
(same formula as above!), and we say that `f` *descends* to `M/K`.
Similarly, if `L` is a submodule of `N`, there is a natural homomorphism
`g: M \to N/L`, we say that `g` *factors* through `f`. Finally, if now `L`
contains the image of `f`, then there is a natural homomorphism `g: M \to L`
(defined, again, by the same formula), and we say `g` is obtained from `f`
by restriction of codomain. Observe also that each of these four operations
is reversible, in the sense that given `g`, one can always (non-uniquely)
find `f` such that `g` is obtained from `f` in the above way.

Note that all modules implemented in AGCA are obtained from free modules by
taking a succession of submodules and quotients. Hence, in order to explain
how to define a homomorphism between arbitrary modules, in light of the above,
we need only explain how to define homomorphisms of free modules. But,
essentially by the definition of free module, a homomorphism from a free module
`A^n` to any module `M` is precisely the same as giving `n` elements of `M`
(the images of the standard basis), and giving an element of a free module
`A^m` is precisely the same as giving `m` elements of `A`. Hence a
homomorphism of free modules `A^n \to A^m` can be specified via a matrix,
entirely analogously to the case of vector spaces.

The functions ``restrict_domain`` etc. of the class ``Homomorphism`` can be
used
to carry out the operations described above, and homomorphisms of free modules
can in principle be instantiated by hand. Since these operations are so common,
there is a convenience function ``homomorphism`` to define a homomorphism
between arbitrary modules via the method outlined above. It is essentially
the only way homomorphisms need ever be created by the user.

.. currentmodule:: sympy.polys.agca.homomorphisms
.. autofunction:: homomorphism

Finally, here is the detailed reference of the actual homomorphism class:

.. autoclass:: ModuleHomomorphism
   :members:
