==============================
Philosophy behind unit systems
==============================

Dimensions
==========

Introduction
------------

At the root of unit systems are dimension systems, whose structure mainly
determines the one of unit systems. Our definition could seem rough but they
are largely sufficient for our purposes.

A dimension will be defined as a property which is measurable and assigned to
a specific phenomenon. In this sense dimensions are different from pure numbers
because they carry some extra-sense, and for this reason two different
dimensions cannot be added. For example time or length are dimensions, but
also any other things which has some sense for us, like angle, number of
particles (moles...) or information (bits...).

From this point of view the only truly dimensionless quantity are pure
numbers. The idea of being dimensionless is very system-dependent, as can be
seen from the :math:`(c, \hbar, G)`, in which all units appears to be
dimensionless in the usual common sense. This is unavoidable for computability
of generic unit systems (but at the end we can tell the program what is
dimensionless).

Dimensions can be composed together by taking their product or their ratio (to
be defined below). For example the velocity is defined as length divided by
time, or we can see the length as velocity multiplied by time, depending of
what we see as the more fundamental: in general we can select a set of base
dimensions from which we can describe all the others.


Group structure
---------------

After this short introduction whose aim was to introduce the dimensions from
an intuitive perspective, we describe the mathematical structure. A dimension
system with :math:`n` independent dimensions :math:`\{d_i\}_{i=1,\ldots,n}` is
described by a multiplicative group :math:`G`:

- there an identity element :math:`1` corresponding to pure numbers;
- the product :math:`D_3 = D_1 D_2` of two elements :math:`D_1, D_2 \in G`
  is also in :math:`G`;
- any element :math:`D \in G` has an inverse :math:`D^{-1} \in G`.

We denote

.. math::

    D^n = \underbrace{D \times \cdots \times D}_{\text{$n$ times}},

and by definition :math:`D^0 = 1`. The :math:`\{d_i\}_{i=1,\ldots,n}` are
called generators of the group since any element :math:`D \in G` can be
expressed as the product of powers of the generators:

.. math::

    D = \prod_{i=1}^n d_i^{a_i}, \qquad
    a_i \in \mathbf{Z}.

The identity is given for :math:`a_i = 0, \forall i`, while we recover the
generator :math:`d_i` for `a_i = 1, a_j = 0, \forall j \neq i`. This group has
the following properties:

1. abelian, since the generator commutes, :math:`[d_i, d_j] =  0`;
2. countable (infinite but discrete) since the elements are indexed by the
   powers of the generators [#]_.

One can change the dimension basis :math:`\{d'_i\}_{i=1,\ldots,n}` by taking
some combination of the old generators:

.. math::

    d'_i = \prod_{j=1}^n d_j^{P_{ij}}.


Linear space representation
---------------------------

It is possible to use the linear space :math:`\mathbf{Z}^n` as a representation
of the group since the power coefficients :math:`a_i` carry all the
information one needs (we do not distinguish between the element of the group
and its representation):

.. math::

    (d_i)_j = \delta_{ij}, \qquad
    D =
    \begin{pmatrix}
    a_1 \\ \vdots \\ a_n
    \end{pmatrix}.

The change of basis to :math:`d'_i` follows the usual rule of change of basis
for linear space, the matrix being given by the coefficients
:math:`P_{ij}`, which are simply the coefficients of the new vectors in
term of the old basis:

.. math::

    d'_i = P_{ij} d_j.

We will use this last solution in our algorithm.


An example
----------

In order to illustrate all this formalism, we end this section with a specific
example, the MKS system (m, kg, s) with dimensions (L: length, M: mass,
T: time). They are represented as (we will always sort the vectors in
alphabetic order)

.. math::

    L =
    \begin{pmatrix}
    1 \\ 0 \\ 0
    \end{pmatrix}, \qquad
    M =
    \begin{pmatrix}
    0 \\ 1 \\ 0
    \end{pmatrix}, \qquad
    T =
    \begin{pmatrix}
    0 \\ 0 \\ 1
    \end{pmatrix}.

Other dimensions can be derived, for example velocity :math:`V` or action
:math:`A`

.. math::

    V = L T^{-1},  \qquad
    A = M L^2 T^{-2},\\
    V =
    \begin{pmatrix}
    1 \\ 0 \\ -1
    \end{pmatrix}, \qquad
    A =
    \begin{pmatrix}
    2 \\ 1 \\ -2
    \end{pmatrix}.

We can change the basis to go to the natural system :math:`(m, c, \hbar)` with
dimension (L: length, V: velocity, A: action) [#]_. In this basis the
generators are

.. math::

    A =
    \begin{pmatrix}
    1 \\ 0 \\ 0
    \end{pmatrix}, \qquad
    L =
    \begin{pmatrix}
    0 \\ 1 \\ 0
    \end{pmatrix}, \qquad
    V =
    \begin{pmatrix}
    0 \\ 0 \\ 1
    \end{pmatrix},

whereas the mass and time are given by

.. math::

    T = L V^{-1}, \qquad
    M = A V^{-2},\\
    T =
    \begin{pmatrix}
    0 \\ 1 \\ -1
    \end{pmatrix}, \qquad
    M =
    \begin{pmatrix}
    1 \\ 0 \\ -2
    \end{pmatrix}.

Finally the inverse change of basis matrix :math:`P^{-1}` is obtained by
gluing the vectors expressed in the old basis:

.. math::

    P^{-1} =
    \begin{pmatrix}
    2 & 1 & 1 \\
    1 & 0 & 0 \\
    -2 & 0 & -1
    \end{pmatrix}.

To find the change of basis matrix we just have to take the inverse

.. math::

    P =
    \begin{pmatrix}
    0 & 1 & 0 \\
    1 & 0 & 1 \\
    0 & -2 & -1
    \end{pmatrix}.


Quantities
==========

A quantity is defined by its name, dimension and factor to a canonical quantity
of the same dimension. The canonical quantities are an internal reference of
the units module and should not be relevant for end-users. Both units and
physical constants are quantities.

Units
-----

Units, such as meters,
seconds and kilograms, are usually reference quantities chosen by men to refer
to other quantities.

After defining several units of different dimensions we can form a unit system,
which is basically a dimension system with a notion of scale.

Constants
---------

Physical constants are just quantities. They indicate that we used not
to understand that two dimensions are in fact the same. For example, we see
a velocity for the light different from 1 because we do not think that time
is the same as space (which is normal because of our sense; but it is different
at the fundamental level). For example, once there was the "heat constant"
which allowed to convert between joules and calories since people did not know
that heat was energy. As soon as they understood it they fixed this constant to
1 (this is a very schematic story).

We can interpret the fact that now we fix the value of fundamental constants
in the SI as showing that they are units (and we use them to define the other
usual units).


The need for a reference
========================

It is not possible to define from scratch units and unit systems: one needs
to define some references, and then build the rest over them. Said in another
way, we need an origin for the scales of our units (i.e. a unit with factor 1),
and to be sure that all units of a given dimension are defined consistently we
need to use the same origin for all of them. This can happen if we want to use
a derived unit as a base units in another system: we should not define it as
having a scale 1, because, even if it is inconsistent inside the system, we
could not convert to the first system since we have two different units (from
our point of view) of same scale (which means they are equal for the computer).

We will say that the dimensions and scales defined outside systems are
canonical, because we use them for all computations. On the other side the
dimensions and scales obtained with reference to a system are called physical,
because they ultimately carry a sense.

Let's use a concrete (and important) example: the case of the mass units.
We would like to define the gram as the origin. We would like to define the
gram as the canonical origin for the mass, so we assign it a scale 1. Then we
can define a system (e.g. in chemistry) that take it as a base unit. The
MKS system prefers to use the kilogram; a naive choice would be to attribute it
a scale if 1 since it is a base, but we see that we could not convert to the
chemistry system because g and kg have both been given the same factor. So we
need to define kg as 1000 g, and only then use it as a base in MKS. But as soon
as we ask the question "what is the factor of kg in MKS?", we get the answer 1,
since it is a base unit.

Thus we will define all computations without referring to a system, and it is
only at the end that we can plug the result into a system to give the context
we are interested in.


Literature
==========

.. [Page52] C. H. Page, `Classes of units in the SI
    <https://doi.org/10.1119/1.1927482>`_,
    Am. J. of Phys. 20, 1 (1952): 1.

.. [Page78] C. H. Page, `Units and Dimensions in Physics
    <https://pubs.aip.org/aapt/ajp/article-abstract/20/1/1/1034555/Units-and-Dimensions-in-Physics>`_,
    Am. J. of Phys. 46, 1 (1978): 78.

.. [deBoer79] J. de Boer, `Group properties of quantities and units
    <https://aapt.scitation.org/doi/10.1119/1.11703>`_,
    Am. J. of Phys. 47, 9 (1979): 818.

.. [LevyLeblond77] J.-M. Lévy-Leblond, `On the Conceptual Nature of the
    Physical Constants
    <https://link.springer.com/article/10.1007/BF02748049>`_,
    La Rivista Del Nuovo Cimento 7, no. 2 (1977): 187-214.

.. [NIST] `NIST reference on constants, units and uncertainties
    <https://physics.nist.gov/cuu/Units/introduction.html>`_.


.. rubric:: Footnotes

.. [#] In general we will consider only dimensions with a maximum coefficient,
    so we can only a truncation of the group; but this is not useful for the
    algorithm.
.. [#] We anticipate a little by considering :math:`c` and :math:`\hbar` as
    units and not as physical constants.
