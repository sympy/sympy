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
a specific phenomenom. In this sense dimensions are different from pure numbers
because they carry some extra-sense and for this reason two different
dimensions can not be added. For example time or length are dimensions, but
also any other things which has some sense for us, like angle, number of
particles (moles...) or information (bits...).

From this point of view the only trully dimensionless quantity are pure
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
information one need (we do not distinguish between the element of the group
and their representation):

.. math::

    (d_i)_j = \delta_{ij}, \qquad
    D =
    \begin{pmatrix}
    a_1 \\ \vdots \\ a_n
    \end{pmatrix}.

The change of basis to :math:`d'_i` follows the usual rule of change of basis
for linear space, the matrix being given by the coefficient
:math:`P_{ij}`, which are simplfy the coefficient of the new vectors in
term of the old basis:

.. math::

    d'_i = P_{ij} d_j.

We will use this last solution in our algorithm.


An example
----------

In order to illustrate all this formalism, we end this section with a specific
example, the MKS system (m, kg, s) with dimensions (L: length, M: mass,
T: time). Their are represented as (we will always sort the vectors in
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

Literature
==========

.. [Page52] C. H. Page, `Classes of units in the SI
    <http://ajp.aapt.org/resource/1/ajpias/v46/i1/p78_s1>`_,
    Am. J. of Phys. 20, 1 (1952): 1.

.. [Page78] C. H. Page, `Units and Dimensions in Physics
    <http://ajp.aapt.org/resource/1/ajpias/v20/i1/p1_s1>`_,
    Am. J. of Phys. 46, 1 (1978): 78.

.. [deBoer79] J. de Boer, `Group properties of quantities and units
    <http://ajp.aapt.org/resource/1/ajpias/v47/i9/p818_s1>`_,
    Am. J. of Phys. 47, 9 (1979): 818.

.. [LevyLeblond77] J.-M. Lévy-Leblond, `On the Conceptual Nature of the
    Physical Constants
    <http://link.springer.com/article/10.1007%2FBF02748049>`_,
    La Rivista Del Nuovo Cimento 7, no. 2 (1977): 187-214.

.. [NIST] `NIST reference on constants, units and uncertainties
    <http://physics.nist.gov/cuu/Units/introduction.html>`_.

.. rubric:: Footnotes

.. [#] In general we will consider only dimensions with a maximum coefficient,
    so we can only a truncation of the group; but this is not useful for the
    algorithm.
.. [#] We anticipate a little by considering :math:`c` and :math:`\hbar` as
    units and not as physical constants.
