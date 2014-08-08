===================================
Introduction to Vectors and Dyadics
===================================

Vector
======

A vector is a geometric object that has a magnitude (or length) and a
direction. Vectors in 3-space are often represented on paper as:

.. image:: vec_rep.*
   :height: 175
   :width: 350
   :align: center

Vector Algebra
==============

Vector algebra is the first topic to be discussed.

Two vectors are said to be equal if and only if (iff) they have the same 
magnitude and orientation.

Vector Operations
-----------------
Multiple algebraic operations can be done with vectors: addition between
vectors, scalar multiplication, and vector multiplication.

Vector addition as based on the parallelogram law.

.. image:: vec_add.*
   :height: 200
   :width: 200
   :align: center

Vector addition is also commutative:

.. math::
  \mathbf{a} + \mathbf{b} &= \mathbf{b} + \mathbf{a} \\
  (\mathbf{a} + \mathbf{b}) + \mathbf{c} &= \mathbf{a} + (\mathbf{b} +
  \mathbf{c})

Scalar multiplication is the product of a vector and a scalar; the result is a
vector with the same orientation but whose magnitude is scaled by the scalar.
Note that multiplication by -1 is equivalent to rotating the vector by 180
degrees about an arbitrary axis in the plane perpendicular to the vector.

.. image:: vec_mul.*
   :height: 150
   :width: 200
   :align: center

A unit vector is simply a vector whose magnitude is equal to 1.  Given any
vector :math:`\mathbf{v}` we can define a unit vector as:

.. math::
  \mathbf{\hat{n}_v} = \frac{\mathbf{v}}{\Vert \mathbf{v} \Vert}

Note that every vector can be written as the product of a scalar and unit
vector.

Three vector products are implemented in :mod:`vector`: the dot product, the
cross product, and the outer product.

The dot product operation maps two vectors to a scalar.  It is defined as:

.. math::
  \mathbf{a} \cdot \mathbf{b} = \Vert \mathbf{a} \Vert \Vert \mathbf{b}
  \Vert \cos(\theta)\\

where :math:`\theta` is the angle between :math:`\mathbf{a}` and
:math:`\mathbf{b}`.

The dot product of two unit vectors represent the magnitude of the common
direction; for other vectors, it is the product of the magnitude of the common
direction and the two vectors' magnitudes. The dot product of two perpendicular
is zero. The figure below shows some examples:

.. image:: vec_dot.*
   :height: 250
   :width: 450
   :align: center

The dot product is commutative:

.. math::
  \mathbf{a} \cdot \mathbf{b} = \mathbf{b} \cdot \mathbf{a}

The cross product vector multiplication operation of two vectors returns a
vector:

.. math::
  \mathbf{a} \times \mathbf{b} = \mathbf{c}

The vector :math:`\mathbf{c}` has the following properties: it's orientation is
perpendicular to both :math:`\mathbf{a}` and :math:`\mathbf{b}`, it's magnitude
is defined as :math:`\Vert \mathbf{c} \Vert = \Vert \mathbf{a} \Vert \Vert
\mathbf{b} \Vert \sin(\theta)` (where :math:`\theta` is the angle between
:math:`\mathbf{a}` and :math:`\mathbf{b}`), and has a sense defined by using
the right hand rule between :math:`\Vert \mathbf{a} \Vert \Vert \mathbf{b}
\Vert`. The figure below shows this:

.. image:: vec_cross.*
   :height: 350
   :width: 700
   :align: center

The cross product has the following properties:

It is not commutative:

.. math::
  \mathbf{a} \times \mathbf{b} &\neq \mathbf{b} \times \mathbf{a} \\
  \mathbf{a} \times \mathbf{b} &= - \mathbf{b} \times \mathbf{a}

and not associative:

.. math::
  (\mathbf{a} \times \mathbf{b} ) \times \mathbf{c} \neq \mathbf{a} \times
  (\mathbf{b} \times \mathbf{c})

Two parallel vectors will have a zero cross product.

The outer product between two vectors results in the formation of a dyadic
tensor. 

A dyadic tensor is a second order tensor formed by the juxtaposition of a
pair of vectors. There are various operations defined with respect to dyadics,
which have been implemented in :mod:`vector` in the form of class
:mod:`Dyadic`. Similar to :mod:`Vector`, :mod:`Dyadic` also has all
operations such as addition, multiplication, dotting and crossing defined
with respect to it. To know more, refer to the :mod:`Dyadic` and :mod:`Vector`
class APIs.
