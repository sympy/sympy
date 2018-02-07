============
Introduction
============

This page gives a brief conceptual overview of the functionality present in
:mod:`sympy.vector`.

Vectors and Scalars
===================

In vector math, we deal with two kinds of quantities – scalars and vectors.

A **scalar** is an entity which only has a magnitude – no direction. Examples of
scalar quantities include mass, electric charge, temperature, distance, etc.

A **vector**, on the other hand, is an entity that is characterized by a
magnitude and a direction. Examples of vector quantities are displacement,
velocity, magnetic field, etc.

A scalar can be depicted just by a number, for e.g. a temperature of 300 K.
On the other hand, vectorial quantities like acceleration are usually denoted
by a vector. Given a vector :math:`\mathbf{V}`, the magnitude of the
corresponding quantity can be calculated as the magnitude of the vector
itself :math:`\Vert \mathbf{V} \Vert`, while the direction would be specified
by a unit vector in the direction of the original vector,
:math:`\mathbf{\hat{V}} = \frac{\mathbf{V}}{\Vert \mathbf{V} \Vert}`.

For example, consider a displacement of
:math:`(3\mathbf{\hat{i}} + 4\mathbf{\hat{j}} + 5\mathbf{\hat{k}})` m,
where , as per standard convention, :math:`\mathbf{\hat{i}}`,
:math:`\mathbf{\hat{j}}` and :math:`\mathbf{\hat{k}}` represent unit vectors
along the :math:`\mathbf{X}`, :math:`\mathbf{Y}` and :math:`\mathbf{Z}`
axes respectively. Therefore, it can be concluded that the distance
traveled is
:math:`\Vert 3\mathbf{\hat{i}} + 4\mathbf{\hat{j}} + 5\mathbf{\hat{k}} \Vert`
m = :math:`5\sqrt{2}` m. The direction of travel is given by the unit vector
:math:`\frac{3}{5\sqrt{2}}\mathbf{\hat{i}} +
\frac{4}{5\sqrt{2}}\mathbf{\hat{j}} + \frac{5}{5\sqrt{2}}\mathbf{\hat{k}}`.

Coordinate Systems
==================

A **coordinate system** is an abstract mathematical entity used to define
the notion of directions and locations in n-dimensional spaces. This
module deals with 3-dimensional spaces, with the conventional :math:`X`,
:math:`Y` and :math:`Z` axes defined with respect
to each coordinate system.

Each coordinate system also has a special reference point called the
'origin' defined for it. This point is used either while referring to
locations in 3D space, or while calculating the coordinates of
pre-defined points with respect to the system.

It is a pretty well-known concept that there is no absolute notion
of location or orientation in space. Any given coordinate system
defines a unique 'perspective' of quantifying positions and directions.
Therefore, even if we assume that all systems deal with the same
units of measurement, the expression of vectorial and scalar quantities
differs according to the coordinate system a certain observer deals with.

Consider two points :math:`P` and :math:`Q` in space. Assuming units to
be common throughtout, the distance between these points remains
the same regardless of the coordinate system in which the measurements are
being made. However, the 3-D coordinates of each of the two points, as well
as the position vector of any of the points with respect to the other,
do not.
In fact, these two quantities don't make sense at all, unless they are being
measured keeping in mind a certain location and orientation of the measurer
(essentially the coordinate system).

Therefore, it is quite clear that the orientation and location (of the origin)
of a coordinate system define the way different quantities will be expressed
with respect to it.  Neither of the two properties can be measured on an
absolute scale, but rather with respect to another coordinate system. The
orientation of one system with respect to another is measured using the
the rotation matrix, while the relative position can be quantified via
the position vector of one system's origin with respect to the other.

Fields
======

A **field** is a vector or scalar quantity that can be
specified everywhere in space as a function of position (Note that in general
a field may also be dependent on time and other custom variables). Since we
only deal with 3D spaces in this module, a field is defined as a function of
the :math:`x`, :math:`y` and :math:`z` coordinates corresponding
to a location in the coordinate system. Here, :math:`x`, :math:`y` and
:math:`z` act as scalar variables defining the position of a general point.

For example, temperature in 3 dimensional space (a temperature field) can be
written as :math:`T(x, y, z)` – a scalar function of the position.
An example of a scalar field in electromagnetism is the electric potential.

In a similar manner, a vector field can be defined as a vectorial function
of the location :math:`(x, y, z)` of any point in space.

For instance, every point on the earth may be considered to be in the
gravitational force field of the earth. We may specify the field by the
magnitude and the direction of acceleration due to gravity
(i.e. force per unit mass ) :math:`\vec g(x, y, z)` at every point in
space.

To give an example from electromagnetism, consider an electric potential
of form :math:`2{x}^{2}y`, a scalar field in 3D space. The corresponding
conservative electric field can be computed as the gradient of the electric
potential function, and expressed as :math:`4xy\mathbf{\hat{i}} +
2{x}^{2}\mathbf{\hat{j}}`.
The magnitude of this electric field can in turn be expressed
as a scalar field of the form
:math:`\sqrt{4{x}^{4} + 16{x}^{2}{y}^{2}}`.
