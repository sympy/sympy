=============
Vector Module
=============

Introduction
============

Vectors and Scalars
-------------------

In vector math, we deal with two kinds of quantities – scalars and vectors.

A scalar is an entity which only has a magnitude – no direction. Examples of
scalar quantities include mass, electric charge, temperature, distance, etc.

A vector, on the other hand, is an entity that is characterized by a
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
in the :math:`\mathbf{X}`, :math:`\mathbf{Y}` and :math:`\mathbf{Z}`
directions respectively. Therefore, it can be concluded that the distance
traveled is
:math:`\Vert 3\mathbf{\hat{i}} + 4\mathbf{\hat{j}} + 5\mathbf{\hat{k}} \Vert`
m = :math:`5\sqrt{2}` m. The direction of travel is given by the unit vector
:math:`\frac{3}{5\sqrt{2}}\mathbf{\hat{i}} +
\frac{4}{5\sqrt{2}}\mathbf{\hat{j}} + \frac{5}{5\sqrt{2}}\mathbf{\hat{k}}`.

Coordinate Systems
------------------

A coordinate system is an abstract mathematical entity used to define
the notion of directions and locations in n-dimensional spaces. This
module deals with 3-dimensional spaces, with the conventional :math:`X`, 
:math:`Y` and :math:`Z` directions (or axes) defined with respect 
to each coordinate system.

Each coordinate system also has a special reference point called the 
'origin' defined for it. This point is used either while referring to 
locations in 3D space, or while calculating the coordinates of 
pre-defined points with respect to the system.


The :mod:`sympy.vector` module provides tools for basic vector math 
and differential calculus with respect to 3D Cartesian coordinate 
systems. This documentation provides an overview of all the 
features offered, and relevant API.

Guide to Vector
===============

.. toctree::
    :maxdepth: 2

    intro.rst
    fields.rst
    coordsys.rst

Vector API
==========

.. toctree::
    :maxdepth: 2

    api/classes.rst
    api/orienterclasses.rst
    api/vectorfunctions.rst

References for Vector
================================

.. [WikiDyadics] http://en.wikipedia.org/wiki/Dyadics
.. [WikiDyadicProducts] http://en.wikipedia.org/wiki/Dyadic_product
.. [WikiDelOperator] http://en.wikipedia.org/wiki/Del
