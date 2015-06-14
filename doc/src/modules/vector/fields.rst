=====================================
Scalar and Vector Field Functionality
=====================================

Implementation in sympy.vector
==============================

Scalar and vector fields
------------------------

In :mod:`sympy.vector`, every ``CoordSysCartesian`` instance is assigned basis
vectors corresponding to the :math:`X`, :math:`Y` and
:math:`Z` axes. These can be accessed using the properties
named ``i``, ``j`` and ``k`` respectively. Hence, to define a vector
:math:`\mathbf{v}` of the form
:math:`3\mathbf{\hat{i}} + 4\mathbf{\hat{j}} + 5\mathbf{\hat{k}}` with
respect to a given frame :math:`\mathbf{R}`, you would do

  >>> from sympy.vector import CoordSysCartesian
  >>> R = CoordSysCartesian('R')
  >>> v = 3*R.i + 4*R.j + 5*R.k

Vector math and basic calculus operations with respect to vectors have
already been elaborated upon in the earlier section of this module's
documentation.

On the other hand, base scalars (or coordinate variables) are implemented
in a special class called ``BaseScalar``, and are assigned to every 
coordinate system, one for each axis from :math:`X`, :math:`Y` and 
:math:`Z`. These coordinate variables are used to form the expressions of
vector or scalar fields in 3D space.
For a system ``R``, the :math:`X`, :math:`Y` and :math:`Z` 
``BaseScalars`` instances can be accessed using the ``R.x``, ``R.y``
and ``R.z`` expressions respectively.

Therefore, to generate the expression for the aforementioned electric
potential field :math:`2{x}^{2}y`, you would have to do

  >>> from sympy.vector import CoordSysCartesian
  >>> R = CoordSysCartesian('R')
  >>> electric_potential = 2*R.x**2*R.y
  >>> electric_potential
  2*R.x**2*R.y

It is to be noted that ``BaseScalar`` instances can be used just
like any other SymPy ``Symbol``, except that they store the information
about the coordinate system and axis they correspond to.

Scalar fields can be treated just as any other SymPy expression,
for any math/calculus functionality. Hence, to differentiate the above
electric potential with respect to :math:`x` (i.e. ``R.x``), you would
use the ``diff`` method.

  >>> from sympy.vector import CoordSysCartesian
  >>> R = CoordSysCartesian('R')
  >>> electric_potential = 2*R.x**2*R.y
  >>> from sympy import diff
  >>> diff(electric_potential, R.x)
  4*R.x*R.y

It is worth noting that having a ``BaseScalar`` in the expression implies
that a 'field' changes with position, in 3D space. Technically speaking, a
simple ``Expr`` with no ``BaseScalar`` s is still a field, though 
constant.

Like scalar fields, vector fields that vary with position can also be 
constructed using ``BaseScalar`` s in the measure-number expressions.

  >>> from sympy.vector import CoordSysCartesian
  >>> R = CoordSysCartesian('R')
  >>> v = R.x**2*R.i + 2*R.x*R.z*R.k

The Del operator
----------------

The Del, or 'Nabla' operator - written as :math:`\mathbf{\nabla}` is
commonly known as the vector differential operator. Depending on its 
usage in a mathematical expression, it may denote the gradient of a
scalar field, or the divergence of a vector field, or the curl of a
vector field.

Essentially, :math:`\mathbf{\nabla}` is not technically an 'operator',
but a convenient mathematical notation to denote any one of the
aforementioned field operations.

In :mod:`sympy.vector`, :math:`\mathbf{\nabla}` has been implemented
as the ``delop`` property of the ``CoordSysCartesian`` class.
Hence, assuming ``C`` is a coordinate system, the 
:math:`\mathbf{\nabla}` operator corresponding to the vector
differentials wrt ``C``'s coordinate variables and basis vectors
would be accessible as ``C.delop``.

Given below is an example of usage of the ``delop`` object.

  >>> from sympy.vector import CoordSysCartesian
  >>> C = CoordSysCartesian('C')
  >>> gradient_field = C.delop(C.x*C.y*C.z)
  >>> gradient_field
  (Derivative(C.x*C.y*C.z, C.x))*C.i + (Derivative(C.x*C.y*C.z, C.y))*C.j + (Derivative(C.x*C.y*C.z, C.z))*C.k

The above expression can be evaluated using the SymPy ``doit()``
routine.

  >>> gradient_field.doit()
  C.y*C.z*C.i + C.x*C.z*C.j + C.x*C.y*C.k

Usage of the :math:`\mathbf{\nabla}` notation in :mod:`sympy.vector`
has been described in greater detail in the subsequent subsections.

Field operators and related functions
=====================================

Here we describe some basic field-related functionality implemented in
:mod:`sympy.vector`.

Curl
----

A curl is a mathematical operator that describes an infinitesimal rotation of a
vector in 3D space. The direction is determined by the right-hand rule (along the
axis of rotation), and the magnitude is given by the magnitude of rotation.

In the 3D Cartesian system, the curl of a 3D vector :math:`\mathbf{F}` ,
denoted by :math:`\nabla \times \mathbf{F}` is given by:

:math:`\nabla \times \mathbf{F} = \left(\frac{\partial F_z}{\partial y}  -
\frac{\partial F_y}{\partial z}\right) \mathbf{\hat{i}} +
\left(\frac{\partial F_x}{\partial z} -
\frac{\partial F_z}{\partial x}\right) \mathbf{\hat{j}} +
\left(\frac{\partial F_y}{\partial x} -
\frac{\partial F_x}{\partial y}\right) \mathbf{\hat{k}}`

where :math:`F_x` denotes the :math:`X` component of vector :math:`\mathbf{F}`.

Computing the curl of a vector field in :mod:`sympy.vector` can be 
accomplished in two ways.

One, by using the ``delop`` property

  >>> from sympy.vector import CoordSysCartesian
  >>> C = CoordSysCartesian('C')
  >>> C.delop.cross(C.x*C.y*C.z*C.i).doit()
  C.x*C.y*C.j + (-C.x*C.z)*C.k
  >>> (C.delop ^ C.x*C.y*C.z*C.i).doit()
  C.x*C.y*C.j + (-C.x*C.z)*C.k

Or by using the dedicated function

  >>> from sympy.vector import curl
  >>> curl(C.x*C.y*C.z*C.i, C)
  C.x*C.y*C.j + (-C.x*C.z)*C.k

Divergence
----------

Divergence is a vector operator that measures the magnitude of a vector field's
source or sink at a given point, in terms of a signed scalar.

The divergence operator always returns a scalar after operating on a vector.

In the 3D Cartesian system, the divergence of a 3D vector :math:`\mathbf{F}`,
denoted by :math:`\nabla\cdot\mathbf{F}` is given by:

:math:`\nabla\cdot\mathbf{F} =\frac{\partial U}{\partial x}
+\frac{\partial V}{\partial y}
+\frac{\partial W}{\partial z
}`

where :math:`U`, :math:`V` and :math:`W` denote the :math:`X`, :math:`Y` and
:math:`Z` components of :math:`\mathbf{F}` respectively.

Computing the divergence of a vector field in :mod:`sympy.vector` can be 
accomplished in two ways.

One, by using the ``delop`` property

  >>> from sympy.vector import CoordSysCartesian
  >>> C = CoordSysCartesian('C')
  >>> C.delop.dot(C.x*C.y*C.z*(C.i + C.j + C.k)).doit()
  C.x*C.y + C.x*C.z + C.y*C.z
  >>> (C.delop & C.x*C.y*C.z*(C.i + C.j + C.k)).doit()
  C.x*C.y + C.x*C.z + C.y*C.z

Or by using the dedicated function

  >>> from sympy.vector import divergence
  >>> divergence(C.x*C.y*C.z*(C.i + C.j + C.k), C)
  C.x*C.y + C.x*C.z + C.y*C.z

Gradient
--------

Consider a scalar field :math:`f(x, y, z)` in 3D space. The gradient of this field
is defined as the vector of the 3 partial derivatives of :math:`f` with respect to
:math:`x`, :math:`y` and :math:`z` in the :math:`X`, :math:`Y` and :math:`Z`
axes respectively.

In the 3D Cartesian system, the divergence of a scalar field :math:`f`,
denoted by :math:`\nabla f` is given by -

:math:`\nabla f = \frac{\partial f}{\partial x} \mathbf{\hat{i}} +
\frac{\partial f}{\partial y}  \mathbf{\hat{j}} +
\frac{\partial f}{\partial z} \mathbf{\hat{k}}`

Computing the divergence of a vector field in :mod:`sympy.vector` can be 
accomplished in two ways.

One, by using the ``delop`` property

  >>> from sympy.vector import CoordSysCartesian
  >>> C = CoordSysCartesian('C')
  >>> C.delop.gradient(C.x*C.y*C.z).doit()
  C.y*C.z*C.i + C.x*C.z*C.j + C.x*C.y*C.k
  >>> C.delop(C.x*C.y*C.z).doit()
  C.y*C.z*C.i + C.x*C.z*C.j + C.x*C.y*C.k

Or by using the dedicated function

  >>> from sympy.vector import gradient
  >>> gradient(C.x*C.y*C.z, C)
  C.y*C.z*C.i + C.x*C.z*C.j + C.x*C.y*C.k

Directional Derivative
----------------------

Apart from the above three common applications of :math:`\mathbf{\nabla}`,
it is also possible to compute the directional derivative of a field wrt
a ``Vector`` in :mod:`sympy.vector`.

By definition, the directional derivative of a field :math:`\mathbf{F}`
along a vector :math:`v` at point :math:`x` represents the instantaneous 
rate of change of :math:`\mathbf{F}` moving through :math:`x` with the
velocity :math:`v`. It is represented mathematically as:
:math:`(\vec v \cdot \nabla) \, \mathbf{F}(x)`.

Directional derivatives of vector and scalar fields can be computed in
:mod:`sympy.vector` using the ``delop`` property of
``CoordSysCartesian``.

  >>> from sympy.vector import CoordSysCartesian
  >>> C = CoordSysCartesian('C')
  >>> vel = C.i + C.j + C.k
  >>> scalar_field = C.x*C.y*C.z
  >>> vector_field = C.x*C.y*C.z*C.i
  >>> (vel.dot(C.delop))(scalar_field)
  C.x*C.y + C.x*C.z + C.y*C.z
  >>> (vel & C.delop)(vector_field)
  (C.x*C.y + C.x*C.z + C.y*C.z)*C.i

Conservative and Solenoidal fields
==================================

In vector calculus, a conservative field is a field that is the gradient of
some scalar field. Conservative fields have the property that their line
integral over any path depends only on the end-points, and is independent
of the path travelled.
A conservative vector field is also said to be 'irrotational', since the
curl of a conservative field is always zero.

In physics, conservative fields represent forces in physical systems where
energy is conserved.

To check if a vector field is conservative in :mod:`sympy.vector`, the 
``is_conservative`` function can be used.

  >>> from sympy.vector import CoordSysCartesian, is_conservative
  >>> R = CoordSysCartesian('R')
  >>> field = R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
  >>> is_conservative(field)
  True
  >>> curl(field, R)
  0

A solenoidal field, on the other hand, is a vector field whose divergence
is zero at all points in space.

To check if a vector field is solenoidal in :mod:`sympy.vector`, the 
``is_solenoidal`` function can be used.

  >>> from sympy.vector import CoordSysCartesian, is_solenoidal
  >>> R = CoordSysCartesian('R')
  >>> field = R.y*R.z*R.i + R.x*R.z*R.j + R.x*R.y*R.k
  >>> is_solenoidal(field)
  True
  >>> divergence(field, R)
  0

Scalar potential functions
==========================

We have previously mentioned that every conservative field can be defined as
the gradient of some scalar field. This scalar field is also called the 'scalar
potential field' corresponding to the aforementioned conservative field.

The ``scalar_potential`` function in :mod:`sympy.vector` calculates the
scalar potential field corresponding to a given conservative vector field in
3D space - minus the extra constant of integration, of course.

Example of usage -

  >>> from sympy.vector import CoordSysCartesian, scalar_potential
  >>> R = CoordSysCartesian('R')
  >>> conservative_field = 4*R.x*R.y*R.z*R.i + 2*R.x**2*R.z*R.j + 2*R.x**2*R.y*R.k
  >>> scalar_potential(conservative_field, R)
  2*R.x**2*R.y*R.z

Providing a non-conservative vector field as an argument to
``scalar_potential`` raises a ``ValueError``.

The scalar potential difference, or simply 'potential difference',
corresponding to a conservative vector field can be defined as the difference
between the values of its scalar potential function at two points in space.
This is useful in calculating a line integral with respect to a conservative
function, since it depends only on the endpoints of the path.

This computation is performed as follows in :mod:`sympy.vector`.

  >>> from sympy.vector import CoordSysCartesian, Point
  >>> from sympy.vector import scalar_potential_difference
  >>> R = CoordSysCartesian('R')
  >>> P = R.origin.locate_new('P', 1*R.i + 2*R.j + 3*R.k)
  >>> vectfield = 4*R.x*R.y*R.i + 2*R.x**2*R.j
  >>> scalar_potential_difference(vectfield, R, R.origin, P)
  4

If provided with a scalar expression instead of a vector field,
``scalar_potential_difference`` returns the difference between the values
of that scalar field at the two given points in space.
