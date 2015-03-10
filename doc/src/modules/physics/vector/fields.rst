=====================================
Scalar and Vector Field Functionality
=====================================

Introduction
============

Vectors and Scalars
-------------------

In physics, we deal with two kinds of quantities – scalars and vectors.

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

Fields
------

In general, a :math:`field` is a vector or scalar quantity that can be
specified everywhere in space as a function of position (Note that in general
a field may also be dependent on time and other custom variables). In this
module, we deal with 3-dimensional spaces only. Hence, a field is defined as
a function of the :math:`x`, :math:`y` and :math:`z` coordinates corresponding
to a location in 3D space.

For example, temperate in 3 dimensional space (a temperature field) can be
written as :math:`T(x, y, z)` – a scalar function of the position.
An example of a scalar field in electromagnetism is the electric potential.

In a similar manner, a vector field can be defined as a vectorial function
of the location :math:`(x, y, z)` of any point in space.

For instance, every point on the earth may be considered to be in the
gravitational force field of the earth. We may specify the field by the
magnitude and the direction of acceleration due to gravity
(i.e. force per unit mass ) :math:`g(x, y, z)` at every point in space.

To give an example from electromagnetism, consider an electric potential
of form :math:`2{x}^{2}y`, a scalar field in 3D space. The corresponding
conservative electric field can be computed as the gradient of the electric
potential function, and expressed as :math:`4xy\mathbf{\hat{i}} +
2{x}^{2}\mathbf{\hat{j}}`.
The magnitude of this electric field can in turn be expressed
as a scalar field of the form
:math:`\sqrt{4{x}^{4} + 16{x}^{2}{y}^{2}}`.

Implementation of fields in sympy.physics.vector
================================================

In sympy.physics.vector, every :mod:`ReferenceFrame` instance is assigned basis
vectors corresponding to the :math:`X`, :math:`Y` and
:math:`Z` directions. These can be accessed using the attributes
named :mod:`x`, :mod:`y` and :mod:`z` respectively. Hence, to define a vector
:math:`\mathbf{v}` of the form
:math:`3\mathbf{\hat{i}} + 4\mathbf{\hat{j}} + 5\mathbf{\hat{k}}` with
respect to a given frame :math:`\mathbf{R}`, you would do

  >>> from sympy.physics.vector import ReferenceFrame
  >>> R = ReferenceFrame('R')
  >>> v = 3*R.x + 4*R.y + 5*R.z

Vector math and basic calculus operations with respect to vectors have
already been elaborated upon in other sections of this module's
documentation.

On the other hand, base scalars (or coordinate variables) are implemented
as special :mod:`SymPy` :mod:`Symbol` s assigned to every frame, one for each
direction from :math:`X`, :math:`Y` and :math:`Z`. For a frame
:mod:`R`, the :math:`X`, :math:`Y` and :math:`Z`
base scalar :mod:`Symbol` s can be accessed using the :mod:`R[0]`, :mod:`R[1]`
and :mod:`R[2]` expressions respectively.

Therefore, to generate the expression for the aforementioned electric
potential field :math:`2{x}^{2}y`, you would have to do

  >>> from sympy.physics.vector import ReferenceFrame
  >>> R = ReferenceFrame('R')
  >>> electric_potential = 2*R[0]**2*R[1]
  >>> electric_potential
  2*R_x**2*R_y

In string representation, :mod:`R_x` denotes the :math:`X` base
scalar assigned to :mod:`ReferenceFrame` :mod:`R`. Essentially, :mod:`R_x` is
the string representation of :mod:`R[0]`.

Scalar fields can be treated just as any other :mod:`SymPy` expression,
for any math/calculus functionality. Hence, to differentiate the above
electric potential with respect to :math:`x` (i.e. :mod:`R[0]`), you would
have to use the :mod:`diff` method.

  >>> from sympy.physics.vector import ReferenceFrame
  >>> R = ReferenceFrame('R')
  >>> electric_potential = 2*R[0]**2*R[1]
  >>> from sympy import diff
  >>> diff(electric_potential, R[0])
  4*R_x*R_y

Like vectors (and vector fields), scalar fields can also be re-expressed in
other frames of reference, apart from the one they were defined in – assuming
that an orientation relationship exists between the concerned frames. This
can be done using the :mod:`express` method, in a way similar to vectors -
but with the :mod:`variables` parameter set to :mod:`True`.

  >>> from sympy.physics.vector import ReferenceFrame
  >>> R = ReferenceFrame('R')
  >>> electric_potential = 2*R[0]**2*R[1]
  >>> from sympy.physics.vector import dynamicsymbols, express
  >>> q = dynamicsymbols('q')
  >>> R1 = R.orientnew('R1', rot_type = 'Axis', amounts = [q, R.z])
  >>> express(electric_potential, R1, variables=True)
  2*(R1_x*sin(q(t)) + R1_y*cos(q(t)))*(R1_x*cos(q(t)) - R1_y*sin(q(t)))**2

Moreover, considering scalars can also be functions of time just as vectors,
differentiation with respect to time is also possible. Depending on the
:mod:`Symbol` s present in the expression and the frame with respect to which
the time differentiation is being done, the output will change/remain the same.

  >>> from sympy.physics.vector import ReferenceFrame
  >>> R = ReferenceFrame('R')
  >>> electric_potential = 2*R[0]**2*R[1]
  >>> q = dynamicsymbols('q')
  >>> R1 = R.orientnew('R1', rot_type = 'Axis', amounts = [q, R.z])
  >>> from sympy.physics.vector import time_derivative
  >>> time_derivative(electric_potential, R)
  0
  >>> time_derivative(electric_potential, R1).simplify()
  (R1_x*cos(q(t)) - R1_y*sin(q(t)))*(3*R1_x**2*cos(2*q(t)) - R1_x**2 -
  6*R1_x*R1_y*sin(2*q(t)) - 3*R1_y**2*cos(2*q(t)) - R1_y**2)*Derivative(q(t), t)

Field operators and other related functions
===========================================

Here we describe some basic field-related functionality implemented in
sympy.physics.vector

Curl
----

A curl is a mathematical operator that describes an infinitesimal rotation of a
vector in 3D space. The direction is determined by the right-hand rule (along the
axis of rotation), and the magnitude is given by the magnitude of rotation.

In the 3D Cartesian system, the curl of a 3D vector :math:`\mathbf{F}` ,
denoted by :math:`\nabla \times \mathbf{F}` is given by -

:math:`\nabla \times \mathbf{F} = \left(\frac{\partial F_z}{\partial y}  -
\frac{\partial F_y}{\partial z}\right) \mathbf{\hat{i}} +
\left(\frac{\partial F_x}{\partial z} -
\frac{\partial F_z}{\partial x}\right) \mathbf{\hat{j}} +
\left(\frac{\partial F_y}{\partial x} -
\frac{\partial F_x}{\partial y}\right) \mathbf{\hat{k}}`

where :math:`F_x` denotes the :math:`X` component of vector :math:`\mathbf{F}`.

To compute the curl of a vector field in :mod:`physics.vector`, you would do

  >>> from sympy.physics.vector import ReferenceFrame
  >>> R = ReferenceFrame('R')
  >>> from sympy.physics.vector import curl
  >>> field = R[0]*R[1]*R[2]*R.x
  >>> curl(field, R)
  R_x*R_y*R.y - R_x*R_z*R.z

Divergence
----------

Divergence is a vector operator that measures the magnitude of a vector field's
source or sink at a given point, in terms of a signed scalar.

The divergence operator always returns a scalar after operating on a vector.

In the 3D Cartesian system, the divergence of a 3D vector :math:`\mathbf{F}`,
denoted by :math:`\nabla\cdot\mathbf{F}` is given by -

:math:`\nabla\cdot\mathbf{F} =\frac{\partial U}{\partial x}
+\frac{\partial V}{\partial y}
+\frac{\partial W}{\partial z
}`

where :math:`U`, :math:`V` and :math:`W` denote the :math:`X`, :math:`Y` and
:math:`Z` components of :math:`\mathbf{F}` respectively.

To compute the divergence of a vector field in :mod:`physics.vector`, you
would do

  >>> from sympy.physics.vector import ReferenceFrame
  >>> R = ReferenceFrame('R')
  >>> from sympy.physics.vector import divergence
  >>> field = R[0]*R[1]*R[2] * (R.x+R.y+R.z)
  >>> divergence(field, R)
  R_x*R_y + R_x*R_z + R_y*R_z

Gradient
--------

Consider a scalar field :math:`f(x, y, z)` in 3D space. The gradient of this field
is defined as the vector of the 3 partial derivatives of :math:`f` with respect to
:math:`x`, :math:`y` and :math:`z` in the :math:`X`, :math:`Y` and :math:`Z`
directions respectively.

In the 3D Cartesian system, the divergence of a scalar field :math:`f`,
denoted by :math:`\nabla f` is given by -

:math:`\nabla f = \frac{\partial f}{\partial x} \mathbf{\hat{i}} +
\frac{\partial f}{\partial y}  \mathbf{\hat{j}} +
\frac{\partial f}{\partial z} \mathbf{\hat{k}}`

To compute the divergence of a vector field in :mod:`physics.vector`, you
would do

  >>> from sympy.physics.vector import ReferenceFrame
  >>> R = ReferenceFrame('R')
  >>> from sympy.physics.vector import gradient
  >>> scalar_field = R[0]*R[1]*R[2]
  >>> gradient(scalar_field, R)
  R_y*R_z*R.x + R_x*R_z*R.y + R_x*R_y*R.z

Conservative and Solenoidal fields
----------------------------------

In vector calculus, a conservative field is a field that is the gradient of
some scalar field. Conservative fields have the property that their line
integral over any path depends only on the end-points, and is independent
of the path between them.
A conservative vector field is also said to be 'irrotational', since the
curl of a conservative field is always zero.

In physics, conservative fields represent forces in physical systems where
energy is conserved.

To check if a vector field is conservative in :mod:`physics.vector`, use
the :mod:`is_conservative` function.

  >>> from sympy.physics.vector import ReferenceFrame, is_conservative
  >>> R = ReferenceFrame('R')
  >>> field = R[1]*R[2]*R.x + R[0]*R[2]*R.y + R[0]*R[1]*R.z
  >>> is_conservative(field)
  True
  >>> curl(field, R)
  0

A solenoidal field, on the other hand, is a vector field whose divergence
is zero at all points in space.

To check if a vector field is solenoidal in :mod:`physics.vector`, use
the :mod:`is_solenoidal` function.

  >>> from sympy.physics.vector import ReferenceFrame, is_solenoidal
  >>> R = ReferenceFrame('R')
  >>> field = R[1]*R[2]*R.x + R[0]*R[2]*R.y + R[0]*R[1]*R.z
  >>> is_solenoidal(field)
  True
  >>> divergence(field, R)
  0

Scalar potential functions
--------------------------

We have previously mentioned that every conservative field can be defined as
the gradient of some scalar field. This scalar field is also called the 'scalar
potential field' corresponding to the aforementioned conservative field.

The :mod:`scalar_potential` function in :mod:`physics.vector` calculates the
scalar potential field corresponding to a given conservative vector field in
3D space - minus the extra constant of integration, of course.

Example of usage -

  >>> from sympy.physics.vector import ReferenceFrame, scalar_potential
  >>> R = ReferenceFrame('R')
  >>> conservative_field = 4*R[0]*R[1]*R[2]*R.x + 2*R[0]**2*R[2]*R.y + 2*R[0]**2*R[1]*R.z
  >>> scalar_potential(conservative_field, R)
  2*R_x**2*R_y*R_z

Providing a non-conservative vector field as an argument to
:mod:`scalar_potential` raises a :mod:`ValueError`.

The scalar potential difference, or simply 'potential difference',
corresponding to a conservative vector field can be defined as the difference
between the values of its scalar potential function at two points in space.
This is useful in calculating a line integral with respect to a conservative
function, since it depends only on the endpoints of the path.

This computation is performed as follows in :mod:`physics.vector`.

  >>> from sympy.physics.vector import ReferenceFrame, Point
  >>> from sympy.physics.vector import scalar_potential_difference
  >>> R = ReferenceFrame('R')
  >>> O = Point('O')
  >>> P = O.locatenew('P', 1*R.x + 2*R.y + 3*R.z)
  >>> vectfield = 4*R[0]*R[1]*R.x + 2*R[0]**2*R.y
  >>> scalar_potential_difference(vectfield, R, O, P, O)
  4

If provided with a scalar expression instead of a vector field,
:mod:`scalar_potential_difference` returns the difference between the values
of that scalar field at the two given points in space.
