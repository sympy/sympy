=====================================
Scalar and Vector Field Functionality
=====================================

Introduction
============

Vectors and Scalars
-------------------

In physics, we deal with two kinds of quantities – scalars and vectors.

A scalar is an entity which only has a magnitude – no direction. Examples of 
scalar quantities include mass, electric charge, temperature, distance etc.

A vector, on the other hand, is an entity that is characterized by a 
magnitude and a direction. Examples of vector quantities are displacement, 
velocity, magnetic field, etc.

A scalar can be depicted just by a numeric, for eg. a temperature of 300 K.
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

  >>> electric_potential = 2*R[0]**2*R[1]
  >>> electric_potential
  2*R_x**2*R_y

In string representation, :mod:`R_x` denotes the :math:`X` base 
scalar assigned to :mod:`ReferenceFrame` :mod:`R`.

Scalar fields can be treated just as any other :mod:`SymPy` expression, 
for any math/calculus functionality. Hence, to differentiate the above 
electric potential with respect to :math:`x` (i.e. :mod:`R[0]`), you would 
have to use the :mod:`diff` method.

  >>> from sympy import diff
  >>> diff(electric_potential, R[0])
  4*R_x*R_y

Like vectors (and vector fields), scalar fields can also be re-expressed in 
other frames of reference, apart from the one they were defined in – assuming 
that an orientation relationship exists between the concerned frames. This 
can be done using the :mod:`express` method, in a way similar to vectors - 
but with the :mod:`variables` parameter set to :mod:`True`.

  >>> q = dynamicsymbols('q')
  >>> R1 = R.orientnew('R1', rot_type = 'Axis', amounts = [q, R.z])
  >>> express(electric_potential, R1, variables=True)
  2*(R1_x*sin(q(t)) + R1_y*cos(q(t)))*(R1_x*cos(q(t)) - R1_y*sin(q(t)))**2

Moreover, considering scalars can also be functions of time just as vectors, 
differentiation with respect to time is also possible. Depending on the 
:mod:`Symbol` s present in the expression and the frame with respect to which 
the time differentiation is being done, the output will change/remain the same.

  >>> time_derivative(electric_potential, R)
  0
  >>> time_derivative(electric_potential, R1)
  (R1_x*cos(q(t)) - R1_y*sin(q(t)))*(3*R1_x**2*cos(2*q(t)) - R1_x**2 - 
  6*R1_x*R1_y*sin(2*q(t)) - 3*R1_y**2*cos(2*q(t)) - R1_y**2)*Derivative(q(t), t)

