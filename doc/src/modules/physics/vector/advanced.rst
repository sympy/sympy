============================================================================
Potential Issues/Advanced Topics/Future Features in Physics/Vector Module
============================================================================

This document will describe some of the more advanced functionality that this
module offers but which is not part of the "official" interface. Here, some of
the features that will be implemented in the future will also be covered, along
with unanswered questions about proper functionality. Also, common problems
will be discussed, along with some solutions.

.. _Dyadic:

Dyadic
======

In :mod:`sympy.physics.mechanics`, dyadics are used to represent inertia ([Kane1985]_,
[WikiDyadics]_, [WikiDyadicProducts]_). A dyadic is a linear polynomial of
component unit dyadics, similar to a vector being a linear polynomial of
component unit vectors. A dyadic is the outer product between two vectors which
returns a new quantity representing the juxtaposition of these two vectors. For
example:

.. math::
  \mathbf{\hat{a}_x} \otimes \mathbf{\hat{a}_x} &= \mathbf{\hat{a}_x}
  \mathbf{\hat{a}_x}\\
  \mathbf{\hat{a}_x} \otimes \mathbf{\hat{a}_y} &= \mathbf{\hat{a}_x}
  \mathbf{\hat{a}_y}\\

Where :math:`\mathbf{\hat{a}_x}\mathbf{\hat{a}_x}` and
`\mathbf{\hat{a}_x}\mathbf{\hat{a}_y}` are the outer products obtained by
multiplying the left side as a column vector by the right side as a row vector.
Note that the order is significant.

Some additional properties of a dyadic are:

.. math::
  (x \mathbf{v}) \otimes \mathbf{w} &= \mathbf{v} \otimes (x \mathbf{w}) = x
  (\mathbf{v} \otimes \mathbf{w})\\
  \mathbf{v} \otimes (\mathbf{w} + \mathbf{u}) &= \mathbf{v} \otimes \mathbf{w}
  + \mathbf{v} \otimes \mathbf{u}\\
  (\mathbf{v} + \mathbf{w}) \otimes \mathbf{u} &= \mathbf{v} \otimes \mathbf{u}
  + \mathbf{w} \otimes \mathbf{u}\\

A vector in a reference frame can be represented as
:math:`\begin{bmatrix}a\\b\\c\end{bmatrix}` or :math:`a \mathbf{\hat{i}} + b
\mathbf{\hat{j}} + c \mathbf{\hat{k}}`. Similarly, a dyadic can be represented
in tensor form:

.. math::
  \begin{bmatrix}
  a_{11} & a_{12} & a_{13} \\
  a_{21} & a_{22} & a_{23} \\
  a_{31} & a_{32} & a_{33}
  \end{bmatrix}\\

or in dyadic form:

.. math::
  a_{11} \mathbf{\hat{a}_x}\mathbf{\hat{a}_x} +
  a_{12} \mathbf{\hat{a}_x}\mathbf{\hat{a}_y} +
  a_{13} \mathbf{\hat{a}_x}\mathbf{\hat{a}_z} +
  a_{21} \mathbf{\hat{a}_y}\mathbf{\hat{a}_x} +
  a_{22} \mathbf{\hat{a}_y}\mathbf{\hat{a}_y} +
  a_{23} \mathbf{\hat{a}_y}\mathbf{\hat{a}_z} +
  a_{31} \mathbf{\hat{a}_z}\mathbf{\hat{a}_x} +
  a_{32} \mathbf{\hat{a}_z}\mathbf{\hat{a}_y} +
  a_{33} \mathbf{\hat{a}_z}\mathbf{\hat{a}_z}\\

Just as with vectors, the later representation makes it possible to keep track
of which frames the dyadic is defined with respect to. Also, the two
components of each term in the dyadic need not be in the same frame. The
following is valid:

.. math::
  \mathbf{\hat{a}_x} \otimes \mathbf{\hat{b}_y} = \mathbf{\hat{a}_x}
  \mathbf{\hat{b}_y}

Dyadics can also be crossed and dotted with vectors; again, order matters:

.. math::
  \mathbf{\hat{a}_x}\mathbf{\hat{a}_x} \cdot \mathbf{\hat{a}_x} &=
  \mathbf{\hat{a}_x}\\
  \mathbf{\hat{a}_y}\mathbf{\hat{a}_x} \cdot \mathbf{\hat{a}_x} &=
  \mathbf{\hat{a}_y}\\
  \mathbf{\hat{a}_x}\mathbf{\hat{a}_y} \cdot \mathbf{\hat{a}_x} &= 0\\
  \mathbf{\hat{a}_x} \cdot \mathbf{\hat{a}_x}\mathbf{\hat{a}_x} &=
  \mathbf{\hat{a}_x}\\
  \mathbf{\hat{a}_x} \cdot \mathbf{\hat{a}_x}\mathbf{\hat{a}_y} &=
  \mathbf{\hat{a}_y}\\
  \mathbf{\hat{a}_x} \cdot \mathbf{\hat{a}_y}\mathbf{\hat{a}_x} &= 0\\
  \mathbf{\hat{a}_x} \times \mathbf{\hat{a}_y}\mathbf{\hat{a}_x} &=
  \mathbf{\hat{a}_z}\mathbf{\hat{a}_x}\\
  \mathbf{\hat{a}_x} \times \mathbf{\hat{a}_x}\mathbf{\hat{a}_x} &= 0\\
  \mathbf{\hat{a}_y}\mathbf{\hat{a}_x} \times \mathbf{\hat{a}_z} &=
  - \mathbf{\hat{a}_y}\mathbf{\hat{a}_y}\\

One can also take the time derivative of dyadics or express them in different
frames, just like with vectors.

Common Issues
=============
Here issues with numerically integrating code, choice of `dynamicsymbols` for
coordinate and speed representation, printing, differentiating, and
substitution will occur.

Printing
--------
The default printing options are to use sorting for ``Vector`` and ``Dyadic``
measure numbers, and have unsorted output from the ``vprint``, ``vpprint``, and
``vlatex`` functions. If you are printing something large, please use one of
those functions, as the sorting can increase printing time from seconds to
minutes.

Substitution
------------
Substitution into large expressions can be slow, and take a few minutes.

Acceleration of Points
----------------------
At a minimum, points need to have their velocities defined, as the acceleration
can be calculated by taking the time derivative of the velocity in the same
frame. If the 1 point or 2 point theorems were used to compute the velocity,
the time derivative of the velocity expression will most likely be more complex
than if you were to use the acceleration level 1 point and 2 point theorems.
Using the acceleration level methods can result in shorted expressions at this
point, which will result in shorter expressions later (such as when forming
Kane's equations).


Advanced Interfaces
===================

Here we will cover advanced options in: ``ReferenceFrame``, ``dynamicsymbols``,
and some associated functionality.

ReferenceFrame
--------------
``ReferenceFrame`` is shown as having a ``.name`` attribute and ``.x``, ``.y``,
and ``.z`` attributes for accessing the basis vectors, as well as a fairly
rigidly defined print output. If you wish to have a different set of indices
defined, there is an option for this. This will also require a different
interface for accessing the basis vectors. ::

  >>> from sympy.physics.vector import ReferenceFrame, vprint, vpprint, vlatex
  >>> N = ReferenceFrame('N', indices=['i', 'j', 'k'])
  >>> N['i']
  N['i']
  >>> N.x
  N['i']
  >>> vlatex(N.x)
  '\\mathbf{\\hat{n}_{i}}'

Also, the latex output can have custom strings; rather than just indices
though, the entirety of each basis vector can be specified. The custom latex
strings can occur without custom indices, and also overwrites the latex string
that would be used if there were custom indices. ::

  >>> from sympy.physics.vector import ReferenceFrame, vlatex
  >>> N = ReferenceFrame('N', latexs=['n1','\\mathbf{n}_2','cat'])
  >>> vlatex(N.x)
  'n1'
  >>> vlatex(N.y)
  '\\mathbf{n}_2'
  >>> vlatex(N.z)
  'cat'



dynamicsymbols
--------------
The ``dynamicsymbols`` function also has 'hidden' functionality; the variable
which is associated with time can be changed, as well as the notation for
printing derivatives. ::

  >>> from sympy import symbols
  >>> from sympy.physics.vector import dynamicsymbols, vprint
  >>> q1 = dynamicsymbols('q1')
  >>> q1
  q1(t)
  >>> dynamicsymbols._t = symbols('T')
  >>> q2 = dynamicsymbols('q2')
  >>> q2
  q2(T)
  >>> q1
  q1(t)
  >>> q1d = dynamicsymbols('q1', 1)
  >>> vprint(q1d)
  q1'
  >>> dynamicsymbols._str = 'd'
  >>> vprint(q1d)
  q1d
  >>> dynamicsymbols._str = '\''
  >>> dynamicsymbols._t = symbols('t')


Note that only dynamic symbols created after the change are different. The same
is not true for the `._str` attribute; this affects the printing output only,
so dynamic symbols created before or after will print the same way.

Also note that ``Vector``'s ``.dt`` method uses the ``._t`` attribute of
``dynamicsymbols``, along with a number of other important functions and
methods. Don't mix and match symbols representing time.



Solving Vectors
---------------
We would be looking at an example to understand it better.

PROBLEM:

A ball is thrown at an angle of 45 degrees to the horizontal with a velocity of
10 m/s in a ground frame G. Another reference frame B is at an angle of 30 degrees
with respect to G. The position of the ball in the B reference frame at time 
t = 2 seconds is at point P with a position vector of 3B.x + 4B.y.
What is the velocity vector of the ball at time t = 2 seconds w.r.t ground frame?

SOLUTION:

We can define two reference frames:
The ground reference frame, which is fixed with respect to the Earth and has its
x-axis pointing towards the east and its y-axis pointing towards the north.
The ``B`` reference frame, which is at an angle of 30 degrees with respect to ``G``. 

We can define two points in the problem:

1) Point O, which is the origin of both reference frames.
2) Point P, where the ball is located at time t = 2 seconds w.r.t ``B``.

We can define the position vector of point P in ``B`` as:
>>> r_b = 3*B.x + 4*B.y

We can also define the velocity of the ball in the ``B`` reference frame at time t = 2 seconds as:
>>> v_b = 10*cos(pi/4)*B.x + 10*sin(pi/4)*B.y - 9.82*B.z

Therefore, the velocity vector of the ball in the ``G`` reference frame at time t = 2 seconds is:
>>> v_a = v_b + omega.cross(r_b)
where ``omega`` is the angular velocity of the ``B`` reference frame w.r.t ground frame, and 
cross represents the vector cross product.

>>> from sympy import *
>>> from sympy.physics.vector import *
>>> G = ReferenceFrame('G')
>>> B = ReferenceFrame('B')
>>> B.orient_axis(G, pi/6, G.z)
>>> r_b = 3*B.x + 4*B.y
>>> v_b = 10*cos(pi/4)*B.x + 10*sin(pi/4)*B.y - 9.8*2*B.z
>>> omega = 2*B.z
>>> v_g = v_b + omega.cross(r_b)
>>> v_gx, v_gy, v_gz = symbols('v_gx, v_gy, v_gz')
>>> eqn = Eq(v_g.to_matrix(G), Matrix([[v_gx], [v_gy], [v_gz]]))
>>> eqn
Eq(Matrix([
[-5*sqrt(2)/2 - 3 + sqrt(3)*(-8 + 5*sqrt(2))/2],
[ -4 + 5*sqrt(2)/2 + sqrt(3)*(6 + 5*sqrt(2))/2],
[                                        -19.6]]), Matrix([
[v_gx],
[v_gy],
[v_gz]]))
>>> solve(eqn, [v_gx, v_gy, v_gz])
{ v_gx: -7.34001277925030,
  v_gy: 10.8554106855973,
  v_gz: -19.6000000000000
}