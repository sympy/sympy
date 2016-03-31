============================================================================
Potential Issues/Advanced Topics/Future Features in Physics/Vector Module
============================================================================

This document will describe some of the more advanced functionality that this
module offers but which is not part of the "official" interface. Here, some of
the features that will be implemented in the future will also be covered, along
with unanswered questions about proper functionality. Also, common problems
will be discussed, along with some solutions.

Inertia (Dyadics)
=================

A dyadic tensor is a second order tensor formed by the juxtaposition of a
pair of vectors. There are various operations defined with respect to dyadics,
which have been implemented in :mod:`vector` in the form of class
:mod:`Dyadic`. To know more, refer to the :mod:`Dyadic` and :mod:`Vector`
class APIs.
Dyadics are used to define the inertia of bodies within :mod:`mechanics`.
Inertia dyadics can be defined explicitly but the ``inertia`` function is
typically much more convenient for the user::

  >>> from sympy.physics.mechanics import ReferenceFrame, inertia
  >>> N = ReferenceFrame('N')

  Supply a reference frame and the moments of inertia if the object
  is symmetrical:

  >>> inertia(N, 1, 2, 3)
  (N.x|N.x) + 2*(N.y|N.y) + 3*(N.z|N.z)

  Supply a reference frame along with the products and moments of inertia
  for a general object:

  >>> inertia(N, 1, 2, 3, 4, 5, 6)
  (N.x|N.x) + 4*(N.x|N.y) + 6*(N.x|N.z) + 4*(N.y|N.x) + 2*(N.y|N.y) + 5*(N.y|N.z) + 6*(N.z|N.x) + 5*(N.z|N.y) + 3*(N.z|N.z)

Notice that the ``inertia`` function returns a dyadic with each component
represented as two unit vectors separated by a ``|``. Refer to the
:ref:`Dyadic` section for more information about dyadics.

Inertia is often expressed in a matrix, or tensor, form, especially for
numerical purposes. Since the matrix form does not contain any information
about the reference frame(s) the inertia dyadic is defined in, you must provide
one or two reference frames to extract the measure numbers from the dyadic.
There is a convenience function to do this::

  >>> inertia(N, 1, 2, 3, 4, 5, 6).to_matrix(N)
  Matrix([
  [1, 4, 6],
  [4, 2, 5],
  [6, 5, 3]])

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
  >>> N = ReferenceFrame('N', latexs=['n1','\mathbf{n}_2','cat'])
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
