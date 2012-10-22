============================================================================
Potential Issues/Advanced Topics/Future Features in Physics/Mechanics Module
============================================================================

This document will describe some of the more advanced functionality that this
module offers but which is not part of the "official" interface. Here, some of
the features that will be implemented in the future will also be covered, along
with unanswered questions about proper functionality. Also, common problems
will be discussed, along with some solutions.

Common Issues
=============
Here issues with numerically integrating code, choice of `dynamicsymbols` for
coordinate and speed representation, printing, differentiating, and
substitution will occur.

Numerically Integrating Code
----------------------------
See Future Features: Code Output

Choice of Coordinates and Speeds
--------------------------------
The Kane object is set up with the assumption that the generalized speeds are
not the same symbol as the time derivatives of the generalized coordinates.
This isn't to say that they can't be the same, just that they have to have a
different symbol. If you did this: ::

  >> KM.coords([q1, q2, q3])
  >> KM.speeds([q1d, q2d, q3d])

Your code would not work. Currently, kinematic differential equations are
required to be provided. It is at this point that we hope the user will
discover they should not attempt the behavior shown in the code above.

This behavior might not be true for other methods of forming the equations of
motion though.

Printing
--------
The default printing options are to use sorting for ``Vector`` and ``Dyad``
measure numbers, and have unsorted output from the ``mprint``, ``mpprint``, and
``mlatex`` functions. If you are printing something large, please use one of
those functions, as the sorting can increase printing time from seconds to
minutes.

Differentiating
---------------
Differentiation of very large expressions can take some time in SymPy; it is
possible for large expressions to take minutes for the derivative to be
evaluated. This will most commonly come up in linearization.

Substitution
------------
Substitution into large expressions can be slow, and take a few minutes.

Linearization
-------------
Currently, the ``Kane`` object's ``linearize`` method doesn't support cases
where there are non-coordinate, non-speed dynamic symbols outside of the
"dynamic equations". It also does not support cases where time derivatives of
these types of dynamic symbols show up. This means if you have kinematic
differential equations which have a non-coordinate, non-speed dynamic symbol,
it will not work. It also means if you have defined a system parameter (say a
length or distance or mass) as a dynamic symbol, its time derivative is likely
to show up in the dynamic equations, and this will prevent linearization.

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

  >>> from sympy.physics.mechanics import ReferenceFrame, mprint, mpprint, mlatex
  >>> N = ReferenceFrame('N', indices=['i', 'j', 'k'])
  >>> N['i']
  N['i']
  >>> N.x
  N['i']
  >>> mlatex(N.x)
  '\\mathbf{\\hat{n}_{i}}'

Also, the latex output can have custom strings; rather than just indices
though, the entirety of each basis vector can be specified. The custom latex
strings can occur without custom indices, and also overwrites the latex string
that would be used if there were custom indices. ::

  >>> from sympy.physics.mechanics import ReferenceFrame, mlatex
  >>> N = ReferenceFrame('N', latexs=['n1','\mathbf{n}_2','cat'])
  >>> mlatex(N.x)
  'n1'
  >>> mlatex(N.y)
  '\\mathbf{n}_2'
  >>> mlatex(N.z)
  'cat'

dynamicsymbols
--------------
The ``dynamicsymbols`` function also has 'hidden' functionality; the variable
which is associated with time can be changed, as well as the notation for
printing derivatives. ::

  >>> from sympy import symbols
  >>> from sympy.physics.mechanics import dynamicsymbols, mprint
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
  >>> mprint(q1d)
  q1'
  >>> dynamicsymbols._str = 'd'
  >>> mprint(q1d)
  q1d
  >>> dynamicsymbols._str = '\''
  >>> dynamicsymbols._t = symbols('t')


Note that only dynamic symbols created after the change are different. The same
is not true for the `._str` attribute; this affects the printing output only,
so dynamic symbols created before or after will print the same way.

Also note that ``Vector``'s ``.dt`` method uses the ``._t`` attribute of
``dynamicsymbols``, along with a number of other important functions and
methods. Don't mix and match symbols representing time.

Advanced Functionality
----------------------
Remember that the ``Kane`` object supports bodies which have time-varying
masses and inertias, although this functionality isn't completely compatible
with the linearization method.

Operators were discussed earlier as a potential way to do mathematical
operations on ``Vector`` and ``Dyad`` objects. The majority of the code in this
module is actually coded with them, as it can (subjectively) result in cleaner,
shorter, more readable code. If using this interface in your code, remember to
take care and use parentheses; the default order of operations in Python
results in addition occurring before some of the vector products, so use
parentheses liberally.


Future Features
===============

This will cover the planned features to be added to this submodule.

Code Output
-----------
A function for generating code output for numerical integration is the highest
priority feature to implement next. There are a number of considerations here.

Code output for C (using the GSL libraries), Fortran 90 (using LSODA), MATLAB,
and SciPy is the goal. Things to be considered include: use of ``cse`` on large
expressions for MATLAB and SciPy, which are interpretive. It is currently unclear
whether compiled languages will benefit from common subexpression elimination,
especially considering that it is a common part of compiler optimization, and
there can be a significant time penalty when calling ``cse``.

Care needs to be taken when constructing the strings for these expressions, as
well as handling of input parameters, and other dynamic symbols. How to deal
with output quantities when integrating also needs to be decided, with the
potential for multiple options being considered.

Additional Options on Initialization of Kane, RigidBody, and Particle
---------------------------------------------------------------------
This would allow a user to specify all relevant information using keyword
arguments when creating these objects. This is fairly clear for ``RigidBody``
and ``Point``. For ``Kane``, everything but the force and body lists will be
able to be entered, as computation of Fr and Fr* can take a while, and produce
an output.

Additional Methods for RigidBody and Particle
---------------------------------------------
For ``RigidBody`` and ``Particle`` (not all methods for ``Particle`` though),
add methods for getting: momentum, angular momentum, and kinetic energy.
Additionally, adding a attribute and method for defining potential energy would
allow for a total energy method/property.

Also possible is including the method which creates a transformation matrix for
3D animations; this would require a "reference orientation" for a camera as
well as a "reference point" for distance to the camera. Development of this
could also be tied into code output.


