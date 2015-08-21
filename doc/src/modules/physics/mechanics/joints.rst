====================================================
Body, Joints and Joint's Method in Physics/Mechanics
====================================================

:mod:`mechanics` provides functionality for deriving equations of motion
using Kane's method [Kane1985]_. This document will describe Kane's method
as used in this module, but not how the equations are actually derived.

Body
====

Bodies are created with the class ``Body`` in :mod:`mechanics`. A Body object
is a common representation of a RigidBody or a Particle in classical mechanics.
A ``Body`` object has an associated name, center of mass, mass and reference
frame.

    >>> from sympy.physics.mechanics import Body
    >>> body = Body('body')

It may of may not have inertia depending upon if it is a common represenation
of a RigidBody or not. It provides common functionalities over RigidBody
and Particle classes like its not necessary to pass all the attributes of
RigidBody. If not passed, it defines defaults and use them. It also associates
a frame with Particle too and provides public functions to apply force and
torque to the ``Body``.

