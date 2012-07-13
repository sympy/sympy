===================
Classical Mechanics
===================

:Authors: Gilbert Gede, Luke Peterson, Angadh Nanjangud

.. topic:: Abstract

   In this documentation many components of the physics/mechanics module will
   be discussed. :mod:`mechanics` has been written to allow for creation of
   symbolic equations of motion for complicated multibody systems.

Mechanics
=========

In physics, mechanics describes conditions of rest or motion; statics or
dynamics. First, an idealized representation of a system is described. Next, we
use physical laws to generate equations that define the system's behaviors.
Then, we solve these equations, in multiple possible ways. Finally, we extract
information from these equations and solutions. Mechanics is currently set up
to create equations of motion, for dynamics. It is also limited to rigid bodies
and particles.


Guide to Mechanics
==================

.. toctree::
    :maxdepth: 2

    vectors.rst
    kinematics.rst
    masses.rst
    kane.rst
    examples.rst
    advanced.rst
    reference.rst

Mechanics API
=============

.. toctree::
    :maxdepth: 2

    api/essential.rst
    api/functions.rst
    api/kinematics.rst
    api/part_bod.rst
    api/kane.rst
    api/printing.rst
