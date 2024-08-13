==========================
Biomechanics API Reference
==========================

.. topic:: Abstract

   In this documentation many components of the physics/biomechanics module will
   be discussed. :mod:`sympy.physics.biomechanics` extends :mod:`sympy.physics.mechanics`
   with a number of function, classes, and other utilities that facilitate the
   creation of biomechanical models.

.. module:: sympy.physics.biomechanics

Introduction
============

Biomechanical models typically involve simultaneously modeling numerous aspects
of anatomy. These can be the skeletal, muscular, and neurological systems. The
skeletal system is typically modeled using multibody dynamics. This can be done
using :mod:`sympy.physics.mechanics`. :mod:`sympy.physics.biomechanics` provides
functionality for modeling musculotendons (the muscular system) and activation
dynamics (the neurological system).

Mechanics
=========

This module acts as an extension module to :mod:`sympy.physics.mechanics`. It is
not intended that this module be used by itself. Rather, a user would import
both :mod:`sympy.physics.mechanics` and :mod:`sympy.physics.biomechanics`, and
use objects from both interchangeably. :mod:`sympy.physics.biomechanics` has
been designed in such a way that its class hierarchies are related to, and
interfaces (e.g. attribute names, call signatures, and return types) mimic,
those of :mod:`sympy.physics.mechanics`. Consequentially, :mod:`sympy.physics.mechanics`
will correctly generate equations of motion for multibody systems that
incorporate biomechanical components.

Guide to Biomechanics
=====================
.. toctree::
   :titlesonly:

   musculotendon.rst
   activation.rst
   curve.rst
