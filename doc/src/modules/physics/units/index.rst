============
Unit Systems
============

This module integrates unit systems into SymPy, allowing a user choose which
system to use when doing their computations and providing utilities to display
and convert units.

Units (like meters, pounds, seconds) and constants (like light years,
Boltzmann's constant) are all considered quantities. A ``Quantity`` object
defines both units and physical constants (though its subclass
``PhysicalConstant`` may be preferred for physical constants).

The relations between quantities are defined by their dimensions and the scale
factor to at least another quantity of the same dimension.  These two types of
relations are usually defined inside ``UnitSystem`` objects, except for
properties valid in every unit system. For example, 1 kilometer is equal to
1000 meters in all unit systems and its dimension is ``length`` in all
dimension systems.  On the other hand, the speed of light is equal to 299792458
meters per second in SI units, while it is equal to 1 (unitless) in natural
units. In both SI and natural units the dimension of the speed of light in
``velocity``, but in the dimension system of natural units ``velocity`` is
dimensionless because ``length`` and ``time`` are equivalent. Similarly, there
are discrepancies in the dimensions and scale factors of electromagnetic
quantities between SI unit system and CGS and gaussian unit systems, as the
last two ones do not consider the ``current`` to be a fundamental dimension.

The advantage of this implementation over the one found in other libraries is
that it handles relations between units differently in different unit systems,
without restrictions to the assumption of relations between units and physical
constants provided by the SI unit system.

Examples
--------

The most important function in the units module is ``convert_to``, it allows
the given quantity to be rewritten as the product of powers of some target
quantities. For example, to represent the speed of light in terms of meters and
seconds:

>>> from sympy.physics.units import speed_of_light, meter, second
>>> from sympy.physics.units import convert_to
>>> convert_to(speed_of_light, [meter, second])
299792458*meter/second

If it is not possible to represent the given quantity in the target units, the
given quantity will be returned unchanged:

>>> convert_to(speed_of_light, [meter])
speed_of_light

The relations between quantities depend on the unit systems. So, ``convert_to``
accepts an optional third parameter representing the unit system, which is
``SI`` by default. The conversion may return different results depending on the
chosen unit system, for example, in the ``cgs_gauss`` unit system the current
is not a fundamental dimension, rather it can be represented as a combination
of length, time and mass:

>>> from sympy.physics.units.systems.si import SI
>>> from sympy.physics.units.systems.cgs import cgs_gauss
>>> from sympy.physics.units import ampere, gram, second
>>> convert_to(ampere, [meter, gram, second], SI)
ampere
>>> convert_to(ampere, [meter, gram, second], cgs_gauss)
149896229*sqrt(gram)*meter**(3/2)/(50*second**2)

Quantities of the same dimension do not get simplified automatically, for
example if you divide meters by kilometers, you will get an object representing
the division between the two units. In order to simplify this kind of
expressions, you can either call the ``.simplify()`` method or import the
``quantity_simplify( )`` function, the last one also accepting a unit system as
optional parameter.

>>> from sympy.physics.units.util import quantity_simplify
>>> from sympy.physics.units import kilometer
>>> meter/kilometer
meter/kilometer
>>> (meter/kilometer).simplify()
1/1000
>>> quantity_simplify(meter/kilometer)
1/1000

More
----

Ideas about future developments can be found on the `Github wiki
<https://github.com/sympy/sympy/wiki/Unit-systems>`_.

.. toctree::
   :titlesonly:

   philosophy.rst
   examples.rst
   dimensions.rst
   prefixes.rst
   unitsystem.rst
   quantities.rst
