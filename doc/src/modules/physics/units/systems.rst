===================================
Unit systems and dimension systems
===================================

This page explains how dimensions, dimension systems, quantities and unit
systems are related, where the information is stored, and how the systems
defined in SymPy are derived from each other.

The objects
===========

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Object
     - Role
   * - ``Dimension``
     - A name, like ``length`` or ``charge``. It does not know how it is
       related to the other dimensions.
   * - ``DimensionSystem``
     - A choice of base dimensions, with the expressions of the derived
       dimensions in terms of the base ones. It does not contain units.
   * - ``Quantity``
     - A unit or a physical constant, like ``meter`` or ``speed_of_light``.
       It is a name as well: dimension and magnitude are given by the unit
       system.
   * - ``UnitSystem``
     - A dimension system, together with the units: the base units and the
       scale factors relating the units of the same dimension.

The relations are:

- every unit system has exactly one dimension system, returned by
  ``get_dimension_system``;
- a dimension system does not depend on unit systems. It can be used alone,
  and unit systems sharing the same base dimensions may share the dimension
  system;
- dimensions and quantities are shared by all systems. Their meaning depends
  on the system.

    >>> from sympy.physics.units.systems.si import SI, dimsys_SI
    >>> SI.get_dimension_system() is dimsys_SI
    True

Dimensions depend on the dimension system
=========================================

The dimension ``charge`` is the same object everywhere, but its relation with
the other dimensions is not the same. In the SI it is the product of current
and time. In Gaussian units there is no base dimension for the electromagnetic
quantities, the charge is expressed in terms of length, mass and time. In
Planck units it is dimensionless:

    >>> from sympy.physics.units import charge, current, time
    >>> from sympy.physics.units.systems.cgs import cgs_gauss
    >>> from sympy.physics.units.systems import planck_units
    >>> for unit_system in [SI, cgs_gauss, planck_units]:
    ...     dimsys = unit_system.get_dimension_system()
    ...     print(dimsys.get_dimensional_dependencies(charge))
    {Dimension(current): 1, Dimension(time): 1}
    {Dimension(mass): 1/2, Dimension(length): 3/2, Dimension(time): -1}
    {}

For this reason the dimensions are compared by the dimension systems. The
equality of two dimensions only compares their names:

    >>> charge == current*time
    False
    >>> dimsys_SI.equivalent_dims(charge, current*time)
    True

A dimension which is not known to a dimension system is considered
independent of the other ones. The dimension system of the MKS unit system
does not contain the electromagnetic dimensions:

    >>> from sympy.physics.units.systems import MKS
    >>> MKS.get_dimension_system().get_dimensional_dependencies(charge)
    {Dimension(charge, Q): 1}

Quantities depend on the unit system
====================================

The dimension of a quantity is given by the unit system:

    >>> from sympy.physics.units import coulomb_constant
    >>> SI.get_quantity_dimension(coulomb_constant)
    Dimension(force*length**2/charge**2)
    >>> cgs_gauss.get_quantity_dimension(coulomb_constant)
    Dimension(1)

The magnitude of a quantity is given by its scale factor. The scale factors
are relative to a reference which is not specified, they are meaningful only
if they are compared to the scale factors of other quantities of the same
dimension, in the same unit system. This is what ``convert_to`` does:

    >>> from sympy.physics.units import coulomb, convert_to
    >>> from sympy.physics.units import centimeter, gram, second
    >>> cgs_gauss.get_quantity_scale_factor(coulomb)
    149896229/50
    >>> convert_to(coulomb, [centimeter, gram, second], cgs_gauss)
    2997924580*centimeter**(3/2)*sqrt(gram)/second

The scale factors are not relative to the base units of the unit system. The
scale factor of the kilogram is 1000 in the SI, in order to be compatible with
the prefix, despite being a base unit:

    >>> from sympy.physics.units import kilogram
    >>> SI.get_quantity_scale_factor(kilogram)
    1000

In the unit systems where some dimensions are equivalent, the scale factors
of the units of these dimensions are related:

    >>> from sympy.physics.units import meter
    >>> from sympy.physics.units.systems import geometrized_units
    >>> geometrized_units.get_quantity_scale_factor(second)
    299792458
    >>> geometrized_units.get_quantity_scale_factor(meter)
    1

.. note::

   The properties ``dimension`` and ``scale_factor`` of the quantities refer
   to the SI. The methods ``get_quantity_dimension`` and
   ``get_quantity_scale_factor`` of the unit system should be used with the
   other unit systems.

Where the information is stored
===============================

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Information
     - Stored by
     - Defined with
   * - Base and derived dimensions, dependencies of the derived dimensions
     - ``DimensionSystem``
     - Constructor, ``DimensionSystem.extend``
   * - Dimension of a quantity
     - ``DimensionSystem``, ``UnitSystem``, or globally
     - ``set_quantity_dimension``, ``Quantity.set_global_dimension``
   * - Scale factor of a quantity
     - ``DimensionSystem``, ``UnitSystem``, or globally
     - ``set_quantity_scale_factor``,
       ``Quantity.set_global_relative_scale_factor``
   * - Base units, units, units of the derived dimensions
     - ``UnitSystem``
     - Constructor, ``UnitSystem.extend``
   * - Defining constants
     - ``UnitSystem``
     - Constructor, ``UnitSystem.contract``

Dimension and scale factor of a quantity may be defined in three places. The
unit system looks for them in this order:

1. the dimension system of the unit system;
2. the unit system;
3. the global definitions, valid in all unit systems. For example, the
   kilometer is 1000 meters in all unit systems, its scale factor is relative
   to the one of the meter.

If no definition is found, the dimension has the name of the quantity and the
scale factor is one:

    >>> from sympy.physics.units import Quantity
    >>> SI.get_quantity_dimension(Quantity("my_quantity"))
    Dimension(my_quantity)

The definitions stored by a dimension system are inherited by the dimension
systems derived from it. The ones stored by a unit system are valid in that
unit system only. For example, the scale factors of the electromagnetic units
are stored by the dimension system of MKSA, they are inherited by the SI. The
Boltzmann constant is stored by the SI unit system.

The available systems
=====================

.. list-table::
   :header-rows: 1
   :widths: 22 28 25 25

   * - Unit system
     - Dimension system
     - Base dimensions
     - Base units
   * - ``MKS``
     - ``dimsys_length_weight_time``
     - length, mass, time
     - meter, kilogram, second
   * - ``MKSA``
     - ``dimsys_MKSA``
     - the ones of ``MKS``, current
     - the ones of ``MKS``, ampere
   * - ``SI``
     - ``dimsys_SI``
     - the ones of ``MKSA``, temperature, amount of substance, luminous
       intensity
     - the ones of ``MKSA``, kelvin, mole, candela
   * - ``cgs_gauss``
     - ``dimsys_cgs``
     - length, mass, time
     - centimeter, gram, second
   * - ``natural``
     - created by ``contract``
     - energy, temperature, amount of substance, luminous intensity
     - electronvolt, kelvin, mole, candela
   * - ``geometrized_units``
     - created by ``contract``
     - length, current, temperature, amount of substance, luminous intensity
     - meter, ampere, kelvin, mole, candela
   * - ``strong_units``
     - created by ``contract``
     - current, temperature, amount of substance, luminous intensity
     - ampere, kelvin, mole, candela
   * - ``stoney_units``, ``schroedinger_units``, ``hartree_atomic_units``,
       ``rydberg_atomic_units``
     - created by ``contract``
     - temperature, amount of substance, luminous intensity
     - kelvin, mole, candela
   * - ``planck_units``
     - created by ``contract``
     - amount of substance, luminous intensity
     - mole, candela

The dimension systems are defined in the modules of
``sympy.physics.units.systems``. There is a further dimension system,
``dimsys_default``, which adds the information to the base dimensions of the
SI. It is not used by any unit system.

``MKS`` and ``cgs_gauss`` have the same base dimensions, but their dimension
systems are different: the one of ``cgs_gauss`` contains the electromagnetic
dimensions as derived dimensions.

How the systems are derived
===========================

There are three ways to create a system from another one.

.. graphviz::

    digraph {
        node [fontsize=10];
        edge [fontsize=9];

        node [shape=ellipse];
        lwt [label="dimsys_length_weight_time"];
        dcgs [label="dimsys_cgs"];
        dmksa [label="dimsys_MKSA"];
        dsi [label="dimsys_SI"];
        dnat [label="created by contract"];

        node [shape=box];
        cgs_gauss; MKS; MKSA; SI;
        nat [label="natural, planck_units, ..."];

        {rank=same; lwt -> MKS [style=dashed, arrowhead=none]}
        {rank=same; cgs_gauss -> dcgs [style=dashed, arrowhead=none]; dmksa -> MKSA [style=dashed, arrowhead=none]; dcgs -> dmksa [style=invis]}
        {rank=same; dsi -> SI [style=dashed, arrowhead=none]}
        {rank=same; dnat -> nat [style=dashed, arrowhead=none]}

        lwt -> dcgs [label="extend"];
        lwt -> dmksa [label="extend"];
        dmksa -> dsi [label="extend"];
        dsi -> dnat [label="contract"];
        MKS -> MKSA [label="extend"];
        MKSA -> SI [label="extend"];
        SI -> nat [label="contract"];
    }

The boxes are unit systems, the ellipses are dimension systems. The dashed
lines connect the unit systems to their dimension systems, the arrows go from
a system to the ones derived from it.

Extension
---------

New base dimensions and new base units are added. ``DimensionSystem.extend``
creates the dimension system, ``UnitSystem.extend`` creates the unit system.
The dimension system has to be passed to ``UnitSystem.extend``:

    >>> from sympy.physics.units.systems.mksa import MKSA, dimsys_MKSA
    >>> from sympy.physics.units import temperature, kelvin
    >>> dimsys = dimsys_MKSA.extend([temperature])
    >>> unit_system = MKSA.extend([kelvin], dimension_system=dimsys)
    >>> unit_system.get_dimension_system().base_dims
    (Dimension(current, I), Dimension(length, L), Dimension(mass, M), Dimension(temperature, T), Dimension(time, T))

The relations among the existing dimensions are not modified. ``MKSA`` and
``SI`` are defined in this way.

New derived dimensions
----------------------

``DimensionSystem.extend`` is also used to add derived dimensions, leaving the
base dimensions as they are. The dimension system of Gaussian units adds the
electromagnetic dimensions to the mechanical ones:

    >>> from sympy.physics.units.systems.cgs import dimsys_cgs
    >>> from sympy.physics.units.systems.length_weight_time import dimsys_length_weight_time
    >>> dimsys_cgs.base_dims == dimsys_length_weight_time.base_dims
    True
    >>> dimsys_cgs.get_dimensional_dependencies(current)
    {Dimension(length): 3/2, Dimension(mass): 1/2, Dimension(time): -2}

Contraction
-----------

``UnitSystem.contract`` creates the unit system where some physical constants
are pure numbers. Every constant removes a base dimension. Both the dimension
system and the scale factors are computed from the ones of the original unit
system, see :doc:`natural_units`.

If the Coulomb constant of MKSA is set to one, the current is not a base
dimension any more. The dimension of the charge is the one of Gaussian units:

    >>> unit_system = MKSA.contract([coulomb_constant])
    >>> dimsys = unit_system.get_dimension_system()
    >>> dimsys.base_dims
    (Dimension(length, L), Dimension(mass, M), Dimension(time, T))
    >>> dimsys.get_dimensional_dependencies(charge)
    {Dimension(length, L): 3/2, Dimension(mass, M): 1/2, Dimension(time, T): -1}

The relations among the dimensions of the original unit system are still
valid. This is the reason why Gaussian units are not defined by a
contraction: their magnetic field has the dimension of the electric field,
while in the SI the ratio of electric and magnetic field is a velocity.

    >>> from sympy.physics.units import voltage, length, magnetic_density, velocity
    >>> dimsys_SI.equivalent_dims(voltage/length, magnetic_density*velocity)
    True
    >>> dimsys.equivalent_dims(voltage/length, magnetic_density*velocity)
    True
    >>> dimsys_cgs.equivalent_dims(voltage/length, magnetic_density)
    True

Relations between two unit systems
==================================

The dimension systems determine whether the equations of a unit system can be
converted to another one with ``convert_unit_system``, see
:doc:`unit_system_conversion`. The base dimensions of one of the unit systems
have to be independent dimensions in the other one:

- ``MKS``, ``MKSA`` and ``SI`` are related by extensions, the equations are
  the same in all of them;
- ``cgs_gauss`` and the unit systems created by ``contract`` have fewer
  independent dimensions than ``SI``. The equations lose some constants in
  the conversion from ``SI``, they get them back in the opposite direction;
- two unit systems created by ``contract`` are related through the unit
  system they are derived from.
