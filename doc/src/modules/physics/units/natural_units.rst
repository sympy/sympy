=============
Natural units
=============

In natural units some physical constants are equal to one. The dimensions
related by these constants cannot be distinguished any more: if the speed of
light is equal to one, length and time have the same dimension, and the
velocity is dimensionless.

Available unit systems
======================

The following unit systems are defined in ``sympy.physics.units.systems``.
The constants are the speed of light `c`, the reduced Planck constant
`\hbar`, the gravitational constant `G`, the Boltzmann constant `k_B`, the
Coulomb constant `k_e = 1/(4 \pi \epsilon_0)`, the vacuum permittivity
`\epsilon_0`, the elementary charge `e`, the masses of electron and proton
`m_e` and `m_p`. The fine-structure constant is `\alpha`.

.. list-table::
   :header-rows: 1

   * - Unit system
     - Name
     - Defining constants
     - Remarks
   * - ``natural``
     - ``"Natural system"``
     - `c = \hbar = \epsilon_0 = 1`
     - Particle physics, with Heaviside-Lorentz units. Energies are measured
       in electronvolt, `e = \sqrt{4 \pi \alpha}`.
   * - ``planck_units``
     - ``"Planck"``
     - `c = G = \hbar = k_B = k_e = 1`
     - `e = \sqrt{\alpha}`
   * - ``stoney_units``
     - ``"Stoney"``
     - `c = G = k_e = e = 1`
     - `\hbar = 1/\alpha`
   * - ``schroedinger_units``
     - ``"Schroedinger"``
     - `\hbar = G = k_e = e = 1`
     - `c = 1/\alpha`
   * - ``geometrized_units``
     - ``"geometrized"``
     - `c = G = 1`
     - General relativity. Masses and times are lengths.
   * - ``hartree_atomic_units``
     - ``"Hartree"``
     - `\hbar = m_e = e = k_e = 1`
     - Atomic physics. `c = 1/\alpha`
   * - ``rydberg_atomic_units``
     - ``"Rydberg"``
     - `\hbar = 2 m_e = e^2/2 = k_e = 1`
     - Atomic physics. `c = 2/\alpha`
   * - ``strong_units``
     - ``"strong"``
     - `c = \hbar = m_p = 1`
     - Nuclear physics.

The relation of these unit systems with the other ones is described in
:doc:`systems`.

All of them are derived from the SI. The dimensions that are not related to
the defining constants are not modified: the temperature is a base dimension
in all these unit systems but Planck units, the current is a base dimension
in geometrized and strong units.

The defining constants are stored by the unit systems:

    >>> from sympy.physics.units.systems import natural, planck_units
    >>> natural.defining_constants
    {hbar: 1, speed_of_light: 1, vacuum_permittivity: 1}
    >>> planck_units.defining_constants
    {boltzmann_constant: 1, coulomb_constant: 1, gravitational_constant: 1, hbar: 1, speed_of_light: 1}

Dimensions and units
====================

In the natural units of particle physics the only mechanical base dimension
is the energy. Lengths and times are inverse energies, masses are energies:

    >>> from sympy.physics.units import length, time, mass, energy, velocity
    >>> dimsys = natural.get_dimension_system()
    >>> dimsys.get_dimensional_dependencies(length)
    {Dimension(energy, E): -1}
    >>> dimsys.equivalent_dims(mass, energy)
    True
    >>> dimsys.is_dimensionless(velocity)
    True

The function ``convert_to`` relates units that are not related in the SI:

    >>> from sympy.physics.units import convert_to, meter, second, kilogram
    >>> from sympy.physics.units import electronvolt, speed_of_light, hbar
    >>> convert_to(speed_of_light, 1, natural)
    1
    >>> convert_to(hbar, 1, natural)
    1
    >>> convert_to(1/meter, electronvolt, natural).n(6)
    1.97327e-7*electronvolt
    >>> convert_to(1/second, electronvolt, natural).n(6)
    6.58212e-16*electronvolt
    >>> convert_to(kilogram, electronvolt, natural).n(6)
    5.60959e+35*electronvolt

In Planck and atomic units all mechanical and electromagnetic quantities are
dimensionless. The value of a quantity in these units is given by the
conversion to ``1``:

    >>> from sympy.physics.units import joule, planck_mass, elementary_charge
    >>> from sympy.physics.units.systems import hartree_atomic_units
    >>> convert_to(speed_of_light, 1, hartree_atomic_units).n(6)
    137.036
    >>> convert_to(joule, 1, hartree_atomic_units).n(6)
    2.29371e+17
    >>> convert_to(kilogram, planck_mass, planck_units).n(6)
    4.59467e+7*planck_mass
    >>> convert_to(elementary_charge, 1, planck_units).n(6)
    0.0854245

In geometrized units masses are lengths:

    >>> from sympy.physics.units.systems import geometrized_units
    >>> convert_to(kilogram, meter, geometrized_units).n(6)
    7.42616e-28*meter

Creating other unit systems
===========================

The method ``contract`` of the unit systems creates the unit system where
the given constants are pure numbers. For example, the Boltzmann constant can
be set to one in the natural unit system, in order to measure the
temperatures in electronvolt:

    >>> from sympy.physics.units import boltzmann_constant, kelvin
    >>> unit_system = natural.contract([boltzmann_constant])
    >>> convert_to(kelvin, electronvolt, unit_system).n(6)
    8.61733e-5*electronvolt

The values of the constants are given by a dictionary. These are reduced
Planck units, where `8 \pi G = 1`:

    >>> from sympy import pi
    >>> from sympy.physics.units import gravitational_constant
    >>> from sympy.physics.units.systems import SI
    >>> reduced_planck_units = SI.contract({speed_of_light: 1, hbar: 1,
    ...     boltzmann_constant: 1, gravitational_constant: 1/(8*pi)})
    >>> convert_to(gravitational_constant, 1, reduced_planck_units)
    1/(8*pi)

Every constant removes a base dimension. The base dimensions of the new unit
system can be selected with units having that dimension, otherwise they are
chosen among the base dimensions of the original unit system:

    >>> unit_system = SI.contract([speed_of_light])
    >>> unit_system.get_dimension_system().equivalent_dims(time, length)
    True
    >>> length in unit_system.get_dimension_system().base_dims
    True
    >>> unit_system = SI.contract([speed_of_light], [second])
    >>> time in unit_system.get_dimension_system().base_dims
    True

The dimensions of the constants have to be independent:

    >>> from sympy.physics.units import vacuum_permittivity, coulomb_constant
    >>> SI.contract([vacuum_permittivity, coulomb_constant])
    Traceback (most recent call last):
    ...
    ValueError: the dimensions of constants and base units are not independent

Conversion of equations
=======================

The equations are converted to natural units and back with
``convert_unit_system``, see :doc:`unit_system_conversion`. The defining
constants disappear from the equations:

    >>> from sympy import symbols, Eq
    >>> from sympy.physics.units import convert_unit_system, momentum
    >>> E, p, m, n = symbols("E p m n")
    >>> dims = {E: energy, p: momentum, m: mass, n: 1}
    >>> eq = Eq(E**2, p**2*speed_of_light**2 + m**2*speed_of_light**4)
    >>> convert_unit_system(eq, dims, SI, natural)
    Eq(E**2, m**2 + p**2)

In the opposite direction they are restored. These are the energy levels of
the hydrogen atom in Hartree atomic units and in the SI:

    >>> eq = Eq(E, -1/(2*n**2))
    >>> convert_unit_system(eq, dims, hartree_atomic_units, SI)
    Eq(E, -coulomb_constant**2*elementary_charge**4*electron_rest_mass/(2*hbar**2*n**2))

The other constants are not replaced by their numerical value:

    >>> eq = Eq(E, m*speed_of_light**2)
    >>> convert_unit_system(eq, dims, SI, hartree_atomic_units)
    Eq(E, speed_of_light**2*m)

Two natural unit systems are in general not related directly, because their
defining constants are different. The conversion passes through the unit
system they are derived from:

    >>> from sympy.physics.units.systems import rydberg_atomic_units
    >>> eq = Eq(E, -1/(2*n**2))
    >>> convert_unit_system(eq, dims, hartree_atomic_units, rydberg_atomic_units)
    Eq(E, -1/n**2)
    >>> convert_unit_system(eq, dims, hartree_atomic_units, natural)
    Eq(E, -elementary_charge**4*electron_rest_mass/(32*pi**2*n**2))

Reference
=========

.. automodule:: sympy.physics.units.systems.natural
