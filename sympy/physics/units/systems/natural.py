r"""
Natural unit systems.

In natural units some physical constants are equal to one, and the dimensions
related by these constants cannot be distinguished any more. For example, if
the speed of light is equal to one, length and time have the same dimension.

The following unit systems are defined:

======================== ======================================================
``natural``              `c = \hbar = \epsilon_0 = 1`
``planck_units``         `c = G = \hbar = k_B = k_e = 1`
``stoney_units``         `c = G = k_e = e = 1`
``schroedinger_units``   `\hbar = G = k_e = e = 1`
``geometrized_units``    `c = G = 1`
``hartree_atomic_units`` `\hbar = m_e = e = k_e = 1`
``rydberg_atomic_units`` `\hbar = 2 m_e = e^2/2 = k_e = 1`
``strong_units``         `c = \hbar = m_p = 1`
======================== ======================================================

They are derived from the SI. The quantities that are not related to the
constants keep their dimension, for example the temperature is a base
dimension in all unit systems but Planck units.

In the natural unit system the electronvolt is the unit of energy, and the
electromagnetic quantities are the ones of Heaviside-Lorentz units. In Planck
units the Coulomb constant is equal to one, in order to have the Planck charge
as unit of charge.

Other ones can be created with
:meth:`~sympy.physics.units.unitsystem.UnitSystem.contract`.
"""
from __future__ import annotations

from sympy.core.singleton import S
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.physics.units.definitions import (
    boltzmann_constant, c, coulomb_constant, electron_rest_mass,
    elementary_charge, eV, gravitational_constant, hbar, meter,
    planck_charge, planck_length, planck_mass, planck_temperature,
    planck_time, proton_rest_mass, vacuum_permittivity)
from sympy.physics.units.prefixes import PREFIXES, prefix_unit
from sympy.physics.units.systems.si import SI


units = prefix_unit(eV, PREFIXES)

natural = SI.contract(
    [c, hbar, vacuum_permittivity],
    base_units=[eV], units=units, name="Natural system",
    description="Natural units of particle physics, with Heaviside-Lorentz "
                "units for the electromagnetic quantities")

planck_units = SI.contract(
    [c, gravitational_constant, hbar, boltzmann_constant, coulomb_constant],
    units=[planck_length, planck_mass, planck_time, planck_temperature, planck_charge],
    name="Planck", description="Planck units")

stoney_units = SI.contract(
    [c, gravitational_constant, coulomb_constant, elementary_charge],
    name="Stoney", description="Stoney units")

schroedinger_units = SI.contract(
    [hbar, gravitational_constant, coulomb_constant, elementary_charge],
    name="Schroedinger", description="Schroedinger units")

geometrized_units = SI.contract(
    [c, gravitational_constant],
    base_units=[meter], name="geometrized", description="Geometrized units")

hartree_atomic_units = SI.contract(
    [hbar, electron_rest_mass, elementary_charge, coulomb_constant],
    name="Hartree", description="Hartree atomic units")

rydberg_atomic_units = SI.contract(
    {hbar: 1, electron_rest_mass: S.Half, elementary_charge: sqrt(2), coulomb_constant: 1},
    name="Rydberg", description="Rydberg atomic units")

strong_units = SI.contract(
    [c, hbar, proton_rest_mass],
    name="strong", description="Strong units")


__all__ = [
    'natural', 'planck_units', 'stoney_units', 'schroedinger_units',
    'geometrized_units', 'hartree_atomic_units', 'rydberg_atomic_units',
    'strong_units', 'units',
]
