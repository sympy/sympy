from __future__ import annotations
from sympy.physics.units import DimensionSystem, joule, second, ampere

from sympy.core.numbers import Float, Rational, pi
from sympy.core.singleton import S
from sympy.core.sympify import sympify
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.physics.units.definitions import (
    c, kg, m, s, boltzmann_constant, coulomb, coulomb_constant,
    electron_rest_mass, electronvolt, elementary_charge, farad,
    gravitational_constant, hbar, henry, kilogram, magnetic_constant, meter,
    ohm, planck_charge, planck_length, planck_mass, planck_temperature,
    planck_time, proton_rest_mass, speed_of_light, tesla, vacuum_impedance,
    vacuum_permittivity, volt)
from sympy.physics.units.definitions.dimension_definitions import (
    charge, energy, length, mass, time, velocity)
from sympy.physics.units.quantities import Quantity
from sympy.physics.units.systems import (
    MKSA, SI, natural, planck_units, stoney_units, schroedinger_units,
    geometrized_units, hartree_atomic_units, rydberg_atomic_units,
    strong_units)
from sympy.physics.units.unitsystem import UnitSystem
from sympy.physics.units.util import convert_to
from sympy.testing.pytest import raises


def test_definition():
    # want to test if the system can have several units of the same dimension
    dm = Quantity("dm")
    base = (m, s)
    # base_dim = (m.dimension, s.dimension)
    ms = UnitSystem(base, (c, dm), "MS", "MS system")
    ms.set_quantity_dimension(dm, length)
    ms.set_quantity_scale_factor(dm, Rational(1, 10))

    assert set(ms._base_units) == set(base)
    assert set(ms._units) == {m, s, c, dm}
    # assert ms._units == DimensionSystem._sort_dims(base + (velocity,))
    assert ms.name == "MS"
    assert ms.descr == "MS system"


def test_str_repr():
    assert str(UnitSystem((m, s), name="MS")) == "MS"
    assert str(UnitSystem((m, s))) == "UnitSystem((meter, second))"

    assert repr(UnitSystem((m, s))) == "<UnitSystem: (%s, %s)>" % (m, s)


def test_convert_to():
    A = Quantity("A")
    A.set_global_relative_scale_factor(S.One, ampere)

    Js = Quantity("Js")
    Js.set_global_relative_scale_factor(S.One, joule*second)

    mksa = UnitSystem((m, kg, s, A), (Js,))
    assert convert_to(Js, mksa._base_units) == m**2*kg*s**-1


def test_extend():
    ms = UnitSystem((m, s), (c,))
    Js = Quantity("Js")
    Js.set_global_relative_scale_factor(1, joule*second)
    mks = ms.extend((kg,), (Js,))

    res = UnitSystem((m, s, kg), (c, Js))
    assert set(mks._base_units) == set(res._base_units)
    assert set(mks._units) == set(res._units)


def test_dim():
    dimsys = UnitSystem((m, kg, s), (c,))
    assert dimsys.dim == 3


def test_is_consistent():
    dimension_system = DimensionSystem([length, time])
    us = UnitSystem([m, s], dimension_system=dimension_system)
    assert us.is_consistent == True


def test_get_units_non_prefixed():
    from sympy.physics.units import volt, ohm
    unit_system = UnitSystem.get_unit_system("SI")
    units = unit_system.get_units_non_prefixed()
    for prefix in ["giga", "tera", "peta", "exa", "zetta", "yotta", "kilo", "hecto", "deca", "deci", "centi", "milli", "micro", "nano", "pico", "femto", "atto", "zepto", "yocto"]:
        for unit in units:
            assert isinstance(unit, Quantity), f"{unit} must be a Quantity, not {type(unit)}"
            assert not unit.is_prefixed, f"{unit} is marked as prefixed"
            assert not unit.is_physical_constant, f"{unit} is marked as physics constant"
            assert not unit.name.name.startswith(prefix), f"Unit {unit.name} has prefix {prefix}"
    assert volt in units
    assert ohm in units

def test_derived_units_must_exist_in_unit_system():
    for unit_system in UnitSystem._unit_systems.values():
        for preferred_unit in unit_system.derived_units.values():
            units = preferred_unit.atoms(Quantity)
            for unit in units:
                assert unit in unit_system._units, f"Unit {unit} is not in unit system {unit_system}"


def test_mksa():
    assert convert_to(volt, [kilogram, meter, second, ampere], MKSA) == \
        kilogram*meter**2/(ampere*second**3)
    assert convert_to(farad*volt, coulomb, MKSA) == coulomb
    for quantity, units in [
            (vacuum_impedance, ohm), (vacuum_permittivity, farad/meter),
            (magnetic_constant, henry/meter), (elementary_charge, coulomb)]:
        assert convert_to(quantity, units, MKSA) == convert_to(quantity, units, SI)
        assert convert_to(quantity, units, MKSA) != quantity


def test_contract():
    unit_system = SI.contract([speed_of_light], [meter])
    dimsys = unit_system.get_dimension_system()
    assert unit_system.defining_constants == {speed_of_light: 1}
    assert unit_system.is_consistent
    assert length in dimsys.base_dims
    assert time not in dimsys.base_dims
    assert dimsys.equivalent_dims(time, length)
    assert dimsys.equivalent_dims(energy, mass)
    assert dimsys.is_dimensionless(velocity)
    assert convert_to(second, meter, unit_system) == 299792458*meter
    assert convert_to(joule, kilogram, unit_system) == kilogram/299792458**2
    assert convert_to(speed_of_light, 1, unit_system) == 1
    assert convert_to(volt, [kilogram, meter, ampere], unit_system) == \
        kilogram/(299792458**3*ampere*meter)

    unit_system2 = unit_system.contract({hbar: 2})
    dimsys2 = unit_system2.get_dimension_system()
    assert unit_system2.defining_constants == {speed_of_light: 1, hbar: 2}
    assert dimsys2.equivalent_dims(energy, 1/length)
    assert convert_to(speed_of_light, 1, unit_system2) == 1
    assert convert_to(hbar, 1, unit_system2) == 2
    assert convert_to(second, meter, unit_system2) == 299792458*meter
    assert abs(convert_to(1/meter, electronvolt, unit_system2)/electronvolt - 9.86634902296511e-8) < 1e-20

    dimsys = MKSA.contract([coulomb_constant]).get_dimension_system()
    assert list(dimsys.base_dims) == [length, mass, time]
    assert dimsys.get_dimensional_dependencies(charge) == {mass: S.Half, length: S(3)/2, time: -1}
    dimsys = SI.contract([speed_of_light]).get_dimension_system()
    assert length in dimsys.base_dims
    assert time not in dimsys.base_dims
    dimsys = SI.contract([speed_of_light], [s]).get_dimension_system()
    assert time in dimsys.base_dims
    assert length not in dimsys.base_dims

    raises(ValueError, lambda: SI.contract([vacuum_permittivity, coulomb_constant]))
    raises(ValueError, lambda: SI.contract([speed_of_light, hbar], [meter, second]))
    raises(ValueError, lambda: SI.contract([speed_of_light], [speed_of_light]))


def _value(quantity, unit_system):
    return convert_to(quantity, 1, unit_system)


def test_natural_unit_systems():
    alpha = 0.0072973525693
    constants = [
        speed_of_light, hbar, gravitational_constant, boltzmann_constant,
        coulomb_constant, vacuum_permittivity, elementary_charge,
        electron_rest_mass, proton_rest_mass]
    for unit_system, values in [
            (natural, [1, 1, None, None, 1/(4*pi), 1, sqrt(4*pi*alpha), None, None]),
            (planck_units, [1, 1, 1, 1, 1, 1/(4*pi), sqrt(alpha), None, None]),
            (stoney_units, [1, 1/alpha, 1, None, 1, 1/(4*pi), 1, None, None]),
            (schroedinger_units, [1/alpha, 1, 1, None, 1, 1/(4*pi), 1, None, None]),
            (geometrized_units, [1, None, 1, None, None, None, None, None, None]),
            (hartree_atomic_units, [1/alpha, 1, None, None, 1, 1/(4*pi), 1, 1, 1836.15267343]),
            (rydberg_atomic_units, [2/alpha, 1, None, None, 1, 1/(4*pi), sqrt(2), S.Half, 1836.15267343/2]),
            (strong_units, [1, 1, None, None, None, None, None, 1/1836.15267343, 1])]:
        assert unit_system.is_consistent
        dimsys = unit_system.get_dimension_system()
        for constant, value in zip(constants, values):
            if value is None:
                continue
            assert dimsys.is_dimensionless(unit_system.get_quantity_dimension(constant))
            if constant in unit_system.defining_constants or not sympify(value).has(Float):
                assert _value(constant, unit_system) == value
            else:
                assert abs(_value(constant, unit_system)/value - 1) < 1e-8

    dimsys = natural.get_dimension_system()
    assert energy in dimsys.base_dims
    assert dimsys.equivalent_dims(mass, energy)
    assert dimsys.equivalent_dims(length, 1/energy)
    assert dimsys.equivalent_dims(time, 1/energy)
    assert dimsys.is_dimensionless(charge)
    assert abs(convert_to(1/meter, electronvolt, natural)/electronvolt - 1.97326980459302e-7) < 1e-20
    assert abs(convert_to(kilogram, electronvolt, natural)/electronvolt - 5.60958860380445e+35) < 1e22
    assert abs(convert_to(tesla, electronvolt**2, natural)/electronvolt**2 - 195.35) < 1e-2

    for unit in [planck_length, planck_mass, planck_time, planck_temperature, planck_charge]:
        assert abs(_value(unit, planck_units) - 1) < 1e-12
    assert abs(convert_to(kilogram, planck_mass, planck_units)/planck_mass - 45946710.3) < 1

    dimsys = geometrized_units.get_dimension_system()
    assert dimsys.equivalent_dims(mass, length)
    assert abs(convert_to(kilogram, meter, geometrized_units)/meter - 7.42616e-28) < 1e-33

    assert abs(_value(joule, hartree_atomic_units)*4.3597447222071e-18 - 1) < 1e-8
    assert abs(_value(meter, hartree_atomic_units)*5.29177210903e-11 - 1) < 1e-8
    assert abs(_value(joule, rydberg_atomic_units)*4.3597447222071e-18 - 2) < 1e-8
    assert abs(_value(meter, rydberg_atomic_units)*5.29177210903e-11 - 1) < 1e-8
