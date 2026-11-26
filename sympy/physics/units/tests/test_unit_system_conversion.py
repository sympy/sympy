from __future__ import annotations

from sympy.core.function import Derivative, Function
from sympy.core.numbers import I, pi
from sympy.core.relational import Eq
from sympy.core.symbol import symbols
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.integrals.integrals import Integral
from sympy.matrices.dense import Matrix
from sympy.physics.units import (
    DimensionSystem, UnitSystem, convert_to, convert_unit_system)
from sympy.physics.units.definitions import (
    ampere, boltzmann_constant, coulomb, coulomb_constant, electron_rest_mass,
    electronvolt, elementary_charge, gravitational_constant, hbar, joule,
    kelvin, kilogram, magnetic_constant, meter, newton, proton_rest_mass,
    second, speed_of_light, tesla, vacuum_impedance, vacuum_permittivity,
    volt)
from sympy.physics.units.definitions.dimension_definitions import (
    capacitance, charge, current, energy, force, impedance, inductance,
    length, magnetic_density, mass, momentum, temperature, time, velocity,
    voltage, volume)
from sympy.physics.units.systems import (
    MKS, MKSA, SI, natural, planck_units, stoney_units, schroedinger_units,
    geometrized_units, hartree_atomic_units, rydberg_atomic_units,
    strong_units)
from sympy.physics.units.systems.cgs import cgs_gauss
from sympy.testing.pytest import raises

F, W, U, R, L, C, T = symbols("F W U R L C T")
q, q1, q2, i, m, p, r, v, x, t, omega, k = symbols("q q1 q2 i m p r v x t omega k")
E, B, J, rho = symbols("E B J rho", cls=Function)

dims = {
    F: force, W: energy, U: voltage, R: impedance, L: inductance,
    C: capacitance, T: temperature, q: charge, q1: charge, q2: charge,
    i: current, m: mass, p: momentum, r: length, v: velocity, x: length,
    t: time, omega: 1/time, k: 1/length, E: voltage/length,
    B: magnetic_density, J: current/length**2, rho: charge/volume,
}

e0 = vacuum_permittivity
u0 = magnetic_constant
k_e = coulomb_constant
c = speed_of_light


def _check(equation_si, equation_gauss, equation_si_back=None):
    if equation_si_back is None:
        equation_si_back = equation_si
    assert convert_unit_system(equation_si, dims, SI, cgs_gauss) == equation_gauss
    assert convert_unit_system(equation_si, dims, MKSA, cgs_gauss) == equation_gauss
    assert convert_unit_system(equation_gauss, dims, cgs_gauss, SI) == equation_si_back
    assert convert_unit_system(equation_gauss, dims, cgs_gauss, MKSA) == equation_si_back


def test_si_gauss():
    _check(Eq(F, q1*q2/(4*pi*e0*r**2)), Eq(F, q1*q2/r**2), Eq(F, k_e*q1*q2/r**2))
    _check(Eq(F, k_e*q1*q2/r**2), Eq(F, q1*q2/r**2))
    _check(Eq(F, q*(E(x, t) + v*B(x, t))), Eq(F, q*(E(x, t) + v*B(x, t)/c)))
    _check(Eq(W, m*c**2), Eq(W, m*c**2))
    _check(Eq(U, R*i), Eq(U, R*i))
    _check(Eq(q, C*U), Eq(q, C*U))
    _check(Eq(W, L*i**2/2), Eq(W, L*i**2/2))
    _check(Eq(U, Integral(E(x, t), (x, 0, r))), Eq(U, Integral(E(x, t), (x, 0, r))))
    _check(
        Eq(W, exp(-q*U/(boltzmann_constant*T))*e0*E(x, t)**2/2 + B(x, t)**2/(2*u0)),
        Eq(W, exp(-q*U/(boltzmann_constant*T))*E(x, t)**2/(8*pi) + B(x, t)**2/(8*pi)),
        Eq(W, exp(-q*U/(boltzmann_constant*T))*E(x, t)**2/(8*pi*k_e) + c**2*B(x, t)**2/(8*pi*k_e)))


def test_si_gauss_maxwell_equations():
    _check(
        Eq(E(x, t).diff(x), rho(x, t)/e0),
        Eq(E(x, t).diff(x), 4*pi*rho(x, t)),
        Eq(E(x, t).diff(x), 4*pi*k_e*rho(x, t)))
    _check(
        Eq(E(x, t).diff(x), -B(x, t).diff(t)),
        Eq(E(x, t).diff(x), -B(x, t).diff(t)/c))
    _check(
        Eq(B(x, t).diff(x), u0*J(x, t) + u0*e0*E(x, t).diff(t)),
        Eq(B(x, t).diff(x), 4*pi*J(x, t)/c + E(x, t).diff(t)/c),
        Eq(B(x, t).diff(x), 4*pi*k_e*J(x, t)/c**2 + E(x, t).diff(t)/c**2))
    _check(
        Eq(0, rho(x, t).diff(t) + J(x, t).diff(x)),
        Eq(0, rho(x, t).diff(t) + J(x, t).diff(x)))


def test_si_gauss_constants():
    eq = Eq(F, q1*q2/r**2)
    assert convert_unit_system(eq, dims, cgs_gauss, SI, constants=[e0, c]) == \
        Eq(F, q1*q2/(4*pi*e0*r**2))
    eq = Eq(B(x, t).diff(x), 4*pi*J(x, t)/c + E(x, t).diff(t)/c)
    assert convert_unit_system(eq, dims, cgs_gauss, SI, constants=[e0, c]) == \
        Eq(B(x, t).diff(x), J(x, t)/(e0*c**2) + E(x, t).diff(t)/c**2)
    eq = Eq(F, q*(E(x, t) + v*B(x, t)))
    assert convert_unit_system(eq, dims, SI, cgs_gauss, constants=[e0, c]) == \
        Eq(F, q*(E(x, t) + v*B(x, t)/c))

    raises(ValueError, lambda: convert_unit_system(eq, dims, SI, cgs_gauss, constants=[u0, c]))
    raises(ValueError, lambda: convert_unit_system(eq, dims, SI, cgs_gauss, constants=[e0, k_e, c]))
    raises(ValueError, lambda: convert_unit_system(eq, dims, SI, cgs_gauss, constants=[e0]))

    alpha = elementary_charge**2/(4*pi*e0*hbar*c)
    assert convert_unit_system(alpha, {}, SI, cgs_gauss) == elementary_charge**2/(hbar*c)
    assert convert_unit_system(elementary_charge**2/(hbar*c), {}, cgs_gauss, SI) == \
        k_e*elementary_charge**2/(hbar*c)
    assert convert_unit_system(vacuum_impedance, {}, SI, cgs_gauss) == 4*pi/c
    assert convert_unit_system(
        Eq(vacuum_impedance, 1/(e0*c)), {}, SI, cgs_gauss, dimension=impedance) == True
    assert convert_unit_system(Eq(W, 5*coulomb*U), dims, SI, cgs_gauss) == Eq(W, 5*coulomb*U)


def test_si_gauss_expressions():
    assert convert_unit_system(q1*q2/(4*pi*e0*r**2), dims, SI, cgs_gauss) == q1*q2/r**2
    assert convert_unit_system(q1*q2/r**2, dims, cgs_gauss, SI) == k_e*q1*q2/r**2

    expr = u0*i/(2*pi*r)
    assert convert_unit_system(expr, dims, SI, cgs_gauss) == 2*i/(c**2*r)
    assert convert_unit_system(expr, dims, SI, cgs_gauss, dimension=magnetic_density) == 2*i/(c*r)
    assert convert_unit_system(2*i/(c*r), dims, cgs_gauss, SI, dimension=magnetic_density) == \
        2*k_e*i/(c**2*r)

    assert convert_unit_system(q/r, dims, cgs_gauss, SI, dimension=voltage) == k_e*q/r
    assert convert_unit_system(q/(4*pi*e0*r), dims, SI, cgs_gauss, dimension=voltage) == q/r

    assert convert_unit_system([Eq(W, q*U), q**2/r], dims, cgs_gauss, SI) == \
        [Eq(W, q*U), k_e*q**2/r]
    assert convert_unit_system(Matrix([q*U, q**2/r]), dims, cgs_gauss, SI) == \
        Matrix([q*U, k_e*q**2/r])
    assert convert_unit_system(Eq(W, q**2/r), dims, "cgs_gauss", "SI") == Eq(W, k_e*q**2/r)

    eq = Eq(W, Piecewise((q**2/r, r > x), (q*U, True)))
    assert convert_unit_system(eq, dims, cgs_gauss, SI) == \
        Eq(W, Piecewise((k_e*q**2/r, r > x), (q*U, True)))
    assert convert_unit_system(Eq(W, Integral(q/r, (q, 0, q1))), dims, cgs_gauss, SI) == \
        Eq(W, k_e*Integral(q/r, (q, 0, q1)))
    assert convert_unit_system(Eq(U, Derivative(W, q)), dims, cgs_gauss, SI) == \
        Eq(U, Derivative(W, q))
    assert convert_unit_system(Eq(U, Derivative(q**2/r, q)), dims, cgs_gauss, SI) == \
        Eq(U, k_e*Derivative(q**2/r, q))

    dims2 = dict(dims)
    del dims2[E]
    dims2[E(x, t)] = voltage/length
    assert convert_unit_system(Eq(E(x, t), q/x**2), dims2, cgs_gauss, SI) == \
        Eq(E(x, t), k_e*q/x**2)


def test_si_gauss_numerical_values():
    for expr, dimension, values, unit in [
            (q1*q2/(4*pi*e0*r**2), force, {q1: 3*coulomb, q2: 5*coulomb, r: 2*meter}, newton),
            (q/(4*pi*e0*r), voltage, {q: 3*coulomb, r: 2*meter}, volt),
            (u0*i/(2*pi*r), magnetic_density, {i: 3*ampere, r: 2*meter}, tesla),
            (e0*U**2/(2*r**2) + U**2/(2*u0*v**2*r**2), energy/volume,
             {U: 3*volt, r: 2*meter, v: 5*meter/second}, newton/meter**2)]:
        expr_gauss = convert_unit_system(expr, dims, SI, cgs_gauss, dimension=dimension)
        assert expr_gauss != expr
        assert convert_to(expr_gauss.subs(values), unit, cgs_gauss) == \
            convert_to(expr.subs(values), unit, SI)
        expr_si = convert_unit_system(expr_gauss, dims, cgs_gauss, SI, dimension=dimension)
        assert convert_to(expr_si.subs(values), unit, SI) == \
            convert_to(expr.subs(values), unit, SI)


def test_equivalent_unit_systems():
    eq = Eq(W, exp(-q*U/(boltzmann_constant*T)) + m*c**2 + e0*E(x, t)**2*r**3)
    for source, target in [(SI, MKSA), (MKSA, SI), (SI, MKS), (MKS, SI)]:
        assert convert_unit_system(eq, dims, source, target) == eq


def test_dimensionless_constants():
    dimsys = DimensionSystem([energy], dimensional_dependencies={
        "mass": {"energy": 1},
        "length": {"energy": -1},
        "time": {"energy": -1},
        "momentum": {"energy": 1},
        "velocity": {},
        "action": {},
    })
    unit_system = UnitSystem(
        base_units=[electronvolt], name="test_dimensionless_constants",
        dimension_system=dimsys)
    for constant in [speed_of_light, hbar]:
        unit_system.set_quantity_dimension(constant, 1)
        unit_system.set_quantity_scale_factor(constant, 1)

    for other in [MKS, SI, cgs_gauss]:
        for eq1, eq2 in [
                (Eq(W**2, p**2*c**2 + m**2*c**4), Eq(W**2, p**2 + m**2)),
                (Eq(W, hbar*omega), Eq(W, omega)),
                (Eq(p, hbar*k), Eq(p, k)),
                (exp(I*(k*x - omega*t)), exp(I*(k*x - omega*t))),
                (Eq(v, p*c**2/W), Eq(v, p/W))]:
            assert convert_unit_system(eq1, dims, other, unit_system) == eq2
            assert convert_unit_system(eq2, dims, unit_system, other) == eq1

    raises(ValueError, lambda: convert_unit_system(Eq(W, q*U), dims, SI, unit_system))


def test_unrelated_unit_systems():
    dimsys = DimensionSystem([mass, time, current], dimensional_dependencies={
        "length": {"time": 1},
        "velocity": {},
    })
    unit_system = UnitSystem(
        base_units=[kilogram, second, ampere], name="test_unrelated_unit_systems",
        dimension_system=dimsys)
    unit_system.set_quantity_dimension(speed_of_light, 1)
    unit_system.set_quantity_scale_factor(speed_of_light, 1)
    assert convert_unit_system(Eq(x, c*t), dims, SI, unit_system) == Eq(x, t)
    assert convert_unit_system(Eq(x, t), dims, unit_system, SI) == Eq(x, c*t)
    raises(NotImplementedError, lambda: convert_unit_system(Eq(x, c*t), dims, cgs_gauss, unit_system))


def test_natural_unit_systems():
    G = gravitational_constant
    e = elementary_charge
    me = electron_rest_mass
    kB = boltzmann_constant
    M, n = symbols("M n")
    psi = Function("psi")
    dims2 = {**dims, M: mass, n: 1, psi: 1/sqrt(volume)}

    energy_momentum = Eq(W**2, p**2*c**2 + m**2*c**4)
    gravity = Eq(W, G*m*M/r + q**2/(4*pi*e0*r) + hbar*omega)
    hawking = Eq(T, hbar*c**3/(8*pi*G*M*kB))
    schroedinger = Eq(
        W*psi(r), -hbar**2/(2*me)*psi(r).diff(r, 2) - e**2/(4*pi*e0*r)*psi(r))
    schroedinger_back = Eq(
        W*psi(r), -hbar**2/(2*me)*psi(r).diff(r, 2) - k_e*e**2/r*psi(r))
    levels = Eq(W, -me*e**4/(2*(4*pi*e0)**2*hbar**2*n**2))
    levels_back = Eq(W, -me*k_e**2*e**4/(2*hbar**2*n**2))
    alpha = e**2/(4*pi*e0*hbar*c)

    for unit_system, eq_si, eq, eq_back in [
            (natural, energy_momentum, Eq(W**2, p**2 + m**2), None),
            (natural, alpha, e**2/(4*pi), None),
            (natural, Eq(W, kB*T), Eq(W, kB*T), None),
            (natural, Eq(F, q*(E(x, t) + v*B(x, t))), Eq(F, q*(E(x, t) + v*B(x, t))), None),
            (planck_units, energy_momentum, Eq(W**2, p**2 + m**2), None),
            (planck_units, hawking, Eq(T, 1/(8*pi*M)), None),
            (planck_units, alpha, e**2, k_e*e**2/(hbar*c)),
            (stoney_units, gravity, Eq(W, m*M/r + q**2/r + hbar*omega),
             Eq(W, G*m*M/r + k_e*q**2/r + hbar*omega)),
            (stoney_units, alpha, 1/hbar, k_e*e**2/(hbar*c)),
            (schroedinger_units, gravity, Eq(W, m*M/r + q**2/r + omega),
             Eq(W, G*m*M/r + k_e*q**2/r + hbar*omega)),
            (schroedinger_units, alpha, 1/c, k_e*e**2/(hbar*c)),
            (geometrized_units, Eq(r, 2*G*M/c**2), Eq(r, 2*M), None),
            (geometrized_units, energy_momentum, Eq(W**2, p**2 + m**2), None),
            (hartree_atomic_units, schroedinger,
             Eq(W*psi(r), -psi(r).diff(r, 2)/2 - psi(r)/r), schroedinger_back),
            (hartree_atomic_units, levels, Eq(W, -1/(2*n**2)), levels_back),
            (hartree_atomic_units, Eq(W, m*c**2), Eq(W, m*c**2), None),
            (hartree_atomic_units, alpha, 1/c, k_e*e**2/(hbar*c)),
            (rydberg_atomic_units, schroedinger,
             Eq(W*psi(r), -psi(r).diff(r, 2) - 2*psi(r)/r), schroedinger_back),
            (rydberg_atomic_units, levels, Eq(W, -1/n**2), levels_back),
            (rydberg_atomic_units, alpha, 2/c, k_e*e**2/(hbar*c)),
            (strong_units, Eq(W, proton_rest_mass*c**2 + hbar*omega), Eq(W, 1 + omega), None),
            ]:
        if eq_back is None:
            eq_back = eq_si
        assert convert_unit_system(eq_si, dims2, SI, unit_system) == eq
        assert convert_unit_system(eq, dims2, unit_system, SI) == eq_back

    assert convert_unit_system(
        Eq(r, 4*pi*e0*hbar**2/(me*e**2)), dims2, SI, hartree_atomic_units) == Eq(r, 1)
    assert convert_unit_system(Eq(r, 1), dims2, hartree_atomic_units, SI) == \
        Eq(r, hbar**2/(k_e*me*e**2))


def test_natural_unit_systems_numerical_values():
    G = gravitational_constant
    M = symbols("M")
    dims2 = {**dims, M: mass}
    values = {
        m: 3*kilogram, M: 5*kilogram, r: 2*meter, q: 7*coulomb,
        omega: 11/second, i: 13*ampere, T: 17*kelvin}
    for unit_system in [
            cgs_gauss, natural, planck_units, stoney_units, schroedinger_units,
            geometrized_units, hartree_atomic_units, rydberg_atomic_units,
            strong_units]:
        for expr in [
                G*m*M/r, q**2/(4*pi*e0*r), hbar*omega, m*c**2, u0*i**2*r,
                boltzmann_constant*T]:
            if unit_system == cgs_gauss and expr.has(T):
                continue
            converted = convert_unit_system(expr, dims2, SI, unit_system)
            value = convert_to(converted.subs(values), joule, unit_system)/joule
            expected = convert_to(expr.subs(values), joule, SI)/joule
            assert abs(value/expected - 1) < 1e-12


def test_natural_unit_systems_intermediate():
    G = gravitational_constant
    M = symbols("M")
    dims2 = {**dims, M: mass}

    eq = Eq(W, G*m*M/r + q**2/r + hbar*omega + m*c**2)
    assert convert_unit_system(eq, dims2, cgs_gauss, planck_units) == \
        Eq(W, m*M/r + q**2/r + omega + m)
    assert convert_unit_system(eq, dims2, cgs_gauss, natural) == \
        Eq(W, G*m*M/r + q**2/(4*pi*r) + omega + m)
    assert convert_unit_system(eq, dims2, cgs_gauss, hartree_atomic_units) == \
        Eq(W, G*m*M/r + q**2/r + omega + m*c**2)
    assert convert_unit_system(eq, dims2, cgs_gauss, geometrized_units) == \
        Eq(W, m*M/r + k_e*q**2/r + hbar*omega + m)

    eq = Eq(W, m*M/r + q**2/r + omega + m)
    assert convert_unit_system(eq, dims2, planck_units, cgs_gauss) == \
        Eq(W, G*m*M/r + q**2/r + hbar*omega + m*c**2)
    assert convert_unit_system(eq, dims2, planck_units, natural) == \
        Eq(W, G*m*M/r + q**2/(4*pi*r) + omega + m)
    assert convert_unit_system(eq, dims2, planck_units, stoney_units) == \
        Eq(W, m*M/r + q**2/r + hbar*omega + m)
    assert convert_unit_system(eq, dims2, planck_units, hartree_atomic_units) == \
        Eq(W, G*m*M/r + q**2/r + omega + m*c**2)
    assert convert_unit_system(eq, dims2, planck_units, rydberg_atomic_units) == \
        Eq(W, G*m*M/r + q**2/r + omega + m*c**2)

    eq = Eq(W, -1/(2*symbols("n")**2))
    assert convert_unit_system(eq, {**dims, symbols("n"): 1}, hartree_atomic_units, rydberg_atomic_units) == \
        Eq(W, -1/symbols("n")**2)

    raises(NotImplementedError, lambda: convert_unit_system(
        eq, {**dims, symbols("n"): 1}, hartree_atomic_units, rydberg_atomic_units,
        constants=[hbar]))


def test_missing_dimensions():
    f = Function("f")
    raises(ValueError, lambda: convert_unit_system(Eq(W, q*U), {W: energy}, SI, cgs_gauss))
    raises(ValueError, lambda: convert_unit_system(Eq(W, f(x)), dims, SI, cgs_gauss))
    raises(TypeError, lambda: convert_unit_system(Eq(W, q*U), {W: energy, q: charge, U: volt}, SI, cgs_gauss))
    n = symbols("n")
    assert convert_unit_system(Eq(W, n*q**2/r), {**dims, n: 1}, cgs_gauss, SI) == Eq(W, k_e*n*q**2/r)
