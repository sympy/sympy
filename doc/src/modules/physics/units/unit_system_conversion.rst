===============================
Conversion between unit systems
===============================

The form of the equations of physics depends on the unit system. The function
``convert_unit_system`` converts expressions and equations from a unit system
to another one. It should not be confused with ``convert_to``, which expresses
a quantity with different units of the same unit system.

Usage
=====

The dimensions of all symbols and undefined functions appearing in the
expression have to be specified with a dictionary. Dimensionless symbols
have dimension ``1``. Physical constants and units are represented by
quantities, their dimensions are known to the unit systems.

    >>> from sympy import symbols, pi, Eq
    >>> from sympy.physics.units import convert_unit_system
    >>> from sympy.physics.units import force, charge, length
    >>> from sympy.physics.units import vacuum_permittivity
    >>> from sympy.physics.units.systems.si import SI
    >>> from sympy.physics.units.systems.cgs import cgs_gauss
    >>> F, q1, q2, r = symbols("F q1 q2 r")
    >>> dims = {F: force, q1: charge, q2: charge, r: length}
    >>> coulomb_law = Eq(F, q1*q2/(4*pi*vacuum_permittivity*r**2))
    >>> coulomb_law_gauss = convert_unit_system(coulomb_law, dims, SI, cgs_gauss)
    >>> coulomb_law_gauss
    Eq(F, q1*q2/r**2)

The unit systems may also be specified by their names:

    >>> convert_unit_system(coulomb_law, dims, "SI", "cgs_gauss")
    Eq(F, q1*q2/r**2)

The SI has the current among its base dimensions. In Gaussian units the
current is a derived dimension, because the Coulomb constant is equal to one.
Converting from Gaussian units to the SI, the Coulomb constant is restored:

    >>> convert_unit_system(coulomb_law_gauss, dims, cgs_gauss, SI)
    Eq(F, coulomb_constant*q1*q2/r**2)

The constants appearing in the result can be chosen with the parameter
``constants``:

    >>> from sympy.physics.units import speed_of_light
    >>> convert_unit_system(coulomb_law_gauss, dims, cgs_gauss, SI,
    ...     constants=[vacuum_permittivity, speed_of_light])
    Eq(F, q1*q2/(4*vacuum_permittivity*pi*r**2))

Fields are represented by undefined functions. The dimension may be given
either to the function or to the applied function. As an example, these are
Maxwell's equations in one space dimension:

    >>> from sympy import Function
    >>> from sympy.physics.units import voltage, magnetic_density, current
    >>> from sympy.physics.units import time, volume, magnetic_constant
    >>> x, t = symbols("x t")
    >>> E, B, J, rho = symbols("E B J rho", cls=Function)
    >>> dims = {x: length, t: time, E: voltage/length, B: magnetic_density,
    ...     J: current/length**2, rho: charge/volume}
    >>> equations = [
    ...     Eq(E(x, t).diff(x), rho(x, t)/vacuum_permittivity),
    ...     Eq(E(x, t).diff(x), -B(x, t).diff(t)),
    ...     Eq(B(x, t).diff(x), magnetic_constant*J(x, t)
    ...         + magnetic_constant*vacuum_permittivity*E(x, t).diff(t))]
    >>> for eq in convert_unit_system(equations, dims, SI, cgs_gauss):
    ...     print(eq)
    Eq(Derivative(E(x, t), x), 4*pi*rho(x, t))
    Eq(Derivative(E(x, t), x), -Derivative(B(x, t), t)/speed_of_light)
    Eq(Derivative(B(x, t), x), 4*pi*J(x, t)/speed_of_light + Derivative(E(x, t), t)/speed_of_light)

Expressions
-----------

An equation relates quantities of the same dimension, so that both sides are
converted in the same way. An expression which is not an equation needs
the dimension of the quantity it represents, unless the quantity is defined in
the same way in both unit systems. The magnetic field is an example of a
quantity with a different definition in the SI and in Gaussian units:

    >>> i = symbols("i")
    >>> dims = {i: current, r: length}
    >>> expr = magnetic_constant*i/(2*pi*r)
    >>> convert_unit_system(expr, dims, SI, cgs_gauss, dimension=magnetic_density)
    2*i/(speed_of_light*r)

Without the parameter ``dimension`` the symbols and the constants are
replaced, but the result is the magnetic field of the SI:

    >>> convert_unit_system(expr, dims, SI, cgs_gauss)
    2*i/(speed_of_light**2*r)

Check of the results
--------------------

The results can be checked replacing the symbols with quantities and
converting to the same unit with ``convert_to``:

    >>> from sympy.physics.units import convert_to, ampere, meter, tesla
    >>> values = {i: 3*ampere, r: 2*meter}
    >>> expr_gauss = convert_unit_system(expr, dims, SI, cgs_gauss,
    ...     dimension=magnetic_density)
    >>> convert_to(expr.subs(values), tesla, SI)
    3*tesla/10000000
    >>> convert_to(expr_gauss.subs(values), tesla, cgs_gauss)
    3*tesla/10000000

How it works
============

The unit system with more independent dimensions is called *full*, the other
one *reduced*. A quantity is represented by `X_f` in the full unit system and
by `X_r` in the reduced one, with

.. math::
    X_f = K X_r

The factor `K` is a product of powers of physical constants of the full unit
system. They are of two kinds:

- the constants that are pure numbers in the reduced unit system, like the
  Coulomb constant in Gaussian units;
- the constants with the same dimension in both unit systems, like the speed
  of light in Gaussian units.

The exponents are determined by the requirement that the dimension of `K` is
the ratio between the dimensions of `X_f` and `X_r`. For example, the
dimensions of the charge are

.. math::
    [q_f] = \mathsf{I} \mathsf{T}, \qquad
    [q_r] = \mathsf{M}^{1/2} \mathsf{L}^{3/2} \mathsf{T}^{-1}

while the dimension of the Coulomb constant in the SI is

.. math::
    [k_e] = \mathsf{M} \mathsf{L}^3 \mathsf{T}^{-4} \mathsf{I}^{-2}

The only possibility is

.. math::
    q_f = \frac{q_r}{\sqrt{k_e}}

The magnetic field requires the speed of light as well:

.. math::
    B_f = \frac{\sqrt{k_e}}{c} B_r

The conversion from the full to the reduced unit system replaces every
symbol with the symbol multiplied by `K`, then the constants that are pure
numbers in the reduced unit system are replaced by their values. The
constants of the full unit system whose value is fixed by the reduced unit
system are replaced as well: in Gaussian units the vacuum permittivity is
`1/(4 \pi)`, the magnetic constant is `4 \pi/c^2`.

The conversion from the reduced to the full unit system replaces every symbol
with the symbol divided by `K`.

The equations are divided by the factor of the left hand side. If the terms
of a sum have different factors, the simplest one is collected.

Limitations
===========

- The conversion is determined by the dimensions only. The definition of some
  quantities differs by a numerical factor, for example the fields
  `\mathbf{D}` and `\mathbf{H}` of Gaussian units contain a factor `4 \pi`
  with respect to the ones of the SI. These factors have to be added
  manually.
- The base dimensions of one of the unit systems have to be independent
  dimensions in the other one.
- There is no check of the dimensional consistency of the expression.

Reference
=========

.. automodule:: sympy.physics.units.unit_system_conversion

.. autofunction:: convert_unit_system
