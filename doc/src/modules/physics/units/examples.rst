=============
More examples
=============

In the following sections we give few examples of what can be done with this
module.


Dimensional analysis
====================

We will start from Newton's second law

.. math::
    m a = F

where :math:`m, a` and :math:`F` are the mass, the acceleration and the force
respectively. Knowing the dimensions of :math:`m` (:math:`M`) and :math:`a`
(:math:`L T^{-2}`), we will determine the dimension of :math:`F`; obviously we
will find that it is a force: :math:`M L T^{-2}`.

From there we will use the expression of the gravitational force between the
particle of mass :math:`m` and the body of mass :math:`M`, at a distance
:math:`r`

.. math::
    F = \frac{G m M}{r^2}

to determine the dimension of the Newton's constant :math:`G`. The result
should be :math:`L^3 M^{-1} T^{-2}`.

    >>> from sympy import symbols
    >>> from sympy.physics.units.systems import SI
    >>> from sympy.physics.units import length, mass, acceleration, force
    >>> from sympy.physics.units import gravitational_constant as G
    >>> from sympy.physics.units.systems.si import dimsys_SI
    >>> F = mass*acceleration
    >>> F
    Dimension(acceleration*mass)
    >>> dimsys_SI.get_dimensional_dependencies(F)
    {Dimension(length): 1, Dimension(mass, M): 1, Dimension(time): -2}
    >>> dimsys_SI.get_dimensional_dependencies(force)
    {Dimension(length): 1, Dimension(mass): 1, Dimension(time): -2}

    Dimensions cannot compared directly, even if in the SI convention they are
    the same:

    >>> F == force
    False

    Dimension system objects provide a way to test the equivalence of
    dimensions:

    >>> dimsys_SI.equivalent_dims(F, force)
    True

    >>> m1, m2, r = symbols("m1 m2 r")
    >>> grav_eq = G * m1 * m2 / r**2
    >>> F2 = grav_eq.subs({m1: mass, m2: mass, r: length, G: G.dimension})
    >>> F2  #doctest: +SKIP
    Dimension(mass*length*time**-2)
    >>> F2.get_dimensional_dependencies()  #doctest: +SKIP
    {'length': 1, 'mass': 1, 'time': -2}

Note that one should first solve the equation, and then substitute with the
dimensions.


Equation with quantities
========================

Using Kepler's third law

.. math::
    \frac{T^2}{a^3} = \frac{4 \pi^2}{GM}

we can find the Venus orbital period using the known values for the other
variables (taken from Wikipedia). The result should be 224.701 days.

    >>> from sympy import solve, symbols, pi, Eq
    >>> from sympy.physics.units import Quantity, length, mass
    >>> from sympy.physics.units import day, gravitational_constant as G
    >>> from sympy.physics.units import meter, kilogram
    >>> T = symbols("T")
    >>> a = Quantity("venus_a")

    Specify the dimension and scale in SI units:

    >>> SI.set_quantity_dimension(a, length)
    >>> SI.set_quantity_scale_factor(a, 108208000e3*meter)

    Add the solar mass as quantity:

    >>> M = Quantity("solar_mass")
    >>> SI.set_quantity_dimension(M, mass)
    >>> SI.set_quantity_scale_factor(M, 1.9891e30*kilogram)

    Now Kepler's law:

    >>> eq = Eq(T**2 / a**3, 4*pi**2 / G / M)
    >>> eq
    Eq(T**2/venus_a**3, 4*pi**2/(gravitational_constant*solar_mass))
    >>> q = solve(eq, T)[1]
    >>> q
    2*pi*venus_a**(3/2)/(sqrt(gravitational_constant)*sqrt(solar_mass))

To convert to days, use the ``convert_to`` function (and possibly approximate
the outcoming result):

    >>> from sympy.physics.units import convert_to
    >>> convert_to(q, day)
    71.5112118495813*pi*day
    >>> convert_to(q, day).n()
    224.659097795948*day

We could also have the solar mass and the day as units coming from the
astrophysical system, but we wanted to show how to create a unit that one needs.

We can see in this example that intermediate dimensions can be ill-defined,
such as sqrt(G), but one should check that the final result - when all
dimensions are combined - is well defined.


Natural units
=============

In natural units some physical constants are equal to one. The unit systems
are defined in ``sympy.physics.units.systems``. In the natural units of
particle physics lengths and times are inverse energies:

    >>> from sympy.physics.units import convert_to, meter, second, electronvolt
    >>> from sympy.physics.units import speed_of_light, hbar
    >>> from sympy.physics.units.systems import natural
    >>> convert_to(speed_of_light, 1, natural)
    1
    >>> convert_to(hbar, 1, natural)
    1
    >>> convert_to(1/meter, electronvolt, natural).n(6)
    1.97327e-7*electronvolt
    >>> convert_to(1/second, electronvolt, natural).n(6)
    6.58212e-16*electronvolt

In Hartree atomic units all mechanical and electromagnetic quantities are
dimensionless, the speed of light is the inverse of the fine-structure
constant:

    >>> from sympy.physics.units.systems import hartree_atomic_units
    >>> convert_to(speed_of_light, 1, hartree_atomic_units).n(6)
    137.036

A unit system where other constants are equal to one is created with the
``contract`` method:

    >>> from sympy.physics.units import boltzmann_constant, kelvin
    >>> unit_system = natural.contract([boltzmann_constant])
    >>> convert_to(kelvin, electronvolt, unit_system).n(6)
    8.61733e-5*electronvolt

The equations can be converted between unit systems with
``convert_unit_system``, provided that the dimensions of the symbols are
specified:

    >>> from sympy import symbols, Eq
    >>> from sympy.physics.units import convert_unit_system, energy, mass, temperature
    >>> from sympy.physics.units.systems import SI
    >>> E, m, T = symbols("E m T")
    >>> dims = {E: energy, m: mass, T: temperature}
    >>> eq = Eq(E, m*speed_of_light**2 + boltzmann_constant*T)
    >>> convert_unit_system(eq, dims, SI, unit_system)
    Eq(E, T + m)
    >>> convert_unit_system(Eq(E, T + m), dims, unit_system, SI)
    Eq(E, boltzmann_constant*T + speed_of_light**2*m)
