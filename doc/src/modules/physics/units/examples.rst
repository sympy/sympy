========
Examples
========

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
    >>> from sympy.physics.units import length, mass, acceleration, force
    >>> from sympy.physics.units import gravitational_constant as G
    >>> from sympy.physics.units.dimensions import dimsys_SI
    >>> F = mass*acceleration
    >>> F
    Dimension(acceleration*mass)
    >>> dimsys_SI.get_dimensional_dependencies(F)
    {'length': 1, 'mass': 1, 'time': -2}
    >>> dimsys_SI.get_dimensional_dependencies(force)
    {'length': 1, 'mass': 1, 'time': -2}

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

    >>> a.set_dimension(length, "SI")
    >>> a.set_scale_factor(108208000e3*meter, "SI")

    Add the solar mass as quantity:

    >>> M = Quantity("solar_mass")
    >>> M.set_dimension(mass, "SI")
    >>> M.set_scale_factor(1.9891e30*kilogram, "SI")

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
    71.5123904642338*pi*day
    >>> convert_to(q, day).n()
    224.662800523082*day

We could also have the solar mass and the day as units coming from the
astrophysical system, but we wanted to show how to create a unit that one needs.

We can see in this example that intermediate dimensions can be ill-defined,
such as sqrt(G), but one should check that the final result - when all
dimensions are combined - is well defined.
