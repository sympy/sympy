========
Examples
========

In the following sections we give few examples of what can be do with this
module.


Dimensional analysis
====================

We will start from the second Newton's law

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
    >>> from sympy.physics.unitsystems import length, mass, acceleration, force
    >>> from sympy.physics.unitsystems import gravitational_constant as G
    >>> F = mass*acceleration
    >>> F
    Dimension(acceleration*mass)
    >>> F.get_dimensional_dependencies()
    {'length': 1, 'mass': 1, 'time': -2}
    >>> force.get_dimensional_dependencies()
    {'length': 1, 'mass': 1, 'time': -2}
    >>> F == force
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
    >>> from sympy.physics.unitsystems import Unit, Quantity, length, mass
    >>> from sympy.physics.unitsystems.systems import mks
    >>> from sympy.physics.unitsystems import day, gravitational_constant as G
    >>> T = symbols("T")
    >>> a = Quantity("venus_a", length, 108208000e3)
    >>> M = Quantity("solar_mass", mass, 1.9891e30)
    >>> Eq(T**2 / a**3, 4*pi**2 / G / M)
    Eq(T**2/venus_a**3, 4*pi**2/(gravitational_constant*solar_mass))
    >>> q = solve(T**2 / a**3 - 4*pi**2 / G / M, T)[1]
    >>> q
    6.28318530717959*sqrt(venus_a**3/(gravitational_constant*solar_mass))
    >>> convert_to(q, day)
    1.59123003109442e-10*sqrt(1993480633430344427209462)*day
    >>> convert_to(q, day).n()
    224.667*day

We could also have the solar mass and the day as units coming from the
astrophysical system, but I wanted to show how to create a unit that one needs.

We can see in this example that intermediate dimensions can be ill-defined,
such as sqrt(G), but one should check that the final result - when all
dimensions are combined - is well defined.
