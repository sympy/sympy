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

    >>> from __future__ import division
    >>> from sympy import solve, Symbol, symbols
    >>> from sympy.physics.unitsystems.systems import mks_dim as mks
    >>> from sympy.physics.unitsystems.simplifiers import dim_simplify
    >>> length, mass, time = mks["length"], mks["mass"], mks["time"]
    >>> acceleration = mks["acceleration"]
    >>> m, a, F = symbols("m a F")
    >>> newton = m*a - F
    >>> sol = solve(newton, F)[0]
    >>> force = dim_simplify(sol.subs({m: mass, a: acceleration}))
    >>> force  #doctest: +SKIP
    {'length': 1, 'mass': 1, 'time': -2}
    >>> force == mks["force"]
    True
    >>> M, r, G = symbols("M r G")
    >>> grav_force = F - G * m * M / r**2
    >>> sol = solve(grav_force, G)[0]
    >>> const = dim_simplify(sol.subs({m: mass, M: mass, r: length, F: force}))
    >>> const  #doctest: +SKIP
    {'length': 3, 'mass': -1, 'time': -2}

Note that one should first solve the equation, and then substitute with the
dimensions.


Equation with quantities
========================

Using Kepler's third law

.. math::
    \frac{T^2}{a^3} = \frac{4 \pi^2}{GM}

we can find the Venus orbital period using the known values for the other
variables (taken from Wikipedia). The result should be 224.701 days.

    >>> from __future__ import division
    >>> from sympy import solve, Symbol, symbols, pi
    >>> from sympy.physics.unitsystems import Unit, Quantity as Q
    >>> from sympy.physics.unitsystems.simplifiers import qsimplify
    >>> from sympy.physics.unitsystems.systems import mks
    >>> m, kg, s = mks["m"], mks["kg"], mks["s"]
    >>> T, a, M, G = symbols("T a M G")
    >>> venus_a = Q(108208000e3, m)
    >>> solar_mass = Q(1.9891e30, kg)
    >>> venus_subs = {"a": venus_a, "M": solar_mass, "G": mks["G"]}
    >>> Tsol = solve(T**2 / a**3 - 4*pi**2 / G / M, T)[1]
    >>> q = qsimplify(Tsol.subs(venus_subs))
    >>> day = Unit(s.dim, abbrev="day", factor=86400)
    >>> print(q.convert_to(qsimplify(day)))
    224.667 day

We could also have the solar mass and the day as units coming from the
astrophysical system, but I wanted to show how to create a unit that one needs.

We can see in this example that intermediate dimensions can be ill-defined,
such as sqrt(G), but one should check that the final result - when all
dimensions are combined - is well defined.
