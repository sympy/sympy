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
    >>> from sympy.physics.unitsystems import dim_simplify
    >>> length, mass, time = mks["length"], mks["mass"], mks["time"]
    >>> acceleration = mks["acceleration"]
    >>> m, a, F = symbols("m a F")
    >>> newton = m*a - F
    >>> sol = solve(newton, F)[0]
    >>> force = dim_simplify(sol.subs({m: mass, a: acceleration}))
    >>> force
    {'length': 1, 'mass': 1, 'time': -2}
    >>> force == mks["force"]
    True
    >>> M, r, G = symbols("m r G")
    >>> grav_force = F - G * m * M / r**2
    >>> sol = solve(grav_force, G)[0]
    >>> const = dim_simplify(sol.subs({m: mass, M: mass, r: length, F: force}))
    >>> const
    {'length': 3, 'mass': -1, 'time': -2}
