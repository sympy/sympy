Operations on holonomic functions
=================================

Addition and Multiplication
---------------------------

Two holonomic functions can be added or multiplied with the result also
a holonomic functions.

    >>> from sympy.holonomic.holonomic import HolonomicFunction, DifferentialOperators
    >>> from sympy.polys.domains import QQ
    >>> from sympy import symbols
    >>> x = symbols('x')
    >>> R, Dx = DifferentialOperators(QQ.old_poly_ring(x),'Dx')
    >>> p = HolonomicFunction(Dx - 1, x, 0, [1])
    >>> q = HolonomicFunction(Dx**2 + 1, x, 0, [0, 1])
    >>> p + q
    HolonomicFunction((-1) + (1)*Dx + (-1)*Dx**2 + (1)*Dx**3, x, 0, [1, 2, 1])
    >>> p * q
    HolonomicFunction((2) + (-2)*Dx + (1)*Dx**2, x, 0, [0, 1])

``p`` and ``q`` here are holonomic representation of :math:`e^x` and
:math:`sin(x)` respectively.

Integration
-----------

Integrates the given holonomic function.

    >>> from sympy.holonomic.holonomic import HolonomicFunction, DifferentialOperators
    >>> from sympy.polys.domains import ZZ, QQ
    >>> from sympy import symbols
    >>> x = symbols('x')
    >>> R, Dx = DifferentialOperators(QQ.old_poly_ring(x),'Dx')
    >>> HolonomicFunction(Dx - 1, x, 0, [1]).integrate((x, 0, x))  # e^x - 1
    HolonomicFunction((-1)*Dx + (1)*Dx**2, x, 0, [0, 1])
    # integrate(cos(x), (x, 0, x)) = sin(x)
    >>> HolonomicFunction(Dx**2 + 1, x, 0, [1, 0]).integrate((x, 0, x))
    HolonomicFunction((1)*Dx + (1)*Dx**3, x, 0, [0, 1, 0])
