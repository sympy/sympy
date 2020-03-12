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

    ``p`` and ``q`` here are holonomic representation of `e^x` and
    `\sin(x)` respectively.

    >>> p = HolonomicFunction(Dx - 1, x, 0, [1])
    >>> q = HolonomicFunction(Dx**2 + 1, x, 0, [0, 1])

    Holonomic representation of `e^x+\sin(x)`

    >>> p + q
    HolonomicFunction((-1) + (1)*Dx + (-1)*Dx**2 + (1)*Dx**3, x, 0, [1, 2, 1])

    Holonomic representation of `e^x \cdot \sin(x)`

    >>> p * q
    HolonomicFunction((2) + (-2)*Dx + (1)*Dx**2, x, 0, [0, 1])

.. currentmodule:: sympy.holonomic.holonomic

Integration and Differentiation
-------------------------------

.. automethod:: HolonomicFunction.integrate

.. automethod:: HolonomicFunction.diff

Composition with polynomials
----------------------------

.. automethod:: HolonomicFunction.composition

Convert to holonomic sequence
-----------------------------

.. automethod:: HolonomicFunction.to_sequence

Series expansion
----------------

.. automethod:: HolonomicFunction.series

Numerical evaluation
--------------------

.. automethod:: HolonomicFunction.evalf

Convert to a linear combination of hypergeometric functions
-----------------------------------------------------------

.. automethod:: HolonomicFunction.to_hyper

Convert to a linear combination of Meijer G-functions
-----------------------------------------------------

.. automethod:: HolonomicFunction.to_meijerg

Convert to expressions
----------------------

.. automethod:: HolonomicFunction.to_expr
