Basic Usage
===========

Here are some basic examples of using the ``solve`` function:

.. code-block:: python

    >>> from sympy import solve, Symbol, Eq
    >>> x = Symbol('x')
    >>> solve(x**2 - 1, x)  # Solve x^2 = 1
    [-1, 1]
    >>> solve(x**3 - x, x)  # Solve x^3 = x
    [-1, 0, 1]
    >>> solve(Eq(x**2, 1), x)  # Can also use Eq objects
    [-1, 1]

Advanced Usage
=============

Here are some more advanced examples:

Systems of Equations
====================

.. code-block:: python

    >>> from sympy import symbols
    >>> x, y = symbols('x y')
    >>> solve([x + y - 2, x - y], [x, y])  # Solve system of linear equations
    {x: 1, y: 1}
    >>> solve([x**2 + y**2 - 1, y - x**2], [x, y])  # Solve system of nonlinear equations
    [(sqrt(2)/2, 1/2), (-sqrt(2)/2, 1/2)]

Equations with Parameters
=========================

.. code-block:: python

    >>> from sympy import Symbol, solve
    >>> a = Symbol('a')
    >>> x = Symbol('x')
    >>> solve(x**2 + a*x + 1, x)  # Solve quadratic equation with parameter
    [-a/2 - sqrt(-4 + a**2)/2, -a/2 + sqrt(-4 + a**2)/2]

Transcendental Equations
========================

.. code-block:: python

    >>> from sympy import exp, log
    >>> solve(exp(x) - 1, x)  # Solve exponential equation
    [0]
    >>> solve(log(x) - 1, x)  # Solve logarithmic equation
    [E]

Trigonometric Equations
=======================

.. code-block:: python

    >>> from sympy import sin, cos
    >>> solve(sin(x), x)  # Find zeros of sine function
    [0, pi]
    >>> solve(cos(x) - 1, x)  # Solve cosine equation
    [0]

Complex Solutions
================

.. code-block:: python

    >>> from sympy import I
    >>> solve(x**2 + 1, x)  # Solve equation with complex solutions
    [-I, I]
    >>> solve(x**2 + 2*x + 2, x)  # Solve quadratic with complex roots
    [-1 - I, -1 + I] 