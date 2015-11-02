=========
 Solvers
=========

    >>> from sympy import *
    >>> x, y, z = symbols('x y z')
    >>> init_printing(use_unicode=True)

A Note about Equations
======================

Recall from the :ref:`gotchas <tutorial_gotchas_equals>` section of this
tutorial that symbolic equations in SymPy are not represented by ``=`` or
``==``, but by ``Eq``.


    >>> Eq(x, y)
    x = y


However, there is an even easier way.  In SymPy, any expression not in an
``Eq`` is automatically assumed to equal 0 by the solving functions.  Since `a
= b` if and only if `a - b = 0`, this means that instead of using ``x == y``,
you can just use ``x - y``.  For example

    >>> solve(Eq(x**2, 1), x)
    [-1, 1]
    >>> solve(Eq(x**2 - 1, 0), x)
    [-1, 1]
    >>> solve(x**2 - 1, x)
    [-1, 1]

This is particularly useful if the equation you wish to solve is already equal
to 0.  Instead of typing ``solve(Eq(expr, 0), x)``, you can just use
``solve(expr, x)``.

Solving Equations Algebraically
===============================

The main function for solving algebraic equations, as we saw above, is
``solve``.  The syntax is ``solve(equations, variables)``, where, as we saw
above, ``equations`` may be in the form of ``Eq`` instances or expressions
that are assumed to be equal to zero.

.. TODO: This is a mess, because solve() has such a complicated interface.

When solving a single equation, the output of ``solve`` is a list of the
solutions.

    >>> solve(x**2 - x, x)
    [0, 1]

If no solutions are found, an empty list is returned, or
``NotImplementedError`` is raised.

    >>> solve(exp(x), x)
    []

.. note::

   If ``solve`` returns ``[]`` or raises ``NotImplementedError``, it doesn't
   mean that the equation has no solutions.  It just means that it couldn't
   find any.  Often this means that the solutions cannot be represented
   symbolically.  For example, the equation `x = \cos(x)` has a solution, but
   it cannot be represented symbolically using standard functions.

   >>> solve(x - cos(x), x)
   Traceback (most recent call last):
   ...
   NotImplementedError: multiple generators [x, exp(I*x)]
   No algorithms are implemented to solve equation exp(I*x)

   In fact, ``solve`` makes *no guarantees whatsoever* about the completeness
   of the solutions it finds.  Much of ``solve`` is heuristics, which may find
   some solutions to an equation or system of equations, but not all of them.

``solve`` can also solve systems of equations.  Pass a list of equations and a
list of variables to solve for.

    >>> solve([x - y + 2, x + y - 3], [x, y])
    {x: 1/2, y: 5/2}
    >>> solve([x*y - 7, x + y - 6], [x, y])
    [(-√2 + 3, √2 + 3), (√2 + 3, -√2 + 3)]

.. note::

   The type of the output of ``solve`` when solving systems of equations
   varies depending on the type of the input.  If you want a consistent
   interface, pass ``dict=True``.

   >>> solve([x - y + 2, x + y - 3], [x, y], dict=True)
   [{x: 1/2, y: 5/2}]
   >>> solve([x*y - 7, x + y - 6], [x, y], dict=True)
   [{x: -√2 + 3, y: √2 + 3}, {x: √2 + 3, y: -√2 + 3}]

.. _tutorial-roots:

``solve`` reports each solution only once.  To get the solutions of a
polynomial including multiplicity use ``roots``.

    >>> solve(x**3 - 6*x**2 + 9*x, x)
    [0, 3]
    >>> roots(x**3 - 6*x**2 + 9*x, x)
    {0: 1, 3: 2}

The output ``{0: 1, 3: 2}`` of ``roots`` means that ``0`` is a root of
multiplicity 1 and ``3`` is a root of multiplicity 2.

.. _tutorial-dsolve:

Solving Differential Equations
==============================

To solve differential equations, use ``dsolve``.  First, create an undefined
function by passing ``cls=Function`` to the ``symbols`` function.


    >>> f, g = symbols('f g', cls=Function)

``f`` and ``g`` are now undefined functions.  We can call ``f(x)``, and it
will represent an unknown function.

    >>> f(x)
    f(x)

Derivatives of ``f(x)`` are unevaluated.

    >>> f(x).diff(x)
    d
    ──(f(x))
    dx

(see the :ref:`Derivatives <tutorial-derivatives>` section for more on
derivatives).

To represent the differential equation `f''(x) - 2f'(x) + f(x) = \sin(x)`, we
would thus use

    >>> diffeq = Eq(f(x).diff(x, x) - 2*f(x).diff(x) + f(x), sin(x))
    >>> diffeq
                          2
             d           d
    f(x) - 2⋅──(f(x)) + ───(f(x)) = sin(x)
             dx           2
                        dx

To solve the ODE, pass it and the function to solve for to ``dsolve``.

    >>> dsolve(diffeq, f(x))
                        x   cos(x)
    f(x) = (C₁ + C₂⋅x)⋅ℯ  + ──────
                              2

``dsolve`` returns an instance of ``Eq``.  This is because in general,
solutions to differential equations cannot be solved explicitly for the
function.

    >>> dsolve(f(x).diff(x)*(1 - sin(f(x))), f(x))
    f(x) + cos(f(x)) = C₁

The arbitrary constants in the solutions from dsolve are symbols of the form
``C1``, ``C2``, ``C3``, and so on.
