sympy.solvers.polysys.solve_poly_system
=====================================

.. currentmodule:: sympy.solvers.polysys

.. autofunction:: solve_poly_system

Examples
--------

Basic Usage
~~~~~~~~~~

Solve a system of polynomial equations:

.. code-block:: pycon

    >>> from sympy import solve_poly_system
    >>> from sympy.abc import x, y
    >>> eq1 = x**2 + y**2 - 4
    >>> eq2 = x - y
    >>> solve_poly_system([eq1, eq2], x, y)
    [(-sqrt(2), -sqrt(2)), (sqrt(2), sqrt(2))]

Solve a system with multiple solutions:

.. code-block:: pycon

    >>> from sympy import solve_poly_system
    >>> from sympy.abc import x, y
    >>> eq1 = x**2 + y**2 - 1
    >>> eq2 = x**2 - y**2
    >>> solve_poly_system([eq1, eq2], x, y)
    [(-1/sqrt(2), -1/sqrt(2)), (-1/sqrt(2), 1/sqrt(2)), (1/sqrt(2), -1/sqrt(2)), (1/sqrt(2), 1/sqrt(2))]

Solve a system with parameters:

.. code-block:: pycon

    >>> from sympy import solve_poly_system
    >>> from sympy.abc import x, y, a
    >>> eq1 = x**2 + y**2 - a**2
    >>> eq2 = x - y
    >>> solve_poly_system([eq1, eq2], x, y)
    [(-a/sqrt(2), -a/sqrt(2)), (a/sqrt(2), a/sqrt(2))]

See Also
--------

- :func:`sympy.solvers.solvers.solve` - General algebraic equation solver
- :func:`sympy.solvers.solveset.solveset` - Set-based equation solver
- :func:`sympy.solvers.solvers.nsolve` - Numerical equation solver
- :func:`sympy.solvers.ode.dsolve` - Differential equation solver

References
----------

- :doc:`/tutorials/intro-tutorial/solvers` - Tutorial on using solvers
- :doc:`/modules/solvers/solvers` - Overview of the solvers module 