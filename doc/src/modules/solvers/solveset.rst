Solveset
========

.. module:: sympy.solvers.solveset

Solveset Documentation
======================


What's wrong with solve():
==========================

SymPy already has a pretty powerful `solve` function. But it has a lot of major
issues

1. It doesn't have a consistent output for various types of solutions
   It needs to return a lot of types of solutions consistently:
   * single solution : ` x == 1`
   * Multiple solutions: `x**2 == 1`
   * No Solution: `x**2 + 1 == 0; x is real`
   * Interval of solution: `floor(x) == 0`
   * Infinitely many solutions: `sin(x) == 0`
   * Multivariate functions with point solutions: `x**2 + y**2 == 0`
   * Multivariate functions with non point solution: `x**2 + y**2 == 1`
   * System of equations: `x + y == 1` and `x - y == 0`
   * Relational: `x > 0`
   * And the most important case "We don't Know"


2. The input API is also a mess, there are a lot of parameters. Many of them
   are not needed and they make it hard for the user and the developers to
   work on solvers.

3. There are cases like finding the maxima and minima of function using
   critical points where it is important to know if it has returned all the
   solutions. `solve` does not guarantee this.

TODO

Why Solveset
============

* `solveset` has a cleaner input and output interface: `solveset` returns a set
  object and a set object take care of all the types of the output. For cases
  where it doesn't "know" all the solutions a `NotImplementedError` is raised.
  For input only takes the equation and the variables for which the equations
  has to be solved.

* `solveset` can return infinitely many solutions. For example solving for
  `sin(x) = 0` returns {2⋅n⋅π | n ∊ ℤ} ∪ {2⋅n⋅π + π | n ∊ ℤ} Whereas `solve`
  only returns [0, π]

* There is a clear code level and interface level separation between solvers
  for equations in complex domain and equations in real domain. For example
  solving `exp(x) = 1` when x is complex returns the set of all solutions that
  is {2⋅n⋅ⅈ⋅π | n ∊ ℤ} . Whereas if x is a real symbol then only {0} is
  returned.

* `solveset` returns a solution only when it can guarantee that it is returning
  all the solutions.

TODO

Design Decision
===============

* There is a code level and interface level separation between solvers
  for equations in complex domain and equations in real domain.
  - `solveset_real()`
  - `solveset_complex()`

* TODO

References
==========

 * https://github.com/sympy/sympy/wiki/GSoC-2015-Ideas#solvers
 * https://github.com/sympy/sympy/wiki/GSoC-2014-Application-Harsh-Gupta:-Solvers


TODO: Remove the line below
-----------------------------------------------------------------------------

Use :func:`solveset` to solve equations or expressions (assumed to be equal to 0) for a single variable.
Solving an equation like x**2 == 1 can be done as follows::

    >>> from sympy.solvers.solveset import *
    >>> from sympy import Symbol, Eq
    >>> x = Symbol('x')
    >>> solveset(Eq(x**2, 1), x)
    {-1, 1}

Or one may manually rewrite the equation as an expression equal to 0::

    >>> solveset(x**2 - 1, x)
    {-1, 1}

The first argument for :func:`solveset` is an expression (equal to zero) or an equation and the second argument
is the symbol that we want to solve the equation for.

.. autofunction:: sympy.solvers.solveset.solveset

.. autofunction:: sympy.solvers.solveset.solveset_real

.. autofunction:: sympy.solvers.solveset.solveset_complex

.. autofunction:: sympy.solvers.solveset.invert_real

.. autofunction:: sympy.solvers.solveset.invert_complex

.. autofunction:: sympy.solvers.solveset.domain_check


linear_eq_to_matrix
-------------------

.. autofunction:: sympy.solvers.solveset.linear_eq_to_matrix


linsolve
--------

.. autofunction:: sympy.solvers.solveset.linsolve


Diophantine Equations (DEs)
---------------------------

See :ref:`diophantine-docs`

Inequalities
------------

See :ref:`inequality-docs`

Ordinary Differential equations (ODEs)
--------------------------------------

See :ref:`ode-docs`.

Partial Differential Equations (PDEs)
-------------------------------------

See :ref:`pde-docs`.