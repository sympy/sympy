Solveset
========

.. module:: sympy.solvers.solveset

This is the official documentation of the `solveset` module in solvers.
It contains the Frequently asked Questions about our new module to solve
equations.

What's wrong with solve():
--------------------------

SymPy already has a pretty powerful `solve` function. But it has a lot of major
issues

1. It doesn't have a consistent output for various types of solutions
   It needs to return a lot of types of solutions consistently:

   * single solution : ` x == 1`
   * Multiple solutions: `x^2 == 1`
   * No Solution: `x^2 + 1 == 0` ; x is real
   * Interval of solution: `floor(x) == 0`
   * Infinitely many solutions: `sin(x) == 0`
   * Multivariate functions with point solutions: `x^2 + y^2 == 0`
   * Multivariate functions with non point solution: `x^2 + y^2 == 1`
   * System of equations: `x + y == 1` and `x - y == 0`
   * Relational: `x > 0`
   * And the most important case "We don't Know"

2. The input API is also a mess, there are a lot of parameters. Many of them
   are not needed and they make it hard for the user and the developers to
   work on solvers.

3. There are cases like finding the maxima and minima of function using
   critical points where it is important to know if it has returned all the
   solutions. `solve` does not guarantee this.


Why Solveset?
-------------

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


Why do we use Sets as an output type?
-------------------------------------

The main reason for using sets as output to solvers is that it can consistently
represent many types of solutions. For single variable case it can represent:

 * No solution (by null set)

 * Finitely many solutions (by FiniteSet)

 * Infinitely many solutions, both countably and uncountably infinite solutions
   (using the ImageSet module) Interval

 * There can also be bizarre solutions to equations like set of all rational number.

No other programmer's object (list, dictionary, python sets) provides the
flexibility of mathematical sets which our sets module try to emulate. The
second reason to use sets is that they are close to the entities which
mathematician's deals with and it makes it easier to reason about them.
Another advantage of using objects closer to mathematical entities is that the
user won't have to "learn" our representation and she can have her expectations
transferred from her mathematical experience.

For multivariate case we are representing solutions as a set of points in a n
dimensional space and a point is represented by a FiniteSet of ordered tuple.
Please Note that, the General FiniteSet is unordered, but a FiniteSet with
a tuple as it's first argument is ordered. For example:

 - FiniteSet(1, 2, 3) : {1, 2, 3}    # Unordered
 - FiniteSet((1, 2, 3)) {(1, 2, 3)}  # Ordered


Why not use Dictionary as Output?

  Dictionary are easy to deal with programatically but mathematically they are not
  very precise and use of them can quickly lead to inconsistency and a lot of
  confusion. For example addition of two solutions is a fairly valid operations
  but dictionaries by themselves have no add operation. Since we are representing
  solutions as dictionary we can define a custom add operation where value of x
  variable is added with value of other x variable and value of y variable
  with the value of other y variable. But what if the user tries to add solution
  of two equations which have different sets of variables? Say one has x and y
  as variables and other has u and v? We cannot allow this addition because
  the keys of dictionaries are unordered but this behavior will be inconsistent
  with the one variable case because it will be absurd to not to allow addition of
  two numbers, say one obtained by solution of x = 1 and other by solving y = 2.


What is this domain argument about?
-----------------------------------

 Solveset is designed to be independent of the assumptions on the
 variable being solved for and instead, uses the `domain` argument to
 decide the solver to dispatch the equation to, namely `solveset\_real`
 or `solveset\_complex`. It's unlike the old `solve` which considers the
 assumption on the variable.

    >>> from sympy.solvers.solveset import solveset
    >>> from sympy.abc import x
    >>> from sympy import S
    >>> solveset(x**2 + 1, x) # domain=S.Complexes is default
    {-I, I}
    >>> solveset(x**2 + 1, x, domain=S.Reals)
    EmptySet()


What are the general methods employed by solveset to solve an equation?
-----------------------------------------------------------------------

 Solveset uses various methods to solve an equation, here is a brief overview
 of the methodology:

 * The `domain` argument is first considered to know the domain in which
   the user is interested to get the solution.

 * If the given function is a relational (`>=`, `<=`, `>`, `<`), and the
   domain is real, then `solve\_univariate\_inequality` and solutions are
   returned. Solving for complex inequalities are not supported yet.

 * Based on the `domain`, the equation is dispatched to one of the two
   functions `solveset\_real` or `solveset\_complex`, which solves the
   given equation in complex and real domain respectively.

 * If the given function (equation) is a product of two or more functions,
   like say `f = g*h`, then the solution to the given equation is the Union
   of the solution of the equations `g = 0` and `h = 0`, if and only if both
   g and h are finite for a finite input. So, the solution is build up
   recursively.

 * The function class is now checked if it's Trigonometric or Hyperbolic, then
   the function `\_solve\_real\_trig` is called.

 * The function is now checked if there is any instance of Piecewise
   expression, if it is, then it's converted to explict expression and
   set pairs and then solved recursively.

 * The respective solvers now tries to invert the equation using the routines
   `invert\_real` and `invert\_complex`.

 * After the invert, the equations are checked for radical or Abs (Modulus),
   then the methods `\_solve\_radical` and `\_solve\_abs` respectively are
   used.

 * If none of the above method is successful, then method of polynomial is
   used as follows:

   - `\_solve\_as\_rational` is called, it's third argument is the
     `solveset\_solver` which can either be `solveset_real` or
     `solveset\_complex` based on these respective poly solvers
     `\_solve\_as\_poly\_real` and `\_solve\_as\_poly\_complex` is called.

   - The underlying method `\_solve\_as\_poly` solves the equation using
     polynomial techniques if it already is a polynomial equation or, with
     a change of variables, can be made so.

 * The Final solution set obtained is taken intersection with the input
   domain, and the resultant solution is returned.


What will you do with the old solve?
------------------------------------

 The `solve` would possibly be deprecated in future versions & we encourage our
 users to use `solveset`. We may proceeds as follows:

 * Replace all internal instances of solve by solveset by next release.
 * Raise a deprecation warning with solve calls possibly from next to
   next release.
 * Possibly remove the `solve` code in the release after the above release.
 * The issues pertaining to old `solve` would be addressed by new issues
   for `solveset`.

References
----------

 * https://github.com/sympy/sympy/wiki/GSoC-2015-Ideas#solvers
 * https://github.com/sympy/sympy/wiki/GSoC-2014-Application-Harsh-Gupta:-Solvers
 * https://github.com/sympy/sympy/wiki/GSoC-2015-Application-AMiT-Kumar--Solvers-:-Extending-Solveset
 * https://github.com/sympy/sympy/pull/9438#issuecomment-109289855


Solveset Module Reference
-------------------------

Use :func:`solveset` to solve equations or expressions (assumed to be equal to 0) for a single variable.
Solving an equation like `x^2 == 1` can be done as follows::

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