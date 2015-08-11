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

   * single solution : `x == 1`
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

* `solveset` has a cleaner input and output interface: `solveset` returns
  a set object and a set object takes care of all types of output. For
  cases where it doesn't "know" all the solutions a `ConditionSet` with partial
  solution is returned. For input it only takes the equation, the variables
  to solve for & the optional argument `domain` in which the equations has to
  be solved.

* `solveset` can return infinitely many solutions. For example solving for
  `sin(x) = 0` returns `\{2 n \pi \| n \in \mathbb{Z}\} \cup \{2 n \pi + \pi \| n \in \mathbb{Z}\}`, whereas `solve`
  only returns [0, π].

* There is a clear code level and interface level separation between solvers
  for equations in complex domain and equations in real domain. For example
  solving `exp(x) == 1` when `x` is to be solved in complex `domain`, returns
  the set of all solutions that is `\{2 n i \pi \| n \in \mathbb{Z}\}`, whereas if `x` is to be
  solved in the Real `domain` then only `\{0\}` is returned.


Why do we use Sets as an output type?
-------------------------------------

.. First give an overview of the Sets module of SymPy

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

For multivariate case we are representing solutions as a set of points in a
n-dimensional space and a point is represented by a FiniteSet of ordered tuple.
Please Note that, the General FiniteSet is unordered, but a FiniteSet with
a tuple as it's first argument is ordered. For example:

 - FiniteSet(1, 2, 3) : {1, 2, 3}    # Unordered
 - FiniteSet((1, 2, 3)) {(1, 2, 3)}  # Ordered


Why not use Dictionary as Output?

  Dictionary are easy to deal with programatically but mathematically they are
  not very precise and use of them can quickly lead to inconsistency and a lot
  of confusion. For example addition of two solutions is a fairly valid
  operation but dictionaries by themselves have no add operation. Since we are
  representing solutions as dictionary we can define a custom add operation
  where value of `x` variable is added with value of other `x` variable and
  value of `y` variable with the value of other `y` variable. But what if the user
  tries to add solution of two equations which have different sets of variables?
  Say one has `x` and `y` as variables and other has `u` and `v`? We cannot allow this
  addition because the keys of dictionaries are unordered but this behavior will
  be inconsistent with the one variable case because it will be absurd to not to
  allow addition of two numbers, say one obtained by solution of `x = 1` and
  other by solving `y = 2`.




What is this domain argument about?
-----------------------------------

.. You first need to explain the input API, maybe you should restructure the
.. question so that it expects us to answer this question.

.. Also explain "why" we added the domain argument

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

.. Remember to change this part when the new solver is implemented.

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
   `g` and `h` are finite for a finite input. So, the solution is build up
   recursively.

 * The function class is now checked if it's Trigonometric or Hyperbolic, then
   the function `\_solve\_real\_trig` is called.
.. Tell what does "solve_real_trig" do, we don't expect the reader of this doc
.. to open up the code and see what it does.

 * The function is now checked if there is any instance of Piecewise
   expression, if it is, then it's converted to explict expression and
   set pairs and then solved recursively.

 * The respective solvers now tries to invert the equation using the routines
   `invert\_real` and `invert\_complex`.

 * After the invert, the equations are checked for radical or Abs (Modulus),
   then the methods `\_solve\_radical` and `\_solve\_abs` respectively are
   used.
.. Same here tell how does _solve_radical and _solve_abs work just don't give
.. the names

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

.. It is not necessary that you give all the steps involved you just need to
.. give a general overview


How do we manipulate and return an infinite solution?
-----------------------------------------------------
.. solution -> solutions probably

.. You should explain  ImageSet, Intergers and other set classes in the set
.. question above

 * In Real Domain, we use our `ImageSet` class in the sets module to
   return infinite solutions. `ImageSet` is an Image of a set under
   a mathematical function. For example, to represent the solution
   of the equation `sin(x) == 0`, we can use the ImageSet as:

   `ImageSet(Lambda(n, 2*n*pi), Integers())`

   Where n is a dummy variable. It is basically the image of the
   Integers set under the function `2*n*pi`.

 * In Complex Domain, we use Complex Sets, which is implemented as
   `ComplexPlane` class in the sets module, to represent infinite
   solution in the argand plane. For example to represent the solution
   of the equation `|z| == 1`, which is a unit circle, we can use
   the ComplexPlane as:

   `ComplexPlane(FiniteSet(1)*Interval(0, 2*pi), polar=True)`

   Where the FiniteSet in the `ProductSet` is the range of the value of `r`,
   which is the radius of the circle and the Interval is the range of theta,
   representing a unit circle in the argand plane.

   Note: We also have non-polar form notation for representing solution
   in rectangular form. For example, to represent first quadrant in argand
   plane, we can write the ComplexPlane as:

   `ComplexPlane(Interval(0, oo)*Interval(0, oo))`

   where the Intervals are the range of `x` and `y` for the set of complex
   number `(x + I*y)`.


How does solveset ensures that it is not returning any wrong solution?
----------------------------------------------------------------------

 Solvers in a Computer Algebra System are based on heuristics algorithms,
 so it's usually very hard to ensure cent percent correctness, in every
 possible case. However there are still a lot of cases where we can ensure
 correctness. Solveset tries to verify correctness wherever it can for
 example:

 Consider the equation `|x| = n`, a naive method to solve this equation
 would return {-n, n} as it's solution, which is not correct since {-n, n}
 can only be it's solution if and only if n is positive. Solveset returns
 this information as well to ensure correctness.

    >>> from sympy.solvers.solveset import solveset
    >>> from sympy import symbols, S
    >>> x, n = symbols('x, n')
    >>> solveset(abs(x) - n, x, domain=S.Reals)
    Intersection([0, oo), {n}) U Intersection((-oo, 0], {-n})

 Though, there still a lot of work needs to be done in this regard.

.. Also mention the search based solver we are trying to implement
.. and the step by step solutions.


How do we deal with cases where only some of the solutions are known?
---------------------------------------------------------------------

 Creating a Universal equation solver, which can solve each and every
 equation we encounter in mathematics is an ideal case for solvers in
 a Computer Algebra System. We always have some cases, which are not
 solved, or solved with incomplete solutions, so it's very important
 to represent that situation. For this type of situation we use
 `ConditionSet` class in the sets module, which acts as an unevaluated
 solveset object.

 ConditionSet is basically a Set of elements which satisfies a given
 condition. For example, to represent the solutions of the Equation
 in Real domain:

 .. math::  (x^2 - 4)*(sin(x) + x)

 We can represent it as:

 `FiniteSet(-2, 2) U ConditionSet(Lambda(x, Eq(sin(x) + x, 0)), S.Reals)`

.. Should we use the mathematical symbols here?


What will you do with the old solve?
------------------------------------

 The (current) `solve` would possibly be deprecated in future versions & we encourage
 our users to use `solveset`. We may proceeds as follows:

.. solveset will be renamed as solve

 * Replace all internal instances of solve by solveset by next release.
 * Raise a deprecation warning with solve calls possibly from next to
   next release.
 * Possibly rename `solve` to `solve\_old`, so that people can easily fix
   their code.
 * The issues pertaining to old `solve` would be addressed by new issues
   for `solveset`.


How are symbolic parameters handled in solveset?
------------------------------------------------

 Solveset is in it's initial phase of development as of now, so the
 symbolic parameters aren't handled well for all the cases, but some
 work has been done in this regard to depict our ideology towards
 symbolic parameters. As an instance the solving of `|x| = n` for `x`
 where `n` is a symbolic parameter.

 Solveset returns the value of `x` considering the domain of the symbolic
 parameter `n` as well, i.e. :

 `Intersection([0, oo), {n}) U Intersection((-oo, 0], {-n})`.

 It simply means `n` is the solution only when it belongs to the `Interval`
 `[0, oo)` and `-n` is the solution only when `n` belongs to the `Interval`
 `(-oo, 0]`.

 There are various other cases as well which needs to be addressed, like
 say, solving of `2**x + (a - 2)` for `x` where `a` is a symbolic parameter.
 As of now, It returns the solution as an intersection with `mathbb{R}`, which
 is trivial, as it doesn't reveals the domain of `a`, in the solution.

.. Also mention the no_empty_in PR by gxyd


References
----------

 * https://github.com/sympy/sympy/wiki/GSoC-2015-Ideas#solvers
 * https://github.com/sympy/sympy/wiki/GSoC-2014-Application-Harsh-Gupta:-Solvers
 * https://github.com/sympy/sympy/wiki/GSoC-2015-Application-AMiT-Kumar--Solvers-:-Extending-Solveset
 * https://github.com/sympy/sympy/pull/9438#issuecomment-109289855
 * http://iamit.in/blog/


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
