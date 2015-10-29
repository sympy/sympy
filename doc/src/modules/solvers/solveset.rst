Solveset
========

.. module:: sympy.solvers.solveset

This is the official documentation of the ``solveset`` module in solvers.
It contains the frequently asked questions about our new module to solve
equations.

What's wrong with solve():
--------------------------

SymPy already has a pretty powerful ``solve`` function. But it has a lot of major
issues

1. It doesn't have a consistent output for various types of solutions
   It needs to return a lot of types of solutions consistently:

   * Single solution : `x = 1`
   * Multiple solutions: `x^2 = 1`
   * No Solution: `x^2 + 1 = 0 ; x \in \mathbb{R}`
   * Interval of solution: `\lfloor x \rfloor = 0`
   * Infinitely many solutions: `sin(x) = 0`
   * Multivariate functions with point solutions: `x^2 + y^2 = 0`
   * Multivariate functions with non-point solution: `x^2 + y^2 = 1`
   * System of equations: `x + y = 1` and `x - y = 0`
   * Relational: `x > 0`
   * And the most important case "We don't Know"

2. The input API is also a mess, there are a lot of parameters. Many of them
   are not needed and they make it hard for the user and the developers to
   work on solvers.

3. There are cases like finding the maxima and minima of function using
   critical points where it is important to know if it has returned all the
   solutions. ``solve`` does not guarantee this.


Why Solveset?
-------------

* ``solveset`` has a cleaner input and output interface: ``solveset`` returns
  a set object and a set object takes care of all types of output. For
  cases where it doesn't "know" all the solutions a ``ConditionSet`` with partial
  solution is returned. For input it only takes the equation, the variables
  to solve for and the optional argument ``domain`` in which the equations has to
  be solved.

* ``solveset`` can return infinitely many solutions. For example solving for
  `\sin{(x)} = 0` returns `\{2 n \pi | n \in \mathbb{Z}\} \cup \{2 n \pi + \pi | n \in \mathbb{Z}\}`,
  whereas ``solve`` only returns `[0, \pi]`.

* There is a clear code level and interface level separation between solvers
  for equations in the complex domain and the real domain. For example
  solving `exp(x) = 1` when `x` is to be solved in the complex domain, returns
  the set of all solutions, that is `\{2 n i \pi | n \in \mathbb{Z}\}`, whereas
  if `x` is to be solved in the real domain then only `\{0\}` is returned.


Why do we use Sets as an output type?
-------------------------------------

SymPy has a well developed sets module, which can represent most of the set
containers in Mathematics such as:


 * ``FiniteSet``

   Represents a finite set of discrete numbers.


 * ``Interval``

   Represents a real interval as a set.


 * ``ProductSet``

   Represents a Cartesian product of sets.


 * ``ImageSet``

   Represents the image of a set under a mathematical function

    >>> from sympy import ImageSet, S, Lambda
    >>> from sympy.abc import x
    >>> squares = ImageSet(Lambda(x, x**2), S.Naturals)  # {x**2 for x in N}
    >>> 4 in squares
    True

 * ``ComplexRegion``

   Represents the set of all complex numbers in a region in the Argand plane.


 * ``ConditionSet``

   Represents the set of elements, which satisfies a given condition.


Also, the predefined set classes such as:

 * ``Naturals`` `\mathbb{N}`

   Represents the natural numbers (or counting numbers), which are all
   positive integers starting from 1.


 * ``Naturals0`` `\mathbb{N_0}`

   Represents the whole numbers, which are all the non-negative integers,
   inclusive of 0.


 * ``Integers`` `\mathbb{Z}`

   Represents all integers: positive, negative and zero.


 * ``Reals`` `\mathbb{R}`

   Represents the set of all real numbers.


 * ``Complexes`` `\mathbb{C}`

   Represents the set of all complex numbers.


 * ``EmptySet`` `\phi`

   Represents the empty set.

 The above six sets are available as Singletons, like say ``S.Integers``.


It is capable of most of the set operations in mathematics:

 * ``Union``
 * ``Intersection``
 * ``Complement``
 * ``SymmetricDifference``

The main reason for using sets as output to solvers is that it can consistently
represent many types of solutions. For the single variable case it can represent:

 * No solution (by the empty set).

 * Finitely many solutions (by ``FiniteSet``).

 * Infinitely many solutions, both countably and uncountably infinite solutions
   (using the ``ImageSet`` module) .

 * ``Interval``

 * There can also be bizarre solutions to equations like set of rational
   numbers.

No other programmer's object (list, dictionary, generator, python sets)
provides the flexibility of mathematical sets which our sets module try to
emulate. The second reason to use sets is that they are close to the entities
which mathematician's deals with and it makes it easier to reason about them.
Set objects conform to Pythonic conventions when possible, i.e., ``x in A`` and
``for i in A`` both work when they can be computed. Another advantage of using
objects closer to mathematical entities is that the user won't have to "learn"
our representation and she can have her expectations transferred from her
mathematical experience.

For the multivariate case we represent solutions as a set of points in a
n-dimensional space and a point is represented by a FiniteSet of ordered
tuples, which is a point in `\mathbb{R}^n` or \mathbb{C}^n\.

Please note that, the general ``FiniteSet`` is unordered, but a ``FiniteSet``
with a tuple as it's only argument becomes ordered, Since a tuple is ordered.
So the order in the tuple is mapped to a pre-defined order of variables,
while returning solutions.

For example:

 >>> from sympy import FiniteSet
 >>> FiniteSet(1, 2, 3)   # Unordered
 {1, 2, 3}
 >>> FiniteSet((1, 2, 3))  # Ordered
 {(1, 2, 3)}


Why not use dicts as output?

  Dictionary are easy to deal with programatically but mathematically they are
  not very precise and use of them can quickly lead to inconsistency and a lot
  of confusion. For example:

  * There are a lot of cases where we don't know the complete solution and we
    may like to output a partial solution, consider the equation `fg = 0`. The
    solution of this equation is the union of the solution of the following
    two equations: `f = 0`, `g = 0`. Let's say that we are able to solve
    `f = 0` but solving `g = 0` isn't supported yet. In this case we cannot
    represent partial solution of the given equation `fg = 0` using dicts.
    This problem is solved with sets using ``ConditionSet`` object:

    `sol_f` U `\{x | x ∊ \mathbb{R} ∧ g = 0\}`, where `sol_f` is the solution
    of the equation `f = 0`.

  * Using dict may lead to surprising results like:

    - ``solve(Eq(x**2, 1), x) != solve(Eq(y**2, 1), y)``
      Though, mathematically, it doesn't make sense. Using ``FiniteSet`` here
      solves the problem.

  * It can also not represent solutions for Equations like `|x| < 1`, which
    is a disk of radius 1 in the Argand Plane. This problem is solved using
    complex sets implemented as ``ComplexRegion``.


Input API of ``solveset``
-------------------------

``solveset`` has a cleaner input API, unlike ``solve``. It takes a maximum
of three arguments:

``solveset(equation, variable=None, domain=S.Complexes)``

* Equation(s)

  The equation(s) to solve.


* Variable(s)

  The variable(s) for which the equation is to be solved.


* Domain

  The domain in which the equation is to be solved.


 ``solveset`` removes the ``flags`` argument of ``solve``, which had made
 the input API messy and output API inconsistent.


What is this domain argument about?
-----------------------------------

 Solveset is designed to be independent of the assumptions on the
 variable being solved for and instead, uses the ``domain`` argument to
 decide the solver to dispatch the equation to, namely ``solveset_real``
 or ``solveset_complex``. It's unlike the old ``solve`` which considers the
 assumption on the variable.

    >>> from sympy import solveset, S
    >>> from sympy.abc import x
    >>> solveset(x**2 + 1, x) # domain=S.Complexes is default
    {-I, I}
    >>> solveset(x**2 + 1, x, domain=S.Reals)
    EmptySet()


What are the general methods employed by solveset to solve an equation?
-----------------------------------------------------------------------

 Solveset uses various methods to solve an equation, here is a brief overview
 of the methodology:

 * The ``domain`` argument is first considered to know the domain in which
   the user is interested to get the solution.


 * If the given function is a relational (``>=``, ``<=``, ``>``, ``<``), and
   the domain is real, then ``solve_univariate_inequality`` and solutions are
   returned. Solving for complex solutions of inequalities, like `x^2 < 0`
   is not yet supported.".


 * Based on the ``domain``, the equation is dispatched to one of the two
   functions ``solveset_real`` or ``solveset_complex``, which solves the
   given equation in the complex or real domain, respectively.


 * If the given expression is a product of two or more functions,
   like say `gh = 0`, then the solution to the given equation is the Union
   of the solution of the equations `g = 0` and `h = 0`, if and only if both
   `g` and `h` are finite for a finite input. So, the solution is built up
   recursively.


 * The function class is now checked if it's trigonometric or hyperbolic, then
   the function ``_solve_real_trig`` is called, which solves it by converting
   it to complex exponential form.


 * The function is now checked if there is any instance of a ``Piecewise``
   expression, if it is, then it's converted to explict expression and
   set pairs and then solved recursively.


 * The respective solver now tries to invert the equation using the routines
   ``invert_real`` and ``invert_complex``. These routines are based on the
   concept of mathematical inverse (though not exactly). It reduces the
   real/complex valued equation ``f(x) = y`` to a set of equations:
   ``{g(x)  = h_1(y), g(x) = h_2(y), ..., g(x) = h_n(y) }`` where ``g(x)`` is a
   simpler function than ``f(x)``. There is some work needed to be done in
   this to find invert of more complex expressions.


 * After the invert, the equations are checked for radical or Abs (Modulus),
   then the method ``_solve_radical`` tries to simplify the radical, by
   removing it using techniques like squarring, cubing etc, and ``_solve_abs``
   solves nested Modulus by considering the positive and negative variants,
   iteratively.


 * If none of the above method is successful, then methods of polynomial is
   used as follows:

   - The method to solve the rational function is called:
     ``_solve_as_rational(f, symbol, solveset_solver, as_poly_solver)``
     its third argument is the ``solveset_solver`` which can either be
     ``solveset_real`` or ``solveset_complex`` based on these, respective
     poly solvers ``_solve_as_poly_real`` or ``_solve_as_poly_complex``
     is called to solve as polynomial.

   - The underlying method ``_solve_as_poly`` solves the equation using
     polynomial techniques if its already a polynomial equation or, with
     a change of variables, can be made so.


 * The final solution set obtained is taken intersection with the input
   domain, and the resulting solution is returned.

.. Remember to change the above part when the new solver is implemented.


How do we manipulate and return an infinite solution?
-----------------------------------------------------

 * In the real domain, we use our ``ImageSet`` class in the sets module to
   return infinite solutions. ``ImageSet`` is an image of a set under
   a mathematical function. For example, to represent the solution
   of the equation `\sin{(x)} = 0`, we can use the ``ImageSet`` as:


   >>> from sympy import ImageSet, Lambda, pi, S, Dummy, pprint
   >>> n = Dummy('n')
   >>> pprint(ImageSet(Lambda(n, 2*pi*n), S.Integers), use_unicode=True)
   {2⋅n⋅π | n ∊ ℤ}


   Where ``n`` is a dummy variable. It is basically the image of the
   set of integers under the function `2\pi n`.

 * In the complex domain, we use complex sets, which are implemented as the
   ``ComplexRegion`` class in the sets module, to represent infinite
   solution in the Argand plane. For example to represent the solution
   of the equation `|z| = 1`, which is a unit circle, we can use
   the ``ComplexRegion`` as:


   >>> from sympy import ComplexRegion, FiniteSet, Interval, pi, pprint
   >>> pprint(ComplexRegion(FiniteSet(1)*Interval(0, 2*pi), polar=True), use_unicode=True)
   {r⋅(ⅈ⋅sin(θ) + cos(θ)) | r, θ ∊ {1} × [0, 2⋅π)}


   Where the ``FiniteSet`` in the ``ProductSet`` is the range of the value
   of `r`, which is the radius of the circle and the ``Interval`` is the range
   of `\theta`, the angle from the `x` axis representing a unit circle in the
   Argand plane.

   Note: We also have non-polar form notation for representing solution
   in rectangular form. For example, to represent first two quadrants in the
   Argand plane, we can write the ``ComplexRegion`` as:


   >>> from sympy import ComplexRegion, Interval, pi, oo, pprint
   >>> pprint(ComplexRegion(Interval(-oo, oo)*Interval(0, oo)), use_unicode=True)
   {x + y⋅ⅈ | x, y ∊ (-∞, ∞) × [0, ∞)}


   where the Intervals are the range of `x` and `y` for the set of complex
   numbers `x + iy`.


How does ``solveset`` ensure that it is not returning any wrong solution?
--------------------------------------------------------------------------

 Solvers in a Computer Algebra System are based on heuristic algorithms,
 so it's usually very hard to ensure 100% percent correctness, in every
 possible case. However there are still a lot of cases where we can ensure
 correctness. Solveset tries to verify correctness wherever it can. For
 example:

 Consider the equation `|x| = n`. A naive method to solve this equation
 would return ``{-n, n}`` as its solution, which is not correct since
 ``{-n, n}`` can be its solution if and only if ``n`` is positive.
 Solveset returns this information as well to ensure correctness.

    >>> from sympy import symbols, S, pprint, solveset
    >>> x, n = symbols('x, n')
    >>> pprint(solveset(abs(x) - n, x, domain=S.Reals), use_unicode=True)
    ([0, ∞) ∩ {n}) ∪ ((-∞, 0] ∩ {-n})

 Though, there still a lot of work needs to be done in this regard.


Search based solver and step-by-step solution
---------------------------------------------

 Note: This is under Development.

 After the introduction of :py:class:`~sympy.sets.conditionset.ConditionSet`, the
 solving of equations can be seen as set transformations. Here is an abstract
 view of the things we can do to solve equations.

 * Apply various set transformations on the given set.
 * Define a metric of the usability of solutions, or a notion of some
   solutions being better than others.
 * Different transformations would be the nodes of a tree.
 * Suitable searching techniques could be applied to get the best solution.

 ``ConditionSet`` gives us the ability to represent unevaluated equations and
 inequalities in forms like `\{x|f(x)=0; x \in S\}` and `\{x|f(x)>0; x \in S\}`
 but a more powerful thing about ``ConditionSet`` is that it allows us to write
 the intermediate steps as set to set transformation. Some of the transformations
 are:

 * Composition: `\{x|f(g(x))=0;x \in S\} \Rightarrow \{x|g(x)=y; x \in S, y \in \{z|f(z)=0; z \in S\}\}`

 * Polynomial Solver: `\{x | P(x) = 0;x \in S\} \Rightarrow  \{x_1,x_2, ... ,x_n\} ∩ S`
                      `\text{ where } `x_i` `\text{ are roots of } P(x)`

 * Invert solver: `\{x|f(x)=0;x \in S\} \Rightarrow  \{g(0)| \text{ all g such that } f(g(x)) = x\}`

 * logcombine: `\{x| log(f(x)) + log(g(x));x \in S\}`
               `\Rightarrow  \{x| log(f(x).g(x)); x \in S\} \text{ if } f(x) > 0 \text{ and } g(x) > 0`
               `\Rightarrow  \{x| log(f(x)) + log(g(x));x \in S\} \text{ otherwise}`

 * product solve: `\{x|f(x)*g(x)=0; x \in S\}`
                  `\Rightarrow  \{x|f(x)=0; x \in S\} U \{x|g(x)=0; x \in S\}`
                  `\text{ given } f(x) \text{ and } g(x) \text{ are bounded.}`
                  `\Rightarrow  \{x|f(x)*g(x)=0; x \in S\}, \text{ otherwise}`

 Since the output type is same as the input type any composition of these
 transformations is also a valid transformation. And our aim is to find
 the right sequence of compositions (given the atoms) which transforms
 the given condition set to a set which is not a condition set i.e.,
 FiniteSet, Interval, Set of Integers and their Union, Intersection,
 Complement or ImageSet. We can assign a cost function to each set,
 such that, the more desirable that form of set is to us, the less the value
 of the cost function. This way our problem is now reduced to finding the path
 from the initial ConditionSet to the lowest valued set on a graph where
 the atomic transformations forms the edges.


How do we deal with cases where only some of the solutions are known?
---------------------------------------------------------------------

 Creating a universal equation solver, which can solve each and every
 equation we encounter in mathematics is an ideal case for solvers in
 a Computer Algebra System. We always have some cases, which are not
 solved, or solved with incomplete solutions, so it's very important
 to represent that situation. For this type of situation we use
 `ConditionSet` class in the sets module, which acts as an unevaluated
 solveset object.

 See `Richardson's theorem <https://en.wikipedia.org/wiki/Richardson%27s_theorem>`_.

 ``ConditionSet`` is basically a Set of elements which satisfy a given
 condition. For example, to represent the solutions of the equation in
 the real domain:

 .. math::  (x^2 - 4)*(\sin(x) + x)

 We can represent it as:

 `\{-2, 2\} ∪ \{x | x \in \mathbb{R} ∧ x + \sin(x) = 0\}`


What will you do with the old solve?
------------------------------------

 There are still a few things ``solveset`` can't do, which the old ``solve``
 can, such as solving non linear multivariate & LambertW type equations.
 Hence, it's not yet a perfect replacement for old ``solve``. The ultimate
 goal is to:

 * Replace ``solve`` with ``solveset``, by the time solveset is
   at least powerful as ``solve``, i.e. ``solveset`` does everything
   that ``solve`` can do currently, and

 * Eventually rename ``solveset`` to ``solve``.


How are symbolic parameters handled in solveset?
------------------------------------------------

 Solveset is in its initial phase of development as of now, so the
 symbolic parameters aren't handled well for all the cases, but some
 work has been done in this regard to depict our ideology towards
 symbolic parameters. As an example the solving of `|x| = n` for `x`
 where `n` is a symbolic parameter. Solveset returns the value of `x`
 considering the domain of the symbolic parameter `n` as well, i.e. :

 `([0, \infty) ∩ \{n\}) ∪ ((-\infty, 0] ∩ \{-n\})`.

 It simply means `n` is the solution only when it belongs to the ``Interval``
 `[0, \infty)` and `-n` is the solution only when `n` belongs to the ``Interval``
 `(- \infty, 0]`.

 There are various other cases as well which needs to be addressed, like
 say, solving of `2^x + (a - 2)` for `x` where `a` is a symbolic parameter.
 As of now, It returns the solution as an intersection with `\mathbb{R}`, which
 is trivial, as it doesn't reveal the domain of `a` in the solution.

 Recently, we have also implemented a function to find the domain of the
 expression in a FiniteSet (Intersection with the interval) in which it is
 not-empty. It is a useful addition for dealing with symbolic parameters.
 For example:

    >>> from sympy import Symbol, FiniteSet, Interval, not_empty_in, sqrt, oo
    >>> from sympy.abc import x
    >>> not_empty_in(FiniteSet(x/2).intersect(Interval(0, 1)), x)
    [0, 2]
    >>> not_empty_in(FiniteSet(x, x**2).intersect(Interval(1, 2)), x)
    [-sqrt(2), -1] U [1, 2]


References
----------

 .. [1] https://github.com/sympy/sympy/wiki/GSoC-2015-Ideas#solvers
 .. [2] https://github.com/sympy/sympy/wiki/GSoC-2014-Application-Harsh-Gupta:-Solvers
 .. [3] https://github.com/sympy/sympy/wiki/GSoC-2015-Application-AMiT-Kumar--Solvers-:-Extending-Solveset
 .. [5] http://iamit.in/blog/
 .. [6] https://github.com/sympy/sympy/pull/2948 : Action Plan for improving solvers.
 .. [7] https://github.com/sympy/sympy/issues/6659 : ``solve()`` is a giant mess
 .. [8] https://github.com/sympy/sympy/pull/7523 : ``solveset`` PR
 .. [9] https://groups.google.com/forum/#!topic/sympy/-SIbX0AFL3Q
 .. [10] https://github.com/sympy/sympy/pull/9696
 .. [11] https://en.wikipedia.org/wiki/Richardson%27s_theorem


Solveset Module Reference
-------------------------

Use :func:`solveset` to solve equations or expressions (assumed to be equal to 0) for a single variable.
Solving an equation like `x^2 == 1` can be done as follows::

    >>> from sympy import solveset
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
