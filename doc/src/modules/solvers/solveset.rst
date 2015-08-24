Solveset
========

.. module:: sympy.solvers.solveset

This is the official documentation of the ``solveset`` module in solvers.
It contains the Frequently asked Questions about our new module to solve
equations.

What's wrong with solve():
--------------------------

SymPy already has a pretty powerful ``solve`` function. But it has a lot of major
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
   solutions. ``solve`` does not guarantee this.


Why Solveset?
-------------

* ``solveset`` has a cleaner input and output interface: ``solveset`` returns
  a set object and a set object takes care of all types of output. For
  cases where it doesn't "know" all the solutions a `ConditionSet` with partial
  solution is returned. For input it only takes the equation, the variables
  to solve for & the optional argument `domain` in which the equations has to
  be solved.

* ``solveset`` can return infinitely many solutions. For example solving for
  `sin(x) = 0` returns `\{2 n \pi \| n \in \mathbb{Z}\} \cup \{2 n \pi + \pi \| n \in \mathbb{Z}\}`, whereas `solve`
  only returns [0, π].

* There is a clear code level and interface level separation between solvers
  for equations in complex domain and equations in real domain. For example
  solving `exp(x) == 1` when `x` is to be solved in complex `domain`, returns
  the set of all solutions that is `\{2 n i \pi \| n \in \mathbb{Z}\}`, whereas if `x` is to be
  solved in the Real `domain` then only `\{0\}` is returned.


Why do we use Sets as an output type?
-------------------------------------

SymPy has a well developed sets module, which can represent most of the set
containers in Mathematics such as:

 * ``FiniteSet``
   Represents a finite set of discrete numbers.

 * ``Interval`` 
   Represents a real interval as a Set.

 * ``ProductSet``
   Represents a Cartesian Product of Sets.

 * ``ImageSet``
   Represents the Image of a set under a mathematical function

    >>> from sympy import ImageSet, S, Lambda
    >>> from sympy.abc import x
    >>> from sympy.sets import ImageSet
    >>> squares = ImageSet(Lambda(x, x**2), S.Naturals)  # {x**2 for x in N}
    >>> 4 in squares
    True

 * ``ComplexPlane``
   Represents the Set of all Complex Numbers

 * ``ConditionSet``
   Represents the Set of elements which satisfies a given condition

Also, the predefined set classes such as:

 * ``Naturals`` `\mathbb{N}`
   Represents the natural numbers (or counting numbers) which are all
   positive integers starting from 1 

 * ``Naturals0`` `\mathbb{W}`
   Represents the whole numbers which are all the non-negative integers,
   inclusive of zero

 * ``Integers`` `\mathbb{Z}`
   Represents all integers: positive, negative and zero. This set is also
   available as the Singleton, ``S.Integers``.

 * ``Reals`` `\mathbb{R}`
   Represents the Set of all Real numbers.

 * ``Complexes`` `\mathbb{C}`
   Represents the Set of all Complex numbers.

It is capable of most of the set operations in Mathematics:

 * ``Union``
 * ``Intersection``
 * ``Complement``
 * ``SymmetricDifference``

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

 - ``FiniteSet(1, 2, 3) : {1, 2, 3}``   # Unordered
 - ``FiniteSet((1, 2, 3)) {(1, 2, 3)}``  # Ordered


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



Input API of `Solveset`
-----------------------

``solveset`` has a cleaner input API, unlike ``solve``, It takes a maximum
of three arguments:

* Equation(s)
  The equation(s) to solve.

* Variable(s)
  The variable(s) for which the equation is to be solved.

* domain
  The domain in which the equation is to be solved.


 ``solveset`` removes the ``flags`` argument of ``solve``, which had made
 the input API messy and output API inconsistent in ``solve``.


What is this domain argument about?
-----------------------------------

 Solveset is designed to be independent of the assumptions on the
 variable being solved for and instead, uses the ``domain`` argument to
 decide the solver to dispatch the equation to, namely ``solveset\_real``
 or ``solveset\_complex``. It's unlike the old ``solve`` which considers the
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

 * The ``domain`` argument is first considered to know the domain in which
   the user is interested to get the solution.

 * If the given function is a relational (`>=`, `<=`, `>`, `<`), and the
   domain is real, then `solve\_univariate\_inequality` and solutions are
   returned. Solving for complex inequalities are not supported yet.


 * Based on the ``domain``, the equation is dispatched to one of the two
   functions ``solveset\_real`` or ``solveset\_complex``, which solves the
   given equation in complex and real domain respectively.

 * If the given function (equation) is a product of two or more functions,
   like say `f = g*h`, then the solution to the given equation is the Union
   of the solution of the equations `g = 0` and `h = 0`, if and only if both
   `g` and `h` are finite for a finite input. So, the solution is build up
   recursively.

 * The function class is now checked if it's Trigonometric or Hyperbolic, then
   the function ``\_solve\_real\_trig`` is called, which solves it by converting
   it in terms of ``exp`` complex exponential form.

 * The function is now checked if there is any instance of Piecewise
   expression, if it is, then it's converted to explict expression and
   set pairs and then solved recursively.

 * The respective solvers now tries to invert the equation using the routines
   ``invert\_real`` and ``invert\_complex``.

 * After the invert, the equations are checked for radical or Abs (Modulus),
   then the method ``\_solve\_radical`` tries to simplify the radical, by
   removing it using techniques like squarring, cubing etc, and ``\_solve\_abs``
   solves nested Modulus by considering the positive and negative variants,
   iteratively.

 * If none of the above method is successful, then method of polynomial is
   used as follows:

   - ``\_solve\_as\_rational`` is called, it's third argument is the
     ``solveset\_solver`` which can either be ``solveset_real`` or
     ``solveset\_complex`` based on these respective poly solvers
     ``\_solve\_as\_poly\_real`` and ``\_solve\_as\_poly\_complex`` is called.

   - The underlying method ``\_solve\_as\_poly`` solves the equation using
     polynomial techniques if it already is a polynomial equation or, with
     a change of variables, can be made so.

 * The Final solution set obtained is taken intersection with the input
   domain, and the resultant solution is returned.

.. Remember to change the above part when the new solver is implemented.


How do we manipulate and return an infinite solution?
-----------------------------------------------------

 * In Real Domain, we use our ``ImageSet`` class in the sets module to
   return infinite solutions. ``ImageSet`` is an Image of a set under
   a mathematical function. For example, to represent the solution
   of the equation `sin(x) == 0`, we can use the ImageSet as:

   ``ImageSet(Lambda(n, 2*n*pi), Integers())``

   Where n is a dummy variable. It is basically the image of the
   Integers set under the function `2*n*pi`.

 * In Complex Domain, we use Complex Sets, which is implemented as
   ``ComplexPlane`` class in the sets module, to represent infinite
   solution in the argand plane. For example to represent the solution
   of the equation `|z| == 1`, which is a unit circle, we can use
   the ComplexPlane as:

   ``ComplexPlane(FiniteSet(1)*Interval(0, 2*pi), polar=True)``

   Where the ``FiniteSet`` in the ``ProductSet`` is the range of the value of `r`,
   which is the radius of the circle and the ``Interval`` is the range of theta,
   representing a unit circle in the argand plane.

   Note: We also have non-polar form notation for representing solution
   in rectangular form. For example, to represent first quadrant in argand
   plane, we can write the ComplexPlane as:

   ``ComplexPlane(Interval(0, oo)*Interval(0, oo))``

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


Search Based solver & Step by Step Solution
-------------------------------------------

 Note: This is under Development.

 After the introduction of ``ConditionSet`` [10], the solving of equations can be
 seen as set transformation, We can do the following things to solve equations
 (abstract View):

 * Apply Various Set Transformations on the given Set.
 * Define a Metric of the usability or define a notion of better solution over others.
 * Different Transformation would be the nodes of the tree.
 * Suitable searching techniques could be applied to get the best solution.

 ConditionSet gives us the ability to represent unevaluated equations and
 inequalities in forms like `{x|f(x)=0; x in S}` and `{x|f(x)>0; x in S}`
 but a more powerful thing about ConditionSet is that it allows us to write
 the intermediate steps as set to set transformation. Some of the transformations
 are:

 * Composition: `{x|f(g(x))=0;x in S} => {x|g(x)=y; x in S,y in {z|f(z)=0; z in S}}`

 * Polynomial Solver: `{x|P(x)=0;x in S} => {x_1,x_2, ... ,x_n}.intersection(S)`
   where `x_i` are roots of P(x)

 * Invert solver: `{x|f(x)=0;x in S} => {g(0)| all g such that f(g(x)) = x}`

 * logcombine: `{x| log(f(x)) + log(g(x));x in S}`
                => `{x| log(f(x)*g(x)); x in S}` if f(x) > 0 and g(x) > 0
                => `{x| log(f(x)) + log(g(x));x in S}` otherwise

 * product solve: {x|f(x)*g(x)=0; x in S}
                => `{x|f(x)=0; x in S} U {x|g(x)=0; x in S}` given f(x) and g(x) are bounded
                => {x|f(x)*g(x)=0; x in S}, otherwise

 Since the output type is same as input type any composition of these
 transformations is also a valid transformation. And our aim is to find
 the right sequence of compositions (given the atoms) which transforms
 the given condition set to a set which is not a condition set i.e.,
 FiniteSet, Interval, Set of Integers and their Union, Intersection,
 Complement or ImageSet. We can assign a cost function to each set,
 such that, more desirable that form of set is to us, the less the value
 of cost function. This way our problem is now reduced to finding the path
 from the initial ConditionSet to the lowest valued set on a graph where
 the atomic transformations forms the edges.


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

 `{-2, 2} ∪ {x | x ∊ ℝ ∧ x + sin(x) = 0}`


What will you do with the old solve?
------------------------------------

 There is still a few things `solveset` can't do, which the old `solve`
 can, such as solving non linear multivariate & LambertW type equations.
 Hence, it's not yet a perfect replacement for old `solve`. The ultimate
 goal is to:

 * Replace ``solve`` with ``solveset``, by the time solveset is
   atleast powerful as ``solve``, i.e. ``solveset`` does everything
   that ``solve`` can do currently, and

 * Eventually rename ``solveset`` as ``solve``. Meanwhile
   issue a deprecation warning for the current behavior of ``solve``.


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

 It simply means `n` is the solution only when it belongs to the ``Interval``
 `[0, oo)` and `-n` is the solution only when `n` belongs to the ``Interval``
 `(-oo, 0]`.

 There are various other cases as well which needs to be addressed, like
 say, solving of `2**x + (a - 2)` for `x` where `a` is a symbolic parameter.
 As of now, It returns the solution as an intersection with `\mathbb{R}`, which
 is trivial, as it doesn't reveals the domain of `a`, in the solution.

 Recently, we have also implemented a function to find the domain of the
 expression in a FiniteSet (Intersection with the interval) in which it is
 not-empty. It is a useful addition for dealing with symbolic parameters.
 For e.g :

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
 .. [4] https://github.com/sympy/sympy/pull/9438#issuecomment-109289855
 .. [5] http://iamit.in/blog/
 .. [6] https://github.com/sympy/sympy/pull/2948 : Action Plan for improving solvers.
 .. [7] https://github.com/sympy/sympy/issues/6659 : ``solve()`` is a giant mess
 .. [8] https://github.com/sympy/sympy/pull/7523 : ``solveset`` PR 
 .. [9] https://groups.google.com/forum/#!topic/sympy/-SIbX0AFL3Q
 .. [10] https://github.com/sympy/sympy/pull/9696


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
