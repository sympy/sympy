
.. _solve_output:

==========================
Understanding Solve Output
==========================

The output of the `solve` function can seem very unwieldy since it may appear to
arbitrarily return one of 6 different types of output (in addition to raising
errors). The reasons for this are historical and are biased toward human
interaction rather than programmatic use. The type of output will depend on the
type of equation(s) (and how they are entered) and the number of symbols that
are provided (and how they are provided).

Since this output is controllable by setting `dict`, `set`, or `tuple` to True
(when not dealing with relational expressions) the following is an explanation for
those that might prefer the default output but wish to understand better the
conditions under which it is given.

    >>> from sympy import sqrt, exp, solve, Symbol, Eq
    >>> from sympy.abc import x, y, a, b

1. An empty list indicates that no solution was found.

    >>> solve(sqrt(x) + 1)
    []

2. A list of values is given when the symbol to solve for was
unambiguous in context because a) the equation was univariate or b) a
single symbol was specified as being of interest.

    >>> solve(x**2 - 4)
    [-2, 2]
    >>> solve(x - y - 1, x)
    [y + 1]

3. A single dictionary with keys being symbols and values being the solutions
for those symbols is the result when one or more equations (passed as a
*list*) has only a single, unique solution (e.g. the system is linear or else is
monotonic) for the symbols specified. Such a system is automatically generated
when ``match=True``.

    >>> solve([x + y - 2, x - y + 2], x, y)
    {x: 0, y: 2}
    >>> solve([exp(x) - y], x)
    {x: log(y)}
    >>> solve(a*x - 2*x + b - 5, {a, b}, match=True)
    {a: 2, b: 5}

4. A list of tuples, each giving a set of values for the symbols in the order
they were given, is given when more than one symbol was given in a
well defined order and either an expression was passed or a list of
expressions with at least one being linear/monotonic. (This is also the
format selected with ``tuple=True``.)

    >>> solve(x - 1, x, y)
    [(1, y)]
    >>> solve([x**2], x)
    [(0,)]
    >>> solve([x**2 - 1], x)
    [(-1,), (1,)]
    >>> solve([x**2 - 1], x, y)
    [(-1, y), (1, y)]

5. A list of dictionaries is returned when the expression was not univariate or
there was a nonlinear/nonmonotonic expression in a list *and* the order of
symbols would otherwise be ambiguous because a) no symbols were passed or b) the
symbols were passed as a set. The dictionary will only contain values that are
distinct from the keys. (This is also the format selected with ``dict=True``.)

    >>> solve(x - y)
    [{x: y}]
    >>> solve([exp(x) - 1, x**2])
    [{x: 0}]
    >>> solve([x + y - 2, x**2 - y + 2], {x, y})
    [{x: -1, y: 3}, {x: 0, y: 2}]

6. A boolean expression is returned when a relational expression other
than an Equality is given as an expression to solve. A single Equality
or a more complicated relational expression might be returned.

    >>> from sympy import Ne
    >>> solve([x**2 - 4, Ne(x, -2)])
    Eq(x, 2)
    >>> solve([x**2 > 4, x > 0])
    (2 < x) & (x < oo)
    >>> b = Symbol('b', positive=True)
    >>> solve([x**2 > b, x > 0], x)
    (0 < x) & (x < oo) & ((sqrt(b) < x) | (x < -sqrt(b)))

    The ``dict=True`` setting is ignored when working with relational
    expressions. When an equality is returned, it can be converted to
    a dictionary as follows:

    >>> soln = Eq(x, 2)
    >>> dict([soln.args]) == {soln.lhs: soln.rhs} == {x: 2}
    True

Note: using relational expressions to filter the solution only works when
solving for a single variable. Simple constraints to limit output to integers,
real numbers, positive values, etc... can be accomplished by setting the
assumptions on the variables for which a solution is being sought. This also
prevents the solution from being returned as a boolean expression.

    >>> p = Symbol('p', positive=True)
    >>> solve(p**2  - 4)
    [2]

For those that used a version of SymPy older than 1.12, `solve` automatically
tried to detect when a single equation might provide a solution for several
variables by matching coefficients on expressions. The output would be a single
dictionary (if the symbols appeared in a linear fashion) or else a list of
tuples if the system could be solved. Otherwise, an algebraic solution was sought
for one or more of the symbols. Newer versions of SymPy require permission for
this coefficient extraction via keyword ``match=True``, otherwise the algebraic
solution is returned. Anyone using an older version of SymPy can avoid the
undetermined coefficients solution by passing the equation as a *list*.
