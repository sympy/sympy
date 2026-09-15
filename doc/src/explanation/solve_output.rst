.. _solve_output:

====================
Solve Output by Type
====================

The output of the :func:`~.solve` function can seem very unwieldy since it may appear to
arbitrarily return one of six different types of output (in addition to raising
errors). The reasons for this are historical and are biased toward human
interaction rather than programmatic use. The type of output will depend on the
type of equation(s) (and how they are entered) and the number of symbols that
are provided (and how they are provided).

    >>> from sympy import sqrt, exp, solve, Symbol, Eq
    >>> from sympy.abc import x, y, z, a, b

    The :func:`~.solve` function attempts to find all values for as many symbols as
    possible that will make each expression given equal to zero. The format
    of the output can be controlled by using the ``dict`` or ``set`` keyword:

    >>> solve(x - 1, dict=True)
    [{x: 1}]
    >>> solve([x**2 - y, x + y - 6], set=True)
    ([x, y], {(-3, 9), (2, 4)})

    The following discussion provides an explanation for the output
    obtained when not using those keywords.

Empty List
----------

    When there is no solution, an empty list is returned.

    >>> solve(sqrt(x) + 1)  # or solve(sqrt(x) + 1, dict=True)
    []
    >>> solve(sqrt(x) + 1, set=True)
    ([x], set())

List Of Values
--------------

    A list of values is given when the symbol to solve for was
    unambiguous in context because a) the equation was univariate or
    b) a single symbol was specified as being of interest.

    >>> solve(x**2 - 4)
    [-2, 2]
    >>> solve(x - y - 1, x)
    [y + 1]

Single Dictionary
-----------------

    A single dictionary with keys being symbols and values being the solutions
    for those symbols is the result when equations are passed as a list and are
    all linear in the symbols given. Note: such a system is automatically generated
    for a single equation (not passed as a list) if there is an
    undetermined-coefficients solution for the symbols specified. If this is not
    what was intended, then pass the expression in a list.

    >>> solve([x + y - 2, x - y + 2], x, y)
    {x: 0, y: 2}
    >>> eq = a*x - 2*x + b - 5
    >>> solve(eq, {a, b})  # undetermined coefficients
    {a: 2, b: 5}
    >>> solve([eq], {a, b})  # algebraic
    {a: -b/x + (2*x + 5)/x}

List of Tuples
--------------

    Each tuple in the list gives a solution for the symbols in the order
    they were given. This format is used when a) a list of equations contains
    at least one nonlinear equation or b) a list of symbols is given in a
    well defined order. (This is also the format for the tuples in the set
    returned when using the flag ``set=True``.)

    >>> solve(x - 1, x, y)  # more than one symbol
    [(1, y)]
    >>> solve([x**2], x)  # list with nonlinear equation
    [(0,)]
    >>> solve([x**2 - 1], x)
    [(-1,), (1,)]
    >>> solve([x**2 - y, x - 3], x, y)  # nonlinear and multiple symbols
    [(3, 9)]

List of Dictionaries
--------------------

    The list of dictionaries is returned when the expression was not
    univariate or there was a nonlinear expression in a list *and* the order of
    symbols would otherwise be ambiguous because a) no symbols were passed or
    b) the symbols were passed as a set. (This is also the format selected with
    ``dict=True``.)

    >>> solve(x - y)
    [{x: y}]
    >>> solve([exp(x) - 1, x*(x - 1)])
    [{x: 0}]
    >>> system = [x + y - z, x**2 - y + z, exp(z) + 1/x + 1/y - 2]
    >>> sol = solve(system[:2]); sol
    [{x: -1, y: z + 1}, {x: 0, y: z}]

    The dictionaries only contain values that are distinct from the keys.
    In the last example above, there is no key for ``z`` in the dictionary
    since only *two* of the three equations were insufficient to determine its
    value. These solutions can be used to eliminate those variables from the
    third equation, however, to give a relationship in a single variable that
    can be solved (perhaps numerically) to obtain a full solution with the
    advantage of only needing to guess a single value instead of three.

        >>> from sympy import nsolve
        >>> [system[-1].subs(s) for s in sol]
        [exp(z) - 3 + 1/(z + 1), exp(z) + zoo + 1/z]
        >>> z_eq = _[0]
        >>> zsol = nsolve(z_eq, 1); zsol
        0.906425478894557
        >>> sol0 = {k: v.subs(z, zsol) for k, v in sol[0].items()}
        >>> sol0[z] = zsol; sol0
        {x: -1, y: 1.90642547889456, z: 0.906425478894557}

Boolean or Relational
---------------------

    A boolean expression is returned when a relational expression other
    than an :class:`~.Equality` is given as an expression to solve. A single `Equality`
    or a more complicated relational expression might be returned. The
    use of :func:`~.solve` here is equivalent to passing the equation set and
    symbols to :func:`~.reduce_inequalities` (and ``dict``, ``set``, and ``check``
    flags are ignored).

    >>> solve([x**2 > 4, x > 0])
    (2 < x) & (x < oo)

    >>> from sympy import Unequality as Ne
    >>> solve([x**2 - 4, Ne(x, -2)])
    Eq(x, 2)

    Any returned `Equality` can be converted to a dictionary:

    >>> {_.lhs: _.rhs}
    {x: 2}
