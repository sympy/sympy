Basic Usage
-----------

The solve_univariate_inequality function is used to solve inequalities in a single variable::

    >>> from sympy import symbols, solve_univariate_inequality
    >>> from sympy import S, Poly
    >>> x = symbols('x')
    >>> solve_univariate_inequality(x**2 >= 4, x)
    ((-oo, -2] | [2, oo))

    >>> solve_univariate_inequality(x**2 < 4, x)
    (-2, 2)

Linear Inequalities
-------------------

Linear inequalities can be solved::

    >>> solve_univariate_inequality(x + 1 > 0, x)
    (-1, oo)

    >>> solve_univariate_inequality(2*x - 3 <= 0, x)
    (-oo, 3/2]

Polynomial Inequalities
-----------------------

Higher-degree polynomial inequalities::

    >>> solve_univariate_inequality(x**3 - 3*x**2 + 2 > 0, x)
    ((-oo, 1) | (2, oo))

    >>> solve_univariate_inequality((x - 1)*(x + 1) < 0, x)
    (-1, 1)

Rational Inequalities
---------------------

Inequalities with rational functions::

    >>> solve_univariate_inequality((x - 1)/(x + 1) >= 0, x)
    ((-oo, -1) | [1, oo))

    >>> solve_univariate_inequality(1/(x - 2) > 0, x)
    (2, oo)

