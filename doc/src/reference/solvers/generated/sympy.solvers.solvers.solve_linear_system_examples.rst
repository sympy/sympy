Basic Usage
-----------

The solve_linear_system function is used to solve systems of linear equations represented by a matrix::

    >>> from sympy import symbols, Matrix, solve_linear_system
    >>> x, y = symbols('x, y')
    >>> system = Matrix(([1, 1, 3], [2, -1, 1]))
    >>> solve_linear_system(system, x, y)
    {x: 1, y: 2}

Systems with Multiple Solutions
-------------------------------

For systems with multiple solutions, free variables are introduced::

    >>> from sympy import symbols, Matrix, solve_linear_system
    >>> x, y, z = symbols('x, y, z')
    >>> system = Matrix(([1, 1, 1, 3], [2, 2, 2, 6]))
    >>> solve_linear_system(system, x, y, z)
    {x: -y - z + 3}

Homogeneous Systems
-------------------

For homogeneous systems::

    >>> from sympy import symbols, Matrix, solve_linear_system
    >>> x, y = symbols('x, y')
    >>> system = Matrix(([1, 1, 0], [2, 2, 0]))
    >>> solve_linear_system(system, x, y)
    {x: -y}

Inconsistent Systems
-------------------

For inconsistent systems, None is returned::

    >>> from sympy import symbols, Matrix, solve_linear_system
    >>> x, y = symbols('x, y')
    >>> system = Matrix(([1, 1, 3], [1, 1, 4]))
    >>> solve_linear_system(system, x, y)
    None 