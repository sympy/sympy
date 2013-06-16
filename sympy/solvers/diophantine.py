from sympy import degree_list
from sympy import igcd
from sympy import symbols
from sympy import Add
from sympy import Integer
from sympy import sign
from sympy import S


def diop_solve(eq):
    """
    Solves diophantine equations. Uses classify_diop() to determine the
    type of eqaution and calls the appropriate solver function.

    **Usage**

        diop_solve(eq) -> Solve diophantine equation, eq.

    **Details**

        ``eq`` should be an expression which is assumed to be zero.

    Examples
    ========

    >>> from sympy.solvers.diophantine import diop_solve
    >>> from sympy.abc import x, y, z, w
    >>> diop_solve(2*x + 3*y - 5)
    {x: 15*t - 5, y: -10*t + 5}
    >>> diop_solve(4*x + 3*y -4*z + 5)
    {x: -15*t + 4*z - 5, y: 20*t - 4*z + 5, z: z}
    >>> diop_solve(x + 3*y - 4*z + w -6)
    {w: 6*t, x: -6*t - 3*y + 4*z + 6, y: y, z: z}
    """
    var, coeff, t = classify_diop(eq)

    if t == "linear":
        return diop_linear(var, coeff)


def classify_diop(eq):
    """
    Helper routine used by diop_solve(). Returns a tuple containing the type of the
    diophantine equation along with the variables(free symbols) and their coefficients.
    Variables are returned as a list and coefficients are returned as a dict with
    the key being the variable name and the constant term is keyed to Integer(1).
    Type is an element in the set {"linear", "quadratic", "pell", "pythogorean", "exponential"}

    **Usage**

        classify_diop(eq) -> Return variables, coefficients and type in order.

    **Details**

        ``eq`` should be an expression which is assumed to be zero.

    Examples
    ========

    >>> from sympy.solvers.diophantine import classify_diop
    >>> from sympy.abc import x, y, z, w, t
    >>> classify_diop(4*x + 6*y - 4)
    ([x, y], {1: -4, x: 4, y: 6}, 'linear')
    >>> classify_diop(x + 3*y -4*z + 5)
    ([x, y, z], {1: 5, x: 1, y: 3, z: -4}, 'linear')
    """
    var = list(eq.free_symbols)
    var.sort()

    if max(degree_list(eq)) == 1:
        diop_type = "linear"
    else:
        raise NotImplementedError("Still not implemented")


    coeff = dict([reversed(t.as_independent(*var)) for t in eq.args])
    for v in coeff:
        if not isinstance(coeff[v], Integer):
            raise TypeError("Coefficients should be Integers")

    return var, coeff, diop_type,


def diop_linear(var, coeff):
    """
    Solves linear diophantine equations.

    **Usage**

        diop_linear(var, coeff) -> var is a list of variables and coeff is a dictionary
        containing coefficients of the symbols.

    **Details**

        ``var`` list of the variables in the equation.
        ``coeff`` dictionary containing coefficients of each variable.
            coefficients are keyed by the variable and the constant term is keyed with
            Integer(1).

    Examples
    ========

    >>> from sympy.solvers.diophantine import diop_linear
    >>> from sympy.abc import x, y, z
    >>> from sympy import Integer
    >>> diop_linear([x, y],{Integer(1): -5, x: 2, y:-3}) #solves equation 2*x - 3*y -5 = 0
    {x: -15*t - 5, y: -10*t - 5}
    >>> diop_linear([x, y, z], {Integer(1): -3, x: 2, y: -3, z: -4}) # 2*x - 3*y - 4*z - 3= 0
    {x: -9*t - 4*z - 3, y: -6*t - 4*z - 3, z: z}
    """
    x = var[0]; y = var[1]
    a = coeff[x]; b = coeff[y]

    if len(var) == len(coeff):
        c = 0
    else:
        c = -coeff[Integer(1)]

    if len(var) == 2:
        sol_x, sol_y = base_solution_linear(c, a, b, True)
        return {x: sol_x, y: sol_y}

    elif len(var) > 2:
        X = []; Y = []

        for v in var[2:]:
            sol_x, sol_y  = base_solution_linear(-coeff[v], a, b)
            X.append(sol_x*v); Y.append(sol_y*v)

        sol_x, sol_y = base_solution_linear(c, a, b, True)
        X.append(sol_x); Y.append(sol_y)

        l = []
        if None not in X and None not in Y:
            l.append((x, Add(*X))); l.append((y, Add(*Y)))
            for v in var[2:]:
                l.append((v, v))
        else:
            for v in var:
                l.append((v, None))

        return dict(l)


def base_solution_linear(c, a, b,param=False):
    """
    Return the base solution for a linear diophantine equation with two
    variables. Called repeatedly by diop_linear().

    **Usage**

        base_solution_linear(c, a, b, param) -> a, b, c are Integers as in a*x + b*y = c
        and param is a boolean value to set parameterized solution.

    **Details**

        ``c`` is the constant term in a*x + b*y = c
        ``a`` is the integer coefficient of x in a*x + b*y = c
        ``b`` is the integer coefficient of y in a*x + b*y = c
        ``param`` is a boolean value used to signal the parameterized solution.
            If set True, returns the general parameterized solution, otherwise return
            a basic solution to the equation.


    Examples
    ========

    >>> from sympy.solvers.diophantine import base_solution_linear
    >>> base_solution_linear(5, 2, 3, False) # equation 2*x + 3*y = 5
    (-5, 5)
    >>> base_solution_linear(0, 5, 7, False) # equation 5*x + 7*y = 0
    (0, 0)
    >>> base_solution_linear(5, 2, 3, True) # equation 2*x + 3*y = 5
    (15*t - 5, -10*t + 5)
    >>> base_solution_linear(0, 5, 7, True) # equation 5*x + 7*y = 0
    (7*t, -5*t)
    """
    t = symbols("t", type = Integer)

    d = igcd(a, igcd(b, c))
    a = a // d; b = b // d; c = c // d

    if c == 0:
        if param:
            return (b*t , -a*t)
        else:
            return (S.Zero, S.Zero)
    else:
        x0, y0, d = extended_euclid(int(abs(a)), int(abs(b)))

        x0 = x0 * sign(a)
        y0 = y0 * sign(b)

        if d == igcd(c, d):
            if param:
                return (c*(x0 + b*t), c*(y0 - a*t))
            else:
                return (Integer(c*x0), Integer(c*y0))
        else:
            return (None, None)


def extended_euclid(a, b):
    """
    For given a, b returns a tuple containing integers x, y and d such that
    a*x + b*y = d. Here d = gcd(a, b).

    **Usage**
        extended_euclid(a, b) -> returns x, y and gcd(a, b)

    **Details**
        ``a`` Any instance of Integer
        ``b`` Any instance of Integer

    Examples
    ========

    >>> from sympy.solvers.diophantine import extended_euclid
    >>> extended_euclid(4, 6)
    (-1, 1, 2)
    >>> extended_euclid(3, 5)
    (2, -1, 1)
    """
    if b == 0:
        return (1, 0, a)

    x0, y0, d = extended_euclid(b, a%b)
    x, y = y0, x0 - (a//b) * y0

    return x, y, d
