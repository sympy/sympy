from __future__ import print_function, division

from sympy import (Poly, igcd, divisors, sign, symbols, S, Integer, Wild, Symbol, factorint,
    Add, Mul, solve, ceiling, floor, sqrt, sympify, Subs, ilcm, Matrix, factor_list, perfect_power,
    isprime, nextprime, integer_nthroot)

from sympy.core.function import _mexpand
from sympy.simplify.radsimp import rad_rationalize
from sympy.utilities import default_sort_key, numbered_symbols
from sympy.core.numbers import igcdex
from sympy.ntheory.residue_ntheory import sqrt_mod
from sympy.core.compatibility import range
from sympy.core.relational import Eq
from sympy.solvers.solvers import check_assumptions

__all__ = ['diophantine', 'diop_solve', 'classify_diop', 'diop_linear', 'base_solution_linear',
'diop_quadratic', 'diop_DN', 'cornacchia', 'diop_bf_DN', 'transformation_to_DN', 'find_DN',
'diop_ternary_quadratic',  'square_factor', 'descent', 'diop_general_pythagorean',
'diop_general_sum_of_squares', 'partition', 'sum_of_three_squares', 'sum_of_four_squares']


def diophantine(eq, param=symbols("t", integer=True)):
    """
    Simplify the solution procedure of diophantine equation ``eq`` by
    converting it into a product of terms which should equal zero.

    For example, when solving, `x^2 - y^2 = 0` this is treated as
    `(x + y)(x - y) = 0` and `x+y = 0` and `x-y = 0` are solved independently
    and combined. Each term is solved by calling ``diop_solve()``.

    Output of ``diophantine()`` is a set of tuples. Each tuple represents a
    solution of the input equation. In a tuple, solution for each variable is
    listed according to the alphabetic order of input variables. i.e. if we have
    an equation with two variables `a` and `b`, first element of the tuple will
    give the solution for `a` and the second element will give the solution for
    `b`.

    Usage
    =====

    ``diophantine(eq, t)``: Solve the diophantine equation ``eq``.
    ``t`` is the parameter to be used by ``diop_solve()``.

    Details
    =======

    ``eq`` should be an expression which is assumed to be zero.
    ``t`` is the parameter to be used in the solution.

    Examples
    ========

    >>> from sympy.solvers.diophantine import diophantine
    >>> from sympy.abc import x, y, z
    >>> diophantine(x**2 - y**2)
    set([(-t, -t), (t, -t)])

    #>>> diophantine(x*(2*x + 3*y - z))
    #set([(0, n1, n2), (3*t - z, -2*t + z, z)])
    #>>> diophantine(x**2 + 3*x*y + 4*x)
    #set([(0, n1), (3*t - 4, -t)])

    See Also
    ========

    diop_solve()
    """
    if isinstance(eq, Eq):
        eq = eq.lhs - eq.rhs

    eq = Poly(eq).as_expr()
    if not eq.is_polynomial() or eq.is_number:
        raise TypeError("Equation input format not supported")

    var = list(eq.expand(force=True).free_symbols)
    var.sort(key=default_sort_key)

    terms = factor_list(eq)[1]

    sols = set([])

    for term in terms:

        base = term[0]

        var_t, jnk, eq_type = classify_diop(base)
        solution = diop_solve(base, param)

        if eq_type in ["linear", "homogeneous_ternary_quadratic", "general_pythagorean"]:
            if merge_solution(var, var_t, solution) != ():
                sols.add(merge_solution(var, var_t, solution))

        elif eq_type in ["binary_quadratic",  "general_sum_of_squares", "univariate"]:
            for sol in solution:
                if merge_solution(var, var_t, sol) != ():
                    sols.add(merge_solution(var, var_t, sol))

    return sols


def merge_solution(var, var_t, solution):
    """
    This is used to construct the full solution from the solutions of sub
    equations.

    For example when solving the equation `(x - y)(x^2 + y^2 - z^2) = 0`,
    solutions for each of the equations `x-y = 0` and `x^2 + y^2 - z^2` are
    found independently. Solutions for `x - y = 0` are `(x, y) = (t, t)`. But
    we should introduce a value for z when we output the solution for the
    original equation. This function converts `(t, t)` into `(t, t, n_{1})`
    where `n_{1}` is an integer parameter.
    """
    l = []

    if None in solution:
        return ()

    solution = iter(solution)
    params = numbered_symbols("n", Integer=True, start=1)
    for v in var:
        if v in var_t:
            l.append(next(solution))
        else:
            l.append(next(params))

    for val, symb in zip(l, var):
        if check_assumptions(val, **symb.assumptions0) is False:
            return tuple()

    return tuple(l)


def diop_solve(eq, param=symbols("t", integer=True)):
    """
    Solves the diophantine equation ``eq``.

    Similar to ``diophantine()`` but doesn't try to factor ``eq`` as latter
    does. Uses ``classify_diop()`` to determine the type of the eqaution and
    calls the appropriate solver function.

    Usage
    =====

    ``diop_solve(eq, t)``: Solve diophantine equation, ``eq`` using ``t``
    as a parameter if needed.

    Details
    =======

    ``eq`` should be an expression which is assumed to be zero.
    ``t`` is a parameter to be used in the solution.

    Examples
    ========

    >>> from sympy.solvers.diophantine import diop_solve
    >>> from sympy.abc import x, y, z, w
    >>> diop_solve(2*x + 3*y - 5)
    (3*t - 5, -2*t + 5)
    >>> diop_solve(4*x + 3*y -4*z + 5)
    (3*t + 4*z - 5, -4*t - 4*z + 5,  z)
    >>> diop_solve(x + 3*y - 4*z + w -6)
    (t, -t - 3*y + 4*z + 6, y, z)
    >>> diop_solve(x**2 + y**2 - 5)
    set([(-2, -1), (-2, 1), (2, -1), (2, 1)])

    See Also
    ========

    diophantine()
    """
    var, coeff, eq_type = classify_diop(eq)

    if eq_type == "linear":
        return _diop_linear(var, coeff, param)

    elif eq_type == "binary_quadratic":
        return _diop_quadratic(var, coeff, param)

    elif eq_type == "homogeneous_ternary_quadratic":
        x_0, y_0, z_0 = _diop_ternary_quadratic(var, coeff)
        return _parametrize_ternary_quadratic((x_0, y_0, z_0), var, coeff)

    elif eq_type == "general_pythagorean":
        return _diop_general_pythagorean(var, coeff, param)

    elif eq_type == "univariate":
        l = solve(eq)
        s = set([])

        for soln in l:
            if isinstance(soln, Integer):
                s.add((soln,))
        return s

    elif eq_type == "general_sum_of_squares":
        return _diop_general_sum_of_squares(var, coeff)


def classify_diop(eq):
    """
    Helper routine used by diop_solve() to find the type of the ``eq`` etc.

    Returns a tuple containing the type of the diophantine equation along with
    the variables(free symbols) and their coefficients. Variables are returned
    as a list and coefficients are returned as a dict with the key being the
    respective term and the constant term is keyed to Integer(1). Type is an
    element in the set {"linear", "binary_quadratic", "general_pythagorean",
    "homogeneous_ternary_quadratic", "univariate", "general_sum_of_squares"}

    Usage
    =====

    ``classify_diop(eq)``: Return variables, coefficients and type of the
    ``eq``.

    Details
    =======

    ``eq`` should be an expression which is assumed to be zero.

    Examples
    ========

    >>> from sympy.solvers.diophantine import classify_diop
    >>> from sympy.abc import x, y, z, w, t
    >>> classify_diop(4*x + 6*y - 4)
    ([x, y], {1: -4, x: 4, y: 6}, 'linear')
    >>> classify_diop(x + 3*y -4*z + 5)
    ([x, y, z], {1: 5, x: 1, y: 3, z: -4}, 'linear')
    >>> classify_diop(x**2 + y**2 - x*y + x + 5)
    ([x, y], {1: 5, x: 1, x**2: 1, y: 0, y**2: 1, x*y: -1}, 'binary_quadratic')
    """
    eq = eq.expand(force=True)
    var = list(eq.free_symbols)
    var.sort(key=default_sort_key)

    coeff = {}
    diop_type = None

    coeff = dict([reversed(t.as_independent(*var)) for t in eq.args])
    for v in coeff:
        if not isinstance(coeff[v], Integer):
            raise TypeError("Coefficients should be Integers")

    if len(var) == 1:
        diop_type = "univariate"
    elif Poly(eq).total_degree() == 1:
        diop_type = "linear"
    elif Poly(eq).total_degree() == 2 and len(var) == 2:
        diop_type = "binary_quadratic"
        x, y = var[:2]

        if isinstance(eq, Mul):
            coeff = {x**2: 0, x*y: eq.args[0], y**2: 0, x: 0, y: 0, Integer(1): 0}
        else:
            for term in [x**2, y**2, x*y, x, y, Integer(1)]:
                if term not in coeff.keys():
                    coeff[term] = Integer(0)

    elif Poly(eq).total_degree() == 2 and len(var) == 3 and Integer(1) not in coeff.keys():
        for v in var:
            if v in coeff.keys():
                diop_type = "inhomogeneous_ternary_quadratic"
                break
        else:
            diop_type = "homogeneous_ternary_quadratic"

            x, y, z = var[:3]

            for term in [x**2, y**2, z**2, x*y, y*z, x*z]:
                if term not in coeff.keys():
                    coeff[term] = Integer(0)

    elif Poly(eq).degree() == 2 and len(var) >= 3:

        for v in var:
            if v in coeff.keys():
                diop_type = "inhomogeneous_general_quadratic"
                break

        else:
            if Integer(1) in coeff.keys():
                constant_term = True
            else:
                constant_term = False

            non_square_degree_2_terms = False
            for v in var:
                for u in var:
                    if u != v and u*v in coeff.keys():
                        non_square_degree_2_terms = True
                        break
                if non_square_degree_2_terms:
                    break

            if constant_term and non_square_degree_2_terms:
                diop_type = "inhomogeneous_general_quadratic"

            elif constant_term and not non_square_degree_2_terms:
                for v in var:
                    if coeff[v**2] != 1:
                        break
                else:
                    diop_type = "general_sum_of_squares"

            elif not constant_term and non_square_degree_2_terms:
                diop_type = "homogeneous_general_quadratic"

            else:
                coeff_sign_sum = 0

                for v in var:
                    if not isinstance(sqrt(abs(Integer(coeff[v**2]))), Integer):
                        break
                    coeff_sign_sum = coeff_sign_sum + sign(coeff[v**2])
                else:
                    if abs(coeff_sign_sum) == len(var) - 2 and not constant_term:
                        diop_type = "general_pythagorean"

    elif Poly(eq).total_degree() == 3 and len(var) == 2:

        x, y = var[:2]
        diop_type = "cubic_thue"

        for term in [x**3, x**2*y, x*y**2, y**3, Integer(1)]:
            if term not in coeff.keys():
                coeff[term] == Integer(0)

    if diop_type is not None:
        return var, coeff, diop_type
    else:
        raise NotImplementedError("Still not implemented")


def diop_linear(eq, param=symbols("t", integer=True)):
    """
    Solves linear diophantine equations.

    A linear diophantine equation is an equation of the form `a_{1}x_{1} +
    a_{2}x_{2} + .. + a_{n}x_{n} = 0` where `a_{1}, a_{2}, ..a_{n}` are
    integer constants and `x_{1}, x_{2}, ..x_{n}` are integer variables.

    Usage
    =====

    ``diop_linear(eq)``: Returns a tuple containing solutions to the
    diophantine equation ``eq``. Values in the tuple is arranged in the same
    order as the sorted variables.

    Details
    =======

    ``eq`` is a linear diophantine equation which is assumed to be zero.
    ``param`` is the parameter to be used in the solution.

    Examples
    ========

    >>> from sympy.solvers.diophantine import diop_linear
    >>> from sympy.abc import x, y, z, t
    >>> from sympy import Integer
    >>> diop_linear(2*x - 3*y - 5) #solves equation 2*x - 3*y -5 = 0
    (-3*t - 5, -2*t - 5)

    Here x = -3*t - 5 and y = -2*t - 5

    >>> diop_linear(2*x - 3*y - 4*z -3)
    (-3*t - 4*z - 3, -2*t - 4*z - 3,  z)

    See Also
    ========

    diop_quadratic(), diop_ternary_quadratic(), diop_general_pythagorean(),
    diop_general_sum_of_squares()
    """
    var, coeff, diop_type = classify_diop(eq)

    if diop_type == "linear":
        return _diop_linear(var, coeff, param)


def _diop_linear(var, coeff, param):
    """
    Break down the multivariate equation into multiple
    bivariate equations of the form:

    ax + by == d

    which can then be solved using base_solution_linear().

    Example:

    a_0*x_0 + a_1*x_1 + a_2*x_2 == c becomes:

    a_0*x_0 + g_0*y_0 == c 

    where g_0 == gcd(a_1, a_2) and
          
          y == a_1*x_1 + a_2*x_2
               ---       ---
               g_0       g_0

    Then we can solve for x_0, y_0 with base_solution_linear().

    x_0 is appended to our return value, while y_0 is used to 
    solve for x_1 and x_2 using base_solution_linear() again.
    """
    
    if len(var) < 2:
        raise ValueError("Fewer than two elements provided for 'var' at _diop_linear()")

    if len(var) == len(coeff):
        c = 0
    else:
        #coeff[] is negated because input is of the form: ax + by - c == 0
        #                                 but is used as: ax + by     == c
        c = -coeff[Integer(1)]

    A = [coeff[v] for v in var]
    B = []
    if len(var) > 2:
        B.append(igcd(A[-2], A[-1]))
        A[-2] = A[-2] // B[0]
        A[-1] = A[-1] // B[0]
        for i in range(len(A) - 3, 0, -1):
            gcd = igcd(B[0], A[i])
            B[0] = B[0] // gcd
            A[i] = A[i] // gcd
            B.insert(0, gcd)
    B.append(A[-1])

    # Some solutions will have multiple free variables in their solutions.
    params = [symbols(str(param) + "_" + str(i)) for i in range(0, len(var))]

    solutions = []
    no_solution = tuple([None] * len(var))
    for i in range(0, len(B)):
        tot_x, tot_y = 0, 0

        if type(c) is Add: 
            # example: 5 + t_0 + 3*t_1
            args = c.args
        else: # c is a Mul, a Symbol, or an Integer
            args = [c]

        for j in range(0, len(args)):
            arg_type = type(args[j])
            if arg_type is Mul:
                # example: 3*t_1 -> k = 3
                k = args[j].as_two_terms()[0]
                param_index = params.index(args[j].as_two_terms()[1]) + 1
            elif arg_type is Symbol: 
                # example: t_0 -> k = 1
                k = 1
                param_index = params.index(args[j]) + 1
            else: #arg_type is Integer
                # example: 5 -> k = 5
                k = args[j]
                param_index = 0

            sol_x, sol_y = base_solution_linear(k, A[i], B[i], params[param_index])
            if arg_type is Mul or arg_type is Symbol:
                if sol_x is None:
                    sol_x = 0
                elif type(sol_x) is Add:
                    sol_x = sol_x.args[0]*params[param_index - 1] + sol_x.args[1]
                elif type(sol_x) is Integer:
                    sol_x = sol_x*params[param_index - 1]
                
                if sol_y is None:
                    sol_y = 0
                elif type(sol_y) is Add:
                    sol_y = sol_y.args[0]*params[param_index - 1] + sol_y.args[1]
                elif type(sol_y) is Integer:
                    sol_y = sol_y*params[param_index - 1]
            
            else:
                if sol_x is None or sol_y is None:
                    return no_solution
            
            tot_x += sol_x
            tot_y += sol_y

        solutions.append(tot_x)
        c = tot_y

    solutions.append(tot_y)

    return tuple(solutions)


def base_solution_linear(c, a, b, t=None):
    """
    Return the base solution for a linear diophantine equation with two
    variables.

    Used by ``diop_linear()`` to find the base solution of a linear
    Diophantine equation. If ``t`` is given then the parametrized solution is
    returned.

    Usage
    =====

    ``base_solution_linear(c, a, b, t)``: ``a``, ``b``, ``c`` are coefficients
    in `ax + by = c` and ``t`` is the parameter to be used in the solution.

    Examples
    ========

    >>> from sympy.solvers.diophantine import base_solution_linear
    >>> from sympy.abc import t
    >>> base_solution_linear(5, 2, 3) # equation 2*x + 3*y = 5
    (-5, 5)
    >>> base_solution_linear(0, 5, 7) # equation 5*x + 7*y = 0
    (0, 0)
    >>> base_solution_linear(5, 2, 3, t) # equation 2*x + 3*y = 5
    (3*t - 5, -2*t + 5)
    >>> base_solution_linear(0, 5, 7, t) # equation 5*x + 7*y = 0
    (7*t, -5*t)
    """
    d = igcd(a, igcd(b, c))
    a = a // d
    b = b // d
    c = c // d

    if c == 0:
        if t != None:
            return (b*t , -a*t)
        else:
            return (S.Zero, S.Zero)
    else:
        x0, y0, d = extended_euclid(int(abs(a)), int(abs(b)))

        x0 = x0 * sign(a)
        y0 = y0 * sign(b)

        if divisible(c, d):
            if t != None:
                return (c*x0 + b*t, c*y0 - a*t)
            else:
                return (Integer(c*x0), Integer(c*y0))
        else:
            return (None, None)


def extended_euclid(a, b):
    """
    For given ``a``, ``b`` returns a tuple containing integers `x`, `y` and `d`
    such that `ax + by = d`. Here `d = gcd(a, b)`.

    Usage
    =====

    ``extended_euclid(a, b)``: returns `x`, `y` and `\gcd(a, b)`.

    Details
    =======

    ``a`` Any instance of Integer.
    ``b`` Any instance of Integer.

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


def divisible(a, b):
    """
    Returns `True` if ``a`` is divisible by ``b`` and `False` otherwise.
    """
    return igcd(int(a), int(b)) == abs(int(b))


def diop_quadratic(eq, param=symbols("t", integer=True)):
    """
    Solves quadratic diophantine equations.

    i.e. equations of the form `Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0`. Returns a
    set containing the tuples `(x, y)` which contains the solutions. If there
    are no solutions then `(None, None)` is returned.

    Usage
    =====

    ``diop_quadratic(eq, param)``: ``eq`` is a quadratic binary diophantine
    equation. ``param`` is used to indicate the parameter to be used in the
    solution.

    Details
    =======

    ``eq`` should be an expression which is assumed to be zero.
    ``param`` is a parameter to be used in the solution.

    Examples
    ========

    >>> from sympy.abc import x, y, t
    >>> from sympy.solvers.diophantine import diop_quadratic
    >>> diop_quadratic(x**2 + y**2 + 2*x + 2*y + 2, t)
    set([(-1, -1)])

    References
    ==========

    .. [1] Methods to solve Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0,[online],
          Available: http://www.alpertron.com.ar/METHODS.HTM
    .. [2] Solving the equation ax^2+ bxy + cy^2 + dx + ey + f= 0, [online],
          Available: http://www.jpr2718.org/ax2p.pdf

    See Also
    ========

    diop_linear(), diop_ternary_quadratic(), diop_general_sum_of_squares(),
    diop_general_pythagorean()
    """
    var, coeff, diop_type = classify_diop(eq)

    if diop_type == "binary_quadratic":
        return _diop_quadratic(var, coeff, param)


def _diop_quadratic(var, coeff, t):

    x, y = var[:2]

    for term in [x**2, y**2, x*y, x, y, Integer(1)]:
        if term not in coeff.keys():
            coeff[term] = Integer(0)

    A = coeff[x**2]
    B = coeff[x*y]
    C = coeff[y**2]
    D = coeff[x]
    E = coeff[y]
    F = coeff[Integer(1)]

    d = igcd(A, igcd(B, igcd(C, igcd(D, igcd(E, F)))))
    A = A // d
    B = B // d
    C = C // d
    D = D // d
    E = E // d
    F = F // d

    # (1) Linear case: A = B = C = 0 ==> considered under linear diophantine equations

    # (2) Simple-Hyperbolic case:A = C = 0, B != 0
    # In this case equation can be converted to (Bx + E)(By + D) = DE - BF
    # We consider two cases; DE - BF = 0 and DE - BF != 0
    # More details, http://www.alpertron.com.ar/METHODS.HTM#SHyperb

    l = set([])

    if A == 0 and C == 0 and B != 0:

        if D*E - B*F == 0:
            if divisible(int(E), int(B)):
                l.add((-E/B, t))
            if divisible(int(D), int(B)):
                l.add((t, -D/B))

        else:
            div = divisors(D*E - B*F)
            div = div + [-term for term in div]

            for d in div:
                if divisible(int(d - E), int(B)):
                    x0  = (d - E) // B
                    if divisible(int(D*E - B*F), int(d)):
                        if divisible(int((D*E - B*F)// d - D), int(B)):
                            y0 = ((D*E - B*F) // d - D) // B
                            l.add((x0, y0))

    # (3) Parabolic case: B**2 - 4*A*C = 0
    # There are two subcases to be considered in this case.
    # sqrt(c)D - sqrt(a)E = 0 and sqrt(c)D - sqrt(a)E != 0
    # More Details, http://www.alpertron.com.ar/METHODS.HTM#Parabol

    elif B**2 - 4*A*C == 0:

        if A == 0:
            s = _diop_quadratic([y, x], coeff, t)
            for soln in s:
                l.add((soln[1], soln[0]))

        else:
            g = igcd(A, C)
            g = abs(g) * sign(A)
            a = A // g
            b = B // g
            c = C // g
            e = sign(B/A)


            if e*sqrt(c)*D - sqrt(a)*E == 0:
                z = symbols("z", real=True)
                roots = solve(sqrt(a)*g*z**2 + D*z + sqrt(a)*F)
                for root in roots:
                    if isinstance(root, Integer):
                        l.add((diop_solve(sqrt(a)*x + e*sqrt(c)*y - root)[0], diop_solve(sqrt(a)*x + e*sqrt(c)*y - root)[1]))

            elif isinstance(e*sqrt(c)*D - sqrt(a)*E, Integer):
                solve_x = lambda u: e*sqrt(c)*g*(sqrt(a)*E - e*sqrt(c)*D)*t**2 - (E + 2*e*sqrt(c)*g*u)*t\
                    - (e*sqrt(c)*g*u**2 + E*u + e*sqrt(c)*F) // (e*sqrt(c)*D - sqrt(a)*E)

                solve_y = lambda u: sqrt(a)*g*(e*sqrt(c)*D - sqrt(a)*E)*t**2 + (D + 2*sqrt(a)*g*u)*t \
                    + (sqrt(a)*g*u**2 + D*u + sqrt(a)*F) // (e*sqrt(c)*D - sqrt(a)*E)

                for z0 in range(0, abs(e*sqrt(c)*D - sqrt(a)*E)):
                    if divisible(sqrt(a)*g*z0**2 + D*z0 + sqrt(a)*F, e*sqrt(c)*D - sqrt(a)*E):
                        l.add((solve_x(z0), solve_y(z0)))

    # (4) Method used when B**2 - 4*A*C is a square, is descibed in p. 6 of the below paper
    # by John P. Robertson.
    # http://www.jpr2718.org/ax2p.pdf

    elif isinstance(sqrt(B**2 - 4*A*C), Integer):
        if A != 0:
            r = sqrt(B**2 - 4*A*C)
            u, v = symbols("u, v", integer=True)
            eq = _mexpand(4*A*r*u*v + 4*A*D*(B*v + r*u + r*v - B*u) + 2*A*4*A*E*(u - v) + 4*A*r*4*A*F)

            sol = diop_solve(eq, t)
            sol = list(sol)

            for solution in sol:
                s0 = solution[0]
                t0 = solution[1]

                x_0 = S(B*t0 + r*s0 + r*t0 - B*s0)/(4*A*r)
                y_0 = S(s0 - t0)/(2*r)

                if isinstance(s0, Symbol) or isinstance(t0, Symbol):
                    if check_param(x_0, y_0, 4*A*r, t) != (None, None):
                        l.add((check_param(x_0, y_0, 4*A*r, t)[0], check_param(x_0, y_0, 4*A*r, t)[1]))

                elif divisible(B*t0 + r*s0 + r*t0 - B*s0, 4*A*r):
                    if divisible(s0 - t0, 2*r):
                        if is_solution_quad(var, coeff, x_0, y_0):
                            l.add((x_0, y_0))
        else:
            _var = var
            _var[0], _var[1] = _var[1], _var[0] # Interchange x and y
            s = _diop_quadratic(_var, coeff, t)

            while len(s) > 0:
                sol = s.pop()
                l.add((sol[1], sol[0]))


    # (5) B**2 - 4*A*C > 0 and B**2 - 4*A*C not a square or B**2 - 4*A*C < 0

    else:

        P, Q = _transformation_to_DN(var, coeff)
        D, N = _find_DN(var, coeff)
        solns_pell = diop_DN(D, N)

        if D < 0:
            for solution in solns_pell:
                for X_i in [-solution[0], solution[0]]:
                    for Y_i in [-solution[1], solution[1]]:
                        x_i, y_i = (P*Matrix([X_i, Y_i]) + Q)[0], (P*Matrix([X_i, Y_i]) + Q)[1]
                        if isinstance(x_i, Integer) and isinstance(y_i, Integer):
                            l.add((x_i, y_i))

        else:
            # In this case equation can be transformed into a Pell equation
            #n = symbols("n", integer=True)

            fund_solns = solns_pell
            solns_pell = set(fund_solns)
            for X, Y in fund_solns:
                solns_pell.add((-X, -Y))

            a = diop_DN(D, 1)
            T = a[0][0]
            U = a[0][1]

            if (isinstance(P[0], Integer) and isinstance(P[1], Integer) and isinstance(P[2], Integer)
                and isinstance(P[3], Integer) and isinstance(Q[0], Integer) and isinstance(Q[1], Integer)):

                for sol in solns_pell:

                    r = sol[0]
                    s = sol[1]
                    x_n = S((r + s*sqrt(D))*(T + U*sqrt(D))**t + (r - s*sqrt(D))*(T - U*sqrt(D))**t)/2
                    y_n = S((r + s*sqrt(D))*(T + U*sqrt(D))**t - (r - s*sqrt(D))*(T - U*sqrt(D))**t)/(2*sqrt(D))

                    x_n = _mexpand(x_n)
                    y_n = _mexpand(y_n)
                    x_n, y_n = (P*Matrix([x_n, y_n]) + Q)[0], (P*Matrix([x_n, y_n]) + Q)[1]

                    l.add((x_n, y_n))

            else:
                L = ilcm(S(P[0]).q, ilcm(S(P[1]).q, ilcm(S(P[2]).q,
                         ilcm(S(P[3]).q, ilcm(S(Q[0]).q, S(Q[1]).q)))))

                k = 1

                T_k = T
                U_k = U

                while (T_k - 1) % L != 0 or U_k % L != 0:
                    T_k, U_k = T_k*T + D*U_k*U, T_k*U + U_k*T
                    k += 1

                for X, Y in solns_pell:

                    for i in range(k):
                        Z = P*Matrix([X, Y]) + Q
                        x, y = Z[0], Z[1]

                        if isinstance(x, Integer) and isinstance(y, Integer):
                            Xt = S((X + sqrt(D)*Y)*(T_k + sqrt(D)*U_k)**t +
                                  (X - sqrt(D)*Y)*(T_k - sqrt(D)*U_k)**t)/ 2
                            Yt = S((X + sqrt(D)*Y)*(T_k + sqrt(D)*U_k)**t -
                                  (X - sqrt(D)*Y)*(T_k - sqrt(D)*U_k)**t)/ (2*sqrt(D))
                            Zt = P*Matrix([Xt, Yt]) + Q
                            l.add((Zt[0], Zt[1]))

                        X, Y = X*T + D*U*Y, X*U + Y*T


    return l


def is_solution_quad(var, coeff, u, v):
    """
    Check whether `(u, v)` is solution to the quadratic binary diophantine
    equation with the variable list ``var`` and coefficient dictionary
    ``coeff``.

    Not intended for use by normal users.
    """
    x, y = var[:2]

    eq = x**2*coeff[x**2] + x*y*coeff[x*y] + y**2*coeff[y**2] + x*coeff[x] + y*coeff[y] + coeff[Integer(1)]

    return _mexpand(Subs(eq, (x, y), (u, v)).doit()) == 0


def diop_DN(D, N, t=symbols("t", integer=True)):
    """
    Solves the equation `x^2 - Dy^2 = N`.

    Mainly concerned in the case `D > 0, D` is not a perfect square, which is
    the same as generalized Pell equation. To solve the generalized Pell
    equation this function Uses LMM algorithm. Refer [1]_ for more details on
    the algorithm.
    Returns one solution for each class of the solutions. Other solutions of
    the class can be constructed according to the values of ``D`` and ``N``.
    Returns a list containing the solution tuples `(x, y)`.

    Usage
    =====

    ``diop_DN(D, N, t)``: D and N are integers as in `x^2 - Dy^2 = N` and
    ``t`` is the parameter to be used in the solutions.

    Details
    =======

    ``D`` and ``N`` correspond to D and N in the equation.
    ``t`` is the parameter to be used in the solutions.

    Examples
    ========

    >>> from sympy.solvers.diophantine import diop_DN
    >>> diop_DN(13, -4) # Solves equation x**2 - 13*y**2 = -4
    [(3, 1), (393, 109), (36, 10)]

    The output can be interpreted as follows: There are three fundamental
    solutions to the equation `x^2 - 13y^2 = -4` given by (3, 1), (393, 109)
    and (36, 10). Each tuple is in the form (x, y), i. e solution (3, 1) means
    that `x = 3` and `y = 1`.

    >>> diop_DN(986, 1) # Solves equation x**2 - 986*y**2 = 1
    [(49299, 1570)]

    See Also
    ========

    find_DN(), diop_bf_DN()

    References
    ==========

    .. [1] Solving the generalized Pell equation x**2 - D*y**2 = N, John P.
        Robertson, July 31, 2004, Pages 16 - 17. [online], Available:
        http://www.jpr2718.org/pell.pdf
    """
    if D < 0:
        if N == 0:
            return [(S.Zero, S.Zero)]
        elif N < 0:
            return []
        elif N > 0:
            d = divisors(square_factor(N))
            sol = []

            for divisor in d:
                sols = cornacchia(1, -D, N // divisor**2)
                if sols:
                    for x, y in sols:
                        sol.append((divisor*x, divisor*y))

            return sol

    elif D == 0:
        if N < 0 or not isinstance(sqrt(N), Integer):
            return []
        if N == 0:
            return [(S.Zero, t)]
        if isinstance(sqrt(N), Integer):
            return [(sqrt(N), t)]

    else: # D > 0
        if isinstance(sqrt(D), Integer):
            r = sqrt(D)

            if N == 0:
                return [(r*t, t)]
            else:
                sol = []

                for y in range(floor(sign(N)*(N - 1)/(2*r)) + 1):
                    if isinstance(sqrt(D*y**2 + N), Integer):
                        sol.append((sqrt(D*y**2 + N), y))

                return sol
        else:
            if N == 0:
                return [(S.Zero, S.Zero)]

            elif abs(N) == 1:

                pqa = PQa(0, 1, D)
                a_0 = floor(sqrt(D))
                l = 0
                G = []
                B = []

                for i in pqa:

                    a = i[2]
                    G.append(i[5])
                    B.append(i[4])

                    if l != 0 and a == 2*a_0:
                        break
                    l = l + 1

                if l % 2 == 1:

                    if N == -1:
                        x = G[l-1]
                        y = B[l-1]
                    else:
                        count = l
                        while count < 2*l - 1:
                            i = next(pqa)
                            G.append(i[5])
                            B.append(i[4])
                            count = count + 1

                        x = G[count]
                        y = B[count]
                else:
                    if N == 1:
                        x = G[l-1]
                        y = B[l-1]
                    else:
                        return []

                return [(x, y)]

            else:

                fs = []
                sol = []
                div = divisors(N)

                for d in div:
                    if divisible(N, d**2):
                        fs.append(d)

                for f in fs:
                    m = N // f**2
                    zs = sqrt_mod(D, abs(m), True)

                    zs = [i for i in zs if i <= abs(m) // 2 ]
                    if abs(m) != 2:
                        zs = zs + [-i for i in zs]
                        if S.Zero in zs:
                            zs.remove(S.Zero) # Remove duplicate zero

                    for z in zs:

                        pqa = PQa(z, abs(m), D)
                        l = 0
                        G = []
                        B = []

                        for i in pqa:

                            a = i[2]
                            G.append(i[5])
                            B.append(i[4])

                            if l != 0 and abs(i[1]) == 1:
                                r = G[l-1]
                                s = B[l-1]

                                if r**2 - D*s**2 == m:
                                    sol.append((f*r, f*s))

                                elif diop_DN(D, -1) != []:
                                    a = diop_DN(D, -1)
                                    sol.append((f*(r*a[0][0] + a[0][1]*s*D), f*(r*a[0][1] + s*a[0][0])))

                                break

                            l = l + 1
                            if l == length(z, abs(m), D):
                                break

                return sol


def cornacchia(a, b, m):
    """
    Solves `ax^2 + by^2 = m` where `\gcd(a, b) = 1 = gcd(a, m)` and `a, b > 0`.

    Uses the algorithm due to Cornacchia. The method only finds primitive
    solutions, i.e. ones with `\gcd(x, y) = 1`. So this method can't be used to
    find the solutions of `x^2 + y^2 = 20` since the only solution to former is
    `(x,y) = (4, 2)` and it is not primitive. When ` a = b = 1`, only the
    solutions with `x \geq y` are found. For more details, see the References.

    Examples
    ========

    >>> from sympy.solvers.diophantine import cornacchia
    >>> cornacchia(2, 3, 35) # equation 2x**2 + 3y**2 = 35
    set([(2, 3), (4, 1)])
    >>> cornacchia(1, 1, 25) # equation x**2 + y**2 = 25
    set([(4, 3)])

    References
    ===========

    .. [1] A. Nitaj, "L'algorithme de Cornacchia"
    .. [2] Solving the diophantine equation ax**2 + by**2 = m by Cornacchia's
        method, [online], Available:
        http://www.numbertheory.org/php/cornacchia.html
    """
    sols = set([])

    a1 = igcdex(a, m)[0]
    v = sqrt_mod(-b*a1, m, True)

    if v is None:
        return None

    if not isinstance(v, list):
        v = [v]

    for t in v:
        if t < m // 2:
            continue

        u, r = t, m

        while True:
            u, r = r, u % r
            if a*r**2 < m:
                break

        m1 = m - a*r**2

        if m1 % b == 0:
            m1 = m1 // b
            if isinstance(sqrt(m1), Integer):
                s = sqrt(m1)
                sols.add((int(r), int(s)))

    return sols


def PQa(P_0, Q_0, D):
    """
    Returns useful information needed to solve the Pell equation.

    There are six sequences of integers defined related to the continued
    fraction representation of `\\frac{P + \sqrt{D}}{Q}`, namely {`P_{i}`},
    {`Q_{i}`}, {`a_{i}`},{`A_{i}`}, {`B_{i}`}, {`G_{i}`}. ``PQa()`` Returns
    these values as a 6-tuple in the same order as mentioned above. Refer [1]_
    for more detailed information.

    Usage
    =====

    ``PQa(P_0, Q_0, D)``: ``P_0``, ``Q_0`` and ``D`` are integers corresponding
    to `P_{0}`, `Q_{0}` and `D` in the continued fraction
    `\\frac{P_{0} + \sqrt{D}}{Q_{0}}`.
    Also it's assumed that `P_{0}^2 == D mod(|Q_{0}|)` and `D` is square free.

    Examples
    ========

    >>> from sympy.solvers.diophantine import PQa
    >>> pqa = PQa(13, 4, 5) # (13 + sqrt(5))/4
    >>> next(pqa) # (P_0, Q_0, a_0, A_0, B_0, G_0)
    (13, 4, 3, 3, 1, -1)
    >>> next(pqa) # (P_1, Q_1, a_1, A_1, B_1, G_1)
    (-1, 1, 1, 4, 1, 3)

    References
    ==========

    .. [1] Solving the generalized Pell equation x^2 - Dy^2 = N, John P.
        Robertson, July 31, 2004, Pages 4 - 8. http://www.jpr2718.org/pell.pdf
    """
    A_i_2 = B_i_1 = 0
    A_i_1 = B_i_2 = 1

    G_i_2 = -P_0
    G_i_1 = Q_0

    P_i = P_0
    Q_i = Q_0

    while(1):

        a_i = floor((P_i + sqrt(D))/Q_i)
        A_i = a_i*A_i_1 + A_i_2
        B_i = a_i*B_i_1 + B_i_2
        G_i = a_i*G_i_1 + G_i_2

        yield P_i, Q_i, a_i, A_i, B_i, G_i

        A_i_1, A_i_2 = A_i, A_i_1
        B_i_1, B_i_2 = B_i, B_i_1
        G_i_1, G_i_2 = G_i, G_i_1

        P_i = a_i*Q_i - P_i
        Q_i = (D - P_i**2)/Q_i


def diop_bf_DN(D, N, t=symbols("t", integer=True)):
    """
    Uses brute force to solve the equation, `x^2 - Dy^2 = N`.

    Mainly concerned with the generalized Pell equation which is the case when
    `D > 0, D` is not a perfect square. For more information on the case refer
    [1]_. Let `(t, u)` be the minimal positive solution of the equation
    `x^2 - Dy^2 = 1`. Then this method requires
    `\sqrt{\\frac{\mid N \mid (t \pm 1)}{2D}}` to be small.

    Usage
    =====

    ``diop_bf_DN(D, N, t)``: ``D`` and ``N`` are coefficients in
    `x^2 - Dy^2 = N` and ``t`` is the parameter to be used in the solutions.

    Details
    =======

    ``D`` and ``N`` correspond to D and N in the equation.
    ``t`` is the parameter to be used in the solutions.

    Examples
    ========

    >>> from sympy.solvers.diophantine import diop_bf_DN
    >>> diop_bf_DN(13, -4)
    [(3, 1), (-3, 1), (36, 10)]
    >>> diop_bf_DN(986, 1)
    [(49299, 1570)]

    See Also
    ========

    diop_DN()

    References
    ==========

    .. [1] Solving the generalized Pell equation x**2 - D*y**2 = N, John P.
        Robertson, July 31, 2004, Page 15. http://www.jpr2718.org/pell.pdf
    """
    sol = []
    a = diop_DN(D, 1)
    u = a[0][0]
    v = a[0][1]


    if abs(N) == 1:
        return diop_DN(D, N)

    elif N > 1:
        L1 = 0
        L2 = floor(sqrt(S(N*(u - 1))/(2*D))) + 1

    elif N < -1:
        L1 = ceiling(sqrt(S(-N)/D))
        L2 = floor(sqrt(S(-N*(u + 1))/(2*D))) + 1

    else:
        if D < 0:
            return [(S.Zero, S.Zero)]
        elif D == 0:
            return [(S.Zero, t)]
        else:
            if isinstance(sqrt(D), Integer):
                return [(sqrt(D)*t, t), (-sqrt(D)*t, t)]
            else:
                return [(S.Zero, S.Zero)]


    for y in range(L1, L2):
        if isinstance(sqrt(N + D*y**2), Integer):
            x = sqrt(N + D*y**2)
            sol.append((x, y))
            if not equivalent(x, y, -x, y, D, N):
                sol.append((-x, y))

    return sol


def equivalent(u, v, r, s, D, N):
    """
    Returns True if two solutions `(u, v)` and `(r, s)` of `x^2 - Dy^2 = N`
    belongs to the same equivalence class and False otherwise.

    Two solutions `(u, v)` and `(r, s)` to the above equation fall to the same
    equivalence class iff both `(ur - Dvs)` and `(us - vr)` are divisible by
    `N`. See reference [1]_. No test is performed to test whether `(u, v)` and
    `(r, s)` are actually solutions to the equation. User should take care of
    this.

    Usage
    =====

    ``equivalent(u, v, r, s, D, N)``: `(u, v)` and `(r, s)` are two solutions
    of the equation `x^2 - Dy^2 = N` and all parameters involved are integers.

    Examples
    ========

    >>> from sympy.solvers.diophantine import equivalent
    >>> equivalent(18, 5, -18, -5, 13, -1)
    True
    >>> equivalent(3, 1, -18, 393, 109, -4)
    False

    References
    ==========

    .. [1] Solving the generalized Pell equation x**2 - D*y**2 = N, John P.
        Robertson, July 31, 2004, Page 12. http://www.jpr2718.org/pell.pdf

    """
    return divisible(u*r - D*v*s, N) and divisible(u*s - v*r, N)


def length(P, Q, D):
    """
    Returns the (length of aperiodic part + length of periodic part) of
    continued fraction representation of `\\frac{P + \sqrt{D}}{Q}`.

    It is important to remember that this does NOT return the length of the
    periodic part but the addition of the legths of the two parts as mentioned
    above.

    Usage
    =====

    ``length(P, Q, D)``: ``P``, ``Q`` and ``D`` are integers corresponding to
    the continued fraction `\\frac{P + \sqrt{D}}{Q}`.

    Details
    =======

    ``P``, ``D`` and ``Q`` corresponds to P, D and Q in the continued fraction,
    `\\frac{P + \sqrt{D}}{Q}`.

    Examples
    ========

    >>> from sympy.solvers.diophantine import length
    >>> length(-2 , 4, 5) # (-2 + sqrt(5))/4
    3
    >>> length(-5, 4, 17) # (-5 + sqrt(17))/4
    4
    """
    x = P + sqrt(D)
    y = Q

    x = sympify(x)
    v, res = [], []
    q = x/y

    if q < 0:
        v.append(q)
        res.append(floor(q))
        q = q - floor(q)
        num, den = rad_rationalize(1, q)
        q = num / den

    while 1:
        v.append(q)
        a = int(q)
        res.append(a)

        if q == a:
            return len(res)

        num, den = rad_rationalize(1,(q - a))
        q = num / den

        if q in v:
            return len(res)


def transformation_to_DN(eq):
    """
    This function transforms general quadratic,
    `ax^2 + bxy + cy^2 + dx + ey + f = 0`
    to more easy to deal with `X^2 - DY^2 = N` form.

    This is used to solve the general quadratic equation by transforming it to
    the latter form. Refer [1]_ for more detailed information on the
    transformation. This function returns a tuple (A, B) where A is a 2 X 2
    matrix and B is a 2 X 1 matrix such that,

    Transpose([x y]) =  A * Transpose([X Y]) + B

    Usage
    =====

    ``transformation_to_DN(eq)``: where ``eq`` is the quadratic to be
    transformed.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.solvers.diophantine import transformation_to_DN
    >>> from sympy.solvers.diophantine import classify_diop
    >>> A, B = transformation_to_DN(x**2 - 3*x*y - y**2 - 2*y + 1)
    >>> A
    Matrix([
    [1/26, 3/26],
    [   0, 1/13]])
    >>> B
    Matrix([
    [-6/13],
    [-4/13]])

    A, B  returned are such that Transpose((x y)) =  A * Transpose((X Y)) + B.
    Substituting these values for `x` and `y` and a bit of simplifying work
    will give an equation of the form `x^2 - Dy^2 = N`.

    >>> from sympy.abc import X, Y
    >>> from sympy import Matrix, simplify, Subs
    >>> u = (A*Matrix([X, Y]) + B)[0] # Transformation for x
    >>> u
    X/26 + 3*Y/26 - 6/13
    >>> v = (A*Matrix([X, Y]) + B)[1] # Transformation for y
    >>> v
    Y/13 - 4/13

    Next we will substitute these formulas for `x` and `y` and do
    ``simplify()``.

    >>> eq = simplify(Subs(x**2 - 3*x*y - y**2 - 2*y + 1, (x, y), (u, v)).doit())
    >>> eq
    X**2/676 - Y**2/52 + 17/13

    By multiplying the denominator appropriately, we can get a Pell equation
    in the standard form.

    >>> eq * 676
    X**2 - 13*Y**2 + 884

    If only the final equation is needed, ``find_DN()`` can be used.

    See Also
    ========

    find_DN()

    References
    ==========

    .. [1] Solving the equation ax^2 + bxy + cy^2 + dx + ey + f = 0,
           John P.Robertson, May 8, 2003, Page 7 - 11.
           http://www.jpr2718.org/ax2p.pdf
    """


    var, coeff, diop_type = classify_diop(eq)
    if diop_type == "binary_quadratic":
        return _transformation_to_DN(var, coeff)


def _transformation_to_DN(var, coeff):

    x, y = var[:2]

    a = coeff[x**2]
    b = coeff[x*y]
    c = coeff[y**2]
    d = coeff[x]
    e = coeff[y]
    f = coeff[Integer(1)]

    g = igcd(a, igcd(b, igcd(c, igcd(d, igcd(e, f)))))
    a = a // g
    b = b // g
    c = c // g
    d = d // g
    e = e // g
    f = f // g

    X, Y = symbols("X, Y", integer=True)

    if b != Integer(0):
        B = (S(2*a)/b).p
        C = (S(2*a)/b).q
        A = (S(a)/B**2).p
        T = (S(a)/B**2).q

        # eq_1 = A*B*X**2 + B*(c*T - A*C**2)*Y**2 + d*T*X + (B*e*T - d*T*C)*Y + f*T*B
        coeff = {X**2: A*B, X*Y: 0, Y**2: B*(c*T - A*C**2), X: d*T, Y: B*e*T - d*T*C, Integer(1): f*T*B}
        A_0, B_0 = _transformation_to_DN([X, Y], coeff)
        return Matrix(2, 2, [S(1)/B, -S(C)/B, 0, 1])*A_0, Matrix(2, 2, [S(1)/B, -S(C)/B, 0, 1])*B_0

    else:
        if d != Integer(0):
            B = (S(2*a)/d).p
            C = (S(2*a)/d).q
            A = (S(a)/B**2).p
            T = (S(a)/B**2).q

            # eq_2 = A*X**2 + c*T*Y**2 + e*T*Y + f*T - A*C**2
            coeff = {X**2: A, X*Y: 0, Y**2: c*T, X: 0, Y: e*T, Integer(1): f*T - A*C**2}
            A_0, B_0 = _transformation_to_DN([X, Y], coeff)
            return Matrix(2, 2, [S(1)/B, 0, 0, 1])*A_0, Matrix(2, 2, [S(1)/B, 0, 0, 1])*B_0 + Matrix([-S(C)/B, 0])

        else:
            if e != Integer(0):
                B = (S(2*c)/e).p
                C = (S(2*c)/e).q
                A = (S(c)/B**2).p
                T = (S(c)/B**2).q

                # eq_3 = a*T*X**2 + A*Y**2 + f*T - A*C**2
                coeff = {X**2: a*T, X*Y: 0, Y**2: A, X: 0, Y: 0, Integer(1): f*T - A*C**2}
                A_0, B_0 = _transformation_to_DN([X, Y], coeff)
                return Matrix(2, 2, [1, 0, 0, S(1)/B])*A_0, Matrix(2, 2, [1, 0, 0, S(1)/B])*B_0 + Matrix([0, -S(C)/B])

            else:
                # TODO: pre-simplification: Not necessary but may simplify
                # the equation.

                return Matrix(2, 2, [S(1)/a, 0, 0, 1]), Matrix([0, 0])


def find_DN(eq):
    """
    This function returns a tuple, `(D, N)` of the simplified form,
    `x^2 - Dy^2 = N`, corresponding to the general quadratic,
    `ax^2 + bxy + cy^2 + dx + ey + f = 0`.

    Solving the general quadratic is then equivalent to solving the equation
    `X^2 - DY^2 = N` and transforming the solutions by using the transformation
    matrices returned by ``transformation_to_DN()``.

    Usage
    =====

    ``find_DN(eq)``: where ``eq`` is the quadratic to be transformed.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.solvers.diophantine import find_DN
    >>> find_DN(x**2 - 3*x*y - y**2 - 2*y + 1)
    (13, -884)

    Interpretation of the output is that we get `X^2 -13Y^2 = -884` after
    transforming `x^2 - 3xy - y^2 - 2y + 1` using the transformation returned
    by ``transformation_to_DN()``.

    See Also
    ========

    transformation_to_DN()

    References
    ==========

    .. [1] Solving the equation ax^2 + bxy + cy^2 + dx + ey + f = 0,
           John P.Robertson, May 8, 2003, Page 7 - 11.
           http://www.jpr2718.org/ax2p.pdf
    """
    var, coeff, diop_type = classify_diop(eq)
    if diop_type == "binary_quadratic":
        return _find_DN(var, coeff)


def _find_DN(var, coeff):

    x, y = var[:2]
    X, Y = symbols("X, Y", integer=True)
    A , B = _transformation_to_DN(var, coeff)

    u = (A*Matrix([X, Y]) + B)[0]
    v = (A*Matrix([X, Y]) + B)[1]
    eq = x**2*coeff[x**2] + x*y*coeff[x*y] + y**2*coeff[y**2] + x*coeff[x] + y*coeff[y] + coeff[Integer(1)]

    simplified = _mexpand(Subs(eq, (x, y), (u, v)).doit())

    coeff = dict([reversed(t.as_independent(*[X, Y])) for t in simplified.args])

    for term in [X**2, Y**2, Integer(1)]:
        if term not in coeff.keys():
            coeff[term] = Integer(0)

    return -coeff[Y**2]/coeff[X**2], -coeff[Integer(1)]/coeff[X**2]


def check_param(x, y, a, t):
    """
    Check if there is a number modulo ``a`` such that ``x`` and ``y`` are both
    integers. If exist, then find a parametric representation for ``x`` and
    ``y``.

    Here ``x`` and ``y`` are functions of ``t``.
    """
    k, m, n = symbols("k, m, n", integer=True)
    p = Wild("p", exclude=[k])
    q = Wild("q", exclude=[k])
    ok = False

    for i in range(a):

        z_x = _mexpand(Subs(x, t, a*k + i).doit()).match(p*k + q)
        z_y = _mexpand(Subs(y, t, a*k + i).doit()).match(p*k + q)

        if (isinstance(z_x[p], Integer) and isinstance(z_x[q], Integer) and
            isinstance(z_y[p], Integer) and isinstance(z_y[q], Integer)):
            ok = True
            break

    if ok == True:

        x_param = x.match(p*t + q)
        y_param = y.match(p*t + q)

        if x_param[p] == 0 or y_param[p] == 0:
            if x_param[p] == 0:
                l1, junk = Poly(y).clear_denoms()
            else:
                l1 = 1

            if y_param[p] == 0:
                l2, junk = Poly(x).clear_denoms()
            else:
                l2 = 1

            return x*ilcm(l1, l2), y*ilcm(l1, l2)

        eq = S(m - x_param[q])/x_param[p] - S(n - y_param[q])/y_param[p]

        lcm_denom, junk = Poly(eq).clear_denoms()
        eq = eq * lcm_denom

        return diop_solve(eq, t)[0], diop_solve(eq, t)[1]
    else:
        return (None, None)


def diop_ternary_quadratic(eq):
    """
    Solves the general quadratic ternary form,
    `ax^2 + by^2 + cz^2 + fxy + gyz + hxz = 0`.

    Returns a tuple `(x, y, z)` which is a base solution for the above
    equation. If there are no solutions, `(None, None, None)` is returned.

    Usage
    =====

    ``diop_ternary_quadratic(eq)``: Return a tuple containing a basic solution
    to ``eq``.

    Details
    =======

    ``eq`` should be an homogeneous expression of degree two in three variables
    and it is assumed to be zero.

    Examples
    ========

    >>> from sympy.abc import x, y, z
    >>> from sympy.solvers.diophantine import diop_ternary_quadratic
    >>> diop_ternary_quadratic(x**2 + 3*y**2 - z**2)
    (1, 0, 1)
    >>> diop_ternary_quadratic(4*x**2 + 5*y**2 - z**2)
    (1, 0, 2)
    >>> diop_ternary_quadratic(45*x**2 - 7*y**2 - 8*x*y - z**2)
    (28, 45, 105)
    >>> diop_ternary_quadratic(x**2 - 49*y**2 - z**2 + 13*z*y -8*x*y)
    (9, 1, 5)
    """
    var, coeff, diop_type = classify_diop(eq)

    if diop_type == "homogeneous_ternary_quadratic":
        return _diop_ternary_quadratic(var, coeff)


def _diop_ternary_quadratic(_var, coeff):

    x, y, z = _var[:3]

    var = [x]*3
    var[0], var[1], var[2] = _var[0], _var[1], _var[2]

    # Equations of the form B*x*y + C*z*x + E*y*z = 0 and At least two of the
    # coefficients A, B, C are non-zero.
    # There are infinitely many solutions for the equation.
    # Ex: (0, 0, t), (0, t, 0), (t, 0, 0)
    # Equation can be re-written as y*(B*x + E*z) = -C*x*z and we can find rather
    # unobviuos solutions. Set y = -C and B*x + E*z = x*z. The latter can be solved by
    # using methods for binary quadratic diophantine equations. Let's select the
    # solution which minimizes |x| + |z|

    if coeff[x**2] == 0 and coeff[y**2] == 0 and coeff[z**2] == 0:
        if coeff[x*z] != 0:
            sols = diophantine(coeff[x*y]*x + coeff[y*z]*z - x*z)
            s = sols.pop()
            min_sum = abs(s[0]) + abs(s[1])

            for r in sols:
                if abs(r[0]) + abs(r[1]) < min_sum:
                    s = r
                    min_sum = abs(s[0]) + abs(s[1])

                x_0, y_0, z_0 = s[0], -coeff[x*z], s[1]

        else:
            var[0], var[1] = _var[1], _var[0]
            y_0, x_0, z_0 = _diop_ternary_quadratic(var, coeff)

        return simplified(x_0, y_0, z_0)

    if coeff[x**2] == 0:
        # If the coefficient of x is zero change the variables
        if coeff[y**2] == 0:
            var[0], var[2] = _var[2], _var[0]
            z_0, y_0, x_0 = _diop_ternary_quadratic(var, coeff)

        else:
            var[0], var[1] = _var[1], _var[0]
            y_0, x_0, z_0 = _diop_ternary_quadratic(var, coeff)

    else:
        if coeff[x*y] != 0 or coeff[x*z] != 0:
        # Apply the transformation x --> X - (B*y + C*z)/(2*A)
            A = coeff[x**2]
            B = coeff[x*y]
            C = coeff[x*z]
            D = coeff[y**2]
            E = coeff[y*z]
            F = coeff[z**2]

            _coeff = dict()

            _coeff[x**2] = 4*A**2
            _coeff[y**2] = 4*A*D - B**2
            _coeff[z**2] = 4*A*F - C**2
            _coeff[y*z] = 4*A*E - 2*B*C
            _coeff[x*y] = 0
            _coeff[x*z] = 0

            X_0, y_0, z_0 = _diop_ternary_quadratic(var, _coeff)

            if X_0 == None:
                return (None, None, None)

            l = (S(B*y_0 + C*z_0)/(2*A)).q
            x_0, y_0, z_0 = X_0*l - (S(B*y_0 + C*z_0)/(2*A)).p, y_0*l, z_0*l

        elif coeff[z*y] != 0:
            if coeff[y**2] == 0:
                if coeff[z**2] == 0:
                    # Equations of the form A*x**2 + E*yz = 0.
                    A = coeff[x**2]
                    E = coeff[y*z]

                    b = (S(-E)/A).p
                    a = (S(-E)/A).q

                    x_0, y_0, z_0 = b, a, b

                else:
                    # Ax**2 + E*y*z + F*z**2  = 0
                    var[0], var[2] = _var[2], _var[0]
                    z_0, y_0, x_0 = _diop_ternary_quadratic(var, coeff)

            else:
                # A*x**2 + D*y**2 + E*y*z + F*z**2 = 0, C may be zero
                var[0], var[1] = _var[1], _var[0]
                y_0, x_0, z_0 = _diop_ternary_quadratic(var, coeff)

        else:
            # Ax**2 + D*y**2 + F*z**2 = 0, C may be zero
            x_0, y_0, z_0 = _diop_ternary_quadratic_normal(var, coeff)

    return simplified(x_0, y_0, z_0)


def transformation_to_normal(eq):
    """
    Returns the transformation Matrix from general ternary quadratic equation
    `eq` to normal form.

    General form of the ternary quadratic equation is `ax^2 + by^2 cz^2 + dxy +
    eyz + fxz`. This function returns a 3X3 transformation Matrix which
    transforms the former equation to the form `ax^2 + by^2 + cz^2 = 0`. This
    is not used in solving ternary quadratics. Only implemented for the sake
    of completeness.
    """
    var, coeff, diop_type = classify_diop(eq)

    if diop_type == "homogeneous_ternary_quadratic":
        return _transformation_to_normal(var, coeff)


def _transformation_to_normal(var, coeff):

    _var = [var[0]]*3
    _var[1], _var[2] = var[1], var[2]

    x, y, z = var[:3]

    if coeff[x**2] == 0:
        # If the coefficient of x is zero change the variables
        if coeff[y**2] == 0:
            _var[0], _var[2] = var[2], var[0]
            T = _transformation_to_normal(_var, coeff)
            T.row_swap(0, 2)
            T.col_swap(0, 2)
            return T

        else:
            _var[0], _var[1] = var[1], var[0]
            T = _transformation_to_normal(_var, coeff)
            T.row_swap(0, 1)
            T.col_swap(0, 1)
            return T

    else:
        # Apply the transformation x --> X - (B*Y + C*Z)/(2*A)
        if coeff[x*y] != 0 or coeff[x*z] != 0:
            A = coeff[x**2]
            B = coeff[x*y]
            C = coeff[x*z]
            D = coeff[y**2]
            E = coeff[y*z]
            F = coeff[z**2]

            _coeff = dict()

            _coeff[x**2] = 4*A**2
            _coeff[y**2] = 4*A*D - B**2
            _coeff[z**2] = 4*A*F - C**2
            _coeff[y*z] = 4*A*E - 2*B*C
            _coeff[x*y] = 0
            _coeff[x*z] = 0

            T_0 = _transformation_to_normal(_var, _coeff)
            return Matrix(3, 3, [1, S(-B)/(2*A), S(-C)/(2*A), 0, 1, 0, 0, 0, 1]) * T_0

        elif coeff[y*z] != 0:
            if coeff[y**2] == 0:
                if coeff[z**2] == 0:
                    # Equations of the form A*x**2 + E*yz = 0.
                    # Apply transformation y -> Y + Z ans z -> Y - Z
                    return Matrix(3, 3, [1, 0, 0, 0, 1, 1, 0, 1, -1])

                else:
                    # Ax**2 + E*y*z + F*z**2  = 0
                    _var[0], _var[2] = var[2], var[0]
                    T = _transformtion_to_normal(_var, coeff)
                    T.row_swap(0, 2)
                    T.col_swap(0, 2)
                    return T

            else:
                # A*x**2 + D*y**2 + E*y*z + F*z**2 = 0, F may be zero
                _var[0], _var[1] = var[1], var[0]
                T = _transformation_to_normal(_var, coeff)
                T.row_swap(0, 1)
                T.col_swap(0, 1)
                return T

        else:
            return Matrix(3, 3, [1, 0, 0, 0, 1, 0, 0, 0, 1])


def simplified(x, y, z):
    """
    Simplify the solution `(x, y, z)`.
    """
    if x == None or y == None or z == None:
        return (x, y, z)

    g = igcd(x, igcd(y, z))

    return x // g, y // g, z // g


def parametrize_ternary_quadratic(eq):
    """
    Returns the parametrized general solution for the ternary quadratic
    equation ``eq`` which has the form
    `ax^2 + by^2 + cz^2 + fxy + gyz + hxz = 0`.

    Examples
    ========

    >>> from sympy.abc import x, y, z
    >>> from sympy.solvers.diophantine import parametrize_ternary_quadratic
    >>> parametrize_ternary_quadratic(x**2 + y**2 - z**2)
    (2*p*q, p**2 - q**2, p**2 + q**2)

    Here `p` and `q` are two co-prime integers.

    >>> parametrize_ternary_quadratic(3*x**2 + 2*y**2 - z**2 - 2*x*y + 5*y*z - 7*y*z)
    (2*p**2 - 2*p*q - q**2, 2*p**2 + 2*p*q - q**2, 2*p**2 - 2*p*q + 3*q**2)
    >>> parametrize_ternary_quadratic(124*x**2 - 30*y**2 - 7729*z**2)
    (-1410*p**2 - 363263*q**2, 2700*p**2 + 30916*p*q - 695610*q**2, -60*p**2 + 5400*p*q + 15458*q**2)

    References
    ==========

    .. [1] The algorithmic resolution of Diophantine equations, Nigel P. Smart,
           London Mathematical Society Student Texts 41, Cambridge University
           Press, Cambridge, 1998.

    """
    var, coeff, diop_type = classify_diop(eq)

    if diop_type == "homogeneous_ternary_quadratic":
        x_0, y_0, z_0 = _diop_ternary_quadratic(var, coeff)
        return _parametrize_ternary_quadratic((x_0, y_0, z_0), var, coeff)


def _parametrize_ternary_quadratic(solution, _var, coeff):

    x, y, z = _var[:3]

    x_0, y_0, z_0 = solution[:3]

    v = [x]*3
    v[0], v[1], v[2] = _var[0], _var[1], _var[2]

    if x_0 == None:
        return (None, None, None)

    if x_0 == 0:
        if y_0 == 0:
            v[0], v[2] = v[2], v[0]
            z_p, y_p, x_p = _parametrize_ternary_quadratic((z_0, y_0, x_0), v, coeff)
            return x_p, y_p, z_p
        else:
            v[0], v[1] = v[1], v[0]
            y_p, x_p, z_p = _parametrize_ternary_quadratic((y_0, x_0, z_0), v, coeff)
            return x_p, y_p, z_p

    x, y, z = v[:3]
    r, p, q = symbols("r, p, q", integer=True)

    eq = x**2*coeff[x**2] + y**2*coeff[y**2] + z**2*coeff[z**2] + x*y*coeff[x*y] + y*z*coeff[y*z] + z*x*coeff[z*x]
    eq_1 = Subs(eq, (x, y, z), (r*x_0, r*y_0 + p, r*z_0 + q)).doit()
    eq_1 = _mexpand(eq_1)
    A, B = eq_1.as_independent(r, as_Add=True)


    x = A*x_0
    y = (A*y_0 - _mexpand(B/r*p))
    z = (A*z_0 - _mexpand(B/r*q))

    return x, y, z


def diop_ternary_quadratic_normal(eq):
    """
    Solves the quadratic ternary diophantine equation,
    `ax^2 + by^2 + cz^2 = 0`.

    Here the coefficients `a`, `b`, and `c` should be non zero. Otherwise the
    equation will be a quadratic binary or univariate equation. If solvable,
    returns a tuple `(x, y, z)` that satisifes the given equation. If the
    equation does not have integer solutions, `(None, None, None)` is returned.

    Usage
    =====

    ``diop_ternary_quadratic_normal(eq)``: where ``eq`` is an equation of the form
    `ax^2 + by^2 + cz^2 = 0`.

    Examples
    ========

    >>> from sympy.abc import x, y, z
    >>> from sympy.solvers.diophantine import diop_ternary_quadratic_normal
    >>> diop_ternary_quadratic_normal(x**2 + 3*y**2 - z**2)
    (1, 0, 1)
    >>> diop_ternary_quadratic_normal(4*x**2 + 5*y**2 - z**2)
    (1, 0, 2)
    >>> diop_ternary_quadratic_normal(34*x**2 - 3*y**2 - 301*z**2)
    (4, 9, 1)
    """
    var, coeff, diop_type = classify_diop(eq)

    if diop_type == "homogeneous_ternary_quadratic":
        return _diop_ternary_quadratic_normal(var, coeff)


def _diop_ternary_quadratic_normal(var, coeff):

    x, y, z = var[:3]

    a = coeff[x**2]
    b = coeff[y**2]
    c = coeff[z**2]

    if a*b*c == 0:
        raise ValueError("Try factoring out you equation or using diophantine()")

    g = igcd(a, igcd(b, c))

    a = a // g
    b = b // g
    c = c // g

    a_0 = square_factor(a)
    b_0 = square_factor(b)
    c_0 = square_factor(c)

    a_1 = a // a_0**2
    b_1 = b // b_0**2
    c_1 = c // c_0**2

    a_2, b_2, c_2 = pairwise_prime(a_1, b_1, c_1)

    A = -a_2*c_2
    B = -b_2*c_2

    # If following two conditions are satisified then there are no solutions
    if A < 0 and B < 0:
        return (None, None, None)

    if (sqrt_mod(-b_2*c_2, a_2) == None or sqrt_mod(-c_2*a_2, b_2) == None or
        sqrt_mod(-a_2*b_2, c_2) == None):
        return (None, None, None)

    z_0, x_0, y_0 = descent(A, B)

    if divisible(z_0, c_2) == True:
        z_0 = z_0 // abs(c_2)
    else:
        x_0 = x_0*(S(z_0)/c_2).q
        y_0 = y_0*(S(z_0)/c_2).q
        z_0 = (S(z_0)/c_2).p

    x_0, y_0, z_0 = simplified(x_0, y_0, z_0)

    # Holzer reduction
    if sign(a) == sign(b):
        x_0, y_0, z_0 = holzer(x_0, y_0, z_0, abs(a_2), abs(b_2), abs(c_2))
    elif sign(a) == sign(c):
        x_0, z_0, y_0 = holzer(x_0, z_0, y_0, abs(a_2), abs(c_2), abs(b_2))
    else:
        y_0, z_0, x_0 = holzer(y_0, z_0, x_0, abs(b_2), abs(c_2), abs(a_2))

    x_0 = reconstruct(b_1, c_1, x_0)
    y_0 = reconstruct(a_1, c_1, y_0)
    z_0 = reconstruct(a_1, b_1, z_0)

    l = ilcm(a_0, ilcm(b_0, c_0))

    x_0 = abs(x_0*l//a_0)
    y_0 = abs(y_0*l//b_0)
    z_0 = abs(z_0*l//c_0)

    return simplified(x_0, y_0, z_0)


def square_factor(a):
    """
    Returns an integer `c` s.t. `a = c^2k, \ c,k \in Z`. Here `k` is square
    free.

    Examples
    ========

    >>> from sympy.solvers.diophantine import square_factor
    >>> square_factor(24)
    2
    >>> square_factor(36)
    6
    >>> square_factor(1)
    1
    """
    f = factorint(abs(a))
    c = 1

    for p, e in f.items():
        c = c * p**(e//2)

    return c


def pairwise_prime(a, b, c):
    """
    Transform `ax^2 + by^2 + cz^2 = 0` into an equivalent equation
    `a'x^2 + b'y^2 + c'z^2 = 0` where `a', b', c'` are pairwise relatively
    prime.

    Returns a tuple containing `a', b', c'`. `\gcd(a, b, c)` should equal `1`
    for this to work. The solutions for `ax^2 + by^2 + cz^2 = 0` can be
    recovered from the solutions of `a'x^2 + b'y^2 + c'z^2 = 0`.

    Examples
    ========

    >>> from sympy.solvers.diophantine import pairwise_prime
    >>> pairwise_prime(6, 15, 10)
    (5, 2, 3)

    See Also
    ========

    make_prime(), reocnstruct()
    """
    a, b, c = make_prime(a, b, c)
    b, c, a = make_prime(b, c, a)
    c, a, b = make_prime(c, a, b)

    return a, b, c


def make_prime(a, b, c):
    """
    Transform the equation `ax^2 + by^2 + cz^2 = 0` to an equivalent equation
    `a'x^2 + b'y^2 + c'z^2 = 0` with `\gcd(a', b') = 1`.

    Returns a tuple `(a', b', c')` which satisfies above conditions. Note that
    in the returned tuple `\gcd(a', c')` and `\gcd(b', c')` can take any value.

    Examples
    ========

    >>> from sympy.solvers.diophantine import make_prime
    >>> make_prime(4, 2, 7)
    (2, 1, 14)

    See Also
    ========

    pairwaise_prime(), reconstruct()
    """
    g = igcd(a, b)

    if g != 1:
        f = factorint(g)
        for p, e in f.items():
            a = a // p**e
            b = b // p**e

            if e % 2 == 1:
                c = p*c

    return a, b, c


def reconstruct(a, b, z):
    """
    Reconstruct the `z` value of an equivalent solution of `ax^2 + by^2 + cz^2`
    from the `z` value of a solution of a transformed version of the above
    equation.
    """
    g = igcd(a, b)

    if g != 1:
        f = factorint(g)
        for p, e in f.items():
            if e %2 == 0:
                z = z*p**(e//2)
            else:
                z = z*p**((e//2)+1)

    return z


def ldescent(A, B):
    """
    Uses Lagrange's method to find a non trivial solution to
    `w^2 = Ax^2 + By^2`.

    Here, `A \\neq 0` and `B \\neq 0` and `A` and `B` are square free. Output a
    tuple `(w_0, x_0, y_0)` which is a solution to the above equation.

    Examples
    ========

    >>> from sympy.solvers.diophantine import ldescent
    >>> ldescent(1, 1) # w^2 = x^2 + y^2
    (1, 1, 0)
    >>> ldescent(4, -7) # w^2 = 4x^2 - 7y^2
    (2, -1, 0)

    This means that `x = -1, y = 0` and `w = 2` is a solution to the equation
    `w^2 = 4x^2 - 7y^2`

    >>> ldescent(5, -1) # w^2 = 5x^2 - y^2
    (2, 1, -1)

    References
    ==========

    .. [1] The algorithmic resolution of Diophantine equations, Nigel P. Smart,
           London Mathematical Society Student Texts 41, Cambridge University
           Press, Cambridge, 1998.
    .. [2] Efficient Solution of Rational Conices, J. E. Cremona and D. Rusin,
           Mathematics of Computation, Volume 00, Number 0.
    """
    if abs(A) > abs(B):
        w, y, x = ldescent(B, A)
        return w, x, y

    if A == 1:
        return (S.One, S.One, 0)

    if B == 1:
        return (S.One, 0, S.One)

    r = sqrt_mod(A, B)

    Q = (r**2 - A) // B

    if Q == 0:
        B_0 = 1
        d = 0
    else:
        div = divisors(Q)
        B_0 = None

        for i in div:
            if isinstance(sqrt(abs(Q) // i), Integer):
                B_0, d = sign(Q)*i, sqrt(abs(Q) // i)
                break

    if B_0 != None:
        W, X, Y = ldescent(A, B_0)
        return simplified((-A*X + r*W), (r*X - W), Y*(B_0*d))
    # In this module Descent will always be called with inputs which have solutions.


def descent(A, B):
    """
    Lagrange's `descent()` with lattice-reduction to find solutions to
    `x^2 = Ay^2 + Bz^2`.

    Here `A` and `B` should be square free and pairwise prime. Always should be
    called with suitable ``A`` and ``B`` so that the above equation has
    solutions.

    This is more faster than the normal Lagrange's descent algorithm because
    the gaussian reduction is used.

    Examples
    ========

    >>> from sympy.solvers.diophantine import descent
    >>> descent(3, 1) # x**2 = 3*y**2 + z**2
    (1, 0, 1)

    `(x, y, z) = (1, 0, 1)` is a solution to the above equation.

    >>> descent(41, -113)
    (-16, -3, 1)

    References
    ==========

    .. [1] Efficient Solution of Rational Conices, J. E. Cremona and D. Rusin,
           Mathematics of Computation, Volume 00, Number 0.
    """
    if abs(A) > abs(B):
        x, y, z = descent(B, A)
        return x, z, y

    if B == 1:
        return (1, 0, 1)
    if A == 1:
        return (1, 1, 0)
    if B == -1:
        return (None, None, None)
    if B == -A:
        return (0, 1, 1)
    if B == A:
        x, z, y = descent(-1, A)
        return (A*y, z, x)

    w = sqrt_mod(A, B)
    x_0, z_0 = gaussian_reduce(w, A, B)

    t = (x_0**2 - A*z_0**2) // B
    t_2 = square_factor(t)
    t_1 = t // t_2**2

    x_1, z_1, y_1 = descent(A, t_1)

    return simplified(x_0*x_1 + A*z_0*z_1, z_0*x_1 + x_0*z_1, t_1*t_2*y_1)


def gaussian_reduce(w, a, b):
    """
    Returns a reduced solution `(x, z)` to the congruence
    `X^2 - aZ^2 \equiv 0 \ (mod \ b)` so that `x^2 + |a|z^2` is minimal.

    Details
    =======

    Here ``w`` is a solution of the congruence `x^2 \equiv a \ (mod \ b)`

    References
    ==========

    .. [1] Gaussian lattice Reduction [online]. Available:
        http://home.ie.cuhk.edu.hk/~wkshum/wordpress/?p=404
    .. [2] Efficient Solution of Rational Conices, J. E. Cremona and D. Rusin,
        Mathematics of Computation, Volume 00, Number 0.
    """
    u = (0, 1)
    v = (1, 0)

    if dot(u, v, w, a, b) < 0:
        v = (-v[0], -v[1])

    if norm(u, w, a, b) < norm(v, w, a, b):
        u, v = v, u

    while norm(u, w, a, b) > norm(v, w, a, b):
        k = dot(u, v, w, a, b) // dot(v, v, w, a, b)
        u, v = v, (u[0]- k*v[0], u[1]- k*v[1])

    u, v = v, u

    if dot(u, v, w, a, b) < dot(v, v, w, a, b)/2 or norm((u[0]-v[0], u[1]-v[1]), w, a, b) > norm(v, w, a, b):
        c = v
    else:
        c = (u[0] - v[0], u[1] - v[1])

    return c[0]*w + b*c[1], c[0]


def dot(u, v, w, a, b):
    """
    Returns a special dot product of the vectors `u = (u_{1}, u_{2})` and
    `v = (v_{1}, v_{2})` which is defined in order to reduce solution of
    the congruence equation `X^2 - aZ^2 \equiv 0 \ (mod \ b)`.
    """
    u_1, u_2 = u[:2]
    v_1, v_2 = v[:2]
    return (w*u_1 + b*u_2)*(w*v_1 + b*v_2) + abs(a)*u_1*v_1


def norm(u, w, a, b):
    """
    Returns the norm of the vector `u = (u_{1}, u_{2})` under the dot product
    defined by `u \cdot v = (wu_{1} + bu_{2})(w*v_{1} + bv_{2}) + |a|*u_{1}*v_{1}`
    where `u = (u_{1}, u_{2})` and `v = (v_{1}, v_{2})`.
    """
    u_1, u_2 = u[:2]
    return sqrt(dot((u_1, u_2), (u_1, u_2), w, a, b))


def holzer(x_0, y_0, z_0, a, b, c):
    """
    Simplify the solution `(x_{0}, y_{0}, z_{0})` of the equation
    `ax^2 + by^2 = cz^2` with `a, b, c > 0` and `z_{0}^2 \geq \mid ab \mid` to
    a new reduced solution `(x, y, z)` such that `z^2 \leq \mid ab \mid`.
    """
    while z_0 > sqrt(a*b):

        if c % 2 == 0:
            k = c // 2
            u_0, v_0 = base_solution_linear(k, y_0, -x_0)

        else:
            k = 2*c
            u_0, v_0 = base_solution_linear(c, y_0, -x_0)

        w = -(a*u_0*x_0 + b*v_0*y_0) // (c*z_0)

        if c % 2 == 1:
            if w % 2 != (a*u_0 + b*v_0) % 2:
                w = w + 1

        x = (x_0*(a*u_0**2 + b*v_0**2 + c*w**2) - 2*u_0*(a*u_0*x_0 + b*v_0*y_0 + c*w*z_0)) // k
        y = (y_0*(a*u_0**2 + b*v_0**2 + c*w**2) - 2*v_0*(a*u_0*x_0 + b*v_0*y_0 + c*w*z_0)) // k
        z = (z_0*(a*u_0**2 + b*v_0**2 + c*w**2) - 2*w*(a*u_0*x_0 + b*v_0*y_0 + c*w*z_0)) // k

        x_0, y_0, z_0 = x, y, z

    return x_0, y_0, z_0


def diop_general_pythagorean(eq, param=symbols("m", integer=True)):
    """
    Solves the general pythagorean equation,
    `a_{1}^2x_{1}^2 + a_{2}^2x_{2}^2 + . . . + a_{n}^2x_{n}^2 - a_{n + 1}^2x_{n + 1}^2 = 0`.

    Returns a tuple which contains a parametrized solution to the equation,
    sorted in the same order as the input variables.

    Usage
    =====

    ``diop_general_pythagorean(eq, param)``: where ``eq`` is a general
    pythagorean equation which is assumed to be zero and ``param`` is the base
    parameter used to construct other parameters by subscripting.

    Examples
    ========

    >>> from sympy.solvers.diophantine import diop_general_pythagorean
    >>> from sympy.abc import a, b, c, d, e
    >>> diop_general_pythagorean(a**2 + b**2 + c**2 - d**2)
    (m1**2 + m2**2 - m3**2, 2*m1*m3, 2*m2*m3, m1**2 + m2**2 + m3**2)
    >>> diop_general_pythagorean(9*a**2 - 4*b**2 + 16*c**2 + 25*d**2 + e**2)
    (10*m1**2  + 10*m2**2  + 10*m3**2 - 10*m4**2, 15*m1**2  + 15*m2**2  + 15*m3**2  + 15*m4**2, 15*m1*m4, 12*m2*m4, 60*m3*m4)
    """
    var, coeff, diop_type  = classify_diop(eq)

    if diop_type == "general_pythagorean":
        return _diop_general_pythagorean(var, coeff, param)


def _diop_general_pythagorean(var, coeff, t):

    if sign(coeff[var[0]**2]) + sign(coeff[var[1]**2]) + sign(coeff[var[2]**2]) < 0:
        for key in coeff.keys():
            coeff[key] = coeff[key] * -1

    n = len(var)
    index = 0

    for i, v in enumerate(var):
        if sign(coeff[v**2]) == -1:
            index = i

    m = symbols(str(t) + "1:" + str(n), integer=True)
    l = []
    ith = 0

    for m_i in m:
        ith = ith + m_i**2

    l.append(ith - 2*m[n - 2]**2)

    for i in range(n - 2):
        l.append(2*m[i]*m[n-2])

    sol = l[:index] + [ith] + l[index:]

    lcm = 1
    for i, v in enumerate(var):
        if i == index or (index > 0 and i == 0) or (index == 0 and i == 1):
            lcm = ilcm(lcm, sqrt(abs(coeff[v**2])))
        else:
            lcm = ilcm(lcm, sqrt(coeff[v**2]) if sqrt(coeff[v**2]) % 2 else sqrt(coeff[v**2]) // 2)

    for i, v in enumerate(var):
        sol[i] = (lcm*sol[i]) / sqrt(abs(coeff[v**2]))

    return tuple(sol)


def diop_general_sum_of_squares(eq, limit=1):
    """
    Solves the equation `x_{1}^2 + x_{2}^2 + . . . + x_{n}^2 - k = 0`.

    Returns at most ``limit`` number of solutions. Currently there is no way to
    set ``limit`` using higher level API's like ``diophantine()`` or
    ``diop_solve()`` but that will be fixed soon.

    Usage
    =====

    ``general_sum_of_squares(eq, limit)`` : Here ``eq`` is an expression which
    is assumed to be zero. Also, ``eq`` should be in the form,
    `x_{1}^2 + x_{2}^2 + . . . + x_{n}^2 - k = 0`. At most ``limit`` number of
    solutions are returned.

    Details
    =======

    When `n = 3` if `k = 4^a(8m + 7)` for some `a, m \in Z` then there will be
    no solutions. Refer [1]_ for more details.

    Examples
    ========

    >>> from sympy.solvers.diophantine import diop_general_sum_of_squares
    >>> from sympy.abc import a, b, c, d, e, f
    >>> diop_general_sum_of_squares(a**2 + b**2 + c**2 + d**2 + e**2 - 2345)
    set([(0, 48, 5, 4, 0)])

    Reference
    =========

    .. [1] Representing an Integer as a sum of three squares, [online],
        Available:
        http://www.proofwiki.org/wiki/Integer_as_Sum_of_Three_Squares
    """
    var, coeff, diop_type = classify_diop(eq)

    if diop_type == "general_sum_of_squares":
        return _diop_general_sum_of_squares(var, coeff, limit)


def _diop_general_sum_of_squares(var, coeff, limit=1):

    n = len(var)
    k = -int(coeff[Integer(1)])
    s = set([])

    if k < 0:
        return set([])

    if n == 3:
        s.add(sum_of_three_squares(k))
    elif n == 4:
        s.add(sum_of_four_squares(k))
    else:

        m = n // 4
        f = partition(k, m, True)

        for j in range(limit):

            soln = []
            try:
                l = next(f)
            except StopIteration:
                break

            for n_i in l:
                a, b, c, d = sum_of_four_squares(n_i)
                soln = soln + [a, b, c, d]

            soln = soln + [0] * (n % 4)

            s.add(tuple(soln))

    return s


## Functions below this comment can be more suitably grouped under an Additive number theory module
## rather than the Diophantine equation module.


def partition(n, k=None, zeros=False):
    """
    Returns a generator that can be used to generate partitions of an integer
    `n`.

    A partition of `n` is a set of positive integers which add upto `n`. For
    example, partitions of 3 are 3 , 1 + 2, 1 + 1+ 1. A partition is returned
    as a tuple. If ``k`` equals None, then all possible partitions are returned
    irrespective of their size, otherwise only the partitions of size ``k`` are
    returned. If there are no partions of `n` with size `k` then an empty tuple
    is returned. If the ``zero`` parameter is set to True then a suitable
    number of zeros are added at the end of every partition of size less than
    ``k``.

    ``zero`` parameter is considered only if ``k`` is not None. When the
    partitions are over, the last `next()` call throws the ``StopIteration``
    exception, so this function should always be used inside a try - except
    block.

    Details
    =======

    ``partition(n, k)``: Here ``n`` is a positive integer and ``k`` is the size
    of the partition which is also positive integer.

    Examples
    ========

    >>> from sympy.solvers.diophantine import partition
    >>> f = partition(5)
    >>> next(f)
    (1, 1, 1, 1, 1)
    >>> next(f)
    (1, 1, 1, 2)
    >>> g = partition(5, 3)
    >>> next(g)
    (3, 1, 1)
    >>> next(g)
    (2, 2, 1)

    Reference
    =========

    .. [1] Generating Integer Partitions, [online],
        Available: http://jeromekelleher.net/partitions.php
    """
    if n < 1:
        yield tuple()

    if k is not None:
        if k < 1:
            yield tuple()

        elif k > n:
            if zeros:
                for i in range(1, n):
                    for t in partition(n, i):
                        yield (t,) + (0,) * (k - i)
            else:
                yield tuple()

        else:
            a = [1 for i in range(k)]
            a[0] = n - k + 1

            yield tuple(a)

            i = 1
            while a[0] >= n // k + 1:
                j = 0

                while j < i and j + 1 < k:
                    a[j] = a[j] - 1
                    a[j + 1] = a[j + 1] + 1

                    yield tuple(a)

                    j = j + 1

                i = i + 1

            if zeros:
                for m in range(1, k):
                    for a in partition(n, m):
                        yield tuple(a) + (0,) * (k - m)

    else:
        a = [0 for i in range(n + 1)]
        l = 1
        y = n - 1

        while l != 0:
            x = a[l - 1] + 1
            l -= 1

            while 2*x <= y:
                a[l] = x
                y -= x
                l += 1

            m = l + 1
            while x <= y:
                a[l] = x
                a[m] = y
                yield tuple(a[:l + 2])
                x += 1
                y -= 1

            a[l] = x + y
            y = x + y - 1
            yield tuple(a[:l + 1])


def prime_as_sum_of_two_squares(p):
    """
    Represent a prime `p` which is congruent to 1 mod 4, as a sum of two
    squares.

    Examples
    ========

    >>> from sympy.solvers.diophantine import prime_as_sum_of_two_squares
    >>> prime_as_sum_of_two_squares(5)
    (2, 1)

    Reference
    =========

    .. [1] Representing a number as a sum of four squares, [online],
        Available: http://www.schorn.ch/howto.html
    """
    if p % 8 == 5:
        b = 2
    else:
        b = 3

        while pow(b, (p - 1) // 2, p) == 1:
            b = nextprime(b)

    b = pow(b, (p - 1) // 4, p)
    a = p

    while b**2 > p:
        a, b = b, a % b

    return (b, a % b)


def sum_of_three_squares(n):
    """
    Returns a 3-tuple `(a, b, c)` such that `a^2 + b^2 + c^2 = n` and
    `a, b, c \geq 0`.

    Returns (None, None, None) if `n = 4^a(8m + 7)` for some `a, m \in Z`. See
    [1]_ for more details.

    Usage
    =====

    ``sum_of_three_squares(n)``: Here ``n`` is a non-negative integer.

    Examples
    ========

    >>> from sympy.solvers.diophantine import sum_of_three_squares
    >>> sum_of_three_squares(44542)
    (207, 37, 18)

    References
    ==========

    .. [1] Representing a number as a sum of three squares, [online],
        Available: http://www.schorn.ch/howto.html
    """
    special = {1:(1, 0, 0), 2:(1, 1, 0), 3:(1, 1, 1), 10: (1, 3, 0), 34: (3, 3, 4), 58:(3, 7, 0),
        85:(6, 7, 0), 130:(3, 11, 0), 214:(3, 6, 13), 226:(8, 9, 9), 370:(8, 9, 15),
        526:(6, 7, 21), 706:(15, 15, 16), 730:(1, 27, 0), 1414:(6, 17, 33), 1906:(13, 21, 36),
        2986: (21, 32, 39), 9634: (56, 57, 57)}

    v = 0

    if n == 0:
        return (0, 0, 0)

    while n % 4 == 0:
        v = v + 1
        n = n // 4

    if n % 8 == 7:
        return (None, None, None)

    if n in special.keys():
        x, y, z = special[n]
        return (2**v*x, 2**v*y, 2**v*z)

    l = int(sqrt(n))

    if n == l**2:
        return (2**v*l, 0, 0)

    x = None

    if n % 8 == 3:
        l = l if l % 2 else l - 1

        for i in range(l, -1, -2):
            if isprime((n - i**2) // 2):
                x = i
                break

        y, z = prime_as_sum_of_two_squares((n - x**2) // 2)
        return (2**v*x, 2**v*(y + z), 2**v*abs(y - z))

    if n % 8 == 2 or n % 8 == 6:
        l = l if l % 2 else l - 1
    else:
        l = l - 1 if l % 2 else l

    for i in range(l, -1, -2):
        if isprime(n - i**2):
            x = i
            break

    y, z = prime_as_sum_of_two_squares(n - x**2)
    return (2**v*x, 2**v*y, 2**v*z)


def sum_of_four_squares(n):
    """
    Returns a 4-tuple `(a, b, c, d)` such that `a^2 + b^2 + c^2 + d^2 = n`.

    Here `a, b, c, d \geq 0`.

    Usage
    =====

    ``sum_of_four_squares(n)``: Here ``n`` is a non-negative integer.

    Examples
    ========

    >>> from sympy.solvers.diophantine import sum_of_four_squares
    >>> sum_of_four_squares(3456)
    (8, 48, 32, 8)
    >>> sum_of_four_squares(1294585930293)
    (0, 1137796, 2161, 1234)

    References
    ==========

    .. [1] Representing a number as a sum of four squares, [online],
        Available: http://www.schorn.ch/howto.html
    """
    if n == 0:
        return (0, 0, 0, 0)

    v = 0
    while n % 4 == 0:
        v = v + 1
        n = n // 4

    if n % 8 == 7:
        d = 2
        n = n - 4
    elif n % 8 == 6 or n % 8 == 2:
        d = 1
        n = n - 1
    else:
        d = 0

    x, y, z = sum_of_three_squares(n)

    return (2**v*d, 2**v*x, 2**v*y, 2**v*z)


def power_representation(n, p, k, zeros=False):
    """
    Returns a generator for finding k-tuples `(n_{1}, n_{2}, . . . n_{k})` such
    that `n = n_{1}^p + n_{2}^p + . . . n_{k}^p`.

    Here `n` is a non-negative integer. StopIteration exception is raised after
    all the solutions are generated, so should always be used within a try-
    catch block.

    Usage
    =====

    ``power_representation(n, p, k, zeros)``: Represent number ``n`` as a sum
    of ``k``, ``p``th powers. If ``zeros`` is true, then the solutions will
    contain zeros.

    Examples
    ========

    >>> from sympy.solvers.diophantine import power_representation
    >>> f = power_representation(1729, 3, 2) # Represent 1729 as a sum of two cubes
    >>> next(f)
    (12, 1)
    >>> next(f)
    (10, 9)
    """
    if p < 1 or k < 1 or n < 1:
        raise ValueError("Expected: n > 0 and k >= 1 and p >= 1")

    if k == 1:
        if perfect_power(n):
            yield (perfect_power(n)[0],)
        else:
            yield tuple()

    elif p == 1:
        for t in partition(n, k, zeros):
            yield t

    else:
        l = []
        a = integer_nthroot(n, p)[0]

        for t in pow_rep_recursive(a, k, n, [], p):
                yield t

        if zeros:
            for i in range(2, k):
                for t in pow_rep_recursive(a, i, n, [], p):
                    yield t + (0,) * (k - i)


def pow_rep_recursive(n_i, k, n_remaining, terms, p):

    if k == 0 and n_remaining == 0:
        yield tuple(terms)
    else:
        if n_i >= 1 and k > 0 and n_remaining >= 0:
            if n_i**p <= n_remaining:
                for t in pow_rep_recursive(n_i, k - 1, n_remaining - n_i**p, terms + [n_i], p):
                    yield t

            for t in pow_rep_recursive(n_i - 1, k, n_remaining, terms, p):
                yield t
