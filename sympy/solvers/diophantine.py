from sympy import degree_list
from sympy import igcd
from sympy import symbols
from sympy import Add
from sympy import Integer
from sympy import sign
from sympy import S
from sympy import Poly
from sympy import divisors
from sympy import solve
from sympy import ceiling, floor
from sympy import sqrt


def diop_solve(eq, param=symbols("t", integer=True)):
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
    var, coeff, eq_type = classify_diop(eq)

    if eq_type == "linear":
        return diop_linear(var, coeff, param)
    elif eq_type == "quadratic":
        return diop_quadratic(var, coeff, param)

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

    if Poly(eq).total_degree() == 1:
        diop_type = "linear"
    elif Poly(eq).total_degree() == 2 and len(var) == 2:
        diop_type = "quadratic"
    else:
        raise NotImplementedError("Still not implemented")


    coeff = dict([reversed(t.as_independent(*var)) for t in eq.args])
    for v in coeff:
        if not isinstance(coeff[v], Integer):
            raise TypeError("Coefficients should be Integers")

    return var, coeff, diop_type


def diop_linear(var, coeff, param):
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
        ``param`` parameter to be used in the solution.

    Examples
    ========

    >>> from sympy.solvers.diophantine import diop_linear
    >>> from sympy.abc import x, y, z, t
    >>> from sympy import Integer
    >>> diop_linear([x, y],{Integer(1): -5, x: 2, y:-3}, t) #solves equation 2*x - 3*y -5 = 0
    {x: -15*t - 5, y: -10*t - 5}
    >>> diop_linear([x, y, z], {Integer(1): -3, x: 2, y: -3, z: -4}, t) # 2*x - 3*y - 4*z - 3= 0
    {x: -9*t - 4*z - 3, y: -6*t - 4*z - 3, z: z}
    """
    x = var[0]; y = var[1]
    a = coeff[x]; b = coeff[y]

    if len(var) == len(coeff):
        c = 0
    else:
        c = -coeff[Integer(1)]

    if len(var) == 2:
        sol_x, sol_y = base_solution_linear(c, a, b, param)
        return {x: sol_x, y: sol_y}

    elif len(var) > 2:
        X = []; Y = []

        for v in var[2:]:
            sol_x, sol_y  = base_solution_linear(-coeff[v], a, b)
            X.append(sol_x*v); Y.append(sol_y*v)

        sol_x, sol_y = base_solution_linear(c, a, b, param)
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


def base_solution_linear(c, a, b, t=None):
    """
    Return the base solution for a linear diophantine equation with two
    variables. Called repeatedly by diop_linear().

    **Usage**

        base_solution_linear(c, a, b, param) -> a, b, c are Integers as in a*x + b*y = c
        and param is the parameter to be used in the solution.

    **Details**

        ``c`` is the constant term in a*x + b*y = c
        ``a`` is the integer coefficient of x in a*x + b*y = c
        ``b`` is the integer coefficient of y in a*x + b*y = c
        ``param`` is the parameter to be used in the solution


    Examples
    ========

    >>> from sympy.solvers.diophantine import base_solution_linear
    >>> from sympy.abc import t
    >>> base_solution_linear(5, 2, 3) # equation 2*x + 3*y = 5
    (-5, 5)
    >>> base_solution_linear(0, 5, 7) # equation 5*x + 7*y = 0
    (0, 0)
    >>> base_solution_linear(5, 2, 3, t) # equation 2*x + 3*y = 5
    (15*t - 5, -10*t + 5)
    >>> base_solution_linear(0, 5, 7, t) # equation 5*x + 7*y = 0
    (7*t, -5*t)
    """
    d = igcd(a, igcd(b, c))
    a = a // d; b = b // d; c = c // d

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


def divisible(a, b):
    return igcd(int(a), int(b)) == abs(int(b))


def diop_quadratic(var, coeff, t):
    """
    Solves quadratic diophantine equations, i.e equations of the form
    Ax**2 + Bxy + Cy**2 + Dx + Ey + F = 0. Returns an set containing
    the tuples (x, y) which contains the solutions.

    **Usage**

        diop_quadratic(var, coeff) -> var is a list of variables and
        coeff is a dictionary containing coefficients of the symbols.

    **Details**

        ``var`` a list which contains two variables x and y.
        ``coeff`` a dict which generally contains six key value pairs.
        The set of keys is {x**2, y**2, x*y, x, y, Integer(1)}.
        ``t`` the parameter to be used in the solution

    Examples
    ========

    >>> from sympy.abc import x, y, t
    >>> from sympy import Integer
    >>> from sympy.solvers.diophantine import diop_quadratic
    >>> diop_quadratic([x, y], {x**2: 1, y**2: 1, x*y: 0, x: 2, y: 2, Integer(1):2}, t)
    set([(-1, -1)])

    References
    ==========

    .. [1] http://www.alpertron.com.ar/METHODS.HTM
    """
    x = var[0]; y = var[1]

    for term in [x**2, y**2, x*y, x, y, Integer(1)]:
        if term not in coeff.keys():
            coeff[term] = Integer(0)

    A = coeff[x**2]; B = coeff[x*y]; C = coeff[y**2]
    D = coeff[x]; E = coeff[y]; F = coeff[Integer(1)]

    d = igcd(A, igcd(B, igcd(C, igcd(D, igcd(E, F)))))
    A = A // d; B = B // d; C = C // d;
    D = D // d; E = E // d; F = F // d;

    # (1) Linear case: A = B = C = 0 -> considered under linear diophantine equations

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

    # (3) Elliptical case: B**2 - 4AC < 0
    # More Details, http://www.alpertron.com.ar/METHODS.HTM#Ellipse
    # In this case x should lie between the roots of
    # (B**2 - 4AC)x**2 + 2(BE - 2CD)x + (E**2 - 4CF) = 0

    elif B**2 - 4*A*C < 0:

        z = symbols("z", real=True)
        roots = solve((B**2 - 4*A*C)*z**2 + 2*(B*E - 2*C*D)*z + E**2 - 4*C*F)

        solve_y = lambda x, e: (-(B*x + E) + e*sqrt((B*x + E)**2 - 4*C*(A*x**2 + D*x + F)))/(2*C)

        if len(roots) == 1 and isinstance(roots[0], Integer):
            x_vals = [roots[0]]
        elif len(roots) == 2:
            x_vals = [i for i in range(ceiling(min(roots)), ceiling(max(roots)))]
        else:
            x_vals = []

        for x0 in x_vals:
            if isinstance(sqrt((B*x0 + E)**2 - 4*C*(A*x0**2 + D*x0 + F)), Integer):
                if isinstance(solve_y(x0, 1), Integer):
                    l.add((Integer(x0), solve_y(x0, 1)))
                if isinstance(solve_y(x0, -1), Integer):
                    l.add((Integer(x0), solve_y(x0, -1)))

    # (4) Parabolic case: B**2 - 4*A*C = 0
    # There are two subcases to be considered in this case.
    # sqrt(c)D - sqrt(a)E = 0 and sqrt(c)D - sqrt(a)E != 0
    # More Details, http://www.alpertron.com.ar/METHODS.HTM#Parabol

    elif B**2 - 4*A*C == 0:

        g = igcd(A, C)
        g = abs(g) * sign(A)
        a = A // g; b = B // g; c = C // g
        e = sign(B/A)


        if e*sqrt(c)*D - sqrt(a)*E == 0:
            z = symbols("z", real=True)
            roots = solve(sqrt(a)*g*z**2 + D*z + sqrt(a)*F)

            for root in roots:
                if isinstance(root, Integer):
                    l.add((diop_solve(sqrt(a)*x + e*sqrt(c)*y - root)[x], diop_solve(sqrt(a)*x + e*sqrt(c)*y - root)[y]))

        elif isinstance(e*sqrt(c)*D - sqrt(a)*E, Integer):
            solve_x = lambda u: e*sqrt(c)*g*(sqrt(a)*E - e*sqrt(c)*D)*t**2 - (E + 2*e*sqrt(c)*g*u)*t\
                - (e*sqrt(c)*g*u**2 + E*u + e*sqrt(c)*F) // (e*sqrt(c)*D - sqrt(a)*E)

            solve_y = lambda u: sqrt(a)*g*(e*sqrt(c)*D - sqrt(a)*E)*t**2 + (D + 2*sqrt(a)*g*u)*t \
                + (sqrt(a)*g*u**2 + D*u + sqrt(a)*F) // (e*sqrt(c)*D - sqrt(a)*E)

            for z0 in range(0, abs(e*sqrt(c)*D - sqrt(a)*E)):
                if divisible(sqrt(a)*g*z0**2 + D*z0 + sqrt(a)*F, e*sqrt(c)*D - sqrt(a)*E):
                    l.add((solve_x(z0), solve_y(z0)))

    return l
