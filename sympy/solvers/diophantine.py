from sympy import (degree_list, Poly, igcd, divisors, sign, symbols, S, Integer, Add, Mul, solve, ceiling, floor, sqrt, sympify, simplify)
from sympy.simplify.simplify import rad_rationalize
from sympy.matrices import Matrix


def diop_solve(eq, param=symbols("t", integer=True)):
    """
    Solves diophantine equations. Uses classify_diop() to determine the
    type of eqaution and calls the appropriate solver function.

    Usage
    =====

        diop_solve(eq, t) -> Solve diophantine equation, eq.

    Details
    =======

        ``eq`` should be an expression which is assumed to be zero.
        ``t`` is a parameter to be used in the solution.

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
    elif eq_type == "univariable":
        return solve(eq)

def classify_diop(eq):
    """
    Helper routine used by diop_solve(). Returns a tuple containing the type of the
    diophantine equation along with the variables(free symbols) and their coefficients.
    Variables are returned as a list and coefficients are returned as a dict with
    the key being the variable name and the constant term is keyed to Integer(1).
    Type is an element in the set {"linear", "quadratic", "pell", "pythogorean", "exponential"}

    Usage
    =====

        classify_diop(eq) -> Return variables, coefficients and type in order.

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
    """
    var = list(eq.free_symbols)
    var.sort()

    coeff = {}
    diop_type = None

    if len(var) == 1:
        diop_type = "univariable"
        return var, coeff, diop_type

    if Poly(eq).total_degree() == 1:
        diop_type = "linear"
    elif Poly(eq).total_degree() == 2:
        diop_type = "quadratic"
        if isinstance(eq, Mul):
            x = var[0]
            y = var[1]
            coeff = {x**2: 0, x*y: eq.args[0], y**2: 0, x: 0, y: 0, Integer(1): 0}
            return var, coeff, diop_type
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

    Usage
    =====

        diop_linear(var, coeff) -> var is a list of variables and coeff is a dictionary
        containing coefficients of the symbols.

    Details
    =======

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
    >>> diop_linear([x, y], {Integer(1): -5, x: 2, y:-3}, t) #solves equation 2*x - 3*y -5 = 0
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

    Usage
    =====

        base_solution_linear(c, a, b, t) -> a, b, c are Integers as in a*x + b*y = c
        and t is the parameter to be used in the solution.

    Details
    =======

        ``c`` is the constant term in a*x + b*y = c
        ``a`` is the integer coefficient of x in a*x + b*y = c
        ``b`` is the integer coefficient of y in a*x + b*y = c
        ``t`` is the parameter to be used in the solution


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
                return (c*(x0 + b*t), c*(y0 - a*t))
            else:
                return (Integer(c*x0), Integer(c*y0))
        else:
            return (None, None)


def extended_euclid(a, b):
    """
    For given a, b returns a tuple containing integers x, y and d such that
    a*x + b*y = d. Here d = gcd(a, b).

    Usage
    =====

        extended_euclid(a, b) -> returns x, y and gcd(a, b).

    Details
    =======

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
    Ax**2 + Bxy + Cy**2 + Dx + Ey + F = 0. Returns a set containing
    the tuples (x, y) which contains the solutions.

    Usage
    =====

        diop_quadratic(var, coeff) -> var is a list of variables and
        coeff is a dictionary containing coefficients of the symbols.

    Details
    =======

        ``var`` a list which contains two variables x and y.
        ``coeff`` a dict which generally contains six key value pairs.
        The set of keys is {x**2, y**2, x*y, x, y, Integer(1)}.
        ``t`` the parameter to be used in the solution.

    Examples
    ========

    >>> from sympy.abc import x, y, t
    >>> from sympy import Integer
    >>> from sympy.solvers.diophantine import diop_quadratic
    >>> diop_quadratic([x, y], {x**2: 1, y**2: 1, x*y: 0, x: 2, y: 2, Integer(1): 2}, t)
    set([(-1, -1)])

    References
    ==========

    .. [1] http://www.alpertron.com.ar/METHODS.HTM
    """
    x = var[0]
    y = var[1]

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
            x_vals = [i for i in range(ceiling(min(roots)), ceiling(max(roots)))] # ceiling = floor +/- 1
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
        a = A // g
        b = B // g
        c = C // g
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

    elif B**2 - 4*A*C > 0:
        # Method used when B**2 - 4*A*C is a square, is descibed in p. 6 of the below paper
        # by John P. Robertson.
        # http://www.jpr2718.org/ax2p.pdf

        if isinstance(sqrt(B**2 - 4*A*C), Integer):
            if A != 0:
                r = sqrt(B**2 - 4*A*C)
                u, v = symbols("u, v", type = Integer)
                eq = simplify(4*A*r*u*v + 4*A*D*(B*v + r*u + r*v - B*u) + 2*A*4*A*E*(u - v) + 4*A*r*4*A*F)
                sol = diop_solve(eq)
                sol = list(sol)

                for solution in sol:
                    s0 = solution[0]
                    t0 = solution[1]
                    if divisible(B*t0 + r*s0 + r*t0 - B*s0, 4*A*r):
                        if divisible(s0 - t0, 2*r):
                                l.add(((B*t0 + r*s0 + r*t0 - B*s0)//(4*A*r), (s0 - t0)//(2*r)))
        else:
            # In this case, equation reduces to the generalized Pell equation, x**2 -D*y**2 = N.
            # An algorithm for the generalized Pell equation has been implemented.
            # Only have to transform this into a Pell equation and solve it.
            # Transformation is described in p. 7 of http://www.jpr2718.org/ax2p.pdf
            # Then recover solutions to the original equation from the solutions to the Pell
            # equation. (p. 13 of the above)
            raise NotImplementedError("Still not implemented")

    return l


def diop_pell(D, N, t=symbols("t", integer=True)):
    """
    Solves the generalized Pell equation x**2 - D*y**2 = N. Uses LMM algorithm.
    Refer [1] for more details on the algorithm. Returns only the fundamental solutions,
    other solutions can be constructed according to the values of D and N.
    Returns a list containing the solution tuples (x, y).

    Usage
    =====

        diop_pell(D, N, t) -> D and N are integers as in x**2 - D*y**2 = N and t is
        the parameter to be used in the solutions.

    Details
    =======

        ``D`` corresponds to the D in the equation
        ``N`` corresponds to the N in the equation
        ``t`` parameter to be used in the solutions

    Examples
    ========

    >>> from sympy.solvers.diophantine import diop_pell
    >>> diop_pell(13, -4) # Solves equation x**2 - 13*y**2 = -4
    [(3, 1), (393, 109), (36, 10)]

    The output can be interpreted as follows: There are three fundamental
    solutions to the equation x**2 - 13*y**2 = -4  given by (3, 1), (393, 109)
    and (36, 10). Each tuple is in the form (x, y), i. e solution (3, 1) means
    that x = 3 and y = 1.

    >>> diop_pell(986, 1) # Solves equation x**2 - 986*y**2 = 1
    [(49299, 1570)]

    References
    ==========

    ..[1] Solving the generalized Pell equation x**2 - D*y**2 = N, John P. Robertson,
          July 31, 2004, Pages 16 - 17.
          http://www.jpr2718.org/pell.pdf
    """
    if D < 0:
        if N == 0:
            return [(S.Zero, S.Zero)]
        elif N < 0:
            return []
        elif N > 0: # Solution method should be improved
            sol = []
            for y in range(floor(sqrt(-S(N)/D)) + 1):
                if isinstance(sqrt(N + D*y**2), Integer):
                    sol.append((sqrt(N + D*y**2), y))
            return sol

    elif D == 0:
        if N < 0 or not isinstance(sqrt(N), Integer):
            return []
        if N == 0:
            return [(S.Zero, t)]
        if isinstance(sqrt(N), Integer):
            return [(sqrt(N), t), (-sqrt(N), t)]

    else: # D > 0
        if isinstance(sqrt(D), Integer):
            r = sqrt(D)

            if N == 0:
                return [(r*t, t), (-r*t, t)]
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
                    zs = []

                    for i in range(floor(S(abs(m))/2) + 1):

                        if (i**2 - D) % abs(m) == 0:
                            zs.append(i)
                            if i < S(abs(m))/2 and i != 0:
                                zs.append(-i)

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

                                elif diop_pell(D, -1) != []:
                                    a = diop_pell(D, -1)
                                    sol.append((f*(r*a[0][0] + a[0][1]*s*D), f*(r*a[0][1] + s*a[0][0])))

                                break

                            l = l + 1
                            if l == length(z, abs(m), D):
                                break

                return sol


def PQa(P_0, Q_0, D):
    """
    Returns the useful information needed to solve the Pell equation.
    There are five sequences of integers defined related to the continued
    fraction representation of (P + sqrt(D))/Q, namely {P_i}, {Q_i}, {a_i},
    {A_i}, {B_i}, {G_i}. Return these values as a 6-tuple in the same order
    mentioned above. Refer [1] for more detailed information.

    Usage
    =====

        PQa(P_0, Q_0, D) -> P_0, Q_0 and D are integers corresponding to the
        continued fraction (P_0 + sqrt(D))/Q_0. Also it's assumed that
        P_0**2 == D mod(|Q_0|) and D is square free.

    Details
    =======

        ``P_0`` corresponds to the P in the continued fraction, (P + sqrt(D))/ Q
        ``D_0`` corresponds to the D in the continued fraction, (P + sqrt(D))/ Q
        ``Q_0`` corresponds to the Q in the continued fraction, (P + sqrt(D))/ Q

    Examples
    ========

    >>> from sympy.solvers.diophantine import PQa
    >>> from sympy.core.compatibility import next
    >>> pqa = PQa(13, 4, 5) # (13 + sqrt(5))/4
    >>> next(pqa) # (P_0, Q_0, a_0, A_0, B_0, G_0)
    (13, 4, 3, 3, 1, -1)
    >>> next(pqa) # (P_1, Q_1, a_1, A_1, B_1, G_1)
    (-1, 1, 1, 4, 1, 3)

    References
    ==========

    .. [1] Solving the generalized Pell equation x**2 - D*y**2 = N, John P. Robertson,
           July 31, 2004, Pages 4 - 8.
           http://www.jpr2718.org/pell.pdf
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


def diop_bf_pell(D, N, t=symbols("t", integer=True)):
    # Implemented mainly to find tests for diop_pell().
    # This method returns all most the same fundamental solutions as Wolfram Alpha.
    # Using these results and the results from diop_pell() in equivalent() will
    # determine whether they belong in the same equivalent class. That way we can
    # indirectly verify our results with Wolfram Alpha.
    """
    Uses brute force to solve the generalized Pell's equation, x**2 - D*y**2 = N.
    For more information refer [1]. Let t, u be the minimal positive solution such that
    t**2 - D*u**2 = 1 (i. e. solutions to the equation x**2 - D*y**2 = 1) then
    this method requires that sqrt(|N|*(t +/- 1) / (2*D)) is not too large.

    Usage
    =====

        diop_bf_pell(D, N, t) -> D and N are integers as in x**2 - D*y**2 = N and t is
        the parameter to be used in the solutions.

    Details
    =======

        ``D`` corresponds to the D in the equation
        ``N`` corresponds to the N in the equation
        ``t`` parameter to be used in the solutions

    Examples
    ========

    >>> from sympy.solvers.diophantine import diop_bf_pell
    >>> diop_bf_pell(13, -4)
    [(3, 1), (-3, 1), (36, 10)]
    >>> diop_bf_pell(986, 1)
    [(49299, 1570)]

    References
    ==========

    .. [1] Solving the generalized Pell equation x**2 - D*y**2 = N, John P. Robertson,
           July 31, 2004, Page 15.
           http://www.jpr2718.org/pell.pdf

    See Also
    ========

    diop_pell()
    """
    sol = []
    a = diop_pell(D, 1)
    u = a[0][0]
    v = a[0][1]


    if abs(N) == 1:
        return diop_pell(D, N)

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
    Returns True is two solutions to the x**2 - D*y**2 = N belongs to the same
    equivalence class and False otherwise. Two solutions (u, v) and (r, s) to
    the above equation falls to the same equivalence class iff both (u*r - D*v*s)
    and (u*s - v*r) are divisible by N. See reference [1]. No check is performed
    to test whether (u, v) and (r, s) are actually solutions to the equation. User
    should take care of this.

    Usage
    =====

        equivalent(u, v, r, s, D, N) -> (u, v) and (r, s) are two solutions of the
        equation x**2 - D*y**2 = N and all parameters involved are integers.

    Examples
    ========

    >>> from sympy.solvers.diophantine import equivalent
    >>> equivalent(18, 5, -18, -5, 13, -1)
    True
    >>> equivalent(3, 1, -18, 393, 109, -4)
    False

    References
    ==========

    .. [1] Solving the generalized Pell equation x**2 - D*y**2 = N, John P. Robertson,
           July 31, 2004, Page 12.
           http://www.jpr2718.org/pell.pdf

    """
    return divisible(u*r - D*v*s, N) and divisible(u*s - v*r, N)


def length(P, Q, D):
    """
    Returns the length of aperiodic part + length of periodic part of
    continued fraction representation of (P + sqrt(D))/Q. It is important
    to remember that this does NOT return the length of the periodic
    part but the addition of the legths of the two parts as mentioned above.

    Usage
    =====

        length(P, Q, D) -> P, Q and D are integers corresponding to the
        continued fraction (P + sqrt(D))/Q.

    Details
    =======

        ``P`` corresponds to the P in the continued fraction, (P + sqrt(D))/ Q
        ``D`` corresponds to the D in the continued fraction, (P + sqrt(D))/ Q
        ``Q`` corresponds to the Q in the continued fraction, (P + sqrt(D))/ Q

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
