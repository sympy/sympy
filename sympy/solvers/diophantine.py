from sympy import (degree_list, Poly, igcd, divisors, sign, symbols, S, Integer, Wild, Symbol)
from sympy import (Add, Mul, solve, ceiling, floor, sqrt, sympify, simplify, Subs, ilcm, Matrix)

from sympy.simplify.simplify import rad_rationalize
from sympy.ntheory.modular import solve_congruence


def diop_solve(eq, param=symbols("t", Integer=True)):
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
    {x: 3*t - 5, y: -2*t + 5}
    >>> diop_solve(4*x + 3*y -4*z + 5)
    {x: 3*t + 4*z - 5, y: -4*t - 4*z + 5, z: z}
    >>> diop_solve(x + 3*y - 4*z + w -6)
    {w: t, x: -t - 3*y + 4*z + 6, y: y, z: z}
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
    eq = eq.expand(force=True)
    var = list(eq.free_symbols)
    var.sort()

    coeff = {}
    diop_type = None

    coeff = dict([reversed(t.as_independent(*var)) for t in eq.args])
    for v in coeff:
        if not isinstance(coeff[v], Integer):
            raise TypeError("Coefficients should be Integers")

    if len(var) == 1:
        diop_type = "univariable"
    elif Poly(eq).total_degree() == 1:
        diop_type = "linear"
    elif Poly(eq).total_degree() == 2:
        diop_type = "quadratic"
        x = var[0]
        y = var[1]

        if isinstance(eq, Mul):
            coeff = {x**2: 0, x*y: eq.args[0], y**2: 0, x: 0, y: 0, Integer(1): 0}
        else:
            for term in [x**2, y**2, x*y, x, y, Integer(1)]:
                if term not in coeff.keys():
                    coeff[term] = Integer(0)
    else:
        raise NotImplementedError("Still not implemented")

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
    {x: -3*t - 5, y: -2*t - 5}
    >>> diop_linear([x, y, z], {Integer(1): -3, x: 2, y: -3, z: -4}, t) # 2*x - 3*y - 4*z - 3 = 0
    {x: -3*t - 4*z - 3, y: -2*t - 4*z - 3, z: z}
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

    # (5) B**2 - 4*A*C > 0

    elif B**2 - 4*A*C > 0:
        # Method used when B**2 - 4*A*C is a square, is descibed in p. 6 of the below paper
        # by John P. Robertson.
        # http://www.jpr2718.org/ax2p.pdf

        if isinstance(sqrt(B**2 - 4*A*C), Integer):
            if A != 0:
                r = sqrt(B**2 - 4*A*C)
                u, v = symbols("u, v", integer=True)
                eq = simplify(4*A*r*u*v + 4*A*D*(B*v + r*u + r*v - B*u) + 2*A*4*A*E*(u - v) + 4*A*r*4*A*F)

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
                var[0], var[1] = var[1], var[0] # Interchange x and y
                s = diop_quadratic(var, coeff, t)

                while len(s) > 0:
                    sol = s.pop()
                    l.add((sol[1], sol[0]))

        else:
            # In this case equation can be transformed into a Pell equation
            A, B = _transformation_to_pell(var, coeff)
            D, N = _find_DN(var, coeff)
            solns_pell = diop_pell(D, N)

            n = symbols("n", integer=True)

            a = diop_pell(D, 1)
            T = a[0][0]
            U = a[0][1]

            if (isinstance(A[0], Integer) and isinstance(A[1], Integer) and isinstance(A[2], Integer)
                and isinstance(A[3], Integer) and isinstance(B[0], Integer) and isinstance(B[1], Integer)):
                for sol in solns_pell:

                    r = sol[0]
                    s = sol[1]
                    x_n = S((r + s*sqrt(D))*(T + U*sqrt(D))**n + (r - s*sqrt(D))*(T - U*sqrt(D))**n)/2
                    y_n = S((r + s*sqrt(D))*(T + U*sqrt(D))**n - (r - s*sqrt(D))*(T - U*sqrt(D))**n)/(2*sqrt(D))

                    x_n, y_n = (A*Matrix([x_n, y_n]) + B)[0], (A*Matrix([x_n, y_n]) + B)[1]

                    l.add((x_n, y_n))

            else:
                L = ilcm(S(A[0]).q, ilcm(S(A[1]).q, ilcm(S(A[2]).q, ilcm(S(A[3]).q, ilcm(S(B[0]).q, S(B[1]).q)))))

                k = 0
                done = False
                T_k = T
                U_k = U

                while not done:
                    k = k + 1
                    if (T_k - 1) % L == 0 and U_k % L == 0:
                        done = True
                    T_k, U_k = T_k*T + D*U_k*U, T_k*U + U_k*T

                for soln in solns_pell:
                    x_0 = soln[0]
                    y_0 = soln[1]

                    x_i = x_0
                    y_i = y_0

                    for i in range(k):

                        X = (A*Matrix([x_i, y_i]) + B)[0]
                        Y = (A*Matrix([x_i, y_i]) + B)[1]

                        if isinstance(X, Integer) and isinstance(Y, Integer):
                            if is_solution_quad(var, coeff, X, Y):
                                x_n = S( (x_i + sqrt(D)*y_i)*(T + sqrt(D)*U)**(n*L) + (x_i - sqrt(D)*y_i)*(T - sqrt(D)*U)**(n*L) )/ 2
                                y_n = S( (x_i + sqrt(D)*y_i)*(T + sqrt(D)*U)**(n*L) - (x_i - sqrt(D)*y_i)*(T - sqrt(D)*U)**(n*L) )/ (2*sqrt(D))

                                x_n, y_n = (A*Matrix([x_n, y_n]) + B)[0], (A*Matrix([x_n, y_n]) + B)[1]
                                l.add((x_n, y_n))

                        x_i = x_i*T + D*U*y_i
                        y_i = x_i*U + y_i*T

    return l


def is_solution_quad(var, coeff, u, v):
    """
    Check whether (u, v) is solution to the quadratic diophantine equation.
    """
    x = var[0]
    y = var[1]

    eq = x**2*coeff[x**2] + x*y*coeff[x*y] + y**2*coeff[y**2] + x*coeff[x] + y*coeff[y] + coeff[Integer(1)]

    return simplify(Subs(eq, (x, y), (u, v)).doit()) == 0


def diop_pell(D, N, t=symbols("t", Integer=True)):
    """
    Solves the generalized Pell equation x**2 - D*y**2 = N. Uses LMM algorithm.
    Refer [1] for more details on the algorithm. Returns one solution for each
    class of the solutions. Other solutions can be constructed according to the
    values of D and N. Returns a list containing the solution tuples (x, y).

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
        elif N > 0: # TODO: Solution method should be improved
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
                        # TODO: efficient algorithm should be used to
                        # solve z^2 = D (mod m)
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


def diop_bf_pell(D, N, t=symbols("t", Integer=True)):
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
    Returns True if two solutions to the x**2 - D*y**2 = N belongs to the same
    equivalence class and False otherwise. Two solutions (u, v) and (r, s) to
    the above equation falls to the same equivalence class iff both (u*r - D*v*s)
    and (u*s - v*r) are divisible by N. See reference [1]. No check is performed to
    test whether (u, v) and (r, s) are actually solutions to the equation.
    User should take care of this.

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


def transformation_to_pell(eq):
    """
    This function transforms general quadratic ax**2 + bxy + cy**2 + dx + ey + f = 0
    to generalized pell equation X**2 - DY**2 = N when delta = b**2 - 4*a*c > 0
    and delta is square free. It can be easily noted that both a and b are
    non zero in this case.

    This can be used to solve the general quadratic with above restrictions
    by transforming it to the Pell equation. Refer [1] for more detailed information
    on the transformation. This function returns a tuple (A, B) where A is 2 * 2
    matrix and B is a 2 * 1 matrix such that,

    Transpose((x y)) =  A * Transpose((X Y)) + B

    Usage
    =====

        transformation_to_pell(eq) -> where eq is the quadratic to be transformed.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.solvers.diophantine import transformation_to_pell
    >>> from sympy.solvers.diophantine import classify_diop
    >>> A, B = transformation_to_pell(x**2 - 3*x*y - y**2 - 2*y + 1)
    >>> A
    Matrix([
    [1/26, 3/26],
    [   0, 1/13]])
    >>> B
    Matrix([
    [-6/13],
    [-4/13]])

    A, B  returned are such that Transpose((x y)) =  A * Transpose((X Y)) + B.
    Substituting these values for x and y and a bit of simplifying work will give
    a pell type equation.

    >>> from sympy.abc import X, Y
    >>> from sympy import Matrix, simplify, Subs
    >>> u = (A*Matrix([X, Y]) + B)[0] # Transformation for x
    >>> u
    X/26 + 3*Y/26 - 6/13
    >>> v = (A*Matrix([X, Y]) + B)[1] # Transformation for y
    >>> v
    Y/13 - 4/13

    Next we will substitute these formulas for x and y and simplify.

    >>> eq = simplify(Subs(x**2 - 3*x*y - y**2 - 2*y + 1, (x, y), (u, v)).doit())
    >>> eq
    X**2/676 - Y**2/52 + 17/13

    By multiplying the denominator appropriately, we can get the Pell equation
    in standard form.

    >>> eq * 676
    X**2 - 13*Y**2 + 884

    If only the final Pell equation is needed, find_DN() can be used.

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
    if diop_type == "quadratic":
        return _transformation_to_pell(var, coeff)


def _transformation_to_pell(var, coeff):

    x = var[0]
    y = var[1]

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
        A_0, B_0 = _transformation_to_pell([X, Y], coeff)
        return Matrix(2, 2, [S(1)/B, -S(C)/B, 0, 1])*A_0, Matrix(2, 2, [S(1)/B, -S(C)/B, 0, 1])*B_0

    else:
        if d != Integer(0):
            B = (S(2*a)/d).p
            C = (S(2*a)/d).q
            A = (S(a)/B**2).p
            T = (S(a)/B**2).q

            # eq_2 = A*X**2 + c*T*Y**2 + e*T*Y + f*T - A*C**2
            coeff = {X**2: A, X*Y: 0, Y**2: c*T, X: 0, Y: e*T, Integer(1): f*T - A*C**2}
            A_0, B_0 = _transformation_to_pell([X, Y], coeff)
            return Matrix(2, 2, [S(1)/B, 0, 0, 1])*A_0, Matrix(2, 2, [S(1)/B, 0, 0, 1])*B_0 + Matrix([-S(C)/B, 0])

        else:
            if e != Integer(0):
                B = (S(2*c)/e).p
                C = (S(2*c)/e).q
                A = (S(c)/B**2).p
                T = (S(c)/B**2).q

                # eq_3 = a*T*X**2 + A*Y**2 + f*T - A*C**2
                coeff = {X**2: a*T, X*Y: 0, Y**2: A, X: 0, Y: 0, Integer(1): f*T - A*C**2}
                A_0, B_0 = _transformation_to_pell([X, Y], coeff)
                return Matrix(2, 2, [1, 0, 0, S(1)/B])*A_0, Matrix(2, 2, [1, 0, 0, S(1)/B])*B_0 + Matrix([0, -S(C)/B])

            else:
                # TODO: pre-simplification: Not necessary but may simplify
                # the equation.

                return Matrix(2, 2, [S(1)/a, 0, 0, 1]), Matrix([0, 0])


def find_DN(eq):
    """
    This function returns a tuple, (D, N) of the Pell equation corresponding to
    the general quadratic ax**2 + bxy + cy**2 + dx + ey + f = 0 when
    delta = b**2 - 4*a*c > 0 and delta is not a perfect square.

    Solving the general quadratic with above restrictons is equivalent to solving
    the generalized Pell equation X**2 - D*Y**2 = N and transforming the solutions
    by using the transformation matrices returned by transformation_to_pell().

    Usage
    =====

        find_DN(eq) -> where eq is the quadratic to be transformed.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.solvers.diophantine import find_DN
    >>> find_DN(x**2 - 3*x*y - y**2 - 2*y + 1)
    (13, -884)

    The result means that after transforming x**2 - 3*x*y - y**2 - 2*y + 1
    using transformation_to_pell(), this can be simplified into X**2 -13*Y**2 = -884.

    See Also
    ========

    transformation_to_pell()

    References
    ==========

    .. [1] Solving the equation ax^2 + bxy + cy^2 + dx + ey + f = 0,
           John P.Robertson, May 8, 2003, Page 7 - 11.
           http://www.jpr2718.org/ax2p.pdf
    """
    var, coeff, diop_type = classify_diop(eq)
    if diop_type == "quadratic":
        return _find_DN(var, coeff)


def _find_DN(var, coeff):

    x = var[0]
    y = var[1]
    X, Y = symbols("X, Y", integer=True)
    A , B = _transformation_to_pell(var, coeff)

    u = (A*Matrix([X, Y]) + B)[0]
    v = (A*Matrix([X, Y]) + B)[1]
    eq = x**2*coeff[x**2] + x*y*coeff[x*y] + y**2*coeff[y**2] + x*coeff[x] + y*coeff[y] + coeff[Integer(1)]

    simplified = simplify(Subs(eq, (x, y), (u, v)).doit())

    coeff = dict([reversed(t.as_independent(*[X, Y])) for t in simplified.args])

    for term in [X**2, Y**2, Integer(1)]:
        if term not in coeff.keys():
            coeff[term] = Integer(0)

    return -coeff[Y**2]/coeff[X**2], -coeff[Integer(1)]/coeff[X**2]


def check_param(x, y, a, t):
    """
    Check if there is a number modulo a such that x and y are both
    integers. If exist, then find a parametric representation for x and y.
    """
    k, m, n = symbols("k, m, n", Integer=True)
    p = Wild("p", exclude=[k])
    q = Wild("q", exclude=[k])
    ok = False

    for i in range(a):

        z_x = simplify(Subs(x, t, a*k + i).doit()).match(p*k + q)
        z_y = simplify(Subs(y, t, a*k + i).doit()).match(p*k + q)

        if (isinstance(z_x[p], Integer) and isinstance(z_x[q], Integer) and
            isinstance(z_y[p], Integer) and isinstance(z_y[q], Integer)):
            ok = True
            break

    if ok == True:

        x_param = x.match(p*t + q)
        y_param = y.match(p*t + q)
        eq = S(m -x_param[q])/x_param[p] - S(n - y_param[q])/y_param[p]

        lcm_denom, junk = Poly(eq).clear_denoms()
        eq = eq * lcm_denom

        return diop_solve(eq, t)[m], diop_solve(eq, t)[n]
    else:
        return (None, None)
