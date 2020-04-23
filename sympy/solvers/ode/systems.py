from sympy import (Derivative, Symbol)
from sympy.core.numbers import I
from sympy.core.relational import Eq
from sympy.core.symbol import Dummy
from sympy.functions import exp, im, cos, sin, re
from sympy.functions.combinatorial.factorials import factorial
from sympy.matrices import zeros, Matrix
from sympy.simplify import simplify, collect
from sympy.solvers.deutils import ode_order
from sympy.solvers.solveset import NonlinearError
from sympy.utilities import numbered_symbols, default_sort_key
from sympy.utilities.iterables import ordered, uniq


def _get_func_order(eqs, funcs):
    return {func: max(ode_order(eq, func) for eq in eqs) for func in funcs}


class ODEOrderError(ValueError):
    """Raised by linear_ode_to_matrix if the system has the wrong order"""
    pass


class ODENonlinearError(NonlinearError):
    """Raised by linear_ode_to_matrix if the system is nonlinear"""
    pass


def linear_ode_to_matrix(eqs, funcs, t, order):
    r"""
    Convert a linear system of ODEs to matrix form

    Explanation
    ===========

    Express a system of linear ordinary differential equations as a single
    matrix differential equation [1]. For example the system $x' = x + y + 1$
    and $y' = x - y$ can be represented as

    .. math:: A_1 X' + A_0 X = b

    where $A_1$ and $A_0$ are $2 \times 2$ matrices and $b$, $X$ and $X'$ are
    $2 \times 1$ matrices with $X = [x, y]^T$.

    Higher-order systems are represented with additional matrices e.g. a
    second-order system would look like

    .. math:: A_2 X'' + A_1 X' + A_0 X = b

    Examples
    ========

    >>> from sympy import (Function, Symbol, Matrix, Eq)
    >>> from sympy.solvers.ode.systems import linear_ode_to_matrix
    >>> t = Symbol('t')
    >>> x = Function('x')
    >>> y = Function('y')

    We can create a system of linear ODEs like

    >>> eqs = [
    ...     Eq(x(t).diff(t), x(t) + y(t) + 1),
    ...     Eq(y(t).diff(t), x(t) - y(t)),
    ... ]
    >>> funcs = [x(t), y(t)]
    >>> order = 1 # 1st order system

    Now ``linear_ode_to_matrix`` can represent this as a matrix
    differential equation.

    >>> (A1, A0), b = linear_ode_to_matrix(eqs, funcs, t, order)
    >>> A1
    Matrix([
    [1, 0],
    [0, 1]])
    >>> A0
    Matrix([
    [-1, -1],
    [-1,  1]])
    >>> b
    Matrix([
    [1],
    [0]])

    The original equations can be recovered from these matrices:

    >>> eqs_mat = Matrix([eq.lhs - eq.rhs for eq in eqs])
    >>> X = Matrix(funcs)
    >>> A1 * X.diff(t) + A0 * X - b == eqs_mat
    True

    If the system of equations has a maximum order greater than the
    order of the system specified, a ODEOrderError exception is raised.

    >>> eqs = [Eq(x(t).diff(t, 2), x(t).diff(t) + x(t)), Eq(y(t).diff(t), y(t) + x(t))]
    >>> linear_ode_to_matrix(eqs, funcs, t, 1)
    Traceback (most recent call last):
    ...
    ODEOrderError: Cannot represent system in 1-order form

    If the system of equations is nonlinear, then ODENonlinearError is
    raised.

    >>> eqs = [Eq(x(t).diff(t), x(t) + y(t)), Eq(y(t).diff(t), y(t)**2 + x(t))]
    >>> linear_ode_to_matrix(eqs, funcs, t, 1)
    Traceback (most recent call last):
    ...
    ODENonlinearError: The system of ODEs is nonlinear.

    Parameters
    ==========

    eqs : list of sympy expressions or equalities
        The equations as expressions (assumed equal to zero).
    funcs : list of applied functions
        The dependent variables of the system of ODEs.
    t : symbol
        The independent variable.
    order : int
        The order of the system of ODEs.

    Returns
    =======

    The tuple ``(As, b)`` where ``As`` is a tuple of matrices and ``b`` is the
    the matrix representing the rhs of the matrix equation.

    Raises
    ======

    ODEOrderError
        When the system of ODEs have an order greater than what was specified
    ODENonlinearError
        When the system of ODEs is nonlinear

    See Also
    ========

    linear_eq_to_matrix: for systems of linear algebraic equations.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Matrix_differential_equation

    """
    from sympy.solvers.solveset import linear_eq_to_matrix

    if any(ode_order(eq, func) > order for eq in eqs for func in funcs):
        msg = "Cannot represent system in {}-order form"
        raise ODEOrderError(msg.format(order))

    As = []

    for o in range(order, -1, -1):
        # Work from the highest derivative down
        funcs_deriv = [func.diff(t, o) for func in funcs]

        # linear_eq_to_matrix expects a proper symbol so substitute e.g.
        # Derivative(x(t), t) for a Dummy.
        rep = {func_deriv: Dummy() for func_deriv in funcs_deriv}
        eqs = [eq.subs(rep) for eq in eqs]
        syms = [rep[func_deriv] for func_deriv in funcs_deriv]

        # Ai is the matrix for X(t).diff(t, o)
        # eqs is minus the remainder of the equations.
        try:
            Ai, b = linear_eq_to_matrix(eqs, syms)
        except NonlinearError:
            raise ODENonlinearError("The system of ODEs is nonlinear.")

        As.append(Ai)
        if o:
            eqs = [-eq for eq in b]
        else:
            rhs = b

    return As, rhs


def matrix_exp(A, t):
    r"""
    Matrix exponential $\exp(A*t)$ for the matrix ``A`` and scalar ``t``.

    Explanation
    ===========

    This functions returns the $\exp(A*t)$ by doing a simple
    matrix multiplication:

    .. math:: \exp(A*t) = P * expJ * P^{-1}

    where $expJ$ is $\exp(J*t)$. $J$ is the Jordan normal
    form of $A$ and $P$ is matrix such that:

    .. math:: A = P * J * P^{-1}

    The matrix exponential $\exp(A*t)$ appears in the solution of linear
    differential equations. For example if $x$ is a vector and $A$ is a matrix
    then the initial value problem

    .. math:: \frac{dx(t)}{dt} = A \times x(t),   x(0) = x0

    has the unique solution

    .. math:: x(t) = \exp(A t) x0

    Examples
    ========

    >>> from sympy import Symbol, Matrix, pprint
    >>> from sympy.solvers.ode.systems import matrix_exp
    >>> t = Symbol('t')

    We will consider a 2x2 matrix for comupting the exponential

    >>> A = Matrix([[2, -5], [2, -4]])
    >>> pprint(A)
    [2  -5]
    [     ]
    [2  -4]

    Now, exp(A*t) is given as follows:

    >>> pprint(matrix_exp(A, t))
    [   -t           -t                    -t              ]
    [3*e  *sin(t) + e  *cos(t)         -5*e  *sin(t)       ]
    [                                                      ]
    [         -t                     -t           -t       ]
    [      2*e  *sin(t)         - 3*e  *sin(t) + e  *cos(t)]

    Parameters
    ==========

    A : Matrix
        The matrix $A$ in the expression $\exp(A*t)$
    t : Symbol
        The independent variable

    See Also
    ========

    matrix_exp_jordan_form: For exponential of Jordan normal form

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Jordan_normal_form
    .. [2] https://en.wikipedia.org/wiki/Matrix_exponential

    """
    P, expJ = matrix_exp_jordan_form(A, t)
    return P * expJ * P.inv()


def matrix_exp_jordan_form(A, t):
    r"""
    Matrix exponential $\exp(A*t)$ for the matrix *A* and scalar *t*.

    Explanation
    ===========

    Returns the Jordan form of the $\exp(A*t)$ along with the matrix $P$ such that:

    .. math::
        \exp(A*t) = P * expJ * P^{-1}

    Examples
    ========

    >>> from sympy import Matrix, Symbol
    >>> from sympy.solvers.ode.systems import matrix_exp, matrix_exp_jordan_form
    >>> t = Symbol('t')

    We will consider a 2x2 defective matrix. This shows that our method
    works even for defective matrices.

    >>> A = Matrix([[1, 1], [0, 1]])

    It can be observed that this function gives us the Jordan normal form
    and the required invertible matrix P.

    >>> P, expJ = matrix_exp_jordan_form(A, t)

    Here, it is shown that P and expJ returned by this function is correct
    as they satisfy the formula: P * expJ * P_inverse = exp(A*t).

    >>> P * expJ * P.inv() == matrix_exp(A, t)
    True

    Parameters
    ==========

    A : Matrix
        The matrix $A$ in the expression $\exp(A*t)$
    t : Symbol
        The independent variable

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Defective_matrix
    .. [2] https://en.wikipedia.org/wiki/Jordan_matrix
    .. [3] https://en.wikipedia.org/wiki/Jordan_normal_form

    """

    N, M = A.shape
    if N != M:
        raise ValueError('Needed square matrix but got shape (%s, %s)' % (N, M))
    elif A.has(t):
        raise ValueError('Matrix A should not depend on t')

    def jordan_chains(A):
        '''Chains from Jordan normal form analogous to M.eigenvects().
        Returns a dict with eignevalues as keys like:
            {e1: [[v111,v112,...], [v121, v122,...]], e2:...}
        where vijk is the kth vector in the jth chain for eigenvalue i.
        '''
        P, blocks = A.jordan_cells()
        basis = [P[:,i] for i in range(P.shape[1])]
        n = 0
        chains = {}
        for b in blocks:
            eigval = b[0, 0]
            size = b.shape[0]
            if eigval not in chains:
                chains[eigval] = []
            chains[eigval].append(basis[n:n+size])
            n += size
        return chains

    eigenchains = jordan_chains(A)

    # Needed for consistency across Python versions:
    eigenchains_iter = sorted(eigenchains.items(), key=default_sort_key)
    isreal = not A.has(I)

    blocks = []
    vectors = []
    seen_conjugate = set()
    for e, chains in eigenchains_iter:
        for chain in chains:
            n = len(chain)
            if isreal and e != e.conjugate() and e.conjugate() in eigenchains:
                if e in seen_conjugate:
                    continue
                seen_conjugate.add(e.conjugate())
                exprt = exp(re(e) * t)
                imrt = im(e) * t
                imblock = Matrix([[cos(imrt), sin(imrt)],
                                  [-sin(imrt), cos(imrt)]])
                expJblock2 = Matrix(n, n, lambda i,j:
                        imblock * t**(j-i) / factorial(j-i) if j >= i
                        else zeros(2, 2))
                expJblock = Matrix(2*n, 2*n, lambda i,j: expJblock2[i//2,j//2][i%2,j%2])

                blocks.append(exprt * expJblock)
                for i in range(n):
                    vectors.append(re(chain[i]))
                    vectors.append(im(chain[i]))
            else:
                vectors.extend(chain)
                fun = lambda i,j: t**(j-i)/factorial(j-i) if j >= i else 0
                expJblock = Matrix(n, n, fun)
                blocks.append(exp(e * t) * expJblock)

    expJ = Matrix.diag(*blocks)
    P = Matrix(N, N, lambda i,j: vectors[j][i])

    return P, expJ


def _neq_linear_first_order_const_coeff_homogeneous(match_):
    r"""
    System of n first-order constant-coefficient linear homogeneous differential equations

    .. math:: y'_k = a_{k1} y_1 + a_{k2} y_2 +...+ a_{kn} y_n; k = 1,2,...,n

    or that can be written as `\vec{y'} = A . \vec{y}`
    where `\vec{y}` is matrix of `y_k` for `k = 1,2,...n` and `A` is a `n \times n` matrix.

    Since these equations are equivalent to a first order homogeneous linear
    differential equation. So the general solution will contain `n` linearly
    independent parts and solution will consist some type of exponential
    functions. Assuming `y = \vec{v} e^{rt}` is a solution of the system where
    `\vec{v}` is a vector of coefficients of `y_1,...,y_n`. Substituting `y` and
    `y' = r v e^{r t}` into the equation `\vec{y'} = A . \vec{y}`, we get

    .. math:: r \vec{v} e^{rt} = A \vec{v} e^{rt}

    .. math:: r \vec{v} = A \vec{v}

    where `r` comes out to be eigenvalue of `A` and vector `\vec{v}` is the eigenvector
    of `A` corresponding to `r`. There are three possibilities of eigenvalues of `A`

    - `n` distinct real eigenvalues
    - complex conjugate eigenvalues
    - eigenvalues with multiplicity `k`

    1. When all eigenvalues `r_1,..,r_n` are distinct with `n` different eigenvectors
    `v_1,...v_n` then the solution is given by

    .. math:: \vec{y} = C_1 e^{r_1 t} \vec{v_1} + C_2 e^{r_2 t} \vec{v_2} +...+ C_n e^{r_n t} \vec{v_n}

    where `C_1,C_2,...,C_n` are arbitrary constants.

    2. When some eigenvalues are complex then in order to make the solution real,
    we take a linear combination: if `r = a + bi` has an eigenvector
    `\vec{v} = \vec{w_1} + i \vec{w_2}` then to obtain real-valued solutions to
    the system, replace the complex-valued solutions `e^{rx} \vec{v}`
    with real-valued solution `e^{ax} (\vec{w_1} \cos(bx) - \vec{w_2} \sin(bx))`
    and for `r = a - bi` replace the solution `e^{-r x} \vec{v}` with
    `e^{ax} (\vec{w_1} \sin(bx) + \vec{w_2} \cos(bx))`

    3. If some eigenvalues are repeated. Then we get fewer than `n` linearly
    independent eigenvectors, we miss some of the solutions and need to
    construct the missing ones. We do this via generalized eigenvectors, vectors
    which are not eigenvectors but are close enough that we can use to write
    down the remaining solutions. For a eigenvalue `r` with eigenvector `\vec{w}`
    we obtain `\vec{w_2},...,\vec{w_k}` using

    .. math:: (A - r I) . \vec{w_2} = \vec{w}

    .. math:: (A - r I) . \vec{w_3} = \vec{w_2}

    .. math:: \vdots

    .. math:: (A - r I) . \vec{w_k} = \vec{w_{k-1}}

    Then the solutions to the system for the eigenspace are `e^{rt} [\vec{w}],
    e^{rt} [t \vec{w} + \vec{w_2}], e^{rt} [\frac{t^2}{2} \vec{w} + t \vec{w_2} + \vec{w_3}],
    ...,e^{rt} [\frac{t^{k-1}}{(k-1)!} \vec{w} + \frac{t^{k-2}}{(k-2)!} \vec{w_2} +...+ t \vec{w_{k-1}}
    + \vec{w_k}]`

    So, If `\vec{y_1},...,\vec{y_n}` are `n` solution of obtained from three
    categories of `A`, then general solution to the system `\vec{y'} = A . \vec{y}`

    .. math:: \vec{y} = C_1 \vec{y_1} + C_2 \vec{y_2} + \cdots + C_n \vec{y_n}

    """
    eq = match_['eq']
    func = match_['func']
    fc = match_['func_coeff']
    n = len(eq)
    t = list(list(eq[0].atoms(Derivative))[0].atoms(Symbol))[0]
    constants = numbered_symbols(prefix='C', cls=Symbol, start=1)

    # This needs to be modified in future so that fc is only of type Matrix
    M = -fc if type(fc) is Matrix else Matrix(n, n, lambda i,j:-fc[i,func[j],0])

    P, J = matrix_exp_jordan_form(M, t)
    P = simplify(P)
    Cvect = Matrix(list(next(constants) for _ in range(n)))
    sol_vector = P * (J * Cvect)

    sol_vector = [collect(s, ordered(J.atoms(exp)), exact=True) for s in sol_vector]

    sol_dict = [Eq(func[i], sol_vector[i]) for i in range(n)]
    return sol_dict


def _matrix_is_constant(M, t):
    """Checks if the matrix M is independent of t or not."""
    return all(coef.as_independent(t, as_Add=True)[1] == 0 for coef in M)


def _canonical_equations(eqs, funcs, t):
    """Helper function that solves for first order derivatives in a system"""
    from sympy.solvers.solvers import solve

    # For now the system of ODEs dealt by this function can have a
    # maximum order of 1.
    if any(ode_order(eq, func) > 1 for eq in eqs for func in funcs):
        msg = "Cannot represent system in {}-order canonical form"
        raise ODEOrderError(msg.format(1))

    canon_eqs = solve(eqs, *[func.diff(t) for func in funcs], dict=True)

    if len(canon_eqs) != 1:
        raise ODENonlinearError("System of ODEs is nonlinear")

    canon_eqs = canon_eqs[0]
    canon_eqs = [Eq(func.diff(t), canon_eqs[func.diff(t)]) for func in funcs]

    return canon_eqs


def neq_nth_linear_constant_coeff_match(eqs, funcs, t):
    r"""
    Returns a dictionary with details of the eqs if every equation is constant coefficient
    and linear else returns None

    Explanation
    ===========

    This function takes the eqs, converts it into a form Ax = b where x is a vector of terms
    containing dependent variables and their derivatives till their maximum order. If it is
    possible to convert eqs into Ax = b, then all the equations in eqs are linear otherwise
    they are non-linear.

    To check if the equations are constant coefficient, we need to check if all the terms in
    A obtained above are constant or not.

    To check if the equations are homogeneous or not, we need to check if b is a zero matrix
    or not.

    Parameters
    ==========

    eqs: List
        List of ODEs
    funcs: List
        List of dependent variables
    t: Symbol
        Independent variable of the equations in eqs

    Returns
    =======

    match = {
        'no_of_equation': len(eqs),
        'eq': eqs,
        'func': funcs,
        'order': order,
        'is_linear': is_linear,
        'is_constant': is_constant,
        'is_homogeneous': is_homogeneous,
    }

    Dict or None
        Dict with values for keys:
            1. no_of_equation: Number of equations
            2. eq: The set of equations
            3. func: List of dependent variables
            4. order: A dictionary that gives the order of the
                      dependent variable in eqs
            5. is_linear: Boolean value indicating if the set of
                          equations are linear or not.
            6. is_constant: Boolean value indicating if the set of
                          equations have constant coefficients or not.
            7. is_homogeneous: Boolean value indicating if the set of
                          equations are homogeneous or not.
        This Dict is the answer returned if the eqs are linear and constant
        coefficient. Otherwise, None is returned.

    """

    # Error for i == 0 can be added but isn't for now

    # Removing the duplicates from the list of funcs
    # meanwhile maintaining the order. This is done
    # since the line in classify_sysode: list(set(funcs)
    # cause some test cases to fail when gives different
    # results in different versions of Python.
    funcs = list(uniq(funcs))

    # Check for len(funcs) == len(eqs)
    if len(funcs) != len(eqs):
        raise ValueError("Number of functions given is not equal to the number of equations %s" % funcs)

    # ValueError when functions have more than one arguments
    for func in funcs:
        if len(func.args) != 1:
            raise ValueError("dsolve() and classify_sysode() work with "
            "functions of one variable only, not %s" % func)

    # Getting the func_dict and order using the helper
    # function
    order = _get_func_order(eqs, funcs)

    if not all(order[func] == 1 for func in funcs):
        return None
    else:

        # TO be changed when this function is updated.
        # This will in future be updated as the maximum
        # order in the system found.
        system_order = 1

    # Not adding the check if the len(func.args) for
    # every func in funcs is 1

    # Linearity check
    try:
        canon_eqs = _canonical_equations(eqs, funcs, t)
        As, b = linear_ode_to_matrix(canon_eqs, funcs, t, system_order)

    # When the system of ODEs is non-linear, an ODENonlinearError is raised.
    # When system has an order greater than what is specified in system_order,
    # ODEOrderError is raised.
    # This function catches these errors and None is returned
    except (ODEOrderError, ODENonlinearError):
        return None

    A = As[1]
    is_linear = True

    # Constant coefficient check
    is_constant = _matrix_is_constant(A, t)

    # Homogeneous check
    is_homogeneous = True if b.is_zero_matrix else False

    match = {
        'no_of_equation': len(eqs),
        'eq': eqs,
        'func': funcs,
        'order': order,
        'is_linear': is_linear,
        'is_constant': is_constant,
        'is_homogeneous': is_homogeneous,
    }

    # The match['is_linear'] check will be added in the future when this
    # function becomes ready to deal with non-linear systems of ODEs
    if match['is_constant']:

        # Converting the equation into canonical form if the
        # equation is first order. There will be a separate
        # function for this in the future.
        if all([order[func] == 1 for func in funcs]) and match['is_homogeneous']:
            match['func_coeff'] = A
            match['type_of_equation'] = "type1"

            return match

    return None
