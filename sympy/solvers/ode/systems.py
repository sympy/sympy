from sympy import (Matrix, Derivative, Symbol)
from sympy.solvers.deutils import ode_order
from sympy.core.numbers import I
from sympy.core.relational import Eq
from sympy.core.symbol import Dummy
from sympy.functions import exp, im, cos, sin, re
from sympy.functions.combinatorial.factorials import factorial
from sympy.matrices import zeros
from sympy.utilities import numbered_symbols, default_sort_key
from sympy.utilities.iterables import uniq
from sympy.simplify import simplify


def _get_func_order(eqs, funcs):
    order = {}

    for func in funcs:
        if not order.get(func, False):
            max_order = -1
            for i, eq in enumerate(eqs):
                max_order = max(max_order, ode_order(eq, func))
            if max_order >= 0:
                order[func] = max_order

    return order


def matrix_exp(A, t):
    '''Matrix exponential exp(A*t) for the matrix A and scalar t.

    The matrix exponential exp(A*t) appears in the solution of linear
    differential equations. For example if x is a vector and A is a matrix
    then the initial value problem

        dx(t)/dt = A x(t),   x(0) = x0

    has the unique solution

        x(t) = exp(A t) x0

    Example:

    >>> from sympy import Symbol, Matrix, pprint
    >>> from sympy.solvers.ode.systems import matrix_exp
    >>> t = Symbol('t')
    >>> A = Matrix([[2, -5], [2, -4]])
    >>> pprint(A)
    [2  -5]
    [     ]
    [2  -4]
    >>> pprint(matrix_exp(A, t))
    [   -t           -t                    -t              ]
    [3*e  *sin(t) + e  *cos(t)         -5*e  *sin(t)       ]
    [                                                      ]
    [         -t                     -t           -t       ]
    [      2*e  *sin(t)         - 3*e  *sin(t) + e  *cos(t)]
    '''
    P, J = matrix_exp_jordan_form(A, t)
    return P * J * P.inv()


def matrix_exp_jordan_form(A, t):
    '''Matrix exponential exp(A*t) for the matrix A and scalar t.

    Returns the Jordan form of the exp(A t).

    >>> from sympy import Matrix, Symbol
    >>> from sympy.solvers.ode.systems import matrix_exp, matrix_exp_jordan_form
    >>> t = Symbol('t')
    >>> A = Matrix([[1, 1], [0, 1]])
    >>> P, J = matrix_exp_jordan_form(A, t)
    >>> P * J * P.inv() == matrix_exp(A, t)
    True
    '''

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
    eigenchains = sorted(eigenchains.items(), key=default_sort_key)
    isreal = not A.has(I)

    blocks = []
    vectors = []
    seen_conjugate = set()
    for e, chains in eigenchains:
        for chain in chains:
            n = len(chain)
            if isreal and im(e).is_nonzero:
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

    M = -fc if type(fc) is Matrix else Matrix(n, n, lambda i,j:-fc[i,func[j],0])
    P, J = matrix_exp_jordan_form(M, t)
    P = simplify(P)
    Cvect = Matrix(list(next(constants) for _ in range(n)))
    sol_vector = P * J * Cvect

    sol_dict = [Eq(func[i], sol_vector[i]) for i in range(n)]
    return sol_dict


def _matrix_is_constant(M, t):
    """Checks if the matrix M is independent of t or not."""
    return all(coef.as_independent(t, as_Add=True)[1] == 0 for coef in M)


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
    from sympy.solvers.solveset import linear_eq_to_matrix
    from sympy.solvers.solvers import solve

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

    # Not adding the check if the len(func.args) for
    # every func in funcs is 1

    rep = {func.diff(t, n): Dummy() for func in order.keys() for n in range(order[func] + 1)}
    eqs_sub = [eq.subs(rep) for eq in eqs]
    # Linearity check
    try:
        A, b = linear_eq_to_matrix(eqs_sub, rep.values())

    except ValueError:
        return None

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
            canon_eqs = solve(eqs, *[func.diff(t) for func in funcs])
            canon_eqs = [func.diff(t) - canon_eqs[func.diff(t)] for func in funcs]
            new_eqs = [canon_eq.subs(rep) for canon_eq in canon_eqs]
            coef = linear_eq_to_matrix(new_eqs, [rep[func] for func in funcs])
            match['func_coeff'] = coef[0]
            match['type_of_equation'] = "type1"

            return match

    return None
