from sympy.core import Add, Mul
from sympy.core.containers import Tuple
from sympy.core.compatibility import iterable
from sympy.core.exprtools import factor_terms
from sympy.core.numbers import I
from sympy.core.relational import Eq, Equality
from sympy.core.symbol import Dummy, Symbol
from sympy.core.function import (expand_mul, expand, Derivative,
                            AppliedUndef)
from sympy.functions import exp, im, cos, sin, re, Piecewise, piecewise_fold
from sympy.functions.combinatorial.factorials import factorial
from sympy.matrices import zeros, Matrix, NonSquareMatrixError, MatrixBase
from sympy.polys import Poly
from sympy.simplify import simplify, collect, powsimp, ratsimp
from sympy.solvers.deutils import ode_order
from sympy.solvers.solveset import NonlinearError
from sympy.utilities import numbered_symbols, default_sort_key
from sympy.utilities.iterables import ordered, uniq
from sympy.utilities.misc import filldedent
from sympy.integrals.integrals import Integral, integrate


def _get_func_order(eqs, funcs):
    return {func: max(ode_order(eq, func) for eq in eqs) for func in funcs}


class ODEOrderError(ValueError):
    """Raised by linear_ode_to_matrix if the system has the wrong order"""
    pass


class ODENonlinearError(NonlinearError):
    """Raised by linear_ode_to_matrix if the system is nonlinear"""
    pass


def _simpsol(soleq):
    lhs = soleq.lhs
    sol = soleq.rhs
    sol = powsimp(sol)
    gens = list(sol.atoms(exp))
    p = Poly(sol, *gens, expand=False)
    gens = [factor_terms(g) for g in gens]
    if not gens:
        gens = p.gens
    syms = [Symbol('C1'), Symbol('C2')]
    terms = []
    for coeff, monom in zip(p.coeffs(), p.monoms()):
        coeff = piecewise_fold(coeff)
        if type(coeff) is Piecewise:
            coeff = Piecewise(*((ratsimp(coef).collect(syms), cond) for coef, cond in coeff.args))
        else:
            coeff = ratsimp(coeff).collect(syms)
        monom = Mul(*(g ** i for g, i in zip(gens, monom)))
        terms.append(coeff * monom)
    return Eq(lhs, Add(*terms))


def _solsimp(e, t):
    no_t, has_t = powsimp(expand_mul(e)).as_independent(t)

    no_t = ratsimp(no_t)
    has_t = has_t.replace(exp, lambda a: exp(factor_terms(a)))

    return no_t + has_t


def linodesolve_type(A, t, b=None):
    r"""
    Helper function that determines the type of the system of ODEs for solving with :obj:`sympy.solvers.ode.systems.linodesolve()`

    Explanation
    ===========

    This function takes in the coefficient matrix and/or the non-homogeneous term
    and returns the type of the equation that can be solved by :obj:`sympy.solvers.ode.systems.linodesolve()`.

    If the system is constant coefficient homogeneous, then "type1" is returned
    If the system is constant coefficient non-homogeneous, then "type2" is returned
    If the system is non-constant coefficient homogeneous, then "type3" is returned
    If the system is non-constnt coefficient non-homogeneous, then "type4" is returned

    Note that, if the system of ODEs is of "type3" or "type4", then along with the type,
    the commutative antiderivative of the coefficient matrix is also returned.

    If the system cannot be solved by :obj:`sympy.solvers.ode.systems.linodesolve()`, then
    NotImplementedError is raised.

    Parameters
    ==========

    A : Matrix
        Coefficient matrix of the system of ODEs
    b : Matrix or None
        Non-homogeneous term of the system. The default value is None.
        If this argument is None, then the system is assumed to be homogeneous.

    Examples
    ========

    >>> from sympy import symbols, Matrix
    >>> from sympy.solvers.ode.systems import linodesolve_type
    >>> t = symbols("t")
    >>> A = Matrix([[1, 1], [2, 3]])
    >>> b = Matrix([t, 1])

    >>> linodesolve_type(A, t)
    {'antiderivative': None, 'type': 'type1'}

    >>> linodesolve_type(A, t, b=b)
    {'antiderivative': None, 'type': 'type2'}

    >>> A_t = Matrix([[1, t], [-t, 1]])

    >>> linodesolve_type(A_t, t)
    {'antiderivative': Matrix([
    [      t, t**2/2],
    [-t**2/2,      t]]), 'type': 'type3'}

    >>> linodesolve_type(A_t, t, b=b)
    {'antiderivative': Matrix([
    [      t, t**2/2],
    [-t**2/2,      t]]), 'type': 'type4'}

    >>> A_non_commutative = Matrix([[1, t], [t, -1]])
    >>> linodesolve_type(A_non_commutative, t)
    Traceback (most recent call last):
    ...
    NotImplementedError:
    The system doesn't have a commutative antiderivative, it can't be
    solved by linodesolve.

    Returns
    =======

    Dict

    Raises
    ======

    NotImplementedError
        When the coefficient matrix doesn't have a commutative antiderivative

    See Also
    ========

    linodesolve: Function for which linodesolve_type gets the information

    """

    is_non_constant = not _matrix_is_constant(A, t)
    is_non_homogeneous = not (b is None or b.is_zero_matrix)
    type = "type{}".format(int("{}{}".format(int(is_non_constant), int(is_non_homogeneous)), 2) + 1)

    B = None
    if is_non_constant:
        B, is_commuting = _is_commutative_anti_derivative(A, t)
        if not is_commuting:
            raise NotImplementedError(filldedent('''
                The system doesn't have a commutative antiderivative, it can't be solved
                by linodesolve.
            '''))

    return {"type": type, "antiderivative": B}


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

        Ai = Ai.applyfunc(expand_mul)

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


def linodesolve(A, t, b=None, B=None, type="auto", doit=False, const_idx=0):
    r"""
    System of n equations linear first-order differential equations

    Explanation
    ===========

    This solver solves the system of ODEs of the follwing form:

    .. math::
        X'(t) = A(t) X(t) +  b(t)

    Here, $A(t)$ is the coefficient matrix, $X(t)$ is the vector of n independent variables,
    $b(t)$ is the non-homogeneous term and $X'(t)$ is the derivative of $X(t)$

    Depending on the properties of $A(t)$ and $b(t)$, this solver evaluates the solution
    differently.

    When $A(t)$ is constant coefficient matrix and $b(t)$ is zero vector i.e. system is homogeneous,
    the solution is:

    .. math::
        X(t) = \exp(A t) C

    Here, $C$ is a vector of constants and $A$ is the constant coefficient matrix.

    When $A(t)$ is constant coefficient matrix and $b(t)$ is non-zero i.e. system is non-homogeneous,
    the solution is:

    .. math::
        X(t) = e^{A t} ( \int e^{- A t} b \,dt + C)

    When $A(t)$ is coefficient matrix such that its commutative with its antiderivative $B(t)$ and
    $b(t)$ is a zero vector i.e. system is homogeneous, the solution is:

    .. math::
        X(t) = \exp(B(t)) C

    When $A(t)$ is commutative with its antiderivative $B(t)$ and $b(t)$ is non-zero i.e. system is
    non-homogeneous, the solution is:

    .. math::
        X(t) =  e^{B(t)} ( \int e^{-B(t)} b(t) \,dt + C)

    The final solution is the general solution for all the four equations since a constant coefficient
    matrix is always commutative with its antidervative.

    Parameters
    ==========

    A : Matrix
        Coefficient matrix of the system of linear first order ODEs.
    t : Symbol
        Independent variable in the system of ODEs.
    b : Matrix or None
        Non-homogeneous term in the system of ODEs. If None is passed,
        a homogeneous system of ODEs is assumed.
    B : Matrix or None
        Antiderivative of the coefficient matrix. If the antiderivative
        is not passed and the solution requires the term, then the solver
        would compute it internally.
    type : String
        Type of the system of ODEs passed. Depending on the type, the
        solution is evaluated. The type values allowed and the corresponding
        system it solves are: "type1" for constant coefficient homogeneous
        "type2" for constant coefficient non-homogeneous, "type3" for non-constant
        coefficient homogeneous and "type4" for non-constant coefficient non-homogeneous.
        The default value is "auto" which will let the solver decide the correct type of
        the system passed.
    doit : Boolean
        Evaluate the solution if True, default value is False

    Examples
    ========

    To solve the system of ODEs using this function directly, several things must be
    done in the right order. Wrong inputs to the function will lead to incorrect results.

    >>> from sympy import symbols, Function, Eq
    >>> from sympy.solvers.ode.systems import canonical_odes, linear_ode_to_matrix, linodesolve, linodesolve_type
    >>> from sympy.solvers.ode.subscheck import checksysodesol
    >>> f, g = symbols("f, g", cls=Function)
    >>> x, a = symbols("x, a")
    >>> funcs = [f(x), g(x)]
    >>> eqs = [Eq(f(x).diff(x) - f(x), a*g(x) + 1), Eq(g(x).diff(x) + g(x), a*f(x))]

    Here, it is important to note that before we derive the coefficient matrix, it is
    important to get the system of ODEs into the desired form. For that we will use
    :obj:`sympy.solvers.ode.systems.canonical_odes()`.

    >>> eqs = canonical_odes(eqs, funcs, x)
    >>> eqs
    [Eq(Derivative(f(x), x), a*g(x) + f(x) + 1), Eq(Derivative(g(x), x), a*f(x) - g(x))]

    Now, we will use :obj:`sympy.solvers.ode.systems.linear_ode_to_matrix()` to get the coefficient matrix and the
    non-homogeneous term if it is there.

    >>> (A1, A0), b = linear_ode_to_matrix(eqs, funcs, x, 1)
    >>> A = -A0

    We have the coefficient matrices and the non-homogeneous term ready. Now, we can use
    :obj:`sympy.solvers.ode.systems.linodesolve_type()` to get the information for the system of ODEs
    to finally pass it to the solver.

    >>> system_info = linodesolve_type(A, x, b=b)
    >>> sol_vector = linodesolve(A, x, b=b, B=system_info['antiderivative'], type=system_info['type'])

    Now, we can prove if the solution is correct or not by using :obj:`sympy.solvers.ode.subscheck.checksysodesol()`

    >>> sol = [Eq(f, s) for f, s in zip(funcs, sol_vector)]
    >>> checksysodesol(eqs, sol)
    (True, [0, 0])

    We can also use the doit method to evaluate the solutions passed by the function.

    >>> sol_vector_evaluated = linodesolve(A, x, b=b, type="type2", doit=True)

    Now, we will look at a system of ODEs which is non-constant.

    >>> eqs = [Eq(f(x).diff(x), f(x) + x*g(x)), Eq(g(x).diff(x), -x*f(x) + g(x))]

    The system defined above is already in the desired form, so we don't have to convert it.

    >>> (A1, A0), b = linear_ode_to_matrix(eqs, funcs, x, 1)
    >>> A = -A0

    A user can also pass the commutative antidervative required for type3 and type4 system of ODEs.
    Passing an incorrect one will lead to incorrect results. If the coefficient matrix is not commutative
    with its antiderivative, then :obj:`sympy.solvers.ode.systems.linodesolve_type()` raises a NotImplementedError.
    If it does have a commutative antiderivative, then the function just returns the information about the system.

    >>> system_info = linodesolve_type(A, x, b=b)

    Now, we can pass the antiderivative as an argument to get the solution. If the system information is not
    passed, then the solver will compute the required arguments internally.

    >>> sol_vector = linodesolve(A, x, b=b)

    Once again, we can verify the solution obtained.

    >>> sol = [Eq(f, s) for f, s in zip(funcs, sol_vector)]
    >>> checksysodesol(eqs, sol)
    (True, [0, 0])

    Returns
    =======

    List

    Raises
    ======

    ValueError
        This error is raised when the coefficient matrix, non-homogeneous term
        or the antiderivative, if passed, aren't a matrix or
        don't have correct dimensions
    NonSquareMatrixError
        When the coefficient matrix or its antiderivative, if passed isn't a square
        matrix
    NotImplementedError
        If the coefficient matrix doesn't have a commutative antiderivative

    See Also
    ========

    linear_ode_to_matrix: Coefficient matrix computation function
    canonical_odes: System of ODEs representation change
    linodesolve_type: Getting information about systems of ODEs to pass in this solver

    """

    if not isinstance(A, MatrixBase):
        raise ValueError(filldedent('''\
            The coefficients of the system of ODEs should be of type Matrix
        '''))

    if not A.is_square:
        raise NonSquareMatrixError(filldedent('''\
            The coefficient matrix must be a square
        '''))

    if b is not None:
        if not isinstance(b, MatrixBase):
            raise ValueError(filldedent('''\
                The non-homogeneous terms of the system of ODEs should be of type Matrix
            '''))

        if A.rows != b.rows:
            raise ValueError(filldedent('''\
                The system of ODEs should have the same number of non-homogeneous terms and the number of
                equations
            '''))

    if B is not None:
        if not isinstance(B, MatrixBase):
            raise ValueError(filldedent('''\
                The antiderivative of coefficients of the system of ODEs should be of type Matrix
            '''))

        if not B.is_square:
            raise NonSquareMatrixError(filldedent('''\
                The antiderivative of the coefficient matrix must be a square
            '''))

        if A.rows != B.rows:
            raise ValueError(filldedent('''\
                        The coefficient matrix and its antiderivative should have same dimensions
                    '''))

    if not any(type == "type{}".format(i) for i in range(1, 5)) and not type == "auto":
        raise ValueError(filldedent('''\
                    The input type should be a valid one
                '''))

    n = A.rows

    constants = numbered_symbols(prefix='C', cls=Symbol, start=const_idx+1)
    Cvect = Matrix(list(next(constants) for _ in range(n)))

    if (type == "type2" or type == "type4") and b is None:
        b = zeros(n, 1)

    if type == "auto":
        system_info = linodesolve_type(A, t, b=b)
        type = system_info["type"]
        B = system_info["antiderivative"]

    if type == "type1" or type == "type2":
        P, J = matrix_exp_jordan_form(A, t)
        P = simplify(P)

        if type == "type1":
            sol_vector = P * (J * Cvect)
        else:
            sol_vector = P * J * ((J.inv() * P.inv() * b).applyfunc(lambda x: Integral(x, t)) + Cvect)

    else:
        if B is None:
            B, _ = _is_commutative_anti_derivative(A, t)

        if type == "type3":
            sol_vector = B.exp() * Cvect
        else:
            sol_vector = B.exp() * (((-B).exp() * b).applyfunc(lambda x: Integral(x, t)) + Cvect)

    gens = sol_vector.atoms(exp)

    if type != "type1":
        sol_vector = [expand_mul(s) for s in sol_vector]

    sol_vector = [collect(s, ordered(gens), exact=True) for s in sol_vector]

    if doit:
        sol_vector = [s.doit() for s in sol_vector]

    return sol_vector


def _matrix_is_constant(M, t):
    """Checks if the matrix M is independent of t or not."""
    return all(coef.as_independent(t, as_Add=True)[1] == 0 for coef in M)


def canonical_odes(eqs, funcs, t):
    r"""
    Function that solves for highest order derivatives in a system

    Explanation
    ===========

    This function inputs a system of ODEs and based on the system,
    the dependent variables and their highest order, returns the system
    in the following form:

    .. math::
        X'(t) = A(t) X(t) + b(t)

    Here, $X(t)$ is the vector of dependent variables of lower order, $A(t)$ is
    the coefficient matrix, $b(t)$ is the non-homogeneous term and $X'(t)$ is the
    vector of dependent variables in their respective highest order. We use the term
    canonical form to imply the system of ODEs which is of the above form.

    If the system passed has a non-linear term with multiple solutions, then a list of
    systems is returned in its canonical form.

    Parameters
    ==========

    eqs : List
        List of the ODEs
    funcs : List
        List of dependent variables
    t : Symbol
        Independent variable

    Examples
    ========

    >>> from sympy import symbols, Function, Eq, Derivative
    >>> from sympy.solvers.ode.systems import canonical_odes
    >>> f, g = symbols("f g", cls=Function)
    >>> x, y = symbols("x y")
    >>> funcs = [f(x), g(x)]
    >>> eqs = [Eq(f(x).diff(x) - 7*f(x), 12*g(x)), Eq(g(x).diff(x) + g(x), 20*f(x))]

    >>> canonical_eqs = canonical_odes(eqs, funcs, x)
    >>> canonical_eqs
    [Eq(Derivative(f(x), x), 7*f(x) + 12*g(x)), Eq(Derivative(g(x), x), 20*f(x) - g(x))]

    >>> system = [Eq(Derivative(f(x), x)**2 - 2*Derivative(f(x), x) + 1, 4), Eq(-y*f(x) + Derivative(g(x), x), 0)]

    >>> canonical_system = canonical_odes(system, funcs, x)
    >>> canonical_system
    [[Eq(Derivative(f(x), x), -1), Eq(Derivative(g(x), x), y*f(x))], [Eq(Derivative(f(x), x), 3), Eq(Derivative(g(x), x), y*f(x))]]

    Returns
    =======

    List

    """
    from sympy.solvers.solvers import solve

    order = _get_func_order(eqs, funcs)

    canon_eqs = solve(eqs, *[func.diff(t, order[func]) for func in funcs], dict=True)

    systems = []
    for eq in canon_eqs:
        system = [Eq(func.diff(t, order[func]), eq[func.diff(t, order[func])]) for func in funcs]
        systems.append(system)

    if len(canon_eqs) == 1:
        systems = systems[0]

    return systems


def _is_commutative_anti_derivative(A, t):
    r"""
    Helper function for determining if the Matrix passed is commutative with its antiderivative

    Explanation
    ===========

    This function checks if the Matrix $A$ passed is commutative with its antiderivative with respect
    to the independent variable $t$.

    .. math::
        B(t) = \int A(t) dt

    The function outputs two values, first one being the antiderivative $B(t)$, second one being a
    boolean value, if True, then the matrix $A(t)$ passed is commutative with $B(t)$, else the matrix
    passed isn't commutative with $B(t)$.

    Parameters
    ==========

    A : Matrix
        The matrix which has to be checked
    t : Symbol
        Independent variable

    Examples
    ========

    >>> from sympy import symbols, Matrix
    >>> from sympy.solvers.ode.systems import _is_commutative_anti_derivative
    >>> t = symbols("t")
    >>> A = Matrix([[1, t], [-t, 1]])

    >>> B, is_commuting = _is_commutative_anti_derivative(A, t)
    >>> is_commuting
    True

    Returns
    =======

    Matrix, Boolean

    """
    B = integrate(A, t)
    is_commuting = (B*A - A*B).applyfunc(expand).applyfunc(factor_terms).is_zero_matrix

    is_commuting = False if is_commuting is None else is_commuting

    return B, is_commuting


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

    Dict or list of Dicts or None
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
            8. commutative_antiderivative: Antiderivative of the coefficient
                          matrix if the coefficient matrix is non-constant
                          and commutative with its antiderivative. This key
                          may or may not exist.
            9. is_general: Boolean value indicating if the system of ODEs is
                           solvable using one of the general case solvers or not.
            10. rhs: rhs of the non-homogeneous system of ODEs in Matrix form. This
                     key may or may not exist.
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
        canon_eqs = canonical_odes(eqs, funcs, t)

        if isinstance(canon_eqs[0], Eq):
            As, b = linear_ode_to_matrix(canon_eqs, funcs, t, system_order)
            A = As[1]
        else:

            match = {
                'is_implicit': True,
                'canon_eqs': canon_eqs
            }

            return match

    # When the system of ODEs is non-linear, an ODENonlinearError is raised.
    # This function catches the error and None is returned.
    except ODENonlinearError:
        return None

    is_linear = True

    # Constant coefficient check
    is_constant = _matrix_is_constant(A, t)

    # Homogeneous check
    is_homogeneous = True if b.is_zero_matrix else False

    # Is general key is used to identify if the system of ODEs can be solved by
    # one of the general case solvers or not.
    match = {
        'no_of_equation': len(eqs),
        'eq': eqs,
        'func': funcs,
        'order': order,
        'is_linear': is_linear,
        'is_constant': is_constant,
        'is_homogeneous': is_homogeneous,
        'is_general': True
    }

    # The match['is_linear'] check will be added in the future when this
    # function becomes ready to deal with non-linear systems of ODEs

    if all([order[func] == 1 for func in funcs]):
        match['func_coeff'] = A

        if not is_homogeneous:
            match['rhs'] = b

        try:
            system_info = linodesolve_type(-A, t, b=b)
        except NotImplementedError:
            return None

        if not is_constant:
            match['commutative_antiderivative'] = system_info["antiderivative"]

        match['type_of_equation'] = system_info["type"]

        return match

    return None


def _preprocess_eqs(eqs):
    processed_eqs = []
    for eq in eqs:
        processed_eqs.append(eq if isinstance(eq, Equality) else Eq(eq, 0))

    return processed_eqs


# For now, this function returns a simple output, later it will
# be capable of dividing the system of ODEs into strongly and
# weakly connected components.
def _component_division(eqs, funcs):
    return [[[eqs, funcs]]]


# Returns: List of equations
def _linear_ode_solver(match):
    t = match['t']
    funcs = match['func']

    rhs = match.get('rhs', None)
    A = -match['func_coeff']
    B = match.get('commutative_antiderivative', None)
    const_idx = match['const_idx']
    type_of_equation = match['type_of_equation']

    sol_vector = linodesolve(A, t, b=rhs, B=B,
                             type=type_of_equation, const_idx=const_idx)

    sol = [Eq(f, s) for f, s in zip(funcs, sol_vector)]

    return sol


# Returns: (List of List of equations, const_idx) or None
def _ode_component_solver(eqs, funcs, t, const_idx=0):
    match = neq_nth_linear_constant_coeff_match(eqs, funcs, t)

    if match:
        match['const_idx'] = const_idx
        match['t'] = t

        if match.get('is_linear', False):
            return [_linear_ode_solver(match)], const_idx + len(eqs)
        if match.get('is_implicit', False):
            canon_eqs = match['canon_eqs']
            sols = []
            temp_const_idx = 0
            for canon_eq in canon_eqs:
                sol = _ode_component_solver(canon_eq, funcs, t, const_idx=const_idx)
                if sol is not None:
                    sol, temp_const_idx = sol
                    sols += sol

            # Note: This const_idx logic has to be
            # verified
            return sols, temp_const_idx

        # To add non-linear case here in future

    return None


def _weak_component_solver(wcc, t, const_idx):
    loop_sol = []
    not_solved_systems = []

    for j, scc in enumerate(wcc):
        eqs, funcs = scc
        scc_sols = []

        # If this is the first SCC to be solved,
        # then looping isn't necessary.
        if not loop_sol:

            # Solving the SCC
            comp_sol = _ode_component_solver(eqs, funcs, t, const_idx=const_idx)
            if comp_sol is not None:

                # scc_sol: List of List of Equations
                # scc_sol is the list of solutions
                scc_sol, const_idx = comp_sol

                # scc_sols: List of List of List of Equations
                # Since only one scc_sol is appended in this
                # block, the scc_sols will look like:
                # scc_sol = [[s1, s2, ...]]
                # where s1, s2, ... are solutions. A solution
                # is a List of Equations.
                scc_sols.append(scc_sol)
        else:

            changed_const_idx = 0
            for sol in loop_sol:

                # If this is a SCC after the first one, then we will
                # loop through loop_sol, which is basically List of
                # solutions, and substitute each solution to solve the
                # SCC. Each time we substitute a different solution from
                # loop_sol, we will get a different subsystem, and different
                # solution.
                comp_eqs = eqs.subs({s.lhs: s.rhs for s in sol})
                comp_sol = _ode_component_solver(comp_eqs, funcs, t, const_idx=const_idx)

                if comp_sol is None:
                    continue

                # scc_sol: List of List of equations
                # scc_sol is List of solutions
                scc_sol, changed_const_idx = comp_sol

                # scc_sols: List of List of List of equations
                # If we got s1, s2 as solutions for first iteration
                # of this loop and s3, s4 for second iteration, then
                # scc_sol = [[s1, s2], [s3, s4]]
                # where each solution is a List of Equations.
                scc_sols.append(scc_sol)

            # Note: This const_idx logic has to be
            # verified
            const_idx = changed_const_idx

        if not scc_sols:

            # All the SCCs which are not solvable
            # If a SCC is found to be unsolvable,
            # then all the other SCCs won't have a
            # solution and we would have to stop
            # the loop.
            not_solved_systems += wcc[j:]
            break

        # Its time to combine the solutions, there are two cases:
        # 1. When loop_sol isn't empty: This would mean we have
        #    solved a SCC before.
        # 2. When loop_sol is empty: This is first SCC we solved.

        # Case 1: loop_sol isn't empty.
        if loop_sol:

            # For each solution in loop_sol, we have to append the
            # appropriate solution from scc_sols
            # For example: loop_sol has two solutions s1 and s2
            # Now, by substituting s1, if we got s3 and s4 and by
            # substituting s2, we got s5 and s6, then the combined
            # solution for the weakly connected component till now
            # will be: [s1+s3, s1+s4, s2+s5, s2+s6]
            loop_sol = [sol + scc_s for sol, scc_sol in zip(loop_sol, scc_sols) for scc_s in scc_sol]
        else:

            # If this is the first SCC to be solved, then
            # scc_sols will be of the form:
            # scc_sols = [[s1, s2, ...]]
            # Hence, our starting loop_sol will be
            # loop_sol = [s1, s2, ...]
            loop_sol = scc_sols[0]

    return loop_sol, const_idx, not_solved_systems


def dsolve_system(eqs, funcs=None, t=None, ics=None):
    r"""
    Solves any(supported) system of Ordinary Differential Equations

    Explanation
    ===========

    This function takes a system of ODEs as an input, determines if the
    it is solvable by this function, and returns the solution if found any.

    This function can handle:
    1. Linear, First Order, Constant coefficient homogeneous system of ODEs
    2. Linear, First Order, Constant coefficient non-homogeneous system of ODEs
    3. Linear, First Order, non-constant coefficient homogeneous system of ODEs
    4. Linear, First Order, non-constant coefficient non-homogeneous system of ODEs
    5. Any implicit system which can be divided into system of ODEs which is of the above 4 forms

    The types of systems described above aren't limited by the number of equations, i.e. this
    function can solve the above types irrespective the number of equations in the system passed.

    This function returns a pair of elements, first being the solutions obtained and second being
    the sub-system of ODEs that it was not able to solve. If the function wasn't able to solve the
    whole system, then just a tuple of two empty lists is returned.

    Parameters
    ==========

    eqs : List
        List of ODEs to be solved
    funcs : List or None
        List of dependent variables that make up the system of ODEs
    t : Symbol
        Independent variable in the system of ODEs
    ics : Dict or None
        Set of initial boundary/conditions for the system of ODEs.

    Examples
    ========

    >>> from sympy import symbols, Eq, Function
    >>> from sympy.solvers.ode.systems import dsolve_system
    >>> f, g = symbols("f g", cls=Function)
    >>> x = symbols("x")

    >>> eqs = [Eq(f(x).diff(x), g(x)), Eq(g(x).diff(x), f(x))]
    >>> dsolve_system(eqs)
    ([[Eq(f(x), -C1*exp(-x) + C2*exp(x)), Eq(g(x), C1*exp(-x) + C2*exp(x))]], [])

    You can also pass the initial conditions for the system of ODEs:
    >>> dsolve_system(eqs, ics={f(0): 1, g(0): 0})
    ([[Eq(f(x), exp(x)/2 + exp(-x)/2), Eq(g(x), exp(x)/2 - exp(-x)/2)]], [])

    Optionally, you can pass the dependent variables and the independent
    variable for which the system is to be solved:
    >>> funcs = [f(x), g(x)]
    >>> dsolve_system(eqs, funcs=funcs, t=x)
    ([[Eq(f(x), -C1*exp(-x) + C2*exp(x)), Eq(g(x), C1*exp(-x) + C2*exp(x))]], [])

    Lets look at an implicit system of ODEs:
    >>> eqs = [Eq(f(x).diff(x)**2, g(x)**2), Eq(g(x).diff(x), g(x))]
    >>> dsolve_system(eqs)
    ([[Eq(f(x), C1 - C2*exp(x)), Eq(g(x), C2*exp(x))], [Eq(f(x), C1 + C2*exp(x)), Eq(g(x), C2*exp(x))]], [])

    Returns
    =======

    Tuple : (List of List of Equations, List of sub-systems not solved)

    """
    from sympy.solvers.ode.ode import solve_ics

    if not iterable(eqs):
        raise ValueError(filldedent('''
            List of equations should be passed. The input is not valid.
        '''))

    eqs = _preprocess_eqs(eqs)

    if funcs is not None and not isinstance(funcs, list):
        raise ValueError(filldedent('''
            Input to the funcs should be a list of functions.
        '''))

    # Note: This is added to solve a major problem encountered.
    # Might be best to make a function for this functions
    # extraction later.
    if funcs is None:
        funcs = []
        for eq in eqs:
            derivs = eq.atoms(Derivative)
            func = set().union(*[d.atoms(AppliedUndef) for d in derivs])
            for func_ in func:
                funcs.append(func_)
        funcs = list(uniq(funcs))

    if len(eqs) != len(funcs):
        raise ValueError(filldedent('''
            Number of equations and number of functions don't match
        '''))

    if t is not None and not isinstance(t, Symbol):
        raise ValueError(filldedent('''
            The indepedent variable must be of type Symbol
        '''))

    if t is None:
        t = list(list(eqs[0].atoms(Derivative))[0].atoms(Symbol))[0]

    sols = []
    const_idx = 0
    components = _component_division(eqs, funcs)

    not_solved_systems = []

    for wcc in components:

        # wcc_sol: List of List of Equations
        wcc_sol, const_idx, wcc_not_solved_systems = _weak_component_solver(wcc, t, const_idx)

        not_solved_systems += wcc_not_solved_systems

        if wcc_sol:

            # If there are solutions for other weakly
            # connected components, then the solution
            # for the combined system will be all the
            # combinations of solutions of other weakly
            # connected components and this one's solutions.
            if sols:
                sols = [s + ls for s in sols for ls in wcc_sol]

            else:

                # If this is the first weakly connected component
                # to be solved, then the solution will start with
                # this weakly connected component.
                sols = wcc_sol

    if ics and sols:
        ics_sols = []

        # This is assuming that all the solutions
        # have the same funcs. This may have to
        # be changed when system division is
        # added.
        funcs = [s.lhs for s in sols[0]]
        for sol in sols:
            constants = Tuple(*sol).free_symbols - Tuple(*eqs).free_symbols
            solved_constants = solve_ics(sol, funcs, constants, ics)
            ics_sols.append([s.subs(solved_constants) for s in sol])

        sols = ics_sols

    if not sols:
        not_solved_systems = []

    # sols: List of List of Equations
    return sols, not_solved_systems
