from sympy.functions.combinatorial.factorials import factorial
from sympy.core.symbol import symbols
from sympy.solvers.solveset import linsolve
from sympy.utilities import public

@public
def pade_approximant(expr, x, x0, n, m=None):
    """
    Return the [n/m] Pade approximation p(x)/q(x) of expr around ``x = x0'',
    where p(x) and q(x) are polynomials of order n and m respectively.

    The polynomials p(x) and q(x) are found by requiring the taylor series of
    expr and p(x)/q(x) to match up to (n + m)th order.

    Parameters
    ----------
    expr : Expression
           The expression whose Pade approximant is to be computed.

    x   : Symbol
        It is the variable of the expression to be calculated.

    x0  : Value
         The value around which ``x`` is calculated. Can be any value
         from ``-oo`` to ``oo``.

    n   : int
        The order of the returned approximating polynomial `p`.

    m   : int, optional
        The order of the returned approximating polynomial `q`. Defaults to n.

    Returns
    -------
    p(x)/q(x) : Rational Function
        The Pade approximation of expr around ``x = x0``.

    Examples
    --------
    >>> import sympy as sp
    >>> x = sp.symbols('x')
    >>> exp_pade = pade_approximant(sp.exp(x), x, 0, 2, 2)
    >>> exp_pade
    (x**2/12 + x/2 + 1)/(x**2/12 - x/2 + 1)

    Compare the [3/4] Pade approximant of sin(x) around x=0 with sin(x)
    >>> sin_pade = pade_approximant(sp.sin(x), x, 0, 3, 4)
    >>> sin_pade
    (-31*x**3/294 + x)/(11*x**4/5880 + 3*x**2/49 + 1)

    >>> sp.N(sin_pade.subs(x, 1))
    0.841465365541513

    >>> sp.N(sp.sin(1))
    0.841470984807897

    """

    if n < 0:
        raise ValueError("Order of p <n> must be greater than 0.")
    if m is None:
        m = n
    elif m < 0:
        raise ValueError("Order of q <m> must be greater than 0.")

    N = m + n
    taylor_coefficients = [
        expr.diff(x, i).subs(x, x0) / factorial(i) for i in range(N + 1)
    ]

    pq_coefficients = symbols(f"p0:{N + 1}")

    A = [[1 if i == j else 0 for j in range(n + 1)] for i in range(N + 1)]
    B = [[0 for _ in range(m)] for _ in range(N + 1)]

    for row in range(1, m + 1):
        B[row][:row] = [-f for f in taylor_coefficients[:row][::-1]]
    for row in range(m + 1, N + 1):
        B[row][:] = [-f for f in taylor_coefficients[row - m : row][::-1]]

    C = [a + b for a, b in zip(A, B)]
    linear_system = [
        sum([c * pq for c, pq in zip(row, pq_coefficients)]) - f
        for row, f in zip(C, taylor_coefficients)
    ]

    pq_solution = list(*linsolve(linear_system, *pq_coefficients))

    # the linear system has no solution iff the leading terms in both p(x) and q(x) are zero
    # in this case, the pade approximation [n/m] is the same as the pade approximation [n-1/m-1]
    if pq_solution == []:
        return pade_approximant(expr, x, x0, n - 1, m - 1)

    p = sum([p * (x - x0) ** i for i, p in enumerate(pq_solution[: n + 1])])
    q = 1 + sum(
        [q * (x - x0) ** (i + 1) for i, q in enumerate(pq_solution[n + 1 :])]
    )

    return p / q
