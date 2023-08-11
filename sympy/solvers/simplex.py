"""Tools for optimizing a linear function for a given simplex.

The minimization of a linear objective, ``f = c*x - d``, with constraints
``A*x >= b`` or b) the maximization of ``f`` with constraints ``A*x <= b``
can be solved with calls to `lpmin` or `lpmax`, in matrix form or
in symbolic form (a function and a list of constraints).

The primal and dual corresponding to a given matrix for the
standard minimization can be generated with `_primal_dual`.

Constraints that are univariate will affect the range of values
returned by the optimization, e.g. ``x <= 3`` will permit negative
values of x unless there is a corresponding ``x >= 0`` condition.
If there are no univariate conditions, only nonnegative
solutions will be found.
"""

from sympy.core import sympify
from sympy.core.exprtools import factor_terms
from sympy.core.relational import Lt, Gt, Le, Ge, Eq
from sympy.core.symbol import Dummy
from sympy.core.singleton import S
from sympy.core.sorting import ordered
from sympy.functions.elementary.complexes import sign
from sympy.matrices.dense import Matrix
from sympy.matrices.matrices import MatrixBase
from sympy.solvers.solveset import linear_eq_to_matrix
from sympy.utilities.iterables import numbered_symbols
from sympy.utilities.misc import filldedent

class UnboundedLPError(Exception):
    """
    A linear programing problem is said to be unbounded if its objective
    function can assume arbitrarily large values.

    Example
    =======

    Suppose you want to maximize
        2x
    subject to
        x >= 0

    There's no upper limit that 2x can take.
    """
    pass


class InfeasibleLPError(Exception):
    """
    A linear programing problem is considered infeasible if its
    constraint set is empty. That is, if the set of all vectors
    satisfying the contraints is empty, then the problem is infeasible.

    Example
    =======

    Suppose you want to maximize
        x
    subject to
        x >= 10
        x <= 9

    No x can satisfy those constraints.
    """
    pass


def _pivot(M, i, j):
    """
    The pivot element `M[i, j]` is inverted and the rest of the matrix
    modified and returned as a new matrix; original is left unmodified.

    Example
    =======

    >>> from sympy.matrices.dense import Matrix
    >>> from sympy.solvers.simplex import _pivot
    >>> from sympy import var
    >>> Matrix(3, 3, var('a:i'))
    Matrix([
    [a, b, c],
    [d, e, f],
    [g, h, i]])
    >>> _pivot(_, 1, 0)
    Matrix([
    [-a/d, -a*e/d + b, -a*f/d + c],
    [ 1/d,        e/d,        f/d],
    [-g/d,  h - e*g/d,  i - f*g/d]])
    """
    Mi, Mj, Mij = M[i,:], M[:,j], M[i,j]
    if Mij == 0:
        raise ZeroDivisionError(
            "Tried to pivot about zero-valued entry.")
    A = M - Mj * (Mi / Mij)
    A[i, :] = Mi / Mij
    A[:, j] = -Mj / Mij
    A[i, j] = 1 / Mij
    return A


def _choose_pivot_row(A, B, candidate_rows, pivot_col, Y):
    # Choose row with smallest ratio
    first_row = candidate_rows[0]
    min_ratio = B[first_row] / A[first_row, pivot_col]
    min_rows = [first_row]
    for i in candidate_rows[1:]:
        ratio = B[i] / A[i, pivot_col]
        if ratio < min_ratio:
            min_ratio = ratio
            min_rows = [i]
        elif ratio == min_ratio:
            min_rows.append(i)

    # If there are ties, pick using Bland's rule
    row = sorted(min_rows, key= lambda r: Y[r])[0]
    return row


def _simplex(A, B, C, D=None, dual=False):
    """
    Return ``(o, x, y)`` obtained from the two-phase simplex method
    using Bland's rule where ``o`` is the minimum value of primal,
    ``Cx - D``, under constraints ``Ax <= B`` (with ``x >= 0``) and
    the maximum of the dual, ``y^{T}B - D``, under constraints
    ``A^{T}*y >= C^{T}`` (with ``y >= 0``). To compute the dual of
    the system, pass `dual=True` and ``(o, y, x)`` will be returned.

    Note: the nonnegative constraints for ``x`` and ``y`` supercede
    any values of ``A`` and ``B`` that are inconsistent with that
    assumption, so if a constraint of ``x >= -1`` is represented
    in ``A`` and ``B``, no value will be obtained that is negative; if
    a constraint of ``x <= -1`` is represented, an error will be
    raised since no solution is possible.

    Examples
    ========

    >>> from sympy.solvers.simplex import _simplex
    >>> from sympy import Matrix

    Consider the simple minimization of ``f = x + y + 1`` under the
    constraint that ``y + 2*x >= 4``. In the nonnegative quadrant,
    this inequality describes a area above a triangle with vertices at
    (0, 4), (0, 0) and (2, 0). The minimum of ``f`` occurs at (2, 0).
    Define A, B, C, D for the standard minimization:

    >>> A, B, C, D = [Matrix(i) for i in [[[2, 1]], [4], [[1, 1]], [-1]]]

    Since `_simplex` will do a minimization for constraints given as
    ``A*x <= B``, the signs of each are negated (``-Ax >= -B``):

    >>> _simplex(-A, -B, C, D)
    (3, [2, 0], [1/2])

    The dual of minimizing ``f`` is maximizing ``F = c*y - d`` for
    ``a*y <= b`` where ``a``, ``b``, ``c``, ``d`` are derived from the
    transpose of the matrix representation of the standard minimization:

    >>> tr = lambda a, b, c, d: [i.T for i in (a, c, b, d)]
    >>> a, b, c, d = tr(A, B, C, D)

    This time ``a*x <= b`` is the expected inequality for the `_simplex`
    method, but to maximize ``F``, the sign of ``c`` and ``d`` must be
    inverted (so that minimizing the negative will give the negative of
    the maximum of ``F``):

    >>> _simplex(a, b, -c, -d)
    (-3, [1/2], [2, 0])

    The negative of ``F`` and the min of ``f`` are the same. The dual
    point `[1/2]` is the value of ``y`` that minimized ``F = c*y - d``
    under constraints a*x <= b``:

    >>> y = Matrix(['y'])
    >>> (c*y - d)[0]
    4*y + 1
    >>> [i <= j for i, j in zip(a*y,b)]
    [2*y <= 1, y <= 1]

    In this 1-dimensional dual system, the more restrictive contraint is
    the first which limits ``y`` between 0 and 1/2 and the maximum of
    ``F`` is attained at the nonzero value, hence is ``4*(1/2) + 1 = 3``.

    In this case the values for ``x`` and ``y`` were the same when the
    dual representation was solved. This is not always the case (though
    the value of the function will be the same).

    >>> l = [[1, 1], [-1, 1], [0, 1], [-1, 0]], [5, 1, 2, -1], [[1, 1]], [-1]
    >>> A, B, C, D = [Matrix(i) for i in l]
    >>> _simplex(A, B, -C, -D)
    (-6, [3, 2], [1, 0, 0, 0])
    >>> _simplex(A, B, -C, -D, dual=True)  # [5, 0] != [3, 2]
    (-6, [1, 0, 0, 0], [5, 0])

    In both cases the function has the same value:

    >>> Matrix(C)*Matrix([3, 2]) == Matrix(C)*Matrix([5, 0])
    True

    See Also
    ========
    _lp - poses min/max problem in form compatible with _simplex
    lpmin - minimization which calls _lp
    lpmax - maximimzation which calls _lp

    References
    ==========

    .. [1] Thomas S. Ferguson, LINEAR PROGRAMMING: A Concise Introduction
           web.tecnico.ulisboa.pt/mcasquilho/acad/or/ftp/FergusonUCLA_lp_args.pdf
    """

    A, B, C, D = [Matrix(i) for i in (A, B, C, D or [0])]
    if dual:
        _o, d, p = _simplex(-A.T, C.T, B.T, -D)
        return -_o, d, p

    if A and B:
        M = Matrix([[A, B], [C, D]])
    else:
        assert not A and not B  # no constraints
        M = Matrix([[C, D]])
    n = M.cols - 1
    m = M.rows - 1

    if not all(i.is_Float or i.is_Rational for i in M):
        # with literal Float and Rational we are guaranteed the
        # ability of determining whether an expression is 0 or not
        raise TypeError(filldedent("""
            Only rationals and floats are allowed in the Simplex method.
            """))

    # x variables have priority over y variables during Bland's rule
    # since False < True
    X = [(False, j) for j in range(n)]
    Y = [(True, i)  for i in range(m)]

    # Phase 1: find a feasible solution or determine none exist
    while True:
        B = M[:-1, -1]
        A = M[:-1, :-1]
        if all(B[i] >= 0 for i in range(B.rows)):
            # We have found a feasible solution
            break

        # Find k: first row with a negative rightmost entry
        for k in range(B.rows):
            if B[k] < 0:
                break  # use current value of k below
        else:
            pass  # XXX is it an error if none was found?

        # Choose pivot column, c
        piv_cols = [_ for _ in range(A.cols) if A[k, _] < 0]
        if not piv_cols:
            raise InfeasibleLPError(filldedent("""
                The constraint set is empty!"""))
        c = sorted(piv_cols, key=lambda _: X[_])[0] # Bland's rule

        # Choose pivot row, r
        piv_rows = [_ for _ in range(A.rows) if A[_, c] > 0 and B[_] > 0]
        piv_rows.append(k)
        r = _choose_pivot_row(A, B, piv_rows, c, Y)

        M = _pivot(M, r, c)
        X[c], Y[r] = Y[r], X[c]

    # Phase 2: from a feasible solution, pivot to optimal
    while True:
        B = M[:-1, -1]
        A = M[:-1, :-1]
        C = M[-1, :-1]

        # Choose a pivot column, c
        piv_cols = []
        piv_cols = [_ for _ in range(n) if C[_] < 0]
        if not piv_cols:
            break
        c = sorted(piv_cols, key=lambda _: X[_])[0] # Bland's rule

        # Choose a pivot row, r
        piv_rows = [_ for _ in range(m) if A[_, c] > 0]
        if not piv_rows:
            raise UnboundedLPError(filldedent("""
                Objective function can assume
                arbitrarily large values!"""))
        r = _choose_pivot_row(A, B, piv_rows, c, Y)

        M = _pivot(M, r, c)
        X[c], Y[r] = Y[r], X[c]

    argmax = [None]*n
    argmin_dual = [None]*m

    for i, (v, n) in enumerate(X):
        if v == False:
            argmax[n] = 0
        else:
            argmin_dual[n] = M[-1, i]

    for i, (v, n) in enumerate(Y):
        if v == True:
            argmin_dual[n] = 0
        else:
            argmax[n] = M[i, -1]

    return -M[-1, -1], argmax, argmin_dual


## routines that use _simplex or support those that do

def _abcd(M, list=False):
    """return parts of M as matrices or lists

    Examples
    ========

    >>> from sympy import Matrix
    >>> from sympy.solvers.simplex import _abcd
    >>> m = Matrix(3, 3, range(9))
    >>> L = _abcd(m, list=True); L
    ([[0, 1], [3, 4]], [2, 5], [[6, 7]], [8])
    >>> _abcd(m)
    (Matrix([
    [0, 1],
    [3, 4]]), Matrix([
    [2],
    [5]]), Matrix([[6, 7]]), Matrix([[8]]))
    >>> assert tuple(Matrix(i) for i in L) == _
    """
    def aslist(i):
        l = i.tolist()
        if len(l[0]) == 1:  # col vector
            return [i[0] for i in l]
        return l
    m = M[:-1, :-1], M[:-1, -1], M[-1, :-1], M[-1:, -1:]
    if not list:
        return m
    return tuple([aslist(i) for i in m])


def _m(a, b, c, d=None):
    """return Matrix([[a, b], [c, d]]) from matrices
    in Matrix or list form.

    Examples
    ========

    >>> from sympy import Matrix
    >>> from sympy.solvers.simplex import _abcd, _m
    >>> m = Matrix(3, 3, range(9))
    >>> L = _abcd(m, list=True); L
    ([[0, 1], [3, 4]], [2, 5], [[6, 7]], [8])
    >>> _abcd(m)
    (Matrix([
    [0, 1],
    [3, 4]]), Matrix([
    [2],
    [5]]), Matrix([[6, 7]]), Matrix([[8]]))
    >>> assert m == _m(*L) == _m(*_)
    """
    a, b, c, d = [Matrix(i) for i in (a, b, c, d or [0])]
    return Matrix([[a, b], [c, d]])


def _primal_dual(M, factor=True):
    """return primal and dual function and constraints
    assuming that ``M = Matrix([[A, b], [c, d]])`` and the
    function ``c*x - d`` is being minimized with ``Ax >= b``
    for nonnegative values of ``x``. The dual and its
    constraints will be for maximizing `b.T*y - d` subject
    to ``A.T*y <= c.T``.

    Examples
    ========

    >>> from sympy.solvers.simplex import _primal_dual, lpmin, lpmax
    >>> from sympy import Matrix

    The following matrix represents the primal task of
    minimizing x + y + 7 for y >= x + 1 and y >= -2*x + 3.
    The dual task seeks to maximize x + 3*y + 7 with
    2*y - x <= 1 and and x + y <= 1:

    >>> M = Matrix([
    ...     [-1, 1,  1],
    ...     [ 2, 1,  3],
    ...     [ 1, 1, -7]])
    >>> p, d = _primal_dual(M)

    The minimum of the primal and maximum of the dual are the same
    (though they occur at different points):

    >>> lpmin(*p)
    (28/3, {x1: 2/3, x2: 5/3})
    >>> lpmax(*d)
    (28/3, {y1: 1/3, y2: 2/3})

    If the equivalent (but canonical) inequalities are
    desired, leave `factor=True`, otherwise the unmodified
    inequalities for M will be returned.

    >>> m = Matrix([
    ... [-3, -2,  4, -2],
    ... [ 2,  0,  0, -2],
    ... [ 0,  1, -3,  0]])

    >>> _primal_dual(m, False)  # last condition is 2*x1 >= -2
    ((x2 - 3*x3,
        [-3*x1 - 2*x2 + 4*x3 >= -2, 2*x1 >= -2]),
    (-2*y1 - 2*y2,
        [-3*y1 + 2*y2 <= 0, -2*y1 <= 1, 4*y1 <= -3]))

    >>> _primal_dual(m)  # condition now x1 >= -1
    ((x2 - 3*x3,
        [-3*x1 - 2*x2 + 4*x3 >= -2, x1 >= -1]),
    (-2*y1 - 2*y2,
        [-3*y1 + 2*y2 <= 0, -2*y1 <= 1, 4*y1 <= -3]))

    If you pass the transpose of the matrix, the primal will be
    identified as the standard minimization problem and the
    dual as the standard maximization:

    >>> _primal_dual(m.T)
    ((-2*x1 - 2*x2,
        [-3*x1 + 2*x2 >= 0, -2*x1 >= 1, 4*x1 >= -3]),
    (y2 - 3*y3,
        [-3*y1 - 2*y2 + 4*y3 <= -2, y1 <= -1]))

    A matrix must have some size or else None will be returned for
    the functions:

    >>> _primal_dual(Matrix([[1, 2]]))
    ((x1 - 2, []), (-2, []))

    >>> _primal_dual(Matrix([]))
    ((None, []), (None, []))

    References
    ==========

    .. [1] David Galvin, Relations between Primal and Dual
           www3.nd.edu/~dgalvin1/30210/30210_F07/presentations/dual_opt.pdf
    """
    if not M:
        return (None, []), (None, [])
    if not hasattr(M, 'shape'):
        if len(M) not in (3, 4):
            raise ValueError('expecting Matrix or 3 or 4 lists')
        M = _m(*M)
    m, n = [i - 1 for i in M.shape]
    A, b, c, d = _abcd(M)
    d = d[0]
    _ = lambda x: numbered_symbols(x, start=1)
    x = Matrix([i for i, j in zip(_('x'), range(n))])
    yT = Matrix([i for i, j in zip(_('y'), range(m))]).T
    def ineq(L, r, op):
        rv = []
        for r in (op(i, j) for i, j in zip(L, r)):
            if r == True:
                continue
            elif r == False:
                return [False]
            if factor:
                f = factor_terms(r)
                if f.lhs.is_Mul and f.rhs % f.lhs.args[0] == 0:
                    assert len(f.lhs.args) == 2, f.lhs
                    k = f.lhs.args[0]
                    r = r.func(sign(k)*f.lhs.args[1], f.rhs//abs(k))
            rv.append(r)
        return rv
    eq = lambda x, d: x[0] - d if x else -d
    F = eq(c*x, d)
    f = eq(yT*b, d)
    return (F, ineq(A*x, b, Ge)), (f, ineq(yT*A, c, Le))


def _rel_as_nonpos(constr):
    """return a (np, d, aux) where np is a list of nonpositive
    expressions that represent the given constraints (possibly rewritten
    in terms of auxilliary variables) expressible with nonnegative
    symbols, and d is a dictionary mapping a given symbols to an
    expression with an auxilliary variable. In some cases a symbol
    will be used as part of the change of variables, e.g. x: x - z1
    instead of x: z1 - z2.

    If any constraint is False/empty, return None. Any constraint
    involving a relationship with oo will create unbounded symbols
    for all symbols present in the expression, e.g. ``x + y < oo``
    will allow ``x`` and ``y`` to take on any value. (This is the
    only type relationship which can be expressed as strict.)

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.solvers.simplex import _rel_as_nonpos
    >>> from sympy.abc import x, y
    >>> _rel_as_nonpos([x >= y])
    ([-x + y], {}, [])
    >>> _rel_as_nonpos([x > -oo])
    ([], {x: _z1 - x}, [_z1])
    >>> _rel_as_nonpos([x >= 3, x <= 5])
    ([_z1 - 2], {x: _z1 + 3}, [_z1])
    >>> _rel_as_nonpos([x <= 5])
    ([], {x: 5 - _z1}, [_z1])
    >>> _rel_as_nonpos([x >= 1])
    ([], {x: _z1 + 1}, [_z1])
    """
    r = {}  # replacements to handle change of variables
    np = []  # nonpositive expressions
    aux = []  # auxilliary symbols added
    ui = numbered_symbols('z', start=1, cls=Dummy) # auxilliary symbols
    univariate = {}  # {x: interval} for univariate constraints
    unbound = set()  # symbols designated as unbound

    # separate out univariates
    for i in constr:
        if i == True:
            continue  # ignore
        if i == False:
            return  # no solution
        i = i.canonical  # +/-oo, if present will be on the rhs
        if isinstance(i, (Le, Ge, Lt, Gt)) and i.rhs.is_infinite:
            if isinstance(i, (Lt, Le)) and i.rhs is S.NegativeInfinity:
                return False
            if isinstance(i, (Gt, Ge)) and i.rhs is S.Infinity:
                return False
            unbound.update(i.lhs.free_symbols) # x < oo, x > -oo
        elif isinstance(i, (Le, Ge)):
            i = i.lts - i.gts
            freei = i.free_symbols
            if len(freei) > 1:
                np.append(i)
            elif freei:
                x = freei.pop()
                if x in unbound:
                    continue  # will handle later
                ivl = Le(i, 0, evaluate=False).as_set()
                if x not in univariate:
                    univariate[x] = ivl
                else:
                    univariate[x] &= ivl
            elif i:
                return False
        else:
            raise TypeError('only Ge or Le is allowed, not %s' % i)

    # introduce auxilliary variables as needed for univariate
    # inequalities
    for x in univariate:
        i = univariate[x]
        if not i:
            return None  # no solution possible
        a, b = i.inf, i.sup
        if a.is_infinite and b.is_infinite:
            unbound.append(x)
        elif a.is_infinite:
            u = next(ui)
            r[x] = b - u
            aux.append(u)
        elif b.is_infinite:
            if a:
                u = next(ui)
                r[x] = a + u
                aux.append(u)
            else:
                # standard nonnegative relationship
                pass
        else:
            u = next(ui)
            aux.append(u)
            # shift so u = x - a => x = u + a
            r[x] = u + a
            # add constraint for u <= b - a
            # since when u = b-a then x = u + a = b - a + a = b:
            # the upper limit for x
            np.append(u - (b - a))

    # make change of variables for unbound variables
    for x in unbound:
        u = next(ui)
        r[x] = u - x  # reusing x
        aux.append(u)

    return np, r, aux


def _lp_matrices(objective, constraints):
    """return A, B, C, D, r, x+X, X for maximizing
    objective = Cx - D with constraints Ax <= B, introducing
    introducing auxilliary variables, X, as necessary to make
    replacements of symbols as given in r, {xi: expression with Xj},
    so all variables in x+X will take on nonnegative values.

    Every univariate condition creates a semi-infinite
    condition, e.g. a single ``x <= 3`` creates the
    interval ``[-oo, 3]`` while ``x <= 3`` and ``x >= 2``
    create an interval ``[2, 3]``. Variables not in a univariate
    expression will take on nonnegative values.
    """

    # sympify input and collect free symbols
    F = sympify(objective)
    np = [sympify(i) for i in constraints]
    syms = set.union(*[i.free_symbols for i in [F] + np], set())

    # change Eq(x, y) to x - y <= 0 and y - x <= 0
    for i in range(len(np)):
        if isinstance(np[i], Eq):
            np[i] = np[i].lhs - np[i].rhs <= 0
            np.append(-np[i].lhs <= 0)

    # convert constraints to nonpositive expressions
    _ = _rel_as_nonpos(np)
    if _ is None:
        raise InfeasibleLPError(filldedent("""
            Inconsistent/False constraint"""))
    np, r, aux = _

    # do change of variables
    F = F.xreplace(r)
    np = [i.xreplace(r) for i in np]

    # convert to matrices
    xx = list(ordered(syms)) + aux
    A, B = linear_eq_to_matrix(np, xx)
    C, D = linear_eq_to_matrix([F], xx)
    return A, B, C, D, r, xx, aux


def _lp(min_max, f, constr):
    """Return the optimization (min or max) of ``f`` with the given
    constraints. Unless constrained, variables will assume only
    nonnegative values in the solution. A constraint like ``x <= 2``
    will permit ``x`` to have values in range ``[-oo, 2]`` while
    ``x < oo`` will allow it to have any real value.

    If `min_max` is 'max' then the results corresponding to the
    maximization of ``f`` will be returned, else the minimization.
    The constraints can be given as Le, Ge or Eq expressions.

    Examples
    ========

    >>> from sympy.solvers.simplex import _lp as lp
    >>> from sympy import Eq
    >>> from sympy.abc import x, y, z
    >>> f = x + y - 2*z
    >>> c = [7*x + 4*y - 7*z <= 3, 3*x - y + 10*z <= 6]
    >>> lp(min, f, c)
    (-6/5, {x: 0, y: 0, z: 3/5})

    By passing max, the maximum value for f under the constraints
    is returned (if possible):

    >>> lp(max, f, c)
    (3/4, {x: 0, y: 3/4, z: 0})

    Constraints that are equalities will require that the solution
    also satisfy them:

    >>> lp(max, f, c + [Eq(y - 9*x, 1)])
    (5/7, {x: 0, y: 1, z: 1/7})

    All symbols are reported, even if they are not in the objective
    function:

    >>> lp(min, x, [y + x >= 3])
    (0, {x: 0, y: 3})
    """
    A, B, C, D, r, xx, aux = _lp_matrices(f, constr)

    how = str(min_max).lower()
    if 'max' in how:
        # _simplex minimizes for Ax <= B so we
        # have to change the sign of the function
        # and negate the optimal value returned
        _o, p, d = _simplex(A, B, -C, -D)
        o = -_o
    elif 'min' in how:
        o, p, d = _simplex(A, B, C, D)
    else:
        raise ValueError('expecting min or max')

    # restore original variables and remove aux from p
    p = dict(zip(xx, p))
    if r:  # p has original symbols and auxilliary symbols
        # if r has x: x - z1 use values from p to update
        r = {k: v.xreplace(p) for k, v in r.items()}
        # then use the actual value of x (= x - z1) in p
        p.update(r)
        # don't show aux
        p = {k: p[k] for k in ordered(p) if k not in aux}

    # not returning dual since there may be extra constraints
    # when a variable has finite bounds
    return o, p


def lpmin(f, constr):
    """return minimum of linear equation ``f`` under
    linear constraints expressed using Ge, Le or Eq.

    All variables will take on nonnegative values unless
    a univariate constraint indicates other ranges for
    a symbol, e.g, ``x <= 3`` will allow for negative
    values of x, while a second constraint of ``x >= 2``
    will constrain ``x`` to the interval ``[2, 3]``.

    Examples
    ========

    >>> from sympy.solvers.simplex import lpmin
    >>> from sympy import Eq
    >>> from sympy.abc import x, y
    >>> lpmin(x, [2*x - 3*y >= -1, Eq(x+ 3*y,2), x <= 2*y])
    (1/3, {x: 1/3, y: 5/9})

    Negative values for variables are permitted unless explicitly
    exluding them.

    >>> lpmin(x, [x >= 0, x <= 3])
    (0, {x: 0})

     Without indicating that ``x`` is nonnegative, there
    is no minimum for this objective:

    >>> lpmin(x, [x <= 3])
    Traceback (most recent call last):
    ...
    sympy.solvers.simplex.UnboundedLPError:
    Objective function can assume arbitrarily large values!

    See Also
    ========
    linprog, lpmax
    """
    return _lp(min, f, constr)


def lpmax(f, constr):
    """return maximum of linear equation ``f`` under
    linear constraints expressed using Ge, Le or Eq.

    All variables will take on nonnegative values unless
    a univariate constraint indicates other ranges for
    a symbol, e.g, ``x <= 3`` will allow for negative
    values of x, while a second constraint of ``x >= 2``
    will constrain ``x`` to the interval ``[2, 3]``.

    Examples
    ========

    >>> from sympy.solvers.simplex import lpmax
    >>> from sympy import Eq
    >>> from sympy.abc import x, y
    >>> lpmax(x, [2*x - 3*y >= -1, Eq(x+ 3*y,2), x <= 2*y])
    (4/5, {x: 4/5, y: 2/5})

    Univariate conditions will allow negative values unless a
    constraint like ``x >= 0`` were included:

    >>> lpmax(x, [x <= -1])
    (-1, {x: -1})
    >>> lpmax(x, [x <= -1, x >= 0])
    Traceback (most recent call last):
    ...
    sympy.solvers.simplex.InfeasibleLPError: inconsistent/False constraint

    See Also
    ========
    linprog, lpmin
    """
    return _lp(max, f, constr)


def linprog(cd, A=None, b=None, A_eq=None, b_eq=None, bounds=None):
    """Return the minimization of ``c*x - d`` with the given
    constraints ``A*x <= b`` and ``A_eq*x = b_eq``. Unless bounds
    are given, variables will have nonnegative values in the solution.

    If ``cd`` is a single matrix then ``d`` will be zero; otherwise
    it should be tuple with two matrices: ``c`` and ``d``.

    By default, all variables will be nonnegative. If ``bounds``
    is given as a single tuple, ``(lo, hi)``, then all variables
    will be constrained to be between ``lo`` and ``hi``. Use
    None for a ``lo`` or ``hi`` if it is unconstrained in the
    negative or positive direction, respectively, e.g.
    ``(None, 0)`` indicates nonpositive values. To set
    individual ranges, pass a list with length equal to the
    number of columns in ``A``, each element being a tuple; if
    only a few variables take on non-default values they can be
    passed as a dictionary with keys giving the corresponding
    column to which the variable is assigned, e.g. ``bounds={2:
    (1, 4)}`` would limit the 3rd variable to have a value in
    range ``[1, 4]``.

    Examples
    ========

    >>> from sympy.solvers.simplex import linprog
    >>> from sympy import symbols, Eq, linear_eq_to_matrix as M, Matrix
    >>> x = x1, x2, x3, x4 = symbols('x1:5')
    >>> X = Matrix(x)
    >>> c, d = cd = M(5*x2 + x3 + 4*x4 - x1, x)
    >>> a, b = M([5*x2 + 2*x3 + 5*x4 - (x1 + 5)], x)
    >>> aeq, beq = M([Eq(3*x2 + x4, 2), Eq(-x1 + x3 + 2*x4, 1)], x)
    >>> constr = [i <= j for i,j in zip(a*X, b)]
    >>> constr += [Eq(i, j) for i,j in zip(aeq*X, beq)]
    >>> linprog(cd, a, b, aeq, beq)
    (9/2, {x1: 0, x2: 1/2, x3: 0, x4: 1/2})
    >>> assert all(i.subs(_[1]) for i in constr)

    See Also
    ========
    lpmin, lpmax
    """

    # convert scipy-like matrices to expressions

    ## the inequalities
    if A is None:
        if b is not None:
            raise ValueError('A and b must both be given')
        A, b = Matrix(0, 0, []), Matrix(0, 0, [])
    else:
        A, b = [Matrix(i) for i in (A, b)]

    ## the equalities
    if A_eq is None:
        if not b_eq is None:
            raise ValueError('A_eq and b_eq must both be given')
    else:
        A_eq, b_eq = [Matrix(i) for i in (A_eq, b_eq)]
        # if x == y then x <= y and x >= y (-x <= -y)
        A = A.col_join(A_eq)
        A = A.col_join(-A_eq)
        b = b.col_join(b_eq)
        b = b.col_join(-b_eq)

    ## the objective
    if isinstance(cd, MatrixBase):
        C, D = cd, Matrix([0])
    elif len(cd) == 2 and all(isinstance(i, (MatrixBase, list)) for i in cd):
        C, D = [Matrix(i) for i in cd]
    elif isinstance(cd, list) and isinstance(cd[0], list) and len(cd[0]) == A.cols:
        C, D = Matrix(cd), Matrix([0])
    else:
        raise ValueError('expecting one or two matrices/lists for cd')

    ## create variables to re-create the symbolic constraints
    ## and objective if necessary, else to return solution
    xit = numbered_symbols('x', start=1)
    x = Matrix([i[1] for i in zip(range(A.cols), xit)])

    ## test for quick exit when all bounds indicate nonnegatives
    if bounds in (None, (0, None)) or set(bounds) == {(0, None)}:
        o, p, d = _simplex(A, b, C, D)
        return o, dict(zip(x, p))

    ## recreate the objective and constraints so aux variables
    ## can be created to handle the bounds
    ### TODO: handle aux variable creation here?
    constr = [i <= j for i, j in zip(A*x, b)]
    f = (C*x - D)[0]

    ## the bounds are interpreted
    if bounds is None:
        bounds = (0, None)
    if type(bounds) is tuple and len(bounds) == 2:
        bounds = [bounds]*A.cols
    elif len(bounds) == A.cols and all(type(i) is tuple and len(i) == 2 for i in bounds):
        pass # individual bounds
    elif type(bounds) is dict and all(type(i) is tuple and len(i) == 2 for i in bounds.values()):
        # sparse bounds
        db = bounds
        bounds = [(0, None)]*A.cols
        while db:
            i, b = db.popitem()
            bounds[i] = b
    else:
        raise ValueError('unexpected bounds %s' % bounds)
    for i, (lo,hi) in enumerate(bounds):
        if lo == 0 and hi is None:
            constr.append(x[i] >= 0)
        elif lo is not None:
            if hi is None:
                constr.append(x[i] >= lo)
            else:
                constr.append(x[i] <= hi)
                constr.append(x[i] >= lo)
        elif lo is None and hi is not None:
            constr.append(x[i] <= hi)
        elif lo is None and hi is None:
            constr.append(x[i] <= S.Infinity)

    ## solve the symbolic system with aux variable handling

    return _lp(min, f, constr)
