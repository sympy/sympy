"""
This module contain solvers for all kinds of equations:

    - algebraic or transcendental, use solve()

    - recurrence, use rsolve()

    - differential, use dsolve()

    - nonlinear (numerically), use nsolve()
      (you will need a good starting point)

"""

from sympy.core.compatibility import iterable, is_sequence
from sympy.utilities.exceptions import SymPyDeprecationWarning
from sympy.core.sympify import sympify
from sympy.core import C, S, Add, Symbol, Wild, Equality, Dummy, Basic, Expr
from sympy.core.function import (expand_mul, expand_multinomial, expand_log,
                          Derivative, AppliedUndef, UndefinedFunction, nfloat,
                          count_ops)
from sympy.core.numbers import ilcm, Float
from sympy.core.relational import Relational
from sympy.logic.boolalg import And, Or
from sympy.core.basic import preorder_traversal

from sympy.functions import (log, exp, LambertW, cos, sin, tan, cot, cosh,
                             sinh, tanh, coth, acos, asin, atan, acot, acosh,
                             asinh, atanh, acoth, Abs)
from sympy.functions.elementary.miscellaneous import real_root
from sympy.simplify import (simplify, collect, powsimp, posify, powdenest,
                            nsimplify)
from sympy.simplify.sqrtdenest import sqrt_depth, _mexpand
from sympy.matrices import Matrix, zeros
from sympy.polys import roots, cancel, Poly, together, factor
from sympy.functions.elementary.piecewise import piecewise_fold, Piecewise

from sympy.utilities.iterables import sift
from sympy.utilities.lambdify import lambdify
from sympy.utilities.misc import default_sort_key, filldedent
from sympy.mpmath import findroot

from sympy.solvers.polysys import solve_poly_system
from sympy.solvers.inequalities import reduce_inequalities

from sympy.core.compatibility import reduce

from sympy.assumptions import Q, ask

from types import GeneratorType
from collections import defaultdict


def _ispow(e):
    """Return True if e is a Pow or is exp."""
    return e.is_Pow or e.func is exp


def denoms(eq, symbols=None):
    """Return (recursively) set of all denominators that appear in eq
    that contain any symbol in iterable ``symbols``; if ``symbols`` is
    None (default) then all denominators with symbols will be returned.

    Examples
    ========

    >>> from sympy.solvers.solvers import denoms
    >>> from sympy.abc import x, y, z

    >>> denoms(x/y)
    set([y])

    >>> denoms(x/(y*z))
    set([y, z])

    >>> denoms(3/x + y/z)
    set([x, z])
    """

    symbols = symbols or eq.free_symbols
    dens = set()
    if not symbols or not eq.has(*symbols):
        return dens
    pt = preorder_traversal(eq)
    for e in pt:
        if _ispow(e):
            n, d = e.as_numer_denom()
            if d in dens:
                pt.skip()
            elif d.has(*symbols):
                dens.add(d.as_base_exp()[0])
    return dens


def checksol(f, symbol, sol=None, **flags):
    """Checks whether sol is a solution of equation f == 0.

    Input can be either a single symbol and corresponding value
    or a dictionary of symbols and values. ``f`` can be a single
    equation or an iterable of equations. A solution must satisfy
    all equations in ``f`` to be considered valid; if a solution
    does not satisfy any equation, False is returned; if one or
    more checks are inconclusive (and none are False) then None
    is returned.

    Examples
    ========

       >>> from sympy import symbols
       >>> from sympy.solvers import checksol
       >>> x, y = symbols('x,y')
       >>> checksol(x**4-1, x, 1)
       True
       >>> checksol(x**4-1, x, 0)
       False
       >>> checksol(x**2 + y**2 - 5**2, {x:3, y: 4})
       True

       None is returned if checksol() could not conclude.

       flags:
           'numerical=True (default)'
               do a fast numerical check if ``f`` has only one symbol.
           'minimal=True (default is False)'
               a very fast, minimal testing.
           'warn=True (default is False)'
               print a warning if checksol() could not conclude.
           'simplify=True (default)'
               simplify solution before substituting into function and
               simplify the function before trying specific simplifications
           'force=True (default is False)'
               make positive all symbols without assumptions regarding sign.

    """

    if sol is not None:
        sol = {symbol: sol}
    elif isinstance(symbol, dict):
        sol = symbol
    else:
        msg = 'Expecting sym, val or {sym: val}, None but got %s, %s'
        raise ValueError(msg % (symbol, sol))

    if iterable(f):
        if not f:
            raise ValueError('no functions to check')
        rv = True
        for fi in f:
            check = checksol(fi, sol, **flags)
            if check:
                continue
            if check is False:
                return False
            rv = None # don't return, wait to see if there's a False
        return rv

    if isinstance(f, Poly):
        f = f.as_expr()
    elif isinstance(f, Equality):
        f = f.lhs - f.rhs

    if not f:
        return True

    if not f.has(*sol.keys()):
        # if f(y) == 0, x=3 does not set f(y) to zero...nor does it not
        return None

    illegal = set([S.NaN,
               S.ComplexInfinity,
               S.Infinity,
               S.NegativeInfinity])
    if any(sympify(v).atoms() & illegal for k, v in sol.iteritems()):
        return False

    was = f
    attempt = -1
    numerical = flags.get('numerical', True)
    while 1:
        attempt += 1
        if attempt == 0:
            val = f.subs(sol)
            if val.atoms() & illegal:
                return False
        elif attempt == 1:
            if val.free_symbols:
                if not val.is_constant(*sol.keys()):
                    return False
                # there are free symbols -- simple expansion might work
                _, val = val.as_content_primitive()
                val = expand_mul(expand_multinomial(val))
        elif attempt == 2:
            if flags.get('minimal', False):
                return
            if flags.get('simplify', True):
                for k in sol:
                    sol[k] = simplify(sol[k])
            # start over without the failed expanded form, possibly
            # with a simplified solution
            val = f.subs(sol)
            if flags.get('force', True):
                val, reps = posify(val)
                # expansion may work now, so try again and check
                exval = expand_mul(expand_multinomial(val))
                if exval.is_number or not exval.free_symbols:
                    # we can decide now
                    val = exval
        elif attempt == 3:
            val = powsimp(val)
        elif attempt == 4:
            val = cancel(val)
        elif attempt == 5:
            val = val.expand()
        elif attempt == 6:
            val = together(val)
        elif attempt == 7:
            val = powsimp(val)
        else:
            # if there are no radicals and no functions then this can't be
            # zero anymore -- can it?
            pot = preorder_traversal(expand_mul(val))
            seen = set()
            saw_pow_func = False
            for p in pot:
                if p in seen:
                    continue
                seen.add(p)
                if p.is_Pow and not p.exp.is_Integer:
                    saw_pow_func = True
                elif p.is_Function:
                    saw_pow_func = True
                elif isinstance(p, UndefinedFunction):
                    saw_pow_func = True
                if saw_pow_func:
                    break
            if saw_pow_func is False:
                return False
            if flags.get('force', True):
                # don't do a zero check with the positive assumptions in place
                val = val.subs(reps)
            nz = val.is_nonzero
            if nz is not None:
                # issue 2574: nz may be True even when False
                # so these are just hacks to keep a false positive
                # from being returned

                # HACK 1: LambertW (issue 2574)
                if val.is_number and val.has(LambertW):
                    # don't eval this to verify solution since if we got here,
                    # numerical must be False
                    return None

                # add other HACKs here if necessary, otherwise we assume
                # the nz value is correct
                return not nz
            break

        if val == was:
            continue
        elif val.is_Rational:
            return val == 0
        if numerical and not val.free_symbols:
            return abs(val.n(chop=True)) < 1e-9
        was = val

    if flags.get('warn', False):
        print("\n\tWarning: could not verify solution %s." % sol)
    # returns None if it can't conclude
    # TODO: improve solution testing


def check_assumptions(expr, **assumptions):
    """Checks whether expression `expr` satisfies all assumptions.

    `assumptions` is a dict of assumptions: {'assumption': True|False, ...}.

    Examples
    ========

       >>> from sympy import Symbol, pi, I, exp
       >>> from sympy.solvers.solvers import check_assumptions

       >>> check_assumptions(-5, integer=True)
       True
       >>> check_assumptions(pi, real=True, integer=False)
       True
       >>> check_assumptions(pi, real=True, negative=True)
       False
       >>> check_assumptions(exp(I*pi/7), real=False)
       True

       >>> x = Symbol('x', real=True, positive=True)
       >>> check_assumptions(2*x + 1, real=True, positive=True)
       True
       >>> check_assumptions(-2*x - 5, real=True, positive=True)
       False

       `None` is returned if check_assumptions() could not conclude.

       >>> check_assumptions(2*x - 1, real=True, positive=True)
       >>> z = Symbol('z')
       >>> check_assumptions(z, real=True)
    """
    expr = sympify(expr)

    result = True
    for key, expected in assumptions.iteritems():
        if expected is None:
            continue
        if expected in [0, 1]:
            expected = bool(expected)
        if not isinstance(expected, bool):
            raise ValueError(_filldendent('''
                A boolean is expected for %s but got %s.''' % (key, expected)))
        if hasattr(Q, key):
            test = ask(getattr(Q, key)(expr))
            if test is expected:
                continue
            elif test is not None:
                return False
        # ask() can't conclude. Try using old assumption system.
        # XXX: remove once transition to new assumption system is finished.
        test = getattr(expr, 'is_' + key, None)
        if test is expected:
            continue
        elif test is not None:
            return False
        result = None # Can't conclude, unless an other test fails.
    return result


def solve(f, *symbols, **flags):
    """
    Algebraically solves equations and systems of equations.

    Currently supported are:
        - univariate polynomial,
        - transcendental
        - piecewise combinations of the above
        - systems of linear and polynomial equations
        - sytems containing relational expressions.

    Input is formed as:

    * f
        - a single Expr or Poly that must be zero,
        - an Equality
        - a Relational expression or boolean
        - iterable of one or more of the above

    * symbols (Symbol, Function or Derivative) specified as
        - none given (all free symbols will be used)
        - single symbol
        - denested list of symbols
          e.g. solve(f, x, y)
        - ordered iterable of symbols
          e.g. solve(f, [x, y])

    * flags
        'dict'=True (default is False)
            return list (perhaps empty) of solution mappings
        'set'=True (default is False)
            return list of symbols and set of tuple(s) of solution(s)
        'exclude=[] (default)'
            don't try to solve for any of the free symbols in exclude;
            if expressions are given, the free symbols in them will
            be extracted automatically.
        'check=True (default)'
            If False, don't do any testing of solutions. This can be
            useful if one wants to include solutions that make any
            denominator zero.
        'numerical=True (default)'
            do a fast numerical check if ``f`` has only one symbol.
        'minimal=True (default is False)'
            a very fast, minimal testing.
        'warning=True (default is False)'
            print a warning if checksol() could not conclude.
        'simplify=True (default)'
            simplify all but cubic and quartic solutions before
            returning them and (if check is not False) use the
            general simplify function on the solutions and the
            expression obtained when they are substituted into the
            function which should be zero
        'force=True (default is False)'
            make positive all symbols without assumptions regarding sign.
        'rational=True (default)'
            recast Floats as Rational; if this option is not used, the
            system containing floats may fail to solve because of issues
            with polys. If rational=None, Floats will be recast as
            rationals but the answer will be recast as Floats. If the
            flag is False then nothing will be done to the Floats.
        'manual=True (default is False)'
            do not use the polys/matrix method to solve a system of
            equations, solve them one at a time as you might "manually".
        'implicit=True (default is False)'
            allows solve to return a solution for a pattern in terms of
            other functions that contain that pattern; this is only
            needed if the pattern is inside of some invertible function
            like cos, exp, ....

    Examples
    ========

    The output varies according to the input and can be seen by example::

        >>> from sympy import solve, Poly, Eq, Function, exp
        >>> from sympy.abc import x, y, z, a, b
        >>> f = Function('f')

    * boolean or univariate Relational

        >>> solve(x < 3)
        And(im(x) == 0, re(x) < 3)

    * to always get a list of solution mappings, use flag dict=True

        >>> solve(x - 3, dict=True)
        [{x: 3}]
        >>> solve([x - 3, y - 1], dict=True)
        [{x: 3, y: 1}]

    * to get a list of symbols and set of solution(s) use flag set=True

        >>> solve([x**2 - 3, y - 1], set=True)
        ([x, y], set([(-sqrt(3), 1), (sqrt(3), 1)]))

    * single expression and single symbol that is in the expression

        >>> solve(x - y, x)
        [y]
        >>> solve(x - 3, x)
        [3]
        >>> solve(Eq(x, 3), x)
        [3]
        >>> solve(Poly(x - 3), x)
        [3]
        >>> set(solve(x**2 - y**2, x))
        set([-y, y])
        >>> set(solve(x**4 - 1, x))
        set([-1, 1, -I, I])

    * single expression with no symbol that is in the expression

        >>> solve(3, x)
        []
        >>> solve(x - 3, y)
        []

    * single expression with no symbol given

          In this case, all free symbols will be selected as potential
          symbols to solve for. If the equation is univariate then a list
          of solutionsis returned; otherwise -- as is the case when symbols are
          given as an iterable of length > 1 -- a list of mappings will be returned.

            >>> solve(x - 3)
            [3]
            >>> solve(x**2 - y**2) # doctest: +SKIP
            [{x: -y}, {x: y}]
            >>> solve(z**2*x**2 - z**2*y**2) # doctest: +SKIP
            [{x: -y}, {x: y}]
            >>> solve(z**2*x - z**2*y**2)
            [{x: y**2}]

    * when a Function or Derivative is given as a symbol, it is
      isolated algebraically and an implicit solution may be obtained;
      to obtain the solution for a function within a derivative, use
      dsolve.

          >>> solve(f(x) - x, f(x))
          [x]
          >>> solve(f(x).diff(x) - f(x) - x, f(x).diff(x))
          [x + f(x)]
          >>> solve(f(x).diff(x) - f(x) - x, f(x))
          [-x + Derivative(f(x), x)]
          >>> set(solve(x + exp(x)**2, exp(x)))
          set([-sqrt(-x), sqrt(-x)])

        * To solve for a *symbol* implicitly, use 'implicit=True':

            >>> solve(x + exp(x), x)
            [-LambertW(1)]
            >>> solve(x + exp(x), x, implicit=True)
            [-exp(x)]

    * single expression and more than 1 symbol

        * when there is a linear solution

            >>> solve(x - y**2, x, y)
            [{x: y**2}]
            >>> solve(x**2 - y, x, y)
            [{y: x**2}]

        * when undetermined coefficients are identified

            * that are linear

                >>> solve((a + b)*x - b + 2, a, b)
                {a: -2, b: 2}

            * that are nonlinear

                >>> set(solve((a + b)*x - b**2 + 2, a, b))
                set([(-sqrt(2), sqrt(2)), (sqrt(2), -sqrt(2))])

        * if there is no linear solution then the first successful
          attempt for a nonlinear solution will be returned

            >>> solve(x**2 - y**2, x, y) # doctest: +SKIP
            [{x: -y}, {x: y}]
            >>> solve(x**2 - y**2/exp(x), x, y)
            [{x: 2*LambertW(y/2)}]
            >>> solve(x**2 - y**2/exp(x), y, x) # doctest: +SKIP
            [{y: -x*exp(x/2)}, {y: x*exp(x/2)}]

    * iterable of one or more of the above

        * involving relationals or bools

            >>> solve([x < 3, x - 2])
            And(re(x) == 2, im(x) == 0)
            >>> solve([x > 3, x - 2])
            False

        * when the system is linear

            * with a solution

                >>> solve([x - 3], x)
                {x: 3}
                >>> solve((x + 5*y - 2, -3*x + 6*y - 15), x, y)
                {x: -3, y: 1}
                >>> solve((x + 5*y - 2, -3*x + 6*y - 15), x, y, z)
                {x: -3, y: 1}
                >>> solve((x + 5*y - 2, -3*x + 6*y - z), z, x, y)
                {x: -5*y + 2, z: 21*y - 6}

            * without a solution

                >>> solve([x + 3, x - 3])
                []

        * when the system is not linear

            >>> set(solve([x**2 + y -2, y**2 - 4], x, y))
            set([(-2, -2), (0, 2), (2, -2)])

        * if no symbols are given, all free symbols will be selected and a list
          of mappings returned

            >>> solve([x - 2, x**2 + y])
            [{x: 2, y: -4}]
            >>> solve([x - 2, x**2 + f(x)], set([f(x), x]))
            [{x: 2, f(x): -4}]

        * if any equation doesn't depend on the symbol(s) given it will be
          eliminated from the equation set and an answer may be given
          implicitly in terms of variables that were not of interest

            >>> solve([x - y, y - 3], x)
            {x: y}

    Notes
    =====

    assumptions aren't checked when `solve()` input involves
    relationals or bools.

    When the solutions are checked, those that make any denominator zero
    are automatically excluded. If you do not want to exclude such solutions
    then use the check=False option:

        >>> from sympy import sin, limit
        >>> solve(sin(x)/x)
        []

    If check=False then a solution to the numerator being zero is found: x = 0.
    In this case, this is a spurious solution since sin(x)/x has the well known
    limit (without dicontinuity) of 1 at x = 0:

        >>> solve(sin(x)/x, check=False)
        [0]

    In the following case, however, the limit exists and is equal to the the
    value of x = 0 that is excluded when check=True:

        >>> eq = x**2*(1/x - z**2/x)
        >>> solve(eq, x)
        []
        >>> solve(eq, x, check=False)
        [0]
        >>> limit(eq, x, 0, '-')
        0
        >>> limit(eq, x, 0, '+')
        0


    See Also
    ========

        - rsolve() for solving recurrence relationships
        - dsolve() for solving differential equations

    """
    # make f and symbols into lists of sympified quantities
    # keeping track of how f was passed since if it is a list
    # a dictionary of results will be returned.
    ###########################################################################

    def _sympified_list(w):
        return map(sympify, w if iterable(w) else [w])
    bare_f = not iterable(f)
    ordered_symbols = (symbols and
                       symbols[0] and
                       (isinstance(symbols[0], Symbol) or
                        is_sequence(symbols[0],
                        include=GeneratorType)
                       )
                      )
    f, symbols = (_sympified_list(w) for w in [f, symbols])

    implicit = flags.get('implicit', False)

    # preprocess equation(s)
    ###########################################################################
    for i, fi in enumerate(f):
        if isinstance(fi, Equality):
            f[i] = fi.lhs - fi.rhs
        elif isinstance(fi, Poly):
            f[i] = fi.as_expr()
        elif isinstance(fi, bool) or fi.is_Relational:
            return reduce_inequalities(f, assume=flags.get('assume'),
                                          symbols=symbols)
        # Any embedded piecewise functions need to be brought out to the
        # top level so that the appropriate strategy gets selected.
        f[i] = piecewise_fold(f[i])

    # preprocess symbol(s)
    ###########################################################################
    if not symbols:
        # get symbols from equations or supply dummy symbols so solve(3)
        # behaves like solve(3, x).
        symbols = reduce(set.union, [fi.free_symbols or set([Dummy()])
                                     for fi in f], set())
        ordered_symbols = False
    elif len(symbols) == 1 and iterable(symbols[0]):
        symbols = symbols[0]

    # remove symbols the user is not interested in
    exclude = flags.pop('exclude', set())
    if exclude:
        if isinstance(exclude, Expr):
            exclude = [exclude]
        exclude = reduce(set.union, [e.free_symbols for e in sympify(exclude)])
    symbols = [s for s in symbols if s not in exclude]

    if not ordered_symbols:
        # we do this to make the results returned canonical in case f
        # contains a system of nonlinear equations; all other cases should
        # be unambiguous
        symbols = sorted(symbols, key=default_sort_key)

    # we can solve for Function and Derivative instances by replacing them
    # with Dummy symbols or functions
    symbols_new = []
    symbol_swapped = False
    funcs = []
    for i, s in enumerate(symbols):
        if s.is_Symbol:
            s_new = s
        elif s.is_Function:
            symbol_swapped = True
            s_new = Dummy('F%d' % i)
            funcs.append(s)
        elif s.is_Derivative:
            symbol_swapped = True
            s_new = Dummy('D%d' % i)
        elif s.is_Pow:
            symbol_swapped = True
            s_new = Dummy('P%d' % i)
        else:
            msg = 'expected Symbol, Function, Power or Derivative but got %s'
            raise TypeError(msg % type(s))
        symbols_new.append(s_new)

    if symbol_swapped:
        swap_sym = zip(symbols, symbols_new)
        f = [fi.subs(swap_sym) for fi in f]
        symbols = symbols_new
        swap_sym = dict([(v, k) for k, v in swap_sym])
    else:
        swap_sym = {}

    # this is needed in the next two events
    symset = set(symbols)

    # get rid of equations that have no symbols of interest; we don't
    # try to solve them because the user didn't ask and they might be
    # hard to solve; this means that solutions may be give in terms
    # of the eliminated equations e.g. solve((x-y, y-3), x) -> {x: y}
    newf = []
    for fi in f:
        # let the solver handle equations that..
        # - have no symbols but are expressions
        # - have symbols of interest
        # - have no symbols of interest but are constant
        # but when an expression is not constant and has no symbols of
        # interest, it can't change what we obtain for a solution from
        # the remaining equations so we don't include it; and if it's
        # zero it can be removed and if it's not zero, there is no
        # solution for the equation set as a whole
        #
        # The reason for doing this filtering is to allow an answer
        # to be obtained to queries like solve((x - y, y), x); without
        # this mod the return value is []
        ok = False
        if fi.has(*symset):
            ok = True
        else:
            free = fi.free_symbols
            if not free:
                if fi.is_Number:
                    if fi.is_zero:
                        continue
                    return []
                ok = True
            else:
              if fi.is_constant():
                ok = True
        if ok:
            newf.append(fi)
    if not newf:
        return []
    f = newf
    del newf

    # mask off any Object that we aren't going to invert: Derivative,
    # Integral, etc... so that solving for anything that they contain will
    # give an implicit solution
    seen = set()
    non_inverts = set()
    for fi in f:
        pot = preorder_traversal(fi)
        for p in pot:
            if isinstance(p, bool) or isinstance(p, Piecewise):
                pass
            elif (isinstance(p, bool) or
                not p.args or
                p in symset or
                p.is_Add or p.is_Mul or
                p.is_Pow and not implicit or
                p.is_Function and not isinstance(p, AppliedUndef) and
                not implicit):
                continue
            elif not p in seen:
                seen.add(p)
                if p.free_symbols & symset:
                    non_inverts.add(p)
                else:
                    continue
            pot.skip()
    del seen
    non_inverts = dict(zip(non_inverts, [Dummy() for d in non_inverts]))
    f = [fi.subs(non_inverts) for fi in f]
    non_inverts = [(v, k.subs(swap_sym)) for k, v in non_inverts.iteritems()]

    # rationalize Floats
    floats = False
    if flags.get('rational', True) is not False:
        for i, fi in enumerate(f):
            if fi.has(Float):
                floats = True
                f[i] = nsimplify(fi, rational=True)

    #
    # try to get a solution
    ###########################################################################
    if bare_f:
        solution = _solve(f[0], *symbols, **flags)
    else:
        solution = _solve_system(f, symbols, **flags)

    #
    # postprocessing
    ###########################################################################
    # Restore masked off derivatives
    if non_inverts:

        def _do_dict(solution):
            return dict([(k, v.subs(non_inverts)) for k, v in
                          solution.iteritems()])
        for i in range(1):
            if type(solution) is dict:
                solution = _do_dict(solution)
                break
            elif solution and type(solution) is list:
                if type(solution[0]) is dict:
                    solution = [_do_dict(s) for s in solution]
                    break
                elif type(solution[0]) is tuple:
                    solution = [tuple([v.subs(non_inverts) for v in s]) for s
                                in solution]
                    break
                else:
                    solution = [v.subs(non_inverts) for v in solution]
                    break
            elif not solution:
                break
        else:
            raise NotImplementedError(filldedent('''
                            no handling of %s was implemented''' % solution))

    # Restore original Functions and Derivatives if a dictionary is returned.
    # This is not necessary for
    #   - the single univariate equation case
    #     since the symbol will have been removed from the solution;
    #   - the nonlinear poly_system since that only supports zero-dimensional
    #     systems and those results come back as a list
    #
    # ** unless there were Derivatives with the symbols, but those were handled
    #    above.
    if symbol_swapped:
        if type(solution) is dict:
            solution = dict([(swap_sym[k], v.subs(swap_sym))
                              for k, v in solution.iteritems()])
        elif solution and type(solution) is list and type(solution[0]) is dict:
            for i, sol in enumerate(solution):
                solution[i] = dict([(swap_sym[k], v.subs(swap_sym))
                              for k, v in sol.iteritems()])
    # undo the dictionary solutions returned when the system was only partially
    # solved with poly-system if all symbols are present
    if (
            solution and
            ordered_symbols and
            type(solution) is not dict and
            type(solution[0]) is dict and
            all(s in solution[0] for s in symbols)
            ):
        solution = [tuple([r[s].subs(r) for s in symbols]) for r in solution]

    # Get assumptions about symbols, to filter solutions.
    # Note that if assumptions about a solution can't be verified, it is still
    # returned.
    check = flags.get('check', True)

    # restore floats
    if floats and solution and flags.get('rational', None) is None:
        solution = nfloat(solution, exponent=False)

    if check and solution:

        warning = flags.get('warn', False)
        got_None = [] # solutions for which one or more symbols gave None
        no_False = [] # solutions for which no symbols gave False
        if type(solution) is list:
            if type(solution[0]) is tuple:
                for sol in solution:
                    for symb, val in zip(symbols, sol):
                        test = check_assumptions(val, **symb.assumptions0)
                        if test is False:
                            break
                        if test is None:
                            got_None.append(sol)
                    else:
                        no_False.append(sol)
            elif type(solution[0]) is dict:
                for sol in solution:
                    a_None = False
                    for symb, val in sol.iteritems():
                        test = check_assumptions(val, **symb.assumptions0)
                        if test:
                            continue
                        if test is False:
                            break
                        a_None = True
                    else:
                        no_False.append(sol)
                        if a_None:
                            got_None.append(sol)
            else: # list of expressions
                for sol in solution:
                    test = check_assumptions(sol, **symbols[0].assumptions0)
                    if test is False:
                        continue
                    no_False.append(sol)
                    if test is None:
                        got_None.append(sol)

        elif type(solution) is dict:
            a_None = False
            for symb, val in solution.iteritems():
                test = check_assumptions(val, **symb.assumptions0)
                if test:
                    continue
                if test is False:
                    no_False = None
                    break
                a_None = True
            else:
                no_False = solution
                if a_None:
                    got_None.append(solution)

        elif isinstance(solution, (Relational, And, Or)):
            assert len(symbols) == 1
            if warning and symbols[0].assumptions0:
                print(filldedent("""
                    \tWarning: assumptions about variable '%s' are
                    not handled currently.""" %symbols[0]))
            # TODO: check also variable assumptions for inequalities

        else:
            raise TypeError('Unrecognized solution') # improve the checker

        solution = no_False
        if warning and got_None:
            print(filldedent("""
                \tWarning: assumptions concerning following solution(s)
                can't be checked:""" + '\n\t' +
                ', '.join(str(s) for s in got_None)))

    #
    # done
    ###########################################################################

    as_dict = flags.get('dict', False)
    as_set = flags.get('set', False)

    if not as_set and isinstance(solution, list):
        # Make sure that a list of solutions is ordered in a canonical way.
        solution.sort(key=default_sort_key)

    if not as_dict and not as_set:
        return solution or []

    # make return a list of mappings or []
    if not solution:
        solution = []
    else:
        if isinstance(solution, dict):
            solution = [solution]
        elif iterable(solution[0]):
            solution = [dict(zip(symbols, s)) for s in solution]
        elif isinstance(solution[0], dict):
            pass
        else:
            assert len(symbols) == 1
            solution = [{symbols[0]: s} for s in solution]
    if as_dict:
        return solution
    assert as_set
    if not solution:
        return [], set()
    k = sorted(solution[0].keys(), key=lambda i: i.sort_key())
    return k, set([tuple([s[ki] for ki in k]) for s in solution])


def _solve(f, *symbols, **flags):
    """Return a checked solution for f in terms of one or more of the
    symbols."""

    if len(symbols) != 1:
        soln = None
        free = f.free_symbols
        ex = free - set(symbols)
        if len(ex) != 1:
            ind, dep = f.as_independent(*symbols)
            ex = ind.free_symbols & dep.free_symbols
        if len(ex) == 1:
            ex = ex.pop()
            try:
                # may come back as dict or list (if non-linear)
                soln = solve_undetermined_coeffs(f, symbols, ex)
            except NotImplementedError:
                pass
        if soln:
            return soln
        # find first successful solution
        failed = []
        for s in symbols:
            n, d = solve_linear(f, symbols=[s])
            if n.is_Symbol:
                # no need to check but we should simplify if desired
                if flags.get('simplify', True):
                    d = simplify(d)
                return [{n: d}]
            elif n and d: # otherwise there was no solution for s
                failed.append(s)
        if not failed:
            return []
        for s in failed:
            try:
                soln = _solve(f, s, **flags)
                return [{s: sol} for sol in soln]
            except NotImplementedError:
                continue
        else:
            msg = "No algorithms are implemented to solve equation %s"
            raise NotImplementedError(msg % f)
    symbol = symbols[0]

    check = flags.get('check', True)
    # build up solutions if f is a Mul
    if f.is_Mul:
        result = set()
        dens = denoms(f, symbols)
        for m in f.args:
            soln = _solve(m, symbol, **flags)
            result.update(set(soln))
        result = list(result)
        if check:
            result = [s for s in result if
                all(not checksol(den, {symbol: s}, **flags) for den in dens)]
        # set flags for quick exit at end
        check = False
        flags['simplify'] = False

    elif f.is_Piecewise:
        result = set()
        for expr, cond in f.args:
            candidates = _solve(expr, *symbols)
            if cond is True:
                # Only include solutions that do not match the condition
                # of any of the other pieces.
                for candidate in candidates:
                    matches_other_piece = False
                    for other_expr, other_cond in f.args:
                        if other_cond is True:
                            continue
                        if bool(other_cond.subs(symbol, candidate)):
                            matches_other_piece = True
                            break
                    if not matches_other_piece:
                        result.add(candidate)
            else:
                for candidate in candidates:
                    if bool(cond.subs(symbol, candidate)):
                        result.add(candidate)
        check = False
    else:
        # first see if it really depends on symbol and whether there
        # is a linear solution
        f_num, sol = solve_linear(f, symbols=symbols)
        if not symbol in f_num.free_symbols:
            return []
        elif f_num.is_Symbol:
            # no need to check but simplify if desired
            if flags.get('simplify', True):
                sol = simplify(sol)
            return [sol]

        result = False # no solution was obtained
        msg = '' # there is no failure message
        dens = denoms(f, symbols) # store these for checking later

        # Poly is generally robust enough to convert anything to
        # a polynomial and tell us the different generators that it
        # contains, so we will inspect the generators identified by
        # polys to figure out what to do.
        poly = Poly(f_num)
        if poly is None:
            raise ValueError('could not convert %s to Poly' % f_num)
        gens = [g for g in poly.gens if g.has(symbol)]

        if len(gens) > 1:
            # If there is more than one generator, it could be that the
            # generators have the same base but different powers, e.g.
            #   >>> Poly(exp(x)+1/exp(x))
            #   Poly(exp(-x) + exp(x), exp(-x), exp(x), domain='ZZ')
            #   >>> Poly(sqrt(x)+sqrt(sqrt(x)))
            #   Poly(sqrt(x) + x**(1/4), sqrt(x), x**(1/4), domain='ZZ')
            # If the exponents are Rational then a change of variables
            # will make this a polynomial equation in a single base.

            def _as_base_q(x):
                """Return (b**e, q) for x = b**(p*e/q) where p/q is the leading
                Rational of the exponent of x, e.g. exp(-2*x/3) -> (exp(x), 3)
                """
                b, e = x.as_base_exp()
                if e.is_Rational:
                    return b, e.q
                if not e.is_Mul:
                    return x, 1
                c, ee = e.as_coeff_Mul()
                if c.is_Rational and c is not S.One: # c could be a Float
                    return b**ee, c.q
                return x, 1

            bases, qs = zip(*[_as_base_q(g) for g in gens])
            bases = set(bases)

            if len(bases) > 1:
                funcs = set(b.func for b in bases if b.is_Function)

                trig = set([cos, sin, tan, cot])
                other = funcs - trig
                if not other and len(funcs.intersection(trig)) > 1:
                    return _solve(f_num.rewrite(tan), symbol, **flags)

                trigh = set([cosh, sinh, tanh, coth])
                other = funcs - trigh
                if not other and len(funcs.intersection(trigh)) > 1:
                    return _solve(f_num.rewrite(tanh), symbol, **flags)

                # just a simple case - see if replacement of single function
                # clears all symbol-dependent functions, e.g.
                # log(x) - log(log(x) - 1) - 3 can be solved even though it has
                # two generators.

                funcs = [f for f in bases if f.is_Function]
                if funcs:
                    funcs.sort(key=count_ops) # put shallowest function first
                    f1 = funcs[0]
                    t = Dummy('t')
                    # perform the substitution
                    ftry = f_num.subs(f1, t)

                    # if no Functions left, we can proceed with usual solve
                    if not ftry.has(symbol):
                        cv_sols = _solve(ftry, t)
                        cv_inv = _solve(t - f1, symbol)[0]
                        sols = list()
                        for sol in cv_sols:
                            sols.append(cv_inv.subs(t, sol))
                        return sols

                msg = 'multiple generators %s' % gens

            elif any(q != 1 for q in qs):
                # e.g. for x**(1/2) + x**(1/4) a change of variables
                # can be made using p**4 to give p**2 + p
                base = bases.pop()
                m = reduce(ilcm, qs)
                p = Dummy('p', positive=True)
                cov = p**m
                fnew = f_num.subs(base, cov)
                poly = Poly(fnew, p) # we now have a single generator, p

                # for cubics and quartics, if the flag wasn't set, DON'T do it
                # by default since the results are quite long. Perhaps one
                # could base this decision on a certain critical length of the
                # roots.
                if poly.degree() > 2:
                    flags['simplify'] = flags.get('simplify', False)

                soln = roots(poly, cubics=True, quartics=True).keys()
                if not soln:
                    soln = poly.all_roots()
                    check = False # RootOf instances can not be checked

                # We now know what the values of p are equal to. Now find out
                # how they are related to the original x, e.g. if p**2 = cos(x)
                # then x = acos(p**2)
                #
                inversion = _solve(cov - base, symbol, **flags)
                result = [i.subs(p, s) for i in inversion for s in soln]

            else: # len(bases) == 1 and all(q == 1 for q in qs):
                # e.g. case where gens are exp(x), exp(-x)
                u = bases.pop()
                t = Dummy('t')
                inv = _solve(u - t, symbol)
                ftry = f_num.subs(u, t)
                if not ftry.has(symbol):
                    soln = _solve(ftry, t)
                    sols = list()
                    for sol in soln:
                        for i in inv:
                            sols.append(i.subs(t, sol))
                    return sols

        elif len(gens) == 1:

            # There is only one generator that we are interested in, but there
            # may have been more than one generator identified by polys (e.g.
            # for symbols other than the one we are interested in) so recast
            # the poly in terms of our generator of interest.

            if len(poly.gens) > 1:
                poly = Poly(poly, gens[0])

            # if we haven't tried tsolve yet, do so now
            if not flags.pop('tsolve', False):
                # for cubics and quartics, if the flag wasn't set, DON'T do it
                # by default since the results are quite long. Perhaps one
                # could base this decision on a certain critical length of the
                # roots.
                if poly.degree() > 2:
                    flags['simplify'] = flags.get('simplify', False)
                soln = roots(poly, cubics=True, quartics=True).keys()
                if not soln:
                    soln = poly.all_roots()
                    check = False # RootOf instances can not be checked
                gen = poly.gen
                if gen != symbol:
                    u = Dummy()
                    flags['tsolve'] = True
                    inversion = _solve(gen - u, symbol, **flags)
                    soln = list(set([i.subs(u, s) for i in
                                inversion for s in soln]))
                result = soln

    # fallback if above fails
    if result is False:
        result = _tsolve(f_num, symbol, **flags)
        if result is None:
            result = False

    if result is False:
        raise NotImplementedError(msg +
        "\nNo algorithms are implemented to solve equation %s" % f)

    if flags.get('simplify', True):
        result = map(simplify, result)
        # we just simplified the solution so we now set the flag to
        # False so the simplification doesn't happen again in checksol()
        flags['simplify'] = False


    if check:
        # reject any result that makes any denom. affirmatively 0;
        # if in doubt, keep it
        result = [s for s in result if
                    all(not checksol(den, {symbol: s}, **flags)
                    for den in dens)]
        # keep only results if the check is not False
        result = [r for r in result if
                  checksol(f_num, {symbol: r}, **flags) is not False]
    return result


def _solve_system(exprs, symbols, **flags):
    check = flags.get('check', True)
    if not exprs:
        return []

    polys = []
    dens = set()
    failed = []
    result = False
    manual = flags.get('manual', False)
    for j, g in enumerate(exprs):
        dens.update(denoms(g, symbols))
        i, d = _invert(g, *symbols)
        g = d - i
        g = exprs[j] = g.as_numer_denom()[0]
        if manual:
            failed.append(g)
            continue

        poly = g.as_poly(*symbols, **{'extension': True})

        if poly is not None:
            polys.append(poly)
        else:
            failed.append(g)

    if not polys:
        solved_syms = []
    else:
        if all(p.is_linear for p in polys):
            n, m = len(polys), len(symbols)
            matrix = zeros(n, m + 1)

            for i, poly in enumerate(polys):
                for monom, coeff in poly.terms():
                    try:
                        j = list(monom).index(1)
                        matrix[i, j] = coeff
                    except ValueError:
                        matrix[i, m] = -coeff

            # returns a dictionary ({symbols: values}) or None
            result = solve_linear_system(matrix, *symbols, **flags)
            if result:
                # it doesn't need to be checked but we need to see
                # that it didn't set any denominators to 0
                if any(checksol(d, result, **flags) for d in dens):
                    result = None
            if failed:
                if result:
                    solved_syms = result.keys()
                else:
                    solved_syms = []

        else:
            if len(symbols) != len(polys):
                from sympy.utilities.iterables import subsets
                from sympy.core.compatibility import set_union

                free = set_union(*[p.free_symbols for p in polys])
                free = list(free.intersection(symbols))
                free.sort(key=default_sort_key)
                for syms in subsets(free, len(polys)):
                    try:
                        # returns [] or list of tuples of solutions for syms
                        result = solve_poly_system(polys, *syms)
                        if result:
                            solved_syms = syms
                            break
                    except NotImplementedError:
                        pass
                else:
                    raise NotImplementedError('no valid subset found')
            else:
                try:
                    result = solve_poly_system(polys, *symbols)
                    solved_syms = symbols
                except NotImplementedError:
                    failed.extend([g.as_expr() for g in polys])
                    solved_syms = []

            if result:
                # we don't know here if the symbols provided were given
                # or not, so let solve resolve that. A list of dictionaries
                # is going to always be returned from here.
                result = [dict(zip(solved_syms, r)) for r in result]

                checked = []
                warning = flags.get('warn', False)
                for r in result:
                    check = checksol(polys, r, **flags)
                    if check is not False:
                        if check is None and warning:
                            print(filldedent("""
                                \tWarning: could not verify solution %s.""" %
                                result))
                        # if it's a solution to any denom then exclude
                        if not dens or not checksol(dens, r, **flags):
                            checked.append(r)
                result = checked

    if failed:
        # For each failed equation, see if we can solve for one of the
        # remaining symbols from that equation. If so, we update the
        # solution set and continue with the next failed equation,
        # repeating until we are done or we get an equation that can't
        # be solved.
        if result:
            if type(result) is dict:
                result = [result]
        else:
            result = [{}]

        def _ok_syms(e, sort=False):
            rv = (e.free_symbols - solved_syms) & legal
            if sort:
                rv = list(rv)
                rv.sort(key=default_sort_key)
            return rv

        solved_syms = set(solved_syms) # set of symbols we have solved for
        legal = set(symbols) # what we are interested in
        simplify_flag = flags.get('simplify', None)
        do_simplify = flags.get('simplify', True)
        # sort so equation with the fewest potential symbols is first;
        # break ties with count_ops and default_sort_key
        short = sift(failed, lambda x: len(_ok_syms(x)))
        failed = []
        for k in sorted(short, key=default_sort_key):
            failed.extend(sorted(sorted(short[k],
                key=lambda x: x.count_ops()),
                key=default_sort_key))
        for eq in failed:
            newresult = []
            got_s = None
            u = Dummy()
            for r in result:
                # update eq with everything that is known so far
                eq2 = eq.subs(r)
                b = checksol(u, u, eq2, minimal=True)
                if b is not None:
                    if b:
                        newresult.append(r)
                    continue
                # search for a symbol amongst those available that
                # can be solved for
                ok_syms = _ok_syms(eq2, sort=True)
                if not ok_syms:
                    break # skip as it's independent of desired symbols
                for s in ok_syms:
                    try:
                        soln = _solve(eq2, s, **flags)
                    except NotImplementedError:
                        continue
                    # put each solution in r and append the now-expanded
                    # result in the new result list; use copy since the
                    # solution for s in being added in-place
                    if do_simplify:
                        flags['simplify'] = False # for checksol's sake
                    for sol in soln:
                        # check that it satisfies *other* equations
                        if check:
                            ok = False
                            for p in polys:
                                if checksol(p, s, sol, **flags) is False:
                                    break
                            else:
                                ok = True
                            if not ok:
                                continue
                            if any(checksol(d, s, sol, **flags) for d in dens):
                                continue
                        # update existing solutions with this new one
                        rnew = r.copy()
                        for k, v in r.iteritems():
                            rnew[k] = v.subs(s, sol)
                        # and add this new solution
                        rnew[s] = sol
                        newresult.append(rnew)
                    if simplify_flag is not None:
                        flags['simplify'] = simplify_flag
                    got_s = s
                    break
                else:
                    raise NotImplementedError('could not solve %s' % eq2)
            if got_s:
                result = newresult
                solved_syms.add(got_s)
        # if there is only one result should we return just the dictionary?
    return result


def solve_linear(lhs, rhs=0, symbols=[], exclude=[]):
    r""" Return a tuple containing derived from f = lhs - rhs that is either:

        (numerator, denominator) of ``f``
            If this comes back as (0, 1) it means
            that ``f`` is independent of the symbols in ``symbols``, e.g::

                y*cos(x)**2 + y*sin(x)**2 - y = y*(0) = 0
                cos(x)**2 + sin(x)**2 = 1

            If it comes back as (0, 0) there is no solution to the equation
            amongst the symbols given.

            If the numerator is not zero then the function is guaranteed
            to be dependent on a symbol in ``symbols``.

        or

        (symbol, solution) where symbol appears linearly in the numerator of
        ``f``, is in ``symbols`` (if given) and is not in ``exclude`` (if given).

        No simplification is done to ``f`` other than and mul=True expansion,
        so the solution will correspond strictly to a unique solution.

    Examples
    ========

    >>> from sympy.solvers.solvers import solve_linear
    >>> from sympy.abc import x, y, z

    These are linear in x and 1/x:

    >>> solve_linear(x + y**2)
    (x, -y**2)
    >>> solve_linear(1/x - y**2)
    (x, y**(-2))

    When not linear in x or y then the numerator and denominator are returned.

    >>> solve_linear(x**2/y**2 - 3)
    (x**2 - 3*y**2, y**2)

    If the numerator is a symbol then (0, 0) is returned if the solution for
    that symbol would have set any denominator to 0:

    >>> solve_linear(1/(1/x - 2))
    (0, 0)
    >>> 1/(1/x) # to SymPy, this looks like x ...
    x
    >>> solve_linear(1/(1/x)) # so a solution is given
    (x, 0)

    If x is allowed to cancel, then this appears linear, but this sort of
    cancellation is not done so the solution will always satisfy the original
    expression without causing a division by zero error.

    >>> solve_linear(x**2*(1/x - z**2/x))
    (x**2*(-z**2 + 1), x)

    You can give a list of what you prefer for x candidates:

    >>> solve_linear(x + y + z, symbols=[y])
    (y, -x - z)

    You can also indicate what variables you don't want to consider:

    >>> solve_linear(x + y + z, exclude=[x, z])
    (y, -x - z)

    If only x was excluded then a solution for y or z might be obtained.

    """
    from sympy import Equality
    if isinstance(lhs, Equality):
        if rhs:
            raise ValueError(filldedent('''
            If lhs is an Equality, rhs must be 0 but was %s''' % rhs))
        rhs = lhs.rhs
        lhs = lhs.lhs
    dens = None
    eq = lhs - rhs
    n, d = eq.as_numer_denom()
    if not n:
        return S.Zero, S.One

    free = n.free_symbols
    if not symbols:
        symbols = free
    else:
        bad = [s for s in symbols if not s.is_Symbol]
        if bad:
            if len(bad) == 1:
                bad = bad[0]
            if len(symbols) == 1:
                eg = 'solve(%s, %s)' % (eq, symbols[0])
            else:
                eg = 'solve(%s, *%s)' % (eq, list(symbols))
            raise ValueError(filldedent('''
                solve_linear only handles symbols, not %s. To isolate
                non-symbols use solve, e.g. >>> %s <<<.
                             ''' % (bad, eg)))
        symbols = free.intersection(symbols)
    symbols = symbols.difference(exclude)

    # derivatives are easy to do but tricky to analyze to see if they are going
    # to disallow a linear solution, so for simplicity we just evaluate the
    # ones that have the symbols of interest
    derivs = defaultdict(list)
    for der in n.atoms(Derivative):
        csym = der.free_symbols & symbols
        for c in csym:
            derivs[c].append(der)

    if symbols:
        all_zero = True
        for xi in symbols:
            # if there are derivatives in this var, calculate them now
            if type(derivs[xi]) is list:
                derivs[xi] = dict([(der, der.doit()) for der in derivs[xi]])
            nn = n.subs(derivs[xi])
            dn = nn.diff(xi)
            if dn:
                all_zero = False
                if not xi in dn.free_symbols:
                    vi = -(nn.subs(xi, 0))/dn
                    if dens is None:
                        dens = denoms(eq, symbols)
                    if not any(checksol(di, {xi: vi}, minimal=True) is True
                              for di in dens):
                        # simplify any trivial integral
                        irep = [(i, i.doit()) for i in vi.atoms(C.Integral) if
                                i.function.is_number]
                        # do a slight bit of simplification
                        vi = expand_mul(vi.subs(irep))
                        return xi, vi

        if all_zero:
            return S.Zero, S.One
    if n.is_Symbol: # there was no valid solution
        n = d = S.Zero
    return n, d # should we cancel now?


def solve_linear_system(system, *symbols, **flags):
    r"""
    Solve system of N linear equations with M variables, which means
    both Cramer and over defined systems are supported. The possible
    number of solutions is zero, one or infinite. Respectively, this
    procedure will return None or dictionary with solutions. In the
    case of over-defined systems all arbitrary parameters are skipped.
    This may cause situation in which an empty dictionary is returned.
    In this case it means all symbols can be assigned arbitrary values.

    Input to this functions is a Nx(M+1) matrix, which means it has
    to be in augmented form. If you prefer to enter N equations and M
    unknowns then use `solve(Neqs, *Msymbols)` instead. Note: a local
    copy of the matrix is made by this routine so the matrix that is
    passed will not be modified.

    The algorithm used here is fraction-free Gaussian elimination,
    which results, after elimination, in an upper-triangular matrix.
    Then solutions are found using back-substitution. This approach
    is more efficient and compact than the Gauss-Jordan method.

    >>> from sympy import Matrix, solve_linear_system
    >>> from sympy.abc import x, y

    Solve the following system::

           x + 4 y ==  2
        -2 x +   y == 14

    >>> system = Matrix(( (1, 4, 2), (-2, 1, 14)))
    >>> solve_linear_system(system, x, y)
    {x: -6, y: 2}

    """
    matrix = system[:, :]
    syms = list(symbols)

    i, m = 0, matrix.cols - 1  # don't count augmentation

    while i < matrix.rows:
        if i == m:
            # an overdetermined system
            if any(matrix[i:, m]):
                return None   # no solutions
            else:
                # remove trailing rows
                matrix = matrix[:i, :]
                break

        if not matrix[i, i]:
            # there is no pivot in current column
            # so try to find one in other columns
            for k in xrange(i + 1, m):
                if matrix[i, k]:
                    break
            else:
                if matrix[i, m]:
                    # we need to know if this is always zero or not. We
                    # assume that if there are free symbols that it is not
                    # identically zero (or that there is more than one way
                    # to make this zero. Otherwise, if there are none, this
                    # is a constant and we assume that it does not simplify
                    # to zero XXX are there better ways to test this?
                    if not matrix[i, m].free_symbols:
                        return None # no solution

                    # zero row with non-zero rhs can only be accepted
                    # if there is another equivalent row, so look for
                    # them and delete them
                    nrows = matrix.rows
                    rowi = matrix.row(i)
                    ip = None
                    j = i + 1
                    while j < matrix.rows:
                        # do we need to see if the rhs of j
                        # is a constant multiple of i's rhs?
                        rowj = matrix.row(j)
                        if rowj == rowi:
                            matrix.row_del(j)
                        elif rowj[:-1] == rowi[:-1]:
                            if ip is None:
                                _, ip = rowi[-1].as_content_primitive()
                            _, jp = rowj[-1].as_content_primitive()
                            if not (simplify(jp - ip) or simplify(jp + ip)):
                                matrix.row_del(j)

                        j += 1

                    if nrows == matrix.rows:
                        # no solution
                        return None
                # zero row or was a linear combination of
                # other rows or was a row with a symbolic
                # expression that matched other rows, e.g. [0, 0, x - y]
                # so now we can safely skip it
                matrix.row_del(i)
                if not matrix:
                    return None
                continue

            # we want to change the order of colums so
            # the order of variables must also change
            syms[i], syms[k] = syms[k], syms[i]
            matrix.col_swap(i, k)

        pivot_inv = S.One/matrix[i, i]

        # divide all elements in the current row by the pivot
        matrix.row(i, lambda x, _: x * pivot_inv)

        for k in xrange(i + 1, matrix.rows):
            if matrix[k, i]:
                coeff = matrix[k, i]

                # subtract from the current row the row containing
                # pivot and multiplied by extracted coefficient
                matrix.row(k, lambda x, j: simplify(x - matrix[i, j]*coeff))

        i += 1

    # if there weren't any problems, augmented matrix is now
    # in row-echelon form so we can check how many solutions
    # there are and extract them using back substitution

    do_simplify = flags.get('simplify', True)

    if len(syms) == matrix.rows:
        # this system is Cramer equivalent so there is
        # exactly one solution to this system of equations
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # run back-substitution for variables
            for j in xrange(k + 1, m):
                content -= matrix[k, j]*solutions[syms[j]]

            if do_simplify:
                solutions[syms[k]] = simplify(content)
            else:
                solutions[syms[k]] = content

            k -= 1

        return solutions
    elif len(syms) > matrix.rows:
        # this system will have infinite number of solutions
        # dependent on exactly len(syms) - i parameters
        k, solutions = i - 1, {}

        while k >= 0:
            content = matrix[k, m]

            # run back-substitution for variables
            for j in xrange(k + 1, i):
                content -= matrix[k, j]*solutions[syms[j]]

            # run back-substitution for parameters
            for j in xrange(i, m):
                content -= matrix[k, j]*syms[j]

            if do_simplify:
                solutions[syms[k]] = simplify(content)
            else:
                solutions[syms[k]] = content

            k -= 1

        return solutions
    else:
        return []   # no solutions


def solve_undetermined_coeffs(equ, coeffs, sym, **flags):
    """Solve equation of a type p(x; a_1, ..., a_k) == q(x) where both
       p, q are univariate polynomials and f depends on k parameters.
       The result of this functions is a dictionary with symbolic
       values of those parameters with respect to coefficients in q.

       This functions accepts both Equations class instances and ordinary
       SymPy expressions. Specification of parameters and variable is
       obligatory for efficiency and simplicity reason.

       >>> from sympy import Eq
       >>> from sympy.abc import a, b, c, x
       >>> from sympy.solvers import solve_undetermined_coeffs

       >>> solve_undetermined_coeffs(Eq(2*a*x + a+b, x), [a, b], x)
       {a: 1/2, b: -1/2}

       >>> solve_undetermined_coeffs(Eq(a*c*x + a+b, x), [a, b], x)
       {a: 1/c, b: -1/c}

    """
    if isinstance(equ, Equality):
        # got equation, so move all the
        # terms to the left hand side
        equ = equ.lhs - equ.rhs

    equ = cancel(equ).as_numer_denom()[0]

    system = collect(equ.expand(), sym, evaluate=False).values()

    if not any(equ.has(sym) for equ in system):
        # consecutive powers in the input expressions have
        # been successfully collected, so solve remaining
        # system using Gaussian elimination algorithm
        return solve(system, *coeffs, **flags)
    else:
        return None # no solutions


def solve_linear_system_LU(matrix, syms):
    """
    Solves the augmented matrix system using LUsolve and returns a dictionary
    in which solutions are keyed to the symbols of syms *as ordered*.

    The matrix must be invertible.

    Examples
    ========

    >>> from sympy import Matrix
    >>> from sympy.abc import x, y, z
    >>> from sympy.solvers.solvers import solve_linear_system_LU

    >>> solve_linear_system_LU(Matrix([
    ... [1, 2, 0, 1],
    ... [3, 2, 2, 1],
    ... [2, 0, 0, 1]]), [x, y, z])
    {x: 1/2, y: 1/4, z: -1/2}

    See Also
    ========

    sympy.matrices.LUsolve

    """
    assert matrix.rows == matrix.cols - 1
    A = matrix[:matrix.rows, :matrix.rows]
    b = matrix[:, matrix.cols - 1:]
    soln = A.LUsolve(b)
    solutions = {}
    for i in range(soln.rows):
        solutions[syms[i]] = soln[i, 0]
    return solutions

_x = Dummy('x')
_a, _b, _c, _d, _e, _f, _g, _h = [Wild(t, exclude=[_x]) for t in 'abcdefgh']
_patterns = None


def _generate_patterns():
    """
    Generates patterns for transcendental equations.

    This is lazily calculated (called) in the tsolve() function and stored in
    the patterns global variable.
    """

    tmp1 = _f ** (_h-(_c*_g/_b))
    tmp2 = (-_e*tmp1/_a)**(1/_d)
    global _patterns
    _patterns = [
        (_a*_x*exp(_b*_x) - _c,
            LambertW(_b*_c/_a)/_b),
        (_a*_x*log(_b*_x) - _c,
            exp(LambertW(_b*_c/_a))/_b),
        (_a*(_b*_x+_c)**_d + _e   ,
            ((-(_e/_a))**(1/_d)-_c)/_b),
        (_b+_c*exp(_d*_x+_e) ,
            (log(-_b/_c)-_e)/_d),
        (_a*_x+_b+_c*exp(_d*_x+_e) ,
            -_b/_a-LambertW(_c*_d*exp(_e-_b*_d/_a)/_a)/_d),
        (_b+_c*_f**(_d*_x+_e) ,
            (log(-_b/_c)-_e*log(_f))/_d/log(_f)),
        (_a*_x+_b+_c*_f**(_d*_x+_e) ,
            -_b/_a-LambertW(_c*_d*_f**(_e-_b*_d/_a)*log(_f)/_a)/_d/log(_f)),
        (_b+_c*log(_d*_x+_e) ,
            (exp(-_b/_c)-_e)/_d),
        (_a*_x+_b+_c*log(_d*_x+_e) ,
            -_e/_d+_c/_a*LambertW(_a/_c/_d*exp(-_b/_c+_a*_e/_c/_d))),
        (_a*(_b*_x+_c)**_d + _e*_f**(_g*_x+_h) ,
            -_c/_b-_d*LambertW(-tmp2*_g*log(_f)/_b/_d)/_g/log(_f))
    ]


def tsolve(eq, sym):
    SymPyDeprecationWarning(
    feature="tsolve()",
    useinstead="solve()"
    ).warn()
    return _tsolve(eq, sym)


def _tsolve(eq, sym, **flags):
    """
    Helper for _solve that solves a transcendental equation with respect
    to the given symbol. Various equations containing powers and logarithms,
    can be solved.

    Only a single solution is returned. This solution is generally
    not unique. In some cases, a complex solution may be returned
    even though a real solution exists.

        >>> from sympy import log
        >>> from sympy.solvers.solvers import _tsolve as tsolve
        >>> from sympy.abc import x

        >>> tsolve(3**(2*x+5)-4, x)
        [(-5*log(3) + 2*log(2))/(2*log(3))]

        >>> tsolve(log(x) + 2*x, x)
        [LambertW(2)/2]

    """
    if _patterns is None:
        _generate_patterns()
    eq2 = eq.subs(sym, _x)
    for p, sol in _patterns:
        m = eq2.match(p)
        if m:
            soln = sol.subs(m).subs(_x, sym)
            if sym not in soln.free_symbols:
                return [soln]

    try:
        u = unrad(eq, sym)
    except ValueError:
        raise NotImplementedError('Radicals cannot be cleared from %s' % eq)
    if u:
        eq, cov, dens = u
        if cov:
            if len(cov) > 1:
                raise NotImplementedError('Not sure how to handle this.')
            isym, ieq = cov[0]
            # since cov is written in terms of positive symbols, set
            # check to False or else 0 would be excluded; _solve will check
            # the results
            flags['check'] = False
            sol = _solve(eq, isym, **flags)
            inv = _solve(ieq, sym, **flags)
            result = []
            for s in sol:
                for i in inv:
                    result.append(i.subs(isym, s))
            return result
        else:
            return _solve(eq, sym, **flags)

    rhs, lhs = _invert(eq, sym)

    if lhs.is_Add:
        # it's time to try factoring
        fac = factor(lhs - rhs)
        if fac.is_Mul:
            return _solve(fac, sym)

    elif lhs.is_Pow:
        if lhs.exp.is_Integer:
            if lhs - rhs != eq:
                return _solve(lhs - rhs, sym)
            elif not rhs:
                return _solve(lhs.base, sym)
        elif sym not in lhs.exp.free_symbols:
            return _solve(lhs.base - rhs**(1/lhs.exp), sym)
        elif not rhs and sym in lhs.exp.free_symbols:
            # f(x)**g(x) only has solutions where f(x) == 0 and g(x) != 0 at
            # the same place
            sol_base = _solve(lhs.base, sym)
            if not sol_base:
                return sol_base
            return list(set(sol_base) - set(_solve(lhs.exp, sym)))
        elif (rhs is not S.Zero and
              lhs.base.is_positive and
              lhs.exp.is_real):
            return _solve(lhs.exp*log(lhs.base) - log(rhs), sym)

    elif lhs.is_Mul and rhs.is_positive:
        llhs = expand_log(log(lhs))
        if llhs.is_Add:
            return _solve(llhs - log(rhs), sym)

    rewrite = lhs.rewrite(exp)
    if rewrite != lhs:
        return _solve(rewrite - rhs, sym)

    if flags.pop('force', True):
        flags['force'] = False
        pos, reps = posify(lhs - rhs)
        for u, s in reps.iteritems():
            if s == sym:
                break
        else:
            u = sym
        try:
            soln = _solve(pos, u, **flags)
        except NotImplementedError:
            return
        return [s.subs(reps) for s in soln]

# TODO: option for calculating J numerically


def nsolve(*args, **kwargs):
    r"""
    Solve a nonlinear equation system numerically::

        nsolve(f, [args,] x0, modules=['mpmath'], **kwargs)

    f is a vector function of symbolic expressions representing the system.
    args are the variables. If there is only one variable, this argument can
    be omitted.
    x0 is a starting vector close to a solution.

    Use the modules keyword to specify which modules should be used to
    evaluate the function and the Jacobian matrix. Make sure to use a module
    that supports matrices. For more information on the syntax, please see the
    docstring of lambdify.

    Overdetermined systems are supported.

    >>> from sympy import Symbol, nsolve
    >>> import sympy
    >>> sympy.mpmath.mp.dps = 15
    >>> x1 = Symbol('x1')
    >>> x2 = Symbol('x2')
    >>> f1 = 3 * x1**2 - 2 * x2**2 - 1
    >>> f2 = x1**2 - 2 * x1 + x2**2 + 2 * x2 - 8
    >>> print nsolve((f1, f2), (x1, x2), (-1, 1))
    [-1.19287309935246]
    [ 1.27844411169911]

    For one-dimensional functions the syntax is simplified:

    >>> from sympy import sin, nsolve
    >>> from sympy.abc import x
    >>> nsolve(sin(x), x, 2)
    3.14159265358979
    >>> nsolve(sin(x), 2)
    3.14159265358979

    mpmath.findroot is used, you can find there more extensive documentation,
    especially concerning keyword parameters and available solvers.
    """
    # interpret arguments
    if len(args) == 3:
        f = args[0]
        fargs = args[1]
        x0 = args[2]
    elif len(args) == 2:
        f = args[0]
        fargs = None
        x0 = args[1]
    elif len(args) < 2:
        raise TypeError('nsolve expected at least 2 arguments, got %i'
                        % len(args))
    else:
        raise TypeError('nsolve expected at most 3 arguments, got %i'
                        % len(args))
    modules = kwargs.get('modules', ['mpmath'])
    if isinstance(f, (list, tuple)):
        f = Matrix(f).T
    if not isinstance(f, Matrix):
        # assume it's a sympy expression
        if isinstance(f, Equality):
            f = f.lhs - f.rhs
        f = f.evalf()
        atoms = f.atoms(Symbol)
        if fargs is None:
            fargs = atoms.copy().pop()
        if not (len(atoms) == 1 and (fargs in atoms or fargs[0] in atoms)):
            raise ValueError(filldedent('''
                expected a one-dimensional and numerical function'''))

        # the function is much better behaved if there is no denominator
        f = f.as_numer_denom()[0]

        f = lambdify(fargs, f, modules)
        return findroot(f, x0, **kwargs)
    if len(fargs) > f.cols:
        raise NotImplementedError(filldedent('''
            need at least as many equations as variables'''))
    verbose = kwargs.get('verbose', False)
    if verbose:
        print 'f(x):'
        print f
    # derive Jacobian
    J = f.jacobian(fargs)
    if verbose:
        print 'J(x):'
        print J
    # create functions
    f = lambdify(fargs, f.T, modules)
    J = lambdify(fargs, J, modules)
    # solve the system numerically
    x = findroot(f, x0, J=J, **kwargs)
    return x


def _invert(eq, *symbols, **kwargs):
    """Return tuple (i, d) where ``i`` is independent of ``symbols`` and ``d``
    contains symbols. ``i`` and ``d`` are obtained after recursively using
    algebraic inversion until an uninvertible ``d`` remains. If there are no
    free symbols then ``d`` will be zero. Some (but not necessarily all)
    solutions to the expression ``i - d`` will be related the solutions of the
    original expression.

    Examples
    ========

    >>> from sympy.solvers.solvers import _invert as invert
    >>> from sympy import sqrt, cos
    >>> from sympy.abc import x, y
    >>> invert(x - 3)
    (3, x)
    >>> invert(3)
    (3, 0)
    >>> invert(2*cos(x) - 1)
    (pi/3, x)
    >>> invert(sqrt(x) - 3)
    (3, sqrt(x))
    >>> invert(sqrt(x) + y, x)
    (-y, sqrt(x))
    >>> invert(sqrt(x) + y, y)
    (-sqrt(x), y)
    >>> invert(sqrt(x) + y, x, y)
    (0, sqrt(x) + y)

    If there is more than one symbol in a power's base and the exponent
    is not an Integer, then the principal root will be used for the
    inversion:

    >>> invert(sqrt(x + y) - 2)
    (4, x + y)
    >>> invert(sqrt(x + y) - 2)
    (4, x + y)

    If the exponent is an integer, setting ``integer_power`` to True
    will force the principal root to be selected:

    >>> invert(x**2 - 4, integer_power=True)
    (2, x)

    """
    eq = sympify(eq)
    free = eq.free_symbols
    if not symbols:
        symbols = free
    if not free & set(symbols):
        return eq, S.Zero

    dointpow = bool(kwargs.get('integer_power', False))

    inverses = {
    asin: sin,
    acos: cos,
    atan: tan,
    acot: cot,
    asinh: sinh,
    acosh: cosh,
    atanh: tanh,
    acoth: coth,
    }

    lhs = eq
    rhs = S.Zero
    while True:
        was = lhs
        while True:
            indep, dep = lhs.as_independent(*symbols)

            # dep + indep == rhs
            if lhs.is_Add:
                # this indicates we have done it all
                if indep is S.Zero:
                    break

                lhs = dep
                rhs -= indep

            # dep * indep == rhs
            else:
                # this indicates we have done it all
                if indep is S.One:
                    break

                lhs = dep
                rhs /= indep

        # collect like-terms in symbols
        if lhs.is_Add:
            terms = {}
            for a in lhs.args:
                i, d = a.as_independent(*symbols)
                terms.setdefault(d, []).append(i)
            if any(len(v) > 1 for v in terms.values()):
                args = []
                for d, i in terms.iteritems():
                    if len(i) > 1:
                        args.append(Add(*i)*d)
                    else:
                        args.append(i[0]*d)
                lhs = Add(*args)

        # if it's a two-term Add with rhs = 0 and two powers we can get the
        # dependent terms together, e.g. 3*f(x) + 2*g(x) -> f(x)/g(x) = -2/3
        if lhs.is_Add and not rhs and len(lhs.args) == 2:
            a, b = lhs.as_two_terms()
            ai, ad = a.as_independent(*symbols)
            bi, bd = b.as_independent(*symbols)
            if any(_ispow(i) for i in (ad, bd)) and \
               ad.as_base_exp()[0] == bd.as_base_exp()[0]:
                # a = -b
                lhs = powsimp(powdenest(ad/bd))
                rhs = -bi/ai

        elif lhs.is_Mul and any(_ispow(a) for a in lhs.args):
            lhs = powsimp(powdenest(lhs))

        #                    -1
        # f(x) = g  ->  x = f  (g)
        elif lhs.is_Function and (lhs.nargs==1 or len(lhs.args) == 1) and \
                                 (hasattr(lhs, 'inverse') or
                                  lhs.func in inverses or
                                  lhs.func is Abs):
            if lhs.func in inverses:
                inv = inverses[lhs.func]
            elif lhs.func is Abs:
                inv = lambda w: w**2
                lhs = Basic(lhs.args[0]**2) # get it ready to remove the args
            else:
                inv = lhs.inverse()
            rhs = inv(rhs)
            lhs = lhs.args[0]

        if rhs and lhs.is_Pow and lhs.exp.is_Integer and lhs.exp < 0:
            lhs = 1/lhs
            rhs = 1/rhs

        # base**a = b -> base = b**(1/a) if
        #    a is an Integer and dointpow=True (this gives real branch of root)
        #    a is not an Integer and the equation is multivariate and the
        #      base has more than 1 symbol in it
        # The rationale for this is that right now the multi-system solvers
        # doesn't try to resolve generators to see, for example, if the whole
        # system is written in terms of sqrt(x + y) so it will just fail, so we
        # do that step here.
        if lhs.is_Pow and (
           lhs.exp.is_Integer and dointpow or not lhs.exp.is_Integer and
           len(symbols) > 1 and len(lhs.base.free_symbols & set(symbols)) > 1):
            rhs = rhs**(1/lhs.exp)
            lhs = lhs.base

        if lhs == was:
            break
    return rhs, lhs


def unrad(eq, *syms, **flags):
    """ Remove radicals with symbolic arguments and return (eq, cov, dens),
    None or raise an error:

    None is returned if there are no radicals to remove.

    ValueError is raised if there are radicals and they cannot be removed.

    Otherwise

        ``eq``, ``cov``
            equation without radicals, perhaps written in terms of
            change variables; the relationship to the original variables
            is given by the expressions in list (``cov``) whose tuples,
            (``v``, ``expr``) give the change variable introduced (``v``)
            and the expression (``expr``) which equates the base of the radical
            to the power of the change variable needed to clear the radical.
            For example, for sqrt(2 - x) the tuple (_p, -_p**2 - x + 2), would
            be obtained.
        ``dens``
            A set containing all denominators encountered while removing
            radicals. This may be of interest since any solution obtained in
            the modified expression should not set any denominator to zero.
        ``syms``
            an iterable of symbols which, if provided, will limit the focus of
            radical removal: only radicals with one or more of the symbols of
            interest will be cleared.

    ``flags`` are used internally for communication during recursive calls.

    Radicals can be removed from an expression if:
        *   all bases of the radicals are the same; a change of variables is
            done in this case.
        *   if all radicals appear in one term of the expression
        *   there are only 4 terms with sqrt() factors or there are less than
            four terms having sqrt() factors

    Examples
    ========

        >>> from sympy.solvers.solvers import unrad
        >>> from sympy.abc import x
        >>> from sympy import sqrt, Rational
        >>> unrad(sqrt(x)*x**Rational(1,3) + 2)
        (x**5 - 64, [], [])
        >>> unrad(sqrt(x) + (x + 1)**Rational(1,3))
        (x**3 - x**2 - 2*x - 1, [], [])
        >>> unrad(sqrt(x) + x**Rational(1,3) + 2)
        (_p**3 + _p**2 + 2, [(_p, -_p**6 + x)], [])

    """
    if eq.is_Atom:
        return
    cov, dens, nwas = [flags.get(k, v) for k, v in
                       sorted(dict(dens=None, cov=None, n=None).items())]

    def _take(d):
        # see if this is a term that has symbols of interest
        # and merits further processing
        free = d.free_symbols
        if not free:
            return False
        return not syms or free.intersection(syms)

    if dens is None:
        dens = set()
    if cov is None:
        cov = []

    eq = powdenest(eq)
    eq, d = eq.as_numer_denom()
    eq = _mexpand(eq)
    if _take(d):
        dens.add(d)

    if not eq.free_symbols:
        return eq, cov, list(dens)

    poly = eq.as_poly()

    rads = set([g for g in poly.gens if _take(g) and
                g.is_Pow and g.exp.as_coeff_mul()[0].q != 1])

    if not rads:
        return

    depth = sqrt_depth(eq)

    # if all the bases are the same or all the radicals are in one
    # term, this is the lcm of the radical's exponent denominators
    lcm = reduce(ilcm, [r.exp.q for r in rads])

    # find the bases of the radicals
    bases = set([r.as_base_exp()[0] for r in rads])

    # get terms together that have common generators
    drad = dict(zip(rads, range(len(rads))))
    rterms = {(): []}
    args = Add.make_args(poly.as_expr())
    for t in args:
        if _take(t):
            common = set(t.as_poly().gens).intersection(rads)
            key = tuple(sorted([drad[i] for i in common]))
        else:
            key = ()
        rterms.setdefault(key, []).append(t)
    args = Add(*rterms.pop(()))
    rterms = [Add(*rterms[k]) for k in rterms.keys()]
    # the output will depend on the order terms are processed, so
    # make it canonical quickly
    rterms.sort(key=default_sort_key)

    # continue handling
    ok = True
    if len(rterms) == 1:
        eq = rterms[0]**lcm - (-args)**lcm

    elif len(rterms) == 2 and not args:
        eq = rterms[0]**lcm - rterms[1]**lcm

    elif log(lcm, 2).is_Integer and (not args and len(rterms) == 4 or len(rterms) < 4):
        if len(rterms) == 4:
            # (r0+r1)**2 - (r2+r3)**2
            t1, t2, t3, t4 = [t**2 for t in rterms]
            eq = t1 + t2 + 2*rterms[0]*rterms[1] - \
                (t3 + t4 + 2*rterms[2]*rterms[3])
        elif len(rterms) == 3:
            # (r0+r1)**2 - (r2+a)**2
            t1, t2, t3 = [t**2 for t in rterms]
            eq = t1 + t2 + 2*rterms[0]*rterms[1] - \
                (t3 + args**2 + 2*args*rterms[2])
        elif len(rterms) == 2:
            t1, t2 = [t**2 for t in rterms[:2]]
            # r0**2 - (r1+a)**2
            eq = t1 - (t2 + args**2 + 2*args*rterms[1])

    elif len(bases) == 1: # change of variables may work
        ok = False
        covwas = len(cov)
        b = bases.pop()
        for p, bexpr in cov:
            pow = (b - bexpr)
            if pow.is_Pow:
                pb, pe = pow.as_base_exp()
                if pe == lcm and pb == p:
                    p = pb
                    break
        else:
            p = Dummy('p', positive=True)
            cov.append((p, b - p**lcm))
        eq = poly.subs(b, p**lcm).as_expr()
        if not eq.free_symbols.intersection(syms):
            ok = True
        else:
            if len(cov) > covwas:
                cov = cov[:-1]
    else:
        ok = False

    new_depth = sqrt_depth(eq)
    if not ok or (nwas is not None and len(rterms) == nwas and new_depth and new_depth == depth):
        # XXX: XFAIL tests indicate other cases that should be handled.
        raise ValueError('Cannot remove all radicals from %s' % eq)

    neq = unrad(eq, *syms, **dict(cov=cov, dens=dens, n=len(rterms)))
    if neq:
        eq = neq[0]
    if eq.could_extract_minus_sign():
        eq = -eq
    return (_mexpand(eq), cov, list(dens))
