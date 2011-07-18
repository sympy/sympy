""" This module contain solvers for all kinds of equations:

    - algebraic, use solve()

    - recurrence, use rsolve()

    - differential, use dsolve()

    - transcendental, use tsolve()

    - nonlinear (numerically), use nsolve()
      (you will need a good starting point)

"""

from sympy.core.compatibility import iterable, is_sequence
from sympy.core.sympify import sympify
from sympy.core import S, Mul, Add, Pow, Symbol, Wild, Equality, Dummy, Basic
from sympy.core.numbers import ilcm
from sympy.core.relational import Relational
from sympy.logic.boolalg import And, Or

from sympy.functions import log, exp, LambertW
from sympy.simplify import simplify, collect, powsimp
from sympy.matrices import Matrix, zeros
from sympy.polys import roots, cancel, Poly, together
from sympy.functions.elementary.piecewise import piecewise_fold

from sympy.utilities.lambdify import lambdify
from sympy.mpmath import findroot

from sympy.solvers.polysys import solve_poly_system
from sympy.solvers.inequalities import reduce_inequalities

from sympy.core.compatibility import reduce

from sympy.assumptions import Q, ask

from warnings import warn
from types import GeneratorType

def denoms(eq, x=None):
    """Return (recursively) set of all denominators that appear in eq
    that contain any symbol in x; if x is None (default) then all
    denominators with symbols will be returned."""
    from sympy.utilities.iterables import preorder_traversal

    if x is None:
        x = eq.free_symbols
    dens = set()
    pt = preorder_traversal(eq)
    for e in pt:
        if e.is_Pow or e.func is exp:
            n, d = e.as_numer_denom()
            if d in dens:
                pt.skip()
            elif d.has(*x):
                dens.add(d.as_base_exp()[0])
    return dens

def checksol(f, symbol, sol=None, **flags):
    """Checks whether sol is a solution of equation f == 0.

    Input can be either a single symbol and corresponding value
    or a dictionary of symbols and values.

    Examples:
    ---------

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
               do a fast numerical check if f has only one symbol.
           'minimal=True (default is False)'
               a very fast, minimal testing.
           'warning=True (default is False)'
               print a warning if checksol() could not conclude.
           'simplified=True (default)'
               solution should be simplified before substituting into function
               and function should be simplified after making substitution.
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

    if hasattr(f, '__iter__') and hasattr(f, '__len__'):
        if not f:
            raise ValueError('no functions to check')
        rv = set()
        for fi in f:
            check = checksol(fi, sol, **flags)
            if check is False:
                return False
            rv.add(check)
        if None in rv: # rv might contain True and/or None
            return None
        assert len(rv) == 1 # True
        return True

    if isinstance(f, Poly):
        f = f.as_expr()
    elif isinstance(f, Equality):
        f = f.lhs - f.rhs

    if not f:
        return True

    if not f.has(*sol.keys()):
        return False

    attempt = -1
    numerical = flags.get('numerical', True)
    while 1:
        attempt += 1
        if attempt == 0:
            val = f.subs(sol)
        elif attempt == 1:
            if not val.atoms(Symbol) and numerical:
                # val is a constant, so a fast numerical test may suffice
                if val not in [S.Infinity, S.NegativeInfinity]:
                    # issue 2088 shows that +/-oo chops to 0
                    val = val.evalf(36).n(30, chop=True)
        elif attempt == 2:
            if flags.get('minimal', False):
                return
            # the flag 'simplified=False' is used in solve to avoid
            # simplifying the solution. So if it is set to False there
            # the simplification will not be attempted here, either. But
            # if the simplification is done here then the flag should be
            # set to False so it isn't done again there.
            # FIXME: this can't work, since `flags` is not passed to
            # `checksol()` as a dict, but as keywords.
            # So, any modification to `flags` here will be lost when returning
            # from `checksol()`.
            if flags.get('simplified', True):
                for k in sol:
                    sol[k] = simplify(sympify(sol[k]))
                flags['simplified'] = False
                val = simplify(f.subs(sol))
            if flags.get('force', False):
                val = posify(val)[0]
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
            break
        if val.is_zero:
            return True
        elif attempt > 0 and numerical and val.is_nonzero:
            return False

    if flags.get('warning', False):
        print("\n\tWarning: could not verify solution %s." % sol)
    # returns None if it can't conclude
    # TODO: improve solution testing

def check_assumptions(expr, **assumptions):
    """Checks whether expression `expr` satisfies all assumptions.

    `assumptions` is a dict of assumptions: {'assumption': True|False, ...}.

    Examples:
    ---------

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
        assert isinstance(expected, bool), 'Argument %s=%s is incorrect. \
                                            A boolean is expected.' %(key, expected)
        if hasattr(Q, key):
            test = ask(getattr(Q, key)(expr))
            if test is expected:
                continue
            elif test is not None:
                return False
        # ask() can't conclude. Try using old assumption system.
        # XXX: remove this once transition to new assumption system is finished.
        test = getattr(expr, 'is_' + key, None)
        if test is expected:
            continue
        elif test is not None:
            return False
        result = None # Can't conclude, unless an other test fails.
    return result

# Codes for guess solve strategy
GS_POLY = 0
GS_RATIONAL = 1
GS_POLY_CV_1 = 2 # can be converted to a polynomial equation via the change of variable y -> x**a, a real
GS_POLY_CV_2 = 3 # can be converted to a polynomial equation multiplying on both sides by x**m
                 # for example, x + 1/x == 0. Multiplying by x yields x**2 + x == 0
GS_RATIONAL_CV_1 = 4 # can be converted to a rational equation via the change of variable y -> x**n
GS_PIECEWISE = 5
GS_TRANSCENDENTAL = 6

def guess_solve_strategy(expr, symbol):
    """
    Tries to guess what approach should be used to solve a specific equation

    Returns
    =======
       - -1: could not guess
       - integer > 0: code representing certain type of equation. See GS_* fields
         on this module for a complete list

    Examples
    ========
    >>> from sympy import Symbol, Rational
    >>> from sympy.solvers.solvers import guess_solve_strategy
    >>> from sympy.abc import x
    >>> guess_solve_strategy(x**2 + 1, x)
    0
    >>> guess_solve_strategy(x**Rational(1,2) + 1, x)
    2

    """
    eq_type = -1
    if expr.is_Add:
        return max([guess_solve_strategy(i, symbol) for i in expr.args])

    elif expr.is_Mul:
        # check for rational functions
        num, denom = expr.as_numer_denom()
        if denom.has(symbol):
            #we have a quotient
            m = max(guess_solve_strategy(num, symbol), guess_solve_strategy(denom, symbol))
            if m == GS_POLY:
                return GS_RATIONAL
            elif m == GS_POLY_CV_1:
                return GS_RATIONAL_CV_1
            else:
                raise NotImplementedError
        else:
            return max([guess_solve_strategy(i, symbol) for i in expr.args])

    elif expr.is_Symbol:
        return GS_POLY

    elif expr.is_Pow:
        if expr.exp.has(symbol):
            return GS_TRANSCENDENTAL
        elif not expr.exp.has(symbol) and expr.base.has(symbol):
            if expr.exp.is_Integer and expr.exp > 0:
                eq_type = max(eq_type, GS_POLY)
            elif expr.exp.is_Integer and expr.exp < 0:
                eq_type = max(eq_type, GS_POLY_CV_2)
            elif expr.exp.is_Rational:
                eq_type = max(eq_type, GS_POLY_CV_1)
            else:
                return GS_TRANSCENDENTAL

    elif expr.is_Piecewise:
        return GS_PIECEWISE

    elif expr.is_Function and expr.has(symbol):
        return GS_TRANSCENDENTAL

    elif not expr.has(symbol):
        return GS_POLY

    return eq_type

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
            f
                - a single Expr or Poly that must be zero,
                - an Equality
                - a Relational expression or boolean
                - iterable of one or more of the above

            symbols (Symbol, Function or Derivative) specified as
                - none given (all free symbols will be used)
                - single symbol
                - denested list of symbols
                  e.g. solve(f, x, y)
                - ordered iterable of symbols
                  e.g. solve(f, [x, y])

            flags
                - ``simplified``, when False, will not simplify solutions
                                 (default=True except for polynomials of
                                  order 3 or greater)
                - ``warning``, when True, will warn every time a solution can
                               not be checked, or assumptions about a variable
                               can't be verified for a solution.

        The output varies according to the input and can be seen by example:

            >>> from sympy import solve, Poly, Eq, Function, exp
            >>> from sympy.abc import x, y, z, a, b

            o boolean or univariate Relational

                >>> solve(x < 3)
                And(im(x) == 0, re(x) < 3)

            o single expression and single symbol that is in the expression

                >>> solve(x - y, x)
                [y]
                >>> solve(x - 3, x)
                [3]
                >>> solve(Eq(x, 3), x)
                [3]
                >>> solve(Poly(x - 3), x)
                [3]
                >>> solve(x**2 - y**2, x)
                [y, -y]
                >>> solve(x**4 - 1, x)
                [1, -1, -I, I]

            o single expression with no symbol that is in the expression

                >>> solve(3, x)
                []
                >>> solve(x - 3, y)
                []

            o when no symbol is given then all free symbols will be used
              and sorted with default_sort_key and the result will be the
              same as above as if those symbols had been supplied

                >>> solve(x - 3)
                [3]
                >>> solve(x**2 - y**2)
                [y, -y]

            o when a Function or Derivative is given as a symbol, it is isolated
              algebraically and an implicit solution may be obtained

                >>> f = Function('f')
                >>> solve(f(x) - x, f(x))
                [x]
                >>> solve(f(x).diff(x) - f(x) - x, f(x).diff(x))
                [x + f(x)]

            o single expression and more than 1 symbol

                when there is a linear solution
                    >>> solve(x - y**2, x, y)
                    {x: y**2}
                    >>> solve(x**2 - y, x, y)
                    {y: x**2}

                when undetermined coefficients are identified
                    that are linear
                        >>> solve((a + b)*x - b + 2, a, b)
                        {a: -2, b: 2}

                    that are nonlinear
                        >>> solve((a + b)*x - b**2 + 2, a, b)
                        [(-2**(1/2), 2**(1/2)), (2**(1/2), -2**(1/2))]

                if there is no linear solution then the first successful
                attempt for a nonlinear solution will be returned
                    >>> solve(x**2 - y**2, x, y)
                    [y, -y]
                    >>> solve(x**2 - y**2/exp(x), x, y)
                    [x*exp(x/2), -x*exp(x/2)]

            o iterable of one or more of the above

                involving relationals or bools
                    >>> solve([x < 3, x - 2])
                    And(im(x) == 0, re(x) == 2)
                    >>> solve([x > 3, x - 2])
                    False

                when the system is linear
                    with a solution
                        >>> solve([x - 3], x)
                        {x: 3}
                        >>> solve((x + 5*y - 2, -3*x + 6*y - 15), x, y)
                        {x: -3, y: 1}
                        >>> solve((x + 5*y - 2, -3*x + 6*y - 15), x, y, z)
                        {x: -3, y: 1}
                        >>> solve((x + 5*y - 2, -3*x + 6*y - z), z, x, y)
                        {x: -5*y + 2, z: 21*y - 6}

                    without a solution
                        >>> solve([x + 3, x - 3])

                when the system is not linear
                    >>> solve([x**2 + y -2, y**2 - 4], x, y)
                    [(-2, -2), (0, 2), (0, 2), (2, -2)]

                Warning: there is a possibility of obtaining ambiguous results
                if no symbols are given for a nonlinear system of equations or
                are given as a set since the symbols are not presently reported
                with the solution. A warning will be issued in this situation.
                    >>> solve([x - 2, x**2 + y])
                    <BLANKLINE>
                        For nonlinear systems of equations, symbols should be
                        given as a list so as to avoid ambiguity in the results.
                        solve sorted the symbols as [x, y]
                    [(2, -4)]

                    >>> solve([x - 2, x**2 + f(x)], set([f(x), x]))
                    <BLANKLINE>
                        For nonlinear systems of equations, symbols should be
                        given as a list so as to avoid ambiguity in the results.
                        solve sorted the symbols as [x, f(x)]
                    [(2, -4)]

                If two variables (or more) don't appear in the result, the assumptions
                can't be checked.
                    >>> solve(z**2*x**2 - z**2*y**2/exp(x), x, y, z, warning=True)
                    <BLANKLINE>
                        Warning: assumptions can't be checked
                        (can't find for which variable equation was solved).
                    [x*exp(x/2), -x*exp(x/2)]

                Presently, assumptions aren't checked either when `solve()` input
                involves relationals or bools.

       See also:
          rsolve() for solving recurrence relationships
          dsolve() for solving differential equations

    """
    # make f and symbols into lists of sympified quantities
    # keeping track of how f was passed since if it is a list
    # a dictionary of results will be returned.
    ###########################################################################
    def sympified_list(w):
        return map(sympify, w if iterable(w) else [w])
    bare_f = not iterable(f)
    ordered_symbols = (symbols and
                       symbols[0] and
                       (isinstance(symbols[0], Symbol) or
                        is_sequence(symbols[0], include=GeneratorType)
                       )
                      )
    f, symbols = (sympified_list(w) for w in [f, symbols])

    # preprocess equation(s)
    ###########################################################################
    for i, fi in enumerate(f):
        if isinstance(fi, Equality):
            f[i] = fi.lhs - fi.rhs
        elif isinstance(fi, Poly):
            f[i] = fi.as_expr()
        elif isinstance(fi, bool) or fi.is_Relational:
            return reduce_inequalities(f, assume=flags.get('assume'))
        # Any embedded piecewise functions need to be brought out to the
        # top level so that the appropriate strategy gets selected.
        f[i] = piecewise_fold(f[i])

    # preprocess symbol(s)
    ###########################################################################
    if not symbols:
        # get symbols from equations or supply dummy symbols so solve(3) behaves
        # like solve(3, x).
        symbols = set([])
        for fi in f:
            symbols |= fi.free_symbols or set([Dummy()])
        ordered_symbols = False
    elif len(symbols) == 1 and iterable(symbols[0]):
        symbols = symbols[0]
    if not ordered_symbols:
        # we do this to make the results returned canonical in case f
        # contains a system of nonlinear equations; all other cases should
        # be unambiguous
        symbols = sorted(symbols, key=lambda i: i.sort_key())

    # we can solve for Function and Derivative instances by replacing them
    # with Dummy symbols
    symbols_new = []
    symbol_swapped = False
    symbols_passed = list(symbols)

    for i, s in enumerate(symbols):
        if s.is_Symbol:
            s_new = s
        elif s.is_Function:
            symbol_swapped = True
            s_new = Dummy('F%d' % i)
        elif s.is_Derivative:
            symbol_swapped = True
            s_new = Dummy('D%d' % i)
        else:
            msg = 'expected Symbol, Function or Derivative but got %s'
            raise TypeError(msg % type(s))
        symbols_new.append(s_new)

    if symbol_swapped:
        swap_back_dict = dict(zip(symbols_new, symbols))
        swap_dict = zip(symbols, symbols_new)
        f = [fi.subs(swap_dict) for fi in f]
        symbols = symbols_new

    #
    # try to get a solution
    ###########################################################################
    if bare_f:
        # pass f the way it was passed to solve; if it wasn't a list then
        # a list of solutions will be returned, otherwise a dictionary is
        # going to be returned
        f = f[0]
    solution = _solve(f, *symbols, **flags)

    #
    # postprocessing
    ###########################################################################
    # Restore original Functions and Derivatives if a dictionary is returned.
    # This is not necessary for
    #   - the single equation, single unknown case
    #     since the symbol will have been removed from the solution;
    #   - the nonlinear poly_system since that only support zero-dimensional
    #     systems and those results come back as a list
    if symbol_swapped and type(solution) is dict:
            solution = dict([(swap_back_dict[k], v.subs(swap_back_dict))
                              for k, v in solution.iteritems()])
    # warn if ambiguous results are being obtained
    # XXX agree on how to make this unambiguous
    # see issue 2405 for logic in how Polys chooses ordering and
    # for discussion of what to return see http://groups.google.com/group/sympy
    #                           Apr 18, 2011 posting 'using results from solve'
    elif (not ordered_symbols and
          len(symbols) > 1 and
          solution and
          is_sequence(solution) and
          is_sequence(solution[0]) and
          any(len(set(s)) > 1 for s in solution)
         ):
        msg = ('\n\tFor nonlinear systems of equations, symbols should be' +
               '\n\tgiven as a list so as to avoid ambiguity in the results.' +
               '\n\tsolve sorted the symbols as %s')
        if symbol_swapped:
            from itertools import izip
            tmp = izip(*swap_dict) # separate for the benefit of 2to3
            print msg % list(tmp.next())
        else:
            print msg % symbols

    # Get assumptions about symbols, to filter solutions.
    # Note that if assumptions about a solution can't be verified, it is still returned.
    # XXX: Currently, there are some cases which are not handled,
    # see issue 2098 comment 13: http://code.google.com/p/sympy/issues/detail?id=2098#c13.
    warn = flags.get('warning', False)
    if type(solution) is list:
        if solution:
            unchecked = []
            filtered = []
            if type(solution[0]) is tuple:
                for sol in solution:
                    full_check = True
                    for symb, val in zip(symbols, sol):
                        test = check_assumptions(val, **symb.assumptions0)
                        if test is None:
                            full_check = False
                        if test is False: # not None nor True
                            break
                    if test is not False:
                        filtered.append(sol)
                    if not full_check:
                        unchecked.append(sol)
                solution = filtered
            else:
                if len(symbols) != 1: # find which one was solved for
                    symbols = list(f.free_symbols - set.union(*(s.free_symbols for s in solution)))
                if len(symbols) == 1:
                    for sol in solution:
                        test = check_assumptions(sol, **symbols[0].assumptions0)
                        if test is None:
                            unchecked.append(sol)
                        if test is not False: # None or True
                            filtered.append(sol)
                    solution = filtered
                else:
                    if warn:
                        print("\n\tWarning: assumptions can't be checked"
                              "\n\t(can't find for which variable equation was solved).")
            if warn and unchecked:
                print("\n\tWarning: assumptions concerning following solution(s) can't be checked:"
                      + '\n\t' + ', '.join(str(s) for s in unchecked))

    elif type(solution) is dict:
        full_check = True
        for symb, val in solution.iteritems():
            test = check_assumptions(val, **symb.assumptions0)
            if test is None:
                full_check = False
            if test is False: # not None nor True
                solution = None
                break

        if warn and not full_check:
            print("\n\tWarning: assumptions concerning solution can't be checked.")
    elif isinstance(solution, (Relational, And, Or)):
        assert len(symbols) == 1
        if warn and symbols[0].assumptions0:
            print("\n\tWarning: assumptions about variable '%s' are not handled currently." %symbols[0])
        # TODO: check also variable assumptions for inequalities
    elif solution is not None:
        raise TypeError('Unrecognized solution') # improve the checker to handle this

    #
    # done
    ###########################################################################

    return solution

def _solve(f, *symbols, **flags):
    """ Return a checked solution for f in terms of one or more of the symbols."""

    if not iterable(f):

        if len(symbols) != 1:
            soln = None
            free = f.free_symbols
            ex = free - set(symbols)
            if len(ex) == 1:
                ex = ex.pop()
                try:
                    # may come back as dict or list (if non-linear)
                    soln = solve_undetermined_coeffs(f, symbols, ex)
                except NotImplementedError:
                    pass
            if not soln is None:
                return soln
            # find first successful solution
            failed = []
            for s in symbols:
                n, d = solve_linear(f, x=[s])
                if n.is_Symbol:
                    soln = {n: cancel(d)}
                    return soln
                failed.append(s)
            for s in failed:
                try:
                    soln = _solve(f, s, **flags)
                    return soln
                except NotImplementedError:
                    pass
            else:
                msg = "No algorithms are implemented to solve equation %s"
                raise NotImplementedError(msg % f)

        symbol = symbols[0]

        # first see if it really depends on symbol and whether there
        # is a linear solution
        f_num, sol = solve_linear(f, x=symbols)
        if not symbol in f_num.free_symbols:
            return []
        elif f_num.is_Symbol:
            return [cancel(sol)]

        strategy = guess_solve_strategy(f, symbol)
        result = False # no solution was obtained

        if strategy == GS_POLY:
            poly = f.as_poly(symbol)
            if poly is None:
                msg = "Cannot solve equation %s for %s" % (f, symbol)
            else:
                # for cubics and quartics, if the flag wasn't set, DON'T do it
                # by default since the results are quite long. Perhaps one could
                # base this decision on a certain critical length of the roots.
                if poly.degree() > 2:
                    flags['simplified'] = flags.get('simplified', False)
                result = roots(poly, cubics=True, quartics=True).keys()

        elif strategy == GS_RATIONAL:
            P, _ = f.as_numer_denom()
            dens = denoms(f, x=symbols)
            try:
                soln = _solve(P, symbol, **flags)
            except NotImplementedError:
                msg = "Cannot solve equation %s for %s" % (P, symbol)
                result = []
            else:
                if dens:
                    # reject any result that makes any denom. affirmatively 0;
                    # if in doubt, keep it
                    result = [s for s in soln if all(not checksol(den, {symbol: s}, **flags) for den in dens)]
                else:
                    result = soln

        elif strategy == GS_POLY_CV_1:
            args = list(f.args)
            if isinstance(f, Pow):
                result = _solve(args[0], symbol, **flags)
            elif isinstance(f, Add):
                # we must search for a suitable change of variables
                # collect exponents
                exponents_denom = list()
                for arg in args:
                    if isinstance(arg, Pow):
                        exponents_denom.append(arg.exp.q)
                    elif isinstance(arg, Mul):
                        for mul_arg in arg.args:
                            if isinstance(mul_arg, Pow):
                                exponents_denom.append(mul_arg.exp.q)
                assert len(exponents_denom) > 0
                if len(exponents_denom) == 1:
                    m = exponents_denom[0]
                else:
                    # get the LCM of the denominators
                    m = reduce(ilcm, exponents_denom)
                # x -> y**m.
                # we assume positive for simplification purposes
                t = Dummy('t', positive=True)
                f_ = f.subs(symbol, t**m)
                if guess_solve_strategy(f_, t) != GS_POLY:
                    msg = "Could not convert to a polynomial equation: %s" % f_
                    result = []
                else:
                    soln = [s**m for s in _solve(f_, t)]
                    # we might have introduced solutions from another branch
                    # when changing variables; check and keep solutions
                    # unless they definitely aren't a solution
                    result = [s for s in soln if checksol(f, {symbol: s}, **flags) is not False]

            elif isinstance(f, Mul):
                result = []
                for m in f.args:
                    result.extend(_solve(m, symbol, **flags) or [])

        elif strategy == GS_POLY_CV_2:
            m = 0
            args = list(f.args)
            if isinstance(f, Add):
                for arg in args:
                    if isinstance(arg, Pow):
                        m = min(m, arg.exp)
                    elif isinstance(arg, Mul):
                        for mul_arg in arg.args:
                            if isinstance(mul_arg, Pow):
                                m = min(m, mul_arg.exp)
            elif isinstance(f, Mul):
                for mul_arg in args:
                    if isinstance(mul_arg, Pow):
                        m = min(m, mul_arg.exp)

            if m and m != 1:
                f_ = simplify(f*symbol**(-m))
                try:
                    sols = _solve(f_, symbol)
                except NotImplementedError:
                    msg = 'Could not solve %s for %s' % (f_, symbol)
                else:
                    # we might have introduced unwanted solutions
                    # when multiplying by x**-m; check and keep solutions
                    # unless they definitely aren't a solution
                    if sols:
                        result = [s for s in sols if checksol(f, {symbol: s}, **flags) is not False]
            else:
                msg = 'CV_2 calculated %d but it should have been other than 0 or 1' % m

        elif strategy == GS_PIECEWISE:
            result = set()
            for expr, cond in f.args:
                candidates = _solve(expr, *symbols)
                if isinstance(cond, bool) or cond.is_Number:
                    if not cond:
                        continue

                    # Only include solutions that do not match the condition
                    # of any of the other pieces.
                    for candidate in candidates:
                        matches_other_piece = False
                        for other_expr, other_cond in f.args:
                            if isinstance(other_cond, bool) \
                               or other_cond.is_Number:
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

            result = list(result)

        elif strategy == -1:
            raise ValueError('Could not parse expression %s' % f)

        # this is the fallback for not getting any other solution
        if result is False or strategy == GS_TRANSCENDENTAL:
            soln = tsolve(f_num, symbol)
            dens = denoms(f, x=symbols)
            if not dens:
                result = soln
            else:
                # reject any result that makes any denom. affirmatively 0;
                # if in doubt, keep it
                result = [s for s in soln if all(not checksol(den, {symbol: s}, **flags) for den in dens)]

        if result is False:
            raise NotImplementedError(msg + "\nNo algorithms are implemented to solve equation %s" % f)

        if flags.get('simplified', True) and strategy != GS_RATIONAL:
            result = map(simplify, result)

        return result
    else:
        if not f:
            return []
        else:

            polys = []

            for g in f:

                poly = g.as_poly(*symbols, **{'extension': True})

                if poly is not None:
                    polys.append(poly)
                else:
                    raise NotImplementedError()

            if all(p.is_linear for p in polys):
                n, m = len(f), len(symbols)
                matrix = zeros((n, m + 1))

                for i, poly in enumerate(polys):
                    for monom, coeff in poly.terms():
                        try:
                            j = list(monom).index(1)
                            matrix[i, j] = coeff
                        except ValueError:
                            matrix[i, m] = -coeff

                # a dictionary of symbols: values or None
                result = solve_linear_system(matrix, *symbols, **flags)
                return result
            else:
                # a list of tuples, T, where T[i] [j] corresponds to the ith solution for symbols[j]
                result = solve_poly_system(polys)
                return result

def solve_linear(lhs, rhs=0, x=[], exclude=[]):
    """ Return a tuple containing derived from f = lhs - rhs that is either:

        (numerator, denominator) of f; if this comes back as (0, 1) it means
            that f is independent of the symbols of x, e.g.
                y*cos(x)**2 + y*sin(x)**2 - y = y*(0) = 0
                cos(x)**2 + sin(x)**2 = 1
            If the numerator is not zero then the function is guaranteed
            to be dependent on a symbol in x.

        or

        (symbol, solution) where symbol appears linearly in the numerator of f,
            is in x (if given) and is not in exclude (if given).

        No simplification is done to f other than and mul=True expansion, so
        the solution will correspond strictly to a unique solution.

    Examples:

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

    If x is allowed to cancel, then this appears linear, but this sort of
    cancellation is not done so the solution will always satisfy the original
    expression without causing a division by zero error.

        >>> solve_linear(x**2*(1/x - z**2/x))
        (x**2*(-x*z**2 + x), x**2)

    You can give a list of what you prefer for x candidates:

        >>> solve_linear(x + y + z, x=[y])
        (y, -x - z)

    You can also indicate what variables you don't want to consider:

        >>> solve_linear(x + y + z, exclude=[x, z])
        (y, -x - z)

    If only x was excluded then a solution for y or z might be obtained.

    """
    from sympy import expand_mul, Equality
    if isinstance(lhs, Equality):
        rhs += lhs.rhs
        lhs = lhs.lhs
    n, d = (lhs - rhs).as_numer_denom()
    ex = expand_mul(n)
    if not ex:
        return ex, S.One

    exclude = set(exclude)
    syms = ex.free_symbols
    if not x:
        x = syms
    else:
        x = syms.intersection(x)
    x = x.difference(exclude)
    d_free = d.free_symbols
    if x:
        all_zero = True
        for xi in x:
            dn = n.diff(xi)
            if dn:
                all_zero = False
                if not xi in dn.free_symbols:
                    vi = -(n.subs(xi, 0))/dn
                    if not checksol(d, {xi: vi}, minimal=True) is True:
                        return xi, vi

        if all_zero:
            return S.Zero, S.One
    return n, d # should we cancel now?

def solve_linear_system(system, *symbols, **flags):
    """Solve system of N linear equations with M variables, which means
       both Cramer and over defined systems are supported. The possible
       number of solutions is zero, one or infinite. Respectively this
       procedure will return None or dictionary with solutions. In the
       case of over-defined systems all arbitrary parameters are skipped.
       This may cause situation in which an empty dictionary is returned.
       In this case it means all symbols can be assigned arbitrary values.

       Input to this functions is a Nx(M+1) matrix, which means it has
       to be in augmented form. If you prefer to enter N equations and M
       unknowns then use 'solve(Neqs, *Msymbols)' instead. Note: a local
       copy of the matrix is made by this routine so the matrix that is
       passed will not be modified.

       The algorithm used here is fraction-free Gaussian elimination,
       which results, after elimination, in an upper-triangular matrix.
       Then solutions are found using back-substitution. This approach
       is more efficient and compact than the Gauss-Jordan method.

       >>> from sympy import Matrix, solve_linear_system
       >>> from sympy.abc import x, y

       Solve the following system:

              x + 4 y ==  2
           -2 x +   y == 14

       >>> system = Matrix(( (1, 4, 2), (-2, 1, 14)))
       >>> solve_linear_system(system, x, y)
       {x: -6, y: 2}

    """
    matrix = system[:,:]
    syms = list(symbols)

    i, m = 0, matrix.cols-1  # don't count augmentation

    while i < matrix.rows:
        if i == m:
            # an overdetermined system
            if any(matrix[i:,m]):
                return None   # no solutions
            else:
                # remove trailing rows
                matrix = matrix[:i,:]
                break

        if not matrix[i, i]:
            # there is no pivot in current column
            # so try to find one in other columns
            for k in xrange(i+1, m):
                if matrix[i, k]:
                    break
            else:
                if matrix[i, m]:
                    return None   # no solutions
                else:
                    # zero row or was a linear combination of
                    # other rows so now we can safely skip it
                    matrix.row_del(i)
                    continue

            # we want to change the order of colums so
            # the order of variables must also change
            syms[i], syms[k] = syms[k], syms[i]
            matrix.col_swap(i, k)

        pivot_inv = S.One / matrix [i, i]

        # divide all elements in the current row by the pivot
        matrix.row(i, lambda x, _: x * pivot_inv)

        for k in xrange(i+1, matrix.rows):
            if matrix[k, i]:
                coeff = matrix[k, i]

                # subtract from the current row the row containing
                # pivot and multiplied by extracted coefficient
                matrix.row(k, lambda x, j: simplify(x - matrix[i, j]*coeff))

        i += 1

    # if there weren't any problems, augmented matrix is now
    # in row-echelon form so we can check how many solutions
    # there are and extract them using back substitution

    simplified = flags.get('simplified', True)

    if len(syms) == matrix.rows:
        # this system is Cramer equivalent so there is
        # exactly one solution to this system of equations
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # run back-substitution for variables
            for j in xrange(k+1, m):
                content -= matrix[k, j]*solutions[syms[j]]

            if simplified:
                solutions[syms[k]] = simplify(content)
            else:
                solutions[syms[k]] = content

            k -= 1

        return solutions
    elif len(syms) > matrix.rows:
        # this system will have infinite number of solutions
        # dependent on exactly len(syms) - i parameters
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # run back-substitution for variables
            for j in xrange(k+1, i):
                content -= matrix[k, j]*solutions[syms[j]]

            # run back-substitution for parameters
            for j in xrange(i, m):
                content -= matrix[k, j]*syms[j]

            if simplified:
                solutions[syms[k]] = simplify(content)
            else:
                solutions[syms[k]] = content

            k -= 1

        return solutions
    else:
        return None   # no solutions

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

    if not any([ equ.has(sym) for equ in system ]):
        # consecutive powers in the input expressions have
        # been successfully collected, so solve remaining
        # system using Gaussian elimination algorithm
        return solve(system, *coeffs, **flags)
    else:
        return None # no solutions

def solve_linear_system_LU(matrix, syms):
    """ LU function works for invertible only """
    assert matrix.rows == matrix.cols-1
    A = matrix[:matrix.rows,:matrix.rows]
    b = matrix[:,matrix.cols-1:]
    soln = A.LUsolve(b)
    solutions = {}
    for i in range(soln.rows):
        solutions[syms[i]] = soln[i,0]
    return solutions

x = Dummy('x')
a,b,c,d,e,f,g,h = [Wild(t, exclude=[x]) for t in 'abcdefgh']
patterns = None

def _generate_patterns():
    """
    Generates patterns for transcendental equations.

    This is lazily calculated (called) in the tsolve() function and stored in
    the patterns global variable.
    """

    tmp1 = f ** (h-(c*g/b))
    tmp2 = (-e*tmp1/a)**(1/d)
    global patterns
    patterns = [
        (a*(b*x+c)**d + e   , ((-(e/a))**(1/d)-c)/b),
        (    b+c*exp(d*x+e) , (log(-b/c)-e)/d),
        (a*x+b+c*exp(d*x+e) , -b/a-LambertW(c*d*exp(e-b*d/a)/a)/d),
        (    b+c*f**(d*x+e) , (log(-b/c)-e*log(f))/d/log(f)),
        (a*x+b+c*f**(d*x+e) , -b/a-LambertW(c*d*f**(e-b*d/a)*log(f)/a)/d/log(f)),
        (    b+c*log(d*x+e) , (exp(-b/c)-e)/d),
        (a*x+b+c*log(d*x+e) , -e/d+c/a*LambertW(a/c/d*exp(-b/c+a*e/c/d))),
        (a*(b*x+c)**d + e*f**(g*x+h) , -c/b-d*LambertW(-tmp2*g*log(f)/b/d)/g/log(f))
    ]

def tsolve(eq, sym):
    """
    Solves a transcendental equation with respect to the given
    symbol. Various equations containing mixed linear terms, powers,
    and logarithms, can be solved.

    Only a single solution is returned. This solution is generally
    not unique. In some cases, a complex solution may be returned
    even though a real solution exists.

        >>> from sympy import tsolve, log
        >>> from sympy.abc import x

        >>> tsolve(3**(2*x+5)-4, x)
        [(-5*log(3) + 2*log(2))/(2*log(3))]

        >>> tsolve(log(x) + 2*x, x)
        [LambertW(2)/2]

    """
    if patterns is None:
        _generate_patterns()
    eq = sympify(eq)
    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs
    sym = sympify(sym)
    eq2 = eq.subs(sym, x)
    # First see if the equation has a linear factor
    # In that case, the other factor can contain x in any way (as long as it
    # is finite), and we have a direct solution to which we add others that
    # may be found for the remaining portion.
    r = Wild('r')
    m = eq2.match((a*x+b)*r)
    if m and m[a]:
        return [(-b/a).subs(m).subs(x, sym)] + solve(m[r], x)
    for p, sol in patterns:
        m = eq2.match(p)
        if m:
            soln = sol.subs(m).subs(x, sym)
            if not(soln is S.NaN or
                   soln.has(S.Infinity) or
                   soln.has(S.NegativeInfinity) or
                   sym in soln.free_symbols):
                return [soln]

    # let's also try to inverse the equation
    lhs = eq
    rhs = S.Zero

    while True:
        indep, dep = lhs.as_independent(sym)

        # dep + indep == rhs
        if lhs.is_Add:
            # this indicates we have done it all
            if indep is S.Zero:
                break

            lhs = dep
            rhs-= indep

        # dep * indep == rhs
        else:
            # this indicates we have done it all
            if indep is S.One:
                break

            lhs = dep
            rhs/= indep

    #                    -1
    # f(x) = g  ->  x = f  (g)
    if lhs.is_Function and lhs.nargs==1 and hasattr(lhs, 'inverse'):
        rhs = lhs.inverse() (rhs)
        lhs = lhs.args[0]

        sol = solve(lhs-rhs, sym)
        return sol

    elif lhs.is_Add:
        # just a simple case - we do variable substitution for first function,
        # and if it removes all functions - let's call solve.
        #      x    -x                   -1
        # UC: e  + e   = y      ->  t + t   = y
        t = Dummy('t')
        terms = lhs.args

        # find first term which is Function
        for f1 in lhs.args:
            if f1.is_Function:
                break
        else:
            raise NotImplementedError("Unable to solve the equation" + \
                "(tsolve: at least one Function expected at this point")

        # perform the substitution
        lhs_ = lhs.subs(f1, t)

        # if no Functions left, we can proceed with usual solve
        if not (lhs_.is_Function or
                any(term.is_Function for term in lhs_.args)):
            cv_sols = solve(lhs_ - rhs, t)
            for sol in cv_sols:
                if sol.has(sym):
                    raise NotImplementedError("Unable to solve the equation")
            cv_inv = solve( t - f1, sym )[0]
            sols = list()
            for sol in cv_sols:
                sols.append(cv_inv.subs(t, sol))
            return sols


    raise NotImplementedError("Unable to solve the equation.")

# TODO: option for calculating J numerically
def nsolve(*args, **kwargs):
    """
    Solve a nonlinear equation system numerically.

    nsolve(f, [args,] x0, modules=['mpmath'], **kwargs)

    f is a vector function of symbolic expressions representing the system.
    args are the variables. If there is only one variable, this argument can be
    omitted.
    x0 is a starting vector close to a solution.

    Use the modules keyword to specify which modules should be used to evaluate
    the function and the Jacobian matrix. Make sure to use a module that
    supports matrices. For more information on the syntax, please see the
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
    if isinstance(f,  (list,  tuple)):
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
            raise ValueError('expected a one-dimensional and numerical function')

        # the function is much better behaved if there is no denominator
        f = f.as_numer_denom()[0]

        f = lambdify(fargs, f, modules)
        return findroot(f, x0, **kwargs)
    if len(fargs) > f.cols:
        raise NotImplementedError('need at least as many equations as variables')
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
