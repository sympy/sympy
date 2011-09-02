""" This module contain solvers for all kinds of equations:

    - algebraic or transcendental, use solve()

    - recurrence, use rsolve()

    - differential, use dsolve()

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
from sympy.simplify import simplify, collect, powsimp, fraction, posify
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

def denoms(eq, symbols=None):
    """Return (recursively) set of all denominators that appear in eq
    that contain any symbol in iterable ``symbols``; if ``symbols`` is
    None (default) then all denominators with symbols will be returned."""
    from sympy.utilities.iterables import preorder_traversal

    symbols = symbols or eq.free_symbols
    dens = set()
    if not symbols or not eq.has(*symbols):
        return dens
    pt = preorder_traversal(eq)
    for e in pt:
        if e.is_Pow or e.func is exp:
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
           'simplify=True (default)'
               simplify solution before substituting into function (which
               will then be simplified regardless of the flag setting)
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

    if is_sequence(f):
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
            # the flag 'simplify=False' is used in solve to avoid
            # simplifying the solution. So if it is set to False there
            # the simplification will not be attempted here, either. But
            # if the simplification is done here then the flag should be
            # set to False so it isn't done again there.
            # FIXME: this can't work, since `flags` is not passed to
            # `checksol()` as a dict, but as keywords.
            # So, any modification to `flags` here will be lost when returning
            # from `checksol()`.
            for k in sol:
                sol[k] = simplify(sympify(sol[k]))
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
                - ``check``, when False, will return all results without checking
                - ``simplify``, when False, will not simplify solutions
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
              and sorted with default_sort_key; the result will the same as
              if the user had provided those symbols. A univariate equation
              will always return a list of solutions; otherwise, a list of
              mappings showing unambiguously the variable that was solved for
              will be returned unless an 'undetermined coefficients' situation
              is detected (see below).

                >>> solve(x - 3)
                [3]
                >>> solve(x**2 - y**2)
                [{x: y}, {x: -y}]
                >>> solve(z**2*x**2 - z**2*y**2)
                [{x: y}, {x: -y}]
                >>> solve(z**2*x - z**2*y**2)
                [{x: y**2}]

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
                    [{x: y**2}]
                    >>> solve(x**2 - y, x, y)
                    [{y: x**2}]

                when undetermined coefficients are identified
                    that are linear
                        >>> solve((a + b)*x - b + 2, a, b)
                        {a: -2, b: 2}

                    that are nonlinear
                        >>> solve((a + b)*x - b**2 + 2, a, b)
                        [(-sqrt(2), sqrt(2)), (sqrt(2), -sqrt(2))]

                if there is no linear solution then the first successful
                attempt for a nonlinear solution will be returned
                    >>> solve(x**2 - y**2, x, y)
                    [{x: y}, {x: -y}]
                    >>> solve(x**2 - y**2/exp(x), x, y)
                    [{y: x*exp(x/2)}, {y: -x*exp(x/2)}]

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

                Warning:
                If no symbols are given for a nonlinear system of equations or
                are given as a set, the solution tuples will contain values for
                the symbols as if the symbols had been sorted with sort_key. In
                the following examples, the solution for x appears first in the
                tuple:
                    >>> solve([x - 2, x**2 + y])
                    [(2, -4)]
                    >>> solve([x - 2, x**2 + f(x)], set([f(x), x]))
                    [(2, -4)]

                Presently, assumptions aren't checked when `solve()` input
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
    #   - the single univariate equation case
    #     since the symbol will have been removed from the solution;
    #   - the nonlinear poly_system since that only support zero-dimensional
    #     systems and those results come back as a list
    if symbol_swapped:
        if type(solution) is dict:
            solution = dict([(swap_back_dict[k], v.subs(swap_back_dict))
                              for k, v in solution.iteritems()])
        elif solution and type(solution) is list and type(solution[0]) is dict:
            for i, sol in enumerate(solution):
                solution[i] = dict([(swap_back_dict[k], v.subs(swap_back_dict))
                              for k, v in sol.iteritems()])

    # Get assumptions about symbols, to filter solutions.
    # Note that if assumptions about a solution can't be verified, it is still returned.
    # XXX: Currently, there are some cases which are not handled,
    # see issue 2098 comment 13: http://code.google.com/p/sympy/issues/detail?id=2098#c13.
    check = flags.get('check', True)
    if not check:
        return solution
    warn = flags.get('warning', False)
    if type(solution) is list:
        if solution:
            unchecked = []
            filtered = []
            if type(solution[0]) is tuple:
                for sol in solution:
                    checked_all_symbols = True
                    for symb, val in zip(symbols, sol):
                        test = check_assumptions(val, **symb.assumptions0)
                        if test is None:
                            checked_all_symbols = False
                        if test is False:
                            break
                    if test is not False:
                        filtered.append(sol)
                    if not checked_all_symbols:
                        unchecked.append(sol)
            elif isinstance(solution[0], dict):
                for s in solution:
                    v = s.values()[0]
                    assumptions = s.keys()[0].assumptions0
                    test = check_assumptions(v, **assumptions)
                    if test is not False:
                        filtered.append(s)
                    if test is None:
                        unchecked.append(s)
            else:
                for sol in solution:
                    test = check_assumptions(sol, **symbols[0].assumptions0)
                    if test is not False:
                        filtered.append(sol)
                    if test is None:
                        unchecked.append(sol)

            solution = filtered
            if warn and unchecked:
                print("\n\tWarning: assumptions concerning following solution(s) can't be checked:"
                      + '\n\t' + ', '.join(str(s) for s in unchecked))

    elif type(solution) is dict:
        checked_all_symbols = True
        for symb, val in solution.iteritems():
            test = check_assumptions(val, **symb.assumptions0)
            if test is None:
                checked_all_symbols = False
            if test is False: # not None nor True
                solution = None
                break
        if warn and not checked_all_symbols:
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

    check  = flags.get('check', True)
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
                n, d = solve_linear(f, symbols=[s])
                if n.is_Symbol:
                    return [{n: cancel(d)}]
                failed.append(s)
            for s in failed:
                try:
                    soln = _solve(f, s, **flags)
                except NotImplementedError:
                    continue
                if soln:
                    return [{s: sol} for sol in soln]
                else:
                    return soln
            else:
                msg = "No algorithms are implemented to solve equation %s"
                raise NotImplementedError(msg % f)

        symbol = symbols[0]

        # build up solutions if f is a Mul
        if f.is_Mul:
            result = set()
            dens = denoms(f, symbols)
            for m in f.args:
                soln = _solve(m, symbol, **flags)
                result.update(set(soln))
            if check:
               result = [s for s in result if all(not checksol(den, {symbol: s}, **flags) for den in dens)]

        elif f.is_Piecewise:
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
            dens = set() # all checking has already been done
        else:
            # first see if it really depends on symbol and whether there
            # is a linear solution
            f_num, sol = solve_linear(f, symbols=symbols)
            if not symbol in f_num.free_symbols:
                return []
            elif f_num.is_Symbol:
                return [cancel(sol)]

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

                def as_base_q(x):
                    """Return (b**e, q) for x = b**(p*e/q) where p/q is the leading
                    Rational of the exponent of x, e.g. exp(-2*x/3) -> (exp(x), 3)
                    """
                    b, e = x.as_base_exp()
                    if e.is_Rational:
                        return b, e.q
                    if not e.is_Mul:
                        return x, 1
                    c, ee = e.as_coeff_Mul()
                    if c.is_Rational and not c is S.One: # c could be a Float
                        return b**ee, c.q
                    return x, 1

                bases, qs = zip(*[as_base_q(g) for g in gens])
                bases = set(bases)
                if len(bases) == 1 and any(q != 1 for q in qs):
                    # e.g. for x**(1/2) + x**(1/4) a change of variables
                    # can be made using p**4 to give p**2 + p
                    base = bases.pop()
                    m = reduce(ilcm, qs)
                    p = Dummy('p', positive=True)
                    cov = p**m
                    fnew = f_num.subs(base, cov)
                    poly = Poly(fnew, p) # we now have a single generator, p

                    # for cubics and quartics, if the flag wasn't set, DON'T do it
                    # by default since the results are quite long. Perhaps one could
                    # base this decision on a certain critical length of the roots.
                    if poly.degree() > 2:
                        flags['simplify'] = flags.get('simplify', False)

                    soln = roots(poly, cubics=True, quartics=True).keys()

                    # We now know what the values of p are equal to. Now find out
                    # how they are related to the original x, e.g. if p**2 = cos(x) then
                    # x = acos(p**2)
                    #
                    inversion = _solve(cov - base, symbol, **flags)
                    result = [i.subs(p, s) for i in inversion for s in soln]
                    if check:
                       result = [r for r in result if checksol(f_num, {symbol: r}, **flags) is not False]
            elif len(gens) == 1:

                # There is only one generator that we are interested in, but there may
                # have been more than one generator identified by polys (e.g. for symbols
                # other than the one we are interested in) so recast the poly in terms
                # of our generator of interest.
                if len(poly.gens) > 1:
                    poly = Poly(poly, gens[0])

                # if we haven't tried tsolve yet, do so now
                if not flags.pop('tsolve', False):
                    # for cubics and quartics, if the flag wasn't set, DON'T do it
                    # by default since the results are quite long. Perhaps one could
                    # base this decision on a certain critical length of the roots.
                    if poly.degree() > 2:
                        flags['simplify'] = flags.get('simplify', False)
                    soln = roots(poly, cubics=True, quartics=True).keys()
                    gen = poly.gen
                    if gen != symbol:
                        u = Dummy()
                        flags['tsolve'] = True
                        inversion = _solve(gen - u, symbol, **flags)
                        soln = list(set([i.subs(u, s) for i in inversion for s in soln]))
                    result = soln
            else:
                msg = 'multiple generators %s' % gens

        # fallback if above fails
        if result is False:
            result = _tsolve(f_num, symbol, **flags) or False

        if result is False:
            raise NotImplementedError(msg +
            "\nNo algorithms are implemented to solve equation %s" % f)

        if flags.get('simplify', True):
            result = map(simplify, result)

        # reject any result that makes any denom. affirmatively 0;
        # if in doubt, keep it
        if check:
           result = [s for s in result if
                     all(not checksol(den, {symbol: s}, **flags)
                     for den in dens)]
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
                matrix = zeros(n, m + 1)

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

def solve_linear(lhs, rhs=0, symbols=[], exclude=[]):
    """ Return a tuple containing derived from f = lhs - rhs that is either:

        (numerator, denominator) of ``f``; if this comes back as (0, 1) it means
            that ``f`` is independent of the symbols in ``symbols``, e.g.
                y*cos(x)**2 + y*sin(x)**2 - y = y*(0) = 0
                cos(x)**2 + sin(x)**2 = 1
            If the numerator is not zero then the function is guaranteed
            to be dependent on a symbol in ``symbols``.

        or

        (symbol, solution) where symbol appears linearly in the numerator of ``f``,
            is in ``symbols`` (if given) and is not in ``exclude`` (if given).

        No simplification is done to ``f`` other than and mul=True expansion, so
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

        >>> solve_linear(x + y + z, symbols=[y])
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
    free = ex.free_symbols
    if not symbols:
        symbols = free
    else:
        symbols = free.intersection(symbols)
    symbols = symbols.difference(exclude)
    d_free = d.free_symbols
    if symbols:
        all_zero = True
        for xi in symbols:
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

    do_simplify = flags.get('simplify', True)

    if len(syms) == matrix.rows:
        # this system is Cramer equivalent so there is
        # exactly one solution to this system of equations
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # run back-substitution for variables
            for j in xrange(k+1, m):
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
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # run back-substitution for variables
            for j in xrange(k+1, i):
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

    if not any(equ.has(sym) for equ in system):
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

_x = Dummy('x')
_a,_b,_c,_d,_e,_f,_g,_h = [Wild(t, exclude=[_x]) for t in 'abcdefgh']
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
    import warnings
    warnings.warn("tsolve is deprecated, use solve.", DeprecationWarning)
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
            if not(soln is S.NaN or
                   soln.has(S.Infinity) or
                   soln.has(S.NegativeInfinity) or
                   sym in soln.free_symbols):
                return [soln]

    # let's also try to invert the equation
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

        # find first term which is a Function
        for f1 in lhs.args:
            if f1.is_Function:
                ok = True
                break
        else:
            ok = False # didn't find a function
        if ok:
            # perform the substitution
            lhs_ = lhs.subs(f1, t)

            # if no Functions left, we can proceed with usual solve
            if not (lhs_.is_Function or
                    any(term.is_Function for term in lhs_.args)):
                cv_sols = solve(lhs_ - rhs, t)
                for sol in cv_sols:
                    if sol.has(sym):
                        # there is more than one function
                        break
                else:
                    cv_inv = solve( t - f1, sym )[0]
                    sols = list()
                    for sol in cv_sols:
                        sols.append(cv_inv.subs(t, sol))
                    return sols

    if flags.pop('posify', True):
        flags['posify'] = False
        pos, reps = posify(lhs)
        u = sym
        for u, s in reps.iteritems():
            if s == sym:
                break
        try:
            soln = _solve(pos - rhs, u, **flags)
        except NotImplementedError:
            return
        return [s.subs(reps) for s in soln]

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
