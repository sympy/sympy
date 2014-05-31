"""
This module contain solvers for all kinds of equations:

    - algebraic or transcendental, use solve()

    - recurrence, use rsolve()

    - differential, use dsolve()

    - nonlinear (numerically), use nsolve()
      (you will need a good starting point)

"""

from __future__ import print_function, division

from sympy.core.compatibility import (iterable, is_sequence, ordered,
    default_sort_key, reduce, xrange)
from sympy.utilities.exceptions import SymPyDeprecationWarning
from sympy.core.sympify import sympify
from sympy.core import (C, S, Add, Symbol, Wild, Equality, Dummy, Basic,
    Expr, Mul, Pow)
from sympy.core.exprtools import factor_terms
from sympy.core.function import (expand_mul, expand_multinomial, expand_log,
                          Derivative, AppliedUndef, UndefinedFunction, nfloat,
                          count_ops, Function, expand_power_exp, Lambda)
from sympy.core.numbers import ilcm, Float
from sympy.core.relational import Relational, Ge
from sympy.logic.boolalg import And, Or
from sympy.core.basic import preorder_traversal

from sympy.functions import (log, exp, LambertW, cos, sin, tan, cot, cosh,
                             sinh, tanh, coth, acos, asin, atan, acot, acosh,
                             asinh, atanh, acoth, Abs, sign, re, im, arg,
                             sqrt, atan2)
from sympy.functions.elementary.miscellaneous import real_root
from sympy.simplify import (simplify, collect, powsimp, posify, powdenest,
                            nsimplify, denom, logcombine)
from sympy.simplify.sqrtdenest import sqrt_depth, _mexpand
from sympy.simplify.fu import TR1, hyper_as_trig
from sympy.matrices import Matrix, zeros
from sympy.polys import (roots, cancel, factor, Poly, together, RootOf,
    degree, PolynomialError)
from sympy.functions.elementary.piecewise import piecewise_fold, Piecewise

from sympy.utilities.lambdify import lambdify
from sympy.utilities.misc import filldedent
from sympy.utilities.iterables import uniq, generate_bell, flatten

from sympy.mpmath import findroot

from sympy.solvers.polysys import solve_poly_system
from sympy.solvers.inequalities import reduce_inequalities

from sympy.assumptions import Q, ask

from types import GeneratorType
from collections import defaultdict
import warnings


def _ispow(e):
    """Return True if e is a Pow or is exp."""
    return isinstance(e, Expr) and (e.is_Pow or e.func is exp)


def denoms(eq, symbols=None):
    """Return (recursively) set of all denominators that appear in eq
    that contain any symbol in iterable ``symbols``; if ``symbols`` is
    None (default) then all denominators will be returned.

    Examples
    ========

    >>> from sympy.solvers.solvers import denoms
    >>> from sympy.abc import x, y, z
    >>> from sympy import sqrt

    >>> denoms(x/y)
    set([y])

    >>> denoms(x/(y*z))
    set([y, z])

    >>> denoms(3/x + y/z)
    set([x, z])

    >>> denoms(x/2 + y/z)
    set([2, z])
    """

    pot = preorder_traversal(eq)
    dens = set()
    for p in pot:
        den =  denom(p)
        if den is S.One:
            continue
        for d in Mul.make_args(den):
            dens.add(d)
    if not symbols:
        return dens
    rv = []
    for d in dens:
        free = d.free_symbols
        if any(s in free for s in symbols):
            rv.append(d)
    return set(rv)


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
    >>> checksol(x**4 - 1, x, 1)
    True
    >>> checksol(x**4 - 1, x, 0)
    False
    >>> checksol(x**2 + y**2 - 5**2, {x: 3, y: 4})
    True

    To check if an expression is zero using checksol, pass it
    as ``f`` and send an empty dictionary for ``symbol``:

    >>> checksol(x**2 + x - x*(x + 1), {})
    True

    None is returned if checksol() could not conclude.

    flags:
        'numerical=True (default)'
           do a fast numerical check if ``f`` has only one symbol.
        'minimal=True (default is False)'
           a very fast, minimal testing.
        'warn=True (default is False)'
           show a warning if checksol() could not conclude.
        'simplify=True (default)'
           simplify solution before substituting into function and
           simplify the function before trying specific simplifications
        'force=True (default is False)'
           make positive all symbols without assumptions regarding sign.

    """
    from sympy.physics.units import Unit

    minimal = flags.get('minimal', False)

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
            rv = None  # don't return, wait to see if there's a False
        return rv

    if isinstance(f, Poly):
        f = f.as_expr()
    elif isinstance(f, Equality):
        f = f.lhs - f.rhs


    if not f:
        return True

    if sol and not f.has(*list(sol.keys())):
        # if f(y) == 0, x=3 does not set f(y) to zero...nor does it not
        return None

    illegal = set([S.NaN,
               S.ComplexInfinity,
               S.Infinity,
               S.NegativeInfinity])
    if any(sympify(v).atoms() & illegal for k, v in sol.items()):
        return False

    was = f
    attempt = -1
    numerical = flags.get('numerical', True)
    while 1:
        attempt += 1
        if attempt == 0:
            val = f.subs(sol)
            if isinstance(val, Mul):
                val = val.as_independent(Unit)[0]
            if val.atoms() & illegal:
                return False
        elif attempt == 1:
            if val.free_symbols:
                if not val.is_constant(*list(sol.keys()), simplify=not minimal):
                    return False
                # there are free symbols -- simple expansion might work
                _, val = val.as_content_primitive()
                val = expand_mul(expand_multinomial(val))
        elif attempt == 2:
            if minimal:
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
                # issue 5673: nz may be True even when False
                # so these are just hacks to keep a false positive
                # from being returned

                # HACK 1: LambertW (issue 5673)
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
            return bool(abs(val.n(18).n(12, chop=True)) < 1e-9)
        was = val

    if flags.get('warn', False):
        warnings.warn("\n\tWarning: could not verify solution %s." % sol)
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
    for key, expected in assumptions.items():
        if expected is None:
            continue
        test = getattr(expr, 'is_' + key, None)
        if test is expected:
            continue
        elif test is not None:
            return False
        result = None  # Can't conclude, unless an other test fails.
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

    * symbols (object(s) to solve for) specified as
        - none given (other non-numeric objects will be used)
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
            show a warning if checksol() could not conclude.
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
        'particular=True (default is False)'
            instructs solve to try to find a particular solution to a linear
            system with as many zeros as possible; this is very expensive
        'quick=True (default is False)'
            when using particular=True, use a fast heuristic instead to find a
            solution with many zeros (instead of using the very slow method
            guaranteed to find the largest number of zeros possible)

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
        >>> solve(x**2 - y**2, x, set=True)
        ([x], set([(-y,), (y,)]))
        >>> solve(x**4 - 1, x, set=True)
        ([x], set([(-1,), (1,), (-I,), (I,)]))

    * single expression with no symbol that is in the expression

        >>> solve(3, x)
        []
        >>> solve(x - 3, y)
        []

    * single expression with no symbol given

          In this case, all free symbols will be selected as potential
          symbols to solve for. If the equation is univariate then a list
          of solutions is returned; otherwise -- as is the case when symbols are
          given as an iterable of length > 1 -- a list of mappings will be returned.

            >>> solve(x - 3)
            [3]
            >>> solve(x**2 - y**2)
            [{x: -y}, {x: y}]
            >>> solve(z**2*x**2 - z**2*y**2)
            [{x: -y}, {x: y}, {z: 0}]
            >>> solve(z**2*x - z**2*y**2)
            [{x: y**2}, {z: 0}]

    * when an object other than a Symbol is given as a symbol, it is
      isolated algebraically and an implicit solution may be obtained.
      This is mostly provided as a convenience to save one from replacing
      the object with a Symbol and solving for that Symbol. It will only
      work if the specified object can be replaced with a Symbol using the
      subs method.

          >>> solve(f(x) - x, f(x))
          [x]
          >>> solve(f(x).diff(x) - f(x) - x, f(x).diff(x))
          [x + f(x)]
          >>> solve(f(x).diff(x) - f(x) - x, f(x))
          [-x + Derivative(f(x), x)]
          >>> solve(x + exp(x)**2, exp(x), set=True)
          ([exp(x)], set([(-sqrt(-x),), (sqrt(-x),)]))

          >>> from sympy import Indexed, IndexedBase, Tuple, sqrt
          >>> A = IndexedBase('A')
          >>> eqs = Tuple(A[1] + A[2] - 3, A[1] - A[2] + 1)
          >>> solve(eqs, eqs.atoms(Indexed))
          {A[1]: 1, A[2]: 2}

        * To solve for a *symbol* implicitly, use 'implicit=True':

            >>> solve(x + exp(x), x)
            [-LambertW(1)]
            >>> solve(x + exp(x), x, implicit=True)
            [-exp(x)]

        * It is possible to solve for anything that can be targeted with
          subs:

            >>> solve(x + 2 + sqrt(3), x + 2)
            [-sqrt(3)]
            >>> solve((x + 2 + sqrt(3), x + 4 + y), y, x + 2)
            {y: -2 + sqrt(3), x + 2: -sqrt(3)}

        * Nothing heroic is done in this implicit solving so you may end up
          with a symbol still in the solution:

            >>> eqs = (x*y + 3*y + sqrt(3), x + 4 + y)
            >>> solve(eqs, y, x + 2)
            {y: -sqrt(3)/(x + 3), x + 2: (-2*x - 6 + sqrt(3))/(x + 3)}
            >>> solve(eqs, y*x, x)
            {x: -y - 4, x*y: -3*y - sqrt(3)}

        * if you attempt to solve for a number remember that the number
          you have obtained does not necessarily mean that the value is
          equivalent to the expression obtained:

            >>> solve(sqrt(2) - 1, 1)
            [sqrt(2)]
            >>> solve(x - y + 1, 1)  # /!\ -1 is targeted, too
            [x/(y - 1)]
            >>> [_.subs(z, -1) for _ in solve((x - y + 1).subs(-1, z), 1)]
            [-x + y]

        * To solve for a function within a derivative, use dsolve.

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

                >>> solve((a + b)*x - b**2 + 2, a, b, set=True)
                ([a, b], set([(-sqrt(2), sqrt(2)), (sqrt(2), -sqrt(2))]))

        * if there is no linear solution then the first successful
          attempt for a nonlinear solution will be returned

            >>> solve(x**2 - y**2, x, y)
            [{x: -y}, {x: y}]
            >>> solve(x**2 - y**2/exp(x), x, y)
            [{x: 2*LambertW(y/2)}]
            >>> solve(x**2 - y**2/exp(x), y, x)
            [{y: -x*exp(x/2)}, {y: x*exp(x/2)}]

    * iterable of one or more of the above

        * involving relationals or bools

            >>> solve([x < 3, x - 2])
            And(im(x) == 0, re(x) == 2)
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

            >>> solve([x**2 + y -2, y**2 - 4], x, y, set=True)
            ([x, y], set([(-2, -2), (0, 2), (2, -2)]))

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
        >>> solve(sin(x)/x)  # 0 is excluded
        [pi]

    If check=False then a solution to the numerator being zero is found: x = 0.
    In this case, this is a spurious solution since sin(x)/x has the well known
    limit (without dicontinuity) of 1 at x = 0:

        >>> solve(sin(x)/x, check=False)
        [0, pi]

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
        return list(map(sympify, w if iterable(w) else [w]))
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
        elif isinstance(fi, (bool, C.BooleanAtom)) or fi.is_Relational:
            return reduce_inequalities(f, assume=flags.get('assume'),
                                       symbols=symbols)

        # if we have a Matrix, we need to iterate over its elements again
        if f[i].is_Matrix:
            bare_f = False
            f.extend(list(f[i]))
            f[i] = S.Zero

        # if we can split it into real and imaginary parts then do so
        freei = f[i].free_symbols
        if freei and all(s.is_real or s.is_imaginary for s in freei):
            fr, fi = f[i].as_real_imag()
            if fr and fi and not any(i.has(re, im, arg) for i in (fr, fi)) \
                    and fr != fi:
                if bare_f:
                    bare_f = False
                f[i: i + 1] = [fr, fi]

    # preprocess symbol(s)
    ###########################################################################
    if not symbols:
        # get symbols from equations
        symbols = reduce(set.union, [fi.free_symbols
                                     for fi in f], set())
        if len(symbols) < len(f):
            for fi in f:
                pot = preorder_traversal(fi)
                for p in pot:
                    if not (p.is_number or p.is_Add or p.is_Mul) or \
                            isinstance(p, AppliedUndef):
                        flags['dict'] = True  # better show symbols
                        symbols.add(p)
                        pot.skip()  # don't go any deeper
        symbols = list(symbols)
        # supply dummy symbols so solve(3) behaves like solve(3, x)
        for i in range(len(f) - len(symbols)):
            symbols.append(Dummy())

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

    # real/imag handling -----------------------------
    w = Dummy('w')
    piece = Lambda(w, Piecewise((w, Ge(w, 0)), (-w, True)))
    for i, fi in enumerate(f):
        # Abs
        reps = []
        for a in fi.atoms(Abs):
            if not a.has(*symbols):
                continue
            if a.args[0].is_real is None:
                raise NotImplementedError('solving %s when the argument '
                    'is not real or imaginary.' % a)
            reps.append((a, piece(a.args[0]) if a.args[0].is_real else \
                piece(a.args[0]*S.ImaginaryUnit)))
        fi = fi.subs(reps)

        # arg
        _arg = [a for a in fi.atoms(arg) if a.has(*symbols)]
        fi = fi.xreplace(dict(list(zip(_arg,
            [atan(im(a.args[0])/re(a.args[0])) for a in _arg]))))

        # save changes
        f[i] = fi

    # see if re(s) or im(s) appear
    irf = []
    for s in symbols:
        if s.is_real or s.is_imaginary:
            continue  # neither re(x) nor im(x) will appear
        # if re(s) or im(s) appear, the auxiliary equation must be present
        if any(fi.has(re(s), im(s)) for fi in f):
            irf.append((s, re(s) + S.ImaginaryUnit*im(s)))
    if irf:
        for s, rhs in irf:
            for i, fi in enumerate(f):
                f[i] = fi.xreplace({s: rhs})
            f.append(s - rhs)
            symbols.extend([re(s), im(s)])
        if bare_f:
            bare_f = False
        flags['dict'] = True
    # end of real/imag handling  -----------------------------

    symbols = list(uniq(symbols))
    if not ordered_symbols:
        # we do this to make the results returned canonical in case f
        # contains a system of nonlinear equations; all other cases should
        # be unambiguous
        symbols = sorted(symbols, key=default_sort_key)

    # we can solve for non-symbol entities by replacing them with Dummy symbols
    symbols_new = []
    symbol_swapped = False
    for i, s in enumerate(symbols):
        if s.is_Symbol:
            s_new = s
        else:
            symbol_swapped = True
            s_new = Dummy('X%d' % i)
        symbols_new.append(s_new)

    if symbol_swapped:
        swap_sym = list(zip(symbols, symbols_new))
        f = [fi.subs(swap_sym) for fi in f]
        symbols = symbols_new
        swap_sym = dict([(v, k) for k, v in swap_sym])
    else:
        swap_sym = {}

    # this is needed in the next two events
    symset = set(symbols)

    # get rid of equations that have no symbols of interest; we don't
    # try to solve them because the user didn't ask and they might be
    # hard to solve; this means that solutions may be given in terms
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
            if not isinstance(p, Expr) or isinstance(p, Piecewise):
                pass
            elif (isinstance(p, bool) or
                    not p.args or
                    p in symset or
                    p.is_Add or p.is_Mul or
                    p.is_Pow and not implicit or
                    p.is_Function and not implicit) and p.func not in (re, im):
                continue
            elif not p in seen:
                seen.add(p)
                if p.free_symbols & symset:
                    non_inverts.add(p)
                else:
                    continue
            pot.skip()
    del seen
    non_inverts = dict(list(zip(non_inverts, [Dummy() for d in non_inverts])))
    f = [fi.subs(non_inverts) for fi in f]

    non_inverts = [(v, k.subs(swap_sym)) for k, v in non_inverts.items()]

    # rationalize Floats
    floats = False
    if flags.get('rational', True) is not False:
        for i, fi in enumerate(f):
            if fi.has(Float):
                floats = True
                f[i] = nsimplify(fi, rational=True)

    # Any embedded piecewise functions need to be brought out to the
    # top level so that the appropriate strategy gets selected.
    # However, this is necessary only if one of the piecewise
    # functions depends on one of the symbols we are solving for.
    def _has_piecewise(e):
        if e.is_Piecewise:
            return e.has(*symbols)
        return any([_has_piecewise(a) for a in e.args])
    for i, fi in enumerate(f):
        if _has_piecewise(fi):
            f[i] = piecewise_fold(fi)

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
    # Restore masked-off objects
    if non_inverts:

        def _do_dict(solution):
            return dict([(k, v.subs(non_inverts)) for k, v in
                         solution.items()])
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

    # Restore original "symbols" if a dictionary is returned.
    # This is not necessary for
    #   - the single univariate equation case
    #     since the symbol will have been removed from the solution;
    #   - the nonlinear poly_system since that only supports zero-dimensional
    #     systems and those results come back as a list
    #
    # ** unless there were Derivatives with the symbols, but those were handled
    #    above.
    if symbol_swapped:
        symbols = [swap_sym[k] for k in symbols]
        if type(solution) is dict:
            solution = dict([(swap_sym[k], v.subs(swap_sym))
                             for k, v in solution.items()])
        elif solution and type(solution) is list and type(solution[0]) is dict:
            for i, sol in enumerate(solution):
                solution[i] = dict([(swap_sym[k], v.subs(swap_sym))
                              for k, v in sol.items()])

    # undo the dictionary solutions returned when the system was only partially
    # solved with poly-system if all symbols are present
    if (
            not flags.get('dict', False) and
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
        got_None = []  # solutions for which one or more symbols gave None
        no_False = []  # solutions for which no symbols gave False
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
                    for symb, val in sol.items():
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
            else:  # list of expressions
                for sol in solution:
                    test = check_assumptions(sol, **symbols[0].assumptions0)
                    if test is False:
                        continue
                    no_False.append(sol)
                    if test is None:
                        got_None.append(sol)

        elif type(solution) is dict:
            a_None = False
            for symb, val in solution.items():
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
            if len(symbols) != 1:
                raise ValueError("Length should be 1")
            if warning and symbols[0].assumptions0:
                warnings.warn(filldedent("""
                    \tWarning: assumptions about variable '%s' are
                    not handled currently.""" % symbols[0]))
            # TODO: check also variable assumptions for inequalities

        else:
            raise TypeError('Unrecognized solution')  # improve the checker

        solution = no_False
        if warning and got_None:
            warnings.warn(filldedent("""
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

    # return a list of mappings or []
    if not solution:
        solution = []
    else:
        if isinstance(solution, dict):
            solution = [solution]
        elif iterable(solution[0]):
            solution = [dict(list(zip(symbols, s))) for s in solution]
        elif isinstance(solution[0], dict):
            pass
        else:
            if len(symbols) != 1:
                raise ValueError("Length should be 1")
            solution = [{symbols[0]: s} for s in solution]
    if as_dict:
        return solution
    assert as_set
    if not solution:
        return [], set()
    k = sorted(list(solution[0].keys()), key=lambda i: i.sort_key())
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
        got_s = set([])
        result = []
        for s in symbols:
            n, d = solve_linear(f, symbols=[s])
            if n.is_Symbol:
                # no need to check but we should simplify if desired
                if flags.get('simplify', True):
                    d = simplify(d)
                if got_s and any([ss in d.free_symbols for ss in got_s]):
                    # sol depends on previously solved symbols: discard it
                    continue
                got_s.add(n)
                result.append({n: d})
            elif n and d:  # otherwise there was no solution for s
                failed.append(s)
        if not failed:
            return result
        for s in failed:
            try:
                soln = _solve(f, s, **flags)
                for sol in soln:
                    if got_s and any([ss in sol.free_symbols for ss in got_s]):
                        # sol depends on previously solved symbols: discard it
                        continue
                    got_s.add(s)
                    result.append({s: sol})
            except NotImplementedError:
                continue
        if got_s:
            return result
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
        for n, (expr, cond) in enumerate(f.args):
            candidates = _solve(expr, *symbols, **flags)

            for candidate in candidates:
                if candidate in result:
                    continue
                cond = (cond == True) or cond.subs(symbol, candidate)
                if cond != False:
                    # Only include solutions that do not match the condition
                    # of any previous pieces.
                    matches_other_piece = False
                    for other_n, (other_expr, other_cond) in enumerate(f.args):
                        if other_n == n:
                            break
                        if other_cond == False:
                            continue
                        if other_cond.subs(symbol, candidate) == True:
                            matches_other_piece = True
                            break
                    if not matches_other_piece:
                        result.add(Piecewise(
                            (candidate, cond == True or cond.doit()),
                            (S.NaN, True)
                        ))
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

        result = False  # no solution was obtained
        msg = ''  # there is no failure message
        dens = denoms(f, symbols)  # store these for checking later

        # Poly is generally robust enough to convert anything to
        # a polynomial and tell us the different generators that it
        # contains, so we will inspect the generators identified by
        # polys to figure out what to do.

        # but first remove radicals as this will help Polys
        if flags.pop('unrad', True):
            try:
                # try remove all...
                u = unrad(f_num)
            except ValueError:
                # ...else hope for the best while letting some remain
                try:
                    u = unrad(f, symbol)
                except ValueError:
                    u = None  # hope for best with original equation
            if u:
                flags['unrad'] = False  # don't unrad next time
                eq, cov, dens2 = u
                dens.update(dens2)
                if cov:
                    if len(cov) > 1:
                        raise NotImplementedError('Not sure how to handle this.')
                    isym, ieq = cov[0]
                    # since cov is written in terms of positive symbols, set
                    # check to False or else 0 would be excluded; the solution
                    # will be checked below
                    absent = Dummy()
                    check = flags.get('check', absent)
                    flags['check'] = False
                    sol = _solve(eq, isym, **flags)
                    inv = _solve(ieq, symbol, **flags)
                    result = []
                    for s in sol:
                        for i in inv:
                            result.append(i.subs(isym, s))
                    if check == absent:
                        flags.pop('check')
                    else:
                        flags['check'] = check
                else:
                    result = _solve(eq, symbol, **flags)

        if result is False:
            # rewrite hyperbolics in terms of exp
            f_num = f_num.replace(lambda w: isinstance(w, C.HyperbolicFunction),
                lambda w: w.rewrite(exp))

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
                    if c.is_Rational and c is not S.One:  # c could be a Float
                        return b**ee, c.q
                    return x, 1

                bases, qs = list(zip(*[_as_base_q(g) for g in gens]))
                bases = set(bases)

                if len(bases) > 1:
                    funcs = set(b for b in bases if b.is_Function)

                    trig = set([_ for _ in funcs if
                        isinstance(_, C.TrigonometricFunction)])
                    other = funcs - trig
                    if not other and len(funcs.intersection(trig)) > 1:
                        newf = TR1(f_num).rewrite(tan)
                        if newf != f_num:
                            return _solve(newf, symbol, **flags)

                    # just a simple case - see if replacement of single function
                    # clears all symbol-dependent functions, e.g.
                    # log(x) - log(log(x) - 1) - 3 can be solved even though it has
                    # two generators.

                    if funcs:
                        funcs = list(ordered(funcs))  # put shallowest function first
                        f1 = funcs[0]
                        t = Dummy()
                        # perform the substitution
                        ftry = f_num.subs(f1, t)

                        # if no Functions left, we can proceed with usual solve
                        if not ftry.has(symbol):
                            cv_sols = _solve(ftry, t, **flags)
                            cv_inv = _solve(t - f1, symbol, **flags)[0]
                            sols = list()
                            for sol in cv_sols:
                                sols.append(cv_inv.subs(t, sol))
                            return list(ordered(sols))

                    msg = 'multiple generators %s' % gens

                else:  # len(bases) == 1 and all(q == 1 for q in qs):
                    # e.g. case where gens are exp(x), exp(-x)
                    u = bases.pop()
                    t = Dummy('t')
                    inv = _solve(u - t, symbol, **flags)
                    if isinstance(u, (Pow, exp)):
                        # this will be resolved by factor in _tsolve but we might
                        # as well try a simple expansion here to get things in
                        # order so something like the following will work now without
                        # having to factor:
                        # >>> eq = (exp(I*(-x-2))+exp(I*(x+2)))
                        # >>> eq.subs(exp(x),y)  # fails
                        # exp(I*(-x - 2)) + exp(I*(x + 2))
                        # >>> eq.expand().subs(exp(x),y)  # works
                        # y**I*exp(2*I) + y**(-I)*exp(-2*I)
                        def _expand(p):
                            b, e = p.as_base_exp()
                            e = expand_mul(e)
                            return expand_power_exp(b**e)
                        ftry = f_num.replace(
                            lambda w: w.is_Pow or isinstance(w, exp),
                            _expand).subs(u, t)
                        if not ftry.has(symbol):
                            soln = _solve(ftry, t, **flags)
                            sols = list()
                            for sol in soln:
                                for i in inv:
                                    sols.append(i.subs(t, sol))
                            return list(ordered(sols))

            elif len(gens) == 1:

                # There is only one generator that we are interested in, but there
                # may have been more than one generator identified by polys (e.g.
                # for symbols other than the one we are interested in) so recast
                # the poly in terms of our generator of interest.

                if len(poly.gens) > 1:
                    poly = Poly(poly, gens[0])

                # if we aren't on the tsolve-pass, use roots
                if not flags.pop('tsolve', False):
                    flags['tsolve'] = True
                    if poly.degree() == 1 and (
                            poly.gen.is_Pow and
                            poly.gen.exp.is_Rational and
                            not poly.gen.exp.is_Integer):
                        pass
                    else:
                        # for cubics and quartics, if the flag wasn't set, DON'T do it
                        # by default since the results are quite long. Perhaps one
                        # could base this decision on a certain critical length of the
                        # roots.
                        deg = poly.degree()
                        if deg > 2:
                            flags['simplify'] = flags.get('simplify', False)

                        # TODO: Just pass composite=True to roots()
                        poly = Poly(poly.as_expr(), poly.gen, composite=True)
                        soln = list(roots(poly, cubics=True, quartics=True,
                                                             quintics=True).keys())

                        if len(soln) < deg:
                            try:
                                # get all_roots if possible
                                soln = list(ordered(uniq(poly.all_roots())))
                            except NotImplementedError:
                                pass
                        gen = poly.gen
                        if gen != symbol:
                            u = Dummy()
                            inversion = _solve(gen - u, symbol, **flags)
                            soln = list(ordered(set([i.subs(u, s) for i in
                                        inversion for s in soln])))
                        result = soln


    # fallback if above fails
    if result is False:

        # allow tsolve to be used on next pass if needed
        flags.pop('tsolve', None)
        try:
            result = _tsolve(f_num, symbol, **flags)
        except PolynomialError:
            result = None
        if result is None:
            result = False

    if result is False:
        raise NotImplementedError(msg +
        "\nNo algorithms are implemented to solve equation %s" % f)

    if flags.get('simplify', True):
        result = list(map(simplify, result))
        # we just simplified the solution so we now set the flag to
        # False so the simplification doesn't happen again in checksol()
        flags['simplify'] = False

    if check:
        # reject any result that makes any denom. affirmatively 0;
        # if in doubt, keep it
        result = [s for s in result if isinstance(s, RootOf) or
                  all(not checksol(den, {symbol: s}, **flags)
                    for den in dens)]
        # keep only results if the check is not False
        result = [r for r in result if isinstance(r, RootOf) or
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

        poly = g.as_poly(*symbols, extension=True)

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
                        j = monom.index(1)
                        matrix[i, j] = coeff
                    except ValueError:
                        matrix[i, m] = -coeff

            # returns a dictionary ({symbols: values}) or None
            if flags.pop('particular', False):
                result = minsolve_linear_system(matrix, *symbols, **flags)
            else:
                result = solve_linear_system(matrix, *symbols, **flags)
            if result:
                # it doesn't need to be checked but we need to see
                # that it didn't set any denominators to 0
                if any(checksol(d, result, **flags) for d in dens):
                    result = None
            if failed:
                if result:
                    solved_syms = list(result.keys())
                else:
                    solved_syms = []

        else:
            if len(symbols) > len(polys):
                from sympy.utilities.iterables import subsets

                free = set.union(*[p.free_symbols for p in polys])
                free = list(free.intersection(symbols))
                free.sort(key=default_sort_key)
                got_s = set([])
                result = []
                for syms in subsets(free, len(polys)):
                    try:
                        # returns [] or list of tuples of solutions for syms
                        res = solve_poly_system(polys, *syms)
                        if res:
                            for r in res:
                                skip = False
                                for r1 in r:
                                    if got_s and any([ss in r1.free_symbols
                                           for ss in got_s]):
                                        # sol depends on previously
                                        # solved symbols: discard it
                                        skip = True
                                if not skip:
                                    got_s.update(syms)
                                    result.extend([dict(list(zip(syms, r)))])
                    except NotImplementedError:
                        pass
                if got_s:
                    solved_syms = list(got_s)
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
                    #
                    # We do not check the solution obtained from polys, either.
                    result = [dict(list(zip(solved_syms, r))) for r in result]

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

        solved_syms = set(solved_syms)  # set of symbols we have solved for
        legal = set(symbols)  # what we are interested in
        simplify_flag = flags.get('simplify', None)
        do_simplify = flags.get('simplify', True)
        # sort so equation with the fewest potential symbols is first
        for eq in ordered(failed, lambda _: len(_ok_syms(_))):
            newresult = []
            bad_results = []
            got_s = set([])
            u = Dummy()
            for r in result:
                # update eq with everything that is known so far
                eq2 = eq.subs(r)
                # if check is True then we see if it satisfies this
                # equation, otherwise we just accept it
                if check and r:
                    b = checksol(u, u, eq2, minimal=True)
                    if b is not None:
                        # this solution is sufficient to know whether
                        # it is valid or not so we either accept or
                        # reject it, then continue
                        if b:
                            newresult.append(r)
                        else:
                            bad_results.append(r)
                        continue
                # search for a symbol amongst those available that
                # can be solved for
                ok_syms = _ok_syms(eq2, sort=True)
                if not ok_syms:
                    if r:
                        newresult.append(r)
                    break  # skip as it's independent of desired symbols
                for s in ok_syms:
                    try:
                        soln = _solve(eq2, s, **flags)
                    except NotImplementedError:
                        continue
                    # put each solution in r and append the now-expanded
                    # result in the new result list; use copy since the
                    # solution for s in being added in-place
                    if do_simplify:
                        flags['simplify'] = False  # for checksol's sake
                    for sol in soln:
                        if got_s and any([ss in sol.free_symbols for ss in got_s]):
                            # sol depends on previously solved symbols: discard it
                            continue
                        if check:
                            # check that it satisfies *other* equations
                            ok = False
                            for p in polys:
                                if checksol(p, s, sol, **flags) is False:
                                    break
                            else:
                                ok = True
                            if not ok:
                                continue
                            # check that it doesn't set any denominator to 0
                            if any(checksol(d, s, sol, **flags) for d in dens):
                                continue
                        # update existing solutions with this new one
                        rnew = r.copy()
                        for k, v in r.items():
                            rnew[k] = v.subs(s, sol)
                        # and add this new solution
                        rnew[s] = sol
                        newresult.append(rnew)
                    if simplify_flag is not None:
                        flags['simplify'] = simplify_flag
                    got_s.add(s)
                if not got_s:
                    raise NotImplementedError('could not solve %s' % eq2)
            if got_s:
                result = newresult
            for b in bad_results:
                result.remove(b)
        # if there is only one result should we return just the dictionary?
    return result


def solve_linear(lhs, rhs=0, symbols=[], exclude=[]):
    r""" Return a tuple derived from f = lhs - rhs that is either:

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
    dfree = d.free_symbols

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
                        if not d.has(xi) or not (d/xi).has(xi):
                            return xi, vi

        if all_zero:
            return S.Zero, S.One
    if n.is_Symbol:  # there was no valid solution
        n = d = S.Zero
    return n, d  # should we cancel now?


def minsolve_linear_system(system, *symbols, **flags):
    r"""
    Find a particular solution to a linear system.

    In particular, try to find a solution with the minimal possible number
    of non-zero variables. This is a very computationally hard prolem.
    If ``quick=True``, a heuristic is used. Otherwise a naive algorithm with
    exponential complexity is used.
    """
    quick = flags.get('quick', False)
    # Check if there are any non-zero solutions at all
    s0 = solve_linear_system(system, *symbols, **flags)
    if not s0 or all(v == 0 for v in s0.values()):
        return s0
    if quick:
        # We just solve the system and try to heuristically find a nice
        # solution.
        s = solve_linear_system(system, *symbols)
        def update(determined, solution):
            delete = []
            for k, v in solution.items():
                solution[k] = v.subs(determined)
                if not solution[k].free_symbols:
                    delete.append(k)
                    determined[k] = solution[k]
            for k in delete:
                del solution[k]
        determined = {}
        update(determined, s)
        while s:
            # NOTE sort by default_sort_key to get deterministic result
            k = max((k for k in s.values()),
                    key=lambda x: (len(x.free_symbols), default_sort_key(x)))
            x = max(k.free_symbols, key=default_sort_key)
            if len(k.free_symbols) != 1:
                determined[x] = S(0)
            else:
                val = solve(k)[0]
                if val == 0 and all(v.subs(x, val) == 0 for v in s.values()):
                    determined[x] = S(1)
                else:
                    determined[x] = val
            update(determined, s)
        return determined
    else:
        # We try to select n variables which we want to be non-zero.
        # All others will be assumed zero. We try to solve the modified system.
        # If there is a non-trivial solution, just set the free variables to
        # one. If we do this for increasing n, trying all combinations of
        # variables, we will find an optimal solution.
        # We speed up slightly by starting at one less than the number of
        # variables the quick method manages.
        from itertools import combinations
        from sympy.utilities.misc import debug
        N = len(symbols)
        bestsol = minsolve_linear_system(system, *symbols, quick=True)
        n0 = len([x for x in bestsol.values() if x != 0])
        for n in range(n0 - 1, 1, -1):
            debug('minsolve: %s' % n)
            thissol = None
            for nonzeros in combinations(list(range(N)), n):
                subm = Matrix([system.col(i).T for i in nonzeros] + [system.col(-1).T]).T
                s = solve_linear_system(subm, *[symbols[i] for i in nonzeros])
                if s and not all(v == 0 for v in s.values()):
                    subs = [(symbols[v], S(1)) for v in nonzeros]
                    for k, v in s.items():
                        s[k] = v.subs(subs)
                    for sym in symbols:
                        if sym not in s:
                            if symbols.index(sym) in nonzeros:
                                s[sym] = S(1)
                            else:
                                s[sym] = S(0)
                    thissol = s
                    break
            if thissol is None:
                break
            bestsol = thissol
        return bestsol


def solve_linear_system(system, *symbols, **flags):
    r"""
    Solve system of N linear equations with M variables, which means
    both under- and overdetermined systems are supported. The possible
    number of solutions is zero, one or infinite. Respectively, this
    procedure will return None or a dictionary with solutions. In the
    case of underdetermined systems, all arbitrary parameters are skipped.
    This may cause a situation in which an empty dictionary is returned.
    In that case, all symbols can be assigned arbitrary values.

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

    A degenerate system returns an empty dictionary.

    >>> system = Matrix(( (0,0,0), (0,0,0) ))
    >>> solve_linear_system(system, x, y)
    {}

    """
    do_simplify = flags.get('simplify', True)

    if system.rows == system.cols - 1 == len(symbols):
        try:
            # well behaved n-equations and n-unknowns
            inv = inv_quick(system[:, :-1])
            rv = dict(zip(symbols, inv*system[:, -1]))
            if do_simplify:
                for k, v in rv.items():
                    rv[k] = simplify(v)
            if not all(i.is_zero for i in rv.values()):
                # non-trivial solution
                return rv
        except ValueError:
            pass

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
                    # We need to know if this is always zero or not. We
                    # assume that if there are free symbols that it is not
                    # identically zero (or that there is more than one way
                    # to make this zero). Otherwise, if there are none, this
                    # is a constant and we assume that it does not simplify
                    # to zero XXX are there better (fast) ways to test this?
                    # The .equals(0) method could be used but that can be
                    # slow; numerical testing is prone to errors of scaling.
                    if not matrix[i, m].free_symbols:
                        return None  # no solution

                    # A row of zeros with a non-zero rhs can only be accepted
                    # if there is another equivalent row. Any such rows will
                    # be deleted.
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
                    # every choice of variable values is a solution
                    # so we return an empty dict instead of None
                    return dict()
                continue

            # we want to change the order of colums so
            # the order of variables must also change
            syms[i], syms[k] = syms[k], syms[i]
            matrix.col_swap(i, k)

        pivot_inv = S.One/matrix[i, i]

        # divide all elements in the current row by the pivot
        matrix.row_op(i, lambda x, _: x * pivot_inv)

        for k in xrange(i + 1, matrix.rows):
            if matrix[k, i]:
                coeff = matrix[k, i]

                # subtract from the current row the row containing
                # pivot and multiplied by extracted coefficient
                matrix.row_op(k, lambda x, j: simplify(x - matrix[i, j]*coeff))

        i += 1

    # if there weren't any problems, augmented matrix is now
    # in row-echelon form so we can check how many solutions
    # there are and extract them using back substitution

    if len(syms) == matrix.rows:
        # this system is Cramer equivalent so there is
        # exactly one solution to this system of equations
        k, solutions = i - 1, {}

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

    system = list(collect(equ.expand(), sym, evaluate=False).values())

    if not any(equ.has(sym) for equ in system):
        # consecutive powers in the input expressions have
        # been successfully collected, so solve remaining
        # system using Gaussian elimination algorithm
        return solve(system, *coeffs, **flags)
    else:
        return None  # no solutions


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
    if matrix.rows != matrix.cols - 1:
        raise ValueError("Rows should be equal to columns - 1")
    A = matrix[:matrix.rows, :matrix.rows]
    b = matrix[:, matrix.cols - 1:]
    soln = A.LUsolve(b)
    solutions = {}
    for i in range(soln.rows):
        solutions[syms[i]] = soln[i, 0]
    return solutions


def det_perm(M):
    """Return the det(``M``) by using permutations to select factors.
    For size larger than 8 the number of permutations becomes prohibitively
    large, or if there are no symbols in the matrix, it is better to use the
    standard determinant routines, e.g. `M.det()`.

    See Also
    ========
    det_minor
    det_quick
    """
    args = []
    s = True
    n = M.rows
    try:
        list = M._mat
    except AttributeError:
        list = flatten(M.tolist())
    for perm in generate_bell(n):
        fac = []
        idx = 0
        for j in perm:
            fac.append(list[idx + j])
            idx += n
        term = Mul(*fac) # disaster with unevaluated Mul -- takes forever for n=7
        args.append(term if s else -term)
        s = not s
    return Add(*args)


def det_minor(M):
    """Return the ``det(M)`` computed from minors without
    introducing new nesting in products.

    See Also
    ========
    det_perm
    det_quick
    """
    n = M.rows
    if n == 2:
        return M[0, 0]*M[1, 1] - M[1, 0]*M[0, 1]
    else:
        return sum([(1, -1)[i % 2]*Add(*[M[0, i]*d for d in
            Add.make_args(det_minor(M.minorMatrix(0, i)))])
            if M[0, i] else S.Zero for i in range(n)])


def det_quick(M, method=None):
    """Return ``det(M)`` assuming that either
    there are lots of zeros or the size of the matrix
    is small. If this assumption is not met, then the normal
    Matrix.det function will be used with method = ``method``.

    See Also
    ========
    det_minor
    det_perm
    """
    if any(i.has(Symbol) for i in M):
        if M.rows < 8 and all(i.has(Symbol) for i in M):
            return det_perm(M)
        return det_minor(M)
    else:
        return M.det(method=method) if method else M.det()


def inv_quick(M):
    """Return the inverse of ``M``, assuming that either
    there are lots of zeros or the size of the matrix
    is small.
    """
    from sympy.matrices import zeros
    if any(i.has(Symbol) for i in M):
        if all(i.has(Symbol) for i in M):
            det = lambda _: det_perm(_)
        else:
            det = lambda _: det_minor(_)
    else:
        return M.inv()
    n = M.rows
    d = det(M)
    if d is S.Zero:
        raise ValueError("Matrix det == 0; not invertible.")
    ret = zeros(n)
    s1 = -1
    for i in range(n):
        s = s1 = -s1
        for j in range(n):
            di = det(M.minorMatrix(i, j))
            ret[j, i] = s*di/d
            s = -s
    return ret


# these are functions that have multiple inverse values per period
multi_inverses = {
    sin: lambda x: (asin(x), S.Pi - asin(x)),
    cos: lambda x: (acos(x), 2*S.Pi - acos(x)),
}


def _tsolve(eq, sym, **flags):
    """
    Helper for _solve that solves a transcendental equation with respect
    to the given symbol. Various equations containing powers and logarithms,
    can be solved.

    There is currently no guarantee that all solutions will be returned or
    that a real solution will be favored over a complex one.

    Examples
    ========

    >>> from sympy import log
    >>> from sympy.solvers.solvers import _tsolve as tsolve
    >>> from sympy.abc import x

    >>> tsolve(3**(2*x + 5) - 4, x)
    [-5/2 + log(2)/log(3), log(-2*sqrt(3)/27)/log(3)]

    >>> tsolve(log(x) + 2*x, x)
    [LambertW(2)/2]

    """
    if 'tsolve_saw' not in flags:
        flags['tsolve_saw'] = []
    if eq in flags['tsolve_saw']:
        return None
    else:
        flags['tsolve_saw'].append(eq)

    rhs, lhs = _invert(eq, sym)

    if lhs == sym:
        return [rhs]
    try:
        if lhs.is_Add:
            # it's time to try factoring; powdenest is used
            # to try get powers in standard form for better factoring
            f = factor(powdenest(lhs - rhs))
            if f.is_Mul:
                return _solve(f, sym, **flags)
            if rhs:
                f = logcombine(lhs, force=flags.get('force', True))
                if f.count(log) != lhs.count(log):
                    if f.func is log:
                        return _solve(f.args[0] - exp(rhs), sym, **flags)
                    return _tsolve(f - rhs, sym)

        elif lhs.is_Pow:
            if lhs.exp.is_Integer:
                if lhs - rhs != eq:
                    return _solve(lhs - rhs, sym, **flags)
            elif sym not in lhs.exp.free_symbols:
                return _solve(lhs.base - rhs**(1/lhs.exp), sym, **flags)
            elif not rhs and sym in lhs.exp.free_symbols:
                # f(x)**g(x) only has solutions where f(x) == 0 and g(x) != 0 at
                # the same place
                sol_base = _solve(lhs.base, sym, **flags)
                if not sol_base:
                    return sol_base  # no solutions to remove so return now
                return list(ordered(set(sol_base) - set(
                    _solve(lhs.exp, sym, **flags))))
            elif (rhs is not S.Zero and
                        lhs.base.is_positive and
                        lhs.exp.is_real):
                return _solve(lhs.exp*log(lhs.base) - log(rhs), sym, **flags)

        elif lhs.is_Mul and rhs.is_positive:
            llhs = expand_log(log(lhs))
            if llhs.is_Add:
                return _solve(llhs - log(rhs), sym, **flags)

        elif lhs.is_Function and len(lhs.args) == 1 and lhs.func in multi_inverses:
            # sin(x) = 1/3 -> x - asin(1/3) & x - (pi - asin(1/3))
            soln = []
            for i in multi_inverses[lhs.func](rhs):
                soln.extend(_solve(lhs.args[0] - i, sym, **flags))
            return list(ordered(soln))

        rewrite = lhs.rewrite(exp)
        if rewrite != lhs:
            return _solve(rewrite - rhs, sym, **flags)
    except NotImplementedError:
        pass

    # maybe it is a lambert pattern
    if flags.pop('bivariate', True):
        # lambert forms may need some help being recognized, e.g. changing
        # 2**(3*x) + x**3*log(2)**3 + 3*x**2*log(2)**2 + 3*x*log(2) + 1
        # to 2**(3*x) + (x*log(2) + 1)**3
        g = _filtered_gens(eq.as_poly(), sym)
        up_or_log = set()
        for gi in g:
            if gi.func is exp or gi.func is log:
                up_or_log.add(gi)
            elif gi.is_Pow:
                gisimp = powdenest(expand_power_exp(gi))
                if gisimp.is_Pow and sym in gisimp.exp.free_symbols:
                    up_or_log.add(gi)
        down = g.difference(up_or_log)
        eq_down = expand_log(expand_power_exp(eq)).subs(
            dict(list(zip(up_or_log, [0]*len(up_or_log)))))
        eq = expand_power_exp(factor(eq_down, deep=True) + (eq - eq_down))
        rhs, lhs = _invert(eq, sym)
        if lhs.has(sym):
            try:
                poly = lhs.as_poly()
                g = _filtered_gens(poly, sym)
                return _solve_lambert(lhs - rhs, sym, g)
            except NotImplementedError:
                # maybe it's a convoluted function
                if len(g) == 2:
                    try:
                        gpu = bivariate_type(lhs - rhs, *g)
                        if gpu is None:
                            raise NotImplementedError
                        g, p, u = gpu
                        flags['bivariate'] = False
                        inversion = _tsolve(g - u, sym, **flags)
                        if inversion:
                            sol = _solve(p, u, **flags)
                            return list(ordered(set([i.subs(u, s)
                                for i in inversion for s in sol])))
                    except NotImplementedError:
                        pass

    if flags.pop('force', True):
        flags['force'] = False
        pos, reps = posify(lhs - rhs)
        for u, s in reps.items():
            if s == sym:
                break
        else:
            u = sym
        try:
            soln = _solve(pos, u, **flags)
        except NotImplementedError:
            return
        return list(ordered([s.subs(reps) for s in soln]))


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
    >>> print(nsolve((f1, f2), (x1, x2), (-1, 1)))
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
    if iterable(f):
        f = list(f)
        for i, fi in enumerate(f):
            if isinstance(fi, Equality):
                f[i] = fi.lhs - fi.rhs
        f = Matrix(f).T
    if not isinstance(f, Matrix):
        # assume it's a sympy expression
        if isinstance(f, Equality):
            f = f.lhs - f.rhs
        f = f.evalf()
        syms = f.free_symbols
        if fargs is None:
            fargs = syms.copy().pop()
        if not (len(syms) == 1 and (fargs in syms or fargs[0] in syms)):
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
        print('f(x):')
        print(f)
    # derive Jacobian
    J = f.jacobian(fargs)
    if verbose:
        print('J(x):')
        print(J)
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
    solutions to the expression ``i - d`` will be related to the solutions of
    the original expression.

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
    (1/2, cos(x))
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
                for d, i in terms.items():
                    if len(i) > 1:
                        args.append(Add(*i)*d)
                    else:
                        args.append(i[0]*d)
                lhs = Add(*args)

        # if it's a two-term Add with rhs = 0 and two powers we can get the
        # dependent terms together, e.g. 3*f(x) + 2*g(x) -> f(x)/g(x) = -2/3
        if lhs.is_Add and not rhs and len(lhs.args) == 2 and \
                not lhs.is_polynomial(*symbols):
            a, b = ordered(lhs.args)
            ai, ad = a.as_independent(*symbols)
            bi, bd = b.as_independent(*symbols)
            if any(_ispow(i) for i in (ad, bd)):
                a_base, a_exp = ad.as_base_exp()
                b_base, b_exp = bd.as_base_exp()
                if a_base == b_base:
                    # a = -b
                    lhs = powsimp(powdenest(ad/bd))
                    rhs = -bi/ai
                else:
                    rat = ad/bd
                    _lhs = powsimp(ad/bd)
                    if _lhs != rat:
                        lhs = _lhs
                        rhs = -bi/ai
            if ai*bi is S.NegativeOne:
                if all(
                        isinstance(i, Function) for i in (ad, bd)) and \
                        ad.func == bd.func and len(ad.args) == len(bd.args):
                    if len(ad.args) == 1:
                        lhs = ad.args[0] - bd.args[0]
                    else:
                        # should be able to solve
                        # f(x, y) == f(2, 3) -> x == 2
                        # f(x, x + y) == f(2, 3) -> x == 2 or x == 3 - y
                        raise NotImplementedError('equal function with more than 1 argument')

        elif lhs.is_Mul and any(_ispow(a) for a in lhs.args):
            lhs = powsimp(powdenest(lhs))

        if lhs.is_Function:
            if hasattr(lhs, 'inverse') and len(lhs.args) == 1:
                #                    -1
                # f(x) = g  ->  x = f  (g)
                #
                # /!\ inverse should not be defined if there are multiple values
                # for the function -- these are handled in _tsolve
                #
                rhs = lhs.inverse()(rhs)
                lhs = lhs.args[0]
            elif lhs.func is atan2:
                y, x = lhs.args
                lhs = 2*atan(y/(sqrt(x**2 + y**2) + x))
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

    Otherwise the tuple, ``(eq, cov, dens)``, is returned where::

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
    Two options are also recognized::

        ``take``, when defined, is interpreted as a single-argument function
        that returns True if a given Pow should be handled.
        ``all``, when True, will signify that an attempt should be made to
        remove all radicals. ``take``, if present, has priority over ``all``.

    Radicals can be removed from an expression if::

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
        >>> unrad(sqrt(x)*x**Rational(1, 3) + 2)
        (x**5 - 64, [], [])
        >>> unrad(sqrt(x) + (x + 1)**Rational(1, 3))
        (x**3 - x**2 - 2*x - 1, [], [])
        >>> unrad(sqrt(x) + x**Rational(1, 3) + 2)
        (_p**3 + _p**2 + 2, [(_p, -_p**6 + x)], [])

    """
    def _canonical(eq):
        # remove constants since these don't change the location of the root
        # and expand the expression
        eq = factor_terms(eq)
        if eq.is_Mul:
            eq = Mul(*[f for f in eq.args if not f.is_number])
        eq = _mexpand(eq)

        # make sign canonical
        free = eq.free_symbols
        if len(free) == 1:
            if (eq.coeff(free.pop()**degree(eq)) < 0) == True:
                eq = -eq
        elif eq.could_extract_minus_sign():
            eq = -eq

        return eq

    if eq.is_Atom:
        return
    cov, dens, nwas, rpt = [flags.get(k, v) for k, v in
        sorted(dict(dens=None, cov=None, n=None, rpt=0).items())]

    if flags.get('take', None):
        _take = flags.pop('take')
    elif flags.pop('all', None):
        _rad = lambda w: w.is_Pow and w.exp.is_Rational and w.exp.q != 1
        def _take(d):
            return _rad(d) or any(_rad(i) for i in d.atoms(Pow))
        if eq.has(S.ImaginaryUnit):
            i = Dummy()
            flags['take'] = _take
            try:
                rv = unrad(eq.xreplace({S.ImaginaryUnit: sqrt(i)}), *syms, **flags)
                rep = {i: S.NegativeOne}
                rv = (_canonical(rv[0].xreplace(rep)),
                      [tuple([j.xreplace(rep) for j in i]) for i in rv[1]],
                      [i.xreplace(rep) for i in rv[2]])
                return rv
            except ValueError as msg:
                raise msg
    else:
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

    eq = powdenest(factor_terms(eq, radical=True))
    eq, d = eq.as_numer_denom()
    eq = _mexpand(eq)
    if _take(d):
        dens.add(d)

    if not eq.free_symbols:
        return eq, cov, list(dens)

    poly = eq.as_poly()

    # if all the bases are the same or all the radicals are in one
    # term, `lcm` will be the lcm of the radical's exponent
    # denominators
    lcm = 1
    rads = set()
    bases = set()
    for g in poly.gens:
        if not _take(g) or not g.is_Pow:
            continue
        ecoeff = g.exp.as_coeff_mul()[0]  # a Rational
        if ecoeff.q != 1:
            rads.add(g)
            lcm = ilcm(lcm, ecoeff.q)
            bases.add(g.base)

    if not rads:
        return

    depth = sqrt_depth(eq)

    # get terms together that have common generators
    drad = dict(list(zip(rads, list(range(len(rads))))))
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
    rterms = list(reversed(list(ordered(rterms))))

    # continue handling
    ok = True
    if len(rterms) == 1:
        eq = rterms[0]**lcm - (-args)**lcm

    elif len(rterms) == 2 and not args:
        eq = rterms[0]**lcm - rterms[1]**lcm

    elif log(lcm, 2).is_Integer and (not args and
            len(rterms) == 4 or len(rterms) < 4):
        def _norm2(a, b):
            return a**2 + b**2 + 2*a*b

        if len(rterms) == 4:
            # (r0+r1)**2 - (r2+r3)**2
            r0, r1, r2, r3 = rterms
            eq = _norm2(r0, r1) - _norm2(r2, r3)
        elif len(rterms) == 3:
            # (r1+r2)**2 - (r0+args)**2
            r0, r1, r2 = rterms
            eq = _norm2(r1, r2) - _norm2(r0, args)
        elif len(rterms) == 2:
            # r0**2 - (r1+args)**2
            r0, r1 = rterms
            eq = r0**2 - _norm2(r1, args)

    elif len(bases) == 1:  # change of variables may work
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
    rpt += 1  # XXX how many repeats with others unchanging is enough?
    if not ok or (
                nwas is not None and len(rterms) == nwas and
                new_depth is not None and new_depth == depth and
                rpt > 3):
        # XXX: XFAIL tests indicate other cases that should be handled.
        raise ValueError('Cannot remove all radicals from %s' % eq)

    neq = unrad(eq, *syms, cov=cov, dens=dens, n=len(rterms), rpt=rpt, take=_take)
    if neq:
        eq = neq[0]

    return (_canonical(eq), cov, list(dens))


from sympy.solvers.bivariate import (
    bivariate_type, _solve_lambert, _filtered_gens)
