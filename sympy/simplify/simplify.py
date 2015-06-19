from __future__ import print_function, division

from collections import defaultdict

from sympy.core import (Basic, S, Add, Mul, Pow,
    Symbol, sympify, expand, expand_mul, expand_func,
    Function, Dummy, Expr, factor_terms,
    FunctionClass, symbols,
    expand_power_exp)
from sympy.core.add import _unevaluated_Add
from sympy.core.compatibility import (combinations_with_replacement, iterable,
    reduce, default_sort_key,
    ordered, range, as_int)
from sympy.core.numbers import Float, I, pi, Rational, Integer
from sympy.core.function import expand_log, count_ops, _mexpand
from sympy.core.mul import prod, _unevaluated_Mul
from sympy.core.rules import Transform
from sympy.core.evaluate import global_evaluate
from sympy.functions import (
    gamma, exp, sqrt, log, root, exp_polar,
    sin, cos, tan, cot, sinh, cosh, tanh, coth, piecewise_fold)
from sympy.functions.elementary.complexes import polar_lift, principal_branch
from sympy.functions.combinatorial.factorials import binomial, CombinatorialFunction, factorial
from sympy.functions.elementary.exponential import ExpBase
from sympy.functions.elementary.hyperbolic import HyperbolicFunction
from sympy.functions.elementary.integers import ceiling
from sympy.functions.elementary.trigonometric import TrigonometricFunction

from sympy.utilities.iterables import has_variety, sift, uniq

from sympy.utilities.misc import debug
from sympy.utilities.timeutils import timethis


from sympy.simplify.radsimp import radsimp, collect_sqrt, fraction
from sympy.simplify.trigsimp import trigsimp, exptrigsimp
from sympy.simplify.besselsimp import besselsimp
from sympy.simplify.powsimp import powsimp
from sympy.simplify.cse_opts import sub_pre, sub_post
from sympy.simplify.sqrtdenest import sqrtdenest

from sympy.polys import (Poly, together, reduced, cancel, factor,
    ComputationFailed, lcm, gcd, parallel_poly_from_expr)
from sympy.polys.monomials import Monomial, monomial_div
from sympy.polys.polyerrors import PolificationFailed, DomainError


import mpmath



def separatevars(expr, symbols=[], dict=False, force=False):
    """
    Separates variables in an expression, if possible.  By
    default, it separates with respect to all symbols in an
    expression and collects constant coefficients that are
    independent of symbols.

    If dict=True then the separated terms will be returned
    in a dictionary keyed to their corresponding symbols.
    By default, all symbols in the expression will appear as
    keys; if symbols are provided, then all those symbols will
    be used as keys, and any terms in the expression containing
    other symbols or non-symbols will be returned keyed to the
    string 'coeff'. (Passing None for symbols will return the
    expression in a dictionary keyed to 'coeff'.)

    If force=True, then bases of powers will be separated regardless
    of assumptions on the symbols involved.

    Notes
    =====
    The order of the factors is determined by Mul, so that the
    separated expressions may not necessarily be grouped together.

    Although factoring is necessary to separate variables in some
    expressions, it is not necessary in all cases, so one should not
    count on the returned factors being factored.

    Examples
    ========

    >>> from sympy.abc import x, y, z, alpha
    >>> from sympy import separatevars, sin
    >>> separatevars((x*y)**y)
    (x*y)**y
    >>> separatevars((x*y)**y, force=True)
    x**y*y**y

    >>> e = 2*x**2*z*sin(y)+2*z*x**2
    >>> separatevars(e)
    2*x**2*z*(sin(y) + 1)
    >>> separatevars(e, symbols=(x, y), dict=True)
    {'coeff': 2*z, x: x**2, y: sin(y) + 1}
    >>> separatevars(e, [x, y, alpha], dict=True)
    {'coeff': 2*z, alpha: 1, x: x**2, y: sin(y) + 1}

    If the expression is not really separable, or is only partially
    separable, separatevars will do the best it can to separate it
    by using factoring.

    >>> separatevars(x + x*y - 3*x**2)
    -x*(3*x - y - 1)

    If the expression is not separable then expr is returned unchanged
    or (if dict=True) then None is returned.

    >>> eq = 2*x + y*sin(x)
    >>> separatevars(eq) == eq
    True
    >>> separatevars(2*x + y*sin(x), symbols=(x, y), dict=True) == None
    True

    """
    expr = sympify(expr)
    if dict:
        return _separatevars_dict(_separatevars(expr, force), symbols)
    else:
        return _separatevars(expr, force)


def _separatevars(expr, force):
    if len(expr.free_symbols) == 1:
        return expr
    # don't destroy a Mul since much of the work may already be done
    if expr.is_Mul:
        args = list(expr.args)
        changed = False
        for i, a in enumerate(args):
            args[i] = separatevars(a, force)
            changed = changed or args[i] != a
        if changed:
            expr = expr.func(*args)
        return expr

    # get a Pow ready for expansion
    if expr.is_Pow:
        expr = Pow(separatevars(expr.base, force=force), expr.exp)

    # First try other expansion methods
    expr = expr.expand(mul=False, multinomial=False, force=force)

    _expr, reps = posify(expr) if force else (expr, {})
    expr = factor(_expr).subs(reps)

    if not expr.is_Add:
        return expr

    # Find any common coefficients to pull out
    args = list(expr.args)
    commonc = args[0].args_cnc(cset=True, warn=False)[0]
    for i in args[1:]:
        commonc &= i.args_cnc(cset=True, warn=False)[0]
    commonc = Mul(*commonc)
    commonc = commonc.as_coeff_Mul()[1]  # ignore constants
    commonc_set = commonc.args_cnc(cset=True, warn=False)[0]

    # remove them
    for i, a in enumerate(args):
        c, nc = a.args_cnc(cset=True, warn=False)
        c = c - commonc_set
        args[i] = Mul(*c)*Mul(*nc)
    nonsepar = Add(*args)

    if len(nonsepar.free_symbols) > 1:
        _expr = nonsepar
        _expr, reps = posify(_expr) if force else (_expr, {})
        _expr = (factor(_expr)).subs(reps)

        if not _expr.is_Add:
            nonsepar = _expr

    return commonc*nonsepar


def _separatevars_dict(expr, symbols):
    if symbols:
        if not all((t.is_Atom for t in symbols)):
            raise ValueError("symbols must be Atoms.")
        symbols = list(symbols)
    elif symbols is None:
        return {'coeff': expr}
    else:
        symbols = list(expr.free_symbols)
        if not symbols:
            return None

    ret = dict(((i, []) for i in symbols + ['coeff']))

    for i in Mul.make_args(expr):
        expsym = i.free_symbols
        intersection = set(symbols).intersection(expsym)
        if len(intersection) > 1:
            return None
        if len(intersection) == 0:
            # There are no symbols, so it is part of the coefficient
            ret['coeff'].append(i)
        else:
            ret[intersection.pop()].append(i)

    # rebuild
    for k, v in ret.items():
        ret[k] = Mul(*v)

    return ret


def ratsimp(expr):
    """
    Put an expression over a common denominator, cancel and reduce.

    Examples
    ========

    >>> from sympy import ratsimp
    >>> from sympy.abc import x, y
    >>> ratsimp(1/x + 1/y)
    (x + y)/(x*y)
    """

    f, g = cancel(expr).as_numer_denom()
    try:
        Q, r = reduced(f, [g], field=True, expand=False)
    except ComputationFailed:
        return f/g

    return Add(*Q) + cancel(r/g)


def ratsimpmodprime(expr, G, *gens, **args):
    """
    Simplifies a rational expression ``expr`` modulo the prime ideal
    generated by ``G``.  ``G`` should be a Groebner basis of the
    ideal.

    >>> from sympy.simplify.simplify import ratsimpmodprime
    >>> from sympy.abc import x, y
    >>> eq = (x + y**5 + y)/(x - y)
    >>> ratsimpmodprime(eq, [x*y**5 - x - y], x, y, order='lex')
    (x**2 + x*y + x + y)/(x**2 - x*y)

    If ``polynomial`` is False, the algorithm computes a rational
    simplification which minimizes the sum of the total degrees of
    the numerator and the denominator.

    If ``polynomial`` is True, this function just brings numerator and
    denominator into a canonical form. This is much faster, but has
    potentially worse results.

    References
    ==========

    M. Monagan, R. Pearce, Rational Simplification Modulo a Polynomial
    Ideal,
    http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.163.6984
    (specifically, the second algorithm)
    """
    from sympy import solve

    quick = args.pop('quick', True)
    polynomial = args.pop('polynomial', False)
    debug('ratsimpmodprime', expr)

    # usual preparation of polynomials:

    num, denom = cancel(expr).as_numer_denom()

    try:
        polys, opt = parallel_poly_from_expr([num, denom] + G, *gens, **args)
    except PolificationFailed:
        return expr

    domain = opt.domain

    if domain.has_assoc_Field:
        opt.domain = domain.get_field()
    else:
        raise DomainError(
            "can't compute rational simplification over %s" % domain)

    # compute only once
    leading_monomials = [g.LM(opt.order) for g in polys[2:]]
    tested = set()

    def staircase(n):
        """
        Compute all monomials with degree less than ``n`` that are
        not divisible by any element of ``leading_monomials``.
        """
        if n == 0:
            return [1]
        S = []
        for mi in combinations_with_replacement(range(len(opt.gens)), n):
            m = [0]*len(opt.gens)
            for i in mi:
                m[i] += 1
            if all([monomial_div(m, lmg) is None for lmg in
                    leading_monomials]):
                S.append(m)

        return [Monomial(s).as_expr(*opt.gens) for s in S] + staircase(n - 1)

    def _ratsimpmodprime(a, b, allsol, N=0, D=0):
        """
        Computes a rational simplification of ``a/b`` which minimizes
        the sum of the total degrees of the numerator and the denominator.

        The algorithm proceeds by looking at ``a * d - b * c`` modulo
        the ideal generated by ``G`` for some ``c`` and ``d`` with degree
        less than ``a`` and ``b`` respectively.
        The coefficients of ``c`` and ``d`` are indeterminates and thus
        the coefficients of the normalform of ``a * d - b * c`` are
        linear polynomials in these indeterminates.
        If these linear polynomials, considered as system of
        equations, have a nontrivial solution, then `\frac{a}{b}
        \equiv \frac{c}{d}` modulo the ideal generated by ``G``. So,
        by construction, the degree of ``c`` and ``d`` is less than
        the degree of ``a`` and ``b``, so a simpler representation
        has been found.
        After a simpler representation has been found, the algorithm
        tries to reduce the degree of the numerator and denominator
        and returns the result afterwards.

        As an extension, if quick=False, we look at all possible degrees such
        that the total degree is less than *or equal to* the best current
        solution. We retain a list of all solutions of minimal degree, and try
        to find the best one at the end.
        """
        c, d = a, b
        steps = 0

        maxdeg = a.total_degree() + b.total_degree()
        if quick:
            bound = maxdeg - 1
        else:
            bound = maxdeg
        while N + D <= bound:
            if (N, D) in tested:
                break
            tested.add((N, D))

            M1 = staircase(N)
            M2 = staircase(D)
            debug('%s / %s: %s, %s' % (N, D, M1, M2))

            Cs = symbols("c:%d" % len(M1), cls=Dummy)
            Ds = symbols("d:%d" % len(M2), cls=Dummy)
            ng = Cs + Ds

            c_hat = Poly(
                sum([Cs[i] * M1[i] for i in range(len(M1))]), opt.gens + ng)
            d_hat = Poly(
                sum([Ds[i] * M2[i] for i in range(len(M2))]), opt.gens + ng)

            r = reduced(a * d_hat - b * c_hat, G, opt.gens + ng,
                        order=opt.order, polys=True)[1]

            S = Poly(r, gens=opt.gens).coeffs()
            sol = solve(S, Cs + Ds, particular=True, quick=True)

            if sol and not all([s == 0 for s in sol.values()]):
                c = c_hat.subs(sol)
                d = d_hat.subs(sol)

                # The "free" variables occuring before as parameters
                # might still be in the substituted c, d, so set them
                # to the value chosen before:
                c = c.subs(dict(list(zip(Cs + Ds, [1] * (len(Cs) + len(Ds))))))
                d = d.subs(dict(list(zip(Cs + Ds, [1] * (len(Cs) + len(Ds))))))

                c = Poly(c, opt.gens)
                d = Poly(d, opt.gens)
                if d == 0:
                    raise ValueError('Ideal not prime?')

                allsol.append((c_hat, d_hat, S, Cs + Ds))
                if N + D != maxdeg:
                    allsol = [allsol[-1]]

                break

            steps += 1
            N += 1
            D += 1

        if steps > 0:
            c, d, allsol = _ratsimpmodprime(c, d, allsol, N, D - steps)
            c, d, allsol = _ratsimpmodprime(c, d, allsol, N - steps, D)

        return c, d, allsol

    # preprocessing. this improves performance a bit when deg(num)
    # and deg(denom) are large:
    num = reduced(num, G, opt.gens, order=opt.order)[1]
    denom = reduced(denom, G, opt.gens, order=opt.order)[1]

    if polynomial:
        return (num/denom).cancel()

    c, d, allsol = _ratsimpmodprime(
        Poly(num, opt.gens), Poly(denom, opt.gens), [])
    if not quick and allsol:
        debug('Looking for best minimal solution. Got: %s' % len(allsol))
        newsol = []
        for c_hat, d_hat, S, ng in allsol:
            sol = solve(S, ng, particular=True, quick=False)
            newsol.append((c_hat.subs(sol), d_hat.subs(sol)))
        c, d = min(newsol, key=lambda x: len(x[0].terms()) + len(x[1].terms()))

    if not domain.has_Field:
        cn, c = c.clear_denoms(convert=True)
        dn, d = d.clear_denoms(convert=True)
    r = Rational(cn, dn)

    return (c*r.q)/(d*r.p)


def _polarify(eq, lift, pause=False):
    from sympy import Integral
    if eq.is_polar:
        return eq
    if eq.is_number and not pause:
        return polar_lift(eq)
    if isinstance(eq, Symbol) and not pause and lift:
        return polar_lift(eq)
    elif eq.is_Atom:
        return eq
    elif eq.is_Add:
        r = eq.func(*[_polarify(arg, lift, pause=True) for arg in eq.args])
        if lift:
            return polar_lift(r)
        return r
    elif eq.is_Function:
        return eq.func(*[_polarify(arg, lift, pause=False) for arg in eq.args])
    elif isinstance(eq, Integral):
        # Don't lift the integration variable
        func = _polarify(eq.function, lift, pause=pause)
        limits = []
        for limit in eq.args[1:]:
            var = _polarify(limit[0], lift=False, pause=pause)
            rest = _polarify(limit[1:], lift=lift, pause=pause)
            limits.append((var,) + rest)
        return Integral(*((func,) + tuple(limits)))
    else:
        return eq.func(*[_polarify(arg, lift, pause=pause)
                         if isinstance(arg, Expr) else arg for arg in eq.args])


def polarify(eq, subs=True, lift=False):
    """
    Turn all numbers in eq into their polar equivalents (under the standard
    choice of argument).

    Note that no attempt is made to guess a formal convention of adding
    polar numbers, expressions like 1 + x will generally not be altered.

    Note also that this function does not promote exp(x) to exp_polar(x).

    If ``subs`` is True, all symbols which are not already polar will be
    substituted for polar dummies; in this case the function behaves much
    like posify.

    If ``lift`` is True, both addition statements and non-polar symbols are
    changed to their polar_lift()ed versions.
    Note that lift=True implies subs=False.

    >>> from sympy import polarify, sin, I
    >>> from sympy.abc import x, y
    >>> expr = (-x)**y
    >>> expr.expand()
    (-x)**y
    >>> polarify(expr)
    ((_x*exp_polar(I*pi))**_y, {_x: x, _y: y})
    >>> polarify(expr)[0].expand()
    _x**_y*exp_polar(_y*I*pi)
    >>> polarify(x, lift=True)
    polar_lift(x)
    >>> polarify(x*(1+y), lift=True)
    polar_lift(x)*polar_lift(y + 1)

    Adds are treated carefully:

    >>> polarify(1 + sin((1 + I)*x))
    (sin(_x*polar_lift(1 + I)) + 1, {_x: x})
    """
    if lift:
        subs = False
    eq = _polarify(sympify(eq), lift)
    if not subs:
        return eq
    reps = dict([(s, Dummy(s.name, polar=True)) for s in eq.free_symbols])
    eq = eq.subs(reps)
    return eq, dict([(r, s) for s, r in reps.items()])


def _unpolarify(eq, exponents_only, pause=False):
    if isinstance(eq, bool) or eq.is_Atom:
        return eq

    if not pause:
        if eq.func is exp_polar:
            return exp(_unpolarify(eq.exp, exponents_only))
        if eq.func is principal_branch and eq.args[1] == 2*pi:
            return _unpolarify(eq.args[0], exponents_only)
        if (
            eq.is_Add or eq.is_Mul or eq.is_Boolean or
            eq.is_Relational and (
                eq.rel_op in ('==', '!=') and 0 in eq.args or
                eq.rel_op not in ('==', '!='))
        ):
            return eq.func(*[_unpolarify(x, exponents_only) for x in eq.args])
        if eq.func is polar_lift:
            return _unpolarify(eq.args[0], exponents_only)

    if eq.is_Pow:
        expo = _unpolarify(eq.exp, exponents_only)
        base = _unpolarify(eq.base, exponents_only,
            not (expo.is_integer and not pause))
        return base**expo

    if eq.is_Function and getattr(eq.func, 'unbranched', False):
        return eq.func(*[_unpolarify(x, exponents_only, exponents_only)
            for x in eq.args])

    return eq.func(*[_unpolarify(x, exponents_only, True) for x in eq.args])


def unpolarify(eq, subs={}, exponents_only=False):
    """
    If p denotes the projection from the Riemann surface of the logarithm to
    the complex line, return a simplified version eq' of `eq` such that
    p(eq') == p(eq).
    Also apply the substitution subs in the end. (This is a convenience, since
    ``unpolarify``, in a certain sense, undoes polarify.)

    >>> from sympy import unpolarify, polar_lift, sin, I
    >>> unpolarify(polar_lift(I + 2))
    2 + I
    >>> unpolarify(sin(polar_lift(I + 7)))
    sin(7 + I)
    """
    if isinstance(eq, bool):
        return eq

    eq = sympify(eq)
    if subs != {}:
        return unpolarify(eq.subs(subs))
    changed = True
    pause = False
    if exponents_only:
        pause = True
    while changed:
        changed = False
        res = _unpolarify(eq, exponents_only, pause)
        if res != eq:
            changed = True
            eq = res
        if isinstance(res, bool):
            return res
    # Finally, replacing Exp(0) by 1 is always correct.
    # So is polar_lift(0) -> 0.
    return res.subs({exp_polar(0): 1, polar_lift(0): 0})


def _is_sum_surds(p):
    args = p.args if p.is_Add else [p]
    for y in args:
        if not ((y**2).is_Rational and y.is_real):
            return False
    return True


def posify(eq):
    """Return eq (with generic symbols made positive) and a restore
    dictionary.

    Any symbol that has positive=None will be replaced with a positive dummy
    symbol having the same name. This replacement will allow more symbolic
    processing of expressions, especially those involving powers and
    logarithms.

    A dictionary that can be sent to subs to restore eq to its original
    symbols is also returned.

    >>> from sympy import posify, Symbol, log
    >>> from sympy.abc import x
    >>> posify(x + Symbol('p', positive=True) + Symbol('n', negative=True))
    (_x + n + p, {_x: x})

    >> log(1/x).expand() # should be log(1/x) but it comes back as -log(x)
    log(1/x)

    >>> log(posify(1/x)[0]).expand() # take [0] and ignore replacements
    -log(_x)
    >>> eq, rep = posify(1/x)
    >>> log(eq).expand().subs(rep)
    -log(x)
    >>> posify([x, 1 + x])
    ([_x, _x + 1], {_x: x})
    """
    eq = sympify(eq)
    if iterable(eq):
        f = type(eq)
        eq = list(eq)
        syms = set()
        for e in eq:
            syms = syms.union(e.atoms(Symbol))
        reps = {}
        for s in syms:
            reps.update(dict((v, k) for k, v in posify(s)[1].items()))
        for i, e in enumerate(eq):
            eq[i] = e.subs(reps)
        return f(eq), dict([(r, s) for s, r in reps.items()])

    reps = dict([(s, Dummy(s.name, positive=True))
                 for s in eq.free_symbols if s.is_positive is None])
    eq = eq.subs(reps)
    return eq, dict([(r, s) for s, r in reps.items()])


def hypersimp(f, k):
    """Given combinatorial term f(k) simplify its consecutive term ratio
       i.e. f(k+1)/f(k).  The input term can be composed of functions and
       integer sequences which have equivalent representation in terms
       of gamma special function.

       The algorithm performs three basic steps:

       1. Rewrite all functions in terms of gamma, if possible.

       2. Rewrite all occurrences of gamma in terms of products
          of gamma and rising factorial with integer,  absolute
          constant exponent.

       3. Perform simplification of nested fractions, powers
          and if the resulting expression is a quotient of
          polynomials, reduce their total degree.

       If f(k) is hypergeometric then as result we arrive with a
       quotient of polynomials of minimal degree. Otherwise None
       is returned.

       For more information on the implemented algorithm refer to:

       1. W. Koepf, Algorithms for m-fold Hypergeometric Summation,
          Journal of Symbolic Computation (1995) 20, 399-417
    """
    f = sympify(f)

    g = f.subs(k, k + 1) / f

    g = g.rewrite(gamma)
    g = expand_func(g)
    g = powsimp(g, deep=True, combine='exp')

    if g.is_rational_function(k):
        return simplify(g, ratio=S.Infinity)
    else:
        return None


def hypersimilar(f, g, k):
    """Returns True if 'f' and 'g' are hyper-similar.

       Similarity in hypergeometric sense means that a quotient of
       f(k) and g(k) is a rational function in k.  This procedure
       is useful in solving recurrence relations.

       For more information see hypersimp().

    """
    f, g = list(map(sympify, (f, g)))

    h = (f/g).rewrite(gamma)
    h = h.expand(func=True, basic=False)

    return h.is_rational_function(k)


class _rf(Function):
    @classmethod
    def eval(cls, a, b):
        if b.is_Integer:
            if not b:
                return S.One

            n, result = int(b), S.One

            if n > 0:
                for i in range(n):
                    result *= a + i

                return result
            elif n < 0:
                for i in range(1, -n + 1):
                    result *= a - i

                return 1/result
        else:
            if b.is_Add:
                c, _b = b.as_coeff_Add()

                if c.is_Integer:
                    if c > 0:
                        return _rf(a, _b)*_rf(a + _b, c)
                    elif c < 0:
                        return _rf(a, _b)/_rf(a + _b + c, -c)

            if a.is_Add:
                c, _a = a.as_coeff_Add()

                if c.is_Integer:
                    if c > 0:
                        return _rf(_a, b)*_rf(_a + b, c)/_rf(_a, c)
                    elif c < 0:
                        return _rf(_a, b)*_rf(_a + c, -c)/_rf(_a + b + c, -c)


@timethis('combsimp')
def combsimp(expr):
    r"""
    Simplify combinatorial expressions.

    This function takes as input an expression containing factorials,
    binomials, Pochhammer symbol and other "combinatorial" functions,
    and tries to minimize the number of those functions and reduce
    the size of their arguments. The result is be given in terms of
    binomials and factorials.

    The algorithm works by rewriting all combinatorial functions as
    expressions involving rising factorials (Pochhammer symbols) and
    applies recurrence relations and other transformations applicable
    to rising factorials, to reduce their arguments, possibly letting
    the resulting rising factorial to cancel. Rising factorials with
    the second argument being an integer are expanded into polynomial
    forms and finally all other rising factorial are rewritten in terms
    more familiar functions. If the initial expression contained any
    combinatorial functions, the result is expressed using binomial
    coefficients and gamma functions. If the initial expression consisted
    of gamma functions alone, the result is expressed in terms of gamma
    functions.

    If the result is expressed using gamma functions, the following three
    additional steps are performed:

    1. Reduce the number of gammas by applying the reflection theorem
       gamma(x)*gamma(1-x) == pi/sin(pi*x).
    2. Reduce the number of gammas by applying the multiplication theorem
       gamma(x)*gamma(x+1/n)*...*gamma(x+(n-1)/n) == C*gamma(n*x).
    3. Reduce the number of prefactors by absorbing them into gammas, where
       possible.

    All transformation rules can be found (or was derived from) here:

    1. http://functions.wolfram.com/GammaBetaErf/Pochhammer/17/01/02/
    2. http://functions.wolfram.com/GammaBetaErf/Pochhammer/27/01/0005/

    Examples
    ========

    >>> from sympy.simplify import combsimp
    >>> from sympy import factorial, binomial
    >>> from sympy.abc import n, k

    >>> combsimp(factorial(n)/factorial(n - 3))
    n*(n - 2)*(n - 1)
    >>> combsimp(binomial(n+1, k+1)/binomial(n, k))
    (n + 1)/(k + 1)

    """

    # as a rule of thumb, if the expression contained gammas initially, it
    # probably makes sense to retain them
    as_gamma = not expr.has(factorial, binomial)


    expr = expr.replace(binomial,
        lambda n, k: _rf((n - k + 1).expand(), k.expand())/_rf(1, k.expand()))
    expr = expr.replace(factorial,
        lambda n: _rf(1, n.expand()))
    expr = expr.rewrite(gamma)
    expr = expr.replace(gamma,
        lambda n: _rf(1, (n - 1).expand()))

    if as_gamma:
        expr = expr.replace(_rf,
            lambda a, b: gamma(a + b)/gamma(a))
    else:
        expr = expr.replace(_rf,
            lambda a, b: binomial(a + b - 1, b)*factorial(b))

    def rule(n, k):
        coeff, rewrite = S.One, False

        cn, _n = n.as_coeff_Add()

        if _n and cn.is_Integer and cn:
            coeff *= _rf(_n + 1, cn)/_rf(_n - k + 1, cn)
            rewrite = True
            n = _n

        # this sort of binomial has already been removed by
        # rising factorials but is left here in case the order
        # of rule application is changed
        if k.is_Add:
            ck, _k = k.as_coeff_Add()
            if _k and ck.is_Integer and ck:
                coeff *= _rf(n - ck - _k + 1, ck)/_rf(_k + 1, ck)
                rewrite = True
                k = _k

        if rewrite:
            return coeff*binomial(n, k)

    expr = expr.replace(binomial, rule)

    def rule_gamma(expr, level=0):
        """ Simplify products of gamma functions further. """

        if expr.is_Atom:
            return expr

        def gamma_rat(x):
            # helper to simplify ratios of gammas
            was = x.count(gamma)
            xx = x.replace(gamma, lambda n: _rf(1, (n - 1).expand()
                ).replace(_rf, lambda a, b: gamma(a + b)/gamma(a)))
            if xx.count(gamma) < was:
                x = xx
            return x

        def gamma_factor(x):
            # return True if there is a gamma factor in shallow args
            if x.func is gamma:
                return True
            if x.is_Add or x.is_Mul:
                return any(gamma_factor(xi) for xi in x.args)
            if x.is_Pow and (x.exp.is_integer or x.base.is_positive):
                return gamma_factor(x.base)
            return False

        # recursion step
        if level == 0:
            expr = expr.func(*[rule_gamma(x, level + 1) for x in expr.args])
            level += 1

        if not expr.is_Mul:
            return expr

        # non-commutative step
        if level == 1:
            args, nc = expr.args_cnc()
            if not args:
                return expr
            if nc:
                return rule_gamma(Mul._from_args(args), level + 1)*Mul._from_args(nc)
            level += 1

        # pure gamma handling, not factor absorbtion
        if level == 2:
            sifted = sift(expr.args, gamma_factor)
            gamma_ind = Mul(*sifted.pop(False, []))
            d = Mul(*sifted.pop(True, []))
            assert not sifted

            nd, dd = d.as_numer_denom()
            for ipass in range(2):
                args = list(ordered(Mul.make_args(nd)))
                for i, ni in enumerate(args):
                    if ni.is_Add:
                        ni, dd = Add(*[
                            rule_gamma(gamma_rat(a/dd), level + 1) for a in ni.args]
                            ).as_numer_denom()
                        args[i] = ni
                        if not dd.has(gamma):
                            break
                nd = Mul(*args)
                if ipass ==  0 and not gamma_factor(nd):
                    break
                nd, dd = dd, nd  # now process in reversed order
            expr = gamma_ind*nd/dd
            if not (expr.is_Mul and (gamma_factor(dd) or gamma_factor(nd))):
                return expr
            level += 1

        # iteration until constant
        if level == 3:
            while True:
                was = expr
                expr = rule_gamma(expr, 4)
                if expr == was:
                    return expr

        numer_gammas = []
        denom_gammas = []
        numer_others = []
        denom_others = []
        def explicate(p):
            if p is S.One:
                return None, []
            b, e = p.as_base_exp()
            if e.is_Integer:
                if b.func is gamma:
                    return True, [b.args[0]]*e
                else:
                    return False, [b]*e
            else:
                return False, [p]

        newargs = list(ordered(expr.args))
        while newargs:
            n, d = newargs.pop().as_numer_denom()
            isg, l = explicate(n)
            if isg:
                numer_gammas.extend(l)
            elif isg is False:
                numer_others.extend(l)
            isg, l = explicate(d)
            if isg:
                denom_gammas.extend(l)
            elif isg is False:
                denom_others.extend(l)

        # =========== level 2 work: pure gamma manipulation =========

        # Try to reduce the number of gamma factors by applying the
        # reflection formula gamma(x)*gamma(1-x) = pi/sin(pi*x)
        for gammas, numer, denom in [(
            numer_gammas, numer_others, denom_others),
                (denom_gammas, denom_others, numer_others)]:
            new = []
            while gammas:
                g1 = gammas.pop()
                if g1.is_integer:
                    new.append(g1)
                    continue
                for i, g2 in enumerate(gammas):
                    n = g1 + g2 - 1
                    if not n.is_Integer:
                        continue
                    numer.append(S.Pi)
                    denom.append(sin(S.Pi*g1))
                    gammas.pop(i)
                    if n > 0:
                        for k in range(n):
                            numer.append(1 - g1 + k)
                    elif n < 0:
                        for k in range(-n):
                            denom.append(-g1 - k)
                    break
                else:
                    new.append(g1)
            # /!\ updating IN PLACE
            gammas[:] = new

        # Try to reduce the number of gammas by using the duplication
        # theorem to cancel an upper and lower: gamma(2*s)/gamma(s) =
        # 2**(2*s + 1)/(4*sqrt(pi))*gamma(s + 1/2). Although this could
        # be done with higher argument ratios like gamma(3*x)/gamma(x),
        # this would not reduce the number of gammas as in this case.
        for ng, dg, no, do in [(numer_gammas, denom_gammas, numer_others,
                                denom_others),
                               (denom_gammas, numer_gammas, denom_others,
                                numer_others)]:

            while True:
                for x in ng:
                    for y in dg:
                        n = x - 2*y
                        if n.is_Integer:
                            break
                    else:
                        continue
                    break
                else:
                    break
                ng.remove(x)
                dg.remove(y)
                if n > 0:
                    for k in range(n):
                        no.append(2*y + k)
                elif n < 0:
                    for k in range(-n):
                        do.append(2*y - 1 - k)
                ng.append(y + S(1)/2)
                no.append(2**(2*y - 1))
                do.append(sqrt(S.Pi))

        # Try to reduce the number of gamma factors by applying the
        # multiplication theorem (used when n gammas with args differing
        # by 1/n mod 1 are encountered).
        #
        # run of 2 with args differing by 1/2
        #
        # >>> combsimp(gamma(x)*gamma(x+S.Half))
        # 2*sqrt(2)*2**(-2*x - 1/2)*sqrt(pi)*gamma(2*x)
        #
        # run of 3 args differing by 1/3 (mod 1)
        #
        # >>> combsimp(gamma(x)*gamma(x+S(1)/3)*gamma(x+S(2)/3))
        # 6*3**(-3*x - 1/2)*pi*gamma(3*x)
        # >>> combsimp(gamma(x)*gamma(x+S(1)/3)*gamma(x+S(5)/3))
        # 2*3**(-3*x - 1/2)*pi*(3*x + 2)*gamma(3*x)
        #
        def _run(coeffs):
            # find runs in coeffs such that the difference in terms (mod 1)
            # of t1, t2, ..., tn is 1/n
            u = list(uniq(coeffs))
            for i in range(len(u)):
                dj = ([((u[j] - u[i]) % 1, j) for j in range(i + 1, len(u))])
                for one, j in dj:
                    if one.p == 1 and one.q != 1:
                        n = one.q
                        got = [i]
                        get = list(range(1, n))
                        for d, j in dj:
                            m = n*d
                            if m.is_Integer and m in get:
                                get.remove(m)
                                got.append(j)
                                if not get:
                                    break
                        else:
                            continue
                        for i, j in enumerate(got):
                            c = u[j]
                            coeffs.remove(c)
                            got[i] = c
                        return one.q, got[0], got[1:]

        def _mult_thm(gammas, numer, denom):
            # pull off and analyze the leading coefficient from each gamma arg
            # looking for runs in those Rationals

            # expr -> coeff + resid -> rats[resid] = coeff
            rats = {}
            for g in gammas:
                c, resid = g.as_coeff_Add()
                rats.setdefault(resid, []).append(c)

            # look for runs in Rationals for each resid
            keys = sorted(rats, key=default_sort_key)
            for resid in keys:
                coeffs = list(sorted(rats[resid]))
                new = []
                while True:
                    run = _run(coeffs)
                    if run is None:
                        break

                    # process the sequence that was found:
                    # 1) convert all the gamma functions to have the right
                    #    argument (could be off by an integer)
                    # 2) append the factors corresponding to the theorem
                    # 3) append the new gamma function

                    n, ui, other = run

                    # (1)
                    for u in other:
                        con = resid + u - 1
                        for k in range(int(u - ui)):
                            numer.append(con - k)

                    con = n*(resid + ui)  # for (2) and (3)

                    # (2)
                    numer.append((2*S.Pi)**(S(n - 1)/2)*
                                 n**(S(1)/2 - con))
                    # (3)
                    new.append(con)

                # restore resid to coeffs
                rats[resid] = [resid + c for c in coeffs] + new

            # rebuild the gamma arguments
            g = []
            for resid in keys:
                g += rats[resid]
            # /!\ updating IN PLACE
            gammas[:] = g

        for l, numer, denom in [(numer_gammas, numer_others, denom_others),
                                (denom_gammas, denom_others, numer_others)]:
            _mult_thm(l, numer, denom)

        # =========== level >= 2 work: factor absorbtion =========

        if level >= 2:
            # Try to absorb factors into the gammas: x*gamma(x) -> gamma(x + 1)
            # and gamma(x)/(x - 1) -> gamma(x - 1)
            # This code (in particular repeated calls to find_fuzzy) can be very
            # slow.
            def find_fuzzy(l, x):
                if not l:
                    return
                S1, T1 = compute_ST(x)
                for y in l:
                    S2, T2 = inv[y]
                    if T1 != T2 or (not S1.intersection(S2) and
                                    (S1 != set() or S2 != set())):
                        continue
                    # XXX we want some simplification (e.g. cancel or
                    # simplify) but no matter what it's slow.
                    a = len(cancel(x/y).free_symbols)
                    b = len(x.free_symbols)
                    c = len(y.free_symbols)
                    # TODO is there a better heuristic?
                    if a == 0 and (b > 0 or c > 0):
                        return y

            # We thus try to avoid expensive calls by building the following
            # "invariants": For every factor or gamma function argument
            #   - the set of free symbols S
            #   - the set of functional components T
            # We will only try to absorb if T1==T2 and (S1 intersect S2 != emptyset
            # or S1 == S2 == emptyset)
            inv = {}

            def compute_ST(expr):
                if expr in inv:
                    return inv[expr]
                return (expr.free_symbols, expr.atoms(Function).union(
                        set(e.exp for e in expr.atoms(Pow))))

            def update_ST(expr):
                inv[expr] = compute_ST(expr)
            for expr in numer_gammas + denom_gammas + numer_others + denom_others:
                update_ST(expr)

            for gammas, numer, denom in [(
                numer_gammas, numer_others, denom_others),
                    (denom_gammas, denom_others, numer_others)]:
                new = []
                while gammas:
                    g = gammas.pop()
                    cont = True
                    while cont:
                        cont = False
                        y = find_fuzzy(numer, g)
                        if y is not None:
                            numer.remove(y)
                            if y != g:
                                numer.append(y/g)
                                update_ST(y/g)
                            g += 1
                            cont = True
                        y = find_fuzzy(denom, g - 1)
                        if y is not None:
                            denom.remove(y)
                            if y != g - 1:
                                numer.append((g - 1)/y)
                                update_ST((g - 1)/y)
                            g -= 1
                            cont = True
                    new.append(g)
                # /!\ updating IN PLACE
                gammas[:] = new

        # =========== rebuild expr ==================================

        return Mul(*[gamma(g) for g in numer_gammas]) \
            / Mul(*[gamma(g) for g in denom_gammas]) \
            * Mul(*numer_others) / Mul(*denom_others)

    # (for some reason we cannot use Basic.replace in this case)
    was = factor(expr)
    expr = rule_gamma(was)
    if expr != was:
        expr = factor(expr)

    return expr


def signsimp(expr, evaluate=None):
    """Make all Add sub-expressions canonical wrt sign.

    If an Add subexpression, ``a``, can have a sign extracted,
    as determined by could_extract_minus_sign, it is replaced
    with Mul(-1, a, evaluate=False). This allows signs to be
    extracted from powers and products.

    Examples
    ========

    >>> from sympy import signsimp, exp, symbols
    >>> from sympy.abc import x, y
    >>> i = symbols('i', odd=True)
    >>> n = -1 + 1/x
    >>> n/x/(-n)**2 - 1/n/x
    (-1 + 1/x)/(x*(1 - 1/x)**2) - 1/(x*(-1 + 1/x))
    >>> signsimp(_)
    0
    >>> x*n + x*-n
    x*(-1 + 1/x) + x*(1 - 1/x)
    >>> signsimp(_)
    0

    Since powers automatically handle leading signs

    >>> (-2)**i
    -2**i

    signsimp can be used to put the base of a power with an integer
    exponent into canonical form:

    >>> n**i
    (-1 + 1/x)**i
    >>> signsimp(_)
    -(1 - 1/x)**i

    By default, signsimp doesn't leave behind any hollow simplification:
    if making an Add canonical wrt sign didn't change the expression, the
    original Add is restored. If this is not desired then the keyword
    ``evaluate`` can be set to False:

    >>> e = exp(y - x)
    >>> signsimp(e) == e
    True
    >>> signsimp(e, evaluate=False)
    exp(-(x - y))

    """
    if evaluate is None:
        evaluate = global_evaluate[0]
    expr = sympify(expr)
    if not isinstance(expr, Expr) or expr.is_Atom:
        return expr
    e = sub_post(sub_pre(expr))
    if not isinstance(e, Expr) or e.is_Atom:
        return e
    if e.is_Add:
        return e.func(*[signsimp(a) for a in e.args])
    if evaluate:
        e = e.xreplace(dict([(m, -(-m)) for m in e.atoms(Mul) if -(-m) != m]))
    return e


def simplify(expr, ratio=1.7, measure=count_ops, fu=False):
    """
    Simplifies the given expression.

    Simplification is not a well defined term and the exact strategies
    this function tries can change in the future versions of SymPy. If
    your algorithm relies on "simplification" (whatever it is), try to
    determine what you need exactly  -  is it powsimp()?, radsimp()?,
    together()?, logcombine()?, or something else? And use this particular
    function directly, because those are well defined and thus your algorithm
    will be robust.

    Nonetheless, especially for interactive use, or when you don't know
    anything about the structure of the expression, simplify() tries to apply
    intelligent heuristics to make the input expression "simpler".  For
    example:

    >>> from sympy import simplify, cos, sin
    >>> from sympy.abc import x, y
    >>> a = (x + x**2)/(x*sin(y)**2 + x*cos(y)**2)
    >>> a
    (x**2 + x)/(x*sin(y)**2 + x*cos(y)**2)
    >>> simplify(a)
    x + 1

    Note that we could have obtained the same result by using specific
    simplification functions:

    >>> from sympy import trigsimp, cancel
    >>> trigsimp(a)
    (x**2 + x)/x
    >>> cancel(_)
    x + 1

    In some cases, applying :func:`simplify` may actually result in some more
    complicated expression. The default ``ratio=1.7`` prevents more extreme
    cases: if (result length)/(input length) > ratio, then input is returned
    unmodified.  The ``measure`` parameter lets you specify the function used
    to determine how complex an expression is.  The function should take a
    single argument as an expression and return a number such that if
    expression ``a`` is more complex than expression ``b``, then
    ``measure(a) > measure(b)``.  The default measure function is
    :func:`count_ops`, which returns the total number of operations in the
    expression.

    For example, if ``ratio=1``, ``simplify`` output can't be longer
    than input.

    ::

        >>> from sympy import sqrt, simplify, count_ops, oo
        >>> root = 1/(sqrt(2)+3)

    Since ``simplify(root)`` would result in a slightly longer expression,
    root is returned unchanged instead::

       >>> simplify(root, ratio=1) == root
       True

    If ``ratio=oo``, simplify will be applied anyway::

        >>> count_ops(simplify(root, ratio=oo)) > count_ops(root)
        True

    Note that the shortest expression is not necessary the simplest, so
    setting ``ratio`` to 1 may not be a good idea.
    Heuristically, the default value ``ratio=1.7`` seems like a reasonable
    choice.

    You can easily define your own measure function based on what you feel
    should represent the "size" or "complexity" of the input expression.  Note
    that some choices, such as ``lambda expr: len(str(expr))`` may appear to be
    good metrics, but have other problems (in this case, the measure function
    may slow down simplify too much for very large expressions).  If you don't
    know what a good metric would be, the default, ``count_ops``, is a good
    one.

    For example:

    >>> from sympy import symbols, log
    >>> a, b = symbols('a b', positive=True)
    >>> g = log(a) + log(b) + log(a)*log(1/b)
    >>> h = simplify(g)
    >>> h
    log(a*b**(-log(a) + 1))
    >>> count_ops(g)
    8
    >>> count_ops(h)
    5

    So you can see that ``h`` is simpler than ``g`` using the count_ops metric.
    However, we may not like how ``simplify`` (in this case, using
    ``logcombine``) has created the ``b**(log(1/a) + 1)`` term.  A simple way
    to reduce this would be to give more weight to powers as operations in
    ``count_ops``.  We can do this by using the ``visual=True`` option:

    >>> print(count_ops(g, visual=True))
    2*ADD + DIV + 4*LOG + MUL
    >>> print(count_ops(h, visual=True))
    2*LOG + MUL + POW + SUB

    >>> from sympy import Symbol, S
    >>> def my_measure(expr):
    ...     POW = Symbol('POW')
    ...     # Discourage powers by giving POW a weight of 10
    ...     count = count_ops(expr, visual=True).subs(POW, 10)
    ...     # Every other operation gets a weight of 1 (the default)
    ...     count = count.replace(Symbol, type(S.One))
    ...     return count
    >>> my_measure(g)
    8
    >>> my_measure(h)
    14
    >>> 15./8 > 1.7 # 1.7 is the default ratio
    True
    >>> simplify(g, measure=my_measure)
    -log(a)*log(b) + log(a) + log(b)

    Note that because ``simplify()`` internally tries many different
    simplification strategies and then compares them using the measure
    function, we get a completely different result that is still different
    from the input expression by doing this.
    """
    expr = sympify(expr)

    try:
        return expr._eval_simplify(ratio=ratio, measure=measure)
    except AttributeError:
        pass

    original_expr = expr = signsimp(expr)

    from sympy.simplify.hyperexpand import hyperexpand
    from sympy.functions.special.bessel import BesselBase
    from sympy import Sum, Product

    if not isinstance(expr, Basic) or not expr.args:  # XXX: temporary hack
        return expr

    if not isinstance(expr, (Add, Mul, Pow, ExpBase)):
        return expr.func(*[simplify(x, ratio=ratio, measure=measure, fu=fu)
                         for x in expr.args])

    # TODO: Apply different strategies, considering expression pattern:
    # is it a purely rational function? Is there any trigonometric function?...
    # See also https://github.com/sympy/sympy/pull/185.

    def shorter(*choices):
        '''Return the choice that has the fewest ops. In case of a tie,
        the expression listed first is selected.'''
        if not has_variety(choices):
            return choices[0]
        return min(choices, key=measure)

    expr = bottom_up(expr, lambda w: w.normal())
    expr = Mul(*powsimp(expr).as_content_primitive())
    _e = cancel(expr)
    expr1 = shorter(_e, _mexpand(_e).cancel())  # issue 6829
    expr2 = shorter(together(expr, deep=True), together(expr1, deep=True))

    if ratio is S.Infinity:
        expr = expr2
    else:
        expr = shorter(expr2, expr1, expr)
    if not isinstance(expr, Basic):  # XXX: temporary hack
        return expr

    expr = factor_terms(expr, sign=False)

    # hyperexpand automatically only works on hypergeometric terms
    expr = hyperexpand(expr)

    expr = piecewise_fold(expr)

    if expr.has(BesselBase):
        expr = besselsimp(expr)

    if expr.has(TrigonometricFunction) and not fu or expr.has(
            HyperbolicFunction):
        expr = trigsimp(expr, deep=True)

    if expr.has(log):
        expr = shorter(expand_log(expr, deep=True), logcombine(expr))

    if expr.has(CombinatorialFunction, gamma):
        expr = combsimp(expr)

    if expr.has(Sum):
        expr = sum_simplify(expr)

    if expr.has(Product):
        expr = product_simplify(expr)

    short = shorter(powsimp(expr, combine='exp', deep=True), powsimp(expr), expr)
    short = shorter(short, factor_terms(short), expand_power_exp(expand_mul(short)))
    if short.has(TrigonometricFunction, HyperbolicFunction, ExpBase):
        short = exptrigsimp(short, simplify=False)

    # get rid of hollow 2-arg Mul factorization
    hollow_mul = Transform(
        lambda x: Mul(*x.args),
        lambda x:
        x.is_Mul and
        len(x.args) == 2 and
        x.args[0].is_Number and
        x.args[1].is_Add and
        x.is_commutative)
    expr = short.xreplace(hollow_mul)

    numer, denom = expr.as_numer_denom()
    if denom.is_Add:
        n, d = fraction(radsimp(1/denom, symbolic=False, max_terms=1))
        if n is not S.One:
            expr = (numer*n).expand()/d

    if expr.could_extract_minus_sign():
        n, d = fraction(expr)
        if d != 0:
            expr = signsimp(-n/(-d))

    if measure(expr) > ratio*measure(original_expr):
        expr = original_expr

    return expr


def sum_simplify(s):
    """Main function for Sum simplification"""
    from sympy.concrete.summations import Sum

    terms = Add.make_args(s)
    s_t = [] # Sum Terms
    o_t = [] # Other Terms

    for term in terms:
        if isinstance(term, Mul):
            constant = 1
            other = 1
            s = 0
            n_sum_terms = 0
            for j in range(len(term.args)):
                if isinstance(term.args[j], Sum):
                    s = term.args[j]
                    n_sum_terms = n_sum_terms + 1
                elif term.args[j].is_number == True:
                    constant = constant * term.args[j]
                else:
                    other = other * term.args[j]
            if other == 1 and n_sum_terms == 1:
                # Insert the constant inside the Sum
                s_t.append(Sum(constant * s.function, *s.limits))
            elif other != 1 and n_sum_terms == 1:
                o_t.append(other * Sum(constant * s.function, *s.limits))
            else:
                o_t.append(term)
        elif isinstance(term, Sum):
            s_t.append(term)
        else:
            o_t.append(term)

    used = [False] * len(s_t)

    for method in range(2):
        for i, s_term1 in enumerate(s_t):
            if not used[i]:
                for j, s_term2 in enumerate(s_t):
                    if not used[j] and i != j:
                        temp = sum_add(s_term1, s_term2, method)
                        if isinstance(temp, Sum):
                            s_t[i] = temp
                            s_term1 = s_t[i]
                            used[j] = True

    result = Add(*o_t)

    for i, s_term in enumerate(s_t):
        if not used[i]:
            result = Add(result, s_term)

    return result


def sum_add(self, other, method=0):
    """Helper function for Sum simplification"""
    from sympy.concrete.summations import Sum

    if type(self) == type(other):
        if method == 0:
            if self.limits == other.limits:
                return Sum(self.function + other.function, *self.limits)
        elif method == 1:
            if simplify(self.function - other.function) == 0:
                if len(self.limits) == len(other.limits) == 1:
                    i = self.limits[0][0]
                    x1 = self.limits[0][1]
                    y1 = self.limits[0][2]
                    j = other.limits[0][0]
                    x2 = other.limits[0][1]
                    y2 = other.limits[0][2]

                    if i == j:
                        if x2 == y1 + 1:
                            return Sum(self.function, (i, x1, y2))
                        elif x1 == y2 + 1:
                            return Sum(self.function, (i, x2, y1))

    return Add(self, other)


def product_simplify(s):
    """Main function for Product simplification"""
    from sympy.concrete.products import Product

    terms = Mul.make_args(s)
    p_t = [] # Product Terms
    o_t = [] # Other Terms

    for term in terms:
        if isinstance(term, Product):
            p_t.append(term)
        else:
            o_t.append(term)

    used = [False] * len(p_t)

    for method in range(2):
        for i, p_term1 in enumerate(p_t):
            if not used[i]:
                for j, p_term2 in enumerate(p_t):
                    if not used[j] and i != j:
                        if isinstance(product_mul(p_term1, p_term2, method), Product):
                            p_t[i] = product_mul(p_term1, p_term2, method)
                            used[j] = True

    result = Mul(*o_t)

    for i, p_term in enumerate(p_t):
        if not used[i]:
            result = Mul(result, p_term)

    return result


def product_mul(self, other, method=0):
    """Helper function for Product simplification"""
    from sympy.concrete.products import Product

    if type(self) == type(other):
        if method == 0:
            if self.limits == other.limits:
                return Product(self.function * other.function, *self.limits)
        elif method == 1:
            if simplify(self.function - other.function) == 0:
                if len(self.limits) == len(other.limits) == 1:
                    i = self.limits[0][0]
                    x1 = self.limits[0][1]
                    y1 = self.limits[0][2]
                    j = other.limits[0][0]
                    x2 = other.limits[0][1]
                    y2 = other.limits[0][2]

                    if i == j:
                        if x2 == y1 + 1:
                            return Product(self.function, (i, x1, y2))
                        elif x1 == y2 + 1:
                            return Product(self.function, (i, x2, y1))

    return Mul(self, other)


def _nthroot_solve(p, n, prec):
    """
     helper function for ``nthroot``
     It denests ``p**Rational(1, n)`` using its minimal polynomial
    """
    from sympy.polys.numberfields import _minimal_polynomial_sq
    from sympy.solvers import solve
    while n % 2 == 0:
        p = sqrtdenest(sqrt(p))
        n = n // 2
    if n == 1:
        return p
    pn = p**Rational(1, n)
    x = Symbol('x')
    f = _minimal_polynomial_sq(p, n, x)
    if f is None:
        return None
    sols = solve(f, x)
    for sol in sols:
        if abs(sol - pn).n() < 1./10**prec:
            sol = sqrtdenest(sol)
            if _mexpand(sol**n) == p:
                return sol


def logcombine(expr, force=False):
    """
    Takes logarithms and combines them using the following rules:

    - log(x) + log(y) == log(x*y) if both are not negative
    - a*log(x) == log(x**a) if x is positive and a is real

    If ``force`` is True then the assumptions above will be assumed to hold if
    there is no assumption already in place on a quantity. For example, if
    ``a`` is imaginary or the argument negative, force will not perform a
    combination but if ``a`` is a symbol with no assumptions the change will
    take place.

    Examples
    ========

    >>> from sympy import Symbol, symbols, log, logcombine, I
    >>> from sympy.abc import a, x, y, z
    >>> logcombine(a*log(x) + log(y) - log(z))
    a*log(x) + log(y) - log(z)
    >>> logcombine(a*log(x) + log(y) - log(z), force=True)
    log(x**a*y/z)
    >>> x,y,z = symbols('x,y,z', positive=True)
    >>> a = Symbol('a', real=True)
    >>> logcombine(a*log(x) + log(y) - log(z))
    log(x**a*y/z)

    The transformation is limited to factors and/or terms that
    contain logs, so the result depends on the initial state of
    expansion:

    >>> eq = (2 + 3*I)*log(x)
    >>> logcombine(eq, force=True) == eq
    True
    >>> logcombine(eq.expand(), force=True)
    log(x**2) + I*log(x**3)

    See Also
    ========
    posify: replace all symbols with symbols having positive assumptions

    """

    def f(rv):
        if not (rv.is_Add or rv.is_Mul):
            return rv

        def gooda(a):
            # bool to tell whether the leading ``a`` in ``a*log(x)``
            # could appear as log(x**a)
            return (a is not S.NegativeOne and  # -1 *could* go, but we disallow
                (a.is_real or force and a.is_real is not False))

        def goodlog(l):
            # bool to tell whether log ``l``'s argument can combine with others
            a = l.args[0]
            return a.is_positive or force and a.is_nonpositive is not False

        other = []
        logs = []
        log1 = defaultdict(list)
        for a in Add.make_args(rv):
            if a.func is log and goodlog(a):
                log1[()].append(([], a))
            elif not a.is_Mul:
                other.append(a)
            else:
                ot = []
                co = []
                lo = []
                for ai in a.args:
                    if ai.is_Rational and ai < 0:
                        ot.append(S.NegativeOne)
                        co.append(-ai)
                    elif ai.func is log and goodlog(ai):
                        lo.append(ai)
                    elif gooda(ai):
                        co.append(ai)
                    else:
                        ot.append(ai)
                if len(lo) > 1:
                    logs.append((ot, co, lo))
                elif lo:
                    log1[tuple(ot)].append((co, lo[0]))
                else:
                    other.append(a)

        # if there is only one log at each coefficient and none have
        # an exponent to place inside the log then there is nothing to do
        if not logs and all(len(log1[k]) == 1 and log1[k][0] == [] for k in log1):
            return rv

        # collapse multi-logs as far as possible in a canonical way
        # TODO: see if x*log(a)+x*log(a)*log(b) -> x*log(a)*(1+log(b))?
        # -- in this case, it's unambiguous, but if it were were a log(c) in
        # each term then it's arbitrary whether they are grouped by log(a) or
        # by log(c). So for now, just leave this alone; it's probably better to
        # let the user decide
        for o, e, l in logs:
            l = list(ordered(l))
            e = log(l.pop(0).args[0]**Mul(*e))
            while l:
                li = l.pop(0)
                e = log(li.args[0]**e)
            c, l = Mul(*o), e
            if l.func is log:  # it should be, but check to be sure
                log1[(c,)].append(([], l))
            else:
                other.append(c*l)

        # logs that have the same coefficient can multiply
        for k in list(log1.keys()):
            log1[Mul(*k)] = log(logcombine(Mul(*[
                l.args[0]**Mul(*c) for c, l in log1.pop(k)]),
                force=force))

        # logs that have oppositely signed coefficients can divide
        for k in ordered(list(log1.keys())):
            if not k in log1:  # already popped as -k
                continue
            if -k in log1:
                # figure out which has the minus sign; the one with
                # more op counts should be the one
                num, den = k, -k
                if num.count_ops() > den.count_ops():
                    num, den = den, num
                other.append(num*log(log1.pop(num).args[0]/log1.pop(den).args[0]))
            else:
                other.append(k*log1.pop(k))

        return Add(*other)

    return bottom_up(expr, f)


def bottom_up(rv, F, atoms=False, nonbasic=False):
    """Apply ``F`` to all expressions in an expression tree from the
    bottom up. If ``atoms`` is True, apply ``F`` even if there are no args;
    if ``nonbasic`` is True, try to apply ``F`` to non-Basic objects.
    """
    try:
        if rv.args:
            args = tuple([bottom_up(a, F, atoms, nonbasic)
                for a in rv.args])
            if args != rv.args:
                rv = rv.func(*args)
            rv = F(rv)
        elif atoms:
            rv = F(rv)
    except AttributeError:
        if nonbasic:
            try:
                rv = F(rv)
            except TypeError:
                pass

    return rv


def nthroot(expr, n, max_len=4, prec=15):
    """
    compute a real nth-root of a sum of surds

    Parameters
    ==========

    expr : sum of surds
    n : integer
    max_len : maximum number of surds passed as constants to ``nsimplify``

    Algorithm
    =========

    First ``nsimplify`` is used to get a candidate root; if it is not a
    root the minimal polynomial is computed; the answer is one of its
    roots.

    Examples
    ========

    >>> from sympy.simplify.simplify import nthroot
    >>> from sympy import Rational, sqrt
    >>> nthroot(90 + 34*sqrt(7), 3)
    sqrt(7) + 3

    """
    expr = sympify(expr)
    n = sympify(n)
    p = expr**Rational(1, n)
    if not n.is_integer:
        return p
    if not _is_sum_surds(expr):
        return p
    surds = []
    coeff_muls = [x.as_coeff_Mul() for x in expr.args]
    for x, y in coeff_muls:
        if not x.is_rational:
            return p
        if y is S.One:
            continue
        if not (y.is_Pow and y.exp == S.Half and y.base.is_integer):
            return p
        surds.append(y)
    surds.sort()
    surds = surds[:max_len]
    if expr < 0 and n % 2 == 1:
        p = (-expr)**Rational(1, n)
        a = nsimplify(p, constants=surds)
        res = a if _mexpand(a**n) == _mexpand(-expr) else p
        return -res
    a = nsimplify(p, constants=surds)
    if _mexpand(a) is not _mexpand(p) and _mexpand(a**n) == _mexpand(expr):
        return _mexpand(a)
    expr = _nthroot_solve(expr, n, prec)
    if expr is None:
        return p
    return expr


def nsimplify(expr, constants=[], tolerance=None, full=False, rational=None):
    """
    Find a simple representation for a number or, if there are free symbols or
    if rational=True, then replace Floats with their Rational equivalents. If
    no change is made and rational is not False then Floats will at least be
    converted to Rationals.

    For numerical expressions, a simple formula that numerically matches the
    given numerical expression is sought (and the input should be possible
    to evalf to a precision of at least 30 digits).

    Optionally, a list of (rationally independent) constants to
    include in the formula may be given.

    A lower tolerance may be set to find less exact matches. If no tolerance
    is given then the least precise value will set the tolerance (e.g. Floats
    default to 15 digits of precision, so would be tolerance=10**-15).

    With full=True, a more extensive search is performed
    (this is useful to find simpler numbers when the tolerance
    is set low).

    Examples
    ========

    >>> from sympy import nsimplify, sqrt, GoldenRatio, exp, I, exp, pi
    >>> nsimplify(4/(1+sqrt(5)), [GoldenRatio])
    -2 + 2*GoldenRatio
    >>> nsimplify((1/(exp(3*pi*I/5)+1)))
    1/2 - I*sqrt(sqrt(5)/10 + 1/4)
    >>> nsimplify(I**I, [pi])
    exp(-pi/2)
    >>> nsimplify(pi, tolerance=0.01)
    22/7

    See Also
    ========
    sympy.core.function.nfloat

    """
    try:
        return sympify(as_int(expr))
    except (TypeError, ValueError):
        pass
    expr = sympify(expr)
    if rational or expr.free_symbols:
        return _real_to_rational(expr, tolerance)

    # SymPy's default tolerance for Rationals is 15; other numbers may have
    # lower tolerances set, so use them to pick the largest tolerance if None
    # was given
    if tolerance is None:
        tolerance = 10**-min([15] +
             [mpmath.libmp.libmpf.prec_to_dps(n._prec)
             for n in expr.atoms(Float)])
    # XXX should prec be set independent of tolerance or should it be computed
    # from tolerance?
    prec = 30
    bprec = int(prec*3.33)

    constants_dict = {}
    for constant in constants:
        constant = sympify(constant)
        v = constant.evalf(prec)
        if not v.is_Float:
            raise ValueError("constants must be real-valued")
        constants_dict[str(constant)] = v._to_mpmath(bprec)

    exprval = expr.evalf(prec, chop=True)
    re, im = exprval.as_real_imag()

    # safety check to make sure that this evaluated to a number
    if not (re.is_Number and im.is_Number):
        return expr

    def nsimplify_real(x):
        orig = mpmath.mp.dps
        xv = x._to_mpmath(bprec)
        try:
            # We'll be happy with low precision if a simple fraction
            if not (tolerance or full):
                mpmath.mp.dps = 15
                rat = mpmath.findpoly(xv, 1)
                if rat is not None:
                    return Rational(-int(rat[1]), int(rat[0]))
            mpmath.mp.dps = prec
            newexpr = mpmath.identify(xv, constants=constants_dict,
                tol=tolerance, full=full)
            if not newexpr:
                raise ValueError
            if full:
                newexpr = newexpr[0]
            expr = sympify(newexpr)
            if x and not expr:  # don't let x become 0
                raise ValueError
            if expr.is_finite is False and not xv in [mpmath.inf, mpmath.ninf]:
                raise ValueError
            return expr
        finally:
            # even though there are returns above, this is executed
            # before leaving
            mpmath.mp.dps = orig
    try:
        if re:
            re = nsimplify_real(re)
        if im:
            im = nsimplify_real(im)
    except ValueError:
        if rational is None:
            return _real_to_rational(expr)
        return expr

    rv = re + im*S.ImaginaryUnit
    # if there was a change or rational is explicitly not wanted
    # return the value, else return the Rational representation
    if rv != expr or rational is False:
        return rv
    return _real_to_rational(expr)


def _real_to_rational(expr, tolerance=None):
    """
    Replace all reals in expr with rationals.

    >>> from sympy import nsimplify
    >>> from sympy.abc import x

    >>> nsimplify(.76 + .1*x**.5, rational=True)
    sqrt(x)/10 + 19/25

    """
    p = expr
    reps = {}
    reduce_num = None
    if tolerance is not None and tolerance < 1:
        reduce_num = ceiling(1/tolerance)
    for float in p.atoms(Float):
        key = float
        if reduce_num is not None:
            r = Rational(float).limit_denominator(reduce_num)
        elif (tolerance is not None and tolerance >= 1 and
                float.is_Integer is False):
            r = Rational(tolerance*round(float/tolerance)
                ).limit_denominator(int(tolerance))
        else:
            r = nsimplify(float, rational=False)
            # e.g. log(3).n() -> log(3) instead of a Rational
            if float and not r:
                r = Rational(float)
            elif not r.is_Rational:
                if float < 0:
                    float = -float
                    d = Pow(10, int((mpmath.log(float)/mpmath.log(10))))
                    r = -Rational(str(float/d))*d
                elif float > 0:
                    d = Pow(10, int((mpmath.log(float)/mpmath.log(10))))
                    r = Rational(str(float/d))*d
                else:
                    r = Integer(0)
        reps[key] = r
    return p.subs(reps, simultaneous=True)
