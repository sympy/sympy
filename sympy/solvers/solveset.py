"""
This module contains functions to:

    - solve a single equation for a single variable, in any domain either real or complex.

    - solve a system of linear equations with N variables and M equations.

    - solve a system of Non Linear Equations with N variables and M equations
"""
from __future__ import print_function, division

from sympy.core.sympify import sympify
from sympy.core import S, Pow, Dummy, pi, Expr, Wild, Mul, Equality
from sympy.core.numbers import I, Number, Rational, oo
from sympy.core.function import (Lambda, expand, expand_complex)
from sympy.core.relational import Eq
from sympy.simplify.simplify import simplify, fraction, trigsimp
from sympy.core.symbol import Symbol
from sympy.functions import (log, Abs, tan, cot, sin, cos, sec, csc, exp,
                             acos, asin, acsc, asec, arg,
                             piecewise_fold)
from sympy.functions.elementary.trigonometric import (TrigonometricFunction,
                                                      HyperbolicFunction)
from sympy.functions.elementary.miscellaneous import real_root
from sympy.sets import (FiniteSet, EmptySet, imageset, Interval, Intersection,
                        Union, ConditionSet, ImageSet)
from sympy.matrices import Matrix
from sympy.polys import (roots, Poly, degree, together, PolynomialError,
                         RootOf)
from sympy.solvers.solvers import checksol, denoms, unrad
from sympy.solvers.polysys import solve_poly_system
from sympy.solvers.inequalities import solve_univariate_inequality
from sympy.utilities import filldedent


def _invert(f_x, y, x, domain=S.Complexes):
    """
    Reduce the complex valued equation ``f(x) = y`` to a set of equations
    ``{g(x) = h_1(y), g(x) = h_2(y), ..., g(x) = h_n(y) }`` where ``g(x)`` is
    a simpler function than ``f(x)``.  The return value is a tuple ``(g(x),
    set_h)``, where ``g(x)`` is a function of ``x`` and ``set_h`` is
    the set of function ``{h_1(y), h_2(y), ..., h_n(y)}``.
    Here, ``y`` is not necessarily a symbol.

    The ``set_h`` contains the functions along with the information
    about their domain in which they are valid, through set
    operations. For instance, if ``y = Abs(x) - n``, is inverted
    in the real domain, then, the ``set_h`` doesn't simply return
    `{-n, n}`, as the nature of `n` is unknown; rather it will return:
    `Intersection([0, oo) {n}) U Intersection((-oo, 0], {-n})`

    By default, the complex domain is used but note that inverting even
    seemingly simple functions like ``exp(x)`` can give very different
    result in the complex domain than are obtained in the real domain.
    (In the case of ``exp(x)``, the inversion via ``log`` is multi-valued
    in the complex domain, having infinitely many branches.)

    If you are working with real values only (or you are not sure which
    function to use) you should probably use set the domain to
    ``S.Reals`` (or use `invert\_real` which does that automatically).


    Examples
    ========

    >>> from sympy.solvers.solveset import invert_complex, invert_real
    >>> from sympy.abc import x, y
    >>> from sympy import exp, log

    When does exp(x) == y?

    >>> invert_complex(exp(x), y, x)
    (x, ImageSet(Lambda(_n, I*(2*_n*pi + arg(y)) + log(Abs(y))), Integers()))
    >>> invert_real(exp(x), y, x)
    (x, Intersection((-oo, oo), {log(y)}))

    When does exp(x) == 1?

    >>> invert_complex(exp(x), 1, x)
    (x, ImageSet(Lambda(_n, 2*_n*I*pi), Integers()))
    >>> invert_real(exp(x), 1, x)
    (x, {0})

    See Also
    ========
    invert_real, invert_complex
    """
    x = sympify(x)
    if not x.is_Symbol:
        raise ValueError("x must be a symbol")
    f_x = sympify(f_x)
    if not f_x.has(x):
        raise ValueError("Inverse of constant function doesn't exist")
    y = sympify(y)
    if y.has(x):
        raise ValueError("y should be independent of x ")

    if domain.is_subset(S.Reals):
        x, s = _invert_real(f_x, FiniteSet(y), x)
    else:
        x, s = _invert_complex(f_x, FiniteSet(y), x)
    return x, s.intersection(domain) if isinstance(s, FiniteSet) else s


invert_complex = _invert


def invert_real(f_x, y, x, domain=S.Reals):
    """
    Inverts a real-valued function. Same as _invert, but sets
    the domain to ``S.Reals`` before inverting.
    """
    return _invert(f_x, y, x, domain)


def _invert_real(f, g_ys, symbol):
    """Helper function for _invert."""

    if f == symbol:
        return (f, g_ys)

    n = Dummy('n', real=True)

    if hasattr(f, 'inverse') and not isinstance(f, (
            TrigonometricFunction,
            HyperbolicFunction,
            )):
        if len(f.args) > 1:
            raise ValueError("Only functions with one argument are supported.")
        return _invert_real(f.args[0],
                            imageset(Lambda(n, f.inverse()(n)), g_ys),
                            symbol)

    if isinstance(f, Abs):
        pos = Interval(0, S.Infinity)
        neg = Interval(S.NegativeInfinity, 0)
        return _invert_real(f.args[0],
                    Union(imageset(Lambda(n, n), g_ys).intersect(pos),
                          imageset(Lambda(n, -n), g_ys).intersect(neg)), symbol)

    if f.is_Add:
        # f = g + h
        g, h = f.as_independent(symbol)
        if g is not S.Zero:
            return _invert_real(h, imageset(Lambda(n, n - g), g_ys), symbol)

    if f.is_Mul:
        # f = g*h
        g, h = f.as_independent(symbol)

        if g is not S.One:
            return _invert_real(h, imageset(Lambda(n, n/g), g_ys), symbol)

    if f.is_Pow:
        base, expo = f.args
        base_has_sym = base.has(symbol)
        expo_has_sym = expo.has(symbol)

        if not expo_has_sym:
            res = imageset(Lambda(n, real_root(n, expo)), g_ys)
            if expo.is_rational:
                numer, denom = expo.as_numer_denom()
                if numer == S.One or numer == - S.One:
                    return _invert_real(base, res, symbol)
                else:
                    if numer % 2 == 0:
                        n = Dummy('n')
                        neg_res = imageset(Lambda(n, -n), res)
                        return _invert_real(base, res + neg_res, symbol)
                    else:
                        return _invert_real(base, res, symbol)
            else:
                if not base.is_positive:
                    raise ValueError("x**w where w is irrational is not "
                                     "defined for negative x")
                return _invert_real(base, res, symbol)

        if not base_has_sym:
            return _invert_real(expo,
                imageset(Lambda(n, log(n)/log(base)), g_ys), symbol)

    if isinstance(f, TrigonometricFunction):
        if isinstance(g_ys, FiniteSet):
            def inv(trig):
                if isinstance(f, (sin, csc)):
                    F = asin if isinstance(f, sin) else acsc
                    return (lambda a: n*pi + (-1)**n*F(a),)
                if isinstance(f, (cos, sec)):
                    F = acos if isinstance(f, cos) else asec
                    return (
                        lambda a: 2*n*pi + F(a),
                        lambda a: 2*n*pi - F(a),)
                if isinstance(f, (tan, cot)):
                    return (lambda a: n*pi + f.inverse()(a),)

            n = Dummy('n', integer=True)
            invs = S.EmptySet
            for L in inv(f):
                invs += Union(*[imageset(Lambda(n, L(g)), S.Integers) for g in g_ys])
            return _invert_real(f.args[0], invs, symbol)

    return (f, g_ys)


def _invert_complex(f, g_ys, symbol):
    """Helper function for _invert."""

    if f == symbol:
        return (f, g_ys)

    n = Dummy('n')

    if f.is_Add:
        # f = g + h
        g, h = f.as_independent(symbol)
        if g is not S.Zero:
            return _invert_complex(h, imageset(Lambda(n, n - g), g_ys), symbol)

    if f.is_Mul:
        # f = g*h
        g, h = f.as_independent(symbol)

        if g is not S.One:
            return _invert_complex(h, imageset(Lambda(n, n/g), g_ys), symbol)

    if hasattr(f, 'inverse') and \
       not isinstance(f, TrigonometricFunction) and \
       not isinstance(f, exp):
        if len(f.args) > 1:
            raise ValueError("Only functions with one argument are supported.")
        return _invert_complex(f.args[0],
                               imageset(Lambda(n, f.inverse()(n)), g_ys), symbol)

    if isinstance(f, exp):
        if isinstance(g_ys, FiniteSet):
            exp_invs = Union(*[imageset(Lambda(n, I*(2*n*pi + arg(g_y)) +
                                               log(Abs(g_y))), S.Integers)
                               for g_y in g_ys if g_y != 0])
            return _invert_complex(f.args[0], exp_invs, symbol)

    return (f, g_ys)


def domain_check(f, symbol, p):
    """Returns False if point p is infinite or any subexpression of f
    is infinite or becomes so after replacing symbol with p. If none of
    these conditions is met then True will be returned.

    Examples
    ========

    >>> from sympy import Mul, oo
    >>> from sympy.abc import x
    >>> from sympy.solvers.solveset import domain_check
    >>> g = 1/(1 + (1/(x + 1))**2)
    >>> domain_check(g, x, -1)
    False
    >>> domain_check(x**2, x, 0)
    True
    >>> domain_check(1/x, x, oo)
    False

    * The function relies on the assumption that the original form
      of the equation has not been changed by automatic simplification.

    >>> domain_check(x/x, x, 0) # x/x is automatically simplified to 1
    True

    * To deal with automatic evaluations use evaluate=False:

    >>> domain_check(Mul(x, 1/x, evaluate=False), x, 0)
    False
    """
    f, p = sympify(f), sympify(p)
    if p.is_infinite:
        return False
    return _domain_check(f, symbol, p)


def _domain_check(f, symbol, p):
    # helper for domain check
    if f.is_Atom and f.is_finite:
        return True
    elif f.subs(symbol, p).is_infinite:
        return False
    else:
        return all([_domain_check(g, symbol, p)
                    for g in f.args])


def _is_finite_with_finite_vars(f, domain=S.Complexes):
    """
    Return True if the given expression is finite. For symbols that
    don't assign a value for `complex` and/or `real`, the domain will
    be used to assign a value; symbols that don't assign a value
    for `finite` will be made finite. All other assumptions are
    left unmodified.
    """
    def assumptions(s):
        A = s.assumptions0
        if A.get('finite', None) is None:
            A['finite'] = True
        A.setdefault('complex', True)
        A.setdefault('real', domain.is_subset(S.Reals))
        return A

    reps = {s: Dummy(**assumptions(s)) for s in f.free_symbols}
    return f.xreplace(reps).is_finite


def _is_function_class_equation(func_class, f, symbol):
    """ Tests whether the equation is an equation of the given function class.

    The given equation belongs to the given function class if it is
    comprised of functions of the function class which are multiplied by
    or added to expressions independent of the symbol. In addition, the
    arguments of all such functions must be linear in the symbol as well.

    Examples
    ========

    >>> from sympy.solvers.solveset import _is_function_class_equation
    >>> from sympy import tan, sin, tanh, sinh, exp
    >>> from sympy.abc import x
    >>> from sympy.functions.elementary.trigonometric import (TrigonometricFunction,
    ... HyperbolicFunction)
    >>> _is_function_class_equation(TrigonometricFunction, exp(x) + tan(x), x)
    False
    >>> _is_function_class_equation(TrigonometricFunction, tan(x) + sin(x), x)
    True
    >>> _is_function_class_equation(TrigonometricFunction, tan(x**2), x)
    False
    >>> _is_function_class_equation(TrigonometricFunction, tan(x + 2), x)
    True
    >>> _is_function_class_equation(HyperbolicFunction, tanh(x) + sinh(x), x)
    True
    """
    if f.is_Mul or f.is_Add:
        return all(_is_function_class_equation(func_class, arg, symbol)
                   for arg in f.args)

    if f.is_Pow:
        if not f.exp.has(symbol):
            return _is_function_class_equation(func_class, f.base, symbol)
        else:
            return False

    if not f.has(symbol):
        return True

    if isinstance(f, func_class):
        try:
            g = Poly(f.args[0], symbol)
            return g.degree() <= 1
        except PolynomialError:
            return False
    else:
        return False


def _solve_as_rational(f, symbol, domain):
    """ solve rational functions"""
    f = together(f, deep=True)
    g, h = fraction(f)
    if not h.has(symbol):
        return _solve_as_poly(g, symbol, domain)
    else:
        valid_solns = _solveset(g, symbol, domain)
        invalid_solns = _solveset(h, symbol, domain)
        return valid_solns - invalid_solns


def _solve_trig(f, symbol, domain):
    """ Helper to solve trigonometric equations """
    f = trigsimp(f)
    f_original = f
    f = f.rewrite(exp)
    f = together(f)
    g, h = fraction(f)
    y = Dummy('y')
    g, h = g.expand(), h.expand()
    g, h = g.subs(exp(I*symbol), y), h.subs(exp(I*symbol), y)
    if g.has(symbol) or h.has(symbol):
        return ConditionSet(symbol, Eq(f, 0), S.Reals)

    solns = solveset_complex(g, y) - solveset_complex(h, y)

    if isinstance(solns, FiniteSet):
        result = Union(*[invert_complex(exp(I*symbol), s, symbol)[1]
                       for s in solns])
        return Intersection(result, domain)
    elif solns is S.EmptySet:
        return S.EmptySet
    else:
        return ConditionSet(symbol, Eq(f_original, 0), S.Reals)


def _solve_as_poly(f, symbol, domain=S.Complexes):
    """
    Solve the equation using polynomial techniques if it already is a
    polynomial equation or, with a change of variables, can be made so.
    """
    result = None
    if f.is_polynomial(symbol):
        solns = roots(f, symbol, cubics=True, quartics=True,
                      quintics=True, domain='EX')
        num_roots = sum(solns.values())
        if degree(f, symbol) <= num_roots:
            result = FiniteSet(*solns.keys())
        else:
            poly = Poly(f, symbol)
            solns = poly.all_roots()
            if poly.degree() <= len(solns):
                result = FiniteSet(*solns)
            else:
                result = ConditionSet(symbol, Eq(f, 0), domain)
    else:
        poly = Poly(f)
        if poly is None:
            result = ConditionSet(symbol, Eq(f, 0), domain)
        gens = [g for g in poly.gens if g.has(symbol)]

        if len(gens) == 1:
            poly = Poly(poly, gens[0])
            gen = poly.gen
            deg = poly.degree()
            poly = Poly(poly.as_expr(), poly.gen, composite=True)
            poly_solns = FiniteSet(*roots(poly, cubics=True, quartics=True,
                                          quintics=True).keys())

            if len(poly_solns) < deg:
                result = ConditionSet(symbol, Eq(f, 0), domain)

            if gen != symbol:
                y = Dummy('y')
                inverter = invert_real if domain.is_subset(S.Reals) else invert_complex
                lhs, rhs_s = inverter(gen, y, symbol)
                if lhs == symbol:
                    result = Union(*[rhs_s.subs(y, s) for s in poly_solns])
                else:
                    result = ConditionSet(symbol, Eq(f, 0), domain)
        else:
            result = ConditionSet(symbol, Eq(f, 0), domain)

    if result is not None:
        if isinstance(result, FiniteSet):
            # this is to simplify solutions like -sqrt(-I) to sqrt(2)/2
            # - sqrt(2)*I/2. We are not expanding for solution with free
            # variables because that makes the solution more complicated. For
            # example expand_complex(a) returns re(a) + I*im(a)
            if all([s.free_symbols == set() and not isinstance(s, RootOf)
                    for s in result]):
                s = Dummy('s')
                result = imageset(Lambda(s, expand_complex(s)), result)
        if isinstance(result, FiniteSet):
            result = result.intersection(domain)
        return result
    else:
        return ConditionSet(symbol, Eq(f, 0), domain)


def _has_rational_power(expr, symbol):
    """
    Returns (bool, den) where bool is True if the term has a
    non-integer rational power and den is the denominator of the
    expression's exponent.

    Examples
    ========

    >>> from sympy.solvers.solveset import _has_rational_power
    >>> from sympy import sqrt
    >>> from sympy.abc import x
    >>> _has_rational_power(sqrt(x), x)
    (True, 2)
    >>> _has_rational_power(x**2, x)
    (False, 1)
    """
    a, p, q = Wild('a'), Wild('p'), Wild('q')
    pattern_match = expr.match(a*p**q) or {}
    if pattern_match.get(a, S.Zero) is S.Zero:
        return (False, S.One)
    elif p not in pattern_match.keys():
        return (False, S.One)
    elif isinstance(pattern_match[q], Rational) \
            and pattern_match[p].has(symbol):
        if not pattern_match[q].q == S.One:
            return (True, pattern_match[q].q)

    if not isinstance(pattern_match[a], Pow) \
            or isinstance(pattern_match[a], Mul):
        return (False, S.One)
    else:
        return _has_rational_power(pattern_match[a], symbol)


def _solve_radical(f, symbol, solveset_solver):
    """ Helper function to solve equations with radicals """
    eq, cov = unrad(f)
    if not cov:
        result = solveset_solver(eq, symbol) - \
            Union(*[solveset_solver(g, symbol) for g in denoms(f, [symbol])])
    else:
        y, yeq = cov
        if not solveset_solver(y - I, y):
            yreal = Dummy('yreal', real=True)
            yeq = yeq.xreplace({y: yreal})
            eq = eq.xreplace({y: yreal})
            y = yreal
        g_y_s = solveset_solver(yeq, symbol)
        f_y_sols = solveset_solver(eq, y)
        result = Union(*[imageset(Lambda(y, g_y), f_y_sols)
                         for g_y in g_y_s])

    return FiniteSet(*[s for s in result if checksol(f, symbol, s) is True])


def _solve_abs(f, symbol, domain):
    """ Helper function to solve equation involving absolute value function """
    if not domain.is_subset(S.Reals):
        raise ValueError(filldedent('''
            Absolute values cannot be inverted in the
            complex domain.'''))
    p, q, r = Wild('p'), Wild('q'), Wild('r')
    pattern_match = f.match(p*Abs(q) + r) or {}
    if not pattern_match.get(p, S.Zero).is_zero:
        f_p, f_q, f_r = pattern_match[p], pattern_match[q], pattern_match[r]
        q_pos_cond = solve_univariate_inequality(f_q >= 0, symbol,
                                                 relational=False)
        q_neg_cond = solve_univariate_inequality(f_q < 0, symbol,
                                                 relational=False)

        sols_q_pos = solveset_real(f_p*f_q + f_r,
                                           symbol).intersect(q_pos_cond)
        sols_q_neg = solveset_real(f_p*(-f_q) + f_r,
                                           symbol).intersect(q_neg_cond)
        return Union(sols_q_pos, sols_q_neg)
    else:
        return ConditionSet(symbol, Eq(f, 0), domain)


def solve_decomposition(f, symbol, domain):
    """
    Function to solve equations via the principle of "Decomposition
    and Rewriting".

    Examples
    ========
    >>> from sympy import exp, sin, Symbol, pprint, S
    >>> from sympy.solvers.solveset import solve_decomposition as sd
    >>> x = Symbol('x')
    >>> f1 = exp(2*x) - 3*exp(x) + 2
    >>> sd(f1, x, S.Reals)
    {0, log(2)}
    >>> f2 = sin(x)**2 + 2*sin(x) + 1
    >>> pprint(sd(f2, x, S.Reals), use_unicode=False)
              3*pi
    {2*n*pi + ---- | n in Integers()}
               2
    >>> f3 = sin(x + 2)
    >>> pprint(sd(f3, x, S.Reals), use_unicode=False)
    {2*n*pi - 2 | n in Integers()} U {pi*(2*n + 1) - 2 | n in Integers()}

    """
    from sympy.solvers.decompogen import decompogen
    from sympy.calculus.util import function_range
    # decompose the given function
    g_s = decompogen(f, symbol)
    # `y_s` represents the set of values for which the function `g` is to be
    # solved.
    # `solutions` represent the solutions of the equations `g = y_s` or
    # `g = 0` depending on the type of `y_s`.
    # As we are interested in solving the equation: f = 0
    y_s = FiniteSet(0)
    for g in g_s:
        frange = function_range(g, symbol, domain)
        y_s = Intersection(frange, y_s)
        result = S.EmptySet
        if isinstance(y_s, FiniteSet):
            for y in y_s:
                solutions = solveset(Eq(g, y), symbol, domain)
                if not isinstance(solutions, ConditionSet):
                    result += solutions

        else:
            if isinstance(y_s, ImageSet):
                iter_iset = (y_s,)

            elif isinstance(y_s, Union):
                iter_iset = y_s.args

            for iset in iter_iset:
                new_solutions = solveset(Eq(iset.lamda.expr, g), symbol, domain)
                dummy_var = tuple(iset.lamda.expr.free_symbols)[0]
                base_set = iset.base_set
                if isinstance(new_solutions, FiniteSet):
                    new_exprs = new_solutions

                elif isinstance(new_solutions, Intersection):
                    if isinstance(new_solutions.args[1], FiniteSet):
                        new_exprs = new_solutions.args[1]

                for new_expr in new_exprs:
                    result += ImageSet(Lambda(dummy_var, new_expr), base_set)

        if result is S.EmptySet:
            return ConditionSet(symbol, Eq(f, 0), domain)

        y_s = result

    return y_s


def _solveset(f, symbol, domain, _check=False):
    """Helper for solveset to return a result from an expression
    that has already been sympify'ed and is known to contain the
    given symbol."""
    # _check controls whether the answer is checked or not

    from sympy.simplify.simplify import signsimp
    orig_f = f
    f = together(f)
    if f.is_Mul:
        _, f = f.as_independent(symbol, as_Add=False)
    if f.is_Add:
        a, h = f.as_independent(symbol)
        m, h = h.as_independent(symbol, as_Add=False)
        f = a/m + h  # XXX condition `m != 0` should be added to soln
    f = piecewise_fold(f)

    # assign the solvers to use
    solver = lambda f, x, domain=domain: _solveset(f, x, domain)
    if domain.is_subset(S.Reals):
        inverter_func = invert_real
    else:
        inverter_func = invert_complex
    inverter = lambda f, rhs, symbol: inverter_func(f, rhs, symbol, domain)

    result = EmptySet()

    if f.expand().is_zero:
        return domain
    elif not f.has(symbol):
        return EmptySet()
    elif f.is_Mul and all(_is_finite_with_finite_vars(m, domain)
            for m in f.args):
        # if f(x) and g(x) are both finite we can say that the solution of
        # f(x)*g(x) == 0 is same as Union(f(x) == 0, g(x) == 0) is not true in
        # general. g(x) can grow to infinitely large for the values where
        # f(x) == 0. To be sure that we are not silently allowing any
        # wrong solutions we are using this technique only if both f and g are
        # finite for a finite input.
        result = Union(*[solver(m, symbol) for m in f.args])
    elif _is_function_class_equation(TrigonometricFunction, f, symbol) or \
            _is_function_class_equation(HyperbolicFunction, f, symbol):
        result = _solve_trig(f, symbol, domain)
    elif f.is_Piecewise:
        dom = domain
        result = EmptySet()
        expr_set_pairs = f.as_expr_set_pairs()
        for (expr, in_set) in expr_set_pairs:
            if in_set.is_Relational:
                in_set = in_set.as_set()
            if in_set.is_Interval:
                dom -= in_set
            solns = solver(expr, symbol, in_set)
            result += solns
    else:
        lhs, rhs_s = inverter(f, 0, symbol)
        if lhs == symbol:
            # do some very minimal simplification since
            # repeated inversion may have left the result
            # in a state that other solvers (e.g. poly)
            # would have simplified; this is done here
            # rather than in the inverter since here it
            # is only done once whereas there it would
            # be repeated for each step of the inversion
            if isinstance(rhs_s, FiniteSet):
                rhs_s = FiniteSet(*[Mul(*
                    signsimp(i).as_content_primitive())
                    for i in rhs_s])
            result = rhs_s
        elif isinstance(rhs_s, FiniteSet):
            for equation in [lhs - rhs for rhs in rhs_s]:
                if equation == f:
                    if any(_has_rational_power(g, symbol)[0]
                           for g in equation.args) or _has_rational_power(
                           equation, symbol)[0]:
                        result += _solve_radical(equation,
                                                 symbol,
                                                 solver)
                    elif equation.has(Abs):
                        result += _solve_abs(f, symbol, domain)
                    else:
                        result += _solve_as_rational(equation, symbol, domain)
                else:
                    result += solver(equation, symbol)
        else:
            result = ConditionSet(symbol, Eq(f, 0), domain)

    if _check:
        if isinstance(result, ConditionSet):
            # it wasn't solved or has enumerated all conditions
            # -- leave it alone
            return result

        # whittle away all but the symbol-containing core
        # to use this for testing
        fx = orig_f.as_independent(symbol, as_Add=True)[1]
        fx = fx.as_independent(symbol, as_Add=False)[1]

        if isinstance(result, FiniteSet):
            # check the result for invalid solutions
            result = FiniteSet(*[s for s in result
                      if isinstance(s, RootOf)
                      or domain_check(fx, symbol, s)])

    return result


def solveset(f, symbol=None, domain=S.Complexes):
    """Solves a given inequality or equation with set as output

    Parameters
    ==========

    f : Expr or a relational.
        The target equation or inequality
    symbol : Symbol
        The variable for which the equation is solved
    domain : Set
        The domain over which the equation is solved

    Returns
    =======

    Set
        A set of values for `symbol` for which `f` is True or is equal to
        zero. An `EmptySet` is returned if `f` is False or nonzero.
        A `ConditionSet` is returned as unsolved object if algorithms
        to evaluatee complete solution are not yet implemented.

    `solveset` claims to be complete in the solution set that it returns.

    Raises
    ======

    NotImplementedError
        The algorithms to solve inequalities in complex domain  are
        not yet implemented.
    ValueError
        The input is not valid.
    RuntimeError
        It is a bug, please report to the github issue tracker.


    Notes
    =====

    Python interprets 0 and 1 as False and True, respectively, but
    in this function they refer to solutions of an expression. So 0 and 1
    return the Domain and EmptySet, respectively, while True and False
    return the opposite (as they are assumed to be solutions of relational
    expressions).


    See Also
    ========

    solveset_real: solver for real domain
    solveset_complex: solver for complex domain

    Examples
    ========

    >>> from sympy import exp, sin, Symbol, pprint, S
    >>> from sympy.solvers.solveset import solveset, solveset_real

    * The default domain is complex. Not specifying a domain will lead
      to the solving of the equation in the complex domain (and this
      is not affected by the assumptions on the symbol):

    >>> x = Symbol('x')
    >>> pprint(solveset(exp(x) - 1, x), use_unicode=False)
    {2*n*I*pi | n in Integers()}

    >>> x = Symbol('x', real=True)
    >>> pprint(solveset(exp(x) - 1, x), use_unicode=False)
    {2*n*I*pi | n in Integers()}

    * If you want to use `solveset` to solve the equation in the
      real domain, provide a real domain. (Using `solveset\_real`
      does this automatically.)

    >>> R = S.Reals
    >>> x = Symbol('x')
    >>> solveset(exp(x) - 1, x, R)
    {0}
    >>> solveset_real(exp(x) - 1, x)
    {0}

    The solution is mostly unaffected by assumptions on the symbol,
    but there may be some slight difference:

    >>> pprint(solveset(sin(x)/x,x), use_unicode=False)
    ({2*n*pi | n in Integers()} \ {0}) U ({2*n*pi + pi | n in Integers()} \ {0})

    >>> p = Symbol('p', positive=True)
    >>> pprint(solveset(sin(p)/p, p), use_unicode=False)
    {2*n*pi | n in Integers()} U {2*n*pi + pi | n in Integers()}

    * Inequalities can be solved over the real domain only. Use of a complex
      domain leads to a NotImplementedError.

    >>> solveset(exp(x) > 1, x, R)
    (0, oo)

    """
    f = sympify(f)

    if f is S.true:
        return domain

    if f is S.false:
        return S.EmptySet

    if not isinstance(f, (Expr, Number)):
        raise ValueError("%s is not a valid SymPy expression" % (f))

    free_symbols = f.free_symbols

    if not free_symbols:
        b = Eq(f, 0)
        if b is S.true:
            return domain
        elif b is S.false:
            return S.EmptySet
        else:
            raise NotImplementedError(filldedent('''
                relationship between value and 0 is unknown: %s''' % b))

    if symbol is None:
        if len(free_symbols) == 1:
            symbol = free_symbols.pop()
        else:
            raise ValueError(filldedent('''
                The independent variable must be specified for a
                multivariate equation.'''))
    elif not getattr(symbol, 'is_Symbol', False):
        raise ValueError('A Symbol must be given, not type %s: %s' %
            (type(symbol), symbol))

    if isinstance(f, Eq):
        from sympy.core import Add
        f = Add(f.lhs, - f.rhs, evaluate=False)
    elif f.is_Relational:
        if not domain.is_subset(S.Reals):
            raise NotImplementedError(filldedent('''
                Inequalities in the complex domain are
                not supported. Try the real domain by
                setting domain=S.Reals'''))
        try:
            result = solve_univariate_inequality(
            f, symbol, relational=False) - _invalid_solutions(
            f, symbol, domain)
        except NotImplementedError:
            result = ConditionSet(symbol, f, domain)
        return result

    return _solveset(f, symbol, domain, _check=True)


def _invalid_solutions(f, symbol, domain):
    bad = S.EmptySet
    for d in denoms(f):
        bad += _solveset(d, symbol, domain, _check=False)
    return bad


def solveset_real(f, symbol):
    return solveset(f, symbol, S.Reals)


def solveset_complex(f, symbol):
    return solveset(f, symbol, S.Complexes)


###############################################################################
################################ LINSOLVE #####################################
###############################################################################


def linear_eq_to_matrix(equations, *symbols):
    r"""
    Converts a given System of Equations into Matrix form.
    Here `equations` must be a linear system of equations in
    `symbols`. The order of symbols in input `symbols` will
    determine the order of coefficients in the returned
    Matrix.

    The Matrix form corresponds to the augmented matrix form.
    For example:

    .. math:: 4x + 2y + 3z  = 1
    .. math:: 3x +  y +  z  = -6
    .. math:: 2x + 4y + 9z  = 2

    This system would return `A` & `b` as given below:

    ::

         [ 4  2  3 ]          [ 1 ]
     A = [ 3  1  1 ]   b  =   [-6 ]
         [ 2  4  9 ]          [ 2 ]

    Examples
    ========

    >>> from sympy import linear_eq_to_matrix, symbols
    >>> x, y, z = symbols('x, y, z')
    >>> eqns = [x + 2*y + 3*z - 1, 3*x + y + z + 6, 2*x + 4*y + 9*z - 2]
    >>> A, b = linear_eq_to_matrix(eqns, [x, y, z])
    >>> A
    Matrix([
    [1, 2, 3],
    [3, 1, 1],
    [2, 4, 9]])
    >>> b
    Matrix([
    [ 1],
    [-6],
    [ 2]])
    >>> eqns = [x + z - 1, y + z, x - y]
    >>> A, b = linear_eq_to_matrix(eqns, [x, y, z])
    >>> A
    Matrix([
    [1,  0, 1],
    [0,  1, 1],
    [1, -1, 0]])
    >>> b
    Matrix([
    [1],
    [0],
    [0]])

    * Symbolic coefficients are also supported

    >>> a, b, c, d, e, f = symbols('a, b, c, d, e, f')
    >>> eqns = [a*x + b*y - c, d*x + e*y - f]
    >>> A, B = linear_eq_to_matrix(eqns, x, y)
    >>> A
    Matrix([
    [a, b],
    [d, e]])
    >>> B
    Matrix([
    [c],
    [f]])

    """

    if not symbols:
        raise ValueError('Symbols must be given, for which coefficients \
                         are to be found.')

    if hasattr(symbols[0], '__iter__'):
        symbols = symbols[0]

    M = Matrix([symbols])
    # initialise Matrix with symbols + 1 columns
    M = M.col_insert(len(symbols), Matrix([1]))
    row_no = 1

    for equation in equations:
        f = sympify(equation)
        if isinstance(f, Equality):
            f = f.lhs - f.rhs

        # Extract coeff of symbols
        coeff_list = []
        for symbol in symbols:
            coeff_list.append(f.coeff(symbol))

        # append constant term (term free from symbols)
        coeff_list.append(-f.as_coeff_add(*symbols)[0])

        # insert equations coeff's into rows
        M = M.row_insert(row_no, Matrix([coeff_list]))
        row_no += 1

    # delete the initialised (Ist) trivial row
    M.row_del(0)
    A, b = M[:, :-1], M[:, -1:]
    return A, b


def linsolve(system, *symbols):
    r"""
    Solve system of N linear equations with M variables, which
    means both under - and overdetermined systems are supported.
    The possible number of solutions is zero, one or infinite.
    Zero solutions throws a ValueError, where as infinite
    solutions are represented parametrically in terms of given
    symbols. For unique solution a FiniteSet of ordered tuple
    is returned.

    All Standard input formats are supported:
    For the given set of Equations, the respective input types
    are given below:

    .. math:: 3x + 2y -   z = 1
    .. math:: 2x - 2y + 4z = -2
    .. math:: 2x -   y + 2z = 0

    * Augmented Matrix Form, `system` given below:

    ::

              [3   2  -1  1]
     system = [2  -2   4 -2]
              [2  -1   2  0]

    * List Of Equations Form

    `system  =  [3x + 2y - z - 1, 2x - 2y + 4z + 2, 2x - y + 2z]`

    * Input A & b Matrix Form (from Ax = b) are given as below:

    ::

         [3   2  -1 ]         [  1 ]
     A = [2  -2   4 ]    b =  [ -2 ]
         [2  -1   2 ]         [  0 ]

    `system = (A, b)`

    Symbols to solve for should be given as input in all the
    cases either in an iterable or as comma separated arguments.
    This is done to maintain consistency in returning solutions
    in the form of variable input by the user.

    The algorithm used here is Gauss-Jordan elimination, which
    results, after elimination, in an row echelon form matrix.

    Returns
    =======

    A FiniteSet of ordered tuple of values of `symbols` for which
    the `system` has solution.

    Please note that general FiniteSet is unordered, the solution
    returned here is not simply a FiniteSet of solutions, rather
    it is a FiniteSet of ordered tuple, i.e. the first & only
    argument to FiniteSet is a tuple of solutions, which is ordered,
    & hence the returned solution is ordered.

    Also note that solution could also have been returned as an
    ordered tuple, FiniteSet is just a wrapper `{}` around
    the tuple. It has no other significance except for
    the fact it is just used to maintain a consistent output
    format throughout the solveset.

    Returns EmptySet(), if the linear system is inconsistent.

    Raises
    ======

    ValueError
        The input is not valid.
        The symbols are not given.

    Examples
    ========

    >>> from sympy import Matrix, S, linsolve, symbols
    >>> x, y, z = symbols("x, y, z")
    >>> A = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
    >>> b = Matrix([3, 6, 9])
    >>> A
    Matrix([
    [1, 2,  3],
    [4, 5,  6],
    [7, 8, 10]])
    >>> b
    Matrix([
    [3],
    [6],
    [9]])
    >>> linsolve((A, b), [x, y, z])
    {(-1, 2, 0)}

    * Parametric Solution: In case the system is under determined, the function
      will return parametric solution in terms of the given symbols.
      Free symbols in the system are returned as it is. For e.g. in the system
      below, `z` is returned as the solution for variable z, which means z is a
      free symbol, i.e. it can take arbitrary values.

    >>> A = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> b = Matrix([3, 6, 9])
    >>> linsolve((A, b), [x, y, z])
    {(z - 1, -2*z + 2, z)}

    * List of Equations as input

    >>> Eqns = [3*x + 2*y - z - 1, 2*x - 2*y + 4*z + 2, - x + S(1)/2*y - z]
    >>> linsolve(Eqns, x, y, z)
    {(1, -2, -2)}

    * Augmented Matrix as input

    >>> aug = Matrix([[2, 1, 3, 1], [2, 6, 8, 3], [6, 8, 18, 5]])
    >>> aug
    Matrix([
    [2, 1,  3, 1],
    [2, 6,  8, 3],
    [6, 8, 18, 5]])
    >>> linsolve(aug, x, y, z)
    {(3/10, 2/5, 0)}

    * Solve for symbolic coefficients

    >>> a, b, c, d, e, f = symbols('a, b, c, d, e, f')
    >>> eqns = [a*x + b*y - c, d*x + e*y - f]
    >>> linsolve(eqns, x, y)
    {((-b*f + c*e)/(a*e - b*d), (a*f - c*d)/(a*e - b*d))}

    * A degenerate system returns solution as set of given
      symbols.

    >>> system = Matrix(([0,0,0], [0,0,0], [0,0,0]))
    >>> linsolve(system, x, y)
    {(x, y)}

    * For an empty system linsolve returns empty set

    >>> linsolve([ ], x)
    EmptySet()

    """

    if not system:
        return S.EmptySet

    if not symbols:
        raise ValueError('Symbols must be given, for which solution of the '
                         'system is to be found.')

    if hasattr(symbols[0], '__iter__'):
        symbols = symbols[0]

    try:
        sym = symbols[0].is_Symbol
    except AttributeError:
        sym = False

    if not sym:
        raise ValueError('Symbols or iterable of symbols must be given as '
                         'second argument, not type %s: %s' % (type(symbols[0]), symbols[0]))

    # 1). Augmented Matrix input Form
    if isinstance(system, Matrix):
        A, b = system[:, :-1], system[:, -1:]

    elif hasattr(system, '__iter__'):

        # 2). A & b as input Form
        if len(system) == 2 and system[0].is_Matrix:
            A, b = system[0], system[1]

        # 3). List of equations Form
        if not system[0].is_Matrix:
            A, b = linear_eq_to_matrix(system, symbols)

    else:
        raise ValueError("Invalid arguments")

    # Solve using Gauss-Jordan elimination
    try:
        sol, params, free_syms = A.gauss_jordan_solve(b, freevar=True)
    except ValueError:
        # No solution
        return EmptySet()

    # Replace free parameters with free symbols
    solution = []
    if params:
        for s in sol:
            for k, v in enumerate(params):
                s = s.xreplace({v: symbols[free_syms[k]]})
            solution.append(simplify(s))
    else:
        for s in sol:
            solution.append(simplify(s))

    # Return solutions
    solution = FiniteSet(tuple(solution))
    return solution


###############################################################################
################################ NLINSOLVE ####################################
###############################################################################


def substitution(system, symbols, result, known_symbols, all_symbols):
    r"""
    Solves the `system` using substitution method.

    Parameters
    ==========

    system : list of equations

    symbols : list of unsolved symbols.

    result : If non of the symbols are solved then it is empty list, otherwise
    list of already solved_syms symbols( dict symbol : value )

    known_symbols : list of symbols already solved (might be in terms of other
    symbols, that will be solved).

    all_symbols : known_symbols + unsolved symbols.

    Returns
    =======
    Returns Finiteset .

    Examples
    ========

    >>> from sympy.core.symbol import symbols
    >>> x, y = symbols('x, y', real = True)
    >>> from sympy.solvers.solveset import substitution
    >>> substitution([x +y], [x], [{y : 1}], [y], [x, y])
    {{x: -1, y: 1}}

    """
    # TODO: known_symbol is not needed, remove this variable.
    # It is equal to keys of result.
    from sympy.core.compatibility import ordered, default_sort_key
    from sympy import Complement
    from sympy.core.containers import Dict
    # sort so equation with the fewest potential symbols is first

    def _ok_syms(e, sort=False):
            rv = (e.free_symbols - set(known_symbols)) & set(all_symbols)
            if sort:
                rv = list(rv)
                rv.sort(key=default_sort_key)
            return rv

    for eq in ordered(system, lambda _: len(_ok_syms(_))):
        u = Dummy()  # used in solution checking
        newresult = [] # TODO: Make it FiniteSet
        bad_results = []
        got_s = set()
        hit = False
        complements = {}
        for r in result:
            # update eq with everything that is known so far
            eq2 = eq.subs(r)
            # if check is True then we see if it satisfies this
            # equation, otherwise we just accept it
            if r:
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
                not_solvable = False
                soln_imageset = None
                try:
                    soln = solveset_real(eq2, s)
                    # Not sure to add the complex solution or not
                    # soln = soln + solveset_complex(eq2, s)
                    if not soln:
                        soln = solveset_complex(eq2, s)
                except NotImplementedError:
                    continue
                # put each solution in r and append the now-expanded
                # result in the new result list; use copy since the
                # solution for s in being added in-place
                if isinstance(soln, ImageSet):
                    soln_imageset = soln
                    soln = FiniteSet(soln.lamda.expr)
                elif isinstance(soln, ConditionSet):
                    soln = FiniteSet()
                    not_solvable = True
                elif isinstance(soln, Complement):
                    # extract solution and complement
                    complements[s] = list(soln.args[1])[0]
                    soln = soln.args[0]
                    # complement should be added at the end
                elif isinstance(soln, Intersection):
                    # sometimes solveset returns Intersection with S.Real
                    soln = soln.args[1]
                for sol in soln:
                    if got_s and any([ss in sol.free_symbols for ss in got_s]):
                        # sol depends on previously solved symbols: discard it
                        continue
                    rnew = dict(r.copy())
                    for k, v in r.items():
                        if isinstance(v, Expr):
                            # if any unsolved symbol is present
                            # Then subs known value
                            rnew[k] = v.subs(s, sol)
                    # and add this new solution
                    rnew[s] = sol if not soln_imageset else soln_imageset
                    newresult.append(rnew)
                hit = True
                if not not_solvable:
                    got_s.add(s)
            if not hit:
                raise NotImplementedError('could not solve %s' % eq2)
        else:
            result = newresult
            for b in bad_results:
                if b in result:
                    result.remove(b)

    result_finiteset = FiniteSet()
    infinite_soln = 0
    # If soln have general soln and some finite soln,
    # then general soln return(when infinte_soln == 1).
    for r in result:
        if not r:
            # if {None : None} is present.
            # No solution then first will be {None : None}
            return S.EmptySet

        # If length < len(all_symbols) means infinite soln.
        # Some or all the soln is dependent on 1 symbol.
        # eg. {x: y+2} then final soln is {x: y+2, y: y}
        if len(r) < len(all_symbols):
            solved_symbols = r.keys()
            unsolved = list(filter(lambda x: x not in solved_symbols,
                                           all_symbols))
            rcopy = dict(r.copy())
            # If `r` is SymPy Dict. Convert it to python dict
            for us in unsolved:
                rcopy[us] = us
            # if symbol is not in `v` , means it can take any value.
            # eg. if k, v => y, exp(x) then above lines will add {x: x}
            # if we have another symbol `z` then add {z: z} using above line.
            r = rcopy
            infinite_soln += 1
        result_finiteset = result_finiteset +  FiniteSet(Dict(r))

    # if infinte_soln == 1 means we have general soln
    # eg : {{x: -1, y : 1}, {x : -y , y: y}} then
    # return {{x : -y, y : y}} only which is last element always
    result_finiteset = result_finiteset \
    if not infinite_soln == 1 else FiniteSet(list(result_finiteset)[-1])

    if complements:
        # If solveset have returned some complements for any symbol
        result = FiniteSet()
        for res in result_finiteset:
            res_copy = dict(res)
            for k_res, v_res in res.items():
                for k_c, v_c in complements.items():
                    if k_c == k_res:
                        res_copy[k_res] = FiniteSet(v_res) - FiniteSet(v_c)
            result = result + FiniteSet(Dict(res_copy))
        result_finiteset = result
    return result_finiteset

def nlinsolve(system, symbols):
    r"""
    Solve system of N non linear equations with M variables, which
    means both under - and overdetermined systems are supported.
    Positive dimensional system is also supported (Infinite solution).
    in Positive dimensional system solution will be dependent on at
    least one symbol.
    If system doesn't have real solution then complex solution will
    be returned in general form ( in terms of `_n`, where `_n`
    is any real number).
    The possible number of solutions is zero, one or infinite.

    Parameters
    ==========

    system : list of equations
        The target system of equations
    symbols : list of Symbol

    * Note : List of symbol is used here.

    Returns
    =======

    A FiniteSet of ordered tuple of values of `symbols` for which
    the `system` has solution, when system is zero dimensional
    system. For positive dimensional system A Finiteset of Dict
    ( key = symbol and value = symbol solution). Dict is defined
    at sympy.core.containers .

    Please note that general FiniteSet is unordered, the solution
    returned here is not simply a FiniteSet of solutions, rather
    it is a FiniteSet of ordered tuple, i.e. the first & only
    argument to FiniteSet is a tuple of solutions, which is ordered,
    & hence the returned solution is ordered.

    Also note that solution could also have been returned as an
    ordered tuple, FiniteSet is just a wrapper `{}` around
    the tuple. It has no other significance except for
    the fact it is just used to maintain a consistent output
    format throughout the solveset.

    For the given set of Equations, the respective input types
    are given below:

    .. math:: x*y - 1 = 0
    .. math:: 4*x**2 + y**2 - 5 = 0

    `system  = [x*y - 1, 4*x**2 + y**2 - 5]`
    `symbols = [x, y]`

    >>> from sympy.core.symbol import symbols
    >>> from sympy.solvers.solveset import nlinsolve
    >>> x, y = symbols('x, y', real = True)
    >>> nlinsolve([x*y - 1, 4*x**2 + y**2 - 5], [x, y])
    {(-1, -1), (-1/2, -2), (1/2, 2), (1, 1)}


    * Positive dimensional system :

    Examples
    ========

    * If solveset (substitution method) returns complement
    for any symbol in then that will also be considered in the
    final solution. Following example is that type.

    >>> from sympy.core.symbol import symbols
    >>> from sympy import pprint
    >>> from sympy.polys.polytools import is_zero_dimensional
    >>> from sympy.solvers.solveset import nlinsolve
    >>> a, b, c, d = symbols('a, b, c, d', real = True)
    >>> foo =  a + b + c + d
    >>> bar = a*b + b*c + c*d + d*a
    >>> foo_bar = a*b*c + b*c*d + c*d*a + d*a*b
    >>> bar_foo = a*b*c*d -1
    >>> system = [foo, bar, foo_bar, bar_foo]
    >>> is_zero_dimensional(system)
    False
    >>> pprint(nlinsolve(system, [a, b, c, d]))
         -1             1
    {{a: ---, b: -d, c: -, d: {d} \ {0}}}
          d             d


    >>> x, y = symbols('x, y', real = True)
    >>> nlinsolve([(x+y)**2 - 4, x + y - 2], [x, y])
    {{x: -y + 2, y: y}}

    * Note: You have to take assumption `real = True`, otherwise solveset returns
    solution with Intersection S.Reals.

    * Non linear system having non polynomial equation, if there is
    complex solution for particular symbol and no real solution, then `nlinsolve`
    returns it's general solution.

    Example
    =======

    >>> from sympy.core.symbol import symbols
    >>> from sympy import sqrt, exp, sin
    >>> from sympy.solvers.solveset import nlinsolve
    >>> x, y, _n = symbols('x, y, _n')
    >>> nlinsolve([exp(x) - sin(y), y**2 - 4], [x, y])
    {{x: log(sin(2)), y: 2}, {x: ImageSet(Lambda(_n, I*(2*_n*pi + pi) + log(sin(2))), Integers()), y: -2}}


    * Non linear system having all the equations polynomial, then it
    returns both real and complex solutions.

    Example
    =======

    >>> from sympy.core.symbol import symbols
    >>> from sympy import sqrt, exp, sin
    >>> from sympy.solvers.solveset import nlinsolve
    >>> x, y = symbols('x, y')
    >>> nlinsolve([x**2 - 2*y**2 -2, x*y - 2], [x, y])
    {(-2, -1), (2, 1), (-sqrt(2)*I, sqrt(2)*I), (sqrt(2)*I, -sqrt(2)*I)}


    * If system is linear positive dimensional system, `nlinsolve` can
    solve this system also, since it is using `groebner` method.

    Examples
    ========

    >>> from sympy.core.symbol import symbols
    >>> from sympy import pprint
    >>> from sympy.solvers.solveset import nlinsolve
    >>> x, y, z = symbols('x, y, z', real = True)
    >>> pprint(nlinsolve([x + 2*y -z - 3, x - y - 4*z + 9 , y + z - 4], [x, y, z]))
    {{x: 3*z - 5, y: -z + 4, z: z}}


    More examples:
    =============

    >>> from sympy.core.symbol import symbols
    >>> from sympy import sqrt
    >>> from sympy.solvers.solveset import nlinsolve
    >>> x, y, z = symbols('x, y, z', real = True)
    >>> e1 = sqrt(x**2 + y**2) - 10
    >>> e2 = sqrt(y**2 + (-x + 10)**2) - 3
    >>> nlinsolve((e1, e2), (x, y))
    {(191/20, -3*sqrt(391)/20), (191/20, 3*sqrt(391)/20)}
    >>> nlinsolve([x**2 + 2/y - 2, x + y - 3], [x, y])
    {{x: 1, y: 2}, {x: 1 + sqrt(5), y: -sqrt(5) + 2}, {x: -sqrt(5) + 1, y: 2 + sqrt(5)}}

    The last example is zero dimensional but is_zero_dimensional returned false
    thats why solution is coming from substitution method

    Note :
    =======

    1. If system if zero dimensional system (Finite solution,
    solvable using `solve_poly_system`) then if symbols = [y, z, x]
    then solution is Finiteset((y_solution, z_solution, x_solution)),
    means in the same order.
    2. If infinite solution (substitution method is used ), then solution
    will be in Finiteset Dict and always ordered . If symbols = [y, z, x]
    then final solution = Finiteset({x : x_solution, y : y_solution, z :
        z_solution})
    3. `Dict` is defined in `SymPy`. `Finiteset` can't be used with Python `dict`.

    """
    from sympy.solvers.solvers import _invert as _invert_solver
    from sympy.utilities.iterables import subsets
    from sympy.polys.polytools import is_zero_dimensional, groebner
    from sympy.core.containers import Dict

    if not system:
        return S.EmptySet

    if not symbols:
        raise ValueError('Symbols must be given, for which solution of the '
                         'system is to be found.')
    polys = []
    nonpolys = []
    for j, g in enumerate(system):
        # TODO : solveset `_invert, improve for more than one symbols
        # move all the terms, having any `symbols` in lhs and make it 'g',
        # currently using old solver's _invert
        i, d = _invert_solver(g, *symbols)
        g = d - i
        g = g.as_numer_denom()[0]

        poly = g.as_poly(*symbols, extension=True)
        if poly is not None:
            polys.append(poly)
        else:
            nonpolys.append(g)

    # If none of the equation is Poly
    if not polys:
        solved_syms = []

    result = None
    if len(symbols) == len(polys):
        if is_zero_dimensional(system):
            try:
                # try to solve, when all the equations are poly
                result = solve_poly_system(polys, *symbols)
                return FiniteSet(*[s for s in result])
            except NotImplementedError:
                # Right now We don't know the failed case
                pass
        else:
            # positive dimensional system
            # Do substitution method with groebner basis of the system
            basis = groebner(polys, symbols, polys=True)
            new_system = []
            for p in basis:
                new_system.append(p.as_expr())
            # solved_symbols = []
            result = [{}]
            result = substitution(new_system, symbols, result, [], symbols)

    else:
        result = FiniteSet()
        # all the equations are not Polynomial
        depend_soln = {}
        solved_syms = []
        # first solve the polynomial equations if present
        if polys:
            combinations = list(subsets(symbols, len(polys)))
            for new_symbols in combinations:
                try:
                    # solution for new_symbols in terms of other symbols
                    res = solve_poly_system(polys, new_symbols)
                    for r in res:
                        skip = False
                        # check tuples of res
                        for r1 in r:
                            if depend_soln and any([k in r1.free_symbols
                                for k, v in depend_soln.items()]):
                                # sol depends on previously
                                # solved symbols: discard it
                                # eg : depend_soln=> {x : -y} and
                                # r1 is function of x like r1 = -x for y
                                skip = True
                        if not skip:
                            for ns in new_symbols:
                                for r1 in r:
                                    depend_soln[ns] = r1
                            # need to maintain dict <symbol : value>
                            dlist = dict(list(zip(new_symbols, r)))
                            result += FiniteSet(Dict(dlist))
                            # Finiteset accept the sympy Dict
                    result_update = FiniteSet()
                    for r in result:
                        # If length < len(symbols) means infinite soln.
                        # Some or all the soln is dependent on 1 symbol
                        if len(r) < len(symbols):
                            unsolved = None
                            rcopy = dict(r.copy())
                            # if SymPy Dict then convert it to Python dict
                            # so that we can update using index
                            for k, v in r.items():
                                if isinstance(v, Expr):
                                    unsolved = list(filter(lambda x: x in
                                                           v.free_symbols,
                                                           symbols))
                                    # there will be only one symbol in v if
                                    # unsolved in not Empty
                                    if unsolved:
                                        unsolved = unsolved[0]
                            if unsolved:
                                rcopy[unsolved] = unsolved
                            result_update += FiniteSet(Dict(rcopy))
                        else:
                            result_update += FiniteSet(Dict(r))
                    result = result_update
                except NotImplementedError:
                    pass
            if depend_soln:
                # TODO no need of solved_syms. Remove it from entire method
                solved_syms = list(depend_soln)
        # Polynomail is done. Now use substition method to get the solution
        # for unsolved symbols if nonpolys list is not None.
        if nonpolys:
            if not result:
                result = [{}]
            # if non polynomail equation is present
            unsolved_syms = list(filter(lambda x: x not in solved_syms, symbols))
            # one by one solve for unsolved after substitution of solved symbols values.
            result = substitution(nonpolys, unsolved_syms, list(result), solved_syms, symbols)

    if result:
        return result
    else:
        return S.EmptySet
