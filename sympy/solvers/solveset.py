"""
This module contains functions to solve a single equation for a single variable.
"""
from __future__ import print_function, division

from sympy.core.sympify import sympify
from sympy.core import S, Pow, Dummy, pi, Expr, Wild, Mul, Equality
from sympy.core.numbers import I, Number, Rational
from sympy.core.function import (Lambda, expand, expand_complex)
from sympy.core.relational import Eq
from sympy.simplify.simplify import fraction, trigsimp
from sympy.functions import (log, Abs, tan, cot, exp,
                             arg, Piecewise, piecewise_fold)
from sympy.functions.elementary.trigonometric import (TrigonometricFunction,
                                                      HyperbolicFunction)
from sympy.sets import FiniteSet, EmptySet, imageset, Union
from sympy.matrices import Matrix
from sympy.polys import (roots, Poly, degree, together, PolynomialError,
                         RootOf)
from sympy.solvers.solvers import checksol, denoms
from sympy.utilities import filldedent

import warnings


def invert_real(f_x, y, x):
    """ Inverts a real valued function

    Reduces the real valued equation ``f(x) = y`` to a set of equations ``{g(x)
    = h_1(y), g(x) = h_2(y), ..., g(x) = h_n(y) }`` where ``g(x)`` is a simpler
    function than ``f(x)``.  The return value is a tuple ``(g(x), set_h)``,
    where ``g(x)`` is a function of ``x`` and ``set_h`` is the set of
    functions ``{h_1(y), h_2(y), ..., h_n(y)}``.
    Here, ``y`` is not necessarily a symbol.

    Examples
    ========

    >>> from sympy.solvers.solveset import invert_real
    >>> from sympy import tan, Abs, exp
    >>> from sympy.abc import x, y, n
    >>> invert_real(exp(Abs(x)), y, x)
    (x, {-log(y), log(y)})
    >>> invert_real(exp(x), 1, x)
    (x, {0})
    >>> invert_real(Abs(x**31 + x), y, x)
    (x**31 + x, {-y, y})
    >>> invert_real(tan(x), y, x)
    (x, ImageSet(Lambda(_n, _n*pi + atan(y)), Integers()))

    See Also
    ========
    invert_complex
    """
    y = sympify(y)
    if not y.has(x):
        return _invert_real(f_x, FiniteSet(y), x)
    else:
        raise ValueError(" y should be independent of x ")


def _invert_real(f, g_ys, symbol):
    """ Helper function for invert_real """

    if not f.has(symbol):
        raise ValueError("Inverse of constant function doesn't exist")

    if f is symbol:
        return (f, g_ys)

    n = Dummy('n')
    if hasattr(f, 'inverse') and not isinstance(f, TrigonometricFunction):
        if len(f.args) > 1:
            raise ValueError("Only functions with one argument are supported.")
        return _invert_real(f.args[0],
                            imageset(Lambda(n, f.inverse()(n)), g_ys), symbol)

    if isinstance(f, Abs):
        return _invert_real(f.args[0],
                            Union(g_ys, imageset(Lambda(n, -n), g_ys)), symbol)

    if f.is_Add:
        # f = g + h
        g, h = f.as_independent(symbol)
        if g != S.Zero:
            return _invert_real(h, imageset(Lambda(n, n - g), g_ys), symbol)

    if f.is_Mul:
        # f = g*h
        g, h = f.as_independent(symbol)

        if g != S.One:
            return _invert_real(h, imageset(Lambda(n, n/g), g_ys), symbol)

    if f.is_Pow:
        base, expo = f.args
        base_has_sym = base.has(symbol)
        expo_has_sym = expo.has(symbol)

        if not expo_has_sym:
            res = imageset(Lambda(n, Pow(n, 1/expo)), g_ys)
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
            return _invert_real(expo, imageset(Lambda(n, log(n)/log(base)),
                                               g_ys), symbol)

    if isinstance(f, tan) or isinstance(f, cot):
        n = Dummy('n')
        if isinstance(g_ys, FiniteSet):
            tan_cot_invs = Union(*[imageset(Lambda(n, n*pi + f.inverse()(g_y)),
                                            S.Integers) for g_y in g_ys])
            return _invert_real(f.args[0], tan_cot_invs, symbol)

    return (f, g_ys)


def invert_complex(f_x, y, x):
    """ Inverts a complex valued function.

    Reduces the complex valued equation ``f(x) = y`` to a set of equations
    ``{g(x) = h_1(y), g(x) = h_2(y), ..., g(x) = h_n(y) }`` where ``g(x)`` is
    a simpler function than ``f(x)``.  The return value is a tuple ``(g(x),
    set_h)``, where ``g(x)`` is a function of ``x`` and ``set_h`` is
    the set of function ``{h_1(y), h_2(y), ..., h_n(y)}``.
    Here, ``y`` is not necessarily a symbol.

    Note that `invert\_complex` and `invert\_real` don't always produce the
    same result even for a seemingly simple function like ``exp(x)`` because
    the complex extension of real valued ``log`` is multivariate in the complex
    system and has infinitely many branches. If you are working with real
    values only or you are not sure with function to use you should use
    `invert\_real`.


    Examples
    ========

    >>> from sympy.solvers.solveset import invert_complex
    >>> from sympy.abc import x, y
    >>> from sympy import exp, log
    >>> invert_complex(log(x), y, x)
    (x, {exp(y)})
    >>> invert_complex(log(x), 0, x)  # Second parameter is not a symbol
    (x, {1})
    >>> invert_complex(exp(x), y, x)
    (x, ImageSet(Lambda(_n, I*(2*_n*pi + arg(y)) + log(Abs(y))), Integers()))

    See Also
    ========
    invert_real
    """
    y = sympify(y)
    if not y.has(x):
        return _invert_complex(f_x, FiniteSet(y), x)
    else:
        raise ValueError(" y should be independent of x ")


def _invert_complex(f, g_ys, symbol):
    """ Helper function for invert_complex """

    if not f.has(symbol):
        raise ValueError("Inverse of constant function doesn't exist")

    if f is symbol:
        return (f, g_ys)

    n = Dummy('n')
    if f.is_Add:
        # f = g + h
        g, h = f.as_independent(symbol)
        if g != S.Zero:
            return _invert_complex(h, imageset(Lambda(n, n - g), g_ys), symbol)

    if f.is_Mul:
        # f = g*h
        g, h = f.as_independent(symbol)

        if g != S.One:
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
                               for g_y in g_ys])
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


def _is_finite_with_finite_vars(f):
    """
    Return True if the given expression is finite when all free symbols
    (that are not already specified as finite) are made finite.
    """
    reps = dict([(s, Dummy(s.name, finite=True, **s.assumptions0))
                for s in f.free_symbols if s.is_finite is None])
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


def solveset_real(f, symbol):
    """ Solves a real valued equation.

    Parameters
    ==========

    f : Expr
        The target equation
    symbol : Symbol
        The variable for which the equation is solved

    Returns
    =======

    Set
        A set of values for `symbol` for which `f` is equal to
        zero. An `EmptySet` is returned if no solution is found.

    `solveset_real` claims to be complete in the set of the solution it
    returns.

    Raises
    ======

    NotImplementedError
        The algorithms for to find the solution of the given equation are
        not yet implemented.
    ValueError
        The input is not valid.
    RuntimeError
        It is a bug, please report to the github issue tracker.


    See Also
    =======

    solveset_complex : solver for complex domain

    Examples
    ========

    >>> from sympy import Symbol, exp, sin, sqrt, I
    >>> from sympy.solvers.solveset import solveset_real
    >>> x = Symbol('x', real=True)
    >>> a = Symbol('a', real=True, finite=True, positive=True)
    >>> solveset_real(x**2 - 1, x)
    {-1, 1}
    >>> solveset_real(sqrt(5*x + 6) - 2 - x, x)
    {-1, 2}
    >>> solveset_real(x - I, x)
    EmptySet()
    >>> solveset_real(x - a, x)
    {a}
    >>> solveset_real(exp(x) - a, x)
    {log(a)}

    * In case the equation has infinitely many solutions an infinitely indexed
      `ImageSet` is returned.

    >>> solveset_real(sin(x) - 1, x)
    ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers())

    * If the equation is true for any arbitrary value of the symbol a `S.Reals`
      set is returned.

    >>> solveset_real(x - x, x)
    (-oo, oo)

    """
    if not symbol.is_Symbol:
        raise ValueError(" %s is not a symbol" % (symbol))

    f = sympify(f)
    if not isinstance(f, (Expr, Number)):
        raise ValueError(" %s is not a valid sympy expression" % (f))

    original_eq = f
    f = together(f)

    if f.has(Piecewise):
        f = piecewise_fold(f)
    result = EmptySet()

    if f.expand().is_zero:
        return S.Reals
    elif not f.has(symbol):
        return EmptySet()
    elif f.is_Mul and all([_is_finite_with_finite_vars(m) for m in f.args]):
        # if f(x) and g(x) are both finite we can say that the solution of
        # f(x)*g(x) == 0 is same as Union(f(x) == 0, g(x) == 0) is not true in
        # general. g(x) can grow to infinitely large for the values where
        # f(x) == 0. To be sure that we not are silently allowing any
        # wrong solutions we are using this technique only if both f and g and
        # finite for a finite input.
        result = Union(*[solveset_real(m, symbol) for m in f.args])
    elif _is_function_class_equation(TrigonometricFunction, f, symbol) or \
            _is_function_class_equation(HyperbolicFunction, f, symbol):
        result = _solve_real_trig(f, symbol)
    elif f.is_Piecewise:
        result = EmptySet()
        expr_set_pairs = f.as_expr_set_pairs()
        for (expr, in_set) in expr_set_pairs:
            solns = solveset_real(expr, symbol).intersect(in_set)
            result = result + solns
    else:
        lhs, rhs_s = invert_real(f, 0, symbol)
        if lhs == symbol:
            result = rhs_s
        elif isinstance(rhs_s, FiniteSet):
            equations = [lhs - rhs for rhs in rhs_s]
            for equation in equations:
                if equation == f:
                    if any(_has_rational_power(g, symbol)[0]
                           for g in equation.args):
                        result += _solve_radical(equation,
                                                 symbol,
                                                 solveset_real)
                    elif equation.has(Abs):
                        result += _solve_abs(f, symbol)
                    else:
                        result += _solve_as_rational(equation, symbol,
                                                     solveset_solver=solveset_real,
                                                     as_poly_solver=_solve_as_poly_real)
                else:
                    result += solveset_real(equation, symbol)
        else:
            raise NotImplementedError

    if isinstance(result, FiniteSet):
        result = [s for s in result
                  if isinstance(s, RootOf)
                  or domain_check(original_eq, symbol, s)]
        return FiniteSet(*result).intersect(S.Reals)
    else:
        return result.intersect(S.Reals)


def _solve_as_rational(f, symbol, solveset_solver, as_poly_solver):
    """ solve rational functions"""
    f = together(f, deep=True)
    g, h = fraction(f)
    if not h.has(symbol):
        return as_poly_solver(g, symbol)
    else:
        valid_solns = solveset_solver(g, symbol)
        invalid_solns = solveset_solver(h, symbol)
        return valid_solns - invalid_solns


def _solve_real_trig(f, symbol):
    """ Helper to solve trigonometric equations """
    f = trigsimp(f)
    f = f.rewrite(exp)
    f = together(f)
    g, h = fraction(f)
    y = Dummy('y')
    g, h = g.expand(), h.expand()
    g, h = g.subs(exp(I*symbol), y), h.subs(exp(I*symbol), y)
    if g.has(symbol) or h.has(symbol):
        raise NotImplementedError

    solns = solveset_complex(g, y) - solveset_complex(h, y)

    if isinstance(solns, FiniteSet):
        return Union(*[invert_complex(exp(I*symbol), s, symbol)[1]
                       for s in solns])
    elif solns is S.EmptySet:
        return S.EmptySet
    else:
        raise NotImplementedError


def _solve_as_poly(f, symbol, solveset_solver, invert_func):
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
                raise NotImplementedError("Couldn't find all roots "
                                          "of the equation %s" % f)
    else:
        poly = Poly(f)
        if poly is None:
            raise NotImplementedError("Could not convert %s to Poly" % f)
        gens = [g for g in poly.gens if g.has(symbol)]

        if len(gens) == 1:
            poly = Poly(poly, gens[0])
            gen = poly.gen
            deg = poly.degree()
            poly = Poly(poly.as_expr(), poly.gen, composite=True)
            poly_solns = FiniteSet(*roots(poly, cubics=True, quartics=True,
                                          quintics=True).keys())

            if len(poly_solns) < deg:
                raise NotImplementedError("Couldn't find all the roots of "
                                          "the equation %s" % f)

            if gen != symbol:
                y = Dummy('y')
                lhs, rhs_s = invert_func(gen, y, symbol)
                if lhs is symbol:
                    result = Union(*[rhs_s.subs(y, s) for s in poly_solns])
                else:
                    raise NotImplementedError(
                        "inversion of %s not handled" % gen)
        else:
            raise NotImplementedError("multiple generators not handled"
                                      " by solveset")

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
        return result
    else:
        raise NotImplementedError


def _solve_as_poly_real(f, symbol):
    """
    Solve real valued equation with methods to solve polynomial
    equations.
    """
    return _solve_as_poly(f, symbol,
                          solveset_solver=solveset_real,
                          invert_func=invert_real)


def _solve_as_poly_complex(f, symbol):
    """
    Solve complex valued equation with methods to solve polynomial
    equations.
    """
    return _solve_as_poly(f, symbol,
                          solveset_solver=solveset_complex,
                          invert_func=invert_complex)


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
    pattern_match = expr.match(a*p**q)
    if pattern_match is None or pattern_match[a] is S.Zero:
        return (False, S.One)
    elif p not in pattern_match.keys() or a not in pattern_match.keys():
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
    from sympy.solvers.solvers import unrad
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


def _solve_abs(f, symbol):
    """ Helper function to solve equation involving absolute value function """
    from sympy.solvers.inequalities import solve_univariate_inequality
    assert f.has(Abs)
    p, q, r = Wild('p'), Wild('q'), Wild('r')
    pattern_match = f.match(p*Abs(q) + r)
    if not pattern_match[p].is_zero:
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
        raise NotImplementedError


def solveset_complex(f, symbol):
    """ Solve a complex valued equation.

    Parameters
    ==========

    f : Expr
        The target equation
    symbol : Symbol
        The variable for which the equation is solved

    Returns
    =======

    Set
        A set of values for `symbol` for which `f` equal to
        zero. An `EmptySet` is returned if no solution is found.

    `solveset_complex` claims to be complete in the solution set that
    it returns.

    Raises
    ======

    NotImplementedError
        The algorithms for to find the solution of the given equation are
        not yet implemented.
    ValueError
        The input is not valid.
    RuntimeError
        It is a bug, please report to the github issue tracker.

    See Also
    ========

    solveset_real: solver for real domain

    Examples
    ========

    >>> from sympy import Symbol, exp
    >>> from sympy.solvers.solveset import solveset_complex
    >>> from sympy.abc import x, a, b, c
    >>> solveset_complex(a*x**2 + b*x +c, x)
    {-b/(2*a) - sqrt(-4*a*c + b**2)/(2*a), -b/(2*a) + sqrt(-4*a*c + b**2)/(2*a)}

    * Due to the fact that complex extension of my real valued functions are
      multivariate even some simple equations can have infinitely many
      solution.

    >>> solveset_complex(exp(x) - 1, x)
    ImageSet(Lambda(_n, 2*_n*I*pi), Integers())

    """
    if not symbol.is_Symbol:
        raise ValueError(" %s is not a symbol" % (symbol))

    f = sympify(f)
    original_eq = f
    if not isinstance(f, (Expr, Number)):
        raise ValueError(" %s is not a valid sympy expression" % (f))

    f = together(f)
    # Without this equations like a + 4*x**2 - E keep oscillating
    # into form  a/4 + x**2 - E/4 and (a + 4*x**2 - E)/4
    if not fraction(f)[1].has(symbol):
        f = expand(f)

    if f.is_zero:
        raise NotImplementedError("S.Complex set is not yet implemented")
    elif not f.has(symbol):
        result = EmptySet()
    elif f.is_Mul and all([_is_finite_with_finite_vars(m) for m in f.args]):
        result = Union(*[solveset_complex(m, symbol) for m in f.args])
    else:
        lhs, rhs_s = invert_complex(f, 0, symbol)
        if lhs == symbol:
            result = rhs_s
        elif isinstance(rhs_s, FiniteSet):
            equations = [lhs - rhs for rhs in rhs_s]
            result = EmptySet()
            for equation in equations:
                if equation == f:
                    if any(_has_rational_power(g, symbol)[0]
                           for g in equation.args):
                        result += _solve_radical(equation,
                                                 symbol,
                                                 solveset_complex)
                    else:
                        result += _solve_as_rational(equation, symbol,
                                                 solveset_solver=solveset_complex,
                                                 as_poly_solver=_solve_as_poly_complex)
                else:
                    result += solveset_complex(equation, symbol)
        else:
            raise NotImplementedError

    if isinstance(result, FiniteSet):
        result = [s for s in result
                  if isinstance(s, RootOf)
                  or domain_check(original_eq, symbol, s)]
        return FiniteSet(*result)
    else:
        return result


def solveset(f, symbol=None):
    """Solves a given inequality or equation with set as output

    Parameters
    ==========

    f : Expr or a relational.
        The target equation or inequality
    symbol : Symbol
        The variable for which the equation is solved

    Returns
    =======

    Set
        A set of values for `symbol` for which `f` is True or is equal to
        zero. An `EmptySet` is returned if no solution is found.

    `solveset` claims to be complete in the solution set that it returns.

    Raises
    ======

    NotImplementedError
        The algorithms for to find the solution of the given equation are
        not yet implemented.
    ValueError
        The input is not valid.
    RuntimeError
        It is a bug, please report to the github issue tracker.


    `solveset` uses two underlying functions `solveset_real` and
    `solveset_complex` to solve equations. They are
    the solvers for real and complex domain respectively. The domain of
    the solver is decided by the assumption on the variable for which the
    equation is being solved.


    See Also
    ========

    solveset_real: solver for real domain
    solveset_complex: solver for complex domain

    Examples
    ========

    >>> from sympy import exp, Symbol, Eq, pprint
    >>> from sympy.solvers.solveset import solveset
    >>> from sympy.abc import x

    * Symbols in Sympy are complex by default. A complex variable
      will lead to the solving of the equation in complex domain.

    >>> pprint(solveset(exp(x) - 1, x), use_unicode=False)
    {2*n*I*pi | n in Integers()}

    * If you want to solve equation in real domain by the `solveset`
      interface, then specify the variable to real. Alternatively use
      `solveset\_real`.

    >>> x = Symbol('x', real=True)
    >>> solveset(exp(x) - 1, x)
    {0}
    >>> solveset(Eq(exp(x), 1), x)
    {0}

    * Inequalities are always solved in the real domain irrespective of
      the assumption on the variable for which the inequality is solved.

    >>> solveset(exp(x) > 1, x)
    (0, oo)

    """

    from sympy.solvers.inequalities import solve_univariate_inequality

    if symbol is None:
        free_symbols = f.free_symbols
        if len(free_symbols) == 1:
            symbol = free_symbols.pop()
        else:
            raise ValueError(filldedent('''
                The independent variable must be specified for a
                multivariate equation.'''))
    elif not symbol.is_Symbol:
        raise ValueError('A Symbol must be given, not type %s: %s' % (type(symbol), symbol))

    real = (symbol.is_real is True)

    f = sympify(f)

    if isinstance(f, Eq):
        from sympy.core import Add
        f = Add(f.lhs, - f.rhs, evaluate=False)

    if f.is_Relational:
        if real is False:
            warnings.warn(filldedent('''
                The variable you are solving for is complex
                but will assumed to be real since solving complex
                inequalities is not supported.
            '''))
        return solve_univariate_inequality(f, symbol, relational=False)

    if isinstance(f, (Expr, Number)):
        if real is True:
            return solveset_real(f, symbol)
        else:
            return solveset_complex(f, symbol)


###############################################################################
################################ LINSOLVE #####################################
###############################################################################


def linear_eq_to_matrix(equations, *symbols):
    """
    Converts a given System of Equations into Matrix form.
    Here `equations` must be a linear system of equations in
    `symbols`. The order of symbols in input `symbols` will
    determine the order of coefficients in the augmented
    Matrix.

    The Matrix form corresponds to the augmented matrix form.
    For example:

      x + 2.y + 3.z  = 1
    3.x +   y +   z  = -6
    2.x + 4.y + 9.z  = 2

    This system would return A & b:

          [ 1  2  3 ]         [ 1 ]
    A  =  [ 3  1  1 ]    b =  [-6 ]
          [ 2  4  9 ]         [ 2 ]


    Examples
    ========

    >>> from sympy.solvers.solveset import linear_eq_to_matrix
    >>> from sympy import symbols
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

    """

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
