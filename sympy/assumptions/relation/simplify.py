"""
Module to define simplify logic for equality and inequalities.
"""

from sympy.core import Add, Expr, S
from sympy.core.add import _unevaluated_Add
from sympy.core.compatibility import ordered
from sympy.polys import Poly, poly, PolynomialError, gcd
from sympy.solvers.solveset import linear_coeffs


def relsimp(func, args, **kwargs):
    lhs, rhs = args

    r = func(lhs.simplify(**kwargs), rhs.simplify(**kwargs))
    if not isinstance(r.lhs, Expr) or not isinstance(r.rhs, Expr):
        return r.canonical
    dif = r.lhs - r.rhs
    r = r.canonical
    # If there is only one symbol in the expression,
    # try to write it on a simplified form
    free = list(filter(lambda x: x.is_real is not False, r.free_symbols))
    if len(free) == 1:
        try:
            x = free.pop()
            dif = r.lhs - r.rhs
            m, b = linear_coeffs(dif, x)
            if m.is_zero is False:
                if m.is_negative:
                    # Dividing with a negative number, so change order of arguments
                    # canonical will put the symbol back on the lhs later
                    r = r.function(-b / m, x)
                else:
                    r = r.function(x, -b / m)
            else:
                r = r.function(b, S.Zero)
        except ValueError:
            # maybe not a linear function, try polynomial
            try:
                p = poly(dif, x)
                c = p.all_coeffs()
                constant = c[-1]
                c[-1] = 0
                scale = gcd(c)
                c = [ctmp / scale for ctmp in c]
                r = r.function(Poly.from_list(c, x).as_expr(), -constant / scale)
            except PolynomialError:
                pass
    elif len(free) >= 2:
        try:
            free = list(ordered(free))
            dif = r.lhs - r.rhs
            m = linear_coeffs(dif, *free)
            constant = m[-1]
            del m[-1]
            scale = gcd(m)
            m = [mtmp / scale for mtmp in m]
            nzm = list(filter(lambda f: f[0] != 0, list(zip(m, free))))
            if scale.is_zero is False:
                if constant != 0:
                    # lhs: expression, rhs: constant
                    newexpr = Add(*[i * j for i, j in nzm])
                    r = r.function(newexpr, -constant / scale)
                else:
                    # keep first term on lhs
                    lhsterm = nzm[0][0] * nzm[0][1]
                    del nzm[0]
                    newexpr = Add(*[i * j for i, j in nzm])
                    r = r.function(lhsterm, -newexpr)

            else:
                r = r.function(constant, S.Zero)
        except ValueError:
            pass
    # Did we get a simplified result?
    r = r.canonical
    rel = func(*args)
    measure = kwargs['measure']
    if measure(r) < kwargs['ratio'] * measure(rel):
        return r
    else:
        return rel


def eqsimp(eq, **kwargs):
    # standard simplify
    eq = relsimp(eq.function, eq.arguments, **kwargs)

    if not isinstance(eq.lhs, Expr) or not isinstance(eq.rhs, Expr):
        return eq
    free = eq.free_symbols
    if len(free) == 1:
        try:
            x = free.pop()
            m, b = linear_coeffs(
                _convert_to_Add(eq, evaluate=False), x)
            if m.is_zero is False:
                eqnew = eq.function(x, -b / m)
            else:
                eqnew = eq.function(m * x, -b)
            measure = kwargs['measure']
            if measure(eqnew) <= kwargs['ratio'] * measure(eq):
                eq = eqnew
        except ValueError:
            pass
    return eq.canonical


def _convert_to_Add(eq, evaluate=True):
    # Mimic Equality._eval_rewrite_as_Add.
    args = eq.lhs, eq.rhs
    L, R = args
    if evaluate:
        # allow cancellation of args
        return L - R
    args = Add.make_args(L) + Add.make_args(-R)
    if evaluate is None:
        # no cancellation, but canonical
        return _unevaluated_Add(*args)
    # no cancellation, not canonical
    return Add._from_args(args)
