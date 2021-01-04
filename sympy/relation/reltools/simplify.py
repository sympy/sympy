"""
Module for equation simplification.
"""

from functools import singledispatch

from sympy.core import Add, Expr, S
from sympy.core.compatibility import ordered
from sympy.polys import Poly, poly, PolynomialError, gcd
from sympy.solvers.solveset import linear_coeffs

from sympy.relation.binrel import BinaryRelation
from sympy.relation.equality import Equal


@singledispatch
def eqnsimp(rel, lhs, rhs, **kwargs):
    raise NotImplementedError("Unknown binary predicate.")


@eqnsimp.register(BinaryRelation)
def _(rel, lhs, rhs, **kwargs):
    eqn = rel(lhs, rhs)

    r = rel(lhs.simplify(**kwargs), rhs.simplify(**kwargs))
    if not isinstance(r.lhs, Expr) or not isinstance(r.rhs, Expr):
        return r
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
    measure = kwargs['measure']
    if measure(r) < kwargs['ratio'] * measure(eqn):
        return r
    else:
        return eqn


@eqnsimp.register(Equal)
def _(rel, lhs, rhs, **kwargs):
    eqn = rel(lhs, rhs)

    # standard simplify
    e = eqnsimp.dispatch(BinaryRelation)(rel, lhs, rhs, **kwargs)

    if not isinstance(e.function, Equal):
        return e
    if not isinstance(e.lhs, Expr) or not isinstance(e.rhs, Expr):
        return e
    free = eqn.free_symbols
    if len(free) == 1:
        try:
            x = free.pop()
            m, b = linear_coeffs(
                _convert_to_Add(e, evaluate=False), x)
            if m.is_zero is False:
                enew = e.function(x, -b / m)
            else:
                enew = e.function(m * x, -b)
            measure = kwargs['measure']
            if measure(enew) <= kwargs['ratio'] * measure(e):
                e = enew
        except ValueError:
            pass
    return e.canonical


def _convert_to_Add(rel, evaluate=True):
    # Mimic Relational._eval_rewrite_as_Add.
    # Binary predicate does not have that method because it does not
    # make sense. Boolean cannot be rewritten to other expression.

    # Perhaps this can be just replaced to `rel-rel.rhs` with no evaluation?
    from sympy.core.add import _unevaluated_Add
    args = rel.lhs, rel.rhs
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
