from sympy.multipledispatch import dispatch, Dispatcher
from sympy.core import Basic, Expr, Function, Add, Mul, Pow, Dummy, Integer
from sympy import Min, Max, Set, sympify, Lambda, symbols, exp, log, S
from sympy.sets import (imageset, Interval, FiniteSet, Union, ImageSet,
    ProductSet, EmptySet, Intersection)
from sympy.core.function import FunctionClass
from sympy.logic.boolalg import And, Or, Not, true, false


_x, _y = symbols("x y")


@dispatch(Set, Set)
def add_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x+_y)), x, y)


@dispatch(Expr, Expr)
def add_sets(x, y):
    return x+y


@dispatch(Interval, Interval)
def add_sets(x, y):
    """
    Additions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    return Interval(x.start + y.start, x.end + y.end,
        x.left_open or y.left_open, x.right_open or y.right_open)


@dispatch(Expr, Expr)
def sub_sets(x, y):
    return x-y


@dispatch(Set, Set)
def sub_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x - _y)), x, y)


@dispatch(Interval, Interval)
def sub_sets(x, y):
    """
    Subtractions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    return Interval(x.start - y.end, x.end - y.start,
        x.left_open or y.right_open, x.right_open or y.left_open)


@dispatch(Set, Set)
def mul_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x * _y)), x, y)


@dispatch(Expr, Expr)
def mul_sets(x, y):
    return x*y


@dispatch(Interval, Interval)
def mul_sets(x, y):
    """
    Multiplications in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    comvals = (
        (x.start * y.start, bool(x.left_open or y.left_open)),
        (x.start * y.end, bool(x.left_open or y.right_open)),
        (x.end * y.start, bool(x.right_open or y.left_open)),
        (x.end * y.end, bool(x.right_open or y.right_open)),
    )
    # TODO: handle symbolic intervals
    minval, minopen = min(comvals)
    maxval, maxopen = max(comvals)
    return Interval(
        minval,
        maxval,
        minopen,
        maxopen
    )
    return SetExpr(Interval(start, end))


@dispatch(Expr, Expr)
def div_sets(x, y):
    return x/y


@dispatch(Set, Set)
def div_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x / _y)), x, y)


@dispatch(Interval, Interval)
def div_sets(x, y):
    """
    Divisions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    if (y.start*y.end).is_negative:
        from sympy import oo
        return Interval(-oo, oo)
    return mul_sets(x, Interval(1/y.end, 1/y.start, y.right_open, y.left_open))


@dispatch(Set, Set)
def pow_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x ** _y)), x, y)


@dispatch(Expr, Expr)
def pow_sets(x, y):
    return x**y


@dispatch(Interval, Integer)
def pow_sets(x, y):
    """
    Powers in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    exponent = sympify(exponent)
    if exponent.is_odd:
        return Interval(x.start**exponent, x.end**exponent, x.left_open, x.right_open)
    if exponent.is_even:
        if (x.start*x.end).is_negative:
            if -x.start > x.end:
                left_limit = x.start
                left_open = x.right_open
            else:
                left_limit = x.end
                left_open = x.left_open
            return Interval(S.Zero, left_limit ** exponent, S.Zero not in x, left_open)
        elif x.start.is_negative and x.end.is_negative:
            return Interval(x.end**exponent, x.start**exponent, x.right_open, x.left_open)
        else:
            return Interval(x.start**exponent, x.end**exponent, x.left_open, x.right_open)


FunctionUnion = (FunctionClass, Lambda)


@dispatch(FunctionUnion, FiniteSet)
def function_sets(f, x):
    return FiniteSet(*map(f, x))

@dispatch(Lambda, Interval)
def function_sets(f, x):
    from sympy.functions.elementary.miscellaneous import Min, Max
    from sympy.solvers.solveset import solveset
    from sympy.core.function import diff, Lambda
    from sympy.series import limit
    from sympy.calculus.singularities import singularities
    from sympy.sets import Complement
    # TODO: handle functions with infinitely many solutions (eg, sin, tan)
    # TODO: handle multivariate functions

    expr = f.expr
    if len(expr.free_symbols) > 1 or len(f.variables) != 1:
        return
    var = f.variables[0]

    if expr.is_Piecewise:
        result = S.EmptySet
        domain_set = x
        for (p_expr, p_cond) in expr.args:
            if p_cond is true:
                intrvl = domain_set
            else:
                intrvl = p_cond.as_set()
                intrvl = Intersection(domain_set, intrvl)

            if p_expr.is_Number:
                image = FiniteSet(p_expr)
            else:
                image = imageset(Lambda(var, p_expr), intrvl)
            result = Union(result, image)

            # remove the part which has been `imaged`
            domain_set = Complement(domain_set, intrvl)
            if domain_set.is_EmptySet:
                break
        return result

    if not x.start.is_comparable or not x.end.is_comparable:
        return

    try:
        sing = [i for i in singularities(expr, var)
            if i.is_real and i in x]
    except NotImplementedError:
        return

    if x.left_open:
        _start = limit(expr, var, x.start, dir="+")
    elif x.start not in sing:
        _start = f(x.start)
    if x.right_open:
        _end = limit(expr, var, x.end, dir="-")
    elif x.end not in sing:
        _end = f(x.end)

    if len(sing) == 0:
        solns = list(solveset(diff(expr, var), var))

        extr = [_start, _end] + [f(i) for i in solns
                                 if i.is_real and i in x]
        start, end = Min(*extr), Max(*extr)

        left_open, right_open = False, False
        if _start <= _end:
            # the minimum or maximum value can occur simultaneously
            # on both the edge of the interval and in some interior
            # point
            if start == _start and start not in solns:
                left_open = x.left_open
            if end == _end and end not in solns:
                right_open = x.right_open
        else:
            if start == _end and start not in solns:
                left_open = x.right_open
            if end == _start and end not in solns:
                right_open = x.left_open

        return Interval(start, end, left_open, right_open)
    else:
        return imageset(f, Interval(x.start, sing[0],
                                    x.left_open, True)) + \
            Union(*[imageset(f, Interval(sing[i], sing[i + 1], True, True))
                    for i in range(0, len(sing) - 1)]) + \
            imageset(f, Interval(sing[-1], x.end, True, x.right_open))

@dispatch(FunctionClass, Interval)
def function_sets(f, x):
    if f == exp:
        return Interval(exp(x.start), exp(x.end), x.left_open, x.right_open)
    elif f == log:
        return Interval(log(x.start), log(x.end), x.left_open, x.right_open)
    return ImageSet(Lambda(_x, f(_x)), x)

@dispatch(FunctionUnion, Union)
def function_sets(f, x):
    return Union(imageset(f, arg) for arg in x.args)

@dispatch(FunctionUnion, Intersection)
def function_sets(f, x):
    # If the function is invertible, intersect the maps of the sets.
    u = symbols("u")
    fdiff = f(u).diff(u)
    # TODO: find a better condition for invertible functions:
    if ((f in (exp, log))  # functions known to be invertible
        or (fdiff > 0) == True or (fdiff < 0) == True  # monotonous funcs
        ):
        return Intersection(imageset(f, arg) for arg in x.args)
    else:
        return ImageSet(Lambda(_x, f(_x)), x)

@dispatch(FunctionUnion, EmptySet)
def function_sets(f, x):
    return x

@dispatch(FunctionUnion, Set)
def function_sets(f, x):
    return ImageSet(Lambda(_x, f(_x)), x)
