from sympy.solvers import solve
from sympy.simplify import simplify
from sympy.solvers.solveset import solveset_real
from sympy.sets.sets import Interval, Intersection, FiniteSet, Union, Complement
from sympy import sympify, diff, limit, S, oo


def singularities(expr, sym):
    """
    Finds singularities for a function.
    Currently supported functions are:
    - univariate real rational functions

    Examples
    ========

    >>> from sympy.calculus.singularities import singularities
    >>> from sympy import Symbol
    >>> x = Symbol('x', real=True)
    >>> singularities(x**2 + x + 1, x)
    ()
    >>> singularities(1/(x + 1), x)
    (-1,)

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Mathematical_singularity

    """
    if not expr.is_rational_function(sym):
        raise NotImplementedError("Algorithms finding singularities for"
                                  " non rational functions are not yet"
                                  " implemented")
    else:
        return tuple(sorted(solve(simplify(1/expr), sym)))


def range_func(func, set_value):
    func = sympify(func)

    if set_value.is_EmptySet:
        return EmptySet()

    free_symbol = func.free_symbols
    if len(free_symbol) == 1:
        symbol = free_symbol.pop()
    elif len(free_symbol) == 0:
        return FiniteSet(func)
    else:
        raise NotImplementedError("more than one variables involved")

    def sing_handle(f, set_val, val):
        val1 = limit(f, symbol, val)
        val2 = limit(f, symbol, val, '-')
        if val in set_val.boundary:
            if val is set_val.start:
                return val1
            else:
                return val2

    # all the singularities of the function
    sing = singularities(func, symbol)

    # this block of code can be replaced by
    # sing = Intersection(FiniteSet(*sing), set_value.closure)
    # after the issue #9706 has been fixed
    if set_value.has(Interval):
        if not oo in set_value.boundary:
            if not S.NegativeInfinity in set_value.boundary:
                sing = Intersection(FiniteSet(*sing), set_value.closure)
            else:
                sing = Intersection(FiniteSet(*sing), Union(set_value, FiniteSet(set_value.end)))
        else:
            if not S.NegativeInfinity in set_value.boundary:
                sing = Intersection(FiniteSet(*sing), Union(set_value, FiniteSet(set_value.start)))
            else:
                sing = Intersection(FiniteSet(*sing), set_value)

    def in_intrvl(f, set_val):
        val1 = (set_val.start, set_val.args[2])
        val2 = (set_val.end, set_val.args[3])
        g1 = simplify(diff(f, symbol))
        expr1 = g1 > 0
        expr2 = g1 < 0
        der_zero = solveset_real(g1, symbol)
        if not (FiniteSet(val1[0], val2[0]).contains(-oo) or FiniteSet(val1[0], val2[0]).contains(oo)):
            der_zero = Intersection(der_zero, set_val.closure)
        g2 = simplify(diff(g1, symbol))
        maxi = set()
        mini = set()

        ans1 = limit(f, symbol, val1[0])
        ans2 = limit(f, symbol, val2[0], '-')
        singl = solveset_real(1/f, symbol)

        # this too can be closed after fixing #9706

        for i in singl:
            if i in sing and not i in set_val.boundary:
                return Union(range_func(f, Interval(val1[0], i, val1[1], True)),
                            range_func(f, Interval(i, val2[0], True, val2[1])))
        if ans1 is S.Infinity or ans2 is S.Infinity:
            maxi = set([(oo, True)])
        elif ans1 is S.NegativeInfinity or ans2 is S.NegativeInfinity:
            mini = set([(-oo, True)])
        if maxi == set():
            more = max(ans1, ans2)
            if ans1 > ans2:
                maxi = set([(ans1, val1[1])])
            elif ans1 < ans2:
                maxi = set([(ans2, val2[1])])
            else:
                maxi = set([(ans1, val1[1] and val2[1])])
        if mini == set():
            if ans1 < ans2:
                mini = set([(ans1, val1[1])])
            elif ans1 > ans2:
                mini = set([(ans2, val2[1])])
            else:
                mini = set([(ans1, val1[1] and val2[1])])
        unk = set()

        for i in der_zero:
            exist = not i in set_val
            if g2.subs({symbol: i}) < 0:
                if not i in sing:
                    maxi.add((f.subs({symbol: i}), exist))
                else:
                    maxi.add((oo, True))
            elif g2.subs({symbol: i}) > 0:
                if not i in sing:
                    mini.add((f.subs({symbol: i}), exist))
                else:
                    mini.add((-oo, True))
            else:
                unk.add(f.subs({symbol: i}))
        ma = (-oo, True)
        mi = (oo, True)

        for i in maxi:
            if i[0] > ma[0]:
                ma = i
            elif i[0] == ma[0]:
                ma = (ma[0], i[1] and ma[1])

        for i in mini:
            if i[0] < mi[0]:
                mi = i
            elif i[0] == mi[0]:
                mi = (mi[0], i[1] and mi[1])

        return Union(Interval(mi[0], ma[0], mi[1], ma[1]), FiniteSet(*unk))

    if isinstance(set_value, Union):
        return Union(*[range_func(func, se) for se in set_value.args])

    if isinstance(set_value, Interval):
        return in_intrvl(func, set_value)

    if isinstance(set_value, FiniteSet):
        set_value = Complement(set_value, FiniteSet(*sing))
        return FiniteSet(*[limit(func, symbol, i) if not i in sing else func.subs({symbol: i}) for i in set_value])
