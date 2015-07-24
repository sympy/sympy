from sympy import sympify, diff, limit, S, oo
from sympy.core import Expr
from sympy.core.numbers import Number
from sympy.calculus.singularities import singularities
from sympy.sets.sets import Interval, Intersection, FiniteSet, Union, Complement, Set
from sympy.simplify import simplify
from sympy.solvers.solveset import solveset_real


def range_func(func, set_value):
    """ Finds the range of a real-valued function

    Parameters
    ==========

    func: Expr
          The expression whose range is to be found
    set_value: Union of Sets
          The domain for the variable involved in function

    Raises
    ======

    NotImplementedError
          The algorithms for to find the solution of the given equation are
          not yet implemented.
    ValueError
          The input is not valid.
    RuntimeError
          It is a bug, please report to the github issue tracker.

    Examples
    ========

    >>> from sympy.calculus.range_func import range_func
    >>> from sympy import Symbol, S, Interval, Union, FiniteSet
    >>> x = Symbol('x', real=True)
    >>> range_func(x**2, Interval(-1, 1))
    [0, 1]
    >>> range_func(x/(x**2 - 4), Union(Interval(-1, 3), FiniteSet(5)))
    (-oo, 1/3] U [3/5, oo)
    >>> range_func(x**2/(x**2 - 4), S.Reals)
    (-oo, 0] U (1, oo)
    """

    func = sympify(func)
    if not isinstance(func, (Expr, Number)):
        raise ValueError(" %s is not a valid sympy expression" % (func))

    if not isinstance(set_value, Set):
        raise ValueError('A Symbol must be given, not type %s: %s' % (type(set_value), set_value))

    if set_value.is_EmptySet:
        return EmptySet()

    free_symbol = func.free_symbols
    if len(free_symbol) == 1:
        symbol = free_symbol.pop()
    elif len(free_symbol) == 0:
        return FiniteSet(func)
    else:
        raise NotImplementedError("more than one variables %s not handeled" % (free_symbol))

    if not func.is_rational_function(symbol):
        raise NotImplementedError("Algorithms finding range for non-rational functions"
                                    "are not yet implemented")

    # this block of code can be replaced by
    # sing = Intersection(FiniteSet(*sing), set_val.closure)
    # after the issue #9706 has been fixed
    def closure_handle(set_val, sing):
        if set_value.has(Interval):
            if not oo in set_value.boundary:
                if not S.NegativeInfinity in set_value.boundary:
                    return Intersection(FiniteSet(*sing), set_value.closure)
                return Intersection(FiniteSet(*sing), Union(set_value,
                                    FiniteSet(max(set_value.boundary))))
            else:
                if not S.NegativeInfinity in set_value.boundary:
                    return Intersection(FiniteSet(*sing), Union(set_value,
                                        FiniteSet(min(set_value.boundary))))
                return Intersection(FiniteSet(*sing), set_value)
        return Intersection(FiniteSet(*sing), set_value)

    # all the singularities of the function
    sing = singularities(func, symbol)
    sing = closure_handle(set_value, sing)

    def in_intrvl(f, set_val):
        val1 = (set_val.start, set_val.args[2])
        val2 = (set_val.end, set_val.args[3])
        g1 = simplify(diff(f, symbol))
        g2 = simplify(diff(g1, symbol))
        expr1 = g1 > 0
        der_zero = solveset_real(g1, symbol)
        der_zero = closure_handle(set_val, der_zero)

        maxi = set()
        mini = set()
        ans1 = limit(f, symbol, val1[0])
        ans2 = limit(f, symbol, val2[0], '-')
        singl = solveset_real(1/f, symbol)

        for i in singl:
            if i in sing and not i in set_val.boundary:
                return Union(range_func(f, Interval(val1[0], i, val1[1], True)),
                            range_func(f, Interval(i, val2[0], True, val2[1])))

        if ans1 is S.Infinity or ans2 is S.Infinity:
            maxi = set([(oo, True)])
        elif ans1 is S.NegativeInfinity or ans2 is S.NegativeInfinity:
            mini = set([(-oo, True)])
        if maxi == set():
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
        return FiniteSet(*[limit(func, symbol, i) if i in FiniteSet(-oo, oo)
                            else func.subs({symbol: i}) for i in set_value])
