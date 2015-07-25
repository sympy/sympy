from sympy import S, sympify, diff, limit, oo
from sympy.calculus.singularities import singularities
from sympy.sets.sets import Interval, Intersection, FiniteSet, Union, Complement, Set
from sympy.solvers.solveset import solveset_real


def codomain(func, domain, *syms):
    """ Finds the range of a real-valued function, for a real-domain

    Parameters
    ==========

    func: Expr
          The expression whose range is to be found
    domain: Union of Sets
          The real-domain for the variable involved in function
    syms: Tuple of symbols
          Symbol whose domain is given

    Raises
    ======

    NotImplementedError
          The algorithms to find the range of the given function are
          not yet implemented.
    ValueError
          The input is not valid.
    RuntimeError
          It is a bug, please report to the github issue tracker.

    Examples
    ========

    >>> from sympy import Symbol, S, Interval, Union, FiniteSet, codomain
    >>> x = Symbol('x', real=True)
    >>> codomain(x**2, Interval(-1, 1), x)
    [0, 1]
    >>> codomain(x/(x**2 - 4), Union(Interval(-1, 3), FiniteSet(5)), x)
    (-oo, 1/3] U [3/5, oo)
    >>> codomain(x**2/(x**2 - 4), S.Reals, x)
    (-oo, 0] U (1, oo)
    """

    func = sympify(func)
    if not isinstance(domain, Set):
        raise ValueError('A Set must be given, not %s: %s' % (type(domain), domain))

    if len(syms) == 0:
        raise ValueError("A Symbol or a tuple of symbols must be given: not %s" % (type(syms)))

    if len(syms) == 1:
        symbol = syms[0]
    else:
        raise NotImplementedError("more than one variables %s not handeled" % (syms))

    if not func.is_rational_function(symbol):
        raise NotImplementedError("Algorithms finding range for non-rational functions"
                                    "are not yet implemented")

    if domain.is_EmptySet:
        return EmptySet()

    if not func.has(symbol):
        return FiniteSet(func)

    # this block of code can be replaced by
    # sing = Intersection(FiniteSet(*sing), domain.closure)
    # after the issue #9706 has been fixed
    def closure_handle(set_im, singul):
        if set_im.has(Interval):
            if not oo in set_im.boundary:
                if not S.NegativeInfinity in set_im.boundary:
                    return Intersection(FiniteSet(*singul), set_im.closure)
                return Intersection(FiniteSet(*singul), Union(set_im,
                                    FiniteSet(max(set_im.boundary))))
            else:
                if not S.NegativeInfinity in set_im.boundary:
                    return Intersection(FiniteSet(*singul), Union(set_im,
                                        FiniteSet(min(set_im.boundary))))
                return Intersection(FiniteSet(*singul), set_im)
        return Intersection(FiniteSet(*singul), set_im)

    # all the singularities of the function
    sing = singularities(func, symbol)
    sing_in_domain = closure_handle(domain, sing)

    def codomain_interval(f, set_val, *sym):
        symb = sym[0]
        df1 = diff(f, symb)
        df2 = diff(df1, symb)
        der_zero = solveset_real(df1, symb)
        der_zero_in_dom = closure_handle(set_val, der_zero)

        local_maxima = set()
        local_minima = set()
        start_val = limit(f, symb, set_val.start)
        end_val = limit(f, symb, set_val.end, '-')

        for i in sing_in_domain:
            if not i in set_val.boundary:
                return Union(codomain(f, Interval(set_val.start, i, set_val.left_open, True), symb),
                            codomain(f, Interval(i, set_val.end, True, set_val.right_open), symb))

        if start_val is S.Infinity or end_val is S.Infinity:
            local_maxima = set([(oo, True)])
        elif start_val is S.NegativeInfinity or end_val is S.NegativeInfinity:
            local_minima = set([(-oo, True)])

        if local_maxima == set():
            if start_val > end_val:
                local_maxima = set([(start_val, set_val.left_open)])
            elif start_val < end_val:
                local_maxima = set([(end_val, set_val.right_open)])
            else:
                local_maxima = set([(start_val, set_val.left_open and set_val.right_open)])

        if local_minima == set():
            if start_val < end_val:
                local_minima = set([(start_val, set_val.left_open)])
            elif start_val > end_val:
                local_minima = set([(end_val, set_val.right_open)])
            else:
                local_minima = set([(start_val, set_val.left_open and set_val.right_open)])

        for i in der_zero_in_dom:
            exist = not i in set_val
            if df2.subs({symb: i}) < 0:
                if not i in sing_in_domain:
                    local_maxima.add((f.subs({symb: i}), exist))
                else:
                    local_maxima.add((oo, True))
            elif df2.subs({symb: i}) > 0:
                if not i in sing_in_domain:
                    local_minima.add((f.subs({symb: i}), exist))
                else:
                    local_minima.add((-oo, True))

        maximum = (-oo, True)
        minimum = (oo, True)

        for i in local_maxima:
            if i[0] > maximum[0]:
                maximum = i
            elif i[0] == maximum[0]:
                maximum = (maximum[0], i[1] and maximum[1])

        for i in local_minima:
            if i[0] < minimum[0]:
                minimum = i
            elif i[0] == minimum[0]:
                minimum = (minimum[0], i[1] and minimum[1])

        return Union(Interval(minimum[0], maximum[0], minimum[1], maximum[1]))

    if isinstance(domain, Union):
        return Union(*[codomain(func, intrvl_or_finset, symbol) for intrvl_or_finset in domain.args])

    if isinstance(domain, Interval):
        return codomain_interval(func, domain, symbol)

    if isinstance(domain, FiniteSet):
        domain = Complement(domain, FiniteSet(*sing_in_domain))
        return FiniteSet(*[limit(func, symbol, i) if i in FiniteSet(-oo, oo)
                            else func.subs({symbol: i}) for i in domain])
