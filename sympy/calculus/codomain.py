from sympy import S, sympify, diff, limit, oo, Dummy
from sympy.calculus.singularities import singularities
from sympy.sets.sets import Interval, Intersection, FiniteSet, Union, Complement, Set, EmptySet
from sympy.solvers.solveset import solveset


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

    # TODO: handle piecewise defined functions
    # TODO: handle transcendental functions
    # TODO: handle multivariate functions
    if len(syms) == 0:
        raise ValueError("A Symbol or a tuple of symbols must be given")

    if len(syms) == 1:
        symbol = syms[0]
    else:
        raise NotImplementedError("more than one variables %s not handled" % (syms,))

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
    sing = solveset(func.as_numer_denom()[1], symbol, domain=S.Reals)
    sing_in_domain = closure_handle(domain, sing)
    domain = Complement(domain, sing_in_domain)

    if domain.is_EmptySet:
        return EmptySet()

    def codomain_interval(f, set_val, *sym):
        symb = sym[0]
        df1 = diff(f, symb)
        df2 = diff(df1, symb)
        der_zero = solveset(df1, symb, domain=S.Reals)
        der_zero_in_dom = closure_handle(set_val, der_zero)

        local_maxima = set()
        local_minima = set()
        start_val = limit(f, symb, set_val.start)
        end_val = limit(f, symb, set_val.end, '-')

        if start_val is S.Infinity or end_val is S.Infinity:
            local_maxima = set([(oo, True)])
        elif start_val is S.NegativeInfinity or end_val is S.NegativeInfinity:
            local_minima = set([(-oo, True)])

        if (not start_val.is_real) or (not end_val.is_real):
            raise ValueError('Function does not contain all points of %s '
                            'as its domain' % (domain))

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
                local_maxima.add((f.subs({symb: i}), exist))
            elif df2.subs({symb: i}) > 0:
                local_minima.add((f.subs({symb: i}), exist))

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
        return FiniteSet(*[limit(func, symbol, i) if i in FiniteSet(-oo, oo)
                            else func.subs({symbol: i}) for i in domain])


def not_empty_in(finset_intersection, *syms):
    """ Finds the domain of the functions in `finite_set` in which the
    `finite_set` is not-empty

    Parameters
    ==========

    finset_intersection: The unevaluated intersection of FiniteSet containing
                        real-valued functions with Union of Sets
    syms: Tuple of symbols
            Symbol for which domain is to be found

    Raises
    ======

    NotImplementedError
        The algorithms to find the non-emptiness of the given FiniteSet are
        not yet implemented.
    ValueError
        The input is not valid.
    RuntimeError
        It is a bug, please report to the github issue tracker.

    Examples
    ========

    >>> from sympy import Symbol, codomain, FiniteSet, Interval, sqrt, oo
    >>> from sympy.calculus.codomain import not_empty_in
    >>> x = Symbol('x')
    >>> not_empty_in(FiniteSet(x/2).intersect(Interval(0, 1)), x)
    [0, 2]
    >>> not_empty_in(FiniteSet(x, x**2).intersect(Interval(1, 2)), x)
    [-sqrt(2), -1] U [1, 2]
    >>> not_empty_in(FiniteSet(x**2/(x + 2)).intersect(Interval(1, oo)), x)
    (-2, -1] U [2, oo)
    """

    # TODO: handle piecewise defined functions
    # TODO: handle transcendental functions
    # TODO: handle multivariate functions
    if len(syms) == 0:
        raise ValueError("A Symbol or a tuple of symbols must be given \
                            as the third parameter")
    if isinstance(finset_intersection, FiniteSet):
        finite_set = finset_intersection
        sets = S.Reals
    elif isinstance(finset_intersection, Union):
        elm_in_sets = finset_intersection.args[0]
        return Union(not_empty_in(finset_intersection.args[1], *syms), elm_in_sets)
    else:
        finite_set = finset_intersection.args[1]
        sets = finset_intersection.args[0]
    if not isinstance(finite_set, FiniteSet):
        raise ValueError('A FiniteSet must be given, not %s: %s' % (type(finite_set), finite_set))
    if len(syms) == 1:
        symbol = syms[0]
    else:
        raise NotImplementedError('more than one variables %s not handled' % (syms,))

    y = Dummy('y')

    def elm_domain(expr, codomain_expr, *sym):
        """
        finds the domain of a single element of `finite_set` with `codomain_expr`
        as its codomain and with given symbols `syms`
        """
        domain_union = S.EmptySet

        # find the inverse of items in the finite_set
        invert_expr = solveset(expr - y, symbol, domain=S.Reals)

        if isinstance(invert_expr, Intersection):
            invert_set = invert_expr.args[1]
        elif isinstance(invert_expr, Complement):
            # remove the elements which are not in `domain` of `syms`
            exclude_sing = invert_expr.args[1]
            codomain_expr = Complement(codomain_expr, exclude_sing)
            if isinstance(invert_expr.args[0], Intersection):
                invert_set = invert_expr.args[0].args[1]
            elif isinstance(invert_expr.args[0], FiniteSet):
                invert_set = invert_expr.args[0]
        else:
            raise NotImplementedError("The algorithm to find the domain for the function %s \
                                        are not yet implemented." % (expr))

        # find the domain of `every` inverse of `expr`
        # since expr can have more than one inverse
        for inverse_val in invert_set:
            domain = codomain(inverse_val, codomain_expr, y)
            domain_union = Union(domain_union, domain)
        return domain_union

    return Union(*[elm_domain(element, sets, symbol) for element in finite_set])
