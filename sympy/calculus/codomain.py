from sympy import S, sympify, diff, limit, oo
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
          It is a bug, please report it to the github issue tracker
          (https://github.com/sympy/sympy/issues).

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
        sym = syms[0]
    else:
        raise NotImplementedError("more than one variables %s not handled" % (syms,))

    if not func.has(sym):
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
    sing = solveset(func.as_numer_denom()[1], sym, domain=S.Reals)
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
        return Union(*[codomain(func, intrvl_or_finset, sym) for intrvl_or_finset in domain.args])

    if isinstance(domain, Interval):
        return codomain_interval(func, domain, sym)

    if isinstance(domain, FiniteSet):
        return FiniteSet(*[limit(func, sym, i) if i in FiniteSet(-oo, oo)
                            else func.subs({sym: i}) for i in domain])


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
        It is a bug, please report it to the github issue tracker
        (https://github.com/sympy/sympy/issues).

    Examples
    ========

    >>> from sympy import Symbol, FiniteSet, Interval, sqrt, oo
    >>> from sympy.calculus.codomain import not_empty_in
    >>> from sympy.abc import x
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

    if finset_intersection.is_EmptySet:
        return EmptySet()

    if isinstance(finset_intersection, Union):
        elm_in_sets = finset_intersection.args[0]
        return Union(not_empty_in(finset_intersection.args[1], *syms), elm_in_sets)

    if isinstance(finset_intersection, FiniteSet):
        finite_set = finset_intersection
        _sets = S.Reals
    else:
        finite_set = finset_intersection.args[1]
        _sets = finset_intersection.args[0]

    if not isinstance(finite_set, FiniteSet):
        raise ValueError('A FiniteSet must be given, not %s: %s' % (type(finite_set), finite_set))

    if len(syms) == 1:
        symb = syms[0]
    else:
        raise NotImplementedError('more than one variables %s not handled' % (syms,))

    def elm_domain(expr, intrvl):
        """ Finds the domain of an expression in any given interval """

        _start = intrvl.start
        _end = intrvl.end
        _sing = solveset(expr.as_numer_denom()[1], symb, domain=S.Reals)

        if intrvl.right_open:
            if _end is S.Infinity:
                _sing1 = solveset(expr.as_numer_denom()[1], symb, domain=S.Reals)
                _domain1 = Complement(S.Reals, _sing1)
            else:
                _domain1 = solveset(expr < _end, symb, domain=S.Reals)
        else:
            _domain1 = solveset(expr <= _end, symb, domain=S.Reals)

        if intrvl.left_open:
            if _start is S.NegativeInfinity:
                _sing2 = solveset(expr.as_numer_denom()[1], symb, domain=S.Reals)
                _domain2 = Complement(S.Reals, _sing2)
            else:
                _domain2 = solveset(expr > _start, symb, domain=S.Reals)
        else:
            _domain2 = solveset(expr >= _start, symb, domain=S.Reals)

        # domain in the interval
        expr_with_sing = Intersection(_domain1, _domain2)
        expr_domain = Complement(expr_with_sing, _sing)
        return expr_domain

    if isinstance(_sets, Interval):
        return Union(*[elm_domain(element, _sets) for element in finite_set])

    if isinstance(_sets, Union):
        _dom_final = S.EmptySet
        for intrvl in _sets.args:
            _domain_element = Union(*[elm_domain(element, intrvl) for element in finite_set])
            _dom_final = Union(_dom_final, _domain_element)
        return _dom_final
