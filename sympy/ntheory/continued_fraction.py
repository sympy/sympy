from sympy.core.numbers import Integer, Rational


def continued_fraction_periodic(p, q, d=0):
    r"""
    Compute the continued fraction expansion of a rational or a
    quadratic irrational number, i.e. `\frac{p + \sqrt{d}}{q}`, where
    `p`, `q` and `d \ge 0` are integers.

    Returns the continued fraction representation (canonical form) as
    a list of integers, optionally ending (for quadratic irrationals)
    with repeating block as the last term of this list.

    Parameters
    ==========

    p : int
        the rational part of the number's numerator
    q : int
        the denominator of the number
    d : int, optional
        the irrational part (discriminator) of the number's numerator

    Examples
    ========

    >>> from sympy.ntheory.continued_fraction import continued_fraction_periodic
    >>> continued_fraction_periodic(3, 2, 7)
    [2, [1, 4, 1, 1]]

    Golden ratio has the simplest continued fraction expansion:

    >>> continued_fraction_periodic(1, 2, 5)
    [[1]]

    If the discriminator is zero or a perfect square then the number will be a
    rational number:

    >>> continued_fraction_periodic(4, 3, 0)
    [1, 3]
    >>> continued_fraction_periodic(4, 3, 49)
    [3, 1, 2]

    See Also
    ========

    continued_fraction_iterator

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Periodic_continued_fraction
    .. [2] K. Rosen. Elementary Number theory and its applications.
           Addison-Wesley, 3 Sub edition, pages 379-381, January 1992.

    """
    from sympy.core.compatibility import as_int
    from sympy.functions import sqrt

    p, q, d = list(map(as_int, [p, q, d]))

    sd = sqrt(d)
    if q == 0:
        raise ValueError("The denominator is zero.")

    if d < 0:
        raise ValueError("Delta supposed to be a non-negative "
                         "integer, got %d" % d)
    elif d == 0 or sd.is_integer:
        # the number is a rational number
        return list(continued_fraction_iterator(Rational(p + sd, q)))

    if (d - p**2)%q:
        d *= q**2
        sd *= q
        p *= abs(q)
        q *= abs(q)

    terms = []
    pq = {}

    while (p, q) not in pq:
        pq[(p, q)] = len(terms)
        terms.append(int((p + sd)/q))
        p = terms[-1]*q - p
        q = (d - p**2)/q

    i = pq[(p, q)]
    return terms[:i] + [terms[i:]]

def continued_fraction_quadratic(cf):
    """
    Compute the rational or quadratic irrational number from its
    terminating or periodic continued fraction expansion.  The continued
    fraction expansion (cf) should be supplied as a simple list (for rational
    numbers, in which case this function is the same as
    continued_fraction_rational), or a list of the non-repeating terms (if
    any) with a list of the repeating terms as the last element of the
    list.  This is the format returned by continued_fraction_periodic.

    Returns the largest solution found, which is generally the one sought,
    if the fraction is in canonical form (all terms positive except
    possibly the first).

    Examples:
    =========

    >>> from sympy.ntheory.continued_fraction import continued_fraction_quadratic
    >>> continued_fraction_quadratic([1, 4, 2, [3, 1]])
    (sqrt(21) + 287)/238
    >>> continued_fraction_quadratic([[1]])
    1/2 + sqrt(5)/2
    >>> from sympy.ntheory.continued_fraction import continued_fraction_periodic
    >>> continued_fraction_quadratic(continued_fraction_periodic(8, 5, 13))
    (sqrt(13) + 8)/5

    """

    from sympy.core.symbol import Dummy
    from sympy.solvers import solve

    if not isinstance(cf[-1], list):
        return continued_fraction_rational(cf)

    x=Dummy('x')
    solns=solve(continued_fraction_rational(cf[-1]+[x])-x, x)
    solns.sort()
    pure=solns[-1]
    if cf[:-1]:
        return continued_fraction_rational(cf[:-1]+[x]).subs(x, pure).radsimp()
    else:
        return pure.radsimp()

def continued_fraction_periodic_iterator(p, q, d=0):
    r"""
    Compute the continued fraction expansion of a rational or a
    quadratic irrational number, i.e. `\frac{p + \sqrt{d}}{q}`, where
    `p`, `q` and `d \ge 0` are integers.

    Returns an iterator which yields the continued fraction expansion one
    term at a time.  Note that if `d` is not zero or a perfect square the
    iterator will never terminate.

    Parameters
    ==========

    p : int
        the rational part of the number's numerator
    q : int
        the denominator of the number
    d : int, optional
        the irrational part (discriminator) of the number's numerator

    Examples:
    =========

    >>> from sympy.ntheory.continued_fraction import continued_fraction_periodic_iterator
    >>> from itertools import islice
    >>> it=continued_fraction_periodic_iterator(4, 3, 0)
    >>> list(it)
    [1, 3]
    >>> it=continued_fraction_periodic_iterator(3, 2, 3)
    >>> list(islice(it, 7))
    [2, 2, 1, 2, 1, 2, 1]

    """

    from itertools import cycle, chain

    expansion = continued_fraction_periodic(p, q, d)
    if not isinstance(expansion[-1], list):
        return iter(expansion)

    return chain(iter(expansion[:-1]), cycle(expansion[-1]))


def continued_fraction_iterator(x):
    """
    Return continued fraction expansion of x as iterator.

    Examples
    ========

    >>> from sympy.core import Rational, pi
    >>> from sympy.ntheory.continued_fraction import continued_fraction_iterator

    >>> list(continued_fraction_iterator(Rational(3, 8)))
    [0, 2, 1, 2]
    >>> list(continued_fraction_iterator(Rational(-3, 8)))
    [-1, 1, 1, 1, 2]

    >>> for i, v in enumerate(continued_fraction_iterator(pi)):
    ...    if i > 7:
    ...        break
    ...    print(v)
    3
    7
    15
    1
    292
    1
    1
    1

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Continued_fraction

    """

    from sympy.functions import floor

    while True:
        i = floor(x)
        yield i
        x -= i
        if not x:
            break
        x = 1/x

def continued_fraction_convergents(cf):
    """Return an iterator over the convergents of a continued fraction (cf).

    The parameter should be an iterable returning successive
    partial quotients of the continued fraction, such as might be
    returned by continued_fraction_iterator.  In computing the
    convergents, the continued fraction need not be strictly in
    canonical form (all integers, all but the first positive).
    Rational and negative elements may be present in the expansion.

    Examples
    ========

    >>> from sympy.core import Rational, pi
    >>> from sympy import S
    >>> from sympy.ntheory.continued_fraction import \
            continued_fraction_convergents, continued_fraction_iterator

    >>> list(continued_fraction_convergents([0, 2, 1, 2]))
    [0, 1/2, 1/3, 3/8]

    >>> list(continued_fraction_convergents([1, S('1/2'), -7, S('1/4')]))
    [1, 3, 19/5, 7]

    >>> it=continued_fraction_convergents(continued_fraction_iterator(pi))
    >>> for n in range(7):
    ...  print(next(it))
    3
    22/7
    333/106
    355/113
    103993/33102
    104348/33215
    208341/66317

    """

    from sympy import S

    p_2, q_2 = 0, 1
    p_1, q_1 = 1, 0
    for a in cf:
        p, q = a*p_1 + p_2, a*q_1 + q_2
        p_2, q_2 = p_1, q_1
        p_1, q_1 = p, q
        yield S(p)/S(q)

def continued_fraction_rational(cf):
    """Return the Rational from a continued fraction iterator (cf).

    Reduce the continued fraction iterator to a Rational object.  This
    is equivalent to list(continued_fraction_convergents(cf))[-1], but
    more efficient.  Be careful about supplying a possibly infinite
    iterator.

    Examples:
    ========

    >>> from sympy.ntheory import continued_fraction_rational
    >>> continued_fraction_rational([1, 2, 3, 4, 5])
    225/157
    >>> continued_fraction_rational([-2, 1, 9, 7, 1, 2])
    -256/233
    >>> continued_fraction_rational([2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8]).n(10)
    2.718281835

    """

    it=continued_fraction_convergents(cf)
    a=Integer(0)
    for a in it:
        pass
    return a
