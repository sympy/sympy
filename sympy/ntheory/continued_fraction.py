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

    If the discriminator is zero then the number will be a rational number:

    >>> continued_fraction_periodic(4, 3, 0)
    [1, 3]

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

    if q == 0:
        raise ValueError("The denominator is zero.")

    if d < 0:
        raise ValueError("Delta supposed to be a non-negative "
                         "integer, got %d" % d)
    elif d == 0:
        # the number is a rational number
        return list(continued_fraction_iterator(Rational(p, q)))

    if (d - p**2)%q:
        d *= q**2
        p *= abs(q)
        q *= abs(q)

    terms = []
    pq = {}
    sd = sqrt(d)

    while (p, q) not in pq:
        pq[(p, q)] = len(terms)
        terms.append(int((p + sd)/q))
        p = terms[-1]*q - p
        q = (d - p**2)/q

    i = pq[(p, q)]
    return terms[:i] + [terms[i:]]


def continued_fraction_iterator(x):
    """
    Return continued fraction expansion of x as iterator.

    Examples
    ========

    >>> from sympy.core import Rational, pi
    >>> from sympy.ntheory.continued_fraction import continued_fraction_iterator

    >>> list(continued_fraction_iterator(Rational(3, 8)))
    [0, 2, 1, 2]

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

    while True:
        i = Integer(x)
        yield i
        x -= i
        if not x:
            break
        x = 1/x
