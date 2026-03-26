from __future__ import annotations
import itertools
from typing import TYPE_CHECKING, Iterable

from sympy.core.exprtools import factor_terms
from sympy.core.numbers import Integer, Rational
from sympy.core.singleton import S
from sympy.core.symbol import Dummy
from sympy.core.sympify import _sympify
from sympy.utilities.misc import as_int

if TYPE_CHECKING:
    from sympy.core.expr import Expr

def continued_fraction(a: object) -> list:
    """Return the continued fraction representation of a Rational or
    quadratic irrational.

    Examples
    ========

    >>> from sympy.ntheory.continued_fraction import continued_fraction
    >>> from sympy import sqrt
    >>> continued_fraction((1 + 2*sqrt(3))/5)
    [0, 1, [8, 3, 34, 3]]

    See Also
    ========
    continued_fraction_periodic, continued_fraction_reduce, continued_fraction_convergents
    """
    e = _sympify(a)
    if all(i.is_Rational for i in e.atoms()):
        if e.is_Integer:
            return continued_fraction_periodic(e, 1, 0)
        elif e.is_Rational:
            # FIX: narrow type to Rational before accessing .p and .q
            assert isinstance(e, Rational)
            return continued_fraction_periodic(e.p, e.q, 0)
        elif e.is_Pow:
            # FIX: narrow type before accessing .exp and .base
            from sympy.core.power import Pow
            assert isinstance(e, Pow)
            if e.exp is S.Half and e.base.is_Integer:
                return continued_fraction_periodic(0, 1, e.base)
        elif e.is_Mul and len(e.args) == 2 and (
                e.args[0].is_Rational and
                e.args[1].is_Pow and
                e.args[1].base.is_Integer and
                e.args[1].exp is S.Half):
            a_arg, b_arg = e.args
            assert isinstance(a_arg, Rational)
            from sympy.core.power import Pow
            assert isinstance(b_arg, Pow)
            return continued_fraction_periodic(0, a_arg.q, b_arg.base, a_arg.p)
        else:
            p, d = e.expand().as_numer_denom()
            if d.is_Integer:
                if p.is_Rational:
                    return continued_fraction_periodic(p, d)
                if p.is_Add and len(p.args) == 2:
                    a_val, bc = p.args
                else:
                    a_val = S.Zero
                    bc = p
                if a_val.is_Integer:
                    b_val: Expr = S.NaN
                    c_val: Expr = S.NaN
                    if bc.is_Mul and len(bc.args) == 2:
                        b_val, c_val = bc.args  # type: ignore[assignment]
                    elif bc.is_Pow:
                        b_val = Integer(1)
                        c_val = bc
                    from sympy.core.power import Pow
                    if b_val.is_Integer and (
                            isinstance(c_val, Pow) and c_val.exp is S.Half and
                            c_val.base.is_Integer):
                        # (a + b*sqrt(c))/d
                        c_base = c_val.base
                        return continued_fraction_periodic(a_val, d, c_base, b_val)
    raise ValueError(f"expecting a rational or quadratic irrational, not {e}")


def continued_fraction_periodic(p: int | Expr, q: int | Expr, d: int | Expr = 0, s: int | Expr = 1) -> list:
    r"""
    Find the periodic continued fraction expansion of a quadratic irrational.

    Compute the continued fraction expansion of a rational or a
    quadratic irrational number, i.e. `\frac{p + s\sqrt{d}}{q}`, where
    `p`, `q \ne 0` and `d \ge 0` are integers.

    Returns the continued fraction representation (canonical form) as
    a list of integers, optionally ending (for quadratic irrationals)
    with list of integers representing the repeating digits.

    Parameters
    ==========

    p : int
        the rational part of the number's numerator
    q : int
        the denominator of the number
    d : int, optional
        the irrational part (discriminator) of the number's numerator
    s : int, optional
        the coefficient of the irrational part

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

    continued_fraction_iterator, continued_fraction_reduce

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Periodic_continued_fraction
    .. [2] K. Rosen. Elementary Number theory and its applications.
           Addison-Wesley, 3 Sub edition, pages 379-381, January 1992.

    """
    from sympy.functions import sqrt, floor

    p, q, d, s = int(p), int(q), int(d), int(s)

    if d < 0:
        raise ValueError(f"expected non-negative for `d` but got {d}")

    if q == 0:
        raise ValueError("The denominator cannot be 0.")

    if not s:
        d = 0

    # check for rational case
    sd = sqrt(d)
    if sd.is_Integer:
        return list(continued_fraction_iterator(Rational(p + s*sd, q)))

    # irrational case with sd != Integer
    if q < 0:
        p, q, s = -p, -q, -s

    n = (p + s*sd)/q
    if n.is_negative:
        w = floor(-n)
        f = -n - w
        one_f = continued_fraction(1 - f)  # 1-f < 1 so cf is [0 ... [...]]
        one_f[0] -= w + 1
        return one_f

    d *= s**2
    sd *= s

    if (d - p**2)%q:
        d *= q**2
        sd *= q
        p *= q
        q *= q

    terms: list[int] = []  # FIX: use list[int] instead of list[Expr]
    pq: dict[tuple[int, int], int] = {}

    while (p, q) not in pq:
        pq[(p, q)] = len(terms)
        terms.append(int((p + sd)//q))  # FIX: cast to int explicitly
        p = terms[-1]*q - p
        q = (d - p**2)//q

    i = pq[(p, q)]
    return terms[:i] + [terms[i:]]  # type: ignore


def continued_fraction_reduce(cf: object) -> Expr:
    """
    Reduce a continued fraction to a rational or quadratic irrational.

    Compute the rational or quadratic irrational number from its
    terminating or periodic continued fraction expansion.  The
    continued fraction expansion (cf) should be supplied as a
    terminating iterator supplying the terms of the expansion.  For
    terminating continued fractions, this is equivalent to
    ``list(continued_fraction_convergents(cf))[-1]``, only a little more
    efficient.  If the expansion has a repeating part, a list of the
    repeating terms should be returned as the last element from the
    iterator.  This is the format returned by
    continued_fraction_periodic.

    For quadratic irrationals, returns the largest solution found,
    which is generally the one sought, if the fraction is in canonical
    form (all terms positive except possibly the first).

    Examples
    ========

    >>> from sympy.ntheory.continued_fraction import continued_fraction_reduce
    >>> continued_fraction_reduce([1, 2, 3, 4, 5])
    225/157
    >>> continued_fraction_reduce([-2, 1, 9, 7, 1, 2])
    -256/233
    >>> continued_fraction_reduce([2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8]).n(10)
    2.718281835
    >>> continued_fraction_reduce([1, 4, 2, [3, 1]])
    (sqrt(21) + 287)/238
    >>> continued_fraction_reduce([[1]])
    (1 + sqrt(5))/2
    >>> from sympy.ntheory.continued_fraction import continued_fraction_periodic
    >>> continued_fraction_reduce(continued_fraction_periodic(8, 5, 13))
    (sqrt(13) + 8)/5

    See Also
    ========

    continued_fraction_periodic

    """
    from sympy.solvers import solve

    period = []
    x = Dummy('x')

    def untillist(cf):
        for nxt in cf:
            if isinstance(nxt, list):
                period.extend(nxt)
                yield x
                break
            yield nxt

    a = S.Zero
    for a in continued_fraction_convergents(untillist(cf)):
        pass

    if period:
        y = Dummy('y')
        solns = solve(continued_fraction_reduce(period + [y]) - y, y)
        solns.sort()
        pure = solns[-1]
        rv = a.subs(x, pure).radsimp()
    else:
        rv = a
    if rv.is_Add:
        rv = factor_terms(rv)
        if rv.is_Mul and rv.args[0] == -1:
            rv = rv.func(*rv.args)
    return rv


def continued_fraction_iterator(x: Expr):
    """
    Return continued fraction expansion of x as iterator.

    Examples
    ========

    >>> from sympy import Rational, pi
    >>> from sympy.ntheory.continued_fraction import continued_fraction_iterator

    >>> list(continued_fraction_iterator(Rational(3, 8)))
    [0, 2, 1, 2]
    >>> list(continued_fraction_iterator(Rational(-3, 8)))
    [-1, 1, 1, 1, 2]

    >>> for i, v in enumerate(continued_fraction_iterator(pi)):
    ...     if i > 7:
    ...         break
    ...     print(v)
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

    .. [1] https://en.wikipedia.org/wiki/Continued_fraction

    """
    from sympy.functions import floor
    while True:
        i = floor(x)
        yield i
        x -= i
        if not x:
            break
        x = 1/x


def continued_fraction_convergents(cf: Iterable):
    """
    Return an iterator over the convergents of a continued fraction (cf).

    The parameter should be in either of the following to forms:
    - A list of partial quotients, possibly with the last element being a list
    of repeating partial quotients, such as might be returned by
    continued_fraction and continued_fraction_periodic.
    - An iterable returning successive partial quotients of the continued
    fraction, such as might be returned by continued_fraction_iterator.

    In computing the convergents, the continued fraction need not be strictly
    in canonical form (all integers, all but the first positive).
    Rational and negative elements may be present in the expansion.

    Examples
    ========

    >>> from sympy.core import pi
    >>> from sympy import S
    >>> from sympy.ntheory.continued_fraction import \
            continued_fraction_convergents, continued_fraction_iterator

    >>> list(continued_fraction_convergents([0, 2, 1, 2]))
    [0, 1/2, 1/3, 3/8]

    >>> list(continued_fraction_convergents([1, S('1/2'), -7, S('1/4')]))
    [1, 3, 19/5, 7]

    >>> it = continued_fraction_convergents(continued_fraction_iterator(pi))
    >>> for n in range(7):
    ...     print(next(it))
    3
    22/7
    333/106
    355/113
    103993/33102
    104348/33215
    208341/66317

    >>> it = continued_fraction_convergents([1, [1, 2]])  # sqrt(3)
    >>> for n in range(7):
    ...     print(next(it))
    1
    2
    5/3
    7/4
    19/11
    26/15
    71/41

    See Also
    ========

    continued_fraction_iterator, continued_fraction, continued_fraction_periodic

    """
    if isinstance(cf, list) and isinstance(cf[-1], list):
        cf = itertools.chain(cf[:-1], itertools.cycle(cf[-1]))
    p_2: Expr = S.Zero
    q_2: Expr = S.One
    p_1: Expr = S.One
    q_1: Expr = S.Zero
    for a in cf:
        p, q = a*p_1 + p_2, a*q_1 + q_2
        p_2, q_2 = p_1, q_1
        p_1, q_1 = p, q
        yield p/q
