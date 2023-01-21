from __future__ import annotations
from sympy.core.numbers import igcdex


def migcdex(*x: int) -> tuple[tuple[int, ...], int]:
    r"""Compute extended gcd for multiple integers.

    Explanation
    ===========

    Given the integers $x_1, \cdots, x_n$ and
    an extended gcd for multiple arguments are defined as a solution
    $(y_1, \cdots, y_n), g$ for the diophantine equation
    $x_1 y_1 + \cdots + x_n y_n = g$ such that
    $g = \gcd(x_1, \cdots, x_n)$.

    Examples
    ========

    >>> from sympy.functions.elementary._trig_sqrt import igcdex_multi
    >>> igcdex_multi()
    ((), 0)
    >>> igcdex_multi(4)
    ((1,), 4)
    >>> igcdex_multi(4, 6)
    ((-1, 1), 2)
    >>> igcdex_multi(6, 10, 15)
    ((1, 1, -1), 1)
    """
    if not x:
        return (), 0

    if len(x) == 1:
        return (1,), x[0]

    if len(x) == 2:
        u, v, h = igcdex(x[0], x[1])
        return (u, v), h

    y, g = migcdex(*x[1:])
    u, v, h = igcdex(x[0], g)
    return (u, *(v * i for i in y)), h
